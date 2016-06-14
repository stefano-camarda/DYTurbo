#ifndef MERGER
#define MERGER

/**
 * @file merger.C
 * Description of this file
 *
 * @brief A brief description
 *
 * @author cuto <Jakub.Cuth@cern.ch>
 * @author Stefano <Stefano.Camarda@cern.ch>
 * @date 2015-11-18
 */
#include "isnan.h"

#include <vector>
#include <set>
#include <algorithm>
#include <iostream>
using namespace std;

#include <TFile.h>
#include <TH1.h>
#include <TH3.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TList.h>
#include <TClass.h>
#include <TKey.h>
#include <TString.h>
#include <TMath.h>
#include <TObjectTable.h>

//Force to kill
#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>

// CXX option parser: https://raw.githubusercontent.com/jarro2783/cxxopts/master/src/cxxopts.hpp
#include "cxxopts.hpp"
namespace po=cxxopts; // inspired by po = boost::program_options

#define PRINT_HELP(RTN) do {printf("%s\n",opts.help({""}).c_str()); return RTN ;} while (false)

typedef std::string SString;
typedef std::vector<SString> VecSStr;
typedef std::vector<TString> VecTStr;
typedef std::vector<TFile*> VecTFile;
typedef std::vector<TH1*> VecTH1;
typedef std::vector<double> VecDbl;
typedef std::vector<int> VecInt;
//using std::sort;

// Z pt binning from LHC 7TeV measurement
double bins[23] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 60, 70, 80, 100};

class OutlierRemoval{
    public :
        OutlierRemoval() :
            includeUnderOverFlow(true),
            doXsecNormalization(false),
            doRebin(false),
            doNanIterp(false),
            doMedian(false),
            doOutlierRemoval(false),
            do2D(false),
            doTest(false),
            doCuba(false),
            verbose(1) {};
        ~OutlierRemoval(){
            if (verbose >6) {
                printf ( " Destructor:\n");
                gDirectory->ls("-m");
                gObjectTable->Print();
            }
        };

        // Public methods
        // Main function:
        void Merge(){
            // Take first file and create list of object names
            find_all_objects();
            // NOTE: Previously I opened all files and load all objects. This
            // is RAM extensive for larger root files.
            if (doXsecNormalization) read_Xsection();
            if(verbose>0) {printf(" Merging .... "); fflush(stdout);}
            for (auto p_objname : all_obj_names){
                // setup names
                const char* objname = p_objname.Data();
                size_t len = p_objname.Length();
                char proj = p_objname(len-1); //< projection string
                if (verbose>1) printf("objname: %s\n",objname);
                // Get one object from all files
                VecTH1 in_objs;
                for (auto it_fn : infilenames){
                if (verbose>1) printf("filename: %s\n",it_fn.Data());
                    TFile * it_f =TFile::Open(it_fn.Data(),"READ");
                    /// @todo: test if they are all same binning
                    // Create uniq name to avoid "Potential memory leak" warnings.
                    TString retrieved_name(p_objname); 
                    retrieved_name+="__"; 
                    retrieved_name+=in_objs.size();
                    const char* objname_out = retrieved_name.Data();
                    if (verbose>3) printf(" objname out: %s\n",objname_out);
                    // Retrieve object
                    TH1 * o = (TH1*) it_f->Get(objname);
                    if (o!=0) {
                        TH1 * tmp = o;
                        o = (TH1*)  o->Clone(objname_out);
                        delete tmp;
                    } 
                    // If o==0 than we have wrong histogram or retrieved object might by projection. Test it:
                    if (o==0 && len>3 ){
                        // Extract basename: it should end with '_px' where 'x' might be one of: `zyuvn`
                        TString basename = p_objname(0,len-3); // remove "_px" at the end of the name
                        if (verbose>1) printf("basename: %s , proj %c , objname out: %s \n",basename.Data(), proj, objname_out);
                        // Get original, which will be used to create projections
                        TH1* hnd = (TH1*) it_f->Get(basename.Data());
                        if (is_empty(hnd,objname,objname_out)) continue;
                        if ( proj != 'n') {
                            // 2D -> X axis
                            o = make_projection(hnd,proj); // h2d->ProjectionX("dummy"); // ,-1,0,"e");
                            o->SetName(objname_out);
                            if (doXsecNormalization) normalize(o);
                        } else { // if (proj == 'n') 
                            //Rebin: this not projection but rebining to ptz to measuremnt binning
                            // load correct
                            basename = p_objname(0,len-6); // remove "_rebin" at the end of string
                            if(!hnd) delete hnd;
                            hnd = (TH1*) it_f->Get(basename.Data());
                            if (verbose>1) printf("basename: %s , proj %c \n",basename.Data(), proj);
                            // rebin
                            if (doXsecNormalization) normalize(hnd);
                            o=hnd->Rebin(22,objname_out,bins);
                            // divide by bin width
                            for (int ibin =1; ibin<=o->GetNbinsX(); ibin++ ){
                                o->SetBinContent(ibin,o->GetBinContent(ibin)/o->GetBinWidth(ibin));
                            }
                        }
                        delete hnd;
                        if (is_empty(o,objname,objname_out)) continue;
                    } else if (doXsecNormalization && !isProfile(o)) normalize(o);
                    if (o!=0) {
                        if (verbose>2) o->Print();
                        // Check for NaN and adding object to list
                        bool hasNaN= false;
                        if (doNanIterp){ // slower but precise
                            for (auto ibin: loop_bins(o)){
                                if (isnan_ofast(o->GetBinContent(ibin))) {
                                    double val = interpolate_NaN(o,ibin);
                                    printf("Warning: hist %s in file %s  NaN in bin %d was interpolated to %f \n", 
                                            it_fn.Data(), p_objname.Data(), ibin, val);
                                }
                            }
                        } else { // simple fast test: integral
                            hasNaN = isnan_ofast(o->Integral());
                        }
                        if (hasNaN) {
                            printf("Warning: hist %s in file %s contain NaN \n", it_fn.Data(), p_objname.Data());
                            delete o;
                        }
                        else in_objs.push_back(o);
                    } else printf("Warning: hist %s not in file %s \n", it_fn.Data(), p_objname.Data());
		    it_f->Close();
                    if (doTest && in_objs.size()>100) break;
                } //end loop all files
                // DIFFERENT MERGING TYPES
                TString name = p_objname; // this is name of output object (identical to first input object)
                if (doCuba){
                    // CUBA merging:
                    // Since cuba calculation will be run on kinematic
                    // sub-region, we need to merge the ouptut table
                    // (represented by TH2D) properly.
                    //
                    // 1. Per each file you find all used uniq bins and create histogram:
                    TH1* tmp_c = find_all_bins(in_objs);
                    if (tmp_c==0){ continue; }
                    // 2. Sum correct bin content to correct bin (by bin center)
                    sum_bins(tmp_c,in_objs);
                    // Put to output list and go to next object (*no* averaging/median/outlier removal)
                    tmp_c->SetName(name);
                    output_objects.push_back(tmp_c);
                    delete_objs(in_objs); continue;
                }
                // prepare average mediand and profile object
                TH1* o_average = clone_empty(in_objs[0],(name+"_average").Data());
                int dim = o_average->GetDimension();
                TH1* o_median = 0;
                TH1* o_profile =0;
                if (isProfile(o_average)){
                    if (verbose>2) printf(" do profile \n");
                    o_profile = clone_empty(in_objs[0],(name+"_total").Data());
                }
                if (doMedian||doOutlierRemoval){
                    o_median = (dim==1||do2D) ? clone_empty(in_objs[0],(name+"_median").Data()) : 0;
                    // First Loop :  Get average and median
                    if (verbose>1) printf("    First Loop \n");
                    if (dim==1||do2D) {
                        // This part is same for histograms and profiles
                        create_average_obj(o_median,in_objs,"median");
                        if (verbose>3) {printf("  median average "); print_range(o_median);}
                    }
                }
                if (o_profile==0){
                    // This is not for profile!
                    // NOTE: don't calculate average for profiles, profiles should be summed up
                    create_average_obj(o_average,in_objs,"mean");
                    // save average and median
                    o_average->SetName(name);
                    output_objects.push_back(o_average);
                    if(doMedian||doOutlierRemoval) {
                        //
                        if (dim==1||do2D) {
                            o_median->SetName((name+"_median").Data());
                            output_objects.push_back(o_median);
                            if (verbose>2) printf("saving histogram %s with integral %f\n", o_median->GetName(), o_median->Integral());
                        }
                        if(doOutlierRemoval){
                            // Second loop -- outlier removal for histogram
                            o_average = clone_empty(in_objs[0],(name+"_average").Data());
                            // FIXME: Remove outliers with PDF correlated way
                            remove_outliers_from_hist(o_average,o_median,in_objs);
                            o_average->SetName((name+"_outliers").Data());
                            output_objects.push_back(o_average);
                        }
                    }
                    // STOP here for TH1,2,3
                    delete_objs(in_objs);
                } else {
                    // Profiles only -- total profile with a priori uncertainty
                    // NOTE: Adding up profiles would be enough for profiles,
                    // we are extra dividing by number of inputp objects to
                    // retrive correct xsection in denominator.
                    //
                    for (auto ith : in_objs) o_profile->Add(ith,1./in_objs.size());
                    // save profile and median
                    o_profile->SetName(name);
                    output_objects.push_back(o_profile);
                    if (verbose>2) printf("saving histogram %s with integral %f\n",o_profile->GetName(), o_profile->Integral());
                    if(doMedian||doOutlierRemoval){
                        if (dim>1&&!do2D) {delete_objs(in_objs); continue;}
                        o_median->SetName((name+"_median").Data());
                        output_objects.push_back(o_median);
                        // new median of entries
                        o_median = clone_empty(in_objs[0], (name+"_entries_median").Data() );
			cout << "1 " << o_median->ClassName() << endl;
                        create_average_obj(o_median,in_objs,"median_entries");
			cout << " 2 " << o_median->ClassName() << endl;
                        o_median->SetName((name)+"_entries_median");
                        if (verbose>2) printf("saving histogram %s with integral %f\n",o_median->GetName(), o_median->Integral());
                        output_objects.push_back(o_median);
                        if(doOutlierRemoval){
                            // Second Loop -- discard outliers
                            if (verbose>1) printf("    Second Loop \n");
			    cout << o_median->ClassName() << endl;
                            remove_outliers_from_profile(o_average, o_median, in_objs);
                            if (o_average==0){  delete_objs(in_objs); continue;}; // if not implemented 
                            o_average->SetName((name+"_outlier").Data());
                            output_objects.push_back(o_average);
                            if (verbose>2) {printf("  outliers removed \n"); print_range(o_average); print_range(o_profile);}
                        }
                    }
                    // release input histograms
                    delete_objs(in_objs);
                }
                // just for sure
                delete_objs(in_objs);
            }
            if (verbose>1) printf("End of merge, closing files \n");
	    //            close_all_infiles();
        };

        void Write(){
            if (verbose>1) printf("Write \n");
            TFile * fout = new TFile(outfilename.Data(), "RECREATE");
            for (auto ith : output_objects){
                if (verbose>4) printf("  writing pointer %p \n", ith);
                if (verbose>3) ith->Print("range");
                ith->Write();
            }
            if(verbose>0){
                printf("Merged %d files in %s, closing file ...", (int) infilenames.size(), outfilename.Data());
                fflush(stdout);
            }
            fout->Close();
            if(verbose>0) printf(" file closed.\n");
        };

        void SetOutputFile(TString _outf){
            outfilename=_outf;
            if(verbose>1)printf ("merger Target file: %s\n", _outf.Data()); // hadd like comments
        }

        void AddInputFile(TString _inf){
            infilenames.push_back(_inf); 
            if(verbose>2) printf ("merger Source file %d: %s\n", (int) infilenames.size(), _inf.Data()); // hadd like comments
        }

        void Init(){
  	    TH1::AddDirectory(kFALSE);
	    LOBIN= includeUnderOverFlow ? 0 : 1; // start from 0
            HIBIN= includeUnderOverFlow ? 1 : 0; // add 1 to end
            if (verbose >6) {
                printf ( " Init:\n");
                gDirectory->ls("-m");
                gObjectTable->Print();
            }
        }

        // Public data members
        bool includeUnderOverFlow;
        bool doXsecNormalization;
        bool do2D;
        bool doRebin;
        bool doNanIterp;
        bool doMedian;
        bool doOutlierRemoval;
        bool doCuba;
        bool doTest;
        int verbose;

    private :
        // Private methods
        void find_all_objects(){
            // Get first file and create list of object names
            // FIXME: Scan all files and create uniq list of all objects from all files
            TFile *f = TFile::Open(infilenames.begin()->Data(), "read");
            if(f==0||!f->IsOpen()) printf("Cannot Open file: %s \n", infilenames.begin()->Data());
            // get directory
            TIter iKey(f->GetListOfKeys());
            TKey* key=0;
            TDirectory* dir = 0;
            if ((key=(TKey*)iKey())) dir=(TDirectory*)key->ReadObjectAny(TDirectory::Class());
            // get dir name (if not root)
            TString dirname = "";
            if (dir == 0) dir = f;
            else dirname = (TString)dir->GetName() + "/";
            TIter nextkey(dir->GetListOfKeys());
            while (key = (TKey*)nextkey())
            {
                // filter classes
                TClass *cl = TClass::GetClass(key->GetClassName());
                if (!cl->InheritsFrom( "TH1"      )) continue; // profiles and histograms of all dimensions
                TString name = key->GetName();
                if (cl->InheritsFrom( "TH3"      )) {
                    // if is pt,y,M make projection to y,M
                    if (do2D && name.Contains("h_qtVyVQ")) all_obj_names.push_back(dirname+name+"_py");
                    continue; // Dont merge TH3 is too slow
                }
                if (doTest && !name.Contains("p_qtVy_A4")) continue; // for fast merger testing
                all_obj_names.push_back(dirname+name);
                // make projections on TH2 and remove outlier on 1D separatelly
                if (do2D && cl->InheritsFrom( "TH2"      )) {
                    all_obj_names.push_back(dirname+name+"_px");
                    all_obj_names.push_back(dirname+name+"_py");
                    //all_obj_names.push_back(dirname+name+"_pu");
                    //all_obj_names.push_back(dirname+name+"_pv");
                }
                // rebin pt as used ptz measurement
                if (doRebin) {
                    if ( name.EqualTo("pt") || name.EqualTo("h_qt") ) {
                        all_obj_names.push_back(dirname+name+"_rebin");
                    }
                }
                if (doTest && all_obj_names.size()==2) break; // for fast merger testing
            }
            f->Close();
        }

        VecInt loop_bins(TH1*h){
            VecInt bins;
            int dim=h->GetDimension();
            int nbinsx= (dim<1) ? LOBIN : h->GetNbinsX()+HIBIN;
            int nbinsy= (dim<2) ? LOBIN : h->GetNbinsY()+HIBIN;
            int nbinsz= (dim<3) ? LOBIN : h->GetNbinsZ()+HIBIN;
            //if(verbose>6) printf("    LOOP BINS LO %d HI %d nx %d ny %d nz %d \n", LOBIN, HIBIN, nbinsx, nbinsy, nbinsz);
            for (int izbin=LOBIN;izbin<=nbinsz;izbin++){
                for (int iybin=LOBIN;iybin<=nbinsy;iybin++){
                    for (int ixbin=LOBIN;ixbin<=nbinsx;ixbin++){
                        int ibin=h->GetBin(ixbin,iybin,izbin);
                        //if(verbose>6) printf("    bin %d x %d y %d z %d \n",  ibin, ixbin, iybin, izbin);
                        bins.push_back(ibin);
                    }
                }
            }
            //if(verbose>6) printf("    ALL %d \n", (int) bins.size());
            return bins;
        }

        bool is_empty(TObject *o, const char * objname, const char * objname_out){
            // Test if object is empty
            // FIXME: use exceptions
            if (o==0){
                printf(" skipping object because there is none\n with name:  %s\n objname_out: %s\n", objname, objname_out);
                return true;
            }
            return false;
        }

        void push_sorted(VecDbl &vec, const double val){
            // from:  http://stackoverflow.com/a/15048651
            VecDbl::iterator it = std::lower_bound( vec.begin(), vec.end(), val, std::greater<double>() ); // find proper position in descending order
            vec.insert( it, val ); // insert before iterator it
        }

        double median(VecDbl xi) {
            if (xi.size() == 0) {
                printf( "Error, passed empty vector to double median(vector <double> xi)\n");
                return 0;
            }
            if (xi.size() == 1) return xi[0];
            double med = 0;
            sort(xi.begin(), xi.end());
            if (verbose >6){
                printf ("    sorted ");
                for (auto x : xi) printf ("%f  ",x);
                printf ("\n");
            }
            if (xi.size() % 2) //odd
                med = *(xi.begin() + ((xi.size() + 1) / 2) - 1);
            else { //even 
                med += *(xi.begin() + ((xi.size() / 2) - 1));
                med += *(xi.begin() + ((xi.size() / 2)));
                med /=2;
            }
            return med;
        }

        double delta(VecDbl xi, double central, double ConfLevel) {
            double delta = 0;
            VecDbl deltaxi;
            VecDbl::iterator i = xi.begin();
            while (i != xi.end()) {
                deltaxi.push_back(fabs(*i - central));
                ++i;
            }
            sort(deltaxi.begin(), deltaxi.end());
            VecDbl::iterator di = deltaxi.begin();
            while (di != deltaxi.end()) {
                delta = *di;
                int index = di - deltaxi.begin() + 1;
                double prob = (double)index / (double)deltaxi.size();
                //      cout << index << "  " << *di << "  " << prob << endl;
                if (prob > ConfLevel)
                    break;
                ++di;
            }
            return delta;
        }

        double mean(VecDbl xi)
        {
            if (xi.size() == 0) {
                printf( "Error in pdferrors.cc, passed empty vector to double mean(vector <double> xi) \n");
                return 0;
            }
            double avg = 0;
            for (VecDbl::iterator it = xi.begin(); it != xi.end(); it++)
                avg += *it;
            avg /= xi.size();
            return avg;
        }

        double rms(VecDbl xi,double mean) {
            if (xi.size() == 0) {
                printf( "Error in pdferrors.cc, passed empty vector to double rms(vector <double> xi)\n");
                return 0;
            }
            double sum2 = 0;
            for (VecDbl::iterator it = xi.begin(); it != xi.end(); it++)
                sum2 += pow(*it,2);
            sum2 /= xi.size();
            return sqrt(fabs(sum2 - pow(mean,2)));
        }

        double pl(int nsigma) {
            return TMath::Erfc(nsigma/TMath::Sqrt(2));
        }


        double chi2prob(TH1* test, TH1 * ref, bool cumul=false) {
            if (ref == 0 || test == 0) return 0;
            double c2 = 0;
            int dim  = ref->GetDimension();
            VecInt bins = loop_bins(ref);
            for ( auto b : bins ){
                double d = test->GetBinContent(b) - ref->GetBinContent(b);
                double s = ref->GetBinError(b);
                double chi2 = 0;
                if (s != 0)
                    chi2 = d*d/(s*s);
                if (cumul){
                    // cumulative chi2 of all bins
                    c2 += d*d/(s*s);
                } else {
                    // maximum chi2 of all bins
                    c2 = std::max(c2,chi2);
                }
            }
            if (cumul){
                // return the cumulative chi2 probability of all the bins
                return TMath::Prob(c2, bins.size());
            }
            // return the probability of the maximum deviation of all bins,
            // corrected for the look elsewhere effect (corrected assuming bins are independent)
            return TMath::Prob(c2,1) * bins.size();
        }

        void open_all_infiles(){
            for (auto it_fn : infilenames){
                const char * fname = it_fn.Data();
                all_files.push_back(TFile::Open(fname,"READ"));
            }
        }

        void close_all_infiles(){
            for (auto it_f : all_files) if (it_f->IsOpen()) it_f->Close();
        }

        void delete_objs(VecTH1 &in_objs){
            while (!in_objs.empty()){
                delete in_objs.back();
                in_objs.pop_back();
            }
        }

        TH1* find_all_bins(VecTH1 in_objs){
            TH1*obj=0;
            if (in_objs.empty()) return obj;
            TH1* tmp=in_objs[0];
            TAxis *ax = 0;
            int dim = tmp->GetDimension();
            std::set<double> xset,yset,zset;
            // loop all objs, get x,y,z bins
            for( TH1* o : in_objs ){
                //loop all bins
                ax=o->GetXaxis(); for (int i=0; i<=ax->GetNbins(); i++){ xset.insert(ax->GetBinUpEdge(i)); }
                ax=o->GetYaxis(); for (int i=0; i<=ax->GetNbins(); i++){ yset.insert(ax->GetBinUpEdge(i)); }
                ax=o->GetZaxis(); for (int i=0; i<=ax->GetNbins(); i++){ zset.insert(ax->GetBinUpEdge(i)); }
            }
            // create new empty object
            TString title;
            title.Form("%s;%s;%s;%s'",
                    tmp->GetTitle(),
                    tmp->GetXaxis()->GetTitle(), 
                    tmp->GetYaxis()->GetTitle(), 
                    tmp->GetZaxis()->GetTitle()
                    );
            TString name=tmp->GetName(); name+="_allbins";
            // create array from set (by vector constructor)
            VecDbl xbins(xset.begin(), xset.end());
            VecDbl ybins(yset.begin(), yset.end());
            VecDbl zbins(zset.begin(), zset.end());
            if (dim==1){
                obj = (TH1*) new TH1D (name.Data(), title.Data(),
                        xbins.size()-1, &xbins[0]
                        );
            } else if (dim==2){
                obj = (TH1*) new TH2D (name.Data(), title.Data(),
                        xbins.size()-1, &xbins[0],
                        ybins.size()-1, &ybins[0]
                        );
            // FIXME: I am not sure if this is necessary
            //} else if (dim==3){
                //obj = (TH1*) new TH3D (name.Data(), title.Data(),
                        //xbins.size()-1, &xbins[0],
                        //ybins.size()-1, &ybins[0],
                        //zbins.size()-1, &zbins[0]
                        //);
            }
            return obj;
        }

        void sum_bins(TH1*out, VecTH1 v_objs){
            // loop all objs
            for (TH1* o : v_objs) {
                // loop all bins
                for ( int ibin : loop_bins(o)){
                    double val ,err , curr_val , curr_err;
                    // get bin error, get bin val
                    val = o->GetBinContent(ibin);
                    if (val==0) continue;
                    err = o->GetBinError(ibin)/val;
                    // find bin
                    int ixbin = 0;
                    int iybin = 0;
                    int izbin = 0;
                    o->GetBinXYZ(ibin,ixbin,iybin,izbin);
                    int jxbin = out->GetXaxis()->FindBin(o->GetXaxis()->GetBinCenter(ixbin));
                    int jybin = out->GetYaxis()->FindBin(o->GetYaxis()->GetBinCenter(iybin));
                    int jzbin = out->GetZaxis()->FindBin(o->GetZaxis()->GetBinCenter(izbin));
                    int jbin = out->GetBin(jxbin,jybin,jzbin);
                    // sum values
                    curr_val = out->GetBinContent(jbin);
                    if (curr_val!=0) {
                        curr_err = out->GetBinError(jbin)/curr_val;
                        val+=curr_val;
                        err = sqrt(err*err + curr_err*curr_err) * val;
                    }
                    // set new value
                    out->SetBinContent(jbin, val);
                    out->SetBinError  (jbin, err);
                }
            }
            return;
        }

        void read_Xsection(){
            // Get value of X section from integral histogram
            VecDbl xi;
            for (auto it_f : all_files){
                push_sorted(xi, ((TH1*) it_f->Get("qt_y_total"))->Integral());
            }
            Xsection = mean(xi); // median(xi);
            // NOTE: The `qt_y_total` represents x-section output table of
            // dyturbo in form of the 2D histogram. Sometimes it happend that
            // histograms from "ressummed" term integrated with MC method had
            // different value of integral that final xsection. Therefore we
            // use this histogram to rescale "differential" histograms have
            // proper integral.
        }

        void normalize(TH1* h){
            double fac = Xsection/h->Integral();
            h->Scale(fac);
        }

        double interpolate_NaN(TH1* o, int ibin){
            // TODO: NaN in uncertainty treatment
            int dim = o->GetDimension();
            bool isProf = isProfile(o);
            double val=0;
            double ent=0;
            // create graph from neigbour bins
            if (dim==1){
                TGraph* gr = new TGraph();
                gr->SetPoint(0, ibin-1, o->GetBinContent(ibin-1) );
                gr->SetPoint(1, ibin+1, o->GetBinContent(ibin+1) );
                val=gr->Eval(ibin);
                delete gr;
                if (isProf) {
                    TProfile *p = (TProfile* )o;
                    TGraph* gr = new TGraph();
                    gr->SetPoint(0, ibin-1, p->GetBinEntries(ibin-1) );
                    gr->SetPoint(1, ibin+1, p->GetBinEntries(ibin+1) );
                    ent=gr->Eval(ibin);
                    delete gr;
                }
            } else if(dim==2){
                int xbin,ybin,zbin,jbin;
                o->GetBinXYZ(ibin,xbin,ybin,zbin);
                TGraph2D* gr = new TGraph2D();
                jbin=o->GetBin(xbin-1,ybin-1); gr->SetPoint(0, xbin-1, ybin-1,  o->GetBinContent(jbin) );
                jbin=o->GetBin(xbin-1,ybin+1); gr->SetPoint(1, xbin-1, ybin+1,  o->GetBinContent(jbin) );
                jbin=o->GetBin(xbin+1,ybin-1); gr->SetPoint(2, xbin+1, ybin-1,  o->GetBinContent(jbin) );
                jbin=o->GetBin(xbin+1,ybin+1); gr->SetPoint(3, xbin+1, ybin+1,  o->GetBinContent(jbin) );
                val=gr->Interpolate(xbin,ybin);
                delete gr;
                if (isProf) {
                    TProfile2D *p = (TProfile2D* ) o;
                    TGraph2D* gr = new TGraph2D();
                    jbin=p->GetBin(xbin-1,ybin-1); gr->SetPoint(0, xbin-1, ybin-1,  p->GetBinEntries(jbin) );
                    jbin=p->GetBin(xbin-1,ybin+1); gr->SetPoint(1, xbin-1, ybin+1,  p->GetBinEntries(jbin) );
                    jbin=p->GetBin(xbin+1,ybin-1); gr->SetPoint(2, xbin+1, ybin-1,  p->GetBinEntries(jbin) );
                    jbin=p->GetBin(xbin+1,ybin+1); gr->SetPoint(3, xbin+1, ybin+1,  p->GetBinEntries(jbin) );
                    ent=gr->Interpolate(xbin,ybin);
                    delete gr;
                }
            }
            // interpolate
            // set
            if (isProf){
                if      (dim==1) set_profile_bin( (TProfile   *) o,ibin,val,0,ent,0);
                else if (dim==2) set_profile_bin( (TProfile2D *) o,ibin,val,0,ent,0);
            } else o->SetBinContent(ibin, val);
            return val;
        }

        void create_average_obj(TH1* &tmp_m, VecTH1 &in_objs, TString type){
            int dim  = tmp_m->GetDimension();
            bool doEntries = type.CompareTo("median_entries",TString::kIgnoreCase)==0;
            if (doEntries) type = "median";
            bool isProf = isProfile(tmp_m);
            TProfile*   prof   = (dim==1&&isProf&&!doEntries) ? (TProfile   *)tmp_m : 0;
            TProfile2D* prof2D = (dim==2&&isProf&&!doEntries) ? (TProfile2D *)tmp_m : 0;
            if (isProf && doEntries){ // instead of profile value, take profile entries (denominator)
                TString name=tmp_m->GetName();
                if (doEntries) name+="_entries";
                TH1* old=tmp_m;
                if(dim==1){
                    tmp_m = ((TProfile*)tmp_m) -> ProjectionX(name);
                    delete old;
                } else if (dim==2){
                    tmp_m = ((TProfile2D*)tmp_m) -> ProjectionXY(name);
                    delete old;
                }
                tmp_m->SetDirectory(0);
            }
            tmp_m->Reset();
            double sqrtN=sqrt(in_objs.size());
            for ( auto ibin : loop_bins(tmp_m) ){
                VecDbl vals; // value in case of histogram, ratio in case of profile
                VecDbl entrs; // denominator in case of profile
                if ( verbose>6 ) printf ( " looping bins for average %d \n" , ibin);
                for(auto ith : in_objs){
                    if (ith!=0){
                        double entries = 0;
                        double value = ith->GetBinContent(ibin);
                        if (isProf){
                            if(dim==1){
                                entries = ((TProfile*)   ith)->GetBinEntries (ibin);
                                value   = ((TProfile*)   ith)->At            (ibin);
                            } else if (dim==2){
                                entries = ((TProfile2D*) ith)->GetBinEntries (ibin);
                                value   = ((TProfile2D*) ith)->At            (ibin);
                            }
                            if (doEntries){
                                push_sorted(vals,entries);
                            } else {
			        push_sorted(vals,(entries == 0) ? 0 : value/entries); // make median from ratio
                                push_sorted(entrs,entries);
                            }
                        } else push_sorted(vals,value);
                    }
                }
                if ( verbose>6 ) printf ( "      vals size %lu\n" , vals.size() );
                double centr = 0;
                double sigma = 0;
                double wcentr = 0;
                double wsigma = 0;
                // decide what type of moment to calculate
                if (type.CompareTo("median",TString::kIgnoreCase)==0){
                    centr = median(vals);
                    sigma = delta(vals, centr, 0.68) / sqrtN;
                    // if zero sigma take the value of first input object
                    if (sigma==0) sigma = in_objs[0]->GetBinError(ibin);
                    if (isProf&&!doEntries){
                        wcentr = median (entrs);
                        wsigma = delta (entrs, wcentr, 0.68) / sqrtN;
                        if (verbose>5) printf("   calculated median profile nom %f +- %f denom %f +- %f prof %f \n", centr*wcentr, sigma*wcentr, wcentr, wsigma, centr);
                    }
                } else if (type.CompareTo("mean",TString::kIgnoreCase)==0){
                    centr = mean(vals);
                    sigma = rms(vals, centr)/sqrtN;
                    if ( verbose>6 ) printf ( "      centr %f sigma %f  \n" , centr, sigma );
                }
                if (isProf&&!doEntries&&wcentr!=0){
                    // if zero sigma take the value of first input object
                    if (wsigma==0){
                        if (dim==1) wsigma = ((TProfile   *) in_objs[0])->GetBinSumw2()->At(ibin);
                        if (dim==2) wsigma = ((TProfile2D *) in_objs[0])->GetBinSumw2()->At(ibin);
                        wsigma = TMath::Sqrt(wsigma);
                    }
                    if (wsigma!=0 && wcentr!=0){
                        if (prof   !=0) set_profile_bin(prof   ,ibin,centr*wcentr,sigma*wcentr,wcentr,wsigma);
                        if (prof2D !=0) set_profile_bin(prof2D ,ibin,centr*wcentr,sigma*wcentr,wcentr,wsigma);
                    }
                } else {
                    tmp_m->SetBinContent( ibin, centr  );
                    tmp_m->SetBinError  ( ibin, sigma );
                }
            } // bins
            if (prof   !=0)  prof   ->SetErrorOption("i");
            if (prof2D !=0)  prof2D ->SetErrorOption("i");
        }

        void remove_outliers_from_hist(TH1*hist, TH1*ref, VecTH1 in_objs){
            remove_outliers_from_profile(hist,ref,in_objs, false);
        }

        void remove_outliers_from_profile(TH1* out, TH1* med, VecTH1 in_objs,bool isProf = true){
            if (med==0) return;
            int dim = med->GetDimension();
	    cout << out->GetName() << "  " << med->GetName() << "  " << dim << endl;
	    cout << med->ClassName() << endl;
	    TProfile*   prof   = (isProf&&dim==1) ? (TProfile   *)out : 0;
            TProfile2D* prof2D = (isProf&&dim==2) ? (TProfile2D *)out : 0;
	    cout << "check2" << endl;
	    cout << isProfile(med) << "  " << med << endl;
            for ( auto b : loop_bins(med) ){
                VecDbl v_sumw;
                VecDbl v_sumwy;
                // loop over all objects and skip outliers
                for (auto ith : in_objs){
                    double d = med->GetBinContent(b);
                    double s = med->GetBinError(b) * TMath::Sqrt(in_objs.size());
                    double ent = 0;
                    double val = 0;
                    if (isProf && dim==1){
		        //be careful!!! centr->GetBinContent gives the numerator of TProfile::GetBinContent
		        TProfile * centr = (TProfile*) med;
			cout << centr->GetName() << endl;
			d = (centr->GetBinEntries(b) != 0) ? centr->At(b)/centr->GetBinEntries(b) : 0;
		        s = centr->GetBinError(b) * TMath::Sqrt(in_objs.size());
		        TProfile * test = (TProfile*) ith;
                        //val = test->At (b);
                        ent = test->GetBinEntries (b);
                        //d -= (ent != 0) ? val/ent : 0;
                        val = test->GetBinContent(b);
                        d -= val;
                    } else if (isProf && dim==2) {
		        //be careful!!! centr->GetBinContent gives the numerator of TProfile::GetBinContent
		        TProfile2D * centr = (TProfile2D*) med;
			cout << centr->GetName() << endl;
			cout << isProfile(med) << endl;
			d = (centr->GetBinEntries(b) != 0) ? centr->At(b)/centr->GetBinEntries(b) : 0;
		        s = centr->GetBinError(b) * TMath::Sqrt(in_objs.size());
                        TProfile2D * test = (TProfile2D *) ith;
                        //val = test->At (b);
                        ent = test->GetBinEntries (b);
                        //d -= (ent != 0) ? val/ent : 0;
                        val = test->GetBinContent(b);
                        d -= val;
			printf("    b %d in_obj %s d %f s %f bc %f be %f \n", b, ith->GetName(), centr->GetBinContent(b), s, test->GetBinContent(b), test->GetBinError(b));
                    } else {
                        val = ith->GetBinContent(b);
                        d -= val;
                    }
                    double chi2 = 0;
                    if (s != 0) chi2 = d*d/(s*s);
                    if (verbose>1) printf("    b %d in_obj %s d %f s %f chi2 %f \n", b, ith->GetName(), d, s, chi2);
                    if ( TMath::Prob(chi2,1) > pl(7) ){
                        if (verbose>1) printf("   taken up to pl %g chi2 prob %g \n", pl(7), TMath::Prob(chi2,1));
                        if (isProf) push_sorted( v_sumw  , ent );
                        push_sorted( v_sumwy , val );
                    }
                } // objects
                // set new values from average after outlier removal
                if (!v_sumwy.empty()){
                    //
                    if (isProf){
                        // calculate new average
                        double sumw = mean(v_sumw);
                        double sumwy = mean(v_sumwy);
                        double sigmaw = (v_sumw.size() > 0) ? rms(v_sumw, sumw )/TMath::Sqrt(v_sumw.size()) : 0;
                        double sigmawy = (v_sumwy.size() > 0) ? rms(v_sumwy, sumwy )/TMath::Sqrt(v_sumwy.size()) : 0;
                        // set bin content and error
                        if (prof   !=0) set_profile_bin(prof   ,b, sumwy*sumw,sigmawy*sumw, sumw, sigmaw);
                        if (prof2D !=0) set_profile_bin(prof2D ,b, sumwy*sumw,sigmawy*sumw, sumw, sigmaw);
                    } else {
                        // calculate new average
                        double centr = mean(v_sumwy);
                        double sigma = (v_sumwy.size() > 0) ? rms(v_sumwy, centr)/TMath::Sqrt(v_sumwy.size()) : 0;
                        // set bin content and error
                        out->SetBinContent( b, centr  );
                        out->SetBinError  ( b, sigma );
                    }
                }
            } // all bins
            if (prof   !=0)  prof   ->SetErrorOption("i");
            if (prof2D !=0)  prof2D ->SetErrorOption("i");
        }

        template<class T>
        void set_profile_bin(T* prof,int ibin, double centr, double sigma, double  wcentr, double  wsigma, int  N=1){
	  /*
	    double p_centr2 = TMath::Power(centr/wcentr,2);
            /// Propagate error to ratio
            double p_sigma2 = 0;
            p_sigma2+= TMath::Power( sigma  / centr  ,2);
	    p_sigma2+= TMath::Power( wsigma / wcentr ,2);
            p_sigma2*= p_centr2;
	  */

            double p_centr2 = (wcentr != 0) ? TMath::Power(centr/wcentr,2) : 0;
            double p_sigma2 = (wcentr != 0) ? TMath::Power(sigma ,2)/TMath::Power(wcentr ,2) : 0;
	    
            /// Error definition inside TProfile:
            ///
            /// ERR = spread/sqrt(Neff)
            /// spread = sqrt(abs(  sumwyy/sumw - p_centr2 ) )
            /// Neff = sumw2 / sumww
            ///
            /// ERR**2 = abs(sumwyy/sumw - p_centr2 ) * sumww / sumw2
            /// ERR**2 * sumw2 / sumww = abs(sumwyy/sumw - p_centr2)
            ///
            /// let ERR**2 = p_sigma2
            ///
            /// sumwyy = ( p_sigma2*sumw2/sumww +  p_centr2) * sumw
            /// or
            /// sumwyy = ( p_sigma2*sumw2/sumww -  p_centr2) * sumw
            /// let Neff == 1
            /// sumwyy = ( p_sigma2 + p_centr2) * sumw
            /// or
            /// sumwyy = ( p_sigma2 - p_centr2) * sumw
            double new_sumw   = wcentr;
            double new_sumww  = wcentr*wcentr;
            double new_sumwy  = centr;
            int sign= (wcentr<0 ? -1: 1 );
            double new_sumwyy = (p_sigma2+sign*p_centr2)*TMath::Abs(wcentr);
            /// cross check
            double Neff= new_sumww/(new_sumw*new_sumw);
            double spread = TMath::Abs(new_sumwyy/new_sumw - TMath::Power(new_sumwy/new_sumw,2));
            double ERR = spread/Neff;
            ERR = sqrt(ERR);
            if (verbose>5) printf("   filling median profile bin %d cent %f+-%f wcent %f+-%f prof: %f+-%f+-%f \n",
                    ibin,
                    centr, sigma,
                    wcentr, wsigma, 
                    TMath::Sqrt(p_centr2), TMath::Sqrt(p_sigma2),ERR );
            prof-> SetBinEntries        ( ibin       , new_sumw );
            prof-> GetBinSumw2()->SetAt ( new_sumww  , ibin   );
            prof-> SetAt                ( new_sumwy  , ibin   );
            prof-> GetSumw2()->SetAt    ( new_sumwyy , ibin   );
            return;
        }

        bool isProfile(TH1*h){
            return TString(h->ClassName()).Contains("Profile");
        }

        void print_range(TH1*h){
            if(!isProfile(h)) {
                h->Print("range");
                return;
            }
            // assuming its profile
            int dim = h->GetDimension();
            if (dim==1) print_range_profile((TProfile   *) h);
            if (dim==2) print_range_profile((TProfile2D *) h);
        }

        template<class T>
        void print_range_profile(T*p){
            //TProfile * p = (TProfile *) h;
            double stat[10]; p->GetStats(stat);
            printf( " TProfile dump: name=%s, mem=%p, tsumwy=%f, tsumwyy=%f, tsumw=%f, tsumww=%f, tsumwx=%f, tsumwxx=%f, erroropt=%s \n",
                    p->GetName(), p,
                    stat[4], stat[5], stat[0], stat[1], stat[2], stat[3],
                    p->GetErrorOption()
                    );
            if (
                    stat[0]==0 &&
                    stat[1]==0 &&
                    stat[2]==0 &&
                    stat[3]==0 &&
                    stat[4]==0 &&
                    stat[5]==0
                    ) {
                printf("seems empty..\n");
                return;
            }
            for (auto b: loop_bins(p)){
                double prof = p->GetBinContent(b);
                double err = p->GetBinError(b);
                double sumwy = p->At(b);
                double sumww = (p->GetBinSumw2()->GetSize()==0) ? 0 : p->GetBinSumw2()->At(b);
                double sumw = p->GetBinEntries(b);
                double sumwyy = p->GetSumw2()->At(b);
                printf( " bin[%d]=%f+-%f\tsumwy=%f,\tsumwyy=%f,\tsumw=%f,\tsumww=%f \n",
                        b, prof,err,sumwy,sumwyy,sumw,sumww
                        );
            }
            return;
        }


        TH1* clone_empty(TH1* h,TString newname){
            if (verbose>2) printf( " empty clone %s \n", newname.Data());
	    //TH1* tmp = (TH1*) h->Clone(newname.Data());
	    TH1* tmp;
            bool isProf = isProfile(h);
            int dim = h->GetDimension();
	    cout << "cloning " << newname << "  " << isProf << "  " << dim << endl;
	    //if (isProfile(h) && dim == 1)
	    //  tmp = (TH1*) ((TProfile*)h)->Clone(newname.Data());
	    //if (isProfile(h) && dim == 2)
	    //  tmp = (TH1*) ((TProfile2D*)h)->Clone(newname.Data());
	    //else
	      tmp = (TH1*) h->Clone(newname.Data());
            //tmp->SetDirectory(0);
            tmp->Reset();
            if (verbose>2) print_range(tmp);
            return tmp;
        }

        TH1* make_projection(TH1* orig, char proj){
            TH1* o = 0;
            bool isProf = isProfile(orig);
            int dim = orig->GetDimension();
            //
            if (proj == 'x'){
                if (isProf) { // assuming 2D profile
                    TProfile2D* pr2D=(TProfile2D*) orig;
                    o = pr2D->ProfileX("dummy"); // ,-1,0,"e");
                } else { // normal histogram (2D or 3D)
                    if (dim==2){
                        TH2 *h2 = (TH2 *) orig;
                        o = h2->ProjectionX("dummy");
                    } else if (dim==3) {
                        TH3 *h3 = (TH3 *) orig;
                        o = h3->Project3D("zx");
                    }
                }
            } else if (proj == 'y'){
                if (isProf) { // assuming 2D profile
                    TProfile2D* pr2D=(TProfile2D*) orig;
                    o = pr2D->ProfileY("dummy"); // ,-1,0,"e");
                } else { // normal histogram (2D or 3D)
                    if (dim==2){
                        TH2 *h2 = (TH2 *) orig;
                        o = h2->ProjectionY("dummy");
                    } else if (dim==3) {
                        TH3 *h3 = (TH3 *) orig;
                        o = h3->Project3D("zy");
                    }
                }
            } else if (proj == 'u'){ // only for testing
                // 2D -> Y axis, but with "e" -- compute errors
                TH2 *h2 = (TH2 *) orig;
                o = h2->ProjectionY("dummy",0,-1,"e");
            } else if (proj == 'v'){ // only for testing
                // 2D -> Y axis, but with "e" and without underflow and overflow
                TH2 *h2 = (TH2 *) orig;
                int nxbins = h2->GetNbinsX();
                o = h2->ProjectionY("dummy",1,nxbins,"e");
            }
            return o;
        }

        // Data memebers
        TString  outfilename;
        VecTStr  infilenames;
        VecTFile all_files;
        VecTStr  all_obj_names;
        VecTH1   output_objects;
        // bin loop
        int LOBIN;
        int HIBIN;
        double Xsection;

};


/**
 * Description of main program
 *
 */
int main(int argc, char * argv[]){

    // Declare the supported options.
    po::Options opts(argv[0], " outfile infilenames \n\n Program for merging and averaging histograms.");
    // hidden arguments
    opts.add_options("Hidden")
        ("outfile"     , "Name of output file. ", po::value<SString>() )
        ("infilenames" , "List of input files. ", po::value<VecSStr>() )
        ("t,test"      , "Test parser and die."                        )
    ;
    // Program options
    opts.add_options("")
        ("h,help"            , "Print this help and die."                                  )
        ("m,median"          , "Add median (by default only avarage)"                      )
        ("o,outlier"         , "Add median and outlier removal (by default only avarage)"  )
        ("x,x-section"       , "Normalize histograms to Xsection."                         )
        ("p,2D-proj"         , "Make 2d projections and outliers for 2D."                  )
        ("r,rebin"           , "Add pt histograms with Z pt LHC 7TeV rebin."               )
        ("c,cuba"            , "Add all histograms from cubature integration."             )
        ("n,NaN-interpolate" , "Find NaNs and interpolate"                                 )
        ("v,verbose"         , "Increase verbosity (more v the more chaty (max=vvvvvvv))." )
    ;
    // Parse
    try {
        opts.parse_positional( std::vector<SString>({"outfile", "infilenames"}) );
        opts.parse(argc,argv);
    }
    catch (cxxopts::OptionException &e){
        printf("Bad arguments: %s \n",e.what());
        PRINT_HELP(4);
    }
    // test and die
    if (opts.count("test")){
        //
        printf("Testing parser: oufilename\n");
        if (!opts.count("outfile")){
            printf("No outputfile name. \n");
        } else {
            printf("%s\n",opts["outfile"].as<SString>().c_str());
        }
        //
        printf("Testing parser: oufilename\n");
        if (!opts.count("infilenames")){
            printf("No outputfile name. \n");
        } else {
            for(const auto& s: opts["infilenames"].as<VecSStr>()) printf("%s\n",s.c_str());
        }
        //
        printf("Testing parser: counts\n");
        printf ( "help     : %d \n" , opts.count ( "help"            )  ) ;
        printf ( "median   : %d \n" , opts.count ( "median"          )  ) ;
        printf ( "outlier  : %d \n" , opts.count ( "outlier"         )  ) ;
        printf ( "x-section: %d \n" , opts.count ( "x-section"       )  ) ;
        printf ( "2D-proj  : %d \n" , opts.count ( "2D-proj"         )  ) ;
        printf ( "cuba     : %d \n" , opts.count ( "cuba"            )  ) ;
        printf ( "NaN      : %d \n" , opts.count ( "Nan-interpolate" )  ) ;
        printf ( "verbose  : %d \n" , opts.count ( "verbose"         )  ) ;
    }


    // Print Help and die (nicely).
    if (opts.count("help")) PRINT_HELP(0);
    // setup merger
    OutlierRemoval merger;
    //
    if (opts.count("outfile")) merger.SetOutputFile(opts["outfile"].as<SString>().c_str());
    else {
        printf("Please write outputfile name. \n\n");
        PRINT_HELP(1);
    }
    //
    if (opts.count("infilenames")) for(const auto& s: opts["infilenames"].as<VecSStr>()) merger.AddInputFile(s.c_str());
    else {
        printf("Please write at least one inputfile. \n\n");
        PRINT_HELP(2);
    }
    if (opts.count("median"          )) merger.doMedian=true;
    if (opts.count("outlier"         )) merger.doOutlierRemoval=true;
    if (opts.count("x-section"       )) merger.doXsecNormalization=true;
    if (opts.count("2D-proj"         )) merger.do2D=true;
    if (opts.count("rebin"           )) merger.doRebin=true;
    if (opts.count("cuba"            )) {merger.doCuba=true; merger.doRebin=false;}
    if (opts.count("NaN-interpolate" )) {merger.doNanIterp=true;}
    if (opts.count("verbose"         )) merger.verbose=10*opts.count("verbose");
    if (opts.count("test"            )) merger.doTest=true;
    // run merger
    merger.Init();
    merger.Merge();
    merger.Write();

    // Root kill hack: Force termination
    pid_t pid = getpid();
    TString command;
    //command.Form("gdb --batch --eval-command 'call exit(0)' --pid  %i", (int)pid);
    command.Form("kill %i", (int)pid);
    int ret = system(command.Data());
    return ret;
}


#endif // MERGER
