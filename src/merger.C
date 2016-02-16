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

#include <vector>
#include <algorithm>

#include <TFile.h>
#include <TH1.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TList.h>
#include <TClass.h>
#include <TKey.h>
#include <TString.h>
#include <TMath.h>

typedef std::vector<TString> VecTStr;
typedef std::vector<TFile*> VecTFile;
typedef std::vector<TH1*> VecTH1;
typedef std::vector<double> VecDbl;
typedef std::vector<int> VecInt;
//using std::sort;

double bins[23] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 60, 70, 80, 100};

class OutlierRemoval{
    public :
        OutlierRemoval() :
            includeUnderOverFlow(true),
            doXsecNormalization(false),
            doRebin(true),
            do2D(false),
            verbose(1) {};
        ~OutlierRemoval(){};

        // Public methods
        void Merge(){
            find_all_objects();
	    //            open_all_infiles();
            if (doXsecNormalization) read_Xsection();
            // for object in allobjects
            for (auto p_objname : all_obj_names){
                const char* objname = p_objname.Data();
                size_t len = p_objname.Length();
                char proj = p_objname(len-1);
                if (verbose>1) printf("objname: %s\n",objname);
                // get object from all files
                VecTH1 in_objs;
                //                for (auto it_f : all_files){
                for (auto it_fn : infilenames){
                    const char * fname = it_fn.Data();
                    TFile * it_f =TFile::Open(fname,"READ");

                    /// @todo: test if they are all same binning
                    TString tmp(p_objname); tmp+="__"; tmp+=in_objs.size();
                    const char* objname_out = tmp.Data();
                    if (verbose>3) printf(" objname out: %s\n",objname_out);
                    TH1 * o = (TH1*) it_f->Get(objname);
                    if (o!=0) o = (TH1*)  o->Clone(objname_out);
                    if (o==0 && len>3 ){ // not found histogram
                        // test if it not 1D projection and if yes create it
                        TString basename = p_objname(0,len-3);
                        if (verbose>1) printf("basename: %s , proj %c , objname out: %s \n",basename.Data(), proj, objname_out);
                        TH2* h2d = (TH2*) it_f->Get(basename.Data());
                        if ( proj == 'x') {
                            if (is_empty(h2d,objname,objname_out)) continue;
                            o = h2d->ProjectionX("dummy"); // ,-1,0,"e");
                            o->SetName(objname_out);
                            if (doXsecNormalization) normalize(o);
                            delete h2d;
                        } else if (proj == 'y'){
                            if (is_empty(h2d,objname,objname_out)) continue;
                            o = h2d->ProjectionY("dummy"); // ,-1,0,"e");
                            o->SetName(objname_out);
                            if (doXsecNormalization) normalize(o);
                            delete h2d;
                        } else if (proj == 'u'){
                            if (is_empty(h2d,objname,objname_out)) continue;
                            o = h2d->ProjectionY("dummy",0,-1,"e");
                            o->SetName(objname_out);
                            if (doXsecNormalization) normalize(o);
                            delete h2d;
                        } else if (proj == 'v'){
                            if (is_empty(h2d,objname,objname_out)) continue;
                            o = make_projection(h2d);
                            o->SetName(objname_out);
                            if (doXsecNormalization) normalize(o);
                            delete h2d;
                        } else if (proj == 'n') { // its probably pt so rebin to ptz measurement
                            TString basename = p_objname(0,len-6); // remove "_rebin"
                            if (verbose>1) printf("basename: %s , proj %c \n",basename.Data(), proj);
                            TH1* h1 = (TH1*) it_f->Get(basename.Data());
                            if (is_empty(h1,objname,objname_out)) continue;
                            if (doXsecNormalization) normalize(h1);
                            o=h1->Rebin(22,objname_out,bins);
                            // divide by bin width
                            for (int ibin =1; ibin<=o->GetNbinsX(); ibin++ ){
                                o->SetBinContent(ibin,o->GetBinContent(ibin)/o->GetBinWidth(ibin));
                            }
                            delete h1;
                        }
                        if (is_empty(o,objname,objname_out)) continue;
                        for ( auto ibin : loop_bins(o)){
                            if (o->GetBinContent(ibin) != o->GetBinContent(ibin))
                                printf(" NAN bin: %d",ibin);
                        }
                    } else if (doXsecNormalization && !isProfile(o)) normalize(o);
                    if (verbose>2) o->Print();
                    in_objs.push_back(o);
		    it_f->Close();
                } //end loop all files
                // temporary objects
                TString name = p_objname;
                // average
                TH1* tmp_a = clone_empty(in_objs[0],(name+"_average").Data());
                if(verbose>4) {printf("  new obj "); tmp_a ->Print("range"); }
                int dim = tmp_a->GetDimension();
                // median
                TH1* tmp_m = (dim==1||do2D) ? clone_empty(in_objs[0],(name+"_median").Data()) : 0;
                // total profile
                TH1* tmp_p=0;
                if (isProfile(tmp_a)){
                    if (verbose>2) printf(" do profile \n");
                    tmp_p = clone_empty(in_objs[0],(name+"_total").Data());
                }
                // First Loop :  get average
                if (verbose>1) printf("    First Loop \n");
                if (dim==1||do2D) create_average_obj(tmp_m,in_objs,"median");
                if(verbose>4) {printf("  average "); tmp_m ->Print("range"); }
                if (tmp_p==0){ // is not profile
                    create_average_obj(tmp_a,in_objs,"mean"); // dont calculate average for profile
                    // save average and median
                    tmp_a->SetName(name);
                    output_objects.push_back(tmp_a);
                    if (dim==1||do2D) {
                        tmp_m->SetName((name+"_median").Data());
                        output_objects.push_back(tmp_m);
                        if (verbose>2) printf("saving histogram %s with integral %f\n", tmp_m->GetName(), tmp_m->Integral());
                    }

                    // Stop here for TH1,2,3
                    goto releaseobj;
                }
                // Profiles only -- total profile with a priori uncertainty
                for (auto ith : in_objs) tmp_p->Add(ith,1./in_objs.size());
                // save profile and median
                tmp_p->SetName(name);
                output_objects.push_back(tmp_p);
                if (verbose>2) printf("saving histogram %s with integral %f\n",tmp_p->GetName(), tmp_p->Integral());
                if (dim>1&&!do2D) goto releaseobj;
                tmp_m->SetName((name+"_median").Data());
                output_objects.push_back(tmp_m);
                // new median of entries
                tmp_m = clone_empty(in_objs[0], (name+"_entries_median").Data() );
                create_average_obj(tmp_m,in_objs,"median_entries");
                tmp_m->SetName((name)+"_entries_median");
                if (verbose>2) printf("saving histogram %s with integral %f\n",tmp_m->GetName(), tmp_m->Integral());
                output_objects.push_back(tmp_m);
                //
                // Second Loop -- discard outliers
                if (verbose>1) printf("    Second Loop \n");
                remove_outliers_from_profile(tmp_a, tmp_m, in_objs);
                if (tmp_a==0) goto releaseobj; // if not implemented 
                tmp_a->SetName((name+"_outlier").Data());
                output_objects.push_back(tmp_a);
                if (verbose>2) {printf("  outliers removed \n"); tmp_a->Print("range"); tmp_p->Print("range");}
                //
                // release input histograms
                releaseobj:
                    while (!in_objs.empty()){
                        delete in_objs.back();
                        in_objs.pop_back();
                    }
            }
            if (verbose>1) printf("End of merge, closing files \n");
	    //            close_all_infiles();
        };

        void Write(){
            if (verbose>1) printf("Write \n");
            TFile * fout = new TFile(outfilename.Data(), "RECREATE");
            for (auto ith : output_objects){
                if (verbose>4) printf("  writing pointer %p \n", ith);
                if (verbose>3) ith->Print();
                ith->Write();
            }
            if(verbose>0) printf("Merged %d files in %s\n", (int) infilenames.size(), outfilename.Data());
            fout->Close();
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
        }

        // Public data members
        bool includeUnderOverFlow;
        bool doXsecNormalization;
        bool do2D;
        bool doRebin;
        int verbose;

    private :
        // Private methods
        void find_all_objects(){
            TFile *f = TFile::Open(infilenames.begin()->Data(), "read");
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
                TObject* o = key->ReadObj();
                TString name = o->GetName();
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
            }
            f->Close();
        }

        VecInt loop_bins(TH1*h){
            VecInt bins;
            int dim=h->GetDimension();
            int nbinsx= (dim<1) ? LOBIN : h->GetNbinsX()+HIBIN;
            int nbinsy= (dim<2) ? LOBIN : h->GetNbinsY()+HIBIN;
            int nbinsz= (dim<3) ? LOBIN : h->GetNbinsZ()+HIBIN;
            if(verbose>5) printf("    LOOP BINS LO %d HI %d nx %d ny %d nz %d \n", LOBIN, HIBIN, nbinsx, nbinsy, nbinsz);
            for (int izbin=LOBIN;izbin<=nbinsz;izbin++){
                for (int iybin=LOBIN;iybin<=nbinsy;iybin++){
                    for (int ixbin=LOBIN;ixbin<=nbinsx;ixbin++){
                        int ibin=h->GetBin(ixbin,iybin,izbin);
                        if(verbose>5) printf("    bin %d x %d y %d z %d \n",  ibin, ixbin, iybin, izbin);
                        bins.push_back(ibin);
                    }
                }
            }
            if(verbose>5) printf("    ALL %d \n", (int) bins.size());
            return bins;
        }


        TH1D* make_projection(TH2*h2d){
            int nxbins = h2d->GetNbinsX();
            TH1D * o = h2d->ProjectionY("dummy",1,nxbins,"e");
            return o;
        }

        bool is_empty(TObject *o, const char * objname, const char * objname_out){
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
            double med = 0;
            sort(xi.begin(), xi.end());
            if (verbose > 10){
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

        void read_Xsection(){
            VecDbl xi;
            for (auto it_f : all_files){
                push_sorted(xi, ((TH1*) it_f->Get("qt_y_total"))->Integral());
            }
            Xsection = median(xi);
        }

        void normalize(TH1* h){
            double fac = Xsection/h->Integral();
            h->Scale(fac);
        }

        void create_average_obj(TH1* &tmp_m, VecTH1 &in_objs, TString type){
            int dim  = tmp_m->GetDimension();
            bool doEntries = type.CompareTo("median_entries",TString::kIgnoreCase)==0;
            if (doEntries) type = "median";
            bool isProf = isProfile(tmp_m);
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
                tmp_m->Reset();
            }
            double sqrtN=sqrt(in_objs.size());
            for ( auto ibin : loop_bins(tmp_m) ){
                VecDbl vals;
                VecDbl entrs;
                if ( verbose>5 ) printf ( " looping bins for average %d \n" , ibin);
                for(auto ith : in_objs){
                    if (ith!=0){
                        double entries = 0;
                        double value = ith->GetBinContent(ibin);
                        if (isProf){
                            if(dim==1){
                                entries = ((TProfile*)   ith)->GetBinEntries(ibin);
                            } else if (dim==2){
                                entries = ((TProfile2D*) ith)->GetBinEntries(ibin);
                            }
                            if (doEntries){
                                push_sorted(vals,entries);
                            } else {
                                push_sorted(vals,value);
                                push_sorted(entrs,entries);
                            }
                        } else push_sorted(vals,value);
                    }
                }
                double centr = 0;
                double sigma = 0;
                double wcentr = 0;
                double wsigma = 0;
                // decide what type of moment to calculate
                if (type.CompareTo("median",TString::kIgnoreCase)==0){
                    centr = median(vals);
                    sigma = delta(vals, centr, 0.68) / sqrtN;
                    if (isProf&&!doEntries){
                        wcentr = median (entrs);
                        wsigma = delta (entrs, wcentr, 0.68) / sqrtN;
                    }
                } else if (type.CompareTo("mean",TString::kIgnoreCase)==0){
                    centr = mean(vals);
                    sigma = rms(vals, centr)/sqrtN;
                }
                if (isProf&&!doEntries){
                    ((TProfile*)tmp_m)-> SetBinEntries ( ibin  , wcentr );
                    ((TProfile*)tmp_m)-> SetAt         ( centr , ibin   );
                } else {
                    tmp_m->SetBinContent( ibin, centr  );
                    tmp_m->SetBinError  ( ibin, sigma );
                }
            } // bins
        }


        void remove_outliers_from_profile(TH1*prof, TH1*med, VecTH1 in_objs){
            if (med==0) return;
            int dim = med->GetDimension();
            TProfile* out_obj = (TProfile*)prof;
            if (dim!=1) return; /// @todo: 2D
            for ( auto b : loop_bins(med) ){
                VecDbl v_sumw;
                VecDbl v_sumwy;
                // loop over all objects and skip outliers
                for (auto ith : in_objs){
                    TProfile * test = (TProfile*) ith;
                    double d = test->GetBinEntries(b) - med->GetBinContent(b);
                    double s = med->GetBinError(b) * TMath::Sqrt(in_objs.size());
                    double chi2 = 0;
                    if (s != 0) chi2 = d*d/(s*s);
                    if (verbose>1) printf("    b %d in_obj %s d %f s %f chi2 %f \n", b, ith->GetName(), d, s, chi2);
                    if ( TMath::Prob(chi2,1) > pl(7) ){
                        if (verbose>1) printf("   taken up to pl %g chi2 prob %g \n", pl(7), TMath::Prob(chi2,1));
                        push_sorted( v_sumw  ,test->GetBinEntries (b) );
                        push_sorted( v_sumwy ,test->At            (b) );
                    }
                } // objects
                // set new values from average after outlier removal
                if (!v_sumwy.empty()){
                    //
                    double sumw = mean(v_sumw);
                    double sumwy = mean(v_sumwy);
                    //double sigmaw = rms(v_sumw, sumw );
                    //double sigma = rms(v_sumwy, sumwy );
                    //double sumwyy = sumwy*sumwy + sigma*sigma;
                    //
                    out_obj-> SetBinEntries        ( b      , sumw );
                    out_obj-> SetAt                ( sumwy  , b    );
                    //out_obj-> GetBinSumw2()->SetAt ( sumw   , b    );
                    //out_obj-> GetSumw2()->SetAt    ( sumwyy , b    );
                }
            } // all bins
        }

        bool isProfile(TH1*h){
            return TString(h->ClassName()).Contains("Profile");
        }

        TH1* clone_empty(TH1* h,TString newname){
            if (verbose>2) printf( " empty clone %s \n", newname.Data());
            TH1* tmp = (TH1*) h->Clone(newname.Data());
            tmp->SetDirectory(0);
            tmp->Reset();
            if (verbose>2) tmp->Print();
            return tmp;
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






void help(const char * prog){
      printf ("usage: %s [-h] [-v] [-x]  <output> <input list>\n", prog);
      printf (" Please keep separated switches!!! Its on my todolist! \n");
      printf ("   -x    Normalize histograms to Xsection. \n");
      printf ("   -v    Increase verbosity (very chatty). \n");
      printf ("   -2d   Make 2d projections and outliers for 2D. \n");
      printf ("   -h    Print this help message. \n");
      printf ("\n For more info read README.md\n");
}

bool isOpt(const char* test, const char * argvi){
    return TString(argvi).CompareTo(test,TString::kIgnoreCase)==0;
}

/**
 * Description of main program
 *
 */
int main(int argc, const char * argv[]){
    int i=1;
    if (argc < 3)
    {
        if ( isOpt("-h",argv[i]) || isOpt("--help",argv[i]) ) {
            help(argv[0]);
            i++;
            return 0;
        }
        printf("Not enough arguments (at least 1 output and 1 inputs )\n");
        help(argv[0]);
        return 1;
    }
    OutlierRemoval merger;
    // parse argumets
    bool isInputSet=false;
    for (; i < argc; i++){
        // need help ?
        if ( isOpt("-h",argv[i]) || isOpt("--help",argv[i]) ) {
            help(argv[0]);
            return 0;
        }
        // parse x section normalization
        if (isOpt("-x",argv[i])){
            merger.doXsecNormalization=true;
            i++;
        }
        // parse verbosity level
        if (isOpt("-v",argv[i])){
            merger.verbose=100;
            i++;
        }
        // parse verbosity level
        if (isOpt("-2d",argv[i])){
            merger.do2D=true;
            i++;
        }
        //First non-optional argument is the output file
        if(!isInputSet){
            merger.SetOutputFile(argv[i]);
            isInputSet=true;
            i++;
        }
        //Next arguments are input files
        merger.AddInputFile(argv[i]);
    }// parse end
    merger.Init();
    merger.Merge();
    merger.Write();
    return 0;
}


#endif // MERGER
