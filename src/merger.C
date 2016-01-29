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
using std::sort;

double bins[23] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 60, 70, 80, 100};

class OutlierRemoval{
    public :
        OutlierRemoval() :
            includeUnderOverFlow(false),
            doXsecNormalization(false),
            doTh2dProjections(true),
            verbose(1) {};
        ~OutlierRemoval(){};

        // Public methods
        void Merge(){
            find_all_objects();
            open_all_infiles();
            if (doXsecNormalization) read_Xsection();
            // for object in allobjects
            for (auto p_objname : all_obj_names){
                const char* objname = p_objname.Data();
                size_t len = p_objname.Length();
                char proj = p_objname(len-1);
                if (verbose>1) printf("objname: %s\n",objname);
                // get object from all files
                VecTH1 in_objs;
                for (auto it_f : all_files){
                    /// @todo: test if they are all same binning
                    TH1 * o = (TH1*) it_f->Get(objname);
                    if (o==0 && len>3 ){ // not found histogram
                        // test if it not 1D projection and if yes create it
                        TString basename = p_objname(0,len-3);
                        if (verbose>1) printf("basename: %s , proj %c \n",basename.Data(), proj);
                        TH2* h2d = (TH2*) it_f->Get(basename.Data());
                        if ( proj == 'x') {
                            o = h2d->ProjectionX(objname); // ,-1,0,"e");
                        } else if (proj == 'y'){
                            o = h2d->ProjectionY(objname); // ,-1,0,"e");
                        }
                        for (int ibin=0; ibin < o->GetNbinsX()+1;ibin++){
                            if (o->GetBinContent(ibin) != o->GetBinContent(ibin))
                                printf(" NAN bin: %d",ibin);
                        }
                    }
                    if (verbose>2) o->Print();
                    if ( p_objname.EqualTo("pt") || p_objname.EqualTo("h_qt") ) {
                        o->Rebin(22,"zpt",bins);
                    }
                    in_objs.push_back(o);
                }
                // temporary objects
                TString name = p_objname;
                TH1* tmp_m = (TH1*) in_objs[0]->Clone((name+"median").Data());
                TH1* tmp_a = (TH1*) in_objs[0]->Clone((name+"average").Data());
                tmp_m->SetDirectory(0);
                tmp_m->Reset();
                tmp_a->SetDirectory(0);
                tmp_a->Reset();
                if (verbose>2) tmp_m->Print();
                if (verbose>2) printf(" before isprof\n");
                bool isProfile = TString(tmp_m->ClassName()).Contains("Profile");
                if (verbose>2) printf(" after isprof\n");
                if (verbose>2) printf(" before dim\n");
                int dim = tmp_m->GetDimension();
                if (verbose>2) printf(" after dim\n");
                // prepare total profile
                TH1* tmp_p=0;
                if (isProfile){
                    if (verbose>2) printf(" do profile \n");
                    tmp_p=tmp_m;
                    tmp_p->SetName("tot");
                    if(dim==1){
                        tmp_m=((TProfile *)tmp_p)->ProjectionX("average");
                    } else if (dim==2){
                        tmp_m=((TProfile2D *)tmp_p)->ProjectionXY("average");
                    }
                    tmp_m->SetDirectory(0);
                    tmp_p->Reset();
                }
                // First Loop :  get average
                if (verbose>1) printf("    First Loop \n");
                create_average_obj(tmp_m,in_objs,"median");
                create_average_obj(tmp_a,in_objs,"mean");
                // Stop here for TH1,2,3
                if (!isProfile){
                    tmp_m->SetName(name);
                    tmp_a->SetName((name+"_average").Data());
                    if (verbose>2) printf("writing histogram with integral %f\n",tmp_m->Integral());
                    if (doXsecNormalization) normalize(tmp_m);
                    output_objects.push_back(tmp_m);
                    output_objects.push_back(tmp_a);
                    continue;
                }
                // Second Loop -- discard outliers
                if (verbose>1) printf("    Second Loop \n");
                for(auto ith : in_objs ){
                    double p = chi2prob( ith, tmp_m);
                    if (p < pl(7)){
                        ith=0;
                    }
                }
                // Third Loop -- calculate average without outliers
                if (verbose>1) printf("    Third Loop \n");
                // Add non-outlier profiles
                for (auto ith : in_objs) tmp_p->Add(ith);
                tmp_p->SetName(name);
                if (verbose>2) printf("writing profile with integral %f\n",tmp_m->Integral());
                output_objects.push_back(tmp_p);
                // Save average object
                create_average_obj(tmp_m,in_objs,"mean");
                name+="_average";
                tmp_m->SetName(name);
                if (verbose>2) printf("writing histogram with integral %f\n",tmp_m->Integral());
                output_objects.push_back(tmp_m);
            }
            if (verbose>1) printf("Write \n");
            close_all_infiles();
        };

        void Write(){
            TFile * fout = new TFile(outfilename.Data(), "RECREATE");
            for (auto ith : output_objects){
                //ith->Print();
                ith->Write();
            }
            if(verbose>0) printf("Merged %d files in %s\n", infilenames.size(), outfilename.Data());
            fout->Close();
        };

        void SetOutputFile(TString _outf){
            outfilename=_outf;
            if(verbose>0)printf ("merger Target file: %s\n", _outf.Data()); // hadd like comments
        }

        void AddInputFile(TString _inf){
            infilenames.push_back(_inf); 
            if(verbose>0) printf ("merger Source file %d: %s\n", infilenames.size(), _inf.Data()); // hadd like comments
        }

        void Init(){
            LOBIN= includeUnderOverFlow ? 0 : 1;
            HIBIN= includeUnderOverFlow ? 1 : 2;
        }

        // Public data members
        bool includeUnderOverFlow;
        bool doXsecNormalization;
        bool doTh2dProjections;
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
                all_obj_names.push_back(dirname+o->GetName());
                // make projections on TH2 and remove outlier on 1D separatelly
                if (doTh2dProjections && cl->InheritsFrom( "TH2"      )) {
                    all_obj_names.push_back(dirname+o->GetName()+"_px");
                    all_obj_names.push_back(dirname+o->GetName()+"_py");
                }
            }
            f->Close();
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
            switch (nsigma) {
                case 0:
                    return 1.;
                    break;
                case 1:
                    return 1. - 0.682689492137086;
                    break;
                case 2:
                    return 1. - 0.954499736103642;
                    break;
                case 3:
                    return 1. - 0.997300203936740;
                    break;
                case 4:
                    return 1. - 0.999936657516334;
                    break;
                case 5:
                    return 1. - 0.999999426696856;
                    break;
                case 6:
                    return 1. - 0.999999998026825;
                    break;
                case 7:
                    return 1. - 0.999999999997440;
                    break;
                default:
                    printf("Confidence Level interval available only for sigma = 1-7, requested: %d sigma \n", nsigma );
                    return 1;
            }
        }


        double chi2prob(TH1* test, TH1 * ref, bool cumul=false) {
            if (ref == 0 || test == 0) return 0;
            double c2 = 0;
            int dim  = ref->GetDimension();
            int nbinsx= (dim<1) ? 1 : ref->GetNbinsX()+HIBIN;
            int nbinsy= (dim<2) ? 1 : ref->GetNbinsY()+HIBIN;
            int nbinsz= (dim<3) ? 1 : ref->GetNbinsZ()+HIBIN;
            // loop over all bins -- code inspired by TH1::Add()
            for (int izbin=LOBIN;izbin<nbinsz;izbin++){
                for (int iybin=LOBIN;iybin<nbinsy;iybin++){
                    for (int ixbin=LOBIN;ixbin<nbinsx;ixbin++){
                        int b = ixbin+nbinsx*(iybin +nbinsy*izbin);
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
                }
            }
            if (cumul){
                // return the cumulative chi2 probability of all the bins
                return TMath::Prob(c2, ref->GetNbinsX());
            }
            // return the probability of the maximum deviation of all bins,
            // corrected for the look elsewhere effect (corrected assuming bins are independent)
            return TMath::Prob(c2,1) * ref->GetNbinsX();
        }

        void open_all_infiles(){
            for (auto it_fn : infilenames){
                const char * fname = it_fn.Data();
                all_files.push_back(TFile::Open(fname,"READ"));
            }
        }

        void close_all_infiles(){
            for (auto it_f : all_files) it_f->Close();
        }

        void read_Xsection(){
            VecDbl xi;
            for (auto it_f : all_files){
                xi.push_back(((TH1*) it_f->Get("qt_y_total"))->Integral());
            }
            Xsection = median(xi);
        }

        void normalize(TH1* h){
            double fac = Xsection/h->Integral();
            h->Scale(fac);
        }

        void create_average_obj(TH1* tmp_m, VecTH1 &in_objs, TString type){
            int dim  = tmp_m->GetDimension();
            double sqrtN=sqrt(in_objs.size());
            int nbinsx= (dim<1) ? 1+LOBIN : tmp_m->GetNbinsX()+HIBIN;
            int nbinsy= (dim<2) ? 1+LOBIN : tmp_m->GetNbinsY()+HIBIN;
            int nbinsz= (dim<3) ? 1+LOBIN : tmp_m->GetNbinsZ()+HIBIN;
            // loop over all bins -- code inspired by TH1::Add()
            for (int izbin=LOBIN;izbin<nbinsz;izbin++){
                for (int iybin=LOBIN;iybin<nbinsy;iybin++){
                    for (int ixbin=LOBIN;ixbin<nbinsx;ixbin++){
                        int ibin = tmp_m->GetBin(ixbin,iybin,izbin); //ixbin+nbinsx*(iybin +nbinsy*izbin);
                        VecDbl vals;
                        for(auto ith : in_objs){
                            if (ith!=0){
                                //printf("   ibin %d content %f \n", ibin, ith->GetBinContent(ixbin,iybin,izbin));
                                //ith->Print("range");
                                //return;
                                vals.push_back(ith->GetBinContent(ibin));
                            } 
                        }
                        double centr = 0; // can change this in future
                        double sigma = 0;
                        if (type.CompareTo("median",TString::kIgnoreCase)==0){
                            centr = median(vals); 
                            sigma = delta(vals, centr, 0.68) / sqrtN;
                        } else if (type.CompareTo("mean",TString::kIgnoreCase)==0){
                            centr = mean(vals); 
                            sigma = rms(vals, centr);
                        }
                        tmp_m->SetBinContent( ibin, centr  );
                        tmp_m->SetBinError  ( ibin, sigma );
                    }
                }
            }
        }

        // Data memebers
        TString  outfilename;
        VecTStr  infilenames;
        VecTFile all_files;
        VecTStr  all_obj_names;
        VecTH1   output_objects;
        int LOBIN;
        int HIBIN;
        double Xsection;

};






void help(const char * prog){
      printf ("usage: %s [-X] [-v]  <output> <input list>\n");
      printf ("   -X    Normalize histograms to Xsection. \n");
      printf ("   -v    Increase verbosity. \n");
}

/**
 * Description of main program
 *
 */
int main(int argc, const char * argv[]){

    if (argc < 4)
    {
        printf("Not enough arguments (at least 1 output and 2 inputs )\n");
        help(argv[0]);
        return 1;
    }
    OutlierRemoval merger;
    int i= 1;
    // parse x section normalization
    if (TString(argv[i]).CompareTo("-X",TString::kIgnoreCase)==0){
        merger.doXsecNormalization=true;
        i++;
    }
    // parse verbosity level
    if (TString(argv[i]).CompareTo("-v",TString::kIgnoreCase)==0){
        merger.verbose=8;
        i++;
    }
    //First argument is the output file
    merger.SetOutputFile(argv[i]); i++;
    //Next arguments are input files
    for (; i < argc; i++){
        merger.AddInputFile(argv[i]);
    }
    merger.Init();
    merger.Merge();
    merger.Write();

    return 0;
}


#endif // MERGER
