#ifndef HistoObjectsROOT_H
#define HistoObjectsROOT_H
/**
 * @file HistoObjectsROOT.hpp
 * Template implementation of ROOT interface.
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-08-28
 */

#include "src/settings.h"

// ROOT includes
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile.h>
#include <TProfile2D.h>

#include "Kinematics.h"
#include "HistoHandler.h"

//#include <wordexp.h>
#include <ctime>

namespace HistoHandler {

    TFile * outf;

    void OpenOutputFile(int iworker){
        SStream outfname;
        outfname << "tmp_";
        outfname << result_filename;
        outfname << "_";
        int current_pid = getpid();
        if (current_pid!=parent_pid or iworker!=PARENT_PROC){
            // worker always create uniq file
            outfname << iworker << "_" << current_pid << "_" << time(NULL);
        } else {
            // main is always rewriten, because it stays in memory.
            outfname << "main";
        }
        outfname << file_suffix;
        outf = TFile::Open(outfname.str().c_str(),"RECREATE");
        outf->cd();
    }


    void CloseOutputFile(){
        outf->Write();
        outf->Close();
        outf=0;
    }

    void Merge(int iworker){
        // Merging of files in the end.
        if (getpid()==parent_pid || iworker==PARENT_PROC){ // only in main thread
            SStream tmp;
            // merge it with parent file
            tmp.str("");
            tmp << " hadd -f " << result_filename << file_suffix;
            tmp << " tmp_" << result_filename <<"_*" << file_suffix;
            // @todo if verbose print
            tmp << " > /dev/null ";
            // remove temporary files
            tmp << " && rm -f ";
            tmp << " tmp_" << result_filename <<"_*" << file_suffix;
            if (system(tmp.str().c_str())!=0) printf ("Something went wrong during `%s`. Better to throw exception",tmp.str().c_str());
        } // else do nothing, workers should not terminate anything
    }

    // wrapper for common functionality above ROOT histograms
    template<class TH> class HistoROOT : virtual public HistoBase {
        protected :
            // Common initialization
            VecDbl binsX;
            VecDbl binsY;
            VecDbl binsZ;
            double * p_binsX;
            double * p_binsY;
            double * p_binsZ;
            int N_binsX;
            int N_binsY;
            int N_binsZ;

            ~HistoROOT(){
                Delete();
            }

            virtual void Init(String bin_name_X, String bin_name_Y="", String bin_name_Z=""){
                name="";
                title=";";
                // get binning, name and title
                bins.GetBins(bin_name_X,binsX);
                N_binsX=binsX.size()-1;
                p_binsX=&binsX.front();
                name+=bin_name_X;
                title+=bin_name_X;
                if (bin_name_Y.length()!=0){
                    bins.GetBins(bin_name_Y,binsY);
                    N_binsY=binsY.size()-1;
                    p_binsY=&binsY.front();
                    name+="_vs_";
                    name+=bin_name_Y;
                    title+=";";
                    title+=bin_name_Y;
                }
                if (bin_name_Z.length()!=0){
                    bins.GetBins(bin_name_Z,binsZ);
                    N_binsZ=binsZ.size()-1;
                    p_binsZ=&binsZ.front();
                    name+="_vs_";
                    name+=bin_name_Z;
                    title+=";";
                    title+=bin_name_Z;
                }
            }

            // Histogram and PDF
            // histogram and container for all pdf variations
            TH* current=NULL;
            std::vector<TH*> variations;
            size_t current_variation=0;

        public :
            // changing PDF
            virtual void SetVariation(const KeySuffix ivar){
                // If same as current do nothing. This will turn off
                // completelly usage of variation_list vector if your
                // calculation doesnt need it.
                if (current_variation==ivar.index) return;
                // If empty add current hist to 0-th position
                // assuming we always starting from central
                size_t size = variations.size();
                if (size==0) variations.push_back(current);
                // Create new histograms if necessary
                size = variations.size();
                if (ivar.index == size){
                    // clone, rename and reset
                    TH * newmem = (TH*) current->Clone((name+ivar.suffix).c_str());
                    newmem->Reset();
                    variations.push_back(newmem);
                    size = variations.size();
                } else if(ivar.index>size) printf("ERROR: this should not happen ivar.index = %zu size= %zu\n", ivar.index, size);
                // change poiters to correct variation histograms and set last member
                current = variations[ivar.index];
                current_variation = ivar.index;
                return;
            }



            virtual void Save(){
                //printf("Saving %s : current %p variations.size %zu \n", name.c_str(), current, variations.size() );
                outf->cd();
                if (variations.size()==0) current->Write();
                else {
                    variations[0]->Write();
                    outf->mkdir((name+"_var").c_str());
                    outf->cd((name+"_var").c_str());
                    for (size_t ivar = 1; ivar < variations.size(); ++ivar) {
                        variations[ivar]->Write();
                    }
                    outf->cd();
                }
            }

            void Delete(){
                if (variations.size()==0) {
                    if (current!=NULL){
                        delete current;
                    }
                } else {
                    for (size_t ivar = 0; ivar < variations.size(); ++ivar) {
                        if (variations[ivar]!=NULL){
                            delete variations[ivar];
                            variations[ivar]=NULL;
                        }
                    }
                    variations.clear();
                }
                current = NULL;
            };


        protected :
            // Dipole structure:
            // Handling REAL event calculation by saving all points. Later it
            // is decided, which weights are summed to same bin ( for proper
            // uncertainty calculation)
            struct DipPt {
                int ibin;
                double valX;
                double valY;
                double valZ;
                double weight;
            } current_point;
            std::vector<DipPt> dipole_points;

            virtual int CurrentBin()=0;
            virtual void FillPoint(DipPt point)=0;

            virtual void AddPoint(){
                current_point.ibin = CurrentBin();
                current_point.weight = Kinematics::event_weight;
                dipole_points.push_back(current_point);
            }

        public :
            virtual void FillRealEvent(){
                DipPt point;
                // process all calculated contributions
                while (!dipole_points.empty() && isFillMode ){
                    // until you erase all points
                    point = dipole_points.back(); // take the last one
                    dipole_points.pop_back();
                    // remove all point with same bin index
                    // ( they are filled as one event with cumulated weight)
                    for (auto ipoint=dipole_points.begin(); ipoint!=dipole_points.end() ; ){
                        if (point.ibin == ipoint->ibin){
                            // found same bin: remove from list and sum weight
                            point.weight+=ipoint->weight;
                            ipoint = dipole_points.erase(ipoint);
                        } else ++ipoint; // bin index not same, skip it
                    }
                    // fill histograms in ibin with sum of all weights
                    if (point.weight!=0 ){ 
                        FillPoint(point);
                    }
                }
            };


            // Integrator Mode:
            bool isIntegrationSafe=false;
            inline bool IsIntegrationSafe() const {return isIntegrationSafe;}

            void AddToBin(double int_val, double int_err){
                if (!isIntegrationSafe) return;
                int ibin = CurrentBin();
                // if there somethin in bin add it properly
                double val = current->GetBinContent (ibin);
                double err = current->GetBinError   (ibin);
                val += int_val;
                err = sqrt (err*err + int_err*int_err);
                // set new values
                current->SetBinContent (ibin, val);
                current->SetBinError   (ibin, err);
            }

            inline double GetEntries() const {return current->GetEntries();}
            inline const char* GetName() const {return current->GetName();}

            inline void Reset() {
                current->Reset();
                for (size_t ivar = 0; ivar < variations.size(); ++ivar) variations[ivar]->Reset();
            }
    };


    // specialization
    template<class TX> class Histo1D : public HistoROOT<TH1D> {
        private:
            TX varX;

        public:
            Histo1D(const String &bin_name_X){
                // general init
                Init(bin_name_X);
                name = "s_"+name;
                title += ";#sigma[fb]";
                // create new
                current = new TH1D( name.c_str(), title.c_str(), N_binsX, p_binsX );
                isIntegrationSafe = varX.IsIntegrableObservable();
            }

            virtual void FillEvent(){
                current->Fill(varX(),Kinematics::event_weight);
            };

            virtual void FillDipole(){
                current_point.valX = varX();
                AddPoint();
            };

            virtual int CurrentBin(){
                return current->FindBin(varX());
            };

            void FillPoint(DipPt point){
                current->Fill(point.valX,point.weight);
            }
    };

    template<class TX, class TY> class Histo2D : public HistoROOT<TH2D> {
        public :
            Histo2D(const String &bin_name_X, const String &bin_name_Y){
                // general init
                Init(bin_name_X,bin_name_Y);
                name = "s_"+name;
                title += ";#sigma[fb]";
                // create new
                current = new TH2D( name.c_str(), title.c_str(), N_binsX, p_binsX, N_binsY, p_binsY );
                isIntegrationSafe = varX.IsIntegrableObservable() && varY.IsIntegrableObservable();
            }

            virtual void FillEvent(){
                current->Fill(varX(),varY(),Kinematics::event_weight);
            };

            virtual void FillDipole(){
                current_point.valX = varX();
                current_point.valY = varY();
                AddPoint();
            };

            virtual int CurrentBin(){
                return current->FindBin(varX(),varY());
            };

            void FillPoint(DipPt point){
                current->Fill(point.valX, point.valY ,point.weight);
            }


        private:
            TX varX;
            TY varY;
    };

    template<class TX, class TY, class TZ> class Histo3D : public HistoROOT<TH3D> {
        public :
            Histo3D(const String &bin_name_X, const String &bin_name_Y, const String &bin_name_Z){
                // general init
                Init(bin_name_X,bin_name_Y,bin_name_Z);
                name = "s_"+name;
                // create new
                current = new TH3D( name.c_str(), title.c_str(), N_binsX, p_binsX, N_binsY, p_binsY ,  N_binsZ, p_binsZ );
                isIntegrationSafe = varX.IsIntegrableObservable() && varY.IsIntegrableObservable() && varZ.IsIntegrableObservable();
            }

            virtual void FillEvent(){
                current->Fill(varX(),varY(), varZ(), Kinematics::event_weight);
            };

            virtual void FillDipole(){
                current_point.valX = varX();
                current_point.valY = varY();
                current_point.valZ = varZ();
                AddPoint();
            };

            virtual int CurrentBin(){
                return current->FindBin(varX(),varY(),varX());
            };

            void FillPoint(DipPt point){
                current->Fill(point.valX, point.valY, point.valZ, point.weight);
            }


        private:
            TX varX;
            TY varY;
            TZ varZ;
    };


    template<class TX, class TY> class HistoProfile : public HistoROOT<TProfile> {
        public :
            HistoProfile(const String &bin_name_X,const String &bin_name_Y){
                // general init
                Init(bin_name_X);
                // change name, add axis title
                title+=";";
                title+=bin_name_Y;
                name=bin_name_Y+"_"+name;
                // create new
                current = new TProfile( name.c_str(), title.c_str(), N_binsX, p_binsX);
                isIntegrationSafe = varX.IsIntegrableObservable() && varY.IsIntegrableObservable();
            }

            virtual void FillEvent(){
                current->Fill(varX(),varY(),Kinematics::event_weight);
            };

            virtual void FillDipole(){
                // Since the Ai moments are actualy weighted mean we need to do
                // weighted mean per each dipole point. We started by storing
                // the profiled value times weight.
                current_point.valX = varX();
                current_point.valY = varY()*Kinematics::event_weight;
                AddPoint();
            };

            virtual int CurrentBin(){
                return current->FindBin(varX());
            };

            void FillPoint(DipPt point){
                if (point.weight!=0) current->Fill(point.valX, point.valY/point.weight ,point.weight);
            }


        private:
            TX varX;
            TY varY;
    };

    template<class TX, class TY, class TZ> class HistoProfile2D : public HistoROOT<TProfile2D> {
        public :
            HistoProfile2D(const String &bin_name_X,const String &bin_name_Y, const String &bin_name_Z){
                // general init
                Init(bin_name_X, bin_name_Y);
                // change name, add axis title
                title+=";";
                title+=bin_name_Z;
                name=bin_name_Z+"_"+name;
                // create new
                current = new TProfile2D( name.c_str(), title.c_str(), N_binsX, p_binsX, N_binsY, p_binsY);
                isIntegrationSafe = varX.IsIntegrableObservable() && varY.IsIntegrableObservable() && varZ.IsIntegrableObservable();
            }

            virtual void FillEvent(){
                current->Fill(varX(),varY(),varZ(),Kinematics::event_weight);
            };

            virtual void FillDipole(){
                current_point.valX = varX();
                current_point.valY = varY();
                current_point.valZ = varZ()*Kinematics::event_weight;
                AddPoint();
            };

            virtual int CurrentBin(){
                return current->FindBin(varX(),varY());
            };

            void FillPoint(DipPt point){
                if (point.weight!=0) current->Fill(point.valX, point.valY, point.valZ/point.weight ,point.weight);
            }


        private:
            TX varX;
            TY varY;
            TZ varZ;
    };

}

#endif /* ifndef HistoObjectsROOT_H */
