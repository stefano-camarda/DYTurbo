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

#include "Kinematics.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile.h>
#include <TProfile2D.h>

#include <sstream>
using std::ostringstream;

#include <vector>
using std::vector;

#include <string>
using std::string;
#include <wordexp.h>

namespace HistoHandler {

    TFile * outf;

    void OpenOutputFile(int iworker){
        ostringstream outfname;
        outfname << result_filename;
        //if (iworker!=PARENT_PROC) outfname << "_" << iworker; // better to use getpid
        int current_pid = getpid();
        if (current_pid!=parent_pid) outfname << "_" << iworker << "_" << current_pid;
        outfname << ".root";
        outf = TFile::Open(outfname.str().c_str(),"RECREATE");
        outf->cd();
    }


    void CloseOutputFile(){
        outf->Write();
        outf->Close();
        outf=0;
        // Since there could be some parallel files from jobs which has been
        // finished try to merge them and save result hadd file.
        if (getpid()==parent_pid){ // only in main thread
            ostringstream tmp; tmp.str(""); // clear first
            tmp << "ls -la " << result_filename <<"_*.root > /dev/null 2>&1";
            if (system(tmp.str().c_str())==0) { // if there is ${result_filename}_*.root file:
                tmp.str(""); // clear tmp
                // this will include current workers as well as previous hadd
                tmp << "hadd addtmp.root " << result_filename <<"_*.root > /dev/null && ";
                // remove workers and previous hadd
                tmp << "rm " <<  result_filename <<"_*.root && ";
                // save temporary file file for future runs
                tmp << "mv addtmp.root " << result_filename <<"_hadd.root ";
                if (system(tmp.str().c_str())!=0) printf ("Something went wrong during `%s`. Better to throw exception",tmp.str().c_str());
            } // no results_*.root file: parent file was already rewritten
        } // else do nothing.. merge only after all workers are finished (beacause of async, you have to wait)
    }

    void Merge(int iworker){
        // Merging of files in the end.
        if (getpid()==parent_pid){ // only in main thread
            ostringstream tmp;
            tmp.str(""); tmp << "ls " << result_filename <<"_hadd.root > /dev/null 2>&1";
            if (system(tmp.str().c_str())==0) { // if there is hadd file.
                // merge it with parent file
                tmp.str("");
                tmp << "hadd addtmp.root " << result_filename <<".root " << result_filename <<"_hadd.root > /dev/null && ";
                // remove already merged files
                tmp << "rm " << result_filename <<".root " << result_filename <<"_hadd.root && ";
                // rename final merge
                tmp << "mv addtmp.root " << result_filename <<".root ";
                if (system(tmp.str().c_str())!=0) printf ("Something went wrong during `%s`. Better to throw exception",tmp.str().c_str());
            } // else do nothing, it was single thread run
        } // else do nothing, workers should not terminate anything

        // THis is wrong
        // async exitfun: could hadd wrong files
        // if (iworker==PARENT_PROC){
        //     // if there is tmp.root file:
        //         // hadd add.root result.root tmp.root
        //         // mv add.root result.root
        //         // rm tmp.root
        //     // not tmp.root file: do nothing
        // } else {  // is worker
        //     // if there is tmp.root file:
        //         // hadd add.root result_i.root tmp.root
        //         // mv add.root tmp.root
        //         // rm result_i.root
        //     // not temp file 
        //         // create one
        // }
    }

    // wrapper for common functionality above ROOT histograms
    template<class TH> class HistoROOT : virtual public HistoBase {
        protected :
            // Common initialization
            vector<double> binsX;
            vector<double> binsY;
            vector<double> binsZ;
            double * p_binsX;
            double * p_binsY;
            double * p_binsZ;
            int N_binsX;
            int N_binsY;
            int N_binsZ;

            virtual void Init(string bin_name_X, string bin_name_Y="", string bin_name_Z=""){
                // TODO: Implementing to DYTURBO bins -> get binning
                for (int i = 0; i < 101; ++i) {
                    binsX.push_back(i);
                    binsY.push_back(i);
                    binsZ.push_back(i);
                }
                name="";
                title=";";
                // get binning, name and title
                //bins.Get(bin_name_X,binsX);
                N_binsX=binsX.size()-1;
                p_binsX=&binsX.front();
                name+=bin_name_X;
                title+=bin_name_X;
                if (bin_name_Y.length()!=0){
                    //bins.Get(bin_name_Y,binsY);
                    N_binsY=binsY.size()-1;
                    p_binsY=&binsY.front();
                    name+="_vs_";
                    name+=bin_name_Y;
                    title+=";";
                    title+=bin_name_Y;
                }
                if (bin_name_Z.length()!=0){
                    //bins.Get(bin_name_Z,binsZ);
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
            TH* current;
            vector<TH*> variations;
            size_t current_variation;

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
                } else if(ivar.index>size) printf("ERROR: this should not happen\n");
                // change poiters to correct variation histograms and set last member
                current = variations[ivar.index];
                current_variation = ivar.index;
                return;
            }


            // TODO: Destructor:
            //   * Proper clear memory for current or whole variations
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
                if (current!=0) delete current;
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
            vector<DipPt> dipole_points;

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
                        // TODO: point fiducial
                        FillPoint(point);
                    }
                }
            };


            // Integrator Mode:
            void AddToBin(double int_val, double int_err){
                // TODO: if (!isIntegratorSafe) return; // dont fill if all variables are integrator safe
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
            inline void Clear() {return current->Reset();}


    };


    // specialization
    template<class TX> class Histo1D : public HistoROOT<TH1D> {
        public:
            Histo1D(const string &bin_name_X){
                // general init
                Init(bin_name_X);
                name = "s_"+name;
                title += ";#sigma[fb]";
                // create new
                current = new TH1D( name.c_str(), title.c_str(), N_binsX, p_binsX );
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


        private:
            TX varX;
    };

    template<class TX, class TY> class Histo2D : public HistoROOT<TH2D> {
        public :
            Histo2D(const string &bin_name_X, const string &bin_name_Y){
                // general init
                Init(bin_name_X,bin_name_Y);
                name = "s_"+name;
                title += ";#sigma[fb]";
                // create new
                current = new TH2D( name.c_str(), title.c_str(), N_binsX, p_binsX, N_binsY, p_binsY );
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
            Histo3D(const string &bin_name_X, const string &bin_name_Y, const string &bin_name_Z){
                // general init
                Init(bin_name_X,bin_name_Y,bin_name_Z);
                name = "s_"+name;
                // create new
                current = new TH3D( name.c_str(), title.c_str(), N_binsX, p_binsX, N_binsY, p_binsY ,  N_binsZ, p_binsZ );
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
            HistoProfile(const string &bin_name_X,const string &bin_name_Y){
                // general init
                Init(bin_name_X);
                // change name, add axis title
                title+=";";
                title+=bin_name_Y;
                name=bin_name_Y+"_"+name;
                // create new
                current = new TProfile( name.c_str(), title.c_str(), N_binsX, p_binsX);
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
            HistoProfile2D(const string &bin_name_X,const string &bin_name_Y, const string &bin_name_Z){
                // general init
                Init(bin_name_X, bin_name_Y);
                // change name, add axis title
                title+=";";
                title+=bin_name_Z;
                name=bin_name_Z+"_"+name;
                // create new
                current = new TProfile2D( name.c_str(), title.c_str(), N_binsX, p_binsX, N_binsY, p_binsY);
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
            TY varZ;
    };

}

#endif /* ifndef HistoObjectsROOT_H */
