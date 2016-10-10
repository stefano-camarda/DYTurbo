#ifndef HistoBase_H
#define HistoBase_H
/**
 * @file HistoBase.h
 * @brief A brief description
 *
 * Detailed description of file
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-10-05
 */

#include "HistoHandler.h"
#include "Kinematics.h"
#include "src/settings.h"

namespace HistoHandler {
    // Base Histogram class: make possible to store in one container
    class HistoBase{
        public :
            virtual void FillEvent()=0;
            virtual void FillDipole()=0;
            virtual void FillRealEvent()=0;
            virtual void SetVariation(const KeySuffix)=0;
            virtual void Save()=0;
            virtual void Delete()=0;
            virtual void Reset()=0;
            virtual void AddToBin(double int_val,double int_err)=0;
            virtual double GetEntries() const =0; //{ return 666.;};
            virtual const char* GetName() const =0; // {return "Belzeebos";};
            virtual bool IsIntegrationSafe() const =0; //{ return 666.;};
        protected :
            String name;
            String title;
    };


    // Object class: covering common functionality
    template <class TH>
    class HistoObject : virtual public HistoBase {
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
            string titleX;
            string titleY;
            string titleZ;

            ~HistoObject(){
                Delete();
            }

            virtual void Init(String bin_name_X, String bin_name_Y="", String bin_name_Z=""){
                titleX=bin_name_X;
                titleY=bin_name_Y;
                titleZ=bin_name_Z;
                title=";";
                title+=bin_name_X;
                title+=";";
                title+=bin_name_Y;
                title+=";";
                title+=bin_name_Z;
                // get binning, name
                name="";
                bins.GetBins(bin_name_X,binsX);
                N_binsX=binsX.size()-1;
                p_binsX=&binsX.front();
                name+=bin_name_X;
                if (bin_name_Y.length()!=0){
                    bins.GetBins(bin_name_Y,binsY);
                    N_binsY=binsY.size()-1;
                    p_binsY=&binsY.front();
                    name+="_vs_";
                    name+=bin_name_Y;
                }
                if (bin_name_Z.length()!=0){
                    bins.GetBins(bin_name_Z,binsZ);
                    N_binsZ=binsZ.size()-1;
                    p_binsZ=&binsZ.front();
                    name+="_vs_";
                    name+=bin_name_Z;
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

            void SetName(string newname){
                // set new name
                this->name = newname;
                // and rename all histograms
                current->SetName(this->name.c_str());
                if(!variations.empty()) for (size_t i = 0; i < last_index; ++i) {
                    KeySuffix ivar = variation_suffixes[i];
                    if (ivar.index>=variations.size()) continue;
                    variations[ivar.index]->SetName((this->name+ivar.suffix).c_str());
                }
            }


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

            template <class TT> friend TT* New(HistoObject<TT> *h);
    };


    template <class TH> TH* New(HistoObject<TH> *h){
        return NULL;
    };
}



#endif /* ifndef HistoBase_H */
