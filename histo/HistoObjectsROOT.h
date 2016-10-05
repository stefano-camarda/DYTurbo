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


// ROOT includes
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile.h>
#include <TProfile2D.h>

#include "HistoHandler.h"
#include "HistoBase.h"

//#include <wordexp.h>
#include <ctime>

namespace HistoHandler {

    TFile * outf;

    void OpenHistoFile(String name){
        outf = TFile::Open(name.c_str(),"RECREATE");
        outf->cd();
    }

    void CloseOutputFile(){
        outf->Write();
        outf->Close();
        outf=0;
    }

    // wrapper for common functionality above ROOT histograms
    template<class TH> class HistoROOT : virtual public HistoObject<TH> {
        public :

            virtual void Save(){
                //printf("Saving %s : current %p variations.size %zu \n", name.c_str(), current, variations.size() );
                outf->cd();
                if (this->variations.size()==0) this->current->Write();
                else {
                    this->variations[0]->Write();
                    outf->mkdir((this->name+"_var").c_str());
                    outf->cd((this->name+"_var").c_str());
                    for (size_t ivar = 1; ivar < this->variations.size(); ++ivar) {
                        this->variations[ivar]->Write();
                    }
                    outf->cd();
                }
            }
    };

    template <> TH1D* New(HistoObject<TH1D> *parent){
        TH1D* h = new TH1D( parent ->name.c_str(), parent->title.c_str(), parent->N_binsX, parent->p_binsX);
        return h;
    }

}

#endif /* ifndef HistoObjectsROOT_H */
