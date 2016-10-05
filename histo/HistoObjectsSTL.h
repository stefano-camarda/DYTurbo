#ifndef HistoObjectsSTL_H
#define HistoObjectsSTL_H
/**
 * @file HistoObjectsSTL.hpp
 * Template implementation of STL interface.
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-08-28
 */



// TurboHist includes
#include "TurboHist.h"
#include "TurboHist_File.h"
#include "TurboHist_H1.h"
#include "TurboHist_H2.h"
#include "TurboHist_H3.h"
#include "TurboHist_P1.h"
#include "TurboHist_P2.h"

#include "HistoHandler.h"

namespace HistoHandler {

    TurboHist::File fout;

    void OpenHistoFile(string name){
        // recreate files
        fout.Open(name,"RECREATE");
    }

    void CloseOutputFile(){
        // close file
        fout.Close();
    }

    // wrapper for common functionality above STL histograms
    template<class TH> class HistoSTL : virtual public HistoObject<TH> {
        public :
            virtual void Save(){
                if (this->variations.size()==0) fout << *this->current;
                else {
                    for (size_t ivar = 0; ivar < this->variations.size(); ++ivar) {
                        fout << *this->variations[ivar];
                    }
                }
            }

    };

    template <> TurboHist::H1* New(HistoObject<TurboHist::H1> *parent){
        TurboHist::H1* h = new TurboHist::H1();
        h->SetBins(parent->binsX, parent->title);
        return h;
    }
}

#endif /* ifndef HistoObjectsSTL_H */
