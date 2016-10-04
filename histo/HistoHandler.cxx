#ifndef HistoHandler_CXX
#define HistoHandler_CXX 
/**
 * @file HistoHandler.cxx
 * Description of this macro
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-08-26
 */

#include "Kinematics.h"
#include "HistoHandler.h"
#include "HistoObjects.h"


#include <sstream>
using std::ostringstream;

using namespace Kinematics;

// interface:
void disabled_hists_setpdf_(int * npdf){
    HistoHandler::SetVariation(*npdf);
}

void disabled_hists_fill_pdf_(double p3[4], double p4[4], double *weight, int *npdf){
    HistoHandler::SetVariation(*npdf);
    HistoHandler::FillEvent(p3,p4,*weight);
}

void disabled_hists_real_dipole_(double p3[4], double p4[4], double *weight, int * nd){
    HistoHandler::FillDipole(p3,p4,*weight);
}

void disabled_hists_real_dipole_pdf_(double p3[4], double p4[4], double *weight, int * nd, int* npdf){
    HistoHandler::SetVariation(*npdf);
    HistoHandler::FillDipole(p3,p4,*weight);
}

void disabled_hists_real_event_(){
    HistoHandler::FillRealEvent();
}

void disabled_hists_real_event_pdf_(int *npdf){
    /// @attention This is very dangerous !!!
    HistoHandler::SetVariation(*npdf);
    HistoHandler::FillRealEvent();
}

void disabled_hists_finalize_(){
    HistoHandler::Terminate();
}



namespace HistoHandler{

    // declaration of global data members
    HistoContainer histos;
    bool isFillMode=true;
    String result_filename="results";
    int parent_pid=0;

    void Init() {
        Book();
        SetVariation(0); // always start from 0
        parent_pid=getpid();
    }

    void DeleteNonIntegrableHists(){
        for (size_t ih = 0; ih < histos.size(); ++ih) {
            if (!histos[ih]->IsIntegrationSafe()){
                printf ("Warining: Histogram `%s` contains non-integrable observable and will be deleted \n", histos[ih]->GetName());
                histos[ih]->Delete();
                histos.erase(histos.begin()+ih);
                ih--;
            }
        }
    }

    void DeleteHists(){
        for (size_t ih = 0; ih < histos.size(); ++ih) {
            histos[ih]->Delete();
        }
        histos.clear();
    }

    void Reset(){
        for (size_t ih = 0; ih < histos.size(); ++ih) {
            histos[ih]->Reset();
        }
    }

    void FillEvent(double p3[4],double p4[4], double wgt){
        SetKinematics(p3,p4,wgt);
        for (size_t ih = 0; ih < histos.size(); ++ih) {
            histos[ih]->FillEvent();
        }
    }

    void FillResult(double int_val, double int_err){
        SetMiddlePoint();
        for (size_t ih = 0; ih < histos.size(); ++ih) {
            histos[ih]->AddToBin(int_val,int_err);
        }
    }

    void FillDipole(double p3[4],double p4[4], double wgt){
        SetKinematics(p3,p4,wgt);
        for (size_t ih = 0; ih < histos.size(); ++ih) {
            histos[ih]->FillDipole();
        }
    }

    void FillRealEvent(){
        for (size_t ih = 0; ih < histos.size(); ++ih) {
            histos[ih]->FillRealEvent();
        }
    }

#ifdef USEROOT
    String file_suffix = ".root";
#else
    String file_suffix = ".dat";
#endif

    size_t last_index=0;
    MapSuffix variation_suffixes;

    void SetVariation(int imember){
        // create missing index for variation
        if (variation_suffixes.count(imember)==0){
            // NOTE: if in future we find out that there is negative member, or
            // higher than maximum of LHAPDF memebers we can use this to define
            // scale variations
            KeySuffix addkey = {"",0};
            if (imember==0){ // add central
            } else { // PDF case
                /// @todo Check is PDF memeber i.e is positive and smaller than number of memebrs in pdfset
                ostringstream suffix;
                suffix << "_pdf" << imember;
                addkey.suffix=suffix.str();
                addkey.index=last_index;
            }
            // add it and increase index
            variation_suffixes.insert(MapSuffix::value_type(imember,addkey));
            last_index++;
        }
        // Pickup correct histogram
        for (size_t ih = 0; ih < histos.size(); ++ih) {
            histos[ih]->SetVariation(variation_suffixes.at(imember));
        }
    }

    // forward declaration : specialized in implementation header file
    /**
     * @brief Open file to save histograms.
     *
     * Implementation depends on ROOT or STL switch.
     *
     *  @todo always save to tmp file with timestamp. Only Terminate will save to final file.
     *
     * @param iworker Worker id from cuba.
     */
    void OpenOutputFile(int iworker);

    /**
     * @brief Close file with histograms.
     *
     * Implementation depends on ROOT or STL switch.
     *
     */
    void CloseOutputFile();

    /**
     * @brief Close file with histograms.
     *
     * Implementation depends on ROOT or STL switch.
     * Merge files from all cores.
     */
    void Merge(int iworker);

    void Save(int iworker){
        OpenOutputFile(iworker);
        for (size_t ih = 0; ih < histos.size(); ++ih) {
            histos[ih]->Save();
        }
        CloseOutputFile();
    }

    void Terminate(int iworker){
        Save(iworker);
        Merge(iworker);
        DeleteHists();
    }
}


#endif /* ifndef  */

