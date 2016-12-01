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

#include "handy_typdefs.h"

#include "Kinematics.h"
#include "HistoHandler.h"

#include "HistoObjects.h"


using namespace Kinematics;

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
    String hadd_program = "hadd -f";
#else
    String file_suffix = ".dat";
    String hadd_program = "dyturbo-hadd";
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
                SStream suffix;
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

    /**
     * @brief Open file where histograms will be written
     *
     * Implementation depends on ROOT or STL switch.
     *
     * @param name Name of output file.
     */
    void OpenHistoFile(String name);

    /**
     * @brief Open file to save histograms.
     *
     *  @note Always save to tmp file with timestamp. Only Terminate will save to final file.
     *
     * @param iworker Worker id from cuba.
     */
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
        OpenHistoFile(outfname.str());
    }

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
    void Merge(int iworker){
        // Merging of files in the end.
        if (getpid()==parent_pid || iworker==PARENT_PROC){ // only in main thread
            SStream tmp;
            // merge it with parent file
            tmp.str("");
            tmp << " " << hadd_program << " " << result_filename << file_suffix;
            tmp << " tmp_" << result_filename <<"_*" << file_suffix;
            // @todo if verbose print
            tmp << " > /dev/null ";
            // remove temporary files
            tmp << " && rm -f ";
            tmp << " tmp_" << result_filename <<"_*" << file_suffix;
            if (system(tmp.str().c_str())!=0) printf ("Something went wrong during `%s`. Better to throw exception",tmp.str().c_str());
        } // else do nothing, workers should not terminate anything
    }



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

