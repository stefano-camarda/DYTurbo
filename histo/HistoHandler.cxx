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

#include "HistoHandler.h"
#include "Kinematics.h"

#include <map>
using std::map;

#include <sstream>
using std::ostringstream;

using namespace Kinematics;

// interface:
void histo_setpdf(int *imember){
    HistoHandler::SetVariation(*imember);
}
void histo_fill(double *l1,double *l2,double *wgt){
    HistoHandler::FillEvent(l1,l2,*wgt);
}
void histo_filldipole(double *l1,double *l2, double *wgt){
    HistoHandler::FillDipole(l1,l2,*wgt);
}
void histo_fillreal(){
    HistoHandler::FillRealEvent();
}


namespace HistoHandler{

    // Container of all histograms and profiles
    HistoContainer hists;
    // Mode of integration
    bool isFillMode=true;
    // Output file withou suffix
    string result_filename="results";
    // Pid for main thread
    int parent_pid=0;

    void Init() {
        //TODO: Implement to DYTURBO
        //isFillMode = (bins.plotmode == "integrate");
        Book();
        SetVariation(0); // start from 0
        // TODO: Implmentation to dyturbo
        // Initiate as long vector as necessary in beginning : Check if we are
        // doing pdf scanning, if yes then prebook histograms.
        // TODO: Implmentation to dyturbo
        // Check we will use just one variation, then set it from beginning
        parent_pid=getpid();
    }

    void FillEvent(double p3[4],double p4[4], double wgt){
        SetKinematics(p3,p4,wgt);
        for (auto h_it = hists.begin();h_it!=hists.end();h_it++){
            (*h_it)->FillEvent();
        }
    }

    void FillResult(double int_val, double int_err){
        SetMiddlePoint();
        for (auto h_it = hists.begin();h_it!=hists.end();h_it++){
            (*h_it)->AddToBin(int_val,int_err);
        }
    }

    void FillDipole(double p3[4],double p4[4], double wgt){
        SetKinematics(p3,p4,wgt);
        for (auto h_it = hists.begin();h_it!=hists.end();h_it++){
            (*h_it)->FillDipole();
        }
    }

    void FillRealEvent(){
        for (auto h_it = hists.begin();h_it!=hists.end();h_it++){
            (*h_it)->FillRealEvent();
        }
    }

    // Variations treatment.
    // This part is handeling consitent variation treatment of PDF. It is also
    // possible to extend this by other options scale variations in future.
    size_t last_index=0;
    typedef map<int,const KeySuffix> MapSuffix;
    MapSuffix variation_suffixes;

    void SetVariation(int imember){
        // create missing index for variation
        if (variation_suffixes.count(imember)==0){
            // NOTE: if in future we find out that there is negative member, or
            // higher than maximum of LHAPDF memebers we can use this to define
            // scale variations
            KeySuffix addkey = {"",0};
            if (imember==0){ // add central
            } else { // PDF case TODO: Check is PDF memeber i.e is positive and smaller than number of memebrs in pdfset
                ostringstream suffix;
                suffix << "_pdf" << imember;
                addkey.suffix=suffix.str();
                addkey.index=last_index;
                last_index++;
            }
            // add it and increase index
            variation_suffixes.insert(MapSuffix::value_type(imember,addkey));
            last_index++;
        }
        // Pickup correct histogram
        for (auto h_it = hists.begin();h_it!=hists.end();h_it++){
            (*h_it)->SetVariation(variation_suffixes.at(imember));
        }
    }

    // forward declaration : specialized in implementation header file
    void OpenOutputFile(int iworker);
    void CloseOutputFile();
    void Merge(int iworker);

    void Save(int iworker){
        OpenOutputFile(iworker);
        for (auto h_it = hists.begin();h_it!=hists.end();h_it++){
            (*h_it)->Save();
        }
        CloseOutputFile();
    }

    void Terminate(int iworker){
        Merge(iworker);
        // TODO: Properly remove all histograms from your
    }
}


#endif /* ifndef  */

