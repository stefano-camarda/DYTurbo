#ifndef HistoHandler_H
#define HistoHandler_H

/**
 * @file HistoHandler.h
 * Description of this macro
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-08-27
 */

#include <vector>
using std::vector;

// TODO: C++11 ?
//#include <memory>
//using std::unique_ptr;

#include <string>
using std::string;

#include <sys/wait.h>
#include <unistd.h>



#define PARENT_PROC -1

extern "C" {
    // interface to histoHander can go to interface.h
    void histo_setpdf(int *imember);
    void histo_fill(double *l1,double *l2,double *wgt);
    void histo_filldipole(double *l1,double *l2, double *wgt);
    void histo_fillreal();
}

namespace HistoHandler {
    // Typedef and forward declarations
    class HistoBase;
    typedef vector<HistoBase*> HistoContainer;

    // Main functions
    void Init();
    void Clear();
    void Book();
    void DeleteHists();
    void SetVariation(int imember);
    void FillResult(double int_val,double int_err);
    void FillEvent(double *l1,double *l2, double wgt);
    void FillDipole(double *l1,double *l2, double wgt);
    void FillRealEvent();
    void Save(int iworker=PARENT_PROC);
    void Terminate(int iworker=PARENT_PROC);

    // Container of all histograms and profiles
    extern HistoContainer hists;
    // Mode of integration
    extern bool isFillMode;
    // Output file withou suffix
    extern string result_filename;
    // pid for main thread
    extern int parent_pid;

    template<class T>
    void Add(T* newhist){
        //hists.push_back(unique_ptr<HistoBase>(newhist));
        hists.push_back(newhist);
    }

    // Variations treatment.
    struct KeySuffix{
        string suffix;
        size_t index;
    };
}



#endif /* ifndef HistoHandler_H */
