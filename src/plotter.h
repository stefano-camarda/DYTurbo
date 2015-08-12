#ifndef plotter_h
#define plotter_h

#include "config.h"

#ifdef USEROOT
#include <TH1D.h>
#include <TGraphErrors.h>
#endif // USEROOT

class plotter {
    public :
        plotter();
        ~plotter();

        enum TermType { Resum, CT, LO, Real, Virt };

        void Init();
        void FillEvent(double p3[4], double p4[4], double wgt);
        void FillResult(TermType term, double binlo, double binhi, double qt_val, double qt_error, double time);
        void Dump();
        void Finalise();
        static int *gcounter;

    private :
#ifdef USEROOT
        TH1D * h_l1_pt;
        TGraphErrors qt_resum ;
        TGraphErrors qt_ct    ;
        TGraphErrors qt_lo    ;
        TGraphErrors qt_real  ;
        TGraphErrors qt_virt  ;
#endif // USEROOT

};

extern plotter hists;


#endif //plotter_h

