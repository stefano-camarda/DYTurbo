#ifndef plotter_h
#define plotter_h

#include "config.h"

#ifdef USEROOT
#include <TH1D.h>
#include <TH2D.h>
#endif // USEROOT

class plotter {
    public :
        plotter();
        ~plotter();

        enum TermType { Resum, CT, LO, Real, Virt, Total };

        void Init();
        void FillEvent(double p3[4], double p4[4], double wgt);
        void FillResult(TermType term, double int_val, double int_error, double time);
        void Dump();
        void Finalise(double xsection=0);
        static int *gcounter;

    private :
#ifdef USEROOT
        /// @todo: use one object instead
        double N;
        TH1D * h_l1_pt;
        TH1D * h_qt;
        TH2D* qt_y_resum ;
        TH2D* qt_y_ct    ;
        TH2D* qt_y_lo    ;
        TH2D* qt_y_real  ;
        TH2D* qt_y_virt  ;
        TH2D* qt_y_total ;
#endif // USEROOT

};

extern plotter hists;


#endif //plotter_h

