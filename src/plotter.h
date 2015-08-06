#ifndef plotter_h
#define plotter_h

#include "config.h"

#ifdef USEROOT
#include <TH1D.h>
#endif // USEROOT

class plotter {
    public :
        plotter();
        ~plotter();

        void init();
        void fill(double p3[4], double p4[4], double wgt);
        void finalise();
    private :
#ifdef USEROOT
        TH1D * h_l1_pt;
#endif

};

extern plotter hists;


#endif //plotter_h

