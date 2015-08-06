#include "plotter.h"

plotter hists;

#ifdef USEROOT

plotter::plotter() :
    h_l1_pt(0)
{

    return;
}

plotter::~plotter(){
    if (!h_l1_pt) delete h_l1_pt;
    return;
}


void plotter::init(){
    h_l1_pt = new TH1D ("l1_pt", "lep1 pt", 10, 0,100 );
    h_l1_pt->Sumw2();
    return;
}

void plotter::fill(double p3[4], double p4[4], double wgt){
    double l1_pt = sqrt( pow(p3[0],2) + pow(p3[1],2) );
    h_l1_pt->Fill(l1_pt,wgt);
    //printf(" plotter says: filling histogram %p entries %f value %f", h_l1_pt, h_l1_pt->GetEntries(), l1_pt );
    return;
}

void plotter::finalise(){
    printf(" ploter says: h_l1_pt %p entries %f mean %f RMS %f integral %f \n"
            , h_l1_pt
            , h_l1_pt->GetEntries()
            , h_l1_pt->GetMean()
            , h_l1_pt->GetRMS()
            , h_l1_pt->Integral()
            );
    return;
}

#else // not USEROOT

plotter::plotter(){
    return;
}

plotter::~plotter(){
    return;
}

void plotter::init(){
    return;
}

void plotter::fill(double p3[4], double p4[4], double wgt){
    return;
}

void plotter::finalise(){
    return;
}
#endif //USEROOT
