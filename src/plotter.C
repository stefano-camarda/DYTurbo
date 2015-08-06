#include "plotter.h"

#include <sys/mman.h>
#include <sys/types.h>
//#include <sys/wait.h>
//#include <unistd.h>

plotter hists;
int * plotter::gcounter;

#ifdef USEROOT


plotter::plotter() :
    h_l1_pt(0)
{

   // memeory sharring
   void * m_addr = NULL;
   size_t m_size = sizeof *gcounter;
   int m_prot = PROT_READ | PROT_WRITE ;
   int m_flags = MAP_SHARED | MAP_ANONYMOUS;
   int m_fildes = -1; // posix type memory allocation
   off_t m_offset = 0; // to be zero!

    gcounter = (int *) mmap( m_addr , m_size, m_prot, m_flags, m_fildes, m_offset);

    return;
}

plotter::~plotter(){
    if (!h_l1_pt) delete h_l1_pt;
    return;
}


void plotter::init(){

    *gcounter = 0;

    h_l1_pt = new TH1D ("l1_pt", "lep1 pt", 10, 0,100 );
    h_l1_pt->Sumw2();
    return;
}

void plotter::fill(double p3[4], double p4[4], double wgt){
    *gcounter = (*gcounter)+1;
    double l1_pt = sqrt( pow(p3[0],2) + pow(p3[1],2) );
    h_l1_pt->Fill(l1_pt,wgt);
    //printf(" plotter says: filling histogram %p entries %f value %f", h_l1_pt, h_l1_pt->GetEntries(), l1_pt );
    return;
}

void plotter::dump(){
    printf(" ploter says count (%d) : h_l1_pt %p entries %f mean %f RMS %f integral %f \n"
            , *gcounter
            , h_l1_pt
            , h_l1_pt->GetEntries()
            , h_l1_pt->GetMean()
            , h_l1_pt->GetRMS()
            , h_l1_pt->Integral()
            );
    return;
}

void plotter::finalise(){
    munmap(gcounter,sizeof *gcounter); //< HAVETO DO -- otherwise you need to reboot
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

void plotter::dump(){
    return;
}

void plotter::finalise(){
    return;
}
#endif //USEROOT
