#include "plotter.h"

#include <sys/mman.h>
#include <sys/types.h>
//#include <sys/wait.h>
//#include <unistd.h>

plotter hists;
int * plotter::gcounter;

#ifdef USEROOT
#include "TFile.h"

plotter::plotter() :
    h_l1_pt(0)
{

   /// @todo: memory sharing
   // void * m_addr = NULL;
   // size_t m_size = sizeof *gcounter;
   // int m_prot = PROT_READ | PROT_WRITE ;
   // int m_flags = MAP_SHARED | MAP_ANONYMOUS;
   // int m_fildes = -1; // posix type memory allocation
   // off_t m_offset = 0; // to be zero!

   //  gcounter = (int *) mmap( m_addr , m_size, m_prot, m_flags, m_fildes, m_offset);

    return;
}

plotter::~plotter(){
    if (!h_l1_pt) delete h_l1_pt;
    return;
}


void plotter::Init(){

    *gcounter = 0;

    qt_resum .SetName ( "qt_resum" );
    qt_ct    .SetName ( "qt_ct"    );
    qt_lo    .SetName ( "qt_lo"    );
    qt_real  .SetName ( "qt_real"  );
    qt_virt  .SetName ( "qt_virt"  );

    h_l1_pt = new TH1D ("l1_pt", "lep1 pt", 10, 0,100 );
    h_l1_pt->Sumw2();
    return;
}

void plotter::FillEvent(double p3[4], double p4[4], double wgt){
    *gcounter = (*gcounter)+1;
    double l1_pt = sqrt( pow(p3[0],2) + pow(p3[1],2) );
    h_l1_pt->Fill(l1_pt,wgt);
    return;
}

void plotter::FillResult(TermType term, double binlo, double binhi, double int_val, double int_error, double time){
    TGraphErrors * gr;
    switch (term) {
        case Resum : gr = & qt_resum ; break;
        case CT    : gr = & qt_ct    ; break;
        case LO    : gr = & qt_lo    ; break;
        case Real  : gr = & qt_real  ; break;
        case Virt  : gr = & qt_virt  ; break;
    }
    int n = gr->GetN();
    double qt_val   = (binhi+binlo)/2;
    double qt_error = (binhi-binlo)/2;
    gr->SetPoint      ( n , qt_val   , int_val   );
    gr->SetPointError ( n , qt_error , int_error );
    return;
}

void plotter::Dump(){
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

void plotter::Finalise(){
    //munmap(gcounter,sizeof *gcounter); //< HAVETO DO -- otherwise you need to reboot
    
    // create result dir
    // get name from options
    const char * outfname = "results.root";
    // open file
    TFile *outf = TFile::Open(outfname, "recreate");
    // write
    qt_resum .Write();
    qt_ct    .Write();
    qt_lo    .Write();
    qt_real  .Write();
    qt_virt  .Write();
    // close
    outf->Write();
    outf->Close();
    return;
}

#else // not USEROOT

plotter::plotter(){
    return;
}

plotter::~plotter(){
    return;
}

void plotter::Init(){
    return;
}

void plotter::FillEvent(double p3[4], double p4[4], double wgt){
    return;
}

void plotter::FillResult(TermType term, double binlo, double binhi, double qt_val, double qt_error, double time){
    return;
}

void plotter::Dump(){
    return;
}

void plotter::Finalise(){
    return;
}
#endif //USEROOT
