#include "plotter.h"
#include "integr.h"
#include "settings.h"

#include <sys/mman.h>
#include <sys/types.h>
//#include <sys/wait.h>
//#include <unistd.h>

plotter hists;
int * plotter::gcounter;

#ifdef USEROOT
#include "TFile.h"

plotter::plotter() :
    h_l1_pt(0),
    qt_y_resum (0),
    qt_y_ct    (0),
    qt_y_lo    (0),
    qt_y_real  (0),
    qt_y_virt  (0),
    qt_y_total (0)
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
    if(!qt_y_resum ) delete qt_y_resum ;
    if(!qt_y_ct    ) delete qt_y_ct    ;
    if(!qt_y_lo    ) delete qt_y_lo    ;
    if(!qt_y_real  ) delete qt_y_real  ;
    if(!qt_y_virt  ) delete qt_y_virt  ;
    if(!qt_y_total ) delete qt_y_total ;
    return;
}


void plotter::Init(){

    //*gcounter = 0;
    TH1::SetDefaultSumw2(true);

    int qt_n = bins.hist_qt_bins .size() -1 ;
    int y_n  = bins.hist_y_bins  .size() -1 ;
    double * qt_array = &(bins.hist_qt_bins [0]);
    double * y_array  = &(bins.hist_y_bins  [0]);
    qt_y_resum = new TH2D ( "qt_y_resum" ,"qt_y_resum" , qt_n, qt_array, y_n, y_array );
    qt_y_ct    = new TH2D ( "qt_y_ct"    ,"qt_y_ct"    , qt_n, qt_array, y_n, y_array );
    qt_y_lo    = new TH2D ( "qt_y_lo"    ,"qt_y_lo"    , qt_n, qt_array, y_n, y_array );
    qt_y_real  = new TH2D ( "qt_y_real"  ,"qt_y_real"  , qt_n, qt_array, y_n, y_array );
    qt_y_virt  = new TH2D ( "qt_y_virt"  ,"qt_y_virt"  , qt_n, qt_array, y_n, y_array );
    qt_y_total = new TH2D ( "qt_y_total" ,"qt_y_total" , qt_n, qt_array, y_n, y_array );

    h_l1_pt = new TH1D ("l1_pt", "lep1 pt", 10, 0,100 );
    return;
}

void plotter::FillEvent(double p3[4], double p4[4], double wgt){
    //*gcounter = (*gcounter)+1;
    double l1_pt = sqrt( pow(p3[0],2) + pow(p3[1],2) );
    h_l1_pt->Fill(l1_pt,wgt);
    return;
}

void plotter::FillResult(TermType term, double int_val, double int_error, double time){
    TH2D * h = 0 ;
    switch (term) {
        case Resum : h = qt_y_resum ; break;
        case CT    : h = qt_y_ct    ; break;
        case LO    : h = qt_y_lo    ; break;
        case Real  : h = qt_y_real  ; break;
        case Virt  : h = qt_y_virt  ; break;
        case Total : h = qt_y_total ; break;
    }
    double qt_val = ( qtmax + qtmin )/2.;
    double y_val  = ( ymax  + ymin  )/2.;
    int ibin = h->FindBin(qt_val, y_val);
    h->SetBinContent ( ibin , int_val   );
    h->SetBinError   ( ibin , int_error );
    return;
}

void plotter::Dump(){
    printf(" ploter says count (%d) : h_l1_pt %p entries %f mean %f RMS %f integral %f \n"
            , 0 //*gcounter
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
    qt_y_resum -> Write();
    qt_y_ct    -> Write();
    qt_y_lo    -> Write();
    qt_y_real  -> Write();
    qt_y_virt  -> Write();
    qt_y_total -> Write();
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

void plotter::FillResult(TermType term, double int_val, double int_error, double time){
    return;
}

void plotter::Dump(){
    return;
}

void plotter::Finalise(){
    return;
}
#endif //USEROOT
