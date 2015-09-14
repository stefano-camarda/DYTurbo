#include "plotter.h"
#include "integr.h"
#include "interface.h"
#include "settings.h"

#include <sys/mman.h>
#include <sys/types.h>
#include <cmath>
//#include <sys/wait.h>
//#include <unistd.h>

plotter hists;
int * plotter::gcounter;

#ifdef USEROOT
#include "TFile.h"
#include "TString.h"

plotter::plotter() :
    N(0),
    h_l1_pt(0),
    h_qt(0),
    h_y (0),
    h_qtVy (0),
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
    if (!h_l1_pt ) delete h_l1_pt;
    if (!h_qt    ) delete h_qt;
    if (!h_y     ) delete h_y ;
    if (!h_qtVy  ) delete h_qtVy ;
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
    h_qt    = new TH1D ("h_qt", "VB qt", 200, 0,100 );
    h_y     = new TH1D ("h_y" , "VB y" , 25, 0,5 );
    h_qtVy  = new TH2D ("h_qtVy" , "VB qtVy", 200, 0,100  , 25, 0,5 );

    // create shared pointer
    sh_N = make_shared<int>(0);
    //
    /// Trying to calculate final number of vegas calls
    // consider this part as new function: when I will run all terms at once it can be recalculated
    if (opts.doRES  ) N = opts.vegasncallsRES  ;
    if (opts.doCT   ) N = opts.vegasncallsCT   ;
    if (opts.doREAL ) N = opts.vegasncallsREAL ;
    if (opts.doVIRT ) N = opts.vegasncallsVIRT ;
    N = 0;
    return;
}

void plotter::FillEvent(double p3[4], double p4[4], double wgt){
    if (wgt == 0) return;
    //*gcounter = (*gcounter)+1;
    //double l1_pt = sqrt( pow(p3[0],2) + pow(p3[1],2) );
    //h_l1_pt->Fill(l1_pt,wgt);
    // commenting out because of total number of vegas entries
    //N++; 
    if (!N) wgt/=N;
    double qt = sqrt((float)pow(p3[0]+p4[0],2)+pow(p3[1]+p4[1],2));
    double y = 0.5 *log(float((p3[3]+p4[3] +p3[2]+p4[2]) / (p3[3]+p4[3] -p3[2]-p4[2])));
    h_qt   ->Fill(qt,wgt);
    h_y    ->Fill(y,wgt);
    h_qtVy ->Fill(qt,y,wgt);
    return;
}

// fortran interface to Root
void hists_fill_(double p3[4], double p4[4], double weight){
    hists.FillEvent(p3,p4,weight);
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
    printf(" ploter says count %f, shared count %d (use %ld) : h_qt %p entries %f mean %f RMS %f integral %f \n"
            , N
            , *sh_N
            , sh_N.use_count() //*gcounter
            , h_qt
            , h_qt->GetEntries()
            , h_qt->GetMean()
            , h_qt->GetRMS()
            , h_qt->Integral()
            );
    return;
}

void plotter::Merge(){
    // shared pointer from std doesn't work for fork!'
    //std::lock_guard<std::mutex> lock(m);
    //(*sh_N) += N;
    //Dump();
}

void plotter::Finalise(double xsection){
    //shared memory
    //
    //sh_N.reset(); //release from main
    //
    ///@note: Sending the process id number instead of xsection. This is for
    ///       saving separate file name per each worker.
    double intpart; // < this will serve as file index
    bool isworker = ( std::modf(xsection,&intpart) == 0);
    // change the N to serve as normalizing factor: integral of all bins would
    // be normalized to total xsection
    if(!isworker && xsection!=0.) N = h_qt->Integral()/xsection;
    //munmap(gcounter,sizeof *gcounter); //< HAVETO DO -- otherwise you need to reboot
    if (N!=0){
        h_qt   ->Scale(1./N);
        h_y    ->Scale(1./N);
        h_qtVy ->Scale(1./N);
    }
    // correct qt to binwidth
    if(false) for (int ibin=0; ibin<=h_qt->GetNbinsX()+1; ibin++ ){
        double width = h_qt->GetBinWidth   (ibin);
        double val   = h_qt->GetBinContent (ibin);
        double err   = h_qt->GetBinError   (ibin);
        h_qt->SetBinContent (ibin, val/width);
        h_qt->SetBinError   (ibin, err/width);
    }
    // create result dir
    // get name from options
    TString outfname = "results";
    if (isworker) outfname+=int(intpart);
    outfname+= ".root";
    // open file
    TFile *outf = TFile::Open(outfname.Data(), "recreate");
    // write
    h_qt   -> Write();
    h_y    -> Write();
    h_qtVy -> Write();
    // results
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

void plotter::Finalise(double xsection){
    return;
}
#endif //USEROOT
