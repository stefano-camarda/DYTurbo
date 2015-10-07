#include "plotter.h"
#include "integr.h"
#include "interface.h"
#include "settings.h"

#include <sys/mman.h>
#include <sys/types.h>
#include <cmath>
#include <algorithm>
//#include <sys/wait.h>
//#include <unistd.h>

plotter hists;

#ifdef USEROOT
#include "TFile.h"
#include "TString.h"

plotter::plotter() :
    N          (0),
    h_qt       (0),
    h_y        (0),
    h_qtVy     (0),
    p_qtVy_A4  (0),
    qt_y_resum (0),
    qt_y_ct    (0),
    qt_y_lo    (0),
    qt_y_real  (0),
    qt_y_virt  (0),
    qt_y_total (0)
{
    return;
}

plotter::~plotter(){
    if(!h_qt       ) delete h_qt;
    if(!h_y        ) delete h_y        ;
    if(!h_qtVy     ) delete h_qtVy     ;
    if(!p_qtVy_A4  ) delete p_qtVy_A4  ;
    if(!qt_y_resum ) delete qt_y_resum ;
    if(!qt_y_ct    ) delete qt_y_ct    ;
    if(!qt_y_lo    ) delete qt_y_lo    ;
    if(!qt_y_real  ) delete qt_y_real  ;
    if(!qt_y_virt  ) delete qt_y_virt  ;
    if(!qt_y_total ) delete qt_y_total ;
    return;
}


void plotter::Init(){
    // Use sumw2 by default
    TH1::SetDefaultSumw2(true);
    // get qt-y binning: result histograms
    int qt_n = bins.hist_qt_bins .size() -1 ;
    int y_n  = bins.hist_y_bins  .size() -1 ;
    double * qt_array = &(bins.hist_qt_bins [0]);
    double * y_array  = &(bins.hist_y_bins  [0]);
    // histogram for storing result
    qt_y_resum = new TH2D ( "qt_y_resum" ,"qt_y_resum" , qt_n, qt_array, y_n, y_array );
    qt_y_ct    = new TH2D ( "qt_y_ct"    ,"qt_y_ct"    , qt_n, qt_array, y_n, y_array );
    qt_y_lo    = new TH2D ( "qt_y_lo"    ,"qt_y_lo"    , qt_n, qt_array, y_n, y_array );
    qt_y_real  = new TH2D ( "qt_y_real"  ,"qt_y_real"  , qt_n, qt_array, y_n, y_array );
    qt_y_virt  = new TH2D ( "qt_y_virt"  ,"qt_y_virt"  , qt_n, qt_array, y_n, y_array );
    qt_y_total = new TH2D ( "qt_y_total" ,"qt_y_total" , qt_n, qt_array, y_n, y_array );
    // differential qt-y binning
    double bins_qt [] = {200., 0., 100.};
    double bins_y  [] = {25., 0., 5.};
    // differential xsection and ai profiles
    h_qt      = new TH1D        ("h_qt"     , "VB qt"   , bins_qt [0] , bins_qt [1] , bins_qt [2] );
    h_y       = new TH1D        ("h_y"      , "VB y"    , bins_y  [0] , bins_y  [1] , bins_y  [2] );
    h_qtVy    = new TH2D        ("h_qtVy"   , "VB qtVy" , bins_qt [0] , bins_qt [1] , bins_qt [2] ,  bins_y [0] , bins_y[1] , bins_y[2] );
    p_qtVy_A4 = new TProfile2D( "p_qtVy_A4" , "A4 qtVy" , bins_qt [0] , bins_qt [1] , bins_qt [2] ,  bins_y [0] , bins_y[1] , bins_y[2] );

    /// Trying to calculate final number of vegas calls
    // consider this part as new function: when I will run all terms at once it can be recalculated
    if (opts.doRES  ) N = opts.vegasncallsRES  ;
    if (opts.doCT   ) N = opts.vegasncallsCT   ;
    if (opts.doREAL ) N = opts.vegasncallsREAL ;
    if (opts.doVIRT ) N = opts.vegasncallsVIRT ;
    // hacked to one -- cuba gives correct weights (accoring to number of entries);
    N = 1;

    // vector of weights -- for systematics
    v_wgt.clear();
    return;
}

double calcQt(double p[4]){
    return  sqrt((float)   p[0]*p[0]+p[1]*p[1]);
}
double calcY(double p[4]){
    return 0.5 *log(float( (p[3]+p[2]) / (p[3]-p[2]) ));
}
double calcQ2(double p[4]){
    return p[3]*p[3] -p[0]*p[0] -p[1]*p[1] -p[2]*p[2];
}
double sqrt2=1.4142135623730951454746218587388284504413604736328125;
double Vplus(double p[4]){
    return (p[3] + p[2])/sqrt2  ;
}
double Vminus(double p[4]){
    return (p[3] - p[2])/sqrt2  ;
}
void plotter::CalculateKinematics(double p3[4], double p4[4]){
    // calculate VB
    double p[4]; 
    p[0] = p3[0]+p4[0];
    p[1] = p3[1]+p4[1];
    p[2] = p3[2]+p4[2];
    p[3] = p3[3]+p4[3];
    //
    double Q2 = calcQ2(p);
           qt = calcQt(p);
           y  = calcY(p);
    // Collins-Sopper theta
    double costh_CS = 2;
    costh_CS *= sqrt( float((Q2-qt*qt)/Q2) );
    costh_CS *= (Vplus(p3)*Vminus(p4) - Vplus(p4)*Vminus(p3));
    // Collins-Sopper phi
    /// @todo: phi_CS;
    // define ai moments
    /// @todo: rest of moments;
           a4 = costh_CS;
    //
}


void plotter::print_dipole(XsecPoint pt){
    printf("   ibin %d pt %f y %f wgt %f\n", pt.ibin, pt.qt, pt.y, pt.wgt );
}
void plotter::print_dipoleVec(std::vector<XsecPoint> vec ){
    printf("  beg vec: size %d\n",vec.size());
    int i=0;
    for (auto ipt : vec){
        printf ( "   i %d",i);
        print_dipole(ipt);
        i++;
    }
    printf("  end vec: \n");
}


void plotter::FillRealDipole(double p3[4], double p4[4], double wgt, int nd){
    if (nd!=0 && wgt == 0 ) return; // make sure you have at least first one for kinematics
    if (nd == 5 || nd==6 ) { // calculate pt and y for 0..4
        CalculateKinematics(p3,p4);
    } else {
        qt = dipole_points[0] .qt;
        y  = dipole_points[0] .y;
    }
    if (N!=0) wgt/=N;
    // fill Xsec - point
    point.ibin = h_qtVy->FindBin(qt,y);
    point.qt   = qt   ;
    point.y    = y    ;
    point.wgt  = wgt  ;
    dipole_points.push_back(point);
    // debug
    //printf("FillRealDipole: beg\n");
    //print_dipole(point);
    //print_dipoleVec(dipole_points);
    //printf("FillRealDipole: end\n");
    // fill moments
    p_qtVy_A4 -> Fill(qt,y,a4,wgt);
    return;
}

void plotter::FillRealEvent(){
    // process all calculated contributions
    //printf("FillRealEvent: beg\n");
    //print_dipole(point);
    //print_dipoleVec(dipole_points);
    int i=0;
    while (!dipole_points.empty()){
        // go per each affected qt_y bin
        point = dipole_points.back();
        dipole_points.pop_back();
        for (auto ipoint=dipole_points.begin(); ipoint!=dipole_points.end() ; ) if (point.ibin == ipoint->ibin){
            // found same bin: remove from list and sum weight
            point.wgt+=ipoint->wgt;
            ipoint = dipole_points.erase(ipoint);
        } else ++ipoint; // bin index not same, skip it
        //printf (" dipole group %d",i);
        //print_dipole(point);
        //print_dipoleVec(dipole_points);
        // fill histograms
        h_qt   -> Fill(point.qt         ,point.wgt);
        h_y    -> Fill(point.y          ,point.wgt);
        h_qtVy -> Fill(point.qt,point.y ,point.wgt);
        i++;
    }
    //printf("FillRealEvent: end\n");
    //printf("----------\n");
    return;
}

void plotter::FillEvent(double p3[4], double p4[4], double wgt){
    if (wgt == 0 ) return;
    CalculateKinematics(p3,p4);
    // dividing weight
    if (N!=0) wgt/=N;
    h_qt      ->Fill( qt      ,wgt);
    h_y       ->Fill( y       ,wgt);
    h_qtVy    ->Fill( qt,y    ,wgt);
    p_qtVy_A4 ->Fill( qt,y,a4 ,wgt);
    //v_wgt .push_back( wgt );
    return;
}

// fortran interface to Root
void hists_fill_(double p3[4], double p4[4], double *weight){
    //printf("wt %g\np3 %g %g %g %g\np4 %g %g %g %g\n",
            //*weight,
            //p3[0],
            //p3[1],
            //p3[2],
            //p3[3],
            //p4[0],
            //p4[1],
            //p4[2],
            //p4[3]
          //);
    hists.FillEvent(p3,p4,*weight);
    return;
}

void hists_real_dipole_(double p3[4], double p4[4], double *weight, int * nd){
    hists.FillRealDipole(p3,p4,*weight,*nd);
}

void hists_real_event_(){
    hists.FillRealEvent();
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
    // get the middle of bin
    double qt_val = ( qtmax + qtmin )/2.;
    double y_val  = ( ymax  + ymin  )/2.;
    // set bin content
    int ibin = h->FindBin(qt_val, y_val);
    h->SetBinContent ( ibin , int_val   );
    h->SetBinError   ( ibin , int_error );
    return;
}

void plotter::Dump(){
    //printf(" ploter says count %f : h_qt %p entries %f mean %f RMS %f integral %f \n"
            //, N
            //, h_qt
            //, h_qt->GetEntries()
            //, h_qt->GetMean()
            //, h_qt->GetRMS()
            //, h_qt->Integral()
            //);
    //dump weights
    std::vector<double>::iterator v_b = v_wgt.begin ();
    std::vector<double>::iterator v_e = v_wgt.end   ();
    int Nw = v_wgt.size();
    double dNw   = Nw==0 ? 0.0 : Nw;
    double sumw  = Nw==0 ? 0.0 : std::accumulate(v_b, v_e, 0.0);
    double sumw2 = Nw==0 ? 0.0 :  std::inner_product(v_b, v_e, v_b, 0.0);
    double meanw = Nw==0 ? 0.0 : sumw / dNw;
    double rmsw  = Nw==0 ? 0.0 : std::sqrt(sumw2/dNw - meanw*meanw );
    double minw  = Nw==0 ? 0.0 : (* std::min_element(v_b,v_e) ) ;
    double maxw  = Nw==0 ? 0.0 : (* std::max_element(v_b,v_e) ) ;
    printf(" ploter says count %d, sumw %g, sumw2 %g, meanw %g, RMSw %g, minw %g, maxw %g \n",
            Nw,
            sumw,
            sumw2,
            meanw,
            rmsw,
            minw,
            maxw
          );
    return;
}

void plotter::Merge(){
    //Dump();
}

void plotter::Finalise(double xsection){
    //Dump();
    ///@note: Sending the process id number instead of xsection. This is for
    ///       saving separate file name per each worker.
    double intpart; ///< this will work as file index
    bool isworker = ( std::modf(xsection,&intpart) == 0);
    // now we change the N to serve as normalizing factor: integral of all bins
    // would be normalized to total xsection
    if(!isworker && xsection!=0.) N = h_qt->Integral()/xsection;
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
    // tree
    //t_tree->Write();
    
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

void hists_fill_(double p3[4], double p4[4], double *weight){
    return;
}
#endif //USEROOT
