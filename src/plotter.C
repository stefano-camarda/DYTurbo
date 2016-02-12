#include "plotter.h"

plotter hists;

#include "integr.h"
#include "interface.h"
#include "settings.h"

#include <sys/mman.h>
#include <sys/types.h>
#include <cmath>
#include <algorithm>


#ifdef USEROOT
#include <TFile.h>
#include <TH1.h>
#include <TString.h>
#include <TLorentzVector.h>

plotter::plotter() :
    N           (0),
    h_qt        (0),
    h_y         (0),
    h_qtVy      (0),
    qt_y_resum  (0),
    qt_y_ct     (0),
    qt_y_lo     (0),
    qt_y_real   (0),
    qt_y_virt   (0),
    qt_y_total  (0),
    last_npdf   (-1),
    doAiMoments (true),
    //ai_maarten ("ai_maarten"),
    verbose     (false)
{
    return;
}

plotter::~plotter(){
    for (int i=0; i< h_qtVy_PDF.size(); i++){
        SetPDF(i);
        if(!h_qt   ) delete h_qt   ;
        if(!h_y    ) delete h_y    ;
        if(!h_qtVy ) delete h_qtVy ;
        if(doAiMoments) {
            for (auto i=0; i<NMOM; i++){
                delete pa_qtVy .A[i]; pa_qtVy .A[i]=NULL;
                delete pa_qt   .A[i]; pa_qt   .A[i]=NULL;
                delete pa_y    .A[i]; pa_y    .A[i]=NULL;
            }
            if(!h_costh   ) delete h_costh   ;
            if(!h_phi     ) delete h_phi     ;
            if(!h_phi_lep ) delete h_phi_lep ;
        }
    }
    if(!qt_y_resum ) delete qt_y_resum ;
    if(!qt_y_ct    ) delete qt_y_ct    ;
    if(!qt_y_lo    ) delete qt_y_lo    ;
    if(!qt_y_real  ) delete qt_y_real  ;
    if(!qt_y_virt  ) delete qt_y_virt  ;
    if(!qt_y_total ) delete qt_y_total ;
    return;
}


void plotter::Init(){
    if (verbose) printf(" plotter:  histogram initialization\n");
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
    vector<double> bins_prof_qt = { 0, 2.5, 5.0, 8.0, 11.4, 14.9, 18.5, 22.0, 25.5, 29.0, 32.6, 36.4, 40.4, 44.9, 50.2, 56.4, 63.9, 73.4, 85.4, 105.0, 132.0, 173.0, 253.0, 600.0 };

    // differential xsection and ai profiles
    h_qt      = new TH1D        ("h_qt"     , "VB qt"   , bins_qt [0] , bins_qt [1] , bins_qt [2] );
    h_y       = new TH1D        ("h_y"      , "VB y"    , bins_y  [0] , bins_y  [1] , bins_y  [2] );
    h_qtVy    = new TH2D        ("h_qtVy"   , "VB qtVy" , bins_qt [0] , bins_qt [1] , bins_qt [2] ,  bins_y [0] , bins_y[1] , bins_y[2] );
    // turn off aimomens
    if (doAiMoments){
        for (auto i=0; i<NMOM; i++){
            TString name = "p_qtVy_A";
            name += i;
            TString title = "A"; title+= i; title+= " qt vs y;qt[GeV];y;A_{"; title+=i; title+="}";
            pa_qtVy.A[i] = new TProfile2D(name.Data(), title.Data(), bins_qt [0] , bins_qt [1] , bins_qt [2] ,  bins_y [0] , bins_y[1] , bins_y[2] );
            //
            name = "p_qt_A"; name+=i;
            title = "A"; title+= i; title+= " qt;qt[GeV];A_{"; title+=i; title+="}";
            pa_qt.A[i] = new TProfile(name.Data(), title.Data(), bins_prof_qt.size()-1, &bins_prof_qt[0] );
            //
            name = "p_y_A"; name+=i;
            title = "A"; title+= i; title+= " y;y;A_{"; title+=i; title+="}";
            pa_y.A[i] = new TProfile(name.Data(), title.Data(),   bins_y [0] , bins_y[1] , bins_y[2] );
        }
        h_costh       = new TH1D        ("h_costh"      , "VB costh"    , 100, -1,1 );
        h_phi         = new TH1D        ("h_phi"        , "VB phi"      , 100, 0,TMath::TwoPi() );
        h_phi_lep     = new TH1D        ("h_phi_lep"    , "Lep phi CM"  , 100, 0,TMath::TwoPi() );
    }
    // correction factor for Ai moments
    c[0] = 20./3. ;
    c[1] = 5.     ;
    c[2] = 10.    ;
    c[3] = 4.     ;
    c[4] = 4.     ;
    c[5] = 5.     ;
    c[6] = 5.     ;
    c[7] = 4.     ;
    // force to one -- cuba gives correct weights (accoring to number of entries);
    N = 1;
    // vector of weights -- for systematics
    v_wgt.clear();
    SetPDF(0);

    //ai_maarten.Initialize();

    return;
}

bool plotter::IsInitialized(){return (h_qtVy!=0);};

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
    return (p[3] + p[2])        ;
}
double Vminus(double p[4]){
    return (p[3] - p[2])        ;
}
double calcCosThCS(double Q2,double qt,double p3[4],double p4[4]){
    double qz = p3[2]+p4[2];
    if (qz==0) return 0;
    double costh=0;
    costh = (Vplus(p3)*Vminus(p4) - Vplus(p4)*Vminus(p3));
    costh *= qz<0. ? -1 : 1;
    costh /= sqrt(float(Q2*(Q2 + qt*qt)));
    return costh;
}
double calcPhiCS(double p3[4],double p4[4],double &phi_lep){
    // Collins-Sopper phi
    const double pmass=0.938;
    /// @todo: phi_CS without TLorentzVector...
    // Define TLorentz Vectors
    TLorentzVector lep1(p3);
    TLorentzVector VB,p1,p2;
    VB.SetPxPyPzE(p3[0]+p4[0], p3[1]+p4[1], p3[2]+p4[2], p3[3]+p4[3]);
    double sign = VB.Z()>0 ? 1 : -1;
    p1.SetXYZM(0.,0., +sign*opts.sroot/2., pmass);
    p2.SetXYZM(0.,0., -sign*opts.sroot/2., pmass);
    // Boost
    TVector3 VBboost = -VB.BoostVector();
    p1   .Boost(VBboost);
    p2   .Boost(VBboost);
    lep1 .Boost(VBboost);
    phi_lep=TVector2::Phi_0_2pi(-lep1.Phi()+TMath::Pi()/2.); // phi_lep: cos(-phi+pi/2) = sin(phi)
    // Collins Soper transversal plane
    TVector3 hatp1, hatp2, CSAxis, xAxis, yAxis;
    hatp1 = p1.Vect().Unit();
    hatp2 = p2.Vect().Unit();
    CSAxis = ( hatp1 - hatp2      ).Unit();
    yAxis  = ( hatp1.Cross(hatp2) ).Unit();
    xAxis  = ( yAxis.Cross(CSAxis)).Unit();
    // calculate phi
    return TMath::ATan2(float(lep1.Vect()*yAxis), float(lep1.Vect()*xAxis));
}

void plotter::CalculateKinematics(double p3[4], double p4[4]){
    // calculate VB
    double p[4]; 
    p[0] = p3[0]+p4[0];
    p[1] = p3[1]+p4[1];
    p[2] = p3[2]+p4[2];
    p[3] = p3[3]+p4[3];
    //
    Q2 = calcQ2(p);
    qt = calcQt(p);
    y  = calcY(p);
    if (!doAiMoments) return; // turn off ai moments
    // calc Collins-Sopper
    costh = calcCosThCS(Q2,qt,p3,p4);
    phi   = calcPhiCS(p3,p4,phi_lep);
    // angular terms for moments
    double theta     = TMath::ACos(costh);
    double sintheta  = TMath::Sin(theta);
    double sin2theta = TMath::Sin(2*theta);
    double cosphi    = TMath::Cos(phi);
    double cos2phi   = TMath::Cos(2*phi);
    double sinphi    = TMath::Sin(phi);
    double sin2phi   = TMath::Sin(2*phi);
    // define ai moments
    a[0] = c[0] * (0.5-1.5*costh*costh       ) +2./3.;
    a[1] = c[1] * (sin2theta*cosphi          )       ;
    a[2] = c[2] * (sintheta*sintheta*cos2phi )       ;
    a[3] = c[3] * (sintheta*cosphi           )       ;
    a[4] = c[4] * (costh                     )       ;
    a[5] = c[5] * (sintheta*sintheta*sin2phi )       ;
    a[6] = c[6] * (sin2theta*sinphi          )       ;
    a[7] = c[7] * (sintheta*sinphi           )       ;
    a[8] = a[0] - a[2];
}


void plotter::print_dipole(XsecPoint pt){
    printf("   ibin %d pt %f y %f wgt %f\n", pt.ibin, pt.qt, pt.y, pt.wgt );
}
void plotter::print_dipoleVec(std::vector<XsecPoint> vec ){
    printf("  beg vec: size %lu\n",vec.size());
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
    if (nd == 5 || nd==6 ) { // use 0 dipole kinematics
        qt = dipole_points[0] .qt;
        y  = dipole_points[0] .y;
    } else { // calculate pt and y for 0..4
        CalculateKinematics(p3,p4);
    }
    // fill Xsec - point
    point.ibin = h_qtVy->FindBin(qt,y);
    point.qt    = qt    ;
    point.y     = y     ;
    point.wgt   = wgt   ;
    point.costh = costh ;
    point.phi   = phi   ;
    point.fid   = decide_fiducial(p3,p4);
    dipole_points.push_back(point);
    // fill moments
    if(doAiMoments) for (auto i=0; i<NMOM; i++){
        pa_qtVy .A[i]->Fill(qt,y, a[i],wgt);
        pa_qt   .A[i]->Fill(qt,   a[i],wgt);
        pa_y    .A[i]->Fill(qt,   a[i],wgt);
    } 
    return;
}

void plotter::FillRealEvent(plotter::TermType term){
    // process all calculated contributions
    while (!dipole_points.empty()){
        // go per each affected qt_y bin
        point = dipole_points.back(); // take last one and combare ibin
        dipole_points.pop_back();
        for (auto ipoint=dipole_points.begin(); ipoint!=dipole_points.end() ; ) if (point.ibin == ipoint->ibin){
            // found same bin: remove from list and sum weight
            point.wgt+=ipoint->wgt;
            ipoint = dipole_points.erase(ipoint);
        } else ++ipoint; // bin index not same, skip it
        // fill histograms in ibin with sum of all weights
        if (point.fid ){
            h_qt   -> Fill(point.qt         ,point.wgt);
            h_y    -> Fill(point.y          ,point.wgt);
            h_qtVy -> Fill(point.qt,point.y ,point.wgt);
            if(doAiMoments){
                h_costh   ->Fill( point.costh                    ,point.wgt);
                h_phi     ->Fill( TVector2::Phi_0_2pi(point.phi) ,point.wgt);
            }
        }
    }
    return;
}

void plotter::FillEvent(double p3[4], double p4[4], double wgt){
    // save after 10M events
    //if(int(h_qt->GetEntries()+1)%int(1e6)==0) Finalise(0);
    //
    if (wgt == 0 ) return;
    // TLorentzVector lep1(p3);
    // TLorentzVector lep2(p4);
    // if (verbose){
    //     lep1.Print();
    //     lep2.Print();
    // }
    // ai_maarten.Execute(
    //         lep1.Pt(), lep1.Eta(), lep1.Phi(), 0., 13,
    //         lep2.Pt(), lep2.Eta(), lep2.Phi(), 0., -13,
    //         7e3, wgt
    //         );
    CalculateKinematics(p3,p4);
    if (verbose){
        printf(" plotter:  histogram filling: \n");
        printf("         p3: %g %g %g %g \n", p3[0], p3[1], p3[2], p3[3] );
        printf("         p4: %g %g %g %g \n", p4[0], p4[1], p4[2], p4[3] );
        printf("         pV: %g %g       \n", qt, y  );
    }
    // dividing weight
    if (N!=0) wgt/=N;
    h_qt      ->Fill( qt        ,wgt);
    h_y       ->Fill( y         ,wgt);
    h_qtVy    ->Fill( qt,y      ,wgt);
    if(doAiMoments){
        for (int i=0;i<NMOM;i++) pa_qtVy.A[i] ->Fill( qt,y,a[i] ,wgt);
        for (int i=0;i<NMOM;i++) pa_qt.A[i] ->Fill( qt,a[i] ,wgt);
        for (int i=0;i<NMOM;i++) pa_y.A[i] ->Fill( y,a[i] ,wgt);
        h_costh   ->Fill( costh                    ,wgt);
        h_phi     ->Fill( TVector2::Phi_0_2pi(phi) ,wgt);
        h_phi_lep ->Fill( phi_lep                  ,wgt);
    }
    return;
}

// fortran interface to Root
void hists_fill_(double p3[4], double p4[4], double *weight){
    hists.FillEvent(p3,p4,*weight);
    return;
}

void hists_AiTest_(double pjet[4][12], double p4cm[4],double *Q,double *qt,double *y,double* pcosthCS, double* pphiCS, double *pphiVB, double *wt, double *loHst ){
    if ((*wt)*(*loHst)==0) return;
    //printf("EVENT WEIGHT %g %g\n", *wt, *loHst);
    //printf("PJET test -- px -- py -- pz -- E\n");
    //printf("0 : %f\t%f\t%f\t%f \n", pjet[0][0], pjet[1][0], pjet[2][0], pjet[3][0] );
    //printf("1 : %f\t%f\t%f\t%f \n", pjet[0][1], pjet[1][1], pjet[2][1], pjet[3][1] );
    //printf("2 : %f\t%f\t%f\t%f \n", pjet[0][2], pjet[1][2], pjet[2][2], pjet[3][2] );
    //printf("3 : %f\t%f\t%f\t%f \n", pjet[0][3], pjet[1][3], pjet[2][3], pjet[3][3] );
    //printf("4 : %f\t%f\t%f\t%f \n", pjet[0][4], pjet[1][4], pjet[2][4], pjet[3][4] );
    //printf("5 : %f\t%f\t%f\t%f \n", pjet[0][5], pjet[1][5], pjet[2][5], pjet[3][5] );
    //printf("------------------------\n");
    double pvb[4]; for (int i=0; i<4; i++) pvb[i]=pjet[i][2]+pjet[i][3];
    //printf("vb: %f\t%f\t%f\t%f \n", pvb[0], pvb[1],  pvb[2], pvb[3]);
    //printf("q : %f\t%f\t%f\t%f \n", calcQ2(pvb), calcQt(pvb), calcY(pvb), 0 );
    //printf("------------------------\n");
    TLorentzVector vb(pvb);
    TLorentzVector lep; lep.SetXYZT(pjet[0][3],pjet[1][3],pjet[2][3],pjet[3][3]);
    lep.Boost(-vb.BoostVector());
    TLorentzVector ilep(p4cm);

    double i_costhCS = *pcosthCS;
    double i_phiCS   = *pphiCS;
    double i_phiVB   = *pphiVB;
    double i_q       = (*Q);
    double i_qt      = (*qt);
    double i_y       = (*y);
    double i_p4cm    = p4cm[3];

    double p_costhCS = hists.costh;
    double p_phiCS   = hists.phi;
    double p_phiVB   = TLorentzVector(pvb).Phi();
    double p_q       = TMath::Sqrt(hists.Q2);
    double p_qt      = hists.qt;
    double p_y       = hists.y;
    double p_p4cm    = lep.E();

    //
    double pihalf = TMath::Pi()/2.;
    double pe = lep.E();
    double theta = lep.Theta();
    double phi = hists.phi_lep; // i_phiCS;
    double px = pe * TMath::Sin(theta) * TMath::Sin(phi);
    double py = pe * TMath::Sin(theta) * TMath::Cos(phi);
    double pz = pe * TMath::Cos(theta);
    TLorentzVector jlep(px,py,pz,pe);

    printf("AI Test -- integrand -- plotter \n");
    printf("p_costhCS = %f \t %f \t %g \n"        , i_costhCS , p_costhCS , i_costhCS - p_costhCS );
    printf("p_phiCS   = %f \t %f \t %g\t%g\t%g \n", i_phiCS   , p_phiCS   , TVector2::Phi_mpi_pi(i_phiCS   - p_phiCS),  TVector2::Phi_mpi_pi(i_phiCS   - p_phiCS + pihalf),  TVector2::Phi_mpi_pi(i_phiCS   - p_phiCS - pihalf)  );
    printf("p_phiVB   = %f \t %f \t %g \n"        , i_phiVB   , p_phiVB   , TVector2::Phi_mpi_pi(i_phiVB   - p_phiVB)   );
    printf("p_q       = %f \t %f \t %g \n"        , i_q       , p_q       , i_q       - p_q       );
    printf("p_qt      = %f \t %f \t %g \n"        , i_qt      , p_qt      , i_qt      - p_qt      );
    printf("p_y       = %f \t %f \t %g \n"        , i_y       , p_y       , i_y       - p_y       );
    printf("p_p4cm    = %f \t %f \t %g \n"        , i_p4cm    , p_p4cm    , i_p4cm    - p_p4cm    );
    printf("lep      ");lep.Print();
    printf("ilep     ");ilep.Print();
    printf("jlep     ");jlep.Print();
    printf("lep-ilep ");(lep-ilep).Print();
    printf("lep-jlep ");(jlep-ilep).Print();
    printf("--------------------------------\n");
    printf("--------------------------------\n");

    return;
}

void hists_setpdf_(int * npdf){
    hists.SetPDF(*npdf);
}

void hists_fill_pdf_(double p3[4], double p4[4], double *weight, int *npdf){
    hists.SetPDF(*npdf);
    hists_fill_(p3,p4,weight);
}

void hists_real_dipole_(double p3[4], double p4[4], double *weight, int * nd){
    hists.FillRealDipole(p3,p4,*weight,*nd);
}

void hists_real_dipole_pdf_(double p3[4], double p4[4], double *weight, int * nd, int* npdf){
    hists.SetPDF(*npdf);
    hists_real_dipole_(p3,p4,weight,nd);
}

void hists_real_event_(){
    hists.FillRealEvent();
}

void hists_real_event_pdf_(int *npdf){
    hists.SetPDF(*npdf);
    hists_real_event_();
}

void hists_finalize_(){
    hists.Finalise();
}

void plotter::CumulateResult(TermType term, double wgt){
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
    h->Fill ( qt, y, wgt);
    return;
}

void plotter::addToBin(TH1*h, double int_val, double int_error){
    // Get bin
    // get the middle of bin
    double qt_val = ( qtmax + qtmin )/2.;
    double y_val  = ( ymax  + ymin  )/2.;
    int ibin = 0;
    if      (!TString(h->GetName()).CompareTo( "h_qtVy" )) ibin = h->FindBin( qt_val, y_val);
    else if (!TString(h->GetName()).CompareTo( "h_qt"   )) ibin = h->FindBin( qt_val  );
    else if (!TString(h->GetName()).CompareTo( "h_y"    )) ibin = h->FindBin( y_val   );
    // if there somethin in bin add it properly
    double val = h->GetBinContent (ibin);
    double err = h->GetBinError   (ibin);
    val += int_val;
    err = sqrt (err*err + int_error*int_error);
    // set new values
    h->SetBinContent (ibin, val);
    h->SetBinError   (ibin, err);
}

void plotter::FillQuadrature(double int_val, double int_error){
    addToBin( h_qt    , int_val , int_error);
    addToBin( h_y     , int_val , int_error);
    addToBin( h_qtVy  , int_val , int_error);
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

TH1 * plotter::clone_PDF(TH1*h,int npdf){
    TH1 *out = 0;
    TString name = h->GetName();
    name += npdf;
    out = (TH1*) h->Clone(name.Data());
    out->Reset();
    out->SetDirectory(0);
    name = h->GetTitle();
    name += " "; name+=npdf;
    out->SetTitle(name);
    return out;
}


void plotter::SetPDF(int npdf){
    // if current npdf still same don't change anything
    if (last_npdf==npdf) return;
    // if empty add current hist to 0-th position
    // assuming we always starting from central
    int size = h_qtVy_PDF.size();
    if (size==0) {
        h_qt_PDF   .push_back( h_qt   );
        h_y_PDF    .push_back( h_y    );
        h_qtVy_PDF .push_back( h_qtVy );
        if (doAiMoments) {
            h_costh_PDF   .push_back( h_costh   );
            h_phi_PDF     .push_back( h_phi     );
            h_phi_lep_PDF .push_back( h_phi_lep );
            pa_qtVy_PDF   .push_back( pa_qtVy   );
            pa_qt_PDF     .push_back( pa_qt     );
            pa_y_PDF      .push_back( pa_y      );
        }
    }
    // create new histograms
    size = h_qtVy_PDF.size();
    while (npdf >= size){
        // clone, rename and reset
        h_qt_PDF   .push_back( (TH1D *) clone_PDF( h_qt_PDF   [0] , size) );
        h_y_PDF    .push_back( (TH1D *) clone_PDF( h_y_PDF    [0] , size) );
        h_qtVy_PDF .push_back( (TH2D *) clone_PDF( h_qtVy_PDF [0] , size) );
        if (doAiMoments){
            h_costh_PDF   .push_back( (TH1D *) clone_PDF( h_costh_PDF   [0] , size) );
            h_phi_PDF     .push_back( (TH1D *) clone_PDF( h_phi_PDF     [0] , size) );
            h_phi_lep_PDF .push_back( (TH1D *) clone_PDF( h_phi_lep_PDF [0] , size) );
            clone_Array_PDF( pa_qtVy_PDF  , size) ;
            clone_Array_PDF( pa_qt_PDF    , size) ;
            clone_Array_PDF( pa_y_PDF     , size) ;
        }
        // update size value
        size = h_qtVy_PDF.size();
    }
    // change poiters to correct pdf histograms
    h_qt      = h_qt_PDF      [npdf];
    h_y       = h_y_PDF       [npdf];
    h_qtVy    = h_qtVy_PDF    [npdf];
    if (doAiMoments){
        h_costh   = h_qt_PDF      [npdf];
        h_phi     = h_phi_PDF     [npdf];
        h_phi_lep = h_phi_lep_PDF [npdf];
        for (int i=0; i<NMOM; i++){
            pa_qtVy .A[i] = pa_qtVy_PDF .at(npdf) .A[i];
            pa_qt   .A[i] = pa_qt_PDF   .at(npdf) .A[i];
            pa_y    .A[i] = pa_y_PDF    .at(npdf) .A[i];
        }
    }
    // set last pdf
    last_npdf = npdf;
    return;
}

void plotter::Dump(){
    printf(" ploter says count %f : h_qt %p entries %f mean %f RMS %f integral %f \n"
            , N
            , h_qt
            , h_qt->GetEntries()
            , h_qt->GetMean()
            , h_qt->GetRMS()
            , h_qt->Integral()
          );
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
    if (verbose){
        printf(" plotter:  histogram finalize\n");
        Dump();
    }
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
    for (int i=0; i< h_qtVy_PDF.size(); i++){
        SetPDF(i);
        h_qt   -> Write();
        h_y    -> Write();
        h_qtVy -> Write();
    }
    if(doAiMoments){
        for (auto pp : pa_qtVy .A) pp->Write();
        for (auto pp : pa_qt   .A) pp->Write();
        for (auto pp : pa_y    .A) pp->Write();
        h_costh   -> Write();
        h_phi     -> Write();
        h_phi_lep -> Write();
    } 
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
    //ai_maarten.Finalize();
    return;
}



#else // not USEROOT

plotter::plotter(){return;}
plotter::~plotter(){return;}
void plotter::Init(){return;}
bool plotter::IsInitialized(){return false;}
void plotter::FillQuadrature(double int_val, double int_error){return;}
void plotter::FillEvent(double p3[4], double p4[4], double wgt){return;}
void plotter::FillRealDipole(double p3[4], double p4[4], double wgt,int nd){return;}
void plotter::FillRealEvent(TermType term ){return;}
void plotter::FillResult(TermType term, double int_val, double int_error, double time){return;}
void plotter::SetPDF(int npdf){return;}
void plotter::CumulateResult(TermType term, double wgt){return;}
void plotter::Merge(){return;}
void plotter::Dump(){return;}
void plotter::Finalise(double xsection){return;}


void hists_setpdf_(int * npdf){ return; }
void hists_fill_(double p3[4], double p4[4], double *weight){ return; }
void hists_real_dipole_(double p3[4], double p4[4], double *weight,int *nd){ return; }
void hists_real_event_(){ return; }
void hists_fill_pdf_(double p3[4], double p4[4], double *weight, int *npdf){ return; }
void hists_real_dipole_pdf_(double p3[4], double p4[4], double *weight,int *nd, int *npdf){ return; }
void hists_real_event_pdf_(int* npdf){ return; }
yid hists_AiTest_(double* pcosthCS, double* pphiCS) {return};

#endif //USEROOT
