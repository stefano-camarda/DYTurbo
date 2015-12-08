#include "printsettings.h"
#include "settings.h"
#include "interface.h"
#include "coupling.h"

#include <LHAPDF/LHAPDF.h>

#include <iostream>
#include <iomanip>

void printsettings()
{
  cout << "===================  Kinematical settings ==================" << endl;
  cout << endl;
  if ((density_.ih1_==1) && (density_.ih2_ == -1))
    cout << "   proton-antiproton collisions";
  else if((density_.ih1_==1) && (density_.ih2_ == 1))
    cout << "   proton-proton collisions";

  cout << " at sqrt(s) = " << energy_.sroot_ << " GeV " << endl;

  if (nnlo_.order_ == 1)
    cout << "   Computing NLO(+NLL) cross section" << endl;
  else if (nnlo_.order_ == 2)
    cout << "   Computing NNLO(+NNLL) cross section" << endl;

  if (nproc_.nproc_ == 1)
    cout << "   W+ ->  nu(p3)+e+(p4) production " << endl;
  else if(nproc_.nproc_ == 2)
    cout << "   W- ->  e^-(p3)+nu~(p4) production " << endl;
  else if(nproc_.nproc_==3)
    cout << "   Z ->  e-(p3)+e+(p4) production " << endl;
  cout << endl;
  cout << "========================  CKM matrix =======================" << endl;
  cout << endl;
  cout << "     | Vud Vus Vub |   | " << setw(9) << cabib_.Vud_ << setw(9) <<  cabib_.Vus_ << setw(9)  << cabib_.Vub_ << " | " << endl;
  cout << "     | Vcd Vcs Vcb | = | " << setw(9) << cabib_.Vcd_ << setw(9) <<  cabib_.Vcs_ << setw(9)  << cabib_.Vcb_ << " | " << endl;
  cout << "     | Vtd Vts Vtb |   | " << setw(9) <<          0  << setw(9) <<            0 << setw(9)  <<           0 << " | " << endl;
  cout << endl;
  cout << "========================  EW parameters ====================" << endl;
  cout << setprecision(15) << endl;
  cout << setw(12) << "Gf ="      << setw(19) << ewcouple_.Gf_    << setw(7) << "GeV^-2" << endl;
  cout << setw(12) << "1/alpha =" << setw(19) << 1/coupling::aemmz     << setw(7) << endl;
  cout << setw(12) << "sin(thW)^2 ="      << setw(19) << ewcouple_.xw_    << setw(7) << endl;
  cout << setw(12) << "W mass ="  << setw(19) << dymasses_.wmass_ << setw(7) << "GeV" << endl;
  cout << setw(12) << "Z mass ="  << setw(19) << dymasses_.zmass_ << setw(7) << "GeV" <<endl;
  cout << setw(12) << "W width =" << setw(19) << dymasses_.wwidth_ << setw(7) << "GeV" << endl;
  cout << setw(12) << "Z width =" << setw(19) << dymasses_.zwidth_ << setw(7) << "GeV" <<endl;
  cout << setprecision(6) << endl;
  cout << "========================  QCD settings ======================" << endl;
  cout << endl;
  cout << setw(35) << "Renormalization scale: mur ="  << setw(19) << scale_.scale_ << setw(7) << "GeV" << endl;
  cout << setw(35) << "Factorization scale:   muf ="  << setw(19) << facscale_.facscale_ << setw(7) << "GeV" <<endl;
  cout << setw(35) << "Resummation scale:  m_ll/Q ="  << setw(19) << a_param_.a_param_ << endl;
  cout << setw(35) << "alpha_s(MZ) ="                 << setw(19) << qcdcouple_.as_ << endl;
  cout << setw(35) << "alpha_s running order ="       << setw(19) << (LHAPDF::getOrderAlphaS()+1) << "-loop" << endl;
  cout << setw(35) << "PDF set:"                      << setw(19) << opts.LHAPDFset << endl;
  cout << setw(35) << "PDF member: "                  << setw(19) << opts.LHAPDFmember << endl;
  cout << setw(35) << "Non-perturbative form factor: g =" << setw(19) << g_param_.g_param_ << setw(7) << "GeV^2" <<endl;
  cout << endl;
  cout << "============================================================" << endl;
}
