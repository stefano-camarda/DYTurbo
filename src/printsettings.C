#include "printsettings.h"
#include "settings.h"
#include "interface.h"
#include "coupling.h"

#include <LHAPDF/LHAPDF.h>

#include <iostream>
#include <iomanip>

void printsettings()
{
  cout << "===================  Process settings ==================" << endl;
  cout << endl;
  if ((density_.ih1_==1) && (density_.ih2_ == -1))
    cout << "   proton-antiproton collisions";
  else if((density_.ih1_==1) && (density_.ih2_ == 1))
    cout << "   proton-proton collisions";

  cout << " at sqrt(s) = " << energy_.sroot_ << " GeV " << endl;

  if (opts.fixedorder)
    {
      if (nnlo_.order_ == 1)
	cout << "   Computing NLO fixed order cross section" << endl;
      else if (nnlo_.order_ == 2)
	cout << "   Computing NNLO fixed order cross section" << endl;
    }
  else
    {
      if (nnlo_.order_ == 1)
	cout << "   Computing NLO+NLL resummed cross section" << endl;
      else if (nnlo_.order_ == 2)
	cout << "   Computing NNLO+NNLL resummed cross section" << endl;
    }

  if (nproc_.nproc_ == 1)
    cout << "   W+ ->  nu(p3)+e+(p4) production " << endl;
  else if(nproc_.nproc_ == 2)
    cout << "   W- ->  e^-(p3)+nu~(p4) production " << endl;
  else if(nproc_.nproc_== 3 && !opts.useGamma)
    cout << "   Z ->  e-(p3)+e+(p4) production " << endl;
  else if(nproc_.nproc_== 3 && opts.useGamma)
    cout << "   Z/gamma* ->  e-(p3)+e+(p4) production " << endl;
  
  cout << endl;
  if (nproc_.nproc_ == 1 || nproc_.nproc_ == 2)
    {
      cout << "========================  CKM matrix =======================" << endl;
      cout << endl;
      cout << "     | Vud Vus Vub |   | " << setw(9) << cabib_.Vud_ << setw(9) <<  cabib_.Vus_ << setw(9)  << cabib_.Vub_ << " | " << endl;
      cout << "     | Vcd Vcs Vcb | = | " << setw(9) << cabib_.Vcd_ << setw(9) <<  cabib_.Vcs_ << setw(9)  << cabib_.Vcb_ << " | " << endl;
      cout << "     | Vtd Vts Vtb |   | " << setw(9) <<          0  << setw(9) <<            0 << setw(9)  <<           0 << " | " << endl;
      cout << endl;
    }
  cout << "========================  EW parameters ====================" << endl;
  cout << setprecision(15) << endl;
  cout << setw(12) << "Gf ="      << setw(19) << ewcouple_.Gf_    << setw(7) << "GeV^-2" << endl;
  cout << setw(12) << "1/alpha =" << setw(19) << 1/coupling::aemmz     << setw(7) << endl;
  cout << setw(12) << "sin(thW)^2 ="      << setw(19) << ewcouple_.xw_    << setw(7) << endl;
  cout << setw(12) << "W mass ="  << setw(19) << dymasses_.wmass_ << setw(7) << "GeV" << endl;
  cout << setw(12) << "Z mass ="  << setw(19) << dymasses_.zmass_ << setw(7) << "GeV" <<endl;
  cout << setw(12) << "W width =" << setw(19) << dymasses_.wwidth_ << setw(7) << "GeV" << endl;
  cout << setw(12) << "Z width =" << setw(19) << dymasses_.zwidth_ << setw(7) << "GeV" <<endl;
  cout << setw(12) << "Narrow width:" << setw(19) << opts.zerowidth << endl;
  cout << setprecision(6) << endl;
  cout << "========================  QCD settings ======================" << endl;
  cout << endl;
  cout << setw(35) << "Renormalization scale: mur ="  << setw(19) << scale_.scale_ << setw(7) << "GeV" << endl;
  cout << setw(35) << "Factorization scale:   muf ="  << setw(19) << facscale_.facscale_ << setw(7) << "GeV" <<endl;
  cout << setw(35) << "Resummation scale:  m_ll/Q ="  << setw(19) << a_param_.a_param_ << endl;
  cout << setw(35) << "alpha_s(MZ) ="                 << setw(19) << couple_.amz_ << endl;
  cout << setw(35) << "alpha_s(mur) ="                << setw(19) << qcdcouple_.as_ << endl;
  cout << setw(35) << "alpha_s running order ="       << setw(19) << (LHAPDF::getOrderAlphaS()+1) << "-loop" << endl;
  cout << setw(35) << "Non-perturbative form factor: g =" << setw(19) << g_param_.g_param_ << setw(7) << "GeV^2" <<endl;
  cout << endl;
  cout << "========================  PDF settings ======================" << endl;
  cout << endl;
  cout << setw(30) << "PDF set:"                     << setw(20) << opts.LHAPDFset << endl;
  cout << setw(30) << "PDF member:"                  << setw(20) << opts.LHAPDFmember << endl;
  cout << setw(30) << "PDF errors:"                  << setw(20) << (opts.PDFerrors ? "true" : "false") << endl;
  if (opts.doRES)
    {
      if (opts_.approxpdf_ == 1)
	cout << setw(30) << "Mellin moments of PDFs:"     << setw(20) << "Polynomial interpolation" << endl;
      else
	{
	  cout << setw(30) << "Mellin moments of PDFs:"   << setw(20) << "Numerical integration" << endl;
	  cout << setw(30) << "nodes for the quadrature:" << setw(20) << 64 << endl;
	  cout << setw(30) << "Intervals for the quadrature:" << setw(20) << opts_.pdfintervals_ << endl;
	}
    }
  cout << endl;
  cout << "========================  Integration settings ====================" << endl;
  cout << endl;
  cout << setw(25) << "Random seed:"      << setw(30) << opts.rseed << endl;
  cout << setw(25) << "Parallel cores:"   << setw(30) << opts.cubacores << endl;
  if (opts.doRES)
    {
      if (opts.resintvegas)
	cout << setw(25) << "Resummed:"      << setw(30) << "vegas" << setw(12) << "ncalls =" << setw(12) << opts.vegasncallsRES << endl;
      else if (opts.resint3d)
	cout << setw(25) << "Resummed:"      << setw(30) << "cuhre in m,y,pt" << setw(12) << "iter =" << setw(12) << opts.niterRES << endl;
      else if (opts.resint2d)
	{
	  cout << setw(25) << "Resummed:"    << setw(30) << "cuhre in dm,dpt" << setw(12) << "iter =" << setw(12) << opts.niterRES << endl;
	  cout << setw(25) << "Resummed:"    << setw(30) << "gaussian in dy"  << setw(12) << "nodes ="<< setw(12) << opts.yrule << setw(15) << "intervals ="<< setw(5) << opts.yintervals << endl;
	}
      if (opts_.approxpdf_ == 0)
	{
	  cout << setw(25) << "Mellin inverse transform:"    << setw(30) << "gaussian" << setw(12) << "nodes ="<< setw(12) << opts.mellinrule << setw(15) << "intervals ="<< setw(5) << opts.mellinintervals << endl;
	  cout << setw(25) << "zmax in the imaginary axis:"  << setw(30) << opts.zmax << endl;
	}
    }

  if (opts.doVV)
    cout << setw(25) << "Double virtual:"      << setw(30) << "vegas" << setw(12) << "ncalls =" << setw(12) << opts.vegasncallsVV << endl;
  
  if (opts.doCT)
    if (opts.ctintvegas6d)
      cout << setw(25) << "Counter term:"      << setw(30) << "vegas 6d" << setw(12) << "ncalls =" << setw(12) << opts.vegasncallsCT << endl;
    else if (opts.ctintvegas8d)
      cout << setw(25) << "Counter term:"      << setw(30) << "vegas 8d" << setw(12) << "ncalls =" << setw(12) << opts.vegasncallsCT << endl;
    else if (opts.ctint3d)
      cout << setw(25) << "Counter term:"      << setw(30) << "cuhre in dm,dy,dpt" << setw(12) << "iter =" << setw(12) << opts.niterCT << endl;
    else if (opts.ctint2d)
      {
	cout << setw(25) << "Counter term:"      << setw(30) << "cuhre in dm,dy" << setw(12) << "iter =" << setw(12) << opts.niterCT << endl;
	cout << setw(25) << "Counter term:"      << setw(30) << "gaussian in dpt" << setw(12) << "nodes ="<< setw(12) << 20 << setw(15) << "intervals ="<< setw(5) << 1 << endl;
      }

  if (opts.doLO)
    cout << setw(25) << "Z+j LO:"      << setw(30) << "vegas" << setw(12) << "ncalls =" << setw(12) << opts.vegasncallsLO << endl;

  if (opts.doVIRT)
    cout << setw(25) << "Z+j NLO virt:"      << setw(30) << "vegas" << setw(12) << "ncalls =" << setw(12) << opts.vegasncallsVIRT << endl;

  if (opts.doREAL)
    cout << setw(25) << "Z+j NLO real:"      << setw(30) << "vegas" << setw(12) << "ncalls =" << setw(12) << opts.vegasncallsREAL << endl;

  if (opts.doRES && (opts.resint3d || opts.resint2d)
      || opts.doCT && (opts.ctint3d || opts.ctint3d))
    if (opts.cubaint)
      cout << setw(25) << "Angular variables:"      << setw(30) << "Suave in dcosth,dphi" << setw(12) << "ncalls =" << setw(12) << opts.suavepoints << endl;
    else if (opts.quadint)
      {
	cout << setw(25) << "Angular variable costh:"      << setw(30) << "semi-analytical" << setw(12) << "ncstart =" << setw(12) << opts.ncstart << endl;
	cout << setw(25) << "Angular variables phi:"       << setw(30) << "gaussian"        << setw(12) << "intervals =" << setw(12) << opts.quadnphi << endl;
      }
    else if (opts.trapezint)
      {
	cout << setw(25) << "Angular variables costh:"     << setw(30) << "semi-analytical" << setw(12) << "ncstart =" << setw(12) << opts.ncstart << endl;
	cout << setw(25) << "Angular variables phi:"       << setw(30) << "trapezoidal"  << setw(12) << "points =" << setw(12) << opts.nphitrape << endl;
      }
  cout << endl;
  cout << "========================  qt recoil ====================" << endl;
  cout << endl;
  if (opts.qtrec_cs)
    cout << setw(25) << "qt-recoil prescription:"      << setw(30) << "Collins-Soper" << setw(40) << "kt1 = p_V(x)/2; kt2 = p_V(y)/2" << endl;
  else if (opts.qtrec_kt0)
    cout << setw(25) << "qt-recoil prescription:"      << setw(30) << "kt0" << setw(40) << "kt1 = kt2 = 0" << setw(12) << opts.nphitrape << endl;
  else if (opts.qtrec_naive)
    cout << setw(25) << "qt-recoil prescription:"      << setw(30) << "Naive" << setw(40) << "kt = kt_CS * (1+p_V(z)/(m+E_V))" << endl;
  cout << endl;
  cout << "======================== Kinematic settings ====================" << endl;
  cout << endl;
  cout << setw(25) << "y range:"      << setw(12) << opts.ylow << setw(12) << opts.yhigh << endl;
  cout << setw(25) << "m range:"      << setw(12) << opts.mlow << setw(12) << opts.mhigh << endl;
  cout << setw(25) << "apply lepton cuts:"      << setw(20) << (opts.makelepcuts ? "true" : "false") << endl;
  if (opts.makelepcuts)
    if (opts.fiducial == 0)
      {
	cout << setw(25) << "pt_l > " << setw(6) << opts.lptcut << endl;
	cout << setw(25) << "|eta_l| < " << setw(6) << opts.lycut << endl;
	cout << setw(25) << "pt_l(1st) > " << setw(6) << opts.l1ptcut << endl;
	cout << setw(25) << "|eta_l|(1st) < " << setw(6) << opts.l1ycut << endl;
	cout << setw(25) << "pt_l(2nd) > " << setw(6) << opts.l1ptcut << endl;
	cout << setw(25) << "|eta_l|(2nd) < " << setw(6) << opts.l1ycut << endl;
      }
    else
      cout << setw(25) << "fiducial cuts:"      << setw(20) << opts.fiducial << endl;
  if (!opts.fixedorder && (opts.doRES || opts.doCT))
    {
      cout << endl;
      cout << "======================== Resummation damping ====================" << endl;
      cout << endl;
      cout << setw(20) << "damping mode:"  << setw(10) << opts.dampmode << endl;
      cout << setw(20) << "damp above:" << setw(10) << opts.dampk*opts.rmass  << setw(7) << "GeV" << endl;
      cout << setw(20) << "damping width:" << setw(10) << opts.dampdelta*opts.rmass  << setw(7) << "GeV" << endl;
      cout << setw(20) << "resummation cutoff:" << setw(10) << opts.qtcutoff*1000  << setw(7) << "MeV" << endl;
    }
  cout << endl;
  cout << "======================== Debug settings ====================" << endl;
  cout << endl;
  cout << setw(25) << "timeprofile:"      << setw(30) << (opts.timeprofile ? "true" : "false") << endl;
  cout << setw(25) << "verbose:"          << setw(30) << (opts.verbose ? "true" : "false") << endl;
  cout << setw(25) << "cubaverbosity:"    << setw(30) << opts.cubaverbosity << endl;
  cout << endl;
  cout << "============================================================" << endl;
}
