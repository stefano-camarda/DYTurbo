#include "ctint.h"

#include "qtint.h"
#include "mesq.h"
#include "parton.h"
#include "settings.h"
#include "omegaintegr.h"
#include "pdf.h"
#include "phasespace.h"
#include "abint.h"
#include "resconst.h"
#include "dyres_interface.h"
#include "switch.h"
#include "isnan.h"

#include <iostream>
#include <string.h>
#include <math.h>

using namespace std;


//Counterterm to be subtracted from V+j to get a finite cross section at qt->0
void ctint::calc(double costh, double m, double qt, double y, int mode, double f[])
{

  for (int npdf = 0; npdf < opts.totpdf; npdf++)
    f[npdf] = 0.;
      

  ///////////////////////////////////////////////////////
  double m2 = m*m;
  double qt2 = qt*qt;
  double exppy = exp(y);
  double expmy = 1./exppy;
  double tau = sqrt(m2/pow(opts.sroot,2));

  //amplitudes are set in qtint::calc(), no need to calculate them
  /*
  mesq::setpropagators(m);
  double cthmom0, cthmom1, cthmom2;
  if (mode == 0)
    {
      //apply here lepton cuts
      cthmom0 = 1.;
      cthmom1 = costh;
      cthmom2 = pow(costh,2);
    }
  else if (mode == 1)
    {
      //phasespace::set_mqtyphi(m, 0., y);
      phasespace::set_qt(0.);
      omegaintegr::genV4p();
      omegaintegr::cthmoments(cthmom0,cthmom1,cthmom2);
    }

  mesq::setmesq(cthmom0, cthmom1, cthmom2);
  //cout << cthmom0 << "  " << cthmom1 << "  " << cthmom2 << endl;
  */

  //calculate Bjorken x1 x2
  double x1 = tau*exppy;
  double x2 = tau*expmy;

  if (x1 >= 1 || x2 >= 1)
    return;

  //Set factorization scale
  double muf, mur;
  if (opts.dynamicscale)
    {
      muf = m*opts.kmufac;
      mur = m*opts.kmuren;
      double mur2 = mur*mur;
      scaleset_(mur2); //set renormalization and factorization scales, and calculate ason2pi and ason4pi
    }
  else
    {
      muf = opts.rmass*opts.kmufac;
      mur = opts.rmass*opts.kmuren;
    }

  //a-parameter of the resummation scale, set it for the dynamic case
  if (opts.dynamicresscale)
    a_param_.a_param_ = 1./opts.kmures;
  
  //PDFs
  double fx1[2*MAXNF+1],fx2[2*MAXNF+1];
  fdist_(opts.ih1,x1,muf,fx1);
  fdist_(opts.ih2,x2,muf,fx2);
  //////////////////////////////////////////////////////////

  double LR, LF, LQ;
  if (opts.order >= 2)
    LR = log(m2/pow(mur,2));
  LF = log(m2/pow(muf,2));
  LQ = 2.*log(a_param_.a_param_);

  /*
  //Calculate large logs LL1-LL4
  double LL1,LL2,LL3,LL4;
  if (mode == 0 || mode == 1)
    {
      //xmio is used in besselkfast for Itilde
      xmio_.xmio_ = sqrt(qt2/(m2/pow(a_param_.a_param_,2)));
         
      //LL1,LL2,LL3,LL4: large log (squared) corresponding to eq. (136) 
      //In this way normalization is fixed to dsigma/dqt2

      //depends on xmio, which depends on qt2 and m2      
      LL1 = itilde_(one)/pow(m2,2)*pow(a_param_.a_param_,2);
      LL2 = itilde_(two)/pow(m2,2)*pow(a_param_.a_param_,2);
      LL3 = itilde_(three)/pow(m2,2)*pow(a_param_.a_param_,2);
      LL4 = itilde_(four)/pow(m2,2)*pow(a_param_.a_param_,2);
    }
  */

  //Scaled momentum fractions
  double lx1, lx2;
  double pqqintx1, pqqintx2;
  double d0intx1, d0intx2;
  double d1intx1, d1intx2;
  if (opts.order >= 1)
    {
      lx1 = log(x1);
      lx2 = log(x2);
      pqqintx1 = pqqint_(x1);
      pqqintx2 = pqqint_(x2);
    }
  if (opts.order >= 2)
    {
      d0intx1 = d0int_(x1);
      d0intx2 = d0int_(x2);
      d1intx1 = d1int_(x1);
      d1intx2 = d1int_(x2);
    }
  
  //Preliminary loop for caching
  double cz1[abint::abdim];
  double cz2[abint::abdim];
  double oz1[abint::abdim];
  double oz2[abint::abdim];
  // --> cache the logs, and pass them as parameters to the CxxPxx, Pxxxx functions
  double lz1[abint::abdim];
  double lz2[abint::abdim];
  double l1z1[abint::abdim];
  double l1z2[abint::abdim];
  double lzoz1[abint::abdim];
  double lzoz2[abint::abdim];
  if (opts.order >= 1)
    for (int ab = 0; ab < abint::abdim; ab++)
      {
	cz1[ab] = pow(x1,abint::abx[ab]);
	cz2[ab] = pow(x2,abint::abx[ab]);
	oz1[ab] = 1./(1.-cz1[ab]);
	oz2[ab] = 1./(1.-cz2[ab]);
	lz1[ab] = log(cz1[ab]);
	lz2[ab] = log(cz2[ab]);
	l1z1[ab] = log(1.-cz1[ab]);
	l1z2[ab] = log(1.-cz2[ab]);
	lzoz1[ab] = log(1.-cz1[ab])/(1.-cz1[ab]);
	lzoz2[ab] = log(1.-cz2[ab])/(1.-cz2[ab]);
      }

  // skip PDF loop in the preconditioning phase
  int maxpdf=0;
  if (dofill_.doFill_ != 0) maxpdf = opts.totpdf;
      
  // start PDF loop
  for (int npdf = 0; npdf < maxpdf; npdf++)
    {
      dysetpdf_(npdf);
         
      //cache scaled PDFs
      double fx1p[abint::abdim][2*MAXNF+1],fx2p[abint::abdim][2*MAXNF+1];
      if (opts.order >= 1)
	for (int ab = 0; ab < abint::abdim; ab++)
	  {
	    double fx1temp[2*MAXNF+1],fx2temp[2*MAXNF+1];
	    double xx1 = pow(x1,(1-abint::abx[ab]));
	    double xx2 = pow(x2,(1-abint::abx[ab]));
	    fdist_(opts.ih1,xx1,muf,fx1temp);
	    fdist_(opts.ih2,xx2,muf,fx2temp);
	    memcpy(fx1p[ab], fx1temp, (2*MAXNF+1)*sizeof(double));
	    memcpy(fx2p[ab], fx2temp, (2*MAXNF+1)*sizeof(double));
	  }

      double sumfx1p[abint::abdim];
      double sumfx2p[abint::abdim];
      if (opts.order >= 2)
	for (int ab = 0; ab < abint::abdim; ab++)
	  {
	    sumfx1p[ab] = 0.;
	    sumfx2p[ab] = 0.;
	    for (int f = 0; f < MAXNF; f++)
	      {
		sumfx1p[ab] += fx1p[ab][f]+fx1p[ab][parton::charge_conj(parton::pdgid(f))];
		sumfx2p[ab] += fx2p[ab][f]+fx2p[ab][parton::charge_conj(parton::pdgid(f))];
	      }
	  }
	   
      // Start calculation
      double asopi = qcdcouple_.ason2pi_*2.;
  
      //loop on born subprocesses, i.e. born incoming partons ij
      //double lumi[mesq::totpch];
      double sig11[mesq::totpch];
      double sig12[mesq::totpch];
      double sig21[mesq::totpch];
      double sig22[mesq::totpch];
      double sig24[mesq::totpch];
      double sig23[mesq::totpch];
      for (int sp = 0; sp < mesq::totpch; sp++)
	{
	  //simplify notation
	  //double bornmesqij = real(mesq::mesqij[sp]); //born level amplitudes
	  parton::pdgid i = mesq::pid1[sp];         //parton 1
	  parton::pdgid j = mesq::pid2[sp];         //parton 2
	  parton::pdgid g = parton::G;              //gluon
	  parton::pdgid im = parton::charge_conj(i);
	  parton::pdgid jm = parton::charge_conj(j);

	  //LO term (there is no counterterm at LO...)
	  //Simplest term without convolutions
	  double tdelta = fx1[i]*fx2[j];

	  //NLO terms
	  double th1st = 0; //this piece is used only at NNLO?
	  double th1stF = 0;
	  double th1stQ = 0; //this piece is used only at NNLO
	  //H1st delta term
	  th1st += 2*resconst::C1qqdelta*fx1[i]*fx2[j];
      
	  //add resummation scale dependence
	  th1stQ += -(resconst::B1q+resconst::A1q/2.*LQ)*LQ*fx1[i]*fx2[j]; //this piece is used only at NNLO
      
	  //alfa loop (first leg)
	  for (int a = 0; a < abint::abdim; a++)
	    {
	      if (cz1[a] >= 1) continue;
	  
	      //H1st non delta terms
	      th1st += (fx1p[a][i]*cqq_(cz1[a])+fx1p[a][g]*cqg_(cz1[a]))*(-lx1)*fx2[j] * abint::abw[a];

	      //H1st muf dependence, gammaqq and gammaqg:
	      th1stF += (-lx1*((fx1p[a][i]-fx1[i]*cz1[a])*pqq_(cz1[a])+fx1p[a][g]*dypqg_(cz1[a])))*fx2[j] * abint::abw[a];
	    }
	  th1stF += -pqqintx1*fx1[i]*fx2[j];
      
	  //beta loop (second leg)
	  for (int b = 0; b < abint::abdim; b++)
	    {
	      if (cz2[b] >= 1) continue;

	      //H1st non delta terms
	      th1st += (fx2p[b][j]*cqq_(cz2[b])+fx2p[b][g]*cqg_(cz2[b]))*(-lx2)*fx1[i] * abint::abw[b];

	      //H1st muf dependence, gammaqq and gammaqg:
	      th1stF += (-lx2*((fx2p[b][j]-fx2[j]*cz2[b])*pqq_(cz2[b])+fx2p[b][g]*dypqg_(cz2[b])))*fx1[i] * abint::abw[b];
	    }
	  th1stF += -pqqintx2*fx1[i]*fx2[j];
      

	  sig12[sp] = -0.5*resconst::A1q*tdelta;
	  sig11[sp] = -(resconst::B1q+resconst::A1q*LQ)*tdelta - th1stF;

	  if (opts.order == 1) continue;
	  //end NLO

	  //NNLO terms
	  double tcga = 0;
	  double tgamma2 = 0;
	  double tgaga = 0;
      
	  //alfa loop
	  double diffg1f = 0;
	  double diffg10 = 0;
	  double diff1 = 0;
	  double diffc1f = 0;
	  double diffc10 = 0;

	  for (int a = 0; a < abint::abdim; a++)
	    {
	      if (cz1[a] >= 1) continue;
               
	      //(gamma+gamma)*(gamma+gamma) term

	      //First part: one gamma for each leg
	      diffg1f += (-lx1*(fx1p[a][i]-fx1[i]*cz1[a])*pqq_(cz1[a]) - pqqintx1*fx1[i]) * abint::abw[a];
	      diffg10 += -lx1*fx1p[a][g]*dypqg_(cz1[a]) * abint::abw[a];

	      //Second part: gamma*gamma terms
	      //Pij * Pjk = D1ijjk (log(1-z)/(1-z))_+ + D0ijjk/(1-z)_+ 
	      //          + Pijjk(z) + Deltaijjk delta(1-z)
	      //First leg
	      diff1 += (-lx1*((fx1p[a][i]-fx1[i]*cz1[a])
			      *(resconst::D0qqqq*oz1[a]+resconst::D1qqqq*lzoz1[a])
			      +fx1p[a][i]*pqqqq_(cz1[a])+fx1p[a][g]*(pqqqg_(cz1[a])+pqggg_(cz1[a])))
			+(resconst::Deltaqqqq-resconst::D0qqqq*d0intx1-resconst::D1qqqq*d1intx1)
			*fx1[i]) * abint::abw[a];

	      //Include Pqggq
	      diff1 += -lx1*sumfx1p[a]*pqggq_(cz1[a]) * abint::abw[a];
	      //End of (gamma+gamma)*(gamma+gamma) term

	      //Start  (C+C)*(gamma+gamma) term
	      //C first leg, gamma second leg
	      diffc1f += (-lx1*fx1p[a][i]*cqq_(cz1[a])+resconst::C1qqdelta*fx1[i]) * abint::abw[a];
	      diffc10 += -lx1*fx1p[a][g]*cqg_(cz1[a]) * abint::abw[a];
      
	      //C*gamma: first leg (ignore delta term in Cqq: taken into account with th1stF)
	      tcga += (fx1p[a][i]*cqqpqq_(cz1[a])+fx1p[a][g]*(cqqpqg_(cz1[a])+cqgpgg_(cz1[a])))*(-lx1)*fx2[j] * abint::abw[a];

	      //Add Cqg*Pgq contribution
	      tcga += sumfx1p[a]*cqgpgq_(cz1[a])*(-lx1)*fx2[j] * abint::abw[a];

	      //Start 2-loop AP
	      // Gluon + pure singlet
	      //f == gluon piece
	      tgamma2 += fx1p[a][g]*p2qg_(cz1[a])*(-lx1)*fx2[j] * abint::abw[a];

	      //f != gluon piece
	      tgamma2 += sumfx1p[a]*p2qqs_(cz1[a])*(-lx1)*fx2[j] * abint::abw[a];

	      //P2qq non-singlet: regular part
	      tgamma2 += fx1p[a][i]*p2qqv_(cz1[a])*(-lx1)*fx2[j] * abint::abw[a];

	      //P2qq non-singlet: 1/(1-z)_+
	      tgamma2 += 2./3.*resconst::Kappa*(-lx1*(fx1p[a][i]-fx1[i]*cz1[a])*oz1[a]-d0intx1*fx1[i])*fx2[j] * abint::abw[a];

	      //P2qqb non singlet
	      tgamma2 += fx1p[a][im]*p2qqbv_(cz1[a])*(-lx1)*fx2[j] * abint::abw[a];
	    }

	  //beta loop
	  double diffg2f = 0;
	  double diffg20 = 0;
	  double diff2 = 0;
	  double diffc2f = 0;
	  double diffc20 = 0;
	  for (int b = 0; b < abint::abdim; b++)
	    {
	      if (cz2[b] >= 1) continue;
	  
	      //(gamma+gamma)*(gamma+gamma) term

	      //First part: one gamma for each leg
	      diffg2f += (-lx2*(fx2p[b][j]-fx2[j]*cz2[b])*pqq_(cz2[b]) - pqqintx2*fx2[j]) * abint::abw[b];
	      diffg20 += -lx2*fx2p[b][g]*dypqg_(cz2[b]) * abint::abw[b];
	  
	      //Second part: gamma*gamma terms
	      //Pij * Pjk = D1ijjk (log(1-z)/(1-z))_+ + D0ijjk/(1-z)_+ 
	      //          + Pijjk(z) + Deltaijjk delta(1-z)
	      //Second leg
	      diff2 += (-lx2*((fx2p[b][j]-fx2[j]*cz2[b])
			      *(resconst::D0qqqq*oz2[b]+resconst::D1qqqq*lzoz2[b])
			      +fx2p[b][j]*pqqqq_(cz2[b])+fx2p[b][g]*(pqqqg_(cz2[b])+pqggg_(cz2[b])))
			+(resconst::Deltaqqqq-resconst::D0qqqq*d0intx2-resconst::D1qqqq*d1intx2)
			*fx2[j])* abint::abw[b];
	      //Include Pqggq
	      diff2 += -lx2*sumfx2p[b]*pqggq_(cz2[b]) * abint::abw[b];
	      //End of (gamma+gamma)*(gamma+gamma) term

	      //Start  (C+C)*(gamma+gamma) term
	      //gamma first leg, C second leg
	      diffc2f += (-lx2*fx2p[b][j]*cqq_(cz2[b])+resconst::C1qqdelta*fx2[j]) * abint::abw[b];
	      diffc20 += -lx2*fx2p[b][g]*cqg_(cz2[b]) * abint::abw[b];

	      //C*gamma: second leg (ignore delta term in Cqq: taken into account with th1stF)
	      tcga += (fx2p[b][j]*cqqpqq_(cz2[b])+fx2p[b][g]*(cqqpqg_(cz2[b])+cqgpgg_(cz2[b])))*(-lx2)*fx1[i] * abint::abw[b];

	      //Add Cqg*Pgq contribution
	      tcga += sumfx2p[b]*cqgpgq_(cz2[b])*(-lx2)*fx1[i] * abint::abw[b];

	      //Start 2-loop AP
	      // Gluon + pure singlet
	      //f == gluon piece
	      tgamma2 += fx2p[b][g]*p2qg_(cz2[b])*(-lx2)*fx1[i] * abint::abw[b];

	      //f != gluon piece
	      tgamma2 += sumfx2p[b]*p2qqs_(cz2[b])*(-lx2)*fx1[i] * abint::abw[b];

	      //P2qq non-singlet: regular part
	      tgamma2 += fx2p[b][j]*p2qqv_(cz2[b])*(-lx2)*fx1[i] * abint::abw[b];

	      //P2qq non-singlet: 1/(1-z)_+
	      tgamma2 += 2./3.*resconst::Kappa*(-lx2*(fx2p[b][j]-fx2[j]*cz2[b])*oz2[b]-d0intx2*fx2[j])*fx1[i] * abint::abw[b];

	      //P2qqb non singlet
	      tgamma2 += fx2p[b][jm]*p2qqbv_(cz2[b])*(-lx2)*fx1[i] * abint::abw[b];
	    }

	  tgaga=tgaga+2*(diffg10*diffg20+diffg1f*diffg2f+diffg10*diffg2f+diffg1f*diffg20);
	  tgaga += diff1*fx2[j];
	  tgaga += diff2*fx1[i];

	  tcga += (diffc10*diffg20+diffc1f*diffg2f+diffc10*diffg2f+diffc1f*diffg20);
	  tcga += (diffg10*diffc20+diffg1f*diffc2f+diffg10*diffc2f+diffg1f*diffc20);
      
	  sig24[sp] = pow(resconst::A1q,2) / 8. * tdelta;
	  sig23[sp] = -resconst::beta0*resconst::A1q/3.*tdelta-0.5*resconst::A1q*sig11[sp];
	  sig22[sp] = 0.5*(resconst::beta0*resconst::A1q*(LR-LQ)-resconst::A2q)*tdelta
	    -0.5*resconst::A1q*(th1st+th1stQ+(LF-LQ)*th1stF)
	    -0.5*(resconst::B1q+resconst::A1q*LQ-resconst::beta0)*sig11[sp]
	    +0.5*(resconst::B1q+resconst::A1q*LQ)*th1stF
	    +0.5*tgaga;
	  sig21[sp] =
	    -resconst::beta0*(LR-LQ)*sig11[sp]
	    -(resconst::B1q+resconst::A1q*LQ)*(th1st+th1stQ+(LF-LQ)*th1stF)
	    -(LF-LQ)*tgaga
	    -(resconst::B2q+resconst::A2q*LQ)*tdelta
	    +resconst::beta0*th1st
	    +(resconst::B1q+0.5*resconst::A1q*LQ)*LQ*th1stF
	    -tcga-tgamma2;
	    
	  //include missing delta term from C*gamma (no factor 2 here)
	  sig21[sp] += -resconst::C1qqdelta*th1stF;
	  //include missing term from contact term in 2 loop AP
	  sig21[sp] += -2*resconst::Delta2qq*tdelta;
	}

      double xmsq = 0.;
      for (int sp = 0; sp < mesq::totpch; sp++)
	{
	  //as/pi factor
	  double sig1 = (sig12[sp]*qtint::LL2_mesqij[sp]+sig11[sp]*qtint::LL1_mesqij[sp])*asopi;
	  xmsq += -sig1;
	  
	  if (opts.order == 1) continue;
	  //(as/pi)^2 factor
	  double sig2 = (sig24[sp]*qtint::LL4_mesqij[sp]+sig23[sp]*qtint::LL3_mesqij[sp]+sig22[sp]*qtint::LL2_mesqij[sp]+sig21[sp]*qtint::LL1_mesqij[sp])*pow(asopi,2);

	  // sum O(as) and O(as^2) contributions	
	  xmsq += -sig2;
	}	  
      
      double shad = pow(opts.sroot,2);

      //double fbGeV2=0.38937966e12;
      //double fac = M_PI*6./fbGeV2*(2*m2);
      //double flux = fbGeV2/(2.*x1*x2*shad);
      //double ps = 1./shad;
      //double norm = 1./16./M_PI / 2. / M_PI;
      //xmsq = xmsq * ps * fac * flux * norm;

      xmsq = xmsq/shad * 3./8. /2./M_PI;
      //xmsq = xmsq/shad * 3./8. /2./M_PI/2./M_PI;

      xmsq = xmsq * m2; //Is this a Jacobian for dsigma/dq2 -> dsigma/dq?
      xmsq = xmsq * 2*M_PI; //phiV integration
      if (isnan_ofast(xmsq))
	cout << m << " " << y << " " << costh << "  " << xmsq << endl;
  
      //switching --> switching function is inside qtint, do not apply
      //double swtch;
      //if (mode == 0 || mode == 1)
      //  swtch = switching::swtch(qt, m);
      //else if (mode == 2)
      //  swtch=1.; // qt integration already performed

      //if (swtch < 0.01) return 0.;// do not apply this cut to avoid discontinuities. Instead the phase space is limited to qt and m switching limits
      xmsq = xmsq;//*swtch; //switching function is inside qtint

      f[npdf] = xmsq;
    } //end loop on pdf

  return;
}
