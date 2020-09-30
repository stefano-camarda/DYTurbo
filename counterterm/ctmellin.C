#include "ctmellin.h"

#include "qtint.h"
#include "mesq.h"
#include "anomalous.h"
#include "pdfevol.h"
#include "evolnative.h"
#include "parton.h"
#include "settings.h"
#include "omegaintegr.h"
#include "pdf.h"
#include "scales.h"
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

using namespace anomalous;


//Counterterm to be subtracted from V+j to get a finite cross section at qt->0
void ctmellin::calc(double m, double f[])
{

  for (int npdf = 0; npdf < opts.totpdf; npdf++)
    f[npdf] = 0.;
      

  ///////////////////////////////////////////////////////
  double m2 = m*m;
  //  double qt2 = qt*qt;
  //  double exppy = exp(y);
  //  double expmy = 1./exppy;

  //Set scales
  scales::set(m);
  scales::mcfm();
  double muf = scales::fac;
  double mur = scales::ren;

  //update PDFs in Mellin space at the starting scale, if the factorisation scale is proportional to mll
  //!!! create a special function/module to deal with non-evolved PDFs
  if (opts.fmufac > 0)
    {
      evolnative::allocate();
      evolnative::update();
    }
  pdfevol::allocate_fx();
  
  /*
  //Set factorization scale
  double muf, mur;
  if (opts.dynamicscale)
    {
      //Cannot use dynamic scale, throw an error!!!
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
  */

  //a-parameter of the resummation scale, set it for the dynamic case
  if (opts.fmures > 0)
    a_param_.a_param_ = 1./opts.kmures;
  
  //for fixed resummation scale need to recompute a_param
  else
    a_param_.a_param_ = m/scales::res;
  //////////////////////////////////////////////////////////

  double LR, LF, LQ;
  if (opts.order >= 2)
    LR = log(m2/pow(mur,2));
  LF = log(m2/pow(muf,2));
  LQ = 2.*log(a_param_.a_param_);

  // skip PDF loop in the preconditioning phase
  int maxpdf=0;
  if (dofill_.doFill_ != 0) maxpdf = opts.totpdf;
      
  // Start calculation
  double H1q = 0.;

  //factors to compensate the as/2/pi normalization (a factor of 2 for each power of alphas)
  double beta0N = resconst::beta0*2.;
  double beta1N = resconst::beta1*4.;
  double beta2N = resconst::beta2*8.;
  double A1gN = resconst::A1g*2.;
  double A2gN = resconst::A2g*4.;
  double A3gN = resconst::A3g*8.;
  double B1gN = resconst::B1g*2.;
  double B2gN = resconst::B2g*4.;
  double A1qN = resconst::A1q*2.;
  double A2qN = resconst::A2q*4.;
  double A3qN = resconst::A3q*8.;
  double B1qN = resconst::B1q*2.;
  double B2qN = resconst::B2q*4.;
      
  double asopi = qcdcouple_.ason2pi_*2.;

  double bjx = m2/pow(opts.sroot,2);
  double ax = log(bjx);

  //double lumi[mesq::totpch];
  double sig11[mesq::totpch] = {0.};
  double sig12[mesq::totpch] = {0.};
  double sig21[mesq::totpch] = {0.};
  double sig22[mesq::totpch] = {0.};
  double sig23[mesq::totpch] = {0.};
  double sig24[mesq::totpch] = {0.};
  
  //Start mellin loop
  for (int n = 0; n < mellinint::mdim; n++)
    {
      //retrieve PDFs in Mellin space at the factorisation scale
      pdfevol::retrievemuf(n);
      complex<double>* fn1 = pdfevol::fn1;
      complex<double>* fn2 = pdfevol::fn2;

      int ni = anomalous::index(n,mesq::positive);

      //Mellin transform
      complex <double> cexp = exp(-mellinint::Np[n] * ax)/M_PI * mellinint::CCp/complex <double>(0.,1);

      //      cout << mellinint::mdim << "  " << n << "  " << mellinint::Np[n] << "  " << fn1[4] << "  " << fn2[4] << endl;
      
      //loop on born subprocesses, i.e. born incoming partons ij
      complex<double> s11,s12,s21,s22,s23,s24;     
      for (int sp = 0; sp < mesq::totpch; sp++)
	{
	  //simplify notation
	  //double bornmesqij = real(mesq::mesqij[sp]); //born level amplitudes
	  parton::partid i = mesq::pidn1[sp];         //parton 1
	  parton::partid j = mesq::pidn2[sp];         //parton 2
	  parton::partid g = parton::g;              //gluon
	  parton::partid im = parton::charge_conjn(i);
	  parton::partid jm = parton::charge_conjn(j);

      	  //LO term (there is no counterterm at LO...)
	  //Simplest term without convolutions
	  complex <double> tdelta = fn1[i]*fn2[j];
	  
	  complex <double> th1stF = 2.*gamma1qq[ni]*tdelta + gamma1qg[ni]*(fn1[i]*fn2[g]+fn1[g]*fn2[j]);

	  s12 = -0.5*A1qN*tdelta;
	  s11 = -(B1qN+A1qN*LQ)*tdelta - th1stF;
      
	  sig11[sp] += real(s11 * cexp * mellinint::wn[n]);
	  sig12[sp] += real(s12 * cexp * mellinint::wn[n]);
	  //sig11[sp] += real(tdelta * cexp * mellinint::wn[n])*real(mesq::mesqij[sp]);
 
	  //end NLO
	  if (opts.order == 1) continue;
	  
	  
	  complex <double> th1stQ = -(B1qN+A1qN*LQ/2.)*LQ*tdelta;
	  complex <double> th1st = 2.*C1QQ[ni]*tdelta + C1QG[ni]*(fn1[i]*fn2[g]+fn1[g]*fn2[j]);
	  
	  s24 = pow(A1qN,2)/8.*tdelta;
	  s23 = -beta0N*A1qN/3.*tdelta-0.5*A1qN*s11;
	   
	  complex <double> tgaga = 
	    +4.*pow(gamma1qq[ni],2)*tdelta
	    +3.*gamma1qq[ni]*gamma1qg[ni]*(fn1[i]*fn2[g]+fn1[g]*fn2[j])
	    +2.*pow(gamma1qg[ni],2)*(fn1[g]*fn2[g])
	    +gamma1qg[ni]*gamma1gg[ni]*(fn1[i]*fn2[g]+fn1[g]*fn2[j]);

	  for (int k = 0; k < 2*MAXNF+1; k++)
	    if (k != g)
	      tgaga += gamma1qg[ni]*gamma1gq[ni]*(fn1[i]*fn2[k]+fn1[k]*fn2[j]);

	  s22 =
	    0.5*(beta0N*A1qN*(LR-LQ)-A2qN)*tdelta
	    -0.5*A1qN*(H1q*tdelta+th1st+th1stQ+(LF-LQ)*th1stF)
	    -0.5*(B1qN+A1qN*LQ-beta0N)*s11
	    +0.5*(B1qN+A1qN*LQ)*th1stF
	    +0.5*tgaga;

	  //	  //add this piece to tgaga
	  //	  for (int k = 0; k < 2*MAXNF+1; k++)
	  //	    if (k != g)
	  //	      s22 += 0.5*gamma1qg[ni]*gamma1gq[ni]*(fn1[i]*fn2[k]+fn1[k]*fn2[j]);

	  s21 =
	    -beta0N*(LR-LQ)*s11
	    -(B1qN+A1qN*LQ)*(H1q*tdelta+th1stQ+th1st+(LF-LQ)*th1stF)
	    -(LF-LQ)*tgaga
	    -(B2qN+A2qN*LQ)*tdelta
	    +beta0N*th1st
	    +(B1qN+A1qN*LQ/2.)*LQ*th1stF
	    - (C1QG[ni]*(gamma1qq[ni]+gamma1gg[ni])+(H1q+2.*C1QQ[ni])*gamma1qg[ni])*(fn1[i]*fn2[g]+fn1[g]*fn2[j])  //tcga
	    -(2.*(gamma2qqV[ni]+gamma2qqS[ni]))*tdelta //tgamma2
	    -(gamma2qg[ni])*(fn1[i]*fn2[g]+fn1[g]*fn2[j])
	    -(gamma2qqbV[ni]+gamma2qqbS[ni])*(fn1[i]*fn2[jm]+fn1[im]*fn2[j]) //tgamma2
	    -2.*C1QG[ni]*gamma1qg[ni]*(fn1[g]*fn2[g]) //tgamma2
	    -2.*H1q*gamma1qq[ni]*tdelta
	    -4.*C1QQ[ni]*gamma1qq[ni]*tdelta;

	  for (int k = 0; k < 2*MAXNF+1; k++)
	    {
	      if (k != g && k != j && k != jm)
		s21 -= gamma2qqS[ni]*(fn1[i]*fn2[k]);
	      if (k != g && k != i && k != im)
		s21 -= gamma2qqS[ni]*(fn1[k]*fn2[j]);
	    }


	  for (int k = 0; k < 2*MAXNF+1; k++)
	    if (k != g)
	      s21 -= C1QG[ni]*gamma1gq[ni]*(fn1[i]*fn2[k]+fn1[k]*fn2[j]);
	  
      	  sig21[sp] += real(s21 * cexp * mellinint::wn[n]);
	  sig22[sp] += real(s22 * cexp * mellinint::wn[n]);
	  sig24[sp] += real(s24 * cexp * mellinint::wn[n]);
	  sig23[sp] += real(s23 * cexp * mellinint::wn[n]);
	}

      
    }
  
  double xmsq = 0.;
  for (int sp = 0; sp < mesq::totpch; sp++)
    {
      //as/pi factor
      double sig1 = (sig12[sp]*qtint::LL2_mesqij[sp]+sig11[sp]*qtint::LL1_mesqij[sp])*asopi/2.;
      xmsq += -sig1;
	  
      if (opts.order == 1) continue;
      //(as/pi)^2 factor
      double sig2 = (sig24[sp]*qtint::LL4_mesqij[sp]+sig23[sp]*qtint::LL3_mesqij[sp]+sig22[sp]*qtint::LL2_mesqij[sp]+sig21[sp]*qtint::LL1_mesqij[sp])*pow(asopi/2.,2);

      // sum O(as) and O(as^2) contributions	
      xmsq += -sig2;
    }	  

  //Do not fully understand this factor
  xmsq = xmsq * 3./8. /2./M_PI;

  //m^2/s factor
  xmsq = xmsq * m2/pow(opts.sroot,2);

  //phiV integration  
  xmsq = xmsq * 2*M_PI;

  if (isnan_ofast(xmsq))
    cout << m << " " << xmsq << endl;
  
  //switching function is inside qtint, do not apply
  
  f[0] = xmsq;

  //cout << "xmsq " << xmsq << endl;
  if (opts.fmufac > 0)
    evolnative::free();
  
  pdfevol::free_fx();

  return;
}
