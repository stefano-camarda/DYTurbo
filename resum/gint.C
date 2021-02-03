#include "gint.h"
#include "scales.h"
#include "phasespace.h"
#include "resconst.h"
#include "resint.h"
#include "settings.h"
#include "alphas.h"
#include "blim.h"
#include "gaussrules.h"
#include "mellinint.h"
#include "mesq.h"
#include "anomalous.h"
#include "ccoeff.h"
#include <iostream>

complex <double> gint::logS;
complex <double> gint::logasl;
complex <double> gint::logasl_pdf;
complex <double> gint::logasl_expc;
complex <double> *gint::alogqq;
complex <double> *gint::alogqg;
complex <double> *gint::alogqqb;
complex <double> *gint::alogqqp;
complex <double> *gint::alogqqbp;

using namespace resconst;
using namespace alphas;

void gint::init()
{

}

void gint::allocate()
{
  //allocate memory
  alogqq   = new complex <double> [mellinint::mdim*2];
  alogqg   = new complex <double> [mellinint::mdim*2];
  alogqqb  = new complex <double> [mellinint::mdim*2];
  alogqqp  = new complex <double> [mellinint::mdim*2];
  alogqqbp = new complex <double> [mellinint::mdim*2];
}
void gint::reset()
{
  logS = 0;
  logasl = 0;
  logasl_pdf = 0;
  logasl_expc = 0;
  fill(alogqq  , alogqq+mellinint::mdim*2, 0);
  fill(alogqg  , alogqg+mellinint::mdim*2, 0);
  fill(alogqqb , alogqqb+mellinint::mdim*2, 0);
  fill(alogqqp , alogqqp+mellinint::mdim*2, 0);
  fill(alogqqbp, alogqqbp+mellinint::mdim*2, 0);
}
void gint::free()
{
  delete[] alogqq;
  delete[] alogqg;
  delete[] alogqqb;
  delete[] alogqqp;
  delete[] alogqqbp;
}


void gint::intg(complex <double> q, complex <double> jac, double blim)
{
  //Compute: dq^2/q^2 A(as(q^2)) * log(M^2/q^2) + Btilde (as(q^2)) (Eq.19 of https://arxiv.org/pdf/hep-ph/0508068.pdf)
  alphas::calc(q, opts.order+1);
  
  complex <double> res = 0.; 

  //local bstar
  complex <double> qq;
  if (opts.bprescription == 4)
    {
      complex <double> bstartilde = b0/q;
    
      //undo tilde
      double Q = scales::res;
      complex <double> bstar;
      if (opts.modlog)
	bstar = sqrt(bstartilde*bstartilde - pow(b0/Q,2));
      else
	bstar = bstartilde;  //normal sudakov
      
      //undo star
      complex <double> b;
      b = sqrt(1./(1./pow(bstar,2)-1./pow(blim,2)));

      //recompute tilde
      complex <double> btilde;
      if (opts.modlog)
	btilde = sqrt(pow(b,2) + pow(b0/Q,2)); //modified sudakov
      else
	btilde = b;  //normal sudakov

      qq = b0/btilde;
    }
  else
    qq = q;
  
  //complex <double> ll = log(pow(phasespace::m/qq,2)); // --> as in Eq.19 of https://arxiv.org/pdf/hep-ph/0508068.pdf
  complex <double> ll = log(pow(scales::res/qq,2));     // gives results consistent with the analytic Sudakov
  
  //  double LQR = log(pow(scales::res/scales::ren,2));

  double B1qbar = B1q + A1q*resint::LQ;
  double B2qbar = B2q + A2q*resint::LQ;
  double B3qbar = B3q + A3q*resint::LQ;
  double B4qbar = B4q + A4q*resint::LQ;
  
  //Truncate alphas at the exact order
  if (opts.order_sudak == 0) res = A1q*as1_1l * ll;
  if (opts.order_sudak == 1) res = (A1q*as1_2l + A2q*as2_2l) * ll + B1qbar * as1_1l;
  if (opts.order_sudak == 2) res = (A1q*as1_3l + A2q*as2_3l + A3q*as3_3l) * ll + B1qbar * as1_2l + B2qbar * as2_2l;
  if (opts.order_sudak == 3) res = (A1q*as1_4l + A2q*as2_4l + A3q*as3_4l + A4q*as4_4l) * ll + B1qbar*as1_3l + B2qbar*as2_3l + B3qbar*as3_3l;
  if (opts.order_sudak == 4) res = (A1q*as1_5l + A2q*as2_5l + A3q*as3_5l + A4q*as4_5l + A5q*as5_5l) * ll + B1qbar*as1_4l + B2qbar*as2_4l + B3qbar*as3_4l + B4qbar*as4_4l;

  //Truncate alphas approximately
  //if (opts.order_sudak == 0) res = A1q*asLO * ll;
  //if (opts.order_sudak == 1) res = (A1q*asNLO + A2q*pow(asLO,2)) * ll + B1qbar * asLO;
  //if (opts.order_sudak == 2) res = (A1q*asNNLO + A2q*pow(asNLO,2) + A3q*pow(asLO,3)) * ll + B1qbar * asNLO + B2qbar * pow(asLO,2);
  //if (opts.order_sudak == 3) res = (A1q*asNNNLO + A2q*pow(asNNLO,2) + A3q*pow(asNLO,3) + A4q*pow(asLO,4)) * ll + B1qbar * asNNLO + B2qbar * pow(asNLO,2) + B3qbar * pow(asLO,3);

  //Do not truncate alphas
  if (opts.asrgkt)
    {
      if (opts.order_sudak == 0) res = A1q*as * ll;
      if (opts.order_sudak == 1) res = (A1q*as + A2q*pow(as,2)) * ll + B1qbar * as;
      if (opts.order_sudak == 2) res = (A1q*as + A2q*pow(as,2) + A3q*pow(as,3)) * ll + B1qbar * as + B2qbar * pow(as,2);
      if (opts.order_sudak == 3) res = (A1q*as + A2q*pow(as,2) + A3q*pow(as,3) + A4q*pow(as,4)) * ll + B1qbar * as + B2qbar * pow(as,2) + B3qbar * pow(as,3);
      if (opts.order_sudak == 4) res = (A1q*as + A2q*pow(as,2) + A3q*pow(as,3) + A4q*pow(as,4) + A5q*pow(as,5)) * ll + B1qbar*as + B2qbar*pow(as,2) + B3qbar*pow(as,3) + B4qbar*pow(as,4);
    }

  res *= 2./qq*jac;
  //res *= 1./pow(qq,2)*jac;

  logS -= res;
}


//void gint::intbeta(double q, double jac)
void gint::intbeta(complex <double> q, complex <double> jac, double blim)
{
  //Compute: dq^2/q^2 beta(as(q^2)) (Eq.103,105 of https://arxiv.org/pdf/hep-ph/0508068.pdf)
  alphas::calc(q, opts.order+1);
  
  //local bstar
  complex <double> qq;
  if (opts.bprescription == 4)
    {
      complex <double> bstartilde = b0/q;
    
      //undo tilde
      double Q = scales::res;
      complex <double> bstar;
      if (opts.modlog)
	bstar = sqrt(bstartilde*bstartilde - pow(b0/Q,2));
      else
	bstar = bstartilde;  //normal sudakov
      
      //undo star
      complex <double> b;
      b = sqrt(1./(1./pow(bstar,2)-1./pow(blim,2)));

      //recompute tilde
      complex <double> btilde;
      if (opts.modlog)
	btilde = sqrt(pow(b,2) + pow(b0/Q,2)); //modified sudakov
      else
	btilde = b;  //normal sudakov

      qq = b0/btilde;
    }
  else
    qq = q;
  
  //Truncate alphas at the exact order
  complex <double> fac = 2./qq*jac;
  if (!opts.asrgkt)
    {
      if (opts.order == 1) logasl -= fac*(beta0*as1_1l);
      if (opts.order == 2) logasl -= fac*(beta0*as1_2l + beta1*as2_2l);
      if (opts.order == 3) logasl -= fac*(beta0*as1_3l + beta1*as2_3l + beta2*as3_3l);
      if (opts.order == 4) logasl -= fac*(beta0*as1_4l + beta1*as2_4l + beta2*as3_4l + beta3*as4_4l);
    }

  //do not truncate alphas
  if (opts.asrgkt)
    {
      if (opts.order == 1) logasl -= fac*(beta0*as);
      if (opts.order == 2) logasl -= fac*(beta0*as + beta1*pow(as,2));
      if (opts.order == 3) logasl -= fac*(beta0*as + beta1*pow(as,2) + beta2*pow(as,3));
      if (opts.order == 4) logasl -= fac*(beta0*as + beta1*pow(as,2) + beta2*pow(as,3) + beta3*pow(as,4));
    }
  
  if (opts.expc == 0)
    return;
}  

void gint::intexpc(complex <double> q, complex <double> jac, double blim)
{
  //There is a problem in this routine, when b is not on the real axis (minimal prescription),
  //because the positive and negative branches of alogqg, etc... are identical (they should be different when b is complex
  
  //Compute: dq^2/q^2 beta(as(q^2)) dln(C~)/dln(as(q^2)) (Eq.105 of https://arxiv.org/pdf/hep-ph/0508068.pdf)
  //beta(as) dln(C)/dln(as) = beta(as) * as * C'/C = (b0*as+b1*as^2+b2*as^3)*as*(C1+2*as*C2+3*as^2*C3)/(delta_qa+as*C1+as^2*C2+as^3*C3)
  // ~ (b0*as+b1*as^2+b2*as^3)*as*(C1+2*as*C2+3*as^2*C3)*(1-as*C1)
  // ~ +as^2*b0*C1
  //   +as^3*(b1*C1+b0*(2*C2-C1^2)
  //   +as^4*(b2*C1+b1*(2*C2-C1^2)+b0*(3*C3-3*C1*C2+C1^3))
  
  if (!opts.mellin1d)
    {
      cout << "Error, gint::intexpc is not implemented for mellin2d" << endl;
      exit(0);
    }

  alphas::calc(q, opts.order+1);
  
  //local bstar
  complex <double> qq;
  if (opts.bprescription == 4)
    {
      complex <double> bstartilde = b0/q;
    
      //undo tilde
      double Q = scales::res;
      complex <double> bstar;
      if (opts.modlog)
	bstar = sqrt(bstartilde*bstartilde - pow(b0/Q,2));
      else
	bstar = bstartilde;  //normal sudakov
      
      //undo star
      complex <double> b;
      b = sqrt(1./(1./pow(bstar,2)-1./pow(blim,2)));

      //recompute tilde
      complex <double> btilde;
      if (opts.modlog)
	btilde = sqrt(pow(b,2) + pow(b0/Q,2)); //modified sudakov
      else
	btilde = b;  //normal sudakov

      qq = b0/btilde;
    }
  else
    qq = q;
  
  //Truncate alphas at the exact order
  complex <double> fac = 2./qq*jac;
  
  for (int sign = mesq::positive; sign <= mesq::negative; sign++)
    for (int i = 0; i < mellinint::mdim; i++)
      {
	int idx = anomalous::index(i,sign);
	if (opts.order_expc == 2)
	  {
	    if (opts.expc == 1)
		alogqq[idx]   += fac*(as2_2l*beta0*ccoeff::C1qq_delta);
	    else if (opts.expc == 3)
	      {
		alogqq[idx]   += fac*(as2_2l*beta0*ccoeff::C1qq[idx]);
		alogqg[idx]   += fac*(as2_2l*beta0*ccoeff::C2qg[idx]/ccoeff::C1qg[idx]);
		alogqqb[idx]  += fac*(as2_2l*beta0*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]);
		alogqqp[idx]  += fac*(as2_2l*beta0*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]);
		alogqqbp[idx] += fac*(as2_2l*beta0*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]);
	      }
	    else if (opts.expc == 4)
	      alogqq[idx]   += fac*(as2_2l*beta0*ccoeff::C1qq[idx]);
	    else if (opts.expc == 5)
	      {
		//Do not expand the C function at the denominator in the OD and DOD channels
		//2*C2/C1*beta0*as^2/(1+C2/C1*as) + beta0*as/(1+C2/C1*as) - beta0*as = C2/C1*beta0*as^2/(1+C2/C1*as)
		alogqq[idx]   += fac*(as2_2l*beta0*ccoeff::C1qq[idx]);
		alogqg[idx]   += fac*(as2_2l*beta0*ccoeff::C2qg[idx]  /ccoeff::C1qg[idx]  /(1.+as1_1l*ccoeff::C2qg[idx]  /ccoeff::C1qg[idx]  ));
		alogqqb[idx]  += fac*(as2_2l*beta0*ccoeff::C3qqb[idx] /ccoeff::C2qqb[idx] /(1.+as1_1l*ccoeff::C3qqb[idx] /ccoeff::C2qqb[idx] ));
		alogqqp[idx]  += fac*(as2_2l*beta0*ccoeff::C3qqp[idx] /ccoeff::C2qqp[idx] /(1.+as1_1l*ccoeff::C3qqp[idx] /ccoeff::C2qqp[idx] ));
		alogqqbp[idx] += fac*(as2_2l*beta0*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]/(1.+as1_1l*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]));
	      }
	  }
	if (opts.order_expc == 3)
	  {
	    if (opts.expc == 1)
		alogqq[idx] += fac*(
				    +as2_3l*beta0*ccoeff::C1qq_delta
				    +as3_3l*(beta1*ccoeff::C1qq_delta+beta0*(2.*ccoeff::C2qq_delta-pow(ccoeff::C1qq_delta,2)))
				    );
	    else if (opts.expc == 3)
	      {
		alogqq[idx] += fac*(+ as2_3l*beta0*ccoeff::C1qq[idx]
				    + as3_3l*(0.
					      + beta1*ccoeff::C1qq[idx]
					      - beta0*pow(ccoeff::C1qq[idx],2)
					      + 2.*beta0*ccoeff::C2qq[idx]
					      )
				    );
		alogqg[idx] += fac*(+ as2_3l*beta0*ccoeff::C2qg[idx]/ccoeff::C1qg[idx]
				    + as3_3l*(
					      + beta1*ccoeff::C2qg[idx]/ccoeff::C1qg[idx]
					      - beta0*pow(ccoeff::C2qg[idx]/ccoeff::C1qg[idx],2)
					      + 2.*beta0*ccoeff::C3qg[idx]/ccoeff::C1qg[idx]
					      )
				    );
		
		alogqqb[idx]  += fac*(+ as2_3l*beta0*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]
				      + as3_3l*(
						+ beta1*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]
						- beta0*pow(ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx],2)
						//+ 2.*beta0*C4/C2)
						));
		
		alogqqp[idx]  += fac*(+ as2_3l*beta0*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]
				      + as3_3l*(
						+ beta1*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]
						- beta0*pow(ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx],2)
						//+ 2.*beta0*C4/C2)
						));
		
		alogqqbp[idx] += fac*(+ as2_3l*beta0*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]
				      + as3_3l*(
						+ beta1*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]
						- beta0*pow(ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx],2)
						//+ 2.*beta0*C4/C2)
						));
	      }
	    else if (opts.expc == 4)
	      {
		//At order n expand up to terms containing Cn (i.e. expand up to n for D, up to n-1 for OD, and up to n-2 for DOD)
		alogqq[idx] += fac*(+ as2_3l*beta0*ccoeff::C1qq[idx]
				    + as3_3l*(0.
					      + beta1*ccoeff::C1qq[idx]
					      - beta0*pow(ccoeff::C1qq[idx],2)
					      + 2.*beta0*ccoeff::C2qq[idx]
					      )
				    );
		alogqg[idx]   += fac*(as2_2l*beta0*ccoeff::C2qg[idx]  /ccoeff::C1qg[idx]  /(1.+as1_1l*ccoeff::C2qg[idx]  /ccoeff::C1qg[idx]  ));
	      }
	    else if (opts.expc == 5)
	      {
		//Do not expand the C function at the denominator in the OD and DOD channels
		/*
		//extended formulas:
		N2LL =
		+ 2*C2/C1*beta0*asLO^2 /(1 + C2/C1*asLO)
		+ beta0*asLO	     /(1 + C2/C1*asLO)
		- beta0*asLO
		;
	      
		N3LL =
		+ 3*C3/C1*beta0*asLO^3 /(1 + C3/C1*asLO^2)
		+ 2*C2/C1*beta1*asLO^3 /(1 + C2/C1*asLO)
		+ 2*C2/C1*beta0*asLO^2 /(1 + C2/C1*asLO)
		+ beta1*asLO^2	     /(1 + C2/C1*asLO)
		+ beta0*asLO	     /(1 + C2/C1*asLO + C3/C1*asLO^2)
		- beta0*asLO
		- beta1*asLO^2
		- N2LL
		;

		N3LLb =
		+ 2*C2/C1*beta0*(2.*dasNLO*asLO) /(1 + C2/C1*asLO)
		+ beta0*(2*dasNLO)               /(1 + C2/C1*asLO)
		- beta0*(2*dasNLO)
		;

		N3LLmur =
		+ 2*C2/C1*beta0*(2.*dmuR*asLO) /(1 + C2/C1*asLO)
		+ beta0*(2*dmuR)               /(1 + C2/C1*asLO)
		- beta0*(2*dmuR)
		;

		//simplified formulas:
		N2LL = C2/C1*beta0*asLO^2/(1+C2/C1*asLO);
		N3LL = 2*C3/C1*beta0*asLO^3/(1+C2/C1*asLO+C3/C1*asLO^2)
		+C2/C1*beta1*asLO^3 /(1+C2/C1*asLO);
		N3LLb = C2/C1*beta0*2*asLO*dasNLO/(1+C2/C1*asLO);
		N3LLmur = C2/C1*beta0*2*asLO*dmuR/(1+C2/C1*asLO);
		*/
		
		alogqq[idx] += fac*(+ as2_3l*beta0*ccoeff::C1qq[idx]
				    + as3_3l*(0.
					      + beta1*ccoeff::C1qq[idx]
					      - beta0*pow(ccoeff::C1qq[idx],2)
					      + 2.*beta0*ccoeff::C2qq[idx]
					      )
				    );
		
		alogqg[idx] += fac*(
				    + ccoeff::C2qg[idx]/ccoeff::C1qg[idx]*beta0*as2_3l    /(1. + ccoeff::C2qg[idx]/ccoeff::C1qg[idx]*as1_1l)
				    + 2.*ccoeff::C3qg[idx]/ccoeff::C1qg[idx]*beta0*as3_3l /(1. + ccoeff::C2qg[idx]/ccoeff::C1qg[idx]*as1_1l + ccoeff::C3qg[idx]/ccoeff::C1qg[idx]*as2_2l)
				    + ccoeff::C2qg[idx]/ccoeff::C1qg[idx]*beta1*as3_3l    /(1. + ccoeff::C2qg[idx]/ccoeff::C1qg[idx]*as1_1l)
				    );

		//complex <double> N2LL = fac * ccoeff::C2qg[idx]/ccoeff::C1qg[idx]*beta0*as2_2l /(1. + ccoeff::C2qg[idx]/ccoeff::C1qg[idx]*as1_1l);
		//complex <double> N3LL = fac*(
		//			     + 2.*ccoeff::C3qg[idx]/ccoeff::C1qg[idx]*beta0*as3_3l /(1. + ccoeff::C2qg[idx]/ccoeff::C1qg[idx]*as1_1l + ccoeff::C3qg[idx]/ccoeff::C1qg[idx]*as2_2l)
		//			     + ccoeff::C2qg[idx]/ccoeff::C1qg[idx]*beta1*as3_3l    /(1. + ccoeff::C2qg[idx]/ccoeff::C1qg[idx]*as1_1l)
		//			     );
		//complex <double> N3LLb = fac*(ccoeff::C2qg[idx]/ccoeff::C1qg[idx]*beta0*(as2_3l-as2_2l) /(1. + ccoeff::C2qg[idx]/ccoeff::C1qg[idx]*as1_1l));
		//
		//alogqg[idx] += N2LL;
		//alogqg[idx] += N3LL;
		//alogqg[idx] += N3LLb;

		alogqqb[idx] += fac*(
				     + ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]*beta0*as2_3l    /(1. + ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]*as1_1l)
				     //+ 2*C4/C2*beta0*as3/(1+C3/C2*as1+C4/C2*as2)
				     + ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]*beta1*as3_3l    /(1. + ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]*as1_1l)
				     );
		alogqqp[idx] += fac*(
				     + ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]*beta0*as2_3l    /(1. + ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]*as1_1l)
				     //+ 2*C4/C2*beta0*as3/(1+C3/C2*as1+C4/C2*as2)
				     + ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]*beta1*as3_3l    /(1. + ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]*as1_1l)
				     );
		alogqqbp[idx] += fac*(
				     + ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]*beta0*as2_3l    /(1. + ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]*as1_1l)
				     //+ 2*C4/C2*beta0*as3/(1+C3/C2*as1+C4/C2*as2)
				     + ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]*beta1*as3_3l    /(1. + ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]*as1_1l)
				     );
	      }
	  }
	if (opts.order_expc == 4)
	  {
	    if (opts.expc == 1)
		alogqq[idx] += fac*(
				    +as2_4l*beta0*ccoeff::C1qq_delta
				    +as3_4l*(beta1*ccoeff::C1qq_delta+beta0*(2.*ccoeff::C2qq_delta-pow(ccoeff::C1qq_delta,2)))
				    +as4_4l*(0.
					     +beta2*ccoeff::C1qq_delta
					     +beta1*(2.*ccoeff::C2qq_delta-pow(ccoeff::C1qq_delta,2))
					     +beta0*(0.
						     +3.*ccoeff::C3qq_delta
						     -3.*ccoeff::C1qq_delta*ccoeff::C2qq_delta
						     +pow(ccoeff::C1qq_delta,3)
						     )
					     )
				    );
	    else if (opts.expc == 3)
	      {
		alogqq[idx] += fac*(+ as2_4l*beta0*ccoeff::C1qq[idx]
				    + as3_4l*(beta1*ccoeff::C1qq[idx]+beta0*(2.*ccoeff::C2qq[idx]-pow(ccoeff::C1qq[idx],2)))
				    + as4_4l*(
					      +beta2*ccoeff::C1qq[idx]
					      +beta1*(2.*ccoeff::C2qq[idx]-pow(ccoeff::C1qq[idx],2))
					      +beta0*(3.*ccoeff::C3qq[idx]-3.*ccoeff::C1qq[idx]*ccoeff::C2qq[idx]+pow(ccoeff::C1qq[idx],3))
					      )
				    );
		//...
		
	      }
	    else if (opts.expc = 4)
	      {
		alogqq[idx] += fac*(+ as2_4l*beta0*ccoeff::C1qq[idx]
				    + as3_4l*(beta1*ccoeff::C1qq[idx]+beta0*(2.*ccoeff::C2qq[idx]-pow(ccoeff::C1qq[idx],2)))
				    + as4_4l*(
					      +beta2*ccoeff::C1qq[idx]
					      +beta1*(2.*ccoeff::C2qq[idx]-pow(ccoeff::C1qq[idx],2))
					      +beta0*(3.*ccoeff::C3qq[idx]-3.*ccoeff::C1qq[idx]*ccoeff::C2qq[idx]+pow(ccoeff::C1qq[idx],3))
					      )
				    );
		alogqg[idx] += fac*(
				    + ccoeff::C2qg[idx]/ccoeff::C1qg[idx]*beta0*as2_3l    /(1. + ccoeff::C2qg[idx]/ccoeff::C1qg[idx]*as1_1l)
				    + 2.*ccoeff::C3qg[idx]/ccoeff::C1qg[idx]*beta0*as3_3l /(1. + ccoeff::C2qg[idx]/ccoeff::C1qg[idx]*as1_1l + ccoeff::C3qg[idx]/ccoeff::C1qg[idx]*as2_2l)
				    + ccoeff::C2qg[idx]/ccoeff::C1qg[idx]*beta1*as3_3l    /(1. + ccoeff::C2qg[idx]/ccoeff::C1qg[idx]*as1_1l)
				    );
		alogqqb[idx]  += fac*(as2_2l*beta0*ccoeff::C3qqb[idx] /ccoeff::C2qqb[idx] /(1.+as1_1l*ccoeff::C3qqb[idx] /ccoeff::C2qqb[idx] ));
		alogqqp[idx]  += fac*(as2_2l*beta0*ccoeff::C3qqp[idx] /ccoeff::C2qqp[idx] /(1.+as1_1l*ccoeff::C3qqp[idx] /ccoeff::C2qqp[idx] ));
		alogqqbp[idx] += fac*(as2_2l*beta0*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]/(1.+as1_1l*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]));
	      }
	    else if (opts.expc = 5)
	      {
		alogqq[idx] += fac*(+ as2_4l*beta0*ccoeff::C1qq[idx]
				    + as3_4l*(beta1*ccoeff::C1qq[idx]+beta0*(2.*ccoeff::C2qq[idx]-pow(ccoeff::C1qq[idx],2)))
				    + as4_4l*(
					      +beta2*ccoeff::C1qq[idx]
					      +beta1*(2.*ccoeff::C2qq[idx]-pow(ccoeff::C1qq[idx],2))
					      +beta0*(
						      +3.*ccoeff::C3qq[idx]
						      -3.*ccoeff::C1qq[idx]*ccoeff::C2qq[idx]
						      +pow(ccoeff::C1qq[idx],3)
						      )
					      )
				    );
		alogqg[idx] += fac*(
				    + ccoeff::C2qg[idx]/ccoeff::C1qg[idx]*beta0*as2_4l    /(1. + ccoeff::C2qg[idx]/ccoeff::C1qg[idx]*as1_2l + ccoeff::C3qg[idx]/ccoeff::C1qg[idx]*as2_2l)
				    + 2.*ccoeff::C3qg[idx]/ccoeff::C1qg[idx]*beta0*as3_4l /(1. + ccoeff::C2qg[idx]/ccoeff::C1qg[idx]*as1_1l + ccoeff::C3qg[idx]/ccoeff::C1qg[idx]*as2_2l)
				    + ccoeff::C2qg[idx]/ccoeff::C1qg[idx]*beta1*as3_4l    /(1. + ccoeff::C2qg[idx]/ccoeff::C1qg[idx]*as1_1l)
				    + ccoeff::C2qg[idx]/ccoeff::C1qg[idx]*beta2*as4_4l    /(1. + ccoeff::C2qg[idx]/ccoeff::C1qg[idx]*as1_1l)
				    + 2.*ccoeff::C3qg[idx]/ccoeff::C1qg[idx]*beta1*as4_4l /(1. + ccoeff::C2qg[idx]/ccoeff::C1qg[idx]*as1_1l + ccoeff::C3qg[idx]/ccoeff::C1qg[idx]*as2_2l)
				    //+ 3*C4/C1*beta0*as4 /(1+C2/C1*as1+C3/C1*as2+C4/C1*as3)
				    );
		alogqqb[idx] += fac*(
				     + ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]*beta0*as2_4l    /(1. + ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]*as1_2l)// + C4/C2*as2)
				     //+ 2*C4/C2*beta0*as3/(1+C3/C2*as1+C4/C2*as2)
				     + ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]*beta1*as3_4l    /(1. + ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]*as1_1l)
				     + ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]*beta2*as4_4l    /(1. + ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]*as1_1l)
				     //+ 2.*C4/C2*beta1*as4 /(1+C3/C2*as1+C4/C2*as2)
				     //+ 3*C5/C2*beta0*as4 /(1+C3/C2*as1+C4/C2*as2+C5/C2*as3)
				     );
		alogqqp[idx] += fac*(
				     + ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]*beta0*as2_4l    /(1. + ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]*as1_2l)// + C4/C2*as2)
				     //+ 2*C4/C2*beta0*as3/(1+C3/C2*as1+C4/C2*as2)
				     + ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]*beta1*as3_4l    /(1. + ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]*as1_1l)
				     + ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]*beta2*as4_4l    /(1. + ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]*as1_1l)
				     //+ 2.*C4/C2*beta1*as4 /(1+C3/C2*as1+C4/C2*as2)
				     //+ 3*C5/C2*beta0*as4 /(1+C3/C2*as1+C4/C2*as2+C5/C2*as3)
				     );
		alogqqbp[idx] += fac*(
				     + ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]*beta0*as2_4l    /(1. + ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]*as1_2l)// + C4/C2*as2)
				     //+ 2*C4/C2*beta0*as3/(1+C3/C2*as1+C4/C2*as2)
				     + ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]*beta1*as3_4l    /(1. + ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]*as1_1l)
				     + ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]*beta2*as4_4l    /(1. + ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]*as1_1l)
				     //+ 2.*C4/C2*beta1*as4 /(1+C3/C2*as1+C4/C2*as2)
				     //+ 3*C5/C2*beta0*as4 /(1+C3/C2*as1+C4/C2*as2+C5/C2*as3)
				     );
		
	      }
	    
	  }
      }
}

complex <double> btoqmin(complex <double> b, double blim, bool isbstar)
{
  complex <double> bstar;
  if (opts.bprescription == 0 || opts.bprescription == 4 || isbstar)
    bstar = real(b)/sqrt(1.+pow(real(b)/blim,2));
  else
    bstar = b;

  double Q = scales::res;
  complex <double> bstartilde;
  if (opts.modlog)
    bstartilde = sqrt(pow(bstar,2) + pow(b0/Q,2)); //modified sudakov
  else
    bstartilde = bstar;  //normal sudakov

  /*
  complex <double> mubstar = resconst::b0/bstar;
  complex <double> mubstartilde;
  if (opts.modlog)
    mubstartilde = mubstar * Q / sqrt((pow(mubstar,2) + pow(Q,2)));
  //mubstartilde = mubstar / sqrt(1.+pow(mubstar/Q,2));
  else
    mubstartilde = mubstar;

  //also could write
  //double mubstartilde = resconst::b0/bstartilde;
  */
  
  return b0/bstartilde;
}

void gint::calc(complex <double> b)
{
  //cout << "gint: b = " << b << endl;
  //Compute: - int _[(b0/b)^2] ^ [Q^2] dq^2/q^2 A(as(q^2)) * log(M^2/q^2) + Btilde (as(q^2)) (Eq.19 of https://arxiv.org/pdf/hep-ph/0508068.pdf)

  //In the global bstar prescription it is:
  //- int _[(b0/b*(b))^2] ^ [Q^2] dq^2/q^2 A(as(q^2)) * log(M^2/q^2) + Btilde (as(q^2))
  
  //In the local bstar prescription it is:
  //- int _[(b0/b*(b))^2] ^ [Q^2] dq^2/q^2 A(as(q^2)) * log(M^2/q^2) + Btilde (as(q^2))

  reset();

  int rule = 100;

  // boundaries of integration (allow complex values for the minimal prescription)
  complex <double> qmax,qmin,cq,mq;
  if (opts.bprescription == 2)
    qmax = scales::res+complex <double>(0.,1.)*imag(qmin);
  else
    qmax = scales::res;

  if (opts.numsud)
    {
      qmin = btoqmin(b, blim::sudakov, opts.bstar_sudakov);
      cq = (qmax+qmin)/2.;
      mq = (qmax-qmin)/2.;
      for (int i = 0; i < rule; i++)
	{
	  complex <double> q = cq+mq*gr::xxx[rule-1][i];
	  complex <double> jac = mq * gr::www[rule-1][i];
	  intg(q, jac, blim::sudakov);
	}

      qmin = btoqmin(b, blim::pdf, opts.bstar_pdf);
      cq = (qmax+qmin)/2.;
      mq = (qmax-qmin)/2.;
      for (int i = 0; i < rule; i++)
	{
	  complex <double> q = cq+mq*gr::xxx[rule-1][i];
	  complex <double> jac = mq * gr::www[rule-1][i];
	  intbeta(q, jac, blim::pdf);
	}
      logasl_pdf += logasl;
      logasl = 0.;
    }
  
  if (opts.numexpc)
    {
      qmin = btoqmin(b, blim::expc, opts.bstar_expc);
      cq = (qmax+qmin)/2.;
      mq = (qmax-qmin)/2.;
      for (int i = 0; i < rule; i++)
	{
	  complex <double> q = cq+mq*gr::xxx[rule-1][i];
	  complex <double> jac = mq * gr::www[rule-1][i];
	  intbeta(q, jac, blim::expc);
	  intexpc(q, jac, blim::expc);
	}
      logasl_expc += logasl;
      logasl = 0.;
    }

  //now the integration among the vertical line
  if (opts.bprescription == 2)
    {
      qmin = qmax;
      qmax = scales::res;      

      if (opts.numsud)
	{
	  cq = (qmax+qmin)/2.;
	  mq = (qmax-qmin)/2.;
	  for (int i = 0; i < rule; i++)
	    {
	      complex <double> q = cq+mq*gr::xxx[rule-1][i];
	      complex <double> jac = mq * gr::www[rule-1][i];
	      intg(q, jac, blim::sudakov);
	    }
	  cq = (qmax+qmin)/2.;
	  mq = (qmax-qmin)/2.;
	  for (int i = 0; i < rule; i++)
	    {
	      complex <double> q = cq+mq*gr::xxx[rule-1][i];
	      complex <double> jac = mq * gr::www[rule-1][i];
	      intbeta(q, jac, blim::pdf);
	    }
	  logasl_pdf += logasl;
	  logasl = 0.;
	}
      if (opts.numexpc)
	{
	  cq = (qmax+qmin)/2.;
	  mq = (qmax-qmin)/2.;
	  for (int i = 0; i < rule; i++)
	    {
	      complex <double> q = cq+mq*gr::xxx[rule-1][i];
	      complex <double> jac = mq * gr::www[rule-1][i];
	      intbeta(q, jac, blim::expc);
	      intexpc(q, jac, blim::expc);
	    }
	  logasl_expc += logasl;
	  logasl = 0.;
	}
    }

  
  //cout << qmin << " " << logasl << " int(beta) " << endl;
  //cout << qmin << " " << alogqq[0] << " gint::alogqq " << endl;
  //cout << qmin << " " << alogqg[0] << " gint::alogqg " << endl;
  //cout << qmin << " " << alogqqb[0] << " gint::alogqqb " << endl;
  //cout << qmin << " " << alogqqp[0] << " gint::alogqqp " << endl;
  //cout << qmin << " " << alogqqbp[0] << " gint::alogqqbp " << endl;

  //Do it differently, with a change of variable in the integration:
  /*
  double mu0 = b0/blim::sudakov;
  Q = scales::res - mu0;
  double mub = b0/b;
  double mubstartilde = mub * Q / sqrt((pow(mubstar,2) + pow(Q,2)));
  
  qmax = Q;
  qmin = mub;
  */
}  
