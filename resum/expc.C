#include "expc.h"
#include "resconst.h"
#include "mellinint.h"
#include "scales.h"
#include "blim.h"
#include "pdfevol.h"
#include "anomalous.h"
#include "hcoefficients.h"
#include "mesq.h"
#include "resint.h"
#include "settings.h"
#include "sudakovff.h"
#include "ccoeff.h"
#include "alphas.h"
#include "gint.h"
#include "pdf.h"
#include "Li2.hpp"
#include "isnan.h"
#include "psi.h"
#include "constants.h"
#include <iostream>

using namespace polylogarithm;

complex <double> *expc::qqb;
complex <double> *expc::qg;
complex <double> *expc::qg_1;
complex <double> *expc::qg_2;
complex <double> *expc::qq;
complex <double> *expc::qq_1;
complex <double> *expc::qq_2;
complex <double> *expc::qqp;
complex <double> *expc::qqp_1;
complex <double> *expc::qqp_2;
complex <double> *expc::qqbp;
complex <double> *expc::qqbp_1;
complex <double> *expc::qqbp_2;
complex <double> *expc::gg;
complex <double> *expc::qbg;
complex <double> *expc::qbg_1;
complex <double> *expc::qbg_2;
complex <double> *expc::qpg;
complex <double> *expc::qpg_1;
complex <double> *expc::qpg_2;
complex <double> *expc::qbpg;
complex <double> *expc::qbpg_1;
complex <double> *expc::qbpg_2;

//complex <double> expc::aexp;
//complex <double> expc::aexpB;

using namespace anomalous;
using namespace resconst;
using namespace resint;

void expc::allocate()
{
  //if (opts.order == 0)
  //return;

  //allocate memory
  if (opts.mellin1d)
    {
      qqb  = new complex <double> [mellinint::mdim*2];
      qg   = new complex <double> [mellinint::mdim*2];
      qq   = new complex <double> [mellinint::mdim*2];
      qqp  = new complex <double> [mellinint::mdim*2];
      qqbp = new complex <double> [mellinint::mdim*2];
      gg   = new complex <double> [mellinint::mdim*2];
      qbg  = new complex <double> [mellinint::mdim*2];
      qpg  = new complex <double> [mellinint::mdim*2];
      qbpg = new complex <double> [mellinint::mdim*2];
    }
  else
    {
      qqb    = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qg_1   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qg_2   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qq_1   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qq_2   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qqp_1  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qqp_2  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qqbp_1 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qqbp_2 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      gg     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbg_1  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbg_2  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qpg_1  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qpg_2  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbpg_1 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      qbpg_2 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
    }
}

void expc::free()
{
  //if (opts.order == 0)
  //return;

  if (opts.mellin1d)
    {
      delete[] qqb ;
      delete[] qg  ;
      delete[] qq  ;
      delete[] qqp ;
      delete[] qqbp;
      delete[] gg  ;
      delete[] qbg ;
      delete[] qpg ;
      delete[] qbpg;
    }
  else
    {
      delete[] qqb   ;
      delete[] qg_1  ;
      delete[] qg_2  ;
      delete[] qq_1  ;
      delete[] qq_2  ;
      delete[] qqp_1 ;
      delete[] qqp_2 ;
      delete[] qqbp_1;
      delete[] qqbp_2;
      delete[] gg    ;
      delete[] qbg_1 ;
      delete[] qbg_2 ;
      delete[] qpg_1 ;
      delete[] qpg_2 ;
      delete[] qbpg_1;
      delete[] qbpg_2;
    }
}

void expc::reset()
{
  //if (opts.order == 0)
  //return;

  if (opts.mellin1d)
    if (opts.sumlogs)
      {
	fill(qqb , qqb+mellinint::mdim*2, sudakov::S);
	fill(qg  , qg +mellinint::mdim*2, sudakov::S);
	fill(qq  , qq +mellinint::mdim*2, sudakov::S);
	fill(qqp , qqp +mellinint::mdim*2, sudakov::S);
	fill(qqbp, qqbp +mellinint::mdim*2, sudakov::S);
	fill(gg  , gg +mellinint::mdim*2, sudakov::S);
	fill(qbg , qbg +mellinint::mdim*2, sudakov::S);
	fill(qpg , qpg +mellinint::mdim*2, sudakov::S);
	fill(qbpg, qbpg +mellinint::mdim*2, sudakov::S);
      }
    else
      {
	fill(qqb , qqb+mellinint::mdim*2, 1.);
	fill(qg  , qg +mellinint::mdim*2, 1.);
	fill(qq  , qq +mellinint::mdim*2, 1.);
	fill(qqp , qqp +mellinint::mdim*2, 1.);
	fill(qqbp, qqbp +mellinint::mdim*2, 1.);
	fill(gg  , gg +mellinint::mdim*2, 1.);
	fill(qbg , qbg +mellinint::mdim*2, 1.);
	fill(qpg , qpg +mellinint::mdim*2, 1.);
	fill(qbpg, qbpg +mellinint::mdim*2, 1.);
      }
  else
    if (opts.sumlogs)
      {
	fill(qqb,    qqb +mellinint::mdim*mellinint::mdim*2, sudakov::S);
	fill(qg_1,   qg_1+mellinint::mdim*mellinint::mdim*2, sudakov::S);
	fill(qg_2,   qg_2+mellinint::mdim*mellinint::mdim*2, sudakov::S);
	fill(qq_1,   qq_1+mellinint::mdim*mellinint::mdim*2, sudakov::S);
	fill(qq_2,   qq_2+mellinint::mdim*mellinint::mdim*2, sudakov::S);
	fill(qqp_1,  qqp_1+mellinint::mdim*mellinint::mdim*2, sudakov::S);
	fill(qqp_2,  qqp_2+mellinint::mdim*mellinint::mdim*2, sudakov::S);
	fill(qqbp_1, qqbp_1+mellinint::mdim*mellinint::mdim*2, sudakov::S);
	fill(qqbp_2, qqbp_2+mellinint::mdim*mellinint::mdim*2, sudakov::S);
	fill(gg,     gg  +mellinint::mdim*mellinint::mdim*2, sudakov::S);

      }
    else
      {
	fill(qqb,    qqb +mellinint::mdim*mellinint::mdim*2, 1.);
	fill(qg_1,   qg_1+mellinint::mdim*mellinint::mdim*2, 1.);
	fill(qg_2,   qg_2+mellinint::mdim*mellinint::mdim*2, 1.);
	fill(qq_1,   qq_1+mellinint::mdim*mellinint::mdim*2, 1.);
	fill(qq_2,   qq_2+mellinint::mdim*mellinint::mdim*2, 1.);
	fill(qqp_1,  qqp_1+mellinint::mdim*mellinint::mdim*2, 1.);
	fill(qqp_2,  qqp_2+mellinint::mdim*mellinint::mdim*2, 1.);
	fill(qqbp_1, qqbp_1+mellinint::mdim*mellinint::mdim*2, 1.);
	fill(qqbp_2, qqbp_2+mellinint::mdim*mellinint::mdim*2, 1.);
	fill(gg,     gg  +mellinint::mdim*mellinint::mdim*2, 1.);
      }
}

//b-dependent C exponentiation
void expc::calc(complex <double> b)
{
  //Return 1 if b is close to the Landau pole
  double b_L;
  if (opts.modlog == 0)
    b_L = b0/scales::res * exp(1./(2.*aass*beta0));
  else if (opts.modlog == 1)
    b_L = b0/scales::res * sqrt(exp(1./(aass*beta0))-1.);
  else if (opts.modlog == 2)
    b_L = b0/scales::res * pow(sqrt(exp(opts.p/(aass*beta0))-1.),1./opts.p);

  if (fabs(b-b_L) < 1e-8)
    return;

  //  if (opts.expc == 0)
  //    return;

  //LL (no evolution)
  if (opts.order == 0)
    return;

  double blim = blim::expc;  //better use blim::sudakov?
  
//  the limit below implies xlambda<1/2 and then aa2<= 1      
//      blim=b0p*(1/q)*exp(1/(2*as*beta0)) 
//      blim=b0p*(1/q)*exp(1/(4*as*beta0))
// Set a limit to avoid very large values of b (= very small scales ~1/b)
//       blim=b0p*(1/q)*exp(1/(2*aass*beta0)) ! avoid Landau pole     
//     write(*,*) "blim",blim
//     without this additional blim some scale variations will fail (when mures > muren)
//  blim=0.5; // --> allow this to a separate blim in the settings, or set it using muren instead of mures
//      blim=1.1229190d0
//      blim=3d0
  
  //Set b according to bstar or other prescriptions
  complex <double> bstar;
  if (opts.bprescription == 0 || opts.bprescription == 4 || opts.bstar_expc)
    bstar = real(b)/sqrt(1.+pow(real(b)/blim,2));
  else
    bstar = b;

  double Q = scales::res;
  double aass2 = pow(aass,2);
  
  complex <double> blog;
  if (opts.modlog == 0)
    blog = log(pow(Q*bstar/b0,2));   //normal sudakov
  else if (opts.modlog == 1)
    blog = log(pow(Q*bstar/b0,2) + 1.); //modified sudakov
  else if (opts.modlog == 2)
    blog = 1./opts.p*log(pow(Q*bstar/b0,2*opts.p) + 1.); //modified sudakov with exponent p

  
  //if (opts.modlog)
  //  blog = log(pow(Q*bstar/b0,2) + 1.); //modified sudakov
  //else
  //  blog = log(pow(Q*bstar/b0,2));   //normal sudakov

  //lambda as defined in Eq. (25) of hep-ph/0508068.
  complex <double> xlambda = beta0*aass*blog;
  complex <double> log1xlambda = log(1.-xlambda);

  complex <double> logasl;

  //LL (no evolution)
  logasl = 0.;
  
  //NLL (LO evolution)
  if (opts.order >= 1)
    logasl += log1xlambda;

  //NNLL (NLO evolution)
  if (opts.order >= 2)
    logasl += aass* beta1/beta0*log1xlambda/(1.-xlambda);

  //NNNLL (NNLO evolution)
  if (opts.order >= 3)
    logasl += aass2* ((pow(beta1/beta0,2)-beta2/beta0) *xlambda/pow(1.-xlambda,2)
		      + pow(beta1/beta0,2)             *log1xlambda/pow(1.-xlambda,2)
		      - pow(beta1/beta0,2)             *pow(log1xlambda,2)/(2.*pow(1.-xlambda,2)));

  double LQR = real(logq2mur2-2.*loga);
  //QCD coupling scale dependence
  if (opts.order >= 2)
    logasl += aass*LQR*beta0*xlambda/(1.-xlambda);

  //from c3new.txt
  if (opts.order >= 3)
    logasl += aass2*(+LQR*beta1                   *(xlambda-log1xlambda)/pow(1.-xlambda,2)
		     +LQR*beta1                   *xlambda/(1.-xlambda)                      //missing piece
    		     +0.5*pow(LQR,2)*pow(beta0,2) *xlambda*(xlambda-2.)/pow(1.-xlambda,2));  //missing piece

  //from g4 LQR part by replacing B1qbar->beta0 and B2qbar->beta1 (and all An pieces set to zero)
  //if (opts.order >= 3)
  //  logasl += aass2*(LQR*beta1*(xlambda*(2.-xlambda)-log1xlambda)/pow(1.-xlambda,2)
  //		     +0.5*pow(LQR,2)*pow(beta0,2)*xlambda*(xlambda-2.)/pow(1.-xlambda,2));
  
  //cout << logasl << " logasl" << endl;
  //!!! gint::logasl is evaluated with blim::sudakov, not with blim:expc !!!
  if (opts.numexpc)
    logasl = gint::logasl;
  //cout << gint::logasl << " gint logasl" << endl;

  if (isnan_ofast(real(logasl)) || isnan_ofast(imag(logasl)))
    {
      cout << "Warning, logasl = " << logasl << ", b = "  << b << ", lambda = " << xlambda
	   << " log(1-xlambda) " << log1xlambda
	   << endl;
      logasl = 0.;
    }
  
  //For the aexp evolution operator use the same hard scale as for the PDF evolution, and the same blim
  //aexp should always reflect the PDF evolution, in evolmode 2 and 3 it should be log(ASF/ASI), plus eventually the LQR scale variation pieces
  //cout << "aexp " << logasl << " pdfs " << pdfevol::logasl << endl;
  //--> No aexp should be with Q and not muf

  //if (opts.mufevol)
  //complex <double> logasl = pdfevol::logasl;

  complex <double> aexp = exp(-logasl); //aexp is approximately alphas(Q^2)/alphas(mub^2)
  //aexp = exp(-logasl); //aexp is approximately alphas(Q^2)/alphas(mub^2)

  if (isnan_ofast(real(aexp)) || isnan_ofast(imag(aexp)))
    {
      cout << "Warning, aexp = " << aexp << ", b = "  << b << ", lambda = " << xlambda
	   << " log(1-xlambda) " << log1xlambda
	   << endl;
      aexp = 1.;
    }
  
  /*
  //Set b according to bstar or other prescriptions
  //For the C exponentiation always use Q at the resummation scale, and blim::sudakov, to match what is done in the Sudakov
  blim = blim::expc; //better use blim::sudakov;
  complex <double> bstar;
  if (opts.bprescription == 0 || opts.bstar_expc) //better use opts.bstar_sudakov)
    bstar = real(b)/sqrt(1.+pow(real(b)/blim,2));
  else
    bstar = b;

  double Q = scales::res;
  complex <double> blog;
  if (opts.modlog)
    blog = log(pow(Q*bstar/b0,2) + 1.); //modified sudakov
  else
    blog = log(pow(Q*bstar/b0,2));   //normal sudakov

  //lambda as defined in Eq. (25) of hep-ph/0508068.
  complex <double> xlambda = beta0*aass*blog;
  complex <double> log1xlambda = log(1.-xlambda);
  
  Q = scales::res;
  if (opts.modlog)
    blog = log(pow(Q*bstar/b0,2) + 1.); //modified sudakov
  else
    blog = log(pow(Q*bstar/b0,2));   //normal sudakov

  xlambda = beta0*aass*blog;
  log1xlambda = log(1.-xlambda);
  */

  //Multiply lamB, lamC, lamD by the corresponding power of alphas to avoid overflow/uunderflow problems in the exponential
  complex <double> lamB = aass*xlambda/(1.-xlambda);
  complex <double> aexpB = exp(lamB);

  complex <double> lamC = aass2*xlambda*(xlambda-2.)/pow(1.-xlambda,2);
  complex <double> aexpC = exp(lamC);

  complex <double> lamD = aass2*log1xlambda/pow(1.-xlambda,2);
  complex <double> aexpD = exp(lamD);


//  if (opts.bprescription == 0)
//    {
//      //double lambdaqcd         = muren/(exp(1./(2.*aass*beta0))); //--> corresponds to a divergence in alphas
//      double lambdaqcd     = mures/(exp(1./(2.*aass*beta0))); //--> correspond to a divergence in the Sudakov
//      blim = b0/lambdaqcd;
//      bstar = b/sqrt(1.+pow(b/blim,2));
//      if (opts.modlog)
//	blog = log(pow(Q*bstar/b0,2) + 1.); //modified sudakov
//      else
//	blog = log(pow(Q*bstar/b0,2));   //normal sudakov
//
//      xlambda = beta0*aass*blog;
//      log1xlambda = log(1.-xlambda);
//      
//      lamC = xlambda*(xlambda-2.)/pow(1.-xlambda,2);
//      aexpC = exp(lamC);
//    }
//
  
  double aassh = aass/2.;
  double aasshsq = pow(aassh,2);

  complex <double> aexp2 = pow(aexp,2);

  // NLL
  if (opts.order == 1)
    if (opts.mellin1d)
      if (opts.sumlogs)
	fill(qg, qg+mellinint::mdim*2, exp(-logasl+sudakov::logS));
      else
	fill(qg, qg+mellinint::mdim*2, aexp);
    else
      if (opts.sumlogs)
	{
	  fill(qg_1, qg_1+mellinint::mdim*mellinint::mdim*2, exp(-logasl+sudakov::logS));
	  fill(qg_2, qg_2+mellinint::mdim*mellinint::mdim*2, exp(-logasl+sudakov::logS));
	}
      else
	{
	  fill(qg_1, qg_1+mellinint::mdim*mellinint::mdim*2, aexp);
	  fill(qg_2, qg_2+mellinint::mdim*mellinint::mdim*2, aexp);
	}

  //gint::calc(b);
  //cout << endl;
  //cout << " gint::logasl "   << gint::logasl      << " expc " << logasl      << endl;
  // NNLL
  //complex <double> c1delta = pow(aexpb,aassh*(-2.*C1qqn));
  if (opts.order == 2)
    {
      if (opts.mellin1d)
	{
	  //Use pointers to allocate memory on the heap
	  complex <double> *alogqq   = new complex <double> [mellinint::mdim*2];
	  complex <double> *alogqg   = new complex <double> [mellinint::mdim*2];
	  complex <double> *alogqqb  = new complex <double> [mellinint::mdim*2];
	  complex <double> *alogqqp  = new complex <double> [mellinint::mdim*2];
	  complex <double> *alogqqbp = new complex <double> [mellinint::mdim*2];
	  complex <double> *aexpqq   = new complex <double> [mellinint::mdim*2];
	  complex <double> *aexpqg   = new complex <double> [mellinint::mdim*2];
	  complex <double> *aexpqqb  = new complex <double> [mellinint::mdim*2];
	  complex <double> *aexpqqp  = new complex <double> [mellinint::mdim*2];
	  complex <double> *aexpqqbp = new complex <double> [mellinint::mdim*2];
	  //complex <double> alogqq[mellinint::mdim*2];
	  //complex <double> alogqg[mellinint::mdim*2];
	  //complex <double> alogqqb[mellinint::mdim*2];
	  //complex <double> alogqqp[mellinint::mdim*2];
	  //complex <double> alogqqbp[mellinint::mdim*2];
	  //complex <double> aexpqq[mellinint::mdim*2];
	  //complex <double> aexpqg[mellinint::mdim*2];
	  //complex <double> aexpqqb[mellinint::mdim*2];
	  //complex <double> aexpqqp[mellinint::mdim*2];
	  //complex <double> aexpqqbp[mellinint::mdim*2];

	  //No exponentiation of the C coefficients
	  if (opts.expc == 0)
	    if (opts.sumlogs)
	      {
		fill(alogqq, alogqq+mellinint::mdim*2, 0.);
		fill(alogqg, alogqq+mellinint::mdim*2, 0.);
		fill(alogqqb, alogqqb+mellinint::mdim*2, 0.);
		fill(alogqqp, alogqqp+mellinint::mdim*2, 0.);
		fill(alogqqbp, alogqqbp+mellinint::mdim*2, 0.);
	      }
	    else
	      {
		fill(aexpqq, aexpqq+mellinint::mdim*2, 1.);
		fill(aexpqg, aexpqg+mellinint::mdim*2, 1.);
		fill(aexpqqb, aexpqqb+mellinint::mdim*2, 1.);
		fill(aexpqqp, aexpqqp+mellinint::mdim*2, 1.);
		fill(aexpqqbp, aexpqqbp+mellinint::mdim*2, 1.);
	      }
	  //Only delta terms
	  else if (opts.expc == 1)
	    if (opts.sumlogs)
	      {
		fill(alogqg, alogqq+mellinint::mdim*2, 0.);
		fill(alogqqb, alogqqb+mellinint::mdim*2, 0.);
		fill(alogqqp, alogqqp+mellinint::mdim*2, 0.);
		fill(alogqqbp, alogqqbp+mellinint::mdim*2, 0.);
		for (int sign = mesq::positive; sign <= mesq::negative; sign++)
		  for (int i = 0; i < mellinint::mdim; i++)
		    {
		      int idx = anomalous::index(i,sign);
		      alogqq[idx] = lamB*ccoeff::C1qq_delta;
		    }
	      }
	    else
	      {
		fill(aexpqg, aexpqg+mellinint::mdim*2, 1.);
		fill(aexpqqb, aexpqqb+mellinint::mdim*2, 1.);
		fill(aexpqqp, aexpqqp+mellinint::mdim*2, 1.);
		fill(aexpqqbp, aexpqqbp+mellinint::mdim*2, 1.);
		for (int sign = mesq::positive; sign <= mesq::negative; sign++)
		  for (int i = 0; i < mellinint::mdim; i++)
		    {
		      int idx = anomalous::index(i,sign);
		      aexpqq[idx] = pow(aexpB,ccoeff::C1qq_delta);
		    }
	      }
	  //Taylor expansion
	  else if (opts.expc == 2)
	    for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	      for (int i = 0; i < mellinint::mdim; i++)
		{
		  int idx = anomalous::index(i,sign);
		  complex <double> qqtayl = 1.;
		  complex <double> qgtayl = 1.;
		  complex <double> qqbtayl = 1.;
		  complex <double> qqptayl = 1.;
		  complex <double> qqbptayl = 1.;
		  int fact = 1;
		  for (int t = 1; t <= opts.ntaylor; t++)
		    {
		      fact *= t;
		      qqtayl += 1./double(fact) *pow(lamB*(ccoeff::C1qq[idx]-ccoeff::C1qq_delta),t);
		      qgtayl += 1./double(fact) *pow(lamB*ccoeff::C2qg[idx]/ccoeff::C1qg[idx],t);
		      qqbtayl += 1./double(fact) *pow(lamB*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx],t);
		      qqptayl += 1./double(fact) *pow(lamB*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx],t);
		      qqbptayl += 1./double(fact) *pow(lamB*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx],t);
		    }
		  if (opts.sumlogs)
		    {
		      alogqq[idx] = lamB*ccoeff::C1qq_delta + log(qqtayl);
		      alogqg[idx] = log(qgtayl);
		      alogqqb[idx] = log(qqbtayl);
		      alogqqp[idx] = log(qqptayl);
		      alogqqbp[idx] = log(qqbptayl);
		    }
		  else
		    {
		      aexpqq[idx] = pow(aexpB,ccoeff::C1qq_delta)*qqtayl;
		      aexpqg[idx] = qgtayl;
		      aexpqqb[idx] = qqbtayl;
		      aexpqqp[idx] = qqptayl;
		      aexpqqbp[idx] = qqbptayl;
		    }
		}
	  //Standard formula
	  else if (opts.expc == 3)
	    for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	      for (int i = 0; i < mellinint::mdim; i++)
		{
		  int idx = anomalous::index(i,sign);
		  
		  if (opts.sumlogs)
		    {
		      //alogqq = lamB*1./2.*(C1QQ[idx]-2.*C1qqn);
		      //alogqg = lamB*1./2.*(C2qgM[idx]/C1QG[idx]-2.*C1qqn);
		      alogqq[idx] = lamB*ccoeff::C1qq[idx];
		      alogqg[idx] = lamB*ccoeff::C2qg[idx]/ccoeff::C1qg[idx];
		      alogqqb[idx] = lamB*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx];
		      alogqqp[idx] = lamB*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx];
		      alogqqbp[idx] = lamB*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx];
		    }
		  else
		    {
		      //aexpqq = pow(aexpB,1./2.*(C1QQ[idx]-2.*C1qqn));
		      //aexpqg = pow(aexpB,1./2.*(C2qgM[idx]/C1QG[idx]-2.*C1qqn));
		      aexpqq[idx] = pow(aexpB,ccoeff::C1qq[idx]);
		      aexpqg[idx] = pow(aexpB,ccoeff::C2qg[idx]/ccoeff::C1qg[idx]);
		      aexpqqb[idx] = pow(aexpB,ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]);
		      aexpqqp[idx] = pow(aexpB,ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]);
		      aexpqqbp[idx] = pow(aexpB,ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]);
		    }
		}
	  //cout << endl;
	  //cout << aexpqq[0]   << " expc::aexpqq " << endl;
	  //cout << aexpqg[0]   << " expc::aexpqg " << endl;
	  //cout << aexpqqb[0]  << " expc::aexpqqb " << endl;
	  //cout << aexpqqp[0]  << " expc::aexpqqb " << endl;
	  //cout << aexpqqbp[0] << " expc::aexpqqbp " << endl;
	  //cout << endl;
	  //At order n expand up to terms containing Cn (i.e. expand up to n for D, up to n-1 for OD, and up to n-2 for DOD)
	  else if (opts.expc == 4)
	    {
	      if (opts.sumlogs)
		{
		  fill(alogqg, alogqq+mellinint::mdim*2, 0.);
		  fill(alogqqb, alogqqb+mellinint::mdim*2, 0.);
		  fill(alogqqp, alogqqp+mellinint::mdim*2, 0.);
		  fill(alogqqbp, alogqqbp+mellinint::mdim*2, 0.);
		}
	      else
		{
		  fill(aexpqg, aexpqg+mellinint::mdim*2, 1.);
		  fill(aexpqqb, aexpqqb+mellinint::mdim*2, 1.);
		  fill(aexpqqp, aexpqqp+mellinint::mdim*2, 1.);
		  fill(aexpqqbp, aexpqqbp+mellinint::mdim*2, 1.);
		}
	      for (int sign = mesq::positive; sign <= mesq::negative; sign++)
		for (int i = 0; i < mellinint::mdim; i++)
		  {
		    int idx = anomalous::index(i,sign);
		    if (opts.sumlogs)
		      alogqq[idx] = lamB*ccoeff::C1qq[idx];
		    else
		      aexpqq[idx] = pow(aexpB,ccoeff::C1qq[idx]);
		  }
	    }
	  //Do not expand C at the denominator
	  else if (opts.expc == 5)
	    for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	      for (int i = 0; i < mellinint::mdim; i++)
		{
		  int idx = anomalous::index(i,sign);
		  if (opts.sumlogs)
		    {
		      alogqq[idx] = lamB*ccoeff::C1qq[idx];
		      alogqg[idx]   = log(1.-xlambda/(1.+aass*ccoeff::C2qg[idx]/ccoeff::C1qg[idx])) - log(1.-xlambda);
		      alogqqb[idx]  = log(1.-xlambda/(1.+aass*ccoeff::C3qqb[idx] /ccoeff::C2qqb[idx] )) - log(1.-xlambda);
		      alogqqp[idx]  = log(1.-xlambda/(1.+aass*ccoeff::C3qqp[idx] /ccoeff::C2qqp[idx] )) - log(1.-xlambda);
		      alogqqbp[idx] = log(1.-xlambda/(1.+aass*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx])) - log(1.-xlambda);
		      //alternative with as0 at the denominator
		      //alogqg[idx]   = lamB*ccoeff::C2qg[idx]  /ccoeff::C1qg[idx]  /(1.+aass*ccoeff::C2qg[idx]  /ccoeff::C1qg[idx]  );
		      //alogqqb[idx]  = lamB*ccoeff::C3qqb[idx] /ccoeff::C2qqb[idx] /(1.+aass*ccoeff::C3qqb[idx] /ccoeff::C2qqb[idx] );
		      //alogqqp[idx]  = lamB*ccoeff::C3qqp[idx] /ccoeff::C2qqp[idx] /(1.+aass*ccoeff::C3qqp[idx] /ccoeff::C2qqp[idx] );
		      //alogqqbp[idx] = lamB*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]/(1.+aass*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]);
		      
		    }
		  else
		    {
		      aexpqq[idx] = pow(aexpB,ccoeff::C1qq[idx]);
		      aexpqg[idx]   = (1.-xlambda/(1.+aass*ccoeff::C2qg[idx]/ccoeff::C1qg[idx]))/(1.-xlambda);
		      aexpqqb[idx]  = (1.-xlambda/(1.+aass*ccoeff::C3qqb[idx] /ccoeff::C2qqb[idx] ))/(1.-xlambda);
		      aexpqqp[idx]  = (1.-xlambda/(1.+aass*ccoeff::C3qqp[idx] /ccoeff::C2qqp[idx] ))/(1.-xlambda);
		      aexpqqbp[idx] = (1.-xlambda/(1.+aass*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]))/(1.-xlambda);

		    }
		}
	  //gint::calc(b);
	  //cout << endl;
	  //cout << " gint::logasl "   << gint::logasl      << " expc " << logasl      << endl;
	  //cout << " gint::alogqq "   << gint::alogqq[0]   << " expc " << alogqq[0]   << endl;  
	  //cout << " gint::alogqg "   << gint::alogqg[0]   << " expc " << alogqg[0]   << endl;  
	  //cout << " gint::alogqqb "  << gint::alogqqb[0]  << " expc " << alogqqb[0]  << endl; 
	  //cout << " gint::alogqqp "  << gint::alogqqp[0]  << " expc " << alogqqp[0]  << endl; 
	  //cout << " gint::alogqqbp " << gint::alogqqbp[0] << " expc " << alogqqbp[0] << endl;
	  
	  //cout << aexpqq[0]   << " gint::alogqq " << endl;
	  //cout << aexpqg[0]   << " gint::alogqg " << endl;
	  //cout << aexpqqb[0]  << " gint::alogqqb " << endl;
	  //cout << aexpqqp[0]  << " gint::alogqqb " << endl;
	  //cout << aexpqqbp[0] << " gint::alogqqbp " << endl;
	  //cout << endl;

	  //!!! gint::alog are evaluated with blim::sudakov, not with blim:expc !!!
	  if (opts.numexpc)
	    for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	      for (int i = 0; i < mellinint::mdim; i++)
		{
		  int idx = anomalous::index(i,sign);
		  
		  if (opts.sumlogs)
		    {
		      alogqq[idx]   = gint::alogqq[idx];
		      alogqg[idx]   = gint::alogqg[idx];
		      alogqqb[idx]  = gint::alogqqb[idx];
		      alogqqp[idx]  = gint::alogqqp[idx];
		      alogqqbp[idx] = gint::alogqqbp[idx];
		    }
		  else
		    {
		      aexpqq[idx]   = exp(gint::alogqq[idx]);  
		      aexpqg[idx]   = exp(gint::alogqg[idx]);  
		      aexpqqb[idx]  = exp(gint::alogqqb[idx]); 
		      aexpqqp[idx]  = exp(gint::alogqqp[idx]); 
		      aexpqqbp[idx] = exp(gint::alogqqbp[idx]);
		    }
		}
	  
	  for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	    for (int i = 0; i < mellinint::mdim; i++)
	      {
		int idx = anomalous::index(i,sign);
		int ii = hcoefficients::index(i,sign);
		if (opts.sumlogs)
		  {
		    qqb[ii] = exp(alogqq[idx]+alogqq[idx]+sudakov::logS);
		    qg[ii]  = exp(-logasl+alogqg[idx]+alogqq[idx]+sudakov::logS);
		    qq[ii]  = exp(alogqq[idx]+alogqqb[idx]-2.*logasl+sudakov::logS);
		    qqp[ii] = exp(alogqq[idx]+alogqqbp[idx]-2.*logasl+sudakov::logS);
		    qqbp[ii] = exp(alogqq[idx]+alogqqp[idx]-2.*logasl+sudakov::logS);
		    gg[ii]  = exp(-logasl+alogqg[idx]-logasl+alogqg[idx]+sudakov::logS);
		  }
		else
		  {
		    qqb[ii] = aexpqq[idx] * aexpqq[idx];
		    qg[ii]  = aexp * aexpqg[idx]*aexpqq[idx];
		    qq[ii]  = aexpqq[idx] * aexpqqb[idx] * aexp2; //*c1delta
		    qqp[ii] = aexpqq[idx] * aexpqqbp[idx] * aexp2; //*c1delta
		    qqbp[ii] = aexpqq[idx] * aexpqqp[idx] * aexp2; //*c1delta
		    gg[ii]  = aexp*aexpqg[idx] * aexp*aexpqg[idx];
		  }
	      }
	  delete[] alogqq;
	  delete[] alogqg;
	  delete[] alogqqb;
	  delete[] alogqqp;
	  delete[] alogqqbp;
	  delete[] aexpqq;
	  delete[] aexpqg;
	  delete[] aexpqqb;
	  delete[] aexpqqp;
	  delete[] aexpqqbp;
	}
      else //mellin2d
	{
	  complex <double> aexpqq_1[mellinint::mdim*2];
	  complex <double> aexpqg_1[mellinint::mdim*2];
	  complex <double> aexpqqb_1[mellinint::mdim*2];
	  complex <double> aexpqqp_1[mellinint::mdim*2];
	  complex <double> aexpqqbp_1[mellinint::mdim*2];
	  complex <double> aexpqq_2[mellinint::mdim*2];
	  complex <double> aexpqg_2[mellinint::mdim*2];
	  complex <double> aexpqqb_2[mellinint::mdim*2];
	  complex <double> aexpqqp_2[mellinint::mdim*2];
	  complex <double> aexpqqbp_2[mellinint::mdim*2];

	  //case opts.expc == 0
	  fill(aexpqq_1, aexpqq_1+mellinint::mdim*2, 1.);
	  fill(aexpqg_1, aexpqg_1+mellinint::mdim*2, 1.);
	  fill(aexpqqb_1, aexpqqb_1+mellinint::mdim*2, 1.);
	  fill(aexpqqp_1, aexpqqp_1+mellinint::mdim*2, 1.);
	  fill(aexpqqbp_1, aexpqqbp_1+mellinint::mdim*2, 1.);
	  fill(aexpqq_2, aexpqq_2+mellinint::mdim*2, 1.);
	  fill(aexpqg_2, aexpqg_2+mellinint::mdim*2, 1.);
	  fill(aexpqqb_2, aexpqqb_2+mellinint::mdim*2, 1.);
	  fill(aexpqqp_2, aexpqqp_2+mellinint::mdim*2, 1.);
	  fill(aexpqqbp_2, aexpqqbp_2+mellinint::mdim*2, 1.);

	  //Only delta terms
	  if (opts.expc == 1)
	    for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	      for (int i = 0; i < mellinint::mdim; i++)
		{
		  int idx = anomalous::index(i,sign);
		  aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_delta);
		  aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_delta);
		}
	  //Taylor expansion
	  else if (opts.expc == 2)
	    {
	      for (int sign = mesq::positive; sign <= mesq::negative; sign++)
		for (int i = 0; i < mellinint::mdim; i++)
		  {
		    int idx = anomalous::index(i,sign);
		    complex <double> qqtayl_1 = 1.;
		    complex <double> qgtayl_1 = 1.;
		    complex <double> qqbtayl_1 = 1.;
		    complex <double> qqptayl_1 = 1.;
		    complex <double> qqbptayl_1 = 1.;
		    complex <double> qqtayl_2 = 1.;
		    complex <double> qgtayl_2 = 1.;
		    complex <double> qqbtayl_2 = 1.;
		    complex <double> qqptayl_2 = 1.;
		    complex <double> qqbptayl_2 = 1.;
		    int fact = 1;
		    for (int t = 1; t <= opts.ntaylor; t++)
		      {
			fact *= t;
			qqtayl_1 += 1./double(fact) *pow(lamB*(ccoeff::C1qq_1[idx]-ccoeff::C1qq_delta),t);
			qgtayl_1 += 1./double(fact) *pow(lamB*ccoeff::C2qg_1[idx]/ccoeff::C1qg_1[idx],t);
			qqbtayl_1 += 1./double(fact) *pow(lamB*ccoeff::C3qqb_1[idx]/ccoeff::C2qqb_1[idx],t);
			qqptayl_1 += 1./double(fact) *pow(lamB*ccoeff::C3qqp_1[idx]/ccoeff::C2qqp_1[idx],t);
			qqbptayl_1 += 1./double(fact) *pow(lamB*ccoeff::C3qqbp_1[idx]/ccoeff::C2qqbp_1[idx],t);
			qqtayl_2 += 1./double(fact) *pow(lamB*(ccoeff::C1qq_2[idx]-ccoeff::C1qq_delta),t);
			qgtayl_2 += 1./double(fact) *pow(lamB*ccoeff::C2qg_2[idx]/ccoeff::C1qg_2[idx],t);
			qqbtayl_2 += 1./double(fact) *pow(lamB*ccoeff::C3qqb_2[idx]/ccoeff::C2qqb_2[idx],t);
			qqptayl_2 += 1./double(fact) *pow(lamB*ccoeff::C3qqp_2[idx]/ccoeff::C2qqp_2[idx],t);
			qqbptayl_2 += 1./double(fact) *pow(lamB*ccoeff::C3qqbp_2[idx]/ccoeff::C2qqbp_2[idx],t);
		      }
		    aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_delta)*qqtayl_1;
		    aexpqg_1[idx] = qgtayl_1;
		    aexpqqb_1[idx] = qqbtayl_1;
		    aexpqqp_1[idx] = qqptayl_1;
		    aexpqqbp_1[idx] = qqbptayl_1;
		    aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_delta)*qqtayl_2;
		    aexpqg_2[idx] = qgtayl_2;
		    aexpqqb_2[idx] = qqbtayl_2;
		    aexpqqp_2[idx] = qqptayl_2;
		    aexpqqbp_2[idx] = qqbptayl_2;
		  }
	    }
	  //Standard expression
	  else if (opts.expc == 3)
	    for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	      for (int i = 0; i < mellinint::mdim; i++)
		{
		  int idx = anomalous::index(i,sign);
		  aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_1[idx]);
		  aexpqg_1[idx] = pow(aexpB,ccoeff::C2qg_1[idx]/ccoeff::C1qg_1[idx]);
		  aexpqqb_1[idx] = pow(aexpB,ccoeff::C3qqb_1[idx]/ccoeff::C2qqb_1[idx]);
		  aexpqqp_1[idx] = pow(aexpB,ccoeff::C3qqp_1[idx]/ccoeff::C2qqp_1[idx]);
		  aexpqqbp_1[idx] = pow(aexpB,ccoeff::C3qqbp_1[idx]/ccoeff::C2qqbp_1[idx]);
		  aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_2[idx]);
		  aexpqg_2[idx] = pow(aexpB,ccoeff::C2qg_2[idx]/ccoeff::C1qg_2[idx]);
		  aexpqqb_2[idx] = pow(aexpB,ccoeff::C3qqb_2[idx]/ccoeff::C2qqb_2[idx]);
		  aexpqqp_2[idx] = pow(aexpB,ccoeff::C3qqp_2[idx]/ccoeff::C2qqp_2[idx]);
		  aexpqqbp_2[idx] = pow(aexpB,ccoeff::C3qqbp_2[idx]/ccoeff::C2qqbp_2[idx]);
		}
	  //Do not expand C at the denominator in the off diagonal and double off diagonal channels
	  else if (opts.expc == 5)
	    for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	      for (int i = 0; i < mellinint::mdim; i++)
		{
		  int idx = anomalous::index(i,sign);
		  aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_1[idx]);
		  aexpqg_1[idx]   = (1. - xlambda/(1.+aass*ccoeff::C2qg_1[idx]/ccoeff::C1qg_1[idx]))/(1.-xlambda);
		  aexpqqb_1[idx]  = (1. - xlambda/(1.+aass*ccoeff::C3qqb_1[idx] /ccoeff::C2qqb_1[idx] ))/(1.-xlambda);
		  aexpqqp_1[idx]  = (1. - xlambda/(1.+aass*ccoeff::C3qqp_1[idx] /ccoeff::C2qqp_1[idx] ))/(1.-xlambda);
		  aexpqqbp_1[idx] = (1. - xlambda/(1.+aass*ccoeff::C3qqbp_1[idx]/ccoeff::C2qqbp_1[idx]))/(1.-xlambda);
		  aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_2[idx]);
		  aexpqg_2[idx]   = (1. - xlambda/(1.+aass*ccoeff::C2qg_2[idx]/ccoeff::C1qg_2[idx]))/(1.-xlambda);
		  aexpqqb_2[idx]  = (1. - xlambda/(1.+aass*ccoeff::C3qqb_2[idx] /ccoeff::C2qqb_2[idx] ))/(1.-xlambda);
		  aexpqqp_2[idx]  = (1. - xlambda/(1.+aass*ccoeff::C3qqp_2[idx] /ccoeff::C2qqp_2[idx] ))/(1.-xlambda);
		  aexpqqbp_2[idx] = (1. - xlambda/(1.+aass*ccoeff::C3qqbp_2[idx]/ccoeff::C2qqbp_2[idx]))/(1.-xlambda);
		}

	  for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	    for (int i1 = 0; i1 < mellinint::mdim; i1++)
	      for (int i2 = 0; i2 < mellinint::mdim; i2++)
		{
		  int idx1 = anomalous::index(i1,mesq::positive);
		  int idx2 = anomalous::index(i2,sign);
		  int ii12 = hcoefficients::index(i1,i2,sign);
		  qqb[ii12]  = aexpqq_1[idx1] * aexpqq_2[idx2];
		
		  //qg_1[ii12] = aexp * aexpqg_1[idx1]*aexpqq_1[idx1]; //--> bug fix in dyres
		  //qg_2[ii12] = aexp * aexpqg_2[idx2]*aexpqq_2[idx2]; //--> bug fix in dyres
		  qg_1[ii12] = aexp * aexpqg_1[idx1]*aexpqq_2[idx2];
		  qg_2[ii12] = aexp * aexpqg_2[idx2]*aexpqq_1[idx1];
		
		  qq_1[ii12] = aexpqq_1[idx1] * aexpqqb_2[idx2] * aexp2; //*c1delta
		  qq_2[ii12] = aexpqq_2[idx2] * aexpqqb_1[idx1] * aexp2; //*c1delta
		  qqp_1[ii12] = aexpqq_1[idx1] * aexpqqbp_2[idx2] * aexp2; //*c1delta
		  qqp_2[ii12] = aexpqq_2[idx2] * aexpqqbp_1[idx1] * aexp2; //*c1delta
		  qqbp_1[ii12] = aexpqq_1[idx1] * aexpqqp_2[idx2] * aexp2; //*c1delta
		  qqbp_2[ii12] = aexpqq_2[idx2] * aexpqqp_1[idx1] * aexp2; //*c1delta
		  gg[ii12]   = aexp*aexpqg_1[idx1] * aexp*aexpqg_2[idx2];
		}

	}
      //cout << alogqq[0] << " expc::alogqq " << endl;
      //cout << alogqg[0] << " expc::alogqg " << endl;
      //cout << alogqqb[0] << " expc::alogqqb " << endl;
      //cout << endl;

      //cout << "cexp::qqb " << qqb[0] << endl;
      //cout << "cexp::qg " << qg[0]  << endl;
      //cout << "cexp::qq " << qq[0]  << endl;
      
    }
  

  // NNNLL
  //multiple parton case as in Eqs.(100-107) of https://arxiv.org/pdf/hep-ph/0508068.pdf
  //B3tilde additional pieces as in Eq. (50)
  //-2 beta1 C1 + 2 beta0 C1^2 - 4 beta0 C2

  //(B3qbar*(-2. + y)*y)/(2.*beta0*pow(1. - y,2))
  //1./beta0 * 1./2.  y*(y-2.)/pow(1.-y,2) * (-2 beta1 C1 + 2 beta0 C1^2 - 4 beta0 C2)
  //log(aexpC) * (C1^2 - beta1/beta0*C1 - 2*C2)

  //pieces from g4
  //(B3qbar*(-2. + y)*y)/(2.*beta0*pow(1. - y,2))
  //(B2qbar*beta1*((2. - y)*y + 2.*log1y))/(2.*pow(beta0,2)*pow(1. - y,2))
  //((rlogq2mur2-2.*rloga)*(3.*y*(-2.*B2qbar*pow(beta0,3)*(-2. + y)*(-1. + y)))/(6.*pow(beta0,3)*pow(-1. + y,3));

  //(1/2)*(1/beta0)*B3qbar       *y*(y-2)/(1-y)^2
  //-(1/2)*B2qbar*beta1/beta0^2  *y*(y-2)/(1-y)^2
  //B2qbar*beta1/beta0^2        *log1y  /(1-y)^2
  //-B2qbar((rlogq2mur2-2.*rloga)*y*(y-2)/(1-y)^2

  //(-beta1/beta0*C1+C1^2-2*C2)       *y*(y-2)/(1-y)^2
  //C1*beta1/beta0                      *y*(y-2)/(1-y)^2
  //-2*C1*beta1/beta0                    *log1y  /(1-y)^2
  //(2*beta0*C1)*((rlogq2mur2-2.*rloga) *y*(y-2)/(1-y)^2

  //Diagonal form factor
  //Gca_D =
  //  + (aS/Pi)*C1                           *(lam)/(1-lam)
  //  
  //  + (aS^2/Pi^2)*(
  //		   + (C1^2/2 - C2)           *lam*(lam-2)/(1-lam)^2
  //		   - beta1/beta0*C1          *Log[1-lam]/(1-lam)^2
  //		   + beta0*C1*Log[Q2/muR2]   *lam*(lam-2)/(1-lam)^2
  //		   );
	  
  //Off-diagonal form factor
  //Gca_OD =
  //  - Log[1 - lam]
  //
  //  + (aS/Pi)*(
  //	       + (C2/C1)                         *lam/(1-lam)
  //	       - (beta1/beta0                    *Log[1-lam])/(1-lam)
  //	       - (beta0*Log[Q2/muR2])            *lam/(1-lam)
  //	       )
  //
  //  + (aS^2/Pi^2)*(
  //		   - beta1^2/beta0^2                        *lam/(1-lam)^2
  //		   + beta2/beta0                            *lam/(1-lam)^2
  //		   - beta1^2/beta0^2                        *Log[1-lam]/(1-lam)^2
  //		   + 1/2*beta1^2/beta0^2                    *Log[1-lam]^2/(1-lam)^2
  //		   + (C2^2/C1^2/2 - C1*C3/C1^2)             *lam*(lam-2)/(1-lam)^2
  //		   - beta1/beta0*C2/C1                      *Log[1-lam]/(1-lam)^2
  //
  //		   - beta1                  *Log[Q2/muR2]   *lam/(1-lam)^2
  //		   + beta1                  *Log[Q2/muR2]   *Log[1-lam]/(1-lam)^2
  //		   + beta0*C2/C1            *Log[Q2/muR2]   *lam*(lam-2)/(1-lam)^2
  //		   );

  //Doubly off-diagonal form factor
  //Gca_DOD =
  //  - 2*Log[1 - lam]
  //  + (aS/Pi)* (
  //		+ (C3/C2)                      *lam/(1-lam)
  //		- (2*beta1/beta0               *Log[1-lam])/(1-lam)
  //		- (2*beta0*Log[Q2/muR2])       *lam/(1-lam)
  //		)
  //  
  //  + (aS^2/Pi^2)* (
  //		    - 2*beta1^2/beta0^2                    *lam/(1-lam)^2
  //		    + 2*beta2/beta0                        *lam/(1-lam)^2
  //		    - 2*beta1^2/beta0^2                    *Log[1-lam]/(1-lam)^2
  //		    + beta1^2/beta0^2                      *Log[1-lam]^2/(1-lam)^2
  //		    + (C3^2/C2^2/2 - C4/C2)                *lam*(lam-2)/(1-lam)^2
  //		    + beta1/beta0*C3/C2                    *Log[1-lam]/(1-lam)^2
  //		     
  //		    - 2*beta1      *Log[Q2/muR2]          *lam/(1-lam)^2
  //		    + 2*beta1      *Log[Q2/muR2]          *Log[1-lam]/(1-lam)^2
  //		    + beta0*C3/C2  *Log[Q2/muR2]          *lam*(lam-2)/(1-lam)^2
  //		    );

  if (opts.order == 3)
    {
      if (opts.mellin1d)
	{
	  //Use pointers to allocate memory on the heap
	  complex <double> *alogqq   = new complex <double> [mellinint::mdim*2];
	  complex <double> *alogqg   = new complex <double> [mellinint::mdim*2];
	  complex <double> *alogqqb  = new complex <double> [mellinint::mdim*2];
	  complex <double> *alogqqp  = new complex <double> [mellinint::mdim*2];
	  complex <double> *alogqqbp = new complex <double> [mellinint::mdim*2];
	  complex <double> *aexpqq   = new complex <double> [mellinint::mdim*2];
	  complex <double> *aexpqg   = new complex <double> [mellinint::mdim*2];
	  complex <double> *aexpqqb  = new complex <double> [mellinint::mdim*2];
	  complex <double> *aexpqqp  = new complex <double> [mellinint::mdim*2];
	  complex <double> *aexpqqbp = new complex <double> [mellinint::mdim*2];

	  if (opts.expc == 0)
	    if (opts.sumlogs)
	      {
		fill(alogqq, alogqq+mellinint::mdim*2, 0.);
		fill(alogqg, alogqq+mellinint::mdim*2, 0.);
		fill(alogqqb, alogqqb+mellinint::mdim*2, 0.);
		fill(alogqqp, alogqqp+mellinint::mdim*2, 0.);
		fill(alogqqbp, alogqqbp+mellinint::mdim*2, 0.);
	      }
	    else
	      {
		fill(aexpqq, aexpqq+mellinint::mdim*2, 1.);
		fill(aexpqg, aexpqg+mellinint::mdim*2, 1.);
		fill(aexpqqb, aexpqqb+mellinint::mdim*2, 1.);
		fill(aexpqqp, aexpqqp+mellinint::mdim*2, 1.);
		fill(aexpqqbp, aexpqqbp+mellinint::mdim*2, 1.);
	      }
	  else if (opts.expc == 1)
	    {
	      if (opts.sumlogs)
		{
		  fill(alogqg, alogqq+mellinint::mdim*2, 0.);
		  fill(alogqqb, alogqqb+mellinint::mdim*2, 0.);
		  fill(alogqqp, alogqqp+mellinint::mdim*2, 0.);
		  fill(alogqqbp, alogqqbp+mellinint::mdim*2, 0.);
		  for (int sign = mesq::positive; sign <= mesq::negative; sign++)
		    for (int i = 0; i < mellinint::mdim; i++)
		      {
			int idx = anomalous::index(i,sign);
			alogqq[idx] = lamB*ccoeff::C1qq_delta
			  +lamC*(0.5*pow(ccoeff::C1qq_delta,2) - ccoeff::C2qq_delta)
			  -lamD*beta1/beta0*ccoeff::C1qq_delta
			  +lamC*beta0*ccoeff::C1qq_delta*(resint::rlogq2mur2-2.*resint::rloga)
			  ;
		      }
		}
	      else
		{
		  fill(aexpqg, aexpqg+mellinint::mdim*2, 1.);
		  fill(aexpqqb, aexpqqb+mellinint::mdim*2, 1.);
		  fill(aexpqqp, aexpqqp+mellinint::mdim*2, 1.);
		  fill(aexpqqbp, aexpqqbp+mellinint::mdim*2, 1.);
		  for (int sign = mesq::positive; sign <= mesq::negative; sign++)
		    for (int i = 0; i < mellinint::mdim; i++)
		      {
			int idx = anomalous::index(i,sign);
			aexpqq[idx] = pow(aexpB,ccoeff::C1qq_delta)
			  *pow(aexpC,(0.5*pow(ccoeff::C1qq_delta,2) - ccoeff::C2qq_delta))
			  *pow(aexpD,-beta1/beta0*ccoeff::C1qq_delta)
			  *pow(aexpC,beta0*ccoeff::C1qq_delta*(resint::rlogq2mur2-2.*resint::rloga))
			  ;
		      }
		}
	    }
	  else if (opts.expc == 2)
	    {
	      for (int sign = mesq::positive; sign <= mesq::negative; sign++)
		for (int i = 0; i < mellinint::mdim; i++)
		  {
		    int idx = anomalous::index(i,sign);
		    complex <double> qqtayl = 1.;
		    complex <double> qgtayl = 1.;
		    complex <double> qqbtayl = 1.;
		    complex <double> qqptayl = 1.;
		    complex <double> qqbptayl = 1.;
		    int fact = 1;
		    for (int t = 1; t <= opts.ntaylor; t++)
		      {
			fact *= t;
			qqtayl += 1./double(fact)* (pow(lamB*(ccoeff::C1qq[idx]-ccoeff::C1qq_delta),t)
						    +pow(lamC*(0.5*pow((ccoeff::C1qq[idx]-ccoeff::C1qq_delta),2) - (ccoeff::C2qq[idx]-ccoeff::C2qq_delta)),t)
						    +pow(-lamD*beta1/beta0*(ccoeff::C1qq[idx]-ccoeff::C1qq_delta),t)
						    //+pow(lamC*beta0*(ccoeff::C1qq[idx]-ccoeff::C1qq_delta)*(resint::rlogq2mur2-2.*resint::rloga),t)
						    );
						
			qgtayl += 1./double(fact)* (pow(lamB*ccoeff::C2qg[idx]/ccoeff::C1qg[idx],t)
						    +pow(lamC*(pow(ccoeff::C2qg[idx]/ccoeff::C1qg[idx],2)/2. - ccoeff::C3qg[idx]/ccoeff::C1qg[idx]),t)
						    +pow(-lamD*beta1/beta0*ccoeff::C2qg[idx]/ccoeff::C1qg[idx],t)
						    //+pow(lamC*(beta0*ccoeff::C2qg[idx]/ccoeff::C1qg[idx])*(resint::rlogq2mur2-2.*resint::rloga),t)
						    );
		    
			qqbtayl += 1./double(fact) *(pow(lamB*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx],t)
						     +pow(lamC*(pow(ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx],2)/2.),t)
						     -pow(lamD*beta1/beta0*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx],t)
						     //+pow(lamC*beta0*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]*(resint::rlogq2mur2-2.*resint::rloga),t)
						     );

			qqptayl += 1./double(fact) *(pow(lamB*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx],t)
						     +pow(lamC*(pow(ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx],2)/2.),t)
						     -pow(lamD*beta1/beta0*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx],t)
						     //+pow(lamC*beta0*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]*(resint::rlogq2mur2-2.*resint::rloga),t)
						     );

			qqbptayl += 1./double(fact) *(pow(lamB*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx],t)
						      +pow(lamC*(pow(ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx],2)/2.),t)
						      -pow(lamD*beta1/beta0*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx],t)
						      //+pow(lamC*beta0*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]*(resint::rlogq2mur2-2.*resint::rloga),t)
						      );
		      }
		    if (opts.sumlogs)
		      {
			alogqq[idx] = lamB*ccoeff::C1qq_delta
			  +lamC*(0.5*pow(ccoeff::C1qq_delta,2) - ccoeff::C2qq_delta)
			  -lamD*beta1/beta0*ccoeff::C1qq_delta
			  +lamC*beta0*ccoeff::C1qq_delta*(resint::rlogq2mur2-2.*resint::rloga)
			  +log(qqtayl);
			alogqg[idx] = log(qgtayl);
			alogqqb[idx] = log(qqbtayl);
			alogqqp[idx] = log(qqptayl); 
			alogqqbp[idx] = log(qqbptayl);
		      }
		    else
		      {
			aexpqq[idx] = pow(aexpB,ccoeff::C1qq_delta)
			  *pow(aexpC,(0.5*pow(ccoeff::C1qq_delta,2) - ccoeff::C2qq_delta))
			  *pow(aexpD,-beta1/beta0*ccoeff::C1qq_delta)
			  *pow(aexpC,beta0*ccoeff::C1qq_delta*(resint::rlogq2mur2-2.*resint::rloga))
			  *qqtayl;
			aexpqg[idx] = qgtayl;
			aexpqqb[idx] = qqbtayl;
			aexpqqp[idx] = qqptayl;
			aexpqqbp[idx] = qqbptayl;
		      }
		  }
	    }
	  else if (opts.expc == 3)
	    for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	      for (int i = 0; i < mellinint::mdim; i++)
		{
		  int idx = anomalous::index(i,sign);
		  if (opts.sumlogs)
		    {
//		      alogqq[idx] = lamB*1./2.*C1QQ[idx]
//			+lamC*1./4.*(0.5*pow(C1QQ[idx],2) - (C2NSqqM[idx] + C2SqqbM[idx]))
//			-lamD*1./2.*beta1/beta0*C1QQ[idx]
//			+lamC*1./2.*beta0*C1QQ[idx]*(resint::rlogq2mur2-2.*resint::rloga)
//			;

		      alogqq[idx] = lamB*ccoeff::C1qq[idx]
			+lamC*(0.5*pow(ccoeff::C1qq[idx],2) - ccoeff::C2qq[idx])
			-lamD*beta1/beta0*ccoeff::C1qq[idx]
			+lamC*beta0*ccoeff::C1qq[idx]*(LR-LQ)
			;

		      alogqg[idx] = lamB*(ccoeff::C2qg[idx]/ccoeff::C1qg[idx])
			//+lamC*1./4.*(pow(C2qgM[idx]/C1QG[idx],2)/2. - 8.*ccoeff::C3qg[idx]/C1QG[idx])
			//-lamD*1./2.*beta1/beta0*C2qgM[idx]/C1QG[idx]
			//+lamC*1./2.*(beta0*C2qgM[idx]/C1QG[idx])*(resint::rlogq2mur2-2.*resint::rloga)
			+lamC*(pow(ccoeff::C2qg[idx]/ccoeff::C1qg[idx],2)/2. - ccoeff::C3qg[idx]/ccoeff::C1qg[idx])
			-lamD*beta1/beta0*ccoeff::C2qg[idx]/ccoeff::C1qg[idx]
			+lamC*(beta0*ccoeff::C2qg[idx]/ccoeff::C1qg[idx])*(LR-LQ)
			;
		  
		      alogqqb[idx] = lamB*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]
			+lamC*(pow(ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx],2)/2.)// - C4/C2)
			-lamD*beta1/beta0*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]
			+lamC*beta0*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]*(LR-LQ)
			;

		      alogqqp[idx] = lamB*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]
			+lamC*(pow(ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx],2)/2.)// - C4/C2)
			-lamD*beta1/beta0*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]
			+lamC*beta0*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]*(LR-LQ)
			;
		  
		      alogqqbp[idx] = lamB*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]
			+lamC*(pow(ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx],2)/2.)// - C4/C2)
			-lamD*beta1/beta0*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]
			+lamC*beta0*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]*(LR-LQ)
			;
		  
		    }
		  else
		    {
		      aexpqq[idx] = pow(aexpB,ccoeff::C1qq[idx])
			*pow(aexpC,(0.5*pow(ccoeff::C1qq[idx],2) - ccoeff::C2qq[idx]))
			*pow(aexpD,-beta1/beta0*ccoeff::C1qq[idx])
			*pow(aexpC,beta0*ccoeff::C1qq[idx]*(LR-LQ))
			;
		      
		      aexpqg[idx] = pow(aexpB,ccoeff::C2qg[idx]/ccoeff::C1qg[idx])
			*pow(aexpC,pow(ccoeff::C2qg[idx]/ccoeff::C1qg[idx],2)/2. - ccoeff::C3qg[idx]/ccoeff::C1qg[idx])
			*pow(aexpD,-beta1/beta0*ccoeff::C2qg[idx]/ccoeff::C1qg[idx])
			*pow(aexpC,beta0*ccoeff::C2qg[idx]/ccoeff::C1qg[idx]*(LR-LQ))
			;
		
		      aexpqqb[idx] = pow(aexpB,ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx])
			*pow(aexpC,pow(ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx],2)/2.)// - C4/C2)
			*pow(aexpD,-beta1/beta0*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx])
			*pow(aexpC,beta0*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]*(LR-LQ))
			;
		
		      aexpqqp[idx] = pow(aexpB,ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx])
			*pow(aexpC,pow(ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx],2)/2.)// - C4/C2))
			*pow(aexpD,-beta1/beta0*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx])
			*pow(aexpC,beta0*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]*(LR-LQ))
			;
		
		      aexpqqbp[idx] = pow(aexpB,ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx])
			*pow(aexpC,pow(ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx],2)/2.)// - C4/C2))
			*pow(aexpD,-beta1/beta0*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx])
			*pow(aexpC,beta0*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]*(LR-LQ))
			;

		//if (sign == mesq::positive)
		//	cout << mellinint::Np[i] << "  " << aexpqg[idx]
		//	     << " C2/C1 " << pow(ccoeff::C2qg[idx]/ccoeff::C1qg[idx],2)
		//	     << " C3/C1 " << ccoeff::C3qg[idx]/ccoeff::C2qg[idx]
		//	     << " C3/C2 qqb" << ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]
		//	     << " C3/C2 qqp" << ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]
		//	     << " C3/C2 qqbp" << ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]
		//	     << endl;
		//else
		//	cout << mellinint::Nm[i] << "  " << aexpqg[idx]
		//	     << " C2/C1 " << pow(ccoeff::C2qg[idx]/ccoeff::C1qg[idx],2)
		//	     << " C3/C1 " << ccoeff::C3qg[idx]/ccoeff::C2qg[idx]
		//	     << " C3/C2 qqb" << ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]
		//	     << " C3/C2 qqp" << ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]
		//	     << " C3/C2 qqbp" << ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]
		//	     << endl;

		    }
		}

	  //cout << endl;
	  //cout << alogqq[0] << " expc::alogqq " << endl;
	  //cout << alogqg[0] << " expc::alogqg " << endl;
	  //cout << alogqqb[0] << " expc::alogqqb " << endl;
	  //cout << alogqqp[0] << " expc::alogqqb " << endl;
	  //cout << alogqqbp[0] << " expc::alogqqbp " << endl;
	  //cout << endl;
	  
	  //{
	  //  cout << endl;
	  //  int i;
	  //  i = 0;
	  //  //i = opts.mellinrule - 1;
	  //  cout << aexpqq[i]   << " expc::aexpqq " << endl;
	  //  cout << aexpqg[i]   << " expc::aexpqg " << endl;
	  //  cout << aexpqqb[i]  << " expc::aexpqqb " << endl;
	  //  cout << aexpqqp[i]  << " expc::aexpqqb " << endl;
	  //  cout << aexpqqbp[i] << " expc::aexpqqbp " << endl;
	  //  cout << endl;
	  //}
	  else if (opts.expc == 4)
	    {
	      if (opts.sumlogs)
		{
		  fill(alogqqb, alogqqb+mellinint::mdim*2, 0.);
		  fill(alogqqp, alogqqp+mellinint::mdim*2, 0.);
		  fill(alogqqbp, alogqqbp+mellinint::mdim*2, 0.);
		}
	      else
		{
		  fill(aexpqqb, aexpqqb+mellinint::mdim*2, 1.);
		  fill(aexpqqp, aexpqqp+mellinint::mdim*2, 1.);
		  fill(aexpqqbp, aexpqqbp+mellinint::mdim*2, 1.);
		}
	      for (int sign = mesq::positive; sign <= mesq::negative; sign++)
		for (int i = 0; i < mellinint::mdim; i++)
		  {
		    int idx = anomalous::index(i,sign);
		    if (opts.sumlogs)
		      {
			alogqq[idx] = lamB*ccoeff::C1qq[idx]
			  +lamC*(0.5*pow(ccoeff::C1qq[idx],2) - ccoeff::C2qq[idx])
			  -lamD*beta1/beta0*ccoeff::C1qq[idx]
			  +lamC*beta0*ccoeff::C1qq[idx]*(LR-LQ)
			  ;
			alogqg[idx]   = log(1.-xlambda/(1.+aass*ccoeff::C2qg[idx]/ccoeff::C1qg[idx])) - log(1.-xlambda);
		      }
		    else
		      {
			aexpqq[idx] = pow(aexpB,ccoeff::C1qq[idx])
			  *pow(aexpC,(0.5*pow(ccoeff::C1qq[idx],2) - ccoeff::C2qq[idx]))
			  *pow(aexpD,-beta1/beta0*ccoeff::C1qq[idx])
			  *pow(aexpC,beta0*ccoeff::C1qq[idx]*(LR-LQ))
			  ;
			aexpqg[idx]   = (1.-xlambda/(1.+aass*ccoeff::C2qg[idx]/ccoeff::C1qg[idx]))/(1.-xlambda);
		      }
		  }
	    }
	  //for (int i = 0; i < mellinint::mdim; i++)
	  //{
	  //  int i;
	  //  i = 0;
	  //  //i = opts.mellinrule - 1;
	  //  cout << i << endl;
	  //  cout << aexpqq[i]   << " gint::alogqq " << endl;
	  //  cout << aexpqg[i]   << "  " << exp(gint::alogqg[i])  << " gint::alogqg " << endl;
	  //  cout << aexpqqb[i]  << " gint::alogqqb " << endl;
	  //  cout << aexpqqp[i]  << " gint::alogqqb " << endl;
	  //  cout << aexpqqbp[i] << " gint::alogqqbp " << endl;
	  //  cout << endl;
	  //}
	  else if (opts.expc == 5)
	    {
	    for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	      for (int i = 0; i < mellinint::mdim; i++)
		{
		  int idx = anomalous::index(i,sign);
		  if (opts.sumlogs)
		    {
		      if (!opts.numexpc) gint::calc(b);		      
		      alogqq[idx] = lamB*ccoeff::C1qq[idx]
			+lamC*(0.5*pow(ccoeff::C1qq[idx],2) - ccoeff::C2qq[idx])
			-lamD*beta1/beta0*ccoeff::C1qq[idx]
			+lamC*beta0*ccoeff::C1qq[idx]*(LR-LQ)
			;

		      complex <double> C1,C2,C3,rootC,n2llalog;
		      C1 = ccoeff::C1qg[idx];
		      C2 = ccoeff::C2qg[idx];
		      C3 = ccoeff::C3qg[idx];
		      rootC = sqrt(pow(C2/C1,2)-4.*C3/C1);
		      complex <double> lam = xlambda;
		      //lamB = aass*lam/(1.-lam);
		      //lamC = aass2*lam*(lam-2.)/pow(1.-lam,2);
		      //lamD = aass2*log1lam/pow(1.-lam,2);

		      //original formulas
		      complex <double> oldN2LL,oldN3LL,oldN3LLb,oldN3LLmur;
		      oldN2LL    = lamB*(C2/C1);
		      oldN3LL    = (pow(C2/C1,2)/2. - C3/C1) * lamC - beta1/beta0 * C2/C1 * lamC/2.;
		      oldN3LLb   = - beta1/beta0 * C2/C1*lamD + beta1/beta0 * C2/C1 * lamC/2.;
		      oldN3LLmur = beta0*C2/C1*lamC*(LR-LQ);

		      //Formulas without expanding the denominator
		      complex <double> newN2LL,newN3LL,newN3LLb,newN3LLmur;

		      newN2LL = log(1.-lam/(1.+aass*C2/C1)) - log(1.-lam);
		      //newN2LL = log((1.-lam+aass*C2/C1)/(1.-lam)/(1.+aass*C2/C1)); //--> Check

		      /*
		      newN3LL = 0.
			//These pieces are problematic because of as0 at the denominator -> does not converge for bprescription = 2
			//integral of 3*C3/C1*beta0*as3_3l/(1+C2/C1*as0^2)
			//			-3.*C3/C1/(1.+aass2*C3/C1) * lamC/2.
			//integral of 2*C2/C1*beta1*as3_3l/(1+C2/C1*as0) --> not used
			//-2.*beta1/beta0*C2/C1/(1.+aass*C2/C1) * lamC/2.
			
			//integral of  2.*C2/C1*beta1*as3_3l/(1+C2/C1*as1_1l)
			-2.*beta1/beta0*(C1/C2*log((1.+aass*C2/C1/(1.-lam))/(1.+aass*C2/C1))-lamB)

			//This piece is problematic because of as0 at the denominator -> does not converge for bprescription = 2
			//integral of 	3*C3/C1*beta0*as3_3l/(1+C2/C1*as0^2)+2*C2/C1*beta1*as3_3l/(1+C2/C1*as0)
			-(aass*C2/C1*C3/C1*(3.+2.*aass*beta1/beta0)+2.*beta1/beta0*C2/C1+3.*C3/C1) / (1.+aass*C2/C1)/(1.+aass2*C3/C1) * lamC/2.

			//integral of beta0*as1_1l/(1+C2/C1*as1_1l+C3/C1*as2_2l)
			+C2/C1/sqrt(-pow(C2/C1,2)+4.*C3/C1)*(+atan((aass*C2/C1 + 2.*(1.-lam))/(aass*sqrt(-pow(C2/C1,2) + 4.*C3/C1))) 
							     -atan((aass*C2/C1 + 2.)         /(aass*sqrt(-pow(C2/C1,2) + 4.*C3/C1))))
			+0.5*log((1.+aass*C2/C1+aass2*C3/C1)/(pow(lam-1.,2)+aass*C2/C1*(1.-lam)+aass2*C3/C1))

			//integral of beta1*as2_2l/(1+C2/C1*as1_1l)
			+beta1/beta0*C1/C2*log((1.+aass*C2/C1/(1.-lam))/(1.+aass*C2/C1))
			
			//integral of +2*C2/C1*beta0*as2_2l/(1+C2/C1*as1_1l)
			+2.*log((1.+aass*C2/C1/(1.-lam))/(1.+aass*C2/C1))

			//sum of the two pieces above
			//+(2.+beta1/beta0*C1/C2)*log((1.-lam+aass*C2/C1))/(1.-lam)/(1.+aass*C2/C1)) //--> Check

			//integral of -beta0*as1_1l
			+log(1.-lam)

			//integral of -beta1*as2_2l
			-beta1/beta0*lamB
			;
		      */

		      //simplified alternative:
		      newN3LL =

			//integral of 2*C3/C1*beta0*as3_3l/(1+C3/C1*as2_2l) --> this piece is not correct, should be 2*C3/C1*beta0*as3_3l/(1+C2/C1*as1_1l+C3/C1*as2_2l)
			//-log((1.+aass2*C3/C1)/(pow(lam-1.,2)+aass2*C3/C1))
			//-2.*log(1.-lam)
			// --> equivalent to:
			//+log((pow(1.-lam,2)+aass2*C3/C1)/(pow(1.-lam,2)*(1.+aass2*C3/C1)))

			//integral of 2*C3/C1*beta0*as3_3l /(1+C2/C1*as1_1l+C3/C1*as2_2l)
			//C2/C1/sqrt(pow(C2/C1,2)-4.*C3/C1)
			//*log((aass*(C2/C1+sqrt(pow(C2/C1,2)-4.*C3/C1))+2.)/(aass*(C2/C1+sqrt(pow(C2/C1,2)-4.*C3/C1))+2.*(1.-lam))
			//     *(aass*(-C2/C1+sqrt(pow(C2/C1,2)-4.*C3/C1))-2.*(1.-lam))/(aass*(-C2/C1+sqrt(pow(C2/C1,2)-4.*C3/C1))-2.))
			//+log((pow(1.-lam,2)+aass*C2/C1*(1.-lam) + aass2*C3/C1)/pow(1.-lam,2)/(1.+aass*C2/C1+aass2*C3/C1))

			+log(pow((aass*(rootC+C2/C1)+2.)*(aass*(rootC-C2/C1)-2.*(1.-lam))/(aass*(rootC+C2/C1)+2.*(1.-lam))/(aass*(rootC-C2/C1)-2.),C2/C1/rootC)
			     *(pow(1.-lam,2)+aass*C2/C1*(1.-lam)+aass2*C3/C1)/pow(1.-lam,2)/(1.+aass*C2/C1+aass2*C3/C1))
			     //*(1.+aass*C2/C1/(1.-lam)+aass2*C3/C1/pow(1.-lam,2))/(1.+aass*C2/C1+aass2*C3/C1))

			
			//integral of  C2/C1*beta1*as3_3l/(1+C2/C1*as1_1l)
			//-beta1/beta0*(C1/C2*log((1.+aass*C2/C1/(1.-lam))/(1.+aass*C2/C1))-lamB)
			;
		      
		//newN3LLb = 1./(beta0*C2*(lam-1.))
		//	*beta1*(+2.*aass*C2*lam
		//		+log(1.-lam)*(2.*aass*C2-C1*(lam-1.)*(log(1.-lam)-2.*log(1.-(C1*(lam-1.))/(aass*C2))))
		//		+2.*C1*(lam-1.)*(-Li2(-((C1)/(aass*C2))) + Li2((C1*(lam-1.))/(aass*C2))));


		      newN3LLb =
			
			+2.*aass*beta1/beta0*(lam+log(1.-lam))/(lam-1.)
			-C1/C2*beta1/beta0 * log(1.-lam) *(log(1.-lam)+2.*log( aass*C2/C1/(1.-lam + aass*C2/C1)))
			+2.*C1/C2*beta1/beta0*(-Li2(-C1/C2/aass) + Li2(C1/C2*(lam-1.)/aass));
		      
		      
		      newN3LLmur =
			+2.*beta0*(C1/C2*log((1.+aass*C2/C1/(1.-lam))/(1.+aass*C2/C1))-lamB)*(LR-LQ);
			

		      //alternative with as0 at the denominator as much as possible --> Not good, because asLO/as diverges at lam -> 1
		      //newN2LL = lamB*C2/C1/(1.+aass*C2/C1);
		      //newN3LL = 0.5*
		      //	(-(aass*C2/C1*C3/C1*(3.+aass*beta1/beta0)+beta1/beta0*C2/C1+3.*C3/C1) / (1.+aass*C2/C1)/(1.+aass2*C3/C1) * lamC
		      //
		      //	 +2.*C2/C1/sqrt(-pow(C2/C1,2)+4.*C3/C1)*(+atan((aass*C2/C1 + 2.*(1.-lam))/(aass*sqrt(-pow(C2/C1,2) + 4.*C3/C1))) 
		      //						 -atan((aass*C2/C1 + 2.)         /(aass*sqrt(-pow(C2/C1,2) + 4.*C3/C1))))
		      //	 +2.*log(1.-lam)
		      //	 +4.*log((1.+aass*C2/C1/(1.-lam))/(1.+aass*C2/C1))
		      //	 +log((1.+aass*C2/C1+aass2*C3/C1)/(pow(lam-1.,2)+aass*C2/C1*(1.-lam)+aass2*C3/C1)));
		      //newN3LLb   = beta1/beta0*C2/C1*(lamC/2.-lamD)/(1.+aass*C2/C1);
		      //newN3LLmur = beta0*C2/C1*lamC*(LR-LQ)/(1.+aass*C2/C1);
		      

		      /*** test formulas ***
		      
		      complex <double> lam = 0.001;
		      
		      complex <double> lB = aass*lam/(1.-lam);
		      complex <double> lC = aass2*lam*(lam-2.)/pow(1.-lam,2);
		      complex <double> lD = aass2*log(1.-lam)/pow(1.-lam,2);

		      complex <double> C1 = 1.;
		      complex <double> C2 = complex <double> (0.1, 0.1);
		      complex <double> C3 = complex <double> (0.2, 0.2);
		      
		      complex <double> oldN3LLb = - beta1/beta0 * C2/C1*lD + beta1/beta0 * C2/C1 * lC/2.;
		      complex <double> newN3LLb = 1./(beta0*C2*(lam-1.))
			*beta1*(2.*aass*C2*lam+log(1.-lam)*(2.*aass*C2-C1*(lam-1.)
										   *(log(1.-lam)-2.*log(1.-(C1*(lam-1.))/(aass*C2))))
				+2.*C1*(lam-1.)*(- Li2(-((C1)/(aass*C2))) + Li2((C1*(lam-1.) )/(aass*C2))));
		      */
		      
		      /**********************************************************/
		      //alogqg[idx] = newN2LL+newN3LL+newN3LLb+newN3LLmur;

		      /*
		      if (sign == mesq::positive)
			cout << mellinint::Np[i];
		      else
			cout << mellinint::Nm[i];
		      cout << " old " << oldN3LL//+oldN3LL+oldN3LLb+oldN3LLmur
			   << " new " << newN3LL//+newN3LL+newN3LLb+newN3LLmur
			   << " gint " << gint::alogqg[i] << endl;
		      */
		      
//#define rootC(C1,C2,C3) sqrt(pow(C2/C1,2)-4.*C3/C1)
//#define logfinite(C1,C2) log((1.-lam+aass*C2/C1)/(1.-lam)/(1.+aass*C2/C1))
//#define n3llexpc(C1,C2,C3) logfinite(C1,C2)+C2/C1/rootC(C1,C2,C3)*log((aass*(C2/C1+rootC(C1,C2,C3))+2.)/(aass*(C2/C1+rootC(C1,C2,C3))+2.*(1.-lam))*(aass*(-C2/C1+rootC(C1,C2,C3))-2.*(1.-lam))/(aass*(-C2/C1+rootC(C1,C2,C3))-2.))+log((pow(1.-lam,2)+aass*C2/C1*(1.-lam)+aass2*C3/C1)/pow(1.-lam,2)/(1.+aass*C2/C1+aass2*C3/C1))-beta1/beta0*(C1/C2*logfinite(C1,C2)-lamB)+1./(beta0*C2*(lam-1.))*beta1*(+2.*aass*C2*lam+log(1.-lam)*(2.*aass*C2-C1*(lam-1.)*(log(1.-lam)-2.*log(1.-(C1*(lam-1.))/(aass*C2))))+2.*C1*(lam-1.)*(-Li2(-((C1)/(aass*C2))) + Li2((C1*(lam-1.))/(aass*C2))))+2.*beta0*(C1/C2*logfinite(C1,C2)-lamB)*(LR-LQ);
//		      alogqg[idx]   = n3llexpc(ccoeff::C1qg[idx],   ccoeff::C2qg[idx],   ccoeff::C3qg[idx]);
//		      alogqqb[idx]  = n3llexpc(ccoeff::C2qqb[idx],  ccoeff::C3qqb[idx],  0.);
//		      alogqqp[idx]  = n3llexpc(ccoeff::C2qqp[idx],  ccoeff::C3qqp[idx],  0.);
//		      alogqqbp[idx] = n3llexpc(ccoeff::C2qqbp[idx], ccoeff::C3qqbp[idx], 0.);

#define cdeflog(c1,c2,c3) C1 = c1; C2 = c2; C3 = c3; rootC = sqrt(pow(C2/C1,2)-4.*C3/C1); n2llalog = log((1.-lam+aass*C2/C1)/(1.-lam)/(1.+aass*C2/C1));
#define n3llalog() n2llalog\
			+C2/C1/rootC*log((aass*(C2/C1+rootC)+2.)/(aass*(C2/C1+rootC)+2.*(1.-lam))*(aass*(-C2/C1+rootC)-2.*(1.-lam))/(aass*(-C2/C1+rootC)-2.)) \
			+log((pow(1.-lam,2)+aass*C2/C1*(1.-lam)+aass2*C3/C1)/pow(1.-lam,2)/(1.+aass*C2/C1+aass2*C3/C1)) \
			-beta1/beta0*(C1/C2*n2llalog-lamB)		\
			+1./(beta0*C2*(lam-1.))*beta1*(+2.*aass*C2*lam+log(1.-lam)*(2.*aass*C2-C1*(lam-1.)*(log(1.-lam)-2.*log(1.-(C1*(lam-1.))/(aass*C2))))+2.*C1*(lam-1.)*(-Li2(-((C1)/(aass*C2))) + Li2((C1*(lam-1.))/(aass*C2)))) \
			+2.*beta0*(C1/C2*n2llalog-lamB)*(LR-LQ);
		      
		      cdeflog(ccoeff::C1qg[idx],   ccoeff::C2qg[idx],   ccoeff::C3qg[idx]); alogqg[idx]   = n3llalog();
		      //cdeflog(ccoeff::C1qg[idx],   ccoeff::C2qg[idx],   0.); alogqg[idx]   = n3llalog();
		      cdeflog(ccoeff::C2qqb[idx],  ccoeff::C3qqb[idx],  0.               ); alogqqb[idx]  = n3llalog();
		      cdeflog(ccoeff::C2qqp[idx],  ccoeff::C3qqp[idx],  0.               ); alogqqp[idx]  = n3llalog();
		      cdeflog(ccoeff::C2qqbp[idx], ccoeff::C3qqbp[idx], 0.               ); alogqqbp[idx] = n3llalog();
		      
//		      C1 = ccoeff::C1qg[idx]; C2 = ccoeff::C2qg[idx]; C3 = ccoeff::C3qg[idx]; rootC = sqrt(pow(C2/C1,2)-4.*C3/C1); n2llalog = log((1.-lam+aass*C2/C1)/(1.-lam)/(1.+aass*C2/C1));
//		      alogqg[idx]   = n2llalog+C2/C1/rootC*log((aass*(C2/C1+rootC)+2.)/(aass*(C2/C1+rootC)+2.*(1.-lam))*(aass*(-C2/C1+rootC)-2.*(1.-lam))/(aass*(-C2/C1+rootC)-2.))+log((pow(1.-lam,2)+aass*C2/C1*(1.-lam)+aass2*C3/C1)/pow(1.-lam,2)/(1.+aass*C2/C1+aass2*C3/C1))-beta1/beta0*(C1/C2*n2llalog-lamB)+1./(beta0*C2*(lam-1.))*beta1*(+2.*aass*C2*lam+log(1.-lam)*(2.*aass*C2-C1*(lam-1.)*(log(1.-lam)-2.*log(1.-(C1*(lam-1.))/(aass*C2))))+2.*C1*(lam-1.)*(-Li2(-((C1)/(aass*C2))) + Li2((C1*(lam-1.))/(aass*C2))))+2.*beta0*(C1/C2*n2llalog-lamB)*(LR-LQ);
//		      C1 = ccoeff::C2qqb[idx];  C2 = ccoeff::C3qqb[idx];  C3 = 0.; rootC = sqrt(pow(C2/C1,2)-4.*C3/C1); n2llalog = log((1.-lam+aass*C2/C1)/(1.-lam)/(1.+aass*C2/C1));
//		      alogqqb[idx]  = n2llalog+C2/C1/rootC*log((aass*(C2/C1+rootC)+2.)/(aass*(C2/C1+rootC)+2.*(1.-lam))*(aass*(-C2/C1+rootC)-2.*(1.-lam))/(aass*(-C2/C1+rootC)-2.))+log((pow(1.-lam,2)+aass*C2/C1*(1.-lam)+aass2*C3/C1)/pow(1.-lam,2)/(1.+aass*C2/C1+aass2*C3/C1))-beta1/beta0*(C1/C2*n2llalog-lamB)+1./(beta0*C2*(lam-1.))*beta1*(+2.*aass*C2*lam+log(1.-lam)*(2.*aass*C2-C1*(lam-1.)*(log(1.-lam)-2.*log(1.-(C1*(lam-1.))/(aass*C2))))+2.*C1*(lam-1.)*(-Li2(-((C1)/(aass*C2))) + Li2((C1*(lam-1.))/(aass*C2))))+2.*beta0*(C1/C2*n2llalog-lamB)*(LR-LQ);
//		      C1 = ccoeff::C2qqp[idx];  C2 = ccoeff::C3qqp[idx];  C3 = 0.; rootC = sqrt(pow(C2/C1,2)-4.*C3/C1); n2llalog = log((1.-lam+aass*C2/C1)/(1.-lam)/(1.+aass*C2/C1));
//		      alogqqp[idx]  = n2llalog+C2/C1/rootC*log((aass*(C2/C1+rootC)+2.)/(aass*(C2/C1+rootC)+2.*(1.-lam))*(aass*(-C2/C1+rootC)-2.*(1.-lam))/(aass*(-C2/C1+rootC)-2.))+log((pow(1.-lam,2)+aass*C2/C1*(1.-lam)+aass2*C3/C1)/pow(1.-lam,2)/(1.+aass*C2/C1+aass2*C3/C1))-beta1/beta0*(C1/C2*n2llalog-lamB)+1./(beta0*C2*(lam-1.))*beta1*(+2.*aass*C2*lam+log(1.-lam)*(2.*aass*C2-C1*(lam-1.)*(log(1.-lam)-2.*log(1.-(C1*(lam-1.))/(aass*C2))))+2.*C1*(lam-1.)*(-Li2(-((C1)/(aass*C2))) + Li2((C1*(lam-1.))/(aass*C2))))+2.*beta0*(C1/C2*n2llalog-lamB)*(LR-LQ);
//		      C1 = ccoeff::C2qqbp[idx]; C2 = ccoeff::C3qqbp[idx]; C3 = 0.; rootC = sqrt(pow(C2/C1,2)-4.*C3/C1); n2llalog = log((1.-lam+aass*C2/C1)/(1.-lam)/(1.+aass*C2/C1));
//		      alogqqbp[idx] = n2llalog+C2/C1/rootC*log((aass*(C2/C1+rootC)+2.)/(aass*(C2/C1+rootC)+2.*(1.-lam))*(aass*(-C2/C1+rootC)-2.*(1.-lam))/(aass*(-C2/C1+rootC)-2.))+log((pow(1.-lam,2)+aass*C2/C1*(1.-lam)+aass2*C3/C1)/pow(1.-lam,2)/(1.+aass*C2/C1+aass2*C3/C1))-beta1/beta0*(C1/C2*n2llalog-lamB)+1./(beta0*C2*(lam-1.))*beta1*(+2.*aass*C2*lam+log(1.-lam)*(2.*aass*C2-C1*(lam-1.)*(log(1.-lam)-2.*log(1.-(C1*(lam-1.))/(aass*C2))))+2.*C1*(lam-1.)*(-Li2(-((C1)/(aass*C2))) + Li2((C1*(lam-1.))/(aass*C2))))+2.*beta0*(C1/C2*n2llalog-lamB)*(LR-LQ);

		      //alogqqb[idx]  = log(1.-xlambda/(1.+aass*ccoeff::C3qqb[idx] /ccoeff::C2qqb[idx] )) - log(1.-xlambda);
		      //alogqqp[idx]  = log(1.-xlambda/(1.+aass*ccoeff::C3qqp[idx] /ccoeff::C2qqp[idx] )) - log(1.-xlambda);
		      //alogqqbp[idx] = log(1.-xlambda/(1.+aass*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx])) - log(1.-xlambda);
		    }
		  else
		    {
//		      if (opts.numexpc) gint::calc(b);
//		      complex <double> C1,C2,C3,rootC,n2llaexp;
//		      C1 = ccoeff::C1qg[idx];
//		      C2 = ccoeff::C2qg[idx];
//		      C3 = ccoeff::C3qg[idx];
//		      rootC = sqrt(pow(C2/C1,2)-4.*C3/C1);
//		      complex <double> lam = xlambda;
//
//		      //original formulas
//		      complex <double> oldN2LL,oldN3LL,oldN3LLb,oldN3LLmur;
//		      oldN2LL    = exp(lamB*(C2/C1));
//		      oldN3LL    = exp((pow(C2/C1,2)/2. - C3/C1) * lamC - beta1/beta0 * C2/C1 * lamC/2.);
//		      oldN3LLb   = exp(- beta1/beta0 * C2/C1*lamD + beta1/beta0 * C2/C1 * lamC/2.);
//		      oldN3LLmur = exp(beta0*C2/C1*lamC*(LR-LQ));
//		      
//		      //Formulas without expanding the denominator
//		      complex <double> newN2LL,newN3LL,newN3LLb,newN3LLmur;
//
//#define cdefexp(c1,c2,c3) C1 = c1; C2 = c2; C3 = c3; rootC = sqrt(pow(C2/C1,2)-4.*C3/C1); n2llaexp = (1.+aass*C2/C1/(1.-lam))/(1.+aass*C2/C1); //n2llaexp = (1.-lam+aass*C2/C1)/(1.-lam)/(1.+aass*C2/C1);
//		      cdefexp(ccoeff::C1qg[idx],   ccoeff::C2qg[idx],   ccoeff::C3qg[idx]);
//		      newN2LL = (1.+aass*C2/C1/(1.-lam))/(1.+aass*C2/C1);
//		      newN3LL = exp(log((aass*(rootC+C2/C1)+2.)*(aass*(rootC-C2/C1)-2.*(1.-lam))/(aass*(rootC+C2/C1)+2.*(1.-lam))/(aass*(rootC-C2/C1)-2.))*C2/C1/rootC) \
//			*(pow(1.-lam,2)+aass*C2/C1*(1.-lam)+aass2*C3/C1)/pow(1.-lam,2)/(1.+aass*C2/C1+aass2*C3/C1) \
//			* exp(-log(n2llaexp)*beta1/beta0*C1/C2) *pow(aexpB,beta1/beta0);
//		      newN3LLb = exp(+2.*aass*beta1/beta0*(lam+log(1.-lam))/(lam-1.)-C1/C2*beta1/beta0 * log(1.-lam) *(log(1.-lam)+2.*log( aass*C2/C1/(1.-lam + aass*C2/C1)))+2.*C1/C2*beta1/beta0*(-Li2(-C1/C2/aass)+Li2(C1/C2*(lam-1.)/aass)));
//		      newN3LLmur = exp(log(n2llaexp)*2.*beta0*(C1/C2)*(LR-LQ))	\
//			*pow(aexpB,-2.*beta0*(LR-LQ));
//			
//		      if (sign == mesq::positive)
//			cout << mellinint::Np[i];
//		      else
//			cout << mellinint::Nm[i];
//		      cout << " old " << oldN3LLb//+oldN3LL+oldN3LLb+oldN3LLmur
//			   << " new " << newN3LLb//+newN3LL+newN3LLb+newN3LLmur
//			   << " gint " << exp(gint::alogqg[i]) << endl;
		      
		      
		      aexpqq[idx] = pow(aexpB,ccoeff::C1qq[idx])
			*pow(aexpC,(0.5*pow(ccoeff::C1qq[idx],2) - ccoeff::C2qq[idx]))
			*pow(aexpD,-beta1/beta0*ccoeff::C1qq[idx])
			*pow(aexpC,beta0*ccoeff::C1qq[idx]*(LR-LQ))
			;

		      if (isnan_ofast(real(aexpqq[idx])) || isnan_ofast(imag(aexpqq[idx])))
			{
			  cout << "Warning, aexpqq[" << idx << "] = " << aexpqq[idx] << ", b = "  << b << ", lambda = " << xlambda
			       << " aexpB " << aexpB << " aexpC " << aexpC << " aexpD " << aexpD
			       << endl;
			}
		      
		      complex <double> C1,C2,C3,rootC,n2llaexp;
		      complex <double> lam = xlambda;

#define cdefexp(c1,c2,c3) C1 = c1; C2 = c2; C3 = c3; rootC = sqrt(pow(C2/C1,2)-4.*C3/C1); n2llaexp = (1.+aass*C2/C1/(1.-lam))/(1.+aass*C2/C1); //n2llaexp = (1.-lam+aass*C2/C1)/(1.-lam)/(1.+aass*C2/C1);
#define n3llaexp() n2llaexp\
			*exp(log((aass*(rootC+C2/C1)+2.)*(aass*(rootC-C2/C1)-2.*(1.-lam))/(aass*(rootC+C2/C1)+2.*(1.-lam))/(aass*(rootC-C2/C1)-2.))*C2/C1/rootC) \
			*(pow(1.-lam,2)+aass*C2/C1*(1.-lam)+aass2*C3/C1)/pow(1.-lam,2)/(1.+aass*C2/C1+aass2*C3/C1) \
			* exp(-log(n2llaexp)*beta1/beta0*C1/C2) *pow(aexpB,beta1/beta0) \
			*exp(+2.*aass*beta1/beta0*(lam+log(1.-lam))/(lam-1.)-C1/C2*beta1/beta0 * log(1.-lam) *(log(1.-lam)+2.*log( aass*C2/C1/(1.-lam + aass*C2/C1)))+2.*C1/C2*beta1/beta0*(-Li2(-C1/C2/aass)+Li2(C1/C2*(lam-1.)/aass))) \
			*exp(log(n2llaexp)*2.*beta0*(C1/C2)*(LR-LQ))	\
			*pow(aexpB,-2.*beta0*(LR-LQ)) ;

		      cdefexp(ccoeff::C1qg[idx],   ccoeff::C2qg[idx],   ccoeff::C3qg[idx]); aexpqg[idx]   = n3llaexp();
		      //cdefexp(ccoeff::C1qg[idx],   ccoeff::C2qg[idx],   0.); aexpqg[idx]   = n3llaexp();
		      cdefexp(ccoeff::C2qqb[idx],  ccoeff::C3qqb[idx],  0.               ); aexpqqb[idx]  = n3llaexp();
		      cdefexp(ccoeff::C2qqp[idx],  ccoeff::C3qqp[idx],  0.               ); aexpqqp[idx]  = n3llaexp();
		      cdefexp(ccoeff::C2qqbp[idx], ccoeff::C3qqbp[idx], 0.               ); aexpqqbp[idx] = n3llaexp();
		      
		      //aexpqg[idx]   = (1.-xlambda/(1.+aass*ccoeff::C2qg[idx]/ccoeff::C1qg[idx]))/(1.-xlambda);
		      //aexpqqb[idx]  = (1.-xlambda/(1.+aass*ccoeff::C3qqb[idx] /ccoeff::C2qqb[idx] ))/(1.-xlambda);
		      //aexpqqp[idx]  = (1.-xlambda/(1.+aass*ccoeff::C3qqp[idx] /ccoeff::C2qqp[idx] ))/(1.-xlambda);
		      //aexpqqbp[idx] = (1.-xlambda/(1.+aass*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]))/(1.-xlambda);
		    }
		}
	    }
	  /*
	  //for (int i = 0; i < mellinint::mdim; i++)
	  //{
	  if (opts.numexpc) gint::calc(b);
	  int i;
	  i = 0;
	  //  //i = opts.mellinrule - 1;
	  //cout << i << endl;
	  //  cout << aexpqq[i]   << " gint::alogqq " << endl;
	  cout << alogqg[i]   << "  " << gint::alogqg[i]  << " gint::alogqg " << endl;
	  //  cout << aexpqqb[i]  << " gint::alogqqb " << endl;
	  //  cout << aexpqqp[i]  << " gint::alogqqb " << endl;
	  //  cout << aexpqqbp[i] << " gint::alogqqbp " << endl;
	  cout << endl;
	  //}
	  */

	  //!!! gint::alog are evaluated with blim::sudakov, not with blim:expc !!!
	  if (opts.numexpc)
	    for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	      for (int i = 0; i < mellinint::mdim; i++)
		{
		  int idx = anomalous::index(i,sign);
		  
		  if (opts.sumlogs)
		    {
		      alogqq[idx]   = gint::alogqq[idx];
		      alogqg[idx]   = gint::alogqg[idx];
		      alogqqb[idx]  = gint::alogqqb[idx];
		      alogqqp[idx]  = gint::alogqqp[idx];
		      alogqqbp[idx] = gint::alogqqbp[idx];
		    }
		  else
		    {
		      aexpqq[idx]   = exp(gint::alogqq[idx]);  
		      aexpqg[idx]   = exp(gint::alogqg[idx]);  
		      aexpqqb[idx]  = exp(gint::alogqqb[idx]); 
		      aexpqqp[idx]  = exp(gint::alogqqp[idx]); 
		      aexpqqbp[idx] = exp(gint::alogqqbp[idx]);
		    }
		}
	  
	  
	  for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	    for (int i = 0; i < mellinint::mdim; i++)
	      {
		int idx = anomalous::index(i,sign);
		int ii = hcoefficients::index(i,sign);
		if (opts.sumlogs)
		  {
		    qqb[ii] = exp(alogqq[idx]+alogqq[idx]+sudakov::logS);
		    qg[ii]  = exp(-logasl+alogqg[idx]+alogqq[idx]+sudakov::logS);
		    qq[ii]  = exp(alogqq[idx]+alogqqb[idx]-2.*logasl+sudakov::logS);
		    qqp[ii] = exp(alogqq[idx]+alogqqbp[idx]-2.*logasl+sudakov::logS);
		    qqbp[ii] = exp(alogqq[idx]+alogqqp[idx]-2.*logasl+sudakov::logS);
		    gg[ii]   = exp(-logasl+alogqg[idx]-logasl+alogqg[idx]+sudakov::logS);
		    qbg[ii]  = exp(-logasl+alogqg[idx]-2.*logasl+alogqqb[idx]+sudakov::logS);
		    qpg[ii]  = exp(-logasl+alogqg[idx]-2.*logasl+alogqqp[idx]+sudakov::logS);
		    qbpg[ii] = exp(-logasl+alogqg[idx]-2.*logasl+alogqqbp[idx]+sudakov::logS);
		  }
		else
		  {
		    qqb[ii] = aexpqq[idx] * aexpqq[idx];
		    qg[ii]  = aexp * aexpqg[idx]*aexpqq[idx];
		    qq[ii]  = aexpqq[idx] * aexpqqb[idx] * aexp2; //*c1delta
		    qqp[ii] = aexpqq[idx] * aexpqqbp[idx] * aexp2; //*c1delta
		    qqbp[ii] = aexpqq[idx] * aexpqqp[idx] * aexp2; //*c1delta
		    gg[ii]  = aexp*aexpqg[idx] * aexp*aexpqg[idx];
		    qbg[ii] = aexp * aexpqg[idx] *aexp2*aexpqqb[idx];
		    qpg[ii] = aexp * aexpqg[idx] *aexp2*aexpqqp[idx];
		    qbpg[ii] = aexp * aexpqg[idx] *aexp2*aexpqqbp[idx];

		    //		  cout << endl;
		    //		  cout << bstar << "  " << i << "  " << aexp << "  " << aexpqq[idx] << "  " << aexpqg[idx] << endl;
		    //		  cout << "  " << pow(aexpB,1./2.*(C2qgM[idx]/C1QG[idx]))
		    //		       << "  " << pow(aexpC,1./4.*(pow(C2qgM[idx]/C1QG[idx],2)/2. - 8.*ccoeff::C3qg[idx]/C1QG[idx]))
		    //		       << "  " << pow(aexpC,(pow(ccoeff::C2qg[idx]/ccoeff::C1qg[idx],2)/2. - ccoeff::C3qg[idx]/ccoeff::C1qg[idx]))
		    //		       << "  " << pow(aexpD,-1./2.*beta1/beta0*C2qgM[idx]/C1QG[idx])
		    //		       << "  " << pow(aexpC,1./2.*(beta0*C2qgM[idx]/C1QG[idx])*(resint::rlogq2mur2-2.*resint::rloga))
		    //		       << endl;
		  }
		if (isnan_ofast(real(qg[ii])) || isnan_ofast(imag(qg[ii])))
		  {
		    int idx = anomalous::index(i,sign);
		    cout << "Warning, expc::qg[" << i << "] = " << qg[ii] << ", b = "  << b << ", lambda = " << xlambda << " logasl " << logasl
			 << " aexp " << aexp << " aexpqq " << aexpqq[idx] << " aexpqg " << aexpqg[idx]
			 << endl;
		  }
		
	    }
	  delete[] alogqq;
	  delete[] alogqg;
	  delete[] alogqqb;
	  delete[] alogqqp;
	  delete[] alogqqbp;
	  delete[] aexpqq;
	  delete[] aexpqg;
	  delete[] aexpqqb;
	  delete[] aexpqqp;
	  delete[] aexpqqbp;
	}
      else //mellin2d
	{
	  complex <double> aexpqq_1[mellinint::mdim*2];
	  complex <double> aexpqg_1[mellinint::mdim*2];
	  complex <double> aexpqqb_1[mellinint::mdim*2];
	  complex <double> aexpqqp_1[mellinint::mdim*2];
	  complex <double> aexpqqbp_1[mellinint::mdim*2];
	  complex <double> aexpqq_2[mellinint::mdim*2];
	  complex <double> aexpqg_2[mellinint::mdim*2];
	  complex <double> aexpqqb_2[mellinint::mdim*2];
	  complex <double> aexpqqp_2[mellinint::mdim*2];
	  complex <double> aexpqqbp_2[mellinint::mdim*2];

	  //case opts.expc == 0
	  fill(aexpqq_1, aexpqq_1+mellinint::mdim*2, 1.);
	  fill(aexpqg_1, aexpqg_1+mellinint::mdim*2, 1.);
	  fill(aexpqqb_1, aexpqqb_1+mellinint::mdim*2, 1.);
	  fill(aexpqqp_1, aexpqqp_1+mellinint::mdim*2, 1.);
	  fill(aexpqqbp_1, aexpqqbp_1+mellinint::mdim*2, 1.);
	  fill(aexpqq_2, aexpqq_2+mellinint::mdim*2, 1.);
	  fill(aexpqg_2, aexpqg_2+mellinint::mdim*2, 1.);
	  fill(aexpqqb_2, aexpqqb_2+mellinint::mdim*2, 1.);
	  fill(aexpqqp_2, aexpqqp_2+mellinint::mdim*2, 1.);
	  fill(aexpqqbp_2, aexpqqbp_2+mellinint::mdim*2, 1.);
	  

	  if (opts.expc == 1)
	    {
	      for (int sign = mesq::positive; sign <= mesq::negative; sign++)
		for (int i = 0; i < mellinint::mdim; i++)
		  {
		    int idx = anomalous::index(i,sign);
		    aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_delta)
		      *pow(aexpC,(0.5*pow(ccoeff::C1qq_delta,2) - ccoeff::C2qq_delta))
		      *pow(aexpD,-beta1/beta0*ccoeff::C1qq_delta)
		      *pow(aexpC,beta0*ccoeff::C1qq_delta*(resint::rlogq2mur2-2.*resint::rloga)) //*(LR-LQ)
		      ;
		    aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_delta)
		      *pow(aexpC,(0.5*pow(ccoeff::C1qq_delta,2) - ccoeff::C2qq_delta))
		      *pow(aexpD,-beta1/beta0*ccoeff::C1qq_delta)
		      *pow(aexpC,beta0*ccoeff::C1qq_delta*(resint::rlogq2mur2-2.*resint::rloga)) //*(LR-LQ)
		      ;
		  }
	    }
	  /*
	  else if (opts.expc == 2)
	    {
	      for (int sign = mesq::positive; sign <= mesq::negative; sign++)
		for (int i = 0; i < mellinint::mdim; i++)
		  {
		    int idx = anomalous::index(i,sign);
		    complex <double> qqtayl = 1.;
		    complex <double> qgtayl = 1.;
		    complex <double> qqbtayl = 1.;
		    complex <double> qqptayl = 1.;
		    complex <double> qqbptayl = 1.;
		    int fact = 1;
		    for (int t = 1; t <= opts.ntaylor; t++)
		      {
			fact *= t;
			qqtayl += 1./double(fact)* (pow(lamB*(ccoeff::C1qq[idx]-ccoeff::C1qq_delta),t)
						    +pow(lamC*(0.5*pow((ccoeff::C1qq[idx]-ccoeff::C1qq_delta),2) - (ccoeff::C2qq[idx]-ccoeff::C2qq_delta)),t)
						    +pow(-lamD*beta1/beta0*(ccoeff::C1qq[idx]-ccoeff::C1qq_delta),t)
						    //+pow(lamC*beta0*(ccoeff::C1qq[idx]-ccoeff::C1qq_delta)*(resint::rlogq2mur2-2.*resint::rloga),t)
						    );
						
			qgtayl += 1./double(fact)* (pow(lamB*ccoeff::C2qg[idx]/ccoeff::C1qg[idx],t)
						    +pow(lamC*(pow(ccoeff::C2qg[idx]/ccoeff::C1qg[idx],2)/2. - ccoeff::C3qg[idx]/ccoeff::C1qg[idx]),t)
						    +pow(-lamD*beta1/beta0*ccoeff::C2qg[idx]/ccoeff::C1qg[idx],t)
						    //+pow(lamC*(beta0*ccoeff::C2qg[idx]/ccoeff::C1qg[idx])*(resint::rlogq2mur2-2.*resint::rloga),t)
						    );
		    
			qqbtayl += 1./double(fact) *(pow(lamB*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx],t)
						     +pow(lamC*(pow(ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx],2)/2.),t)
						     -pow(lamD*beta1/beta0*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx],t)
						     //+pow(lamC*beta0*ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]*(resint::rlogq2mur2-2.*resint::rloga),t)
						     );

			qqptayl += 1./double(fact) *(pow(lamB*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx],t)
						     +pow(lamC*(pow(ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx],2)/2.),t)
						     -pow(lamD*beta1/beta0*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx],t)
						     //+pow(lamC*beta0*ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]*(resint::rlogq2mur2-2.*resint::rloga),t)
						     );

			qqbptayl += 1./double(fact) *(pow(lamB*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx],t)
						      +pow(lamC*(pow(ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx],2)/2.),t)
						      -pow(lamD*beta1/beta0*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx],t)
						      //+pow(lamC*beta0*ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]*(resint::rlogq2mur2-2.*resint::rloga),t)
						      );
		      }
		    aexpqq[idx] = pow(aexpB,ccoeff::C1qq_delta)
		      *pow(aexpC,(0.5*pow(ccoeff::C1qq_delta,2) - ccoeff::C2qq_delta))
		      *pow(aexpD,-beta1/beta0*ccoeff::C1qq_delta)
		      *pow(aexpC,beta0*ccoeff::C1qq_delta*(resint::rlogq2mur2-2.*resint::rloga))
		      *qqtayl;
		    aexpqg[idx] = qgtayl;
		    aexpqqb[idx] = qqbtayl;
		    aexpqqp[idx] = qqptayl;
		    aexpqqbp[idx] = qqbptayl;
		  }
	    }
	  */
	  else if (opts.expc == 3)
	    for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	      for (int i = 0; i < mellinint::mdim; i++)
		{
		  int idx = anomalous::index(i,sign);
		  
		  aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_1[idx])
		    *pow(aexpC,(0.5*pow(ccoeff::C1qq_1[idx],2) - ccoeff::C2qq_1[idx]))
		    *pow(aexpD,-beta1/beta0*ccoeff::C1qq_1[idx])
		    *pow(aexpC,beta0*ccoeff::C1qq_1[idx]*(LR-LQ))
		    ;
		  aexpqg_1[idx] = pow(aexpB,ccoeff::C2qg_1[idx]/ccoeff::C1qg_1[idx])
		    *pow(aexpC,pow(ccoeff::C2qg_1[idx]/ccoeff::C1qg_1[idx],2)/2. - ccoeff::C3qg_1[idx]/ccoeff::C1qg_1[idx])
		    *pow(aexpD,-beta1/beta0*ccoeff::C2qg_1[idx]/ccoeff::C1qg_1[idx])
		    *pow(aexpC,beta0*ccoeff::C2qg_1[idx]/ccoeff::C1qg_1[idx]*(LR-LQ))
		    ;
		  aexpqqb_1[idx] = pow(aexpB,ccoeff::C3qqb_1[idx]/ccoeff::C2qqb_1[idx])
		    *pow(aexpC,pow(ccoeff::C3qqb_1[idx]/ccoeff::C2qqb_1[idx],2)/2.)// - C4/C2)
		    *pow(aexpD,-beta1/beta0*ccoeff::C3qqb_1[idx]/ccoeff::C2qqb_1[idx])
		    *pow(aexpC,beta0*ccoeff::C3qqb_1[idx]/ccoeff::C2qqb_1[idx]*(LR-LQ))
		    ;
		  aexpqqp_1[idx] = pow(aexpB,ccoeff::C3qqp_1[idx]/ccoeff::C2qqp_1[idx])
		    *pow(aexpC,pow(ccoeff::C3qqp_1[idx]/ccoeff::C2qqp_1[idx],2)/2.)// - C4/C2))
		    *pow(aexpD,-beta1/beta0*ccoeff::C3qqp_1[idx]/ccoeff::C2qqp_1[idx])
		    *pow(aexpC,beta0*ccoeff::C3qqp_1[idx]/ccoeff::C2qqp_1[idx]*(LR-LQ))
		    ;
		  aexpqqbp_1[idx] = pow(aexpB,ccoeff::C3qqbp_1[idx]/ccoeff::C2qqbp_1[idx])
		    *pow(aexpC,pow(ccoeff::C3qqbp_1[idx]/ccoeff::C2qqbp_1[idx],2)/2.)// - C4/C2))
		    *pow(aexpD,-beta1/beta0*ccoeff::C3qqbp_1[idx]/ccoeff::C2qqbp_1[idx])
		    *pow(aexpC,beta0*ccoeff::C3qqbp_1[idx]/ccoeff::C2qqbp_1[idx]*(LR-LQ))
		    ;
		  aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_2[idx])
		    *pow(aexpC,(0.5*pow(ccoeff::C1qq_2[idx],2) - ccoeff::C2qq_2[idx]))
		    *pow(aexpD,-beta1/beta0*ccoeff::C1qq_2[idx])
		    *pow(aexpC,beta0*ccoeff::C1qq_2[idx]*(LR-LQ))
		    ;
		  aexpqg_2[idx] = pow(aexpB,ccoeff::C2qg_2[idx]/ccoeff::C1qg_2[idx])
		    *pow(aexpC,pow(ccoeff::C2qg_2[idx]/ccoeff::C1qg_2[idx],2)/2. - ccoeff::C3qg_2[idx]/ccoeff::C1qg_2[idx])
		    *pow(aexpD,-beta1/beta0*ccoeff::C2qg_2[idx]/ccoeff::C1qg_2[idx])
		    *pow(aexpC,beta0*ccoeff::C2qg_2[idx]/ccoeff::C1qg_2[idx]*(LR-LQ))
		    ;
		  aexpqqb_2[idx] = pow(aexpB,ccoeff::C3qqb_2[idx]/ccoeff::C2qqb_2[idx])
		    *pow(aexpC,pow(ccoeff::C3qqb_2[idx]/ccoeff::C2qqb_2[idx],2)/2.)// - C4/C2)
		    *pow(aexpD,-beta1/beta0*ccoeff::C3qqb_2[idx]/ccoeff::C2qqb_2[idx])
		    *pow(aexpC,beta0*ccoeff::C3qqb_2[idx]/ccoeff::C2qqb_2[idx]*(LR-LQ))
		    ;
		  aexpqqp_2[idx] = pow(aexpB,ccoeff::C3qqp_2[idx]/ccoeff::C2qqp_2[idx])
		    *pow(aexpC,pow(ccoeff::C3qqp_2[idx]/ccoeff::C2qqp_2[idx],2)/2.)// - C4/C2))
		    *pow(aexpD,-beta1/beta0*ccoeff::C3qqp_2[idx]/ccoeff::C2qqp_2[idx])
		    *pow(aexpC,beta0*ccoeff::C3qqp_2[idx]/ccoeff::C2qqp_2[idx]*(LR-LQ))
		    ;
		  aexpqqbp_2[idx] = pow(aexpB,ccoeff::C3qqbp_2[idx]/ccoeff::C2qqbp_2[idx])
		    *pow(aexpC,pow(ccoeff::C3qqbp_2[idx]/ccoeff::C2qqbp_2[idx],2)/2.)// - C4/C2))
		    *pow(aexpD,-beta1/beta0*ccoeff::C3qqbp_2[idx]/ccoeff::C2qqbp_2[idx])
		    *pow(aexpC,beta0*ccoeff::C3qqbp_2[idx]/ccoeff::C2qqbp_2[idx]*(LR-LQ))
		    ;

		  //cout << mellinint::Np_1[i] << "  " << aexpqg_1[idx] << "  " << aexpqg_2[idx]
		  //     << " C1 " << fabs(ccoeff::C1qg_1[idx])
		  //     << " C2 " << fabs(ccoeff::C2qg_1[idx])
		  //     << " C3 " << fabs(ccoeff::C3qg_1[idx])
		  //     << endl;
		}
	  else if (opts.expc == 4)
	    {
	      fill(aexpqqb_1, aexpqqb_1+mellinint::mdim*2, 1.);
	      fill(aexpqqp_1, aexpqqp_1+mellinint::mdim*2, 1.);
	      fill(aexpqqbp_1, aexpqqbp_1+mellinint::mdim*2, 1.);
	      fill(aexpqqb_2, aexpqqb_2+mellinint::mdim*2, 1.);
	      fill(aexpqqp_2, aexpqqp_2+mellinint::mdim*2, 1.);
	      fill(aexpqqbp_2, aexpqqbp_2+mellinint::mdim*2, 1.);
	      for (int sign = mesq::positive; sign <= mesq::negative; sign++)
		for (int i = 0; i < mellinint::mdim; i++)
		  {
		    int idx = anomalous::index(i,sign);
		    aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_1[idx])
		      *pow(aexpC,(0.5*pow(ccoeff::C1qq_1[idx],2) - ccoeff::C2qq_1[idx]))
		      *pow(aexpD,-beta1/beta0*ccoeff::C1qq_1[idx])
		      *pow(aexpC,beta0*ccoeff::C1qq_1[idx]*(LR-LQ))
		      ;
		    aexpqg_1[idx]   = (1.-xlambda/(1.+aass*ccoeff::C2qg_1[idx]/ccoeff::C1qg_1[idx]))/(1.-xlambda);
		    aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_2[idx])
		      *pow(aexpC,(0.5*pow(ccoeff::C1qq_2[idx],2) - ccoeff::C2qq_2[idx]))
		      *pow(aexpD,-beta1/beta0*ccoeff::C1qq_2[idx])
		      *pow(aexpC,beta0*ccoeff::C1qq_2[idx]*(LR-LQ))
		      ;
		    aexpqg_2[idx]   = (1.-xlambda/(1.+aass*ccoeff::C2qg_2[idx]/ccoeff::C1qg_2[idx]))/(1.-xlambda);
		  }
	    }
	  else if (opts.expc == 5)
	    for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	      for (int i = 0; i < mellinint::mdim; i++)
		{
		  int idx = anomalous::index(i,sign);
		  aexpqq_1[idx] = pow(aexpB,ccoeff::C1qq_1[idx])
		    *pow(aexpC,(0.5*pow(ccoeff::C1qq_1[idx],2) - ccoeff::C2qq_1[idx]))
		    *pow(aexpD,-beta1/beta0*ccoeff::C1qq_1[idx])
		    *pow(aexpC,beta0*ccoeff::C1qq_1[idx]*(LR-LQ))
		    ;
		  aexpqq_2[idx] = pow(aexpB,ccoeff::C1qq_2[idx])
		    *pow(aexpC,(0.5*pow(ccoeff::C1qq_2[idx],2) - ccoeff::C2qq_2[idx]))
		    *pow(aexpD,-beta1/beta0*ccoeff::C1qq_2[idx])
		    *pow(aexpC,beta0*ccoeff::C1qq_2[idx]*(LR-LQ))
		    ;

		  complex <double> C1,C2,C3,rootC,n2llaexp;
		  complex <double> lam = xlambda;

#define cdefexp(c1,c2,c3) C1 = c1; C2 = c2; C3 = c3; rootC = sqrt(pow(C2/C1,2)-4.*C3/C1); n2llaexp = (1.+aass*C2/C1/(1.-lam))/(1.+aass*C2/C1); //n2llaexp = (1.-lam+aass*C2/C1)/(1.-lam)/(1.+aass*C2/C1);
#define n3llaexp() n2llaexp						\
		    *exp(log((aass*(rootC+C2/C1)+2.)*(aass*(rootC-C2/C1)-2.*(1.-lam))/(aass*(rootC+C2/C1)+2.*(1.-lam))/(aass*(rootC-C2/C1)-2.))*C2/C1/rootC) \
		    *(pow(1.-lam,2)+aass*C2/C1*(1.-lam)+aass2*C3/C1)/pow(1.-lam,2)/(1.+aass*C2/C1+aass2*C3/C1) \
		    * exp(-log(n2llaexp)*beta1/beta0*C1/C2) *pow(aexpB,beta1/beta0) \
		    *exp(+2.*aass*beta1/beta0*(lam+log(1.-lam))/(lam-1.)-C1/C2*beta1/beta0 * log(1.-lam) *(log(1.-lam)+2.*log( aass*C2/C1/(1.-lam + aass*C2/C1)))+2.*C1/C2*beta1/beta0*(-Li2(-C1/C2/aass)+Li2(C1/C2*(lam-1.)/aass))) \
		    *exp(log(n2llaexp)*2.*beta0*(C1/C2)*(LR-LQ))	\
		    *pow(aexpB,-2.*beta0*(LR-LQ)) ;

		  //cdefexp(ccoeff::C1qg_1[idx],   ccoeff::C2qg_1[idx],   ccoeff::C3qg_1[idx]); aexpqg_1[idx]   = n3llaexp();
		  cdefexp(ccoeff::C1qg_1[idx],   ccoeff::C2qg_1[idx],   0.               ); aexpqg_1[idx]   = n3llaexp();
		  cdefexp(ccoeff::C2qqb_1[idx],  ccoeff::C3qqb_1[idx],  0.               ); aexpqqb_1[idx]  = n3llaexp();
		  cdefexp(ccoeff::C2qqp_1[idx],  ccoeff::C3qqp_1[idx],  0.               ); aexpqqp_1[idx]  = n3llaexp();
		  cdefexp(ccoeff::C2qqbp_1[idx], ccoeff::C3qqbp_1[idx], 0.               ); aexpqqbp_1[idx] = n3llaexp();

		  //cdefexp(ccoeff::C1qg_2[idx],   ccoeff::C2qg_2[idx],   ccoeff::C3qg_2[idx]); aexpqg_2[idx]   = n3llaexp();
		  cdefexp(ccoeff::C1qg_2[idx],   ccoeff::C2qg_2[idx],   0.               ); aexpqg_2[idx]   = n3llaexp();
		  cdefexp(ccoeff::C2qqb_2[idx],  ccoeff::C3qqb_2[idx],  0.               ); aexpqqb_2[idx]  = n3llaexp();
		  cdefexp(ccoeff::C2qqp_2[idx],  ccoeff::C3qqp_2[idx],  0.               ); aexpqqp_2[idx]  = n3llaexp();
		  cdefexp(ccoeff::C2qqbp_2[idx], ccoeff::C3qqbp_2[idx], 0.               ); aexpqqbp_2[idx] = n3llaexp();
		}
	  

	  for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	    for (int i1 = 0; i1 < mellinint::mdim; i1++)
	      for (int i2 = 0; i2 < mellinint::mdim; i2++)
		{
		  int idx1 = anomalous::index(i1,mesq::positive);
		  int idx2 = anomalous::index(i2,sign);
		  int ii12 = hcoefficients::index(i1,i2,sign);
		  qqb[ii12]  = aexpqq_1[idx1] * aexpqq_2[idx2];
		  qg_1[ii12] = aexp * aexpqg_1[idx1]*aexpqq_2[idx2];
		  qg_2[ii12] = aexp * aexpqg_2[idx2]*aexpqq_1[idx1];
		  qq_1[ii12] = aexpqq_1[idx1] * aexpqqb_2[idx2] * aexp2; //*c1delta
		  qq_2[ii12] = aexpqq_2[idx2] * aexpqqb_1[idx1] * aexp2; //*c1delta
		  qqp_1[ii12] = aexpqq_1[idx1] * aexpqqbp_2[idx2] * aexp2; //*c1delta
		  qqp_2[ii12] = aexpqq_2[idx2] * aexpqqbp_1[idx1] * aexp2; //*c1delta
		  qqbp_1[ii12] = aexpqq_1[idx1] * aexpqqp_2[idx2] * aexp2; //*c1delta
		  qqbp_2[ii12] = aexpqq_2[idx2] * aexpqqp_1[idx1] * aexp2; //*c1delta
		  gg[ii12]   = aexp*aexpqg_1[idx1] * aexp*aexpqg_2[idx2];
		  qbg_1[ii12] = aexp * aexpqg_1[idx1] *aexp2*aexpqqb_2[idx2];
		  qbg_2[ii12] = aexp * aexpqg_2[idx2] *aexp2*aexpqqb_1[idx1];
		  qpg_1[ii12] = aexp * aexpqg_1[idx1] *aexp2*aexpqqp_2[idx2];
		  qpg_2[ii12] = aexp * aexpqg_2[idx2] *aexp2*aexpqqp_1[idx1];
		  qbpg_1[ii12] = aexp * aexpqg_1[idx1] *aexp2*aexpqqbp_2[idx2];
		  qbpg_2[ii12] = aexp * aexpqg_2[idx2] *aexp2*aexpqqbp_1[idx1];
		}
	}
      //cout << alogqq[0] << " expc::alogqq " << endl;
      //cout << alogqg[0] << " expc::alogqg " << endl;
      //cout << alogqqb[0] << " expc::alogqqb " << endl;
      //cout << alogqqp[0] << " expc::alogqqb " << endl;
      //cout << alogqqbp[0] << " expc::alogqqbp " << endl;
      //cout << endl;
    }

  if (opts.order >= 4)
    {
      if (opts.mellin1d)
	{
	  complex <double> alogqq[mellinint::mdim*2];
	  complex <double> alogqg[mellinint::mdim*2];
	  complex <double> alogqqb[mellinint::mdim*2];
	  complex <double> alogqqp[mellinint::mdim*2];
	  complex <double> alogqqbp[mellinint::mdim*2];
	  complex <double> aexpqq[mellinint::mdim*2];
	  complex <double> aexpqg[mellinint::mdim*2];
	  complex <double> aexpqqb[mellinint::mdim*2];
	  complex <double> aexpqqp[mellinint::mdim*2];
	  complex <double> aexpqqbp[mellinint::mdim*2];

	  if (opts.expc == 0)
	    if (opts.sumlogs)
	      {
		fill(alogqq, alogqq+mellinint::mdim*2, 0.);
		fill(alogqg, alogqq+mellinint::mdim*2, 0.);
		fill(alogqqb, alogqqb+mellinint::mdim*2, 0.);
		fill(alogqqp, alogqqp+mellinint::mdim*2, 0.);
		fill(alogqqbp, alogqqbp+mellinint::mdim*2, 0.);
	      }
	    else
	      {
		fill(aexpqq, aexpqq+mellinint::mdim*2, 1.);
		fill(aexpqg, aexpqg+mellinint::mdim*2, 1.);
		fill(aexpqqb, aexpqqb+mellinint::mdim*2, 1.);
		fill(aexpqqp, aexpqqp+mellinint::mdim*2, 1.);
		fill(aexpqqbp, aexpqqbp+mellinint::mdim*2, 1.);
	      }

	  //!!! gint::alog are evaluated with blim::sudakov, not with blim:expc !!!
	  if (opts.numexpc)
	    for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	      for (int i = 0; i < mellinint::mdim; i++)
		{
		  int idx = anomalous::index(i,sign);
		  
		  if (opts.sumlogs)
		    {
		      alogqq[idx]   = gint::alogqq[idx];
		      alogqg[idx]   = gint::alogqg[idx];
		      alogqqb[idx]  = gint::alogqqb[idx];
		      alogqqp[idx]  = gint::alogqqp[idx];
		      alogqqbp[idx] = gint::alogqqbp[idx];
		    }
		  else
		    {
		      aexpqq[idx]   = exp(gint::alogqq[idx]);  
		      aexpqg[idx]   = exp(gint::alogqg[idx]);  
		      aexpqqb[idx]  = exp(gint::alogqqb[idx]); 
		      aexpqqp[idx]  = exp(gint::alogqqp[idx]); 
		      aexpqqbp[idx] = exp(gint::alogqqbp[idx]);
		    }
		}
	  
	  for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	    for (int i = 0; i < mellinint::mdim; i++)
	      {
		int idx = anomalous::index(i,sign);
		int ii = hcoefficients::index(i,sign);
		if (opts.sumlogs)
		  {
		    qqb[ii] = exp(alogqq[idx]+alogqq[idx]+sudakov::logS);
		    qg[ii]  = exp(-logasl+alogqg[idx]+alogqq[idx]+sudakov::logS);
		    qq[ii]  = exp(alogqq[idx]+alogqqb[idx]-2.*logasl+sudakov::logS);
		    qqp[ii] = exp(alogqq[idx]+alogqqbp[idx]-2.*logasl+sudakov::logS);
		    qqbp[ii] = exp(alogqq[idx]+alogqqp[idx]-2.*logasl+sudakov::logS);
		    gg[ii]   = exp(-logasl+alogqg[idx]-logasl+alogqg[idx]+sudakov::logS);
		    qbg[ii]  = exp(-logasl+alogqg[idx]-2.*logasl+alogqqb[idx]+sudakov::logS);
		    qpg[ii]  = exp(-logasl+alogqg[idx]-2.*logasl+alogqqp[idx]+sudakov::logS);
		    qbpg[ii] = exp(-logasl+alogqg[idx]-2.*logasl+alogqqbp[idx]+sudakov::logS);
		  }
		else
		  {
		    qqb[ii] = aexpqq[idx] * aexpqq[idx];
		    qg[ii]  = aexp * aexpqg[idx]*aexpqq[idx];
		    qq[ii]  = aexpqq[idx] * aexpqqb[idx] * aexp2;
		    qqp[ii] = aexpqq[idx] * aexpqqbp[idx] * aexp2;
		    qqbp[ii] = aexpqq[idx] * aexpqqp[idx] * aexp2;
		    gg[ii]  = aexp*aexpqg[idx] * aexp*aexpqg[idx];
		    qbg[ii] = aexp * aexpqg[idx] *aexp2*aexpqqb[idx];
		    qpg[ii] = aexp * aexpqg[idx] *aexp2*aexpqqp[idx];
		    qbpg[ii] = aexp * aexpqg[idx] *aexp2*aexpqqbp[idx];
		  }
	      }
	  
	}
      else //mellin2d
	{
	}
    }

  
  //  return;
  if (opts.npff == 3)
    if (opts.mellin1d)
      for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	for (int i = 0; i < mellinint::mdim; i++)
	  {
	    complex <double> N;
	    if (sign == mesq::positive)
	      N = mellinint::Np[i];
	    else
	      N = mellinint::Nm[i];
	    complex <double> S1N = 0.;//cpsi0(N)+constants::euler;
	    complex <double> dmw = exp(-pow(b,2)/2.*(opts.a2*(-2.*S1N-(4.*N+4.)/(N*(N+2.)))));
	    int ii = hcoefficients::index(i,sign);
	    qqb[ii]  *= dmw*dmw;
	    qg[ii]   *= dmw*dmw;
	    qq[ii]   *= dmw*dmw;
	    qqp[ii]  *= dmw*dmw;
	    qqbp[ii] *= dmw*dmw;
	    gg[ii]   *= dmw*dmw;
	    qbg[ii]  *= dmw*dmw;
	    qpg[ii]  *= dmw*dmw;
	    qbpg[ii] *= dmw*dmw;
	  }
    else
      for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	for (int i1 = 0; i1 < mellinint::mdim; i1++)
	  for (int i2 = 0; i2 < mellinint::mdim; i2++)
	    {
	      complex <double> N_1,N_2;
	      if (sign == mesq::positive)
		{
		  N_1 = mellinint::Np_1[i1];
		  N_2 = mellinint::Np_2[i2];
		}
	      else
		{
		  N_1 = mellinint::Nm_1[i1];
		  N_2 = mellinint::Nm_2[i2];
		}
	      complex <double> S1N_1 = cpsi0(N_1+1.)+constants::euler;
	      complex <double> S1N_2 = cpsi0(N_2+1.)+constants::euler;
	      complex <double> dmw_1 = exp(-pow(b,2)/2.*(opts.a2*(-2.*S1N_1-(4.*N_1+4.)/(N_1*(N_1+2.)))));
	      complex <double> dmw_2 = exp(-pow(b,2)/2.*(opts.a2*(-2.*S1N_2-(4.*N_2+4.)/(N_2*(N_2+2.)))));
	      int ii12 = hcoefficients::index(i1,i2,sign);
	      qqb[ii12]    *= dmw_1*dmw_2;
	      qg_1[ii12]   *= dmw_1;
	      qg_2[ii12]   *= dmw_2;
	      qq_1[ii12]   *= 1.;
	      qq_2[ii12]   *= 1.;
	      qqp_1[ii12]  *= 1.;
	      qqp_2[ii12]  *= 1.;
	      qqbp_1[ii12] *= 1.;
	      qqbp_2[ii12] *= 1.;
	      gg[ii12]     *= 1.;
	      qbg_1[ii12]  *= 1.;
	      qbg_2[ii12]  *= 1.;
	      qpg_1[ii12]  *= 1.;
	      qpg_2[ii12]  *= 1.;
	      qbpg_1[ii12] *= 1.;
	      qbpg_2[ii12] *= 1.;
	    }
  
  //check nans
  if (opts.mellin1d)
    for (int sign = mesq::positive; sign <= mesq::negative; sign++)
      for (int i = 0; i < mellinint::mdim; i++)
	{
	  int ii = hcoefficients::index(i,sign);
	  if (isnan_ofast(real(qg[ii])) || isnan_ofast(imag(qg[ii])))
	    {
	      int idx = anomalous::index(i,sign);
	      cout << "Warning, expc::qg[" << i << "] = " << qg[ii] << ", b = "  << b << ", lambda = " << xlambda << " logasl " << logasl
		//<< " aexpqq " << aexpqq[idx] << " aexpqg " << aexpqg[idx]
		   << endl;
	      if (opts.sumlogs)
		qg[ii] = sudakov::S;
	      else
		qg[ii] = 1.;
	    }
	}
}


