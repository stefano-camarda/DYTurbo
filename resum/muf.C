#include "muf.h"
#include "resint.h"
#include "pmom.h"
#include "anomalous.h"
#include "mesq.h"
#include "settings.h"
#include "scales.h"
#include "pdfevol.h"
#include "phasespace.h"
#include "parton.h"
#include "ccoeff.h"
//#include "anomalous.h"
#include "pmom.h"
#include "resconst.h"
#include "expc.h"
#include "hcoeff.h"

#include <iostream>

//LL
complex <double> *muf::qqb;
//NLL
complex <double> *muf::qg;
complex <double> *muf::qg_1;
complex <double> *muf::qg_2;
//NNLL
complex <double> *muf::qq;
complex <double> *muf::qq_1;
complex <double> *muf::qq_2;
complex <double> *muf::qqp;
complex <double> *muf::qqp_1;
complex <double> *muf::qqp_2;
complex <double> *muf::qqbp;
complex <double> *muf::qqbp_1;
complex <double> *muf::qqbp_2;
complex <double> *muf::gg;
//NNNLL
complex <double> *muf::qbg;
complex <double> *muf::qbg_1;
complex <double> *muf::qbg_2;
complex <double> *muf::qpg;
complex <double> *muf::qpg_1;
complex <double> *muf::qpg_2;
complex <double> *muf::qbpg;
complex <double> *muf::qbpg_1;
complex <double> *muf::qbpg_2;

//using namespace ccoeff;
using namespace pmom;
using namespace resconst;
using namespace resint;

  //allocate memory
void muf::allocate()
{
  //if (opts.order == 0)
  //return;

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

void muf::free()
{
  //if (opts.order == 0)
  //return;

  if (opts.mellin1d)
    {
      delete[] qqb;
      delete[] qg;
      delete[] qq;
      delete[] qqp;  
      delete[] qqbp;
      delete[] gg;
      delete[] qbg;
      delete[] qpg;
      delete[] qbpg;
    }
  else
    {
      delete[] qqb;
      delete[] qg_1;
      delete[] qg_2;
      delete[] qq_1;
      delete[] qq_2;
      delete[] qqp_1;
      delete[] qqp_2;
      delete[] qqbp_1;
      delete[] qqbp_2;
      delete[] gg;
      delete[] qbg_1;
      delete[] qbg_2;
      delete[] qpg_1;
      delete[] qpg_2;
      delete[] qbpg_1;
      delete[] qbpg_2;
    }
}

void muf::reset()
{
  //if (opts.order == 0)
  //return;

  if (opts.mellin1d)
    {
      fill(qqb , qqb  + mellinint::mdim*2, 0.);
      fill(qg  , qg   + mellinint::mdim*2, 0.);
      fill(qq  , qq   + mellinint::mdim*2, 0.);
      fill(qqp , qqp  + mellinint::mdim*2, 0.);
      fill(qqbp, qqbp + mellinint::mdim*2, 0.);
      fill(gg  , gg   + mellinint::mdim*2, 0.);
      fill(qbg , qbg  + mellinint::mdim*2, 0.);
      fill(qpg , qpg  + mellinint::mdim*2, 0.);
      fill(qbpg, qbpg + mellinint::mdim*2, 0.);
    }
  else
    {
      fill(qqb   ,qqb    + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qg_1  ,qg_1   + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qg_2  ,qg_2   + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qq_1  ,qq_1   + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qq_2  ,qq_2   + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qqp_1 ,qqp_1  + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qqp_2 ,qqp_2  + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qqbp_1,qqbp_1 + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qqbp_2,qqbp_2 + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(gg    ,gg     + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbg_1 ,qbg_1  + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbg_2 ,qbg_2  + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qpg_1 ,qpg_1  + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qpg_2 ,qpg_2  + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbpg_1,qbpg_1 + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(qbpg_2,qbpg_2 + mellinint::mdim*mellinint::mdim*2, 0.);
    }
}

void muf::calc(complex <double> b)
{
  //LO
  if (opts.order_hcoef == 0)
    return;

  //if (opts.evolmode == 1 && opts.mufevol)
  //return;
  
  double LQ2 = pow(LQ,2);
  double LQ3 = pow(LQ,3);
  double LQ4 = pow(LQ,4);
  double LQ5 = pow(LQ,5);
  double LQ6 = pow(LQ,6);

  double LR2 = pow(LR,2);

  double LF2 = pow(LF,2);
  double LF3 = pow(LF,3);
  double LQF  = LF - LQ;
  double LQF3 = pow(LQF,3);
  /*
  double LQb = LQ;
  double LFb = LF;
  if (opts.mufevol && real(b) > 0)
    {
      double m = phasespace::m;
      double mub = real(pdfevol::mub);
      double muf = scales::fac;
      double Q = scales::res;
      double m2tilde   = pow(m*mub,2)/(pow(m,2) + pow(mub,2));
      double Q2tilde   = pow(Q*mub,2)/(pow(Q,2) + pow(mub,2));
      double muf2tilde = pow(muf*mub,2)/(pow(muf,2) + pow(mub,2));
      LQb = log(m2tilde/Q2tilde);
      LFb = log(m2tilde/muf2tilde);
      
      //LQb = log( pow(phasespace::m,2)/(pow(scales::res*real(pdfevol::mub),2)/(pow(scales::res,2) + pow(real(pdfevol::mub),2))));
      //LFb = log( pow(phasespace::m,2)/(pow(scales::fac*real(pdfevol::mub),2)/(pow(scales::fac,2) + pow(real(pdfevol::mub),2))));
    }
  //LQb=0.;
  //LFb=0.;
  double LQb2 = pow(LQb,2);
  double LQb3 = pow(LQb,3);
  double LQb4 = pow(LQb,3);
  double LQb5 = pow(LQb,3);

  double LFb2 = pow(LFb,2);
  double LFb3 = pow(LFb,3);
  
  double LQFb = LFb - LQb;
  */

  double as = aass;
  double as2 = pow(aass,2);
  double as3 = pow(aass,3);
  
  //cout << "mub " << real(pdfevol::mub) << "  " << LQb << "  " << LFb << "  " << LQFb << endl;

  if (opts.mellin1d) //rapidity integrated
    for (int sign = mesq::positive; sign <= mesq::negative; sign++)
      for (int i = 0; i < mellinint::mdim; i++)
	{
	  int idx = pmom::index(i,sign);

	  //Include exponentiation in the C coefficients
	  complex <double> eC1qq   = ccoeff::C1qq[idx];
	  complex <double> eC2qq   = ccoeff::C2qq[idx];
	  complex <double> eC1qg   = ccoeff::C1qg[idx]*expc::aexpqg[idx]/expc::aexpqq[idx];
	  complex <double> eC2qg   = ccoeff::C2qg[idx]*expc::aexpqg[idx]/expc::aexpqq[idx];
	  complex <double> eC2qqp  = ccoeff::C2qqp[idx]*expc::aexpqqp[idx]/expc::aexpqq[idx];
	  complex <double> eC2qqb  = ccoeff::C2qqb[idx]*expc::aexpqqb[idx]/expc::aexpqq[idx];
	  complex <double> eC2qqbp = ccoeff::C2qqbp[idx]*expc::aexpqqbp[idx]/expc::aexpqq[idx];
	
	  //NLO
	  complex <double> H1st_qqb = 2.*gamma1qq[idx]*LQF;
	  complex <double> H1st_qg  = gamma1qg[idx]*LQF;
	
	  qqb[idx] += as*H1st_qqb;
	  qg[idx]  += as*H1st_qg;

	  if (opts.order_hcoef == 1)
	    continue;

	  //NNLO
	  complex <double> H2st_qqb  = LQF*(2.*eC1qg*gamma1gq[idx] + 2.*gamma2qq[idx] + gamma1gq[idx]*gamma1qg[idx]*LQF + gamma1qq[idx]*(4.*eC1qq + 2.*gamma1qq[idx]*LF - 2.*(B1q + gamma1qq[idx])*LQ - A1q*LQ2 + beta0*(LF + LQ)));
	  complex <double> H2st_qq   = ((2.*eC1qg*gamma1gq[idx] + 2.*gamma2qqb[idx] + gamma1gq[idx]*gamma1qg[idx]*LQF)*LQF)/2.;
	  complex <double> H2st_qqp  = ((2.*eC1qg*gamma1gq[idx] + 2.*gamma2qqp[idx] + gamma1gq[idx]*gamma1qg[idx]*LQF)*LQF)/2.;
	  complex <double> H2st_qqbp = ((2.*eC1qg*gamma1gq[idx] + 2.*gamma2qqp[idx] + gamma1gq[idx]*gamma1qg[idx]*LQF)*LQF)/2.;
	  complex <double> H2st_qg   = (LQF*(4.*eC1qq*gamma1qg[idx] + 2.*eC1qg*(gamma1gg[idx] + gamma1qq[idx]) + 2.*gamma2qg[idx] + beta0*gamma1qg[idx]*LF + gamma1gg[idx]*gamma1qg[idx]*LF + 3.*gamma1qg[idx]*gamma1qq[idx]*LF - 2.*B1q*gamma1qg[idx]*LQ + beta0*gamma1qg[idx]*LQ - gamma1gg[idx]*gamma1qg[idx]*LQ - 3.*gamma1qg[idx]*gamma1qq[idx]*LQ - A1q*gamma1qg[idx]*LQ2))/2.;
	  complex <double> H2st_gg   = gamma1qg[idx]*(2.*eC1qg + gamma1qg[idx]*LQF)*LQF;
	  
	  //Mixed LR*LQF variations
	  H2st_qqb += -beta0*H1st_qqb*LR;
	  H2st_qg  += -beta0*H1st_qg*LR;
	
	  qqb[idx]  += as2*H2st_qqb;
	  qg[idx]   += as2*H2st_qg;
	  qq[idx]   += as2*H2st_qq;
	  qqp[idx]  += as2*H2st_qqp;
	  qqbp[idx] += as2*H2st_qqbp;
	  gg[idx]   += as2*H2st_gg;

	  if (opts.order_hcoef == 2)
	    continue;

	  //NNNLO
	  double nf = 5.;
	  complex <double> H3st_qqb  = (LQF*(24.*eC2qg*gamma1gq[idx] + 24.*pow(eC1qq,2)*gamma1qq[idx] + 48.*eC2qq*gamma1qq[idx] + 48.*eC1qq*gamma2qq[idx] + 24.*gamma3qq[idx] + 24.*eC1qq*gamma1gq[idx]*gamma1qg[idx]*LF + 12.*beta1*gamma1qq[idx]*LF + 24.*beta0*eC1qq*gamma1qq[idx]*LF + 48.*eC1qq*pow(gamma1qq[idx],2)*LF + 12.*gamma1qg[idx]*gamma2gq[idx]*LF + 12.*gamma1gq[idx]*gamma2qg[idx]*LF + 24.*beta0*gamma2qq[idx]*LF + 48.*gamma1qq[idx]*gamma2qq[idx]*LF + 12.*beta0*gamma1gq[idx]*gamma1qg[idx]*LF2 + 4.*gamma1gg[idx]*gamma1gq[idx]*gamma1qg[idx]*LF2 + 8.*pow(beta0,2)*gamma1qq[idx]*LF2 + 20.*gamma1gq[idx]*gamma1qg[idx]*gamma1qq[idx]*LF2 + 24.*beta0*pow(gamma1qq[idx],2)*LF2 + 16.*pow(gamma1qq[idx],3)*LF2 - 24.*eC1qq*gamma1gq[idx]*gamma1qg[idx]*LQ - 24.*B2q*gamma1qq[idx]*LQ + 12.*beta1*gamma1qq[idx]*LQ - 48.*B1q*eC1qq*gamma1qq[idx]*LQ + 72.*beta0*eC1qq*gamma1qq[idx]*LQ - 48.*eC1qq*pow(gamma1qq[idx],2)*LQ - 12.*gamma1qg[idx]*gamma2gq[idx]*LQ - 12.*gamma1gq[idx]*gamma2qg[idx]*LQ - 24.*B1q*gamma2qq[idx]*LQ + 24.*beta0*gamma2qq[idx]*LQ - 48.*gamma1qq[idx]*gamma2qq[idx]*LQ - 12.*B1q*gamma1gq[idx]*gamma1qg[idx]*LF*LQ - 8.*gamma1gg[idx]*gamma1gq[idx]*gamma1qg[idx]*LF*LQ - 12.*B1q*beta0*gamma1qq[idx]*LF*LQ + 8.*pow(beta0,2)*gamma1qq[idx]*LF*LQ - 40.*gamma1gq[idx]*gamma1qg[idx]*gamma1qq[idx]*LF*LQ - 24.*B1q*pow(gamma1qq[idx],2)*LF*LQ - 32.*pow(gamma1qq[idx],3)*LF*LQ + 12.*B1q*gamma1gq[idx]*gamma1qg[idx]*LQ2 - 12.*beta0*gamma1gq[idx]*gamma1qg[idx]*LQ2 + 4.*gamma1gg[idx]*gamma1gq[idx]*gamma1qg[idx]*LQ2 - 12.*A2q*gamma1qq[idx]*LQ2 + 12.*pow(B1q,2)*gamma1qq[idx]*LQ2 - 24.*B1q*beta0*gamma1qq[idx]*LQ2 + 8.*pow(beta0,2)*gamma1qq[idx]*LQ2 - 24.*A1q*eC1qq*gamma1qq[idx]*LQ2 + 20.*gamma1gq[idx]*gamma1qg[idx]*gamma1qq[idx]*LQ2 + 24.*B1q*pow(gamma1qq[idx],2)*LQ2 - 24.*beta0*pow(gamma1qq[idx],2)*LQ2 + 16.*pow(gamma1qq[idx],3)*LQ2 - 12.*A1q*gamma2qq[idx]*LQ2 - 6.*A1q*gamma1gq[idx]*gamma1qg[idx]*LF*LQ2 - 6.*A1q*beta0*gamma1qq[idx]*LF*LQ2 - 12.*A1q*pow(gamma1qq[idx],2)*LF*LQ2 + 6.*A1q*gamma1gq[idx]*gamma1qg[idx]*LQ3 + 12.*A1q*B1q*gamma1qq[idx]*LQ3 - 14.*A1q*beta0*gamma1qq[idx]*LQ3 + 12.*A1q*pow(gamma1qq[idx],2)*LQ3 + 3.*pow(A1q,2)*gamma1qq[idx]*LQ4 + 12.*eC1qg*(2.*eC1qq*gamma1gq[idx] + 2.*gamma2gq[idx] + gamma1gq[idx]*(beta0*LF + gamma1gg[idx]*LF + 3.*gamma1qq[idx]*LF - 2.*B1q*LQ + 3.*beta0*LQ - gamma1gg[idx]*LQ - 3.*gamma1qq[idx]*LQ - A1q*LQ2))))/12.;
	  complex <double> H3st_qq   = (LQF*(12.*(eC2qg*gamma1gq[idx] + 2.*eC2qqb*gamma1qq[idx] + 2.*eC1qq*gamma2qqb[idx] + gamma3qqb[idx]) + 6.*(2.*eC1qq*gamma1gq[idx]*gamma1qg[idx] + gamma1qg[idx]*gamma2gq[idx] + gamma1gq[idx]*gamma2qg[idx] + 2.*beta0*gamma2qqb[idx] + 4.*gamma1qq[idx]*gamma2qqb[idx])*LF + 2.*gamma1gq[idx]*gamma1qg[idx]*(3.*beta0 + gamma1gg[idx] + 5.*gamma1qq[idx])*LF2 - 2.*(6.*eC1qq*gamma1gq[idx]*gamma1qg[idx] + 3.*gamma1qg[idx]*gamma2gq[idx] + 3.*gamma1gq[idx]*gamma2qg[idx] + 6.*(B1q - beta0 + 2.*gamma1qq[idx])*gamma2qqb[idx] + gamma1gq[idx]*gamma1qg[idx]*(3.*B1q + 2.*gamma1gg[idx] + 10.*gamma1qq[idx])*LF)*LQ + (-6.*A1q*gamma2qqb[idx] + gamma1gq[idx]*gamma1qg[idx]*(6.*B1q - 6.*beta0 + 2.*gamma1gg[idx] + 10.*gamma1qq[idx] - 3.*A1q*LF))*LQ2 + 3.*A1q*gamma1gq[idx]*gamma1qg[idx]*LQ3 + 6.*eC1qg*(2.*gamma2gq[idx] + gamma1gq[idx]*(2.*eC1qq + (beta0 + gamma1gg[idx] + 3.*gamma1qq[idx])*LF - (2.*B1q - 3.*beta0 + gamma1gg[idx] + 3.*gamma1qq[idx])*LQ - A1q*LQ2))))/12.;
	  complex <double> H3st_qqp  = (LQF*(12.*eC2qg*gamma1gq[idx] + 24.*eC2qqbp*gamma1qq[idx] + 24.*eC1qq*gamma2qqbp[idx] + 12.*gamma3qqbp[idx] + 12.*eC1qq*gamma1gq[idx]*gamma1qg[idx]*LF + 6.*gamma1qg[idx]*gamma2gq[idx]*LF + 6.*gamma1gq[idx]*gamma2qg[idx]*LF + 12.*beta0*gamma2qqbp[idx]*LF + 24.*gamma1qq[idx]*gamma2qqbp[idx]*LF + 6.*beta0*gamma1gq[idx]*gamma1qg[idx]*LF2 + 2.*gamma1gg[idx]*gamma1gq[idx]*gamma1qg[idx]*LF2 + 10.*gamma1gq[idx]*gamma1qg[idx]*gamma1qq[idx]*LF2 - 12.*eC1qq*gamma1gq[idx]*gamma1qg[idx]*LQ - 6.*gamma1qg[idx]*gamma2gq[idx]*LQ - 6.*gamma1gq[idx]*gamma2qg[idx]*LQ - 12.*B1q*gamma2qqbp[idx]*LQ + 12.*beta0*gamma2qqbp[idx]*LQ - 24.*gamma1qq[idx]*gamma2qqbp[idx]*LQ - 6.*B1q*gamma1gq[idx]*gamma1qg[idx]*LF*LQ - 4.*gamma1gg[idx]*gamma1gq[idx]*gamma1qg[idx]*LF*LQ - 20.*gamma1gq[idx]*gamma1qg[idx]*gamma1qq[idx]*LF*LQ + 6.*B1q*gamma1gq[idx]*gamma1qg[idx]*LQ2 - 6.*beta0*gamma1gq[idx]*gamma1qg[idx]*LQ2 + 2.*gamma1gg[idx]*gamma1gq[idx]*gamma1qg[idx]*LQ2 + 10.*gamma1gq[idx]*gamma1qg[idx]*gamma1qq[idx]*LQ2 - 6.*A1q*gamma2qqbp[idx]*LQ2 - 3.*A1q*gamma1gq[idx]*gamma1qg[idx]*LF*LQ2 + 3.*A1q*gamma1gq[idx]*gamma1qg[idx]*LQ3 + 6.*eC1qg*(2.*eC1qq*gamma1gq[idx] + 2.*gamma2gq[idx] + gamma1gq[idx]*(beta0*LF + gamma1gg[idx]*LF + 3.*gamma1qq[idx]*LF - 2.*B1q*LQ + 3.*beta0*LQ - gamma1gg[idx]*LQ - 3.*gamma1qq[idx]*LQ - A1q*LQ2))))/12.;
	  complex <double> H3st_qqbp = (LQF*(12.*(eC2qg*gamma1gq[idx] + 2.*eC2qqp*gamma1qq[idx] + 2.*eC1qq*gamma2qqp[idx] + gamma3qqp[idx]) + 6.*(2.*eC1qq*gamma1gq[idx]*gamma1qg[idx] + gamma1qg[idx]*gamma2gq[idx] + gamma1gq[idx]*gamma2qg[idx] + 2.*beta0*gamma2qqp[idx] + 4.*gamma1qq[idx]*gamma2qqp[idx])*LF + 2.*gamma1gq[idx]*gamma1qg[idx]*(3.*beta0 + gamma1gg[idx] + 5.*gamma1qq[idx])*LF2 - 2.*(6.*eC1qq*gamma1gq[idx]*gamma1qg[idx] + 3.*gamma1qg[idx]*gamma2gq[idx] + 3.*gamma1gq[idx]*gamma2qg[idx] + 6.*(B1q - beta0 + 2.*gamma1qq[idx])*gamma2qqp[idx] + gamma1gq[idx]*gamma1qg[idx]*(3.*B1q + 2.*gamma1gg[idx] + 10.*gamma1qq[idx])*LF)*LQ + (-6.*A1q*gamma2qqp[idx] + gamma1gq[idx]*gamma1qg[idx]*(6.*B1q - 6.*beta0 + 2.*gamma1gg[idx] + 10.*gamma1qq[idx] - 3.*A1q*LF))*LQ2 + 3.*A1q*gamma1gq[idx]*gamma1qg[idx]*LQ3 + 6.*eC1qg*(2.*gamma2gq[idx] + gamma1gq[idx]*(2.*eC1qq + (beta0 + gamma1gg[idx] + 3.*gamma1qq[idx])*LF - (2.*B1q - 3.*beta0 + gamma1gg[idx] + 3.*gamma1qq[idx])*LQ - A1q*LQ2))))/12.;
	  complex <double> H3st_qg   = (LQF*(24.*pow(eC1qg,2)*gamma1gq[idx] + 24.*eC2qg*(gamma1gg[idx] + gamma1qq[idx]) + 24.*((pow(eC1qq,2) + 2.*eC2qq + eC2qqb - eC2qqbp - eC2qqp)*gamma1qg[idx] + 2.*eC1qq*gamma2qg[idx] + gamma3qg[idx]) + 12.*(beta1 + 2.*beta0*eC1qq)*gamma1qg[idx]*LF + 24.*eC1qq*gamma1gg[idx]*gamma1qg[idx]*LF + 72.*eC1qq*gamma1qg[idx]*gamma1qq[idx]*LF + 12.*gamma1qg[idx]*gamma2gg[idx]*LF + 24.*beta0*gamma2qg[idx]*LF + 12.*gamma1gg[idx]*gamma2qg[idx]*LF + 36.*gamma1qq[idx]*gamma2qg[idx]*LF + 36.*gamma1qg[idx]*gamma2qq[idx]*LF + 12.*gamma1qg[idx]*gamma2qqb[idx]*LF - 12.*gamma1qg[idx]*gamma2qqbp[idx]*LF - 12.*gamma1qg[idx]*gamma2qqp[idx]*LF + 8.*pow(beta0,2)*gamma1qg[idx]*LF2 + 12.*beta0*gamma1gg[idx]*gamma1qg[idx]*LF2 + 4.*pow(gamma1gg[idx],2)*gamma1qg[idx]*LF2 + 12.*gamma1gq[idx]*pow(gamma1qg[idx],2)*LF2 + 36.*beta0*gamma1qg[idx]*gamma1qq[idx]*LF2 + 16.*gamma1gg[idx]*gamma1qg[idx]*gamma1qq[idx]*LF2 + 28.*gamma1qg[idx]*pow(gamma1qq[idx],2)*LF2 - 24.*B2q*gamma1qg[idx]*LQ + 12.*beta1*gamma1qg[idx]*LQ - 48.*B1q*eC1qq*gamma1qg[idx]*LQ + 72.*beta0*eC1qq*gamma1qg[idx]*LQ - 24.*eC1qq*gamma1gg[idx]*gamma1qg[idx]*LQ - 72.*eC1qq*gamma1qg[idx]*gamma1qq[idx]*LQ - 12.*gamma1qg[idx]*gamma2gg[idx]*LQ - 24.*B1q*gamma2qg[idx]*LQ + 24.*beta0*gamma2qg[idx]*LQ - 12.*gamma1gg[idx]*gamma2qg[idx]*LQ - 36.*gamma1qq[idx]*gamma2qg[idx]*LQ - 36.*gamma1qg[idx]*gamma2qq[idx]*LQ - 12.*gamma1qg[idx]*gamma2qqb[idx]*LQ + 12.*gamma1qg[idx]*gamma2qqbp[idx]*LQ + 12.*gamma1qg[idx]*gamma2qqp[idx]*LQ - 12.*B1q*beta0*gamma1qg[idx]*LF*LQ + 8.*pow(beta0,2)*gamma1qg[idx]*LF*LQ - 12.*B1q*gamma1gg[idx]*gamma1qg[idx]*LF*LQ - 8.*pow(gamma1gg[idx],2)*gamma1qg[idx]*LF*LQ - 24.*gamma1gq[idx]*pow(gamma1qg[idx],2)*LF*LQ - 36.*B1q*gamma1qg[idx]*gamma1qq[idx]*LF*LQ - 32.*gamma1gg[idx]*gamma1qg[idx]*gamma1qq[idx]*LF*LQ - 56.*gamma1qg[idx]*pow(gamma1qq[idx],2)*LF*LQ - 12.*A2q*gamma1qg[idx]*LQ2 + 12.*pow(B1q,2)*gamma1qg[idx]*LQ2 - 24.*B1q*beta0*gamma1qg[idx]*LQ2 + 8.*pow(beta0,2)*gamma1qg[idx]*LQ2 - 24.*A1q*eC1qq*gamma1qg[idx]*LQ2 + 12.*B1q*gamma1gg[idx]*gamma1qg[idx]*LQ2 - 12.*beta0*gamma1gg[idx]*gamma1qg[idx]*LQ2 + 4.*pow(gamma1gg[idx],2)*gamma1qg[idx]*LQ2 + 12.*gamma1gq[idx]*pow(gamma1qg[idx],2)*LQ2 + 36.*B1q*gamma1qg[idx]*gamma1qq[idx]*LQ2 - 36.*beta0*gamma1qg[idx]*gamma1qq[idx]*LQ2 + 16.*gamma1gg[idx]*gamma1qg[idx]*gamma1qq[idx]*LQ2 + 28.*gamma1qg[idx]*pow(gamma1qq[idx],2)*LQ2 - 12.*A1q*gamma2qg[idx]*LQ2 - 6.*A1q*beta0*gamma1qg[idx]*LF*LQ2 - 6.*A1q*gamma1gg[idx]*gamma1qg[idx]*LF*LQ2 - 18.*A1q*gamma1qg[idx]*gamma1qq[idx]*LF*LQ2 + 12.*A1q*B1q*gamma1qg[idx]*LQ3 - 14.*A1q*beta0*gamma1qg[idx]*LQ3 + 6.*A1q*gamma1gg[idx]*gamma1qg[idx]*LQ3 + 18.*A1q*gamma1qg[idx]*gamma1qq[idx]*LQ3 + 3.*pow(A1q,2)*gamma1qg[idx]*LQ4 + 4.*gamma1qg[idx]*(6.*eC2qqbp + 6.*eC2qqp + (3.*gamma2qqbp[idx] + 3.*gamma2qqp[idx] + 2.*gamma1gq[idx]*gamma1qg[idx]*LQF)*LQF)*nf + 12.*eC1qg*(2.*eC1qq*(gamma1gg[idx] + gamma1qq[idx]) + 2.*(gamma2gg[idx] + gamma2qq[idx]) + pow(gamma1gg[idx],2)*LF + 3.*gamma1gq[idx]*gamma1qg[idx]*LF + 2.*gamma1gg[idx]*gamma1qq[idx]*LF + pow(gamma1qq[idx],2)*LF - 2.*B1q*gamma1gg[idx]*LQ - pow(gamma1gg[idx],2)*LQ - 3.*gamma1gq[idx]*gamma1qg[idx]*LQ - 2.*B1q*gamma1qq[idx]*LQ - 2.*gamma1gg[idx]*gamma1qq[idx]*LQ - pow(gamma1qq[idx],2)*LQ - A1q*gamma1gg[idx]*LQ2 - A1q*gamma1qq[idx]*LQ2 + beta0*(gamma1gg[idx] + gamma1qq[idx])*(LF + 3.*LQ) + 2.*gamma1gq[idx]*gamma1qg[idx]*LQF*nf)))/24.;
	  complex <double> H3st_qbg  = ((2.*eC2qqb*gamma1qg[idx] + (eC1qg + gamma1qg[idx]*LQF)*(2.*eC1qg*gamma1gq[idx] + 2.*gamma2qqb[idx] + gamma1gq[idx]*gamma1qg[idx]*LQF))*LQF)/2.;
	  complex <double> H3st_qpg  = ((2.*eC2qqp*gamma1qg[idx] + (eC1qg + gamma1qg[idx]*LQF)*(2.*eC1qg*gamma1gq[idx] + 2.*gamma2qqp[idx] + gamma1gq[idx]*gamma1qg[idx]*LQF))*LQF)/2.;
	  complex <double> H3st_qbpg = ((2.*eC2qqbp*gamma1qg[idx] + (eC1qg + gamma1qg[idx]*LQF)*(2.*eC1qg*gamma1gq[idx] + 2.*gamma2qqbp[idx] + gamma1gq[idx]*gamma1qg[idx]*LQF))*LQF)/2.;
	  complex <double> H3st_gg   = (LQF*(4.*pow(eC1qg,2)*gamma1gg[idx] + 4.*eC2qg*gamma1qg[idx] + 2.*eC1qg*(2.*eC1qq*gamma1qg[idx] + 2.*gamma2qg[idx] + gamma1qg[idx]*(beta0*LF + 3.*gamma1gg[idx]*LF + gamma1qq[idx]*LF - 2.*B1q*LQ + 3.*beta0*LQ - 3.*gamma1gg[idx]*LQ - gamma1qq[idx]*LQ - A1q*LQ2)) + gamma1qg[idx]*LQF*(4.*eC1qq*gamma1qg[idx] + 4.*gamma2qg[idx] + gamma1qg[idx]*(2.*gamma1gg[idx]*LF + 2.*gamma1qq[idx]*LF - 2.*B1q*LQ - 2.*gamma1gg[idx]*LQ - 2.*gamma1qq[idx]*LQ - A1q*LQ2 + 2.*beta0*(LF + LQ)))))/2.;
	  
	  //Mixed LR*LQF variations
	  H3st_qqb  += -2*beta0*H2st_qqb*LR + (-beta1*LR - pow(beta0,2)*pow(LR,2))*H1st_qqb;
	  H3st_qg   += -2*beta0*H2st_qg*LR  + (-beta1*LR - pow(beta0,2)*pow(LR,2))*H1st_qg;
	  H3st_qq   += -2*beta0*H2st_qq*LR;
	  H3st_qqp  += -2*beta0*H2st_qqp*LR;
	  H3st_qqbp += -2*beta0*H2st_qqbp*LR;
	  H3st_gg   += -2*beta0*H2st_gg*LR;
	
	  qqb[idx]  += as3*H3st_qqb;
	  qg[idx]   += as3*H3st_qg;
	  qq[idx]   += as3*H3st_qq;
	  qqp[idx]  += as3*H3st_qqp;
	  qqbp[idx] += as3*H3st_qqbp;
	  gg[idx]   += as3*H3st_gg;
	  qbg[idx]  += as3*H3st_qbg;
	  qpg[idx]  += as3*H3st_qpg;
	  qbpg[idx] += as3*H3st_qbpg;
	}
  else //rapidity dependent
    for (int sign = mesq::positive; sign <= mesq::negative; sign++)
      for (int i1 = 0; i1 < mellinint::mdim; i1++)
	for (int i2 = 0; i2 < mellinint::mdim; i2++)
	  {
	    int ii1 = anomalous::index(i1,mesq::positive);
	    int ii2 = anomalous::index(i2,sign);
	    int idx = hcoeff::index(i1,i2,sign);

	    //Include exponentiation in the C coefficients
	    complex <double> eC1qq_1   = ccoeff::C1qq_1[ii1];
	    complex <double> eC2qq_1   = ccoeff::C2qq_1[ii1];
	    complex <double> eC1qg_1   = ccoeff::C1qg_1[ii1]*expc::aexpqg_1[ii1]/expc::aexpqq_1[ii1];
	    complex <double> eC2qg_1   = ccoeff::C2qg_1[ii1]*expc::aexpqg_1[ii1]/expc::aexpqq_1[ii1];
	    complex <double> eC2qqp_1  = ccoeff::C2qqp_1[ii1]*expc::aexpqqp_1[ii1]/expc::aexpqq_1[ii1];
	    complex <double> eC2qqb_1  = ccoeff::C2qqb_1[ii1]*expc::aexpqqb_1[ii1]/expc::aexpqq_1[ii1];
	    complex <double> eC2qqbp_1 = ccoeff::C2qqbp_1[ii1]*expc::aexpqqbp_1[ii1]/expc::aexpqq_1[ii1];
	    complex <double> eC1qq_2   = ccoeff::C1qq_2[ii2];
	    complex <double> eC2qq_2   = ccoeff::C2qq_2[ii2];
	    complex <double> eC1qg_2   = ccoeff::C1qg_2[ii2]*expc::aexpqg_2[ii2]/expc::aexpqq_2[ii2];
	    complex <double> eC2qg_2   = ccoeff::C2qg_2[ii2]*expc::aexpqg_2[ii2]/expc::aexpqq_2[ii2];
	    complex <double> eC2qqp_2  = ccoeff::C2qqp_2[ii2]*expc::aexpqqp_2[ii2]/expc::aexpqq_2[ii2];
	    complex <double> eC2qqb_2  = ccoeff::C2qqb_2[ii2]*expc::aexpqqb_2[ii2]/expc::aexpqq_2[ii2];
	    complex <double> eC2qqbp_2 = ccoeff::C2qqbp_2[ii2]*expc::aexpqqbp_2[ii2]/expc::aexpqq_2[ii2];

	    
	    complex <double> H1st_qqb  = (gamma1qq_1[ii1] + gamma1qq_2[ii2])*LQF ;
	    complex <double> H1st_qg_1 = gamma1qg_1[ii1]*LQF;
	    complex <double> H1st_qg_2 = gamma1qg_2[ii2]*LQF;

	    qqb[idx]   += as*H1st_qqb;
	    qg_1[idx]  += as*H1st_qg_1;
	    qg_2[idx]  += as*H1st_qg_2;
	    
	    if (opts.order_hcoef == 1)
	      continue;

	    complex <double> H2st_qqb    = LQF*(eC1qg_1*gamma1gq_1[ii1] + eC1qq_2*gamma1qq_1[ii1] + eC1qg_2*gamma1gq_2[ii2] + eC1qq_2*gamma1qq_2[ii2] + eC1qq_1*(gamma1qq_1[ii1] + gamma1qq_2[ii2]) + gamma2qq_1[ii1] + gamma2qq_2[ii2] + ((gamma1gq_1[ii1]*gamma1qg_1[ii1] + gamma1gq_2[ii2]*gamma1qg_2[ii2] + (gamma1qq_1[ii1] + gamma1qq_2[ii2])*(beta0 + gamma1qq_1[ii1] + gamma1qq_2[ii2]))*LQF)/ 2. - ((gamma1qq_1[ii1] + gamma1qq_2[ii2])*LQ*(2.*B1q - 2.*beta0 + A1q*LQ))/2.);
	    complex <double> H2st_qg_1   = (LQF*(2.*eC1qg_1*(gamma1qq_2[ii2] + gamma1gg_1[ii1]) + 2.*gamma2qg_1[ii1] + gamma1qg_1[ii1]*(2.*eC1qq_2 + 2.*eC1qq_1 + (beta0 + 2.*gamma1qq_2[ii2] + gamma1gg_1[ii1] + gamma1qq_1[ii1])*LF - (2.*B1q - beta0 + 2.*gamma1qq_2[ii2] + gamma1gg_1[ii1] + gamma1qq_1[ii1])*LQ - A1q*LQ2)))/2.;
	    complex <double> H2st_qg_2   = (LQF*(2.*eC1qg_2*(gamma1qq_1[ii1] + gamma1gg_2[ii2]) + 2.*gamma2qg_2[ii2] + gamma1qg_2[ii2]*(2.*eC1qq_1 + 2.*eC1qq_2 + (beta0 + 2.*gamma1qq_1[ii1] + gamma1gg_2[ii2] + gamma1qq_2[ii2])*LF - (2.*B1q - beta0 + 2.*gamma1qq_1[ii1] + gamma1gg_2[ii2] + gamma1qq_2[ii2])*LQ - A1q*LQ2)))/2.;
	    complex <double> H2st_qq_1   = ((2.*eC1qg_1*gamma1gq_1[ii1] + 2.*gamma2qqb_1[ii1] + gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQF)* LQF)/2.;
	    complex <double> H2st_qq_2   = ((2.*eC1qg_2*gamma1gq_2[ii2] + 2.*gamma2qqb_2[ii2] + gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQF)* LQF)/2.;
	    complex <double> H2st_qqp_1  = ((2.*eC1qg_1*gamma1gq_1[ii1] + 2.*gamma2qqbp_1[ii1] + gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQF)* LQF)/2.;
	    complex <double> H2st_qqp_2  = ((2.*eC1qg_2*gamma1gq_2[ii2] + 2.*gamma2qqbp_2[ii2] + gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQF)* LQF)/2.;
	    complex <double> H2st_qqbp_1 = ((2.*eC1qg_1*gamma1gq_1[ii1] + 2.*gamma2qqp_1[ii1] + gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQF)* LQF)/2.;
	    complex <double> H2st_qqbp_2 = ((2.*eC1qg_2*gamma1gq_2[ii2] + 2.*gamma2qqp_2[ii2] + gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQF)* LQF)/2.;
	    complex <double> H2st_gg     = (eC1qg_2*gamma1qg_1[ii1] + gamma1qg_2[ii2]*(eC1qg_1 + gamma1qg_1[ii1]*LQF))*LQF;

	    //Mixed LR*LQF variations
	    H2st_qqb  += - beta0*H1st_qqb*LR;
	    H2st_qg_1 += - beta0*H1st_qg_1*LR;
	    H2st_qg_2 += - beta0*H1st_qg_2*LR;
	    
	    qqb[idx]    += as2*H2st_qqb;
	    qg_1[idx]   += as2*H2st_qg_1;
	    qg_2[idx]   += as2*H2st_qg_2;
	    qq_1[idx]   += as2*H2st_qq_1;
	    qq_2[idx]   += as2*H2st_qq_2;
	    qqp_1[idx]  += as2*H2st_qqp_1;
	    qqp_2[idx]  += as2*H2st_qqp_2;
	    qqbp_1[idx] += as2*H2st_qqbp_1;
	    qqbp_2[idx] += as2*H2st_qqbp_2;
	    gg[idx]     += as2*H2st_gg;
	    
	    if (opts.order_hcoef == 2)
	      continue;

	    double nf = 5.;
	    complex <double> H3st_qqb    = eC1qg_2*eC1qq_1*gamma1gq_2[ii2]*LF + eC2qg_2*gamma1gq_2[ii2]*LF + eC1qq_2*eC1qq_1*gamma1qq_2[ii2]*LF + eC2qq_2*gamma1qq_2[ii2]*LF + eC2qq_1*gamma1qq_2[ii2]*LF + eC1qq_2*eC1qg_1*gamma1gq_1[ii1]*LF + eC2qg_1*gamma1gq_1[ii1]*LF + eC1qq_2*eC1qq_1*gamma1qq_1[ii1]*LF + eC2qq_2*gamma1qq_1[ii1]*LF + eC2qq_1*gamma1qq_1[ii1]*LF + eC1qg_1*gamma2gq_1[ii1]*LF + eC1qg_2*gamma2gq_2[ii2]*LF + eC1qq_2*gamma2qq_2[ii2]*LF + eC1qq_1*gamma2qq_2[ii2]*LF + eC1qq_2*gamma2qq_1[ii1]*LF + eC1qq_1*gamma2qq_1[ii1]*LF + gamma3qq_2[ii2]*LF + gamma3qq_1[ii1]*LF + (beta0*eC1qg_2*gamma1gq_2[ii2]*LF2)/2. + (eC1qg_2*gamma1gg_2[ii2]*gamma1gq_2[ii2]*LF2)/2. + (eC1qq_2*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2)/ 2. + (eC1qq_1*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2)/2. + (beta1*gamma1qq_2[ii2]*LF2)/2. + (beta0*eC1qq_2*gamma1qq_2[ii2]*LF2)/2. + (beta0*eC1qq_1*gamma1qq_2[ii2]*LF2)/2. + (eC1qg_2*gamma1gq_2[ii2]*gamma1qq_2[ii2]*LF2)/2. + (eC1qq_2*pow(gamma1qq_2[ii2],2)*LF2)/2. + (eC1qq_1*pow(gamma1qq_2[ii2],2)*LF2)/2. + (beta0*eC1qg_1*gamma1gq_1[ii1]*LF2)/2. + eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LF2 + (eC1qq_2*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2)/2. + (beta1*gamma1qq_1[ii1]*LF2)/2. + (beta0*eC1qq_2*gamma1qq_1[ii1]*LF2)/2. + (beta0*eC1qq_1*gamma1qq_1[ii1]*LF2)/2. + eC1qg_2*gamma1gq_2[ii2]*gamma1qq_1[ii1]*LF2 + eC1qq_2*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LF2 + eC1qq_1*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LF2 + (eC1qq_2*pow(gamma1qq_1[ii1],2)*LF2)/2. + gamma1qg_2[ii2]*gamma2gq_2[ii2]*LF2 + beta0*gamma2qq_2[ii2]*LF2 + gamma1qq_2[ii2]*gamma2qq_2[ii2]*LF2 + gamma1qq_1[ii1]*gamma2qq_2[ii2]*LF2 + gamma1qg_1[ii1]*gamma2gq_1[ii1]*LF2 + beta0*gamma2qq_1[ii1]*LF2 + gamma1qq_2[ii2]*gamma2qq_1[ii1]*LF2 + gamma1qq_1[ii1]*gamma2qq_1[ii1]*LF2 + (beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF3)/2. + (pow(beta0,2)*gamma1qq_2[ii2]*LF3)/3. + (beta0*pow(gamma1qq_2[ii2],2)*LF3)/2. + (beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF3)/2. + (gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF3)/2. + (pow(beta0,2)*gamma1qq_1[ii1]*LF3)/3. + (gamma1gq_2[ii2]*gamma1qg_2[ii2]*gamma1qq_1[ii1]*LF3)/2. + beta0*gamma1qq_2[ii2]*gamma1qq_1[ii1]* LF3 + (pow(gamma1qq_2[ii2],2)*gamma1qq_1[ii1]*LF3)/2. + (beta0*pow(gamma1qq_1[ii1],2)*LF3)/2. + (gamma1qq_2[ii2]*pow(gamma1qq_1[ii1],2)*LF3)/2. + ((gamma1gg_1[ii1]*gamma1gq_1[ii1]*gamma1qg_1[ii1] + 2.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*gamma1qq_1[ii1] + pow(gamma1qq_1[ii1],3))*LF3)/6. + ((gamma1gg_2[ii2]*gamma1gq_2[ii2]*gamma1qg_2[ii2] + 2.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*gamma1qq_2[ii2] + pow(gamma1qq_2[ii2],3))*LQF3)/6. - eC1qg_2*eC1qq_1*gamma1gq_2[ii2]*LQ - eC2qg_2*gamma1gq_2[ii2]*LQ - eC1qq_2*eC1qq_1*gamma1qq_2[ii2]*LQ - eC2qq_2*gamma1qq_2[ii2]*LQ - eC2qq_1*gamma1qq_2[ii2]*LQ - eC1qq_2*eC1qg_1*gamma1gq_1[ii1]*LQ - eC2qg_1*gamma1gq_1[ii1]*LQ - eC1qq_2*eC1qq_1*gamma1qq_1[ii1]*LQ - eC2qq_2*gamma1qq_1[ii1]*LQ - eC2qq_1*gamma1qq_1[ii1]*LQ - eC1qg_1*gamma2gq_1[ii1]*LQ - eC1qg_2*gamma2gq_2[ii2]*LQ - eC1qq_2*gamma2qq_2[ii2]*LQ - eC1qq_1*gamma2qq_2[ii2]*LQ - eC1qq_2*gamma2qq_1[ii1]*LQ - eC1qq_1*gamma2qq_1[ii1]*LQ - gamma3qq_2[ii2]*LQ - gamma3qq_1[ii1]*LQ - B1q*eC1qg_2*gamma1gq_2[ii2]*LF*LQ + beta0*eC1qg_2*gamma1gq_2[ii2]*LF*LQ - eC1qg_2*gamma1gg_2[ii2]*gamma1gq_2[ii2]*LF*LQ - eC1qq_2*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ - eC1qq_1*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ - B2q*gamma1qq_2[ii2]*LF*LQ - B1q*eC1qq_2*gamma1qq_2[ii2]*LF*LQ + beta0*eC1qq_2*gamma1qq_2[ii2]*LF*LQ - B1q*eC1qq_1*gamma1qq_2[ii2]*LF*LQ + beta0*eC1qq_1*gamma1qq_2[ii2]*LF*LQ - eC1qg_2*gamma1gq_2[ii2]*gamma1qq_2[ii2]*LF*LQ - eC1qq_2*pow(gamma1qq_2[ii2],2)*LF*LQ - eC1qq_1*pow(gamma1qq_2[ii2],2)*LF*LQ - B1q*eC1qg_1*gamma1gq_1[ii1]*LF*LQ + beta0*eC1qg_1*gamma1gq_1[ii1]*LF*LQ - 2.*eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LF*LQ - eC1qg_1*gamma1gg_1[ii1]*gamma1gq_1[ii1]*LF*LQ - eC1qq_2*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ - eC1qq_1*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ - B2q*gamma1qq_1[ii1]*LF*LQ - B1q*eC1qq_2*gamma1qq_1[ii1]*LF*LQ + beta0*eC1qq_2*gamma1qq_1[ii1]*LF*LQ - B1q*eC1qq_1*gamma1qq_1[ii1]*LF*LQ + beta0*eC1qq_1*gamma1qq_1[ii1]*LF*LQ - 2.*eC1qg_2*gamma1gq_2[ii2]*gamma1qq_1[ii1]*LF*LQ - 2.*eC1qq_2*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LF*LQ - 2.*eC1qq_1*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LF* LQ - eC1qg_1*gamma1gq_1[ii1]*gamma1qq_1[ii1]*LF*LQ - eC1qq_2*pow(gamma1qq_1[ii1],2)*LF*LQ - eC1qq_1*pow(gamma1qq_1[ii1],2)*LF*LQ - gamma1qg_2[ii2]*gamma2gq_2[ii2]*LF*LQ - gamma1gq_2[ii2]*gamma2qg_2[ii2]*LF*LQ - B1q*gamma2qq_2[ii2]*LF*LQ - 2.*gamma1qq_2[ii2]*gamma2qq_2[ii2]*LF*LQ - 2.*gamma1qq_1[ii1]*gamma2qq_2[ii2]*LF*LQ - gamma1qg_1[ii1]*gamma2gq_1[ii1]*LF*LQ - gamma1gq_1[ii1]*gamma2qg_1[ii1]*LF*LQ - B1q*gamma2qq_1[ii1]*LF*LQ - 2.*gamma1qq_2[ii2]*gamma2qq_1[ii1]*LF*LQ - 2.*gamma1qq_1[ii1]*gamma2qq_1[ii1]*LF*LQ - (B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2*LQ)/2. - (beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2*LQ)/2. - (B1q*beta0*gamma1qq_2[ii2]*LF2*LQ)/2. - (B1q*pow(gamma1qq_2[ii2],2)*LF2*LQ)/2. - (beta0*pow(gamma1qq_2[ii2],2)*LF2*LQ)/2. - (B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2*LQ)/2. - (beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2*LQ)/ 2. - (3.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2*LQ)/2. - (B1q*beta0*gamma1qq_1[ii1]*LF2*LQ)/2. - (3.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*gamma1qq_1[ii1]*LF2*LQ)/2. - B1q*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LF2*LQ - beta0*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LF2*LQ - (3.*pow(gamma1qq_2[ii2],2)*gamma1qq_1[ii1]*LF2*LQ)/2. - (B1q*pow(gamma1qq_1[ii1],2)*LF2*LQ)/2. - (beta0*pow(gamma1qq_1[ii1],2)*LF2*LQ)/2. - (3.*gamma1qq_2[ii2]*pow(gamma1qq_1[ii1],2)*LF2*LQ)/2. - ((gamma1gg_1[ii1]*gamma1gq_1[ii1]*gamma1qg_1[ii1] + 2.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*gamma1qq_1[ii1] + pow(gamma1qq_1[ii1],3))*LF2*LQ)/2. + B1q*eC1qg_2*gamma1gq_2[ii2]*LQ2 - (3.*beta0*eC1qg_2*gamma1gq_2[ii2]*LQ2)/2. + (eC1qg_2*gamma1gg_2[ii2]*gamma1gq_2[ii2]*LQ2)/2. + (eC1qq_2*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ2)/2. + (eC1qq_1*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ2)/ 2. + B2q*gamma1qq_2[ii2]*LQ2 - (beta1*gamma1qq_2[ii2]*LQ2)/2. + B1q*eC1qq_2*gamma1qq_2[ii2]*LQ2 - (3.*beta0*eC1qq_2*gamma1qq_2[ii2]*LQ2)/2. + B1q*eC1qq_1*gamma1qq_2[ii2]*LQ2 - (3.*beta0*eC1qq_1*gamma1qq_2[ii2]*LQ2)/2. + (eC1qg_2*gamma1gq_2[ii2]*gamma1qq_2[ii2]*LQ2)/2. + (eC1qq_2*pow(gamma1qq_2[ii2],2)*LQ2)/2. + (eC1qq_1*pow(gamma1qq_2[ii2],2)*LQ2)/2. + B1q*eC1qg_1*gamma1gq_1[ii1]*LQ2 - (3.*beta0*eC1qg_1*gamma1gq_1[ii1]*LQ2)/2. + eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LQ2 + (eC1qq_2*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ2)/2. + B2q*gamma1qq_1[ii1]*LQ2 - (beta1*gamma1qq_1[ii1]*LQ2)/2. + B1q*eC1qq_2*gamma1qq_1[ii1]*LQ2 - (3.*beta0*eC1qq_2*gamma1qq_1[ii1]*LQ2)/2. + B1q*eC1qq_1*gamma1qq_1[ii1]*LQ2 - (3.*beta0*eC1qq_1*gamma1qq_1[ii1]*LQ2)/2. + eC1qg_2*gamma1gq_2[ii2]*gamma1qq_1[ii1]*LQ2 + eC1qq_2*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LQ2 + eC1qq_1*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LQ2 + (eC1qq_2*pow(gamma1qq_1[ii1],2)*LQ2)/2. + gamma1qg_2[ii2]*gamma2gq_2[ii2]*LQ2 + B1q*gamma2qq_2[ii2]*LQ2 - beta0*gamma2qq_2[ii2]*LQ2 + gamma1qq_2[ii2]*gamma2qq_2[ii2]*LQ2 + gamma1qq_1[ii1]*gamma2qq_2[ii2]*LQ2 + gamma1qg_1[ii1]*gamma2gq_1[ii1]*LQ2 + B1q*gamma2qq_1[ii1]*LQ2 - beta0*gamma2qq_1[ii1]*LQ2 + gamma1qq_2[ii2]*gamma2qq_1[ii1]*LQ2 + gamma1qq_1[ii1]*gamma2qq_1[ii1]*LQ2 - (A1q*eC1qg_2*gamma1gq_2[ii2]*LF*LQ2)/2. + B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ2 - (beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ2)/2. - (A2q*gamma1qq_2[ii2]*LF*LQ2)/2. + (pow(B1q,2)*gamma1qq_2[ii2]*LF*LQ2)/2. - (B1q*beta0*gamma1qq_2[ii2]*LF*LQ2)/2. - (A1q*eC1qq_2*gamma1qq_2[ii2]*LF*LQ2)/2. - (A1q*eC1qq_1*gamma1qq_2[ii2]*LF*LQ2)/2. + B1q*pow(gamma1qq_2[ii2],2)*LF*LQ2 - (beta0*pow(gamma1qq_2[ii2],2)*LF*LQ2)/2. - (A1q*eC1qg_1*gamma1gq_1[ii1]*LF*LQ2)/2. + B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ2 - (beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ2)/2. + (3.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ2)/2. - (A2q*gamma1qq_1[ii1]*LF*LQ2)/2. + (pow(B1q,2)*gamma1qq_1[ii1]*LF*LQ2)/2. - (B1q*beta0*gamma1qq_1[ii1]*LF*LQ2)/2. - (A1q*eC1qq_2*gamma1qq_1[ii1]*LF*LQ2)/2. - (A1q*eC1qq_1*gamma1qq_1[ii1]*LF*LQ2)/2. + (3.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*gamma1qq_1[ii1]*LF*LQ2)/2. + 2.*B1q*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LF*LQ2 - beta0*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LF*LQ2 + (3.*pow(gamma1qq_2[ii2],2)*gamma1qq_1[ii1]*LF*LQ2)/2. + B1q*pow(gamma1qq_1[ii1],2)*LF*LQ2 - (beta0*pow(gamma1qq_1[ii1],2)*LF*LQ2)/2. + (3.*gamma1qq_2[ii2]*pow(gamma1qq_1[ii1],2)*LF*LQ2)/2. + ((gamma1gg_1[ii1]*gamma1gq_1[ii1]*gamma1qg_1[ii1] + 2.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*gamma1qq_1[ii1] + pow(gamma1qq_1[ii1],3))*LF*LQ2)/2. - (A1q*gamma2qq_2[ii2]*LF*LQ2)/2. - (A1q*gamma2qq_1[ii1]*LF*LQ2)/2. - (A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2*LQ2)/4. - (A1q*beta0*gamma1qq_2[ii2]*LF2*LQ2)/4. - (A1q*pow(gamma1qq_2[ii2],2)*LF2*LQ2)/4. - (A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2*LQ2)/4. - (A1q*beta0*gamma1qq_1[ii1]*LF2*LQ2)/4. - (A1q*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LF2*LQ2)/2. - (A1q*pow(gamma1qq_1[ii1],2)*LF2*LQ2)/4. + (A1q*eC1qg_2*gamma1gq_2[ii2]*LQ3)/2. - (B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ3)/2. + (beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ3)/2. + (A2q*gamma1qq_2[ii2]*LQ3)/2. - (pow(B1q,2)*gamma1qq_2[ii2]*LQ3)/2. + B1q*beta0*gamma1qq_2[ii2]*LQ3 - (pow(beta0,2)*gamma1qq_2[ii2]*LQ3)/3. + (A1q*eC1qq_2*gamma1qq_2[ii2]*LQ3)/2. + (A1q*eC1qq_1*gamma1qq_2[ii2]*LQ3)/2. - (B1q*pow(gamma1qq_2[ii2],2)*LQ3)/2. + (beta0*pow(gamma1qq_2[ii2],2)*LQ3)/2. + (A1q*eC1qg_1*gamma1gq_1[ii1]*LQ3)/2. - (B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ3)/2. + (beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ3)/2. - (gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ3)/2. + (A2q*gamma1qq_1[ii1]*LQ3)/2. - (pow(B1q,2)*gamma1qq_1[ii1]*LQ3)/2. + B1q*beta0*gamma1qq_1[ii1]*LQ3 - (pow(beta0,2)*gamma1qq_1[ii1]*LQ3)/3. + (A1q*eC1qq_2*gamma1qq_1[ii1]*LQ3)/2. + (A1q*eC1qq_1*gamma1qq_1[ii1]*LQ3)/2. - (gamma1gq_2[ii2]*gamma1qg_2[ii2]*gamma1qq_1[ii1]*LQ3)/2. - B1q*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LQ3 + beta0*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LQ3 - (pow(gamma1qq_2[ii2],2)*gamma1qq_1[ii1]*LQ3)/2. - (B1q*pow(gamma1qq_1[ii1],2)*LQ3)/2. + (beta0*pow(gamma1qq_1[ii1],2)*LQ3)/2. - (gamma1qq_2[ii2]*pow(gamma1qq_1[ii1],2)*LQ3)/2. - ((gamma1gg_1[ii1]*gamma1gq_1[ii1]*gamma1qg_1[ii1] + 2.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*gamma1qq_1[ii1] + pow(gamma1qq_1[ii1],3))*LQ3)/6. + (A1q*gamma2qq_2[ii2]*LQ3)/2. + (A1q*gamma2qq_1[ii1]*LQ3)/2. + (A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ3)/2. + (A1q*B1q*gamma1qq_2[ii2]*LF*LQ3)/2. - (A1q*beta0*gamma1qq_2[ii2]*LF*LQ3)/3. + (A1q*pow(gamma1qq_2[ii2],2)*LF*LQ3)/2. + (A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ3)/2. + (A1q*B1q*gamma1qq_1[ii1]*LF*LQ3)/2. - (A1q*beta0*gamma1qq_1[ii1]*LF*LQ3)/3. + A1q*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LF*LQ3 + (A1q*pow(gamma1qq_1[ii1],2)*LF*LQ3)/2. - (A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ4)/4. - (A1q*B1q*gamma1qq_2[ii2]*LQ4)/2. + (7.*A1q*beta0*gamma1qq_2[ii2]*LQ4)/12. - (A1q*pow(gamma1qq_2[ii2],2)*LQ4)/4. - (A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ4)/4. - (A1q*B1q*gamma1qq_1[ii1]*LQ4)/2. + (7.*A1q*beta0*gamma1qq_1[ii1]*LQ4)/12. - (A1q*gamma1qq_2[ii2]*gamma1qq_1[ii1]*LQ4)/2. - (A1q*pow(gamma1qq_1[ii1],2)*LQ4)/4. + (pow(A1q,2)*gamma1qq_2[ii2]*LF*LQ4)/8. + (pow(A1q,2)*gamma1qq_1[ii1]*LF*LQ4)/8. - (pow(A1q,2)*gamma1qq_2[ii2]*LQ5)/8. - (pow(A1q,2)*gamma1qq_1[ii1]*LQ5)/8. + ((eC1qg_1*gamma1gq_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1]) + eC1qq_1*(gamma1gq_1[ii1]*gamma1qg_1[ii1] + pow(gamma1qq_1[ii1],2)))*(LF2 + LQ2))/2.;
	    complex <double> H3st_qg_1   = eC2qg_1*gamma1qq_2[ii2]*LF + eC1qg_1*(eC1qg_2*gamma1gq_2[ii2] + eC1qq_2*gamma1qq_2[ii2])*LF + eC2qq_2*gamma1qg_1[ii1]*LF + eC1qq_2*(eC1qg_1*gamma1gg_1[ii1] + eC1qq_1*gamma1qg_1[ii1])*LF + eC1qg_1*gamma2qq_2[ii2]*LF + eC1qq_2*gamma2qg_1[ii1]*LF + (eC1qg_1*gamma2gg_1[ii1] + eC1qq_1*gamma2qg_1[ii1])*LF + gamma3qg_1[ii1]*LF + (beta0*eC1qg_1*gamma1qq_2[ii2]*LF2)/2. + (eC1qg_1*(gamma1gq_2[ii2]*gamma1qg_2[ii2] + pow(gamma1qq_2[ii2],2))*LF2)/2. + (beta1*gamma1qg_1[ii1]*LF2)/2. + (beta0*eC1qq_2*gamma1qg_1[ii1]*LF2)/2. + (eC1qg_2*gamma1gq_2[ii2] + eC1qq_2*gamma1qq_2[ii2])*gamma1qg_1[ii1]*LF2 + (beta0*(eC1qg_1*gamma1gg_1[ii1] + eC1qq_1*gamma1qg_1[ii1])*LF2)/2. + gamma1qq_2[ii2]*(eC1qg_1*gamma1gg_1[ii1] + eC1qq_1*gamma1qg_1[ii1])*LF2 + (eC1qq_2*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LF2)/2. + gamma1qg_1[ii1]*gamma2qq_2[ii2]*LF2 + beta0*gamma2qg_1[ii1]*LF2 + gamma1qq_2[ii2]*gamma2qg_1[ii1]*LF2 + (gamma1qg_1[ii1]*gamma2gg_1[ii1] + gamma1qq_1[ii1]*gamma2qg_1[ii1])* LF2 + (pow(beta0,2)*gamma1qg_1[ii1]*LF3)/3. + beta0*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LF3 + ((gamma1gq_2[ii2]*gamma1qg_2[ii2] + pow(gamma1qq_2[ii2],2))*gamma1qg_1[ii1]*LF3)/2. + (beta0*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LF3)/2. + (gamma1qq_2[ii2]*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LF3)/2. - eC2qg_1*gamma1qq_2[ii2]*LQ - eC1qg_1*(eC1qg_2*gamma1gq_2[ii2] + eC1qq_2*gamma1qq_2[ii2])*LQ - eC2qq_2*gamma1qg_1[ii1]*LQ - eC1qq_2*(eC1qg_1*gamma1gg_1[ii1] + eC1qq_1*gamma1qg_1[ii1])*LQ - eC1qg_1*gamma2qq_2[ii2]*LQ - eC1qq_2*gamma2qg_1[ii1]*LQ - (eC1qg_1*gamma2gg_1[ii1] + eC1qq_1*gamma2qg_1[ii1])*LQ - gamma3qg_1[ii1]*LQ - B1q*eC1qg_1*gamma1qq_2[ii2]*LF*LQ + beta0*eC1qg_1*gamma1qq_2[ii2]*LF*LQ - eC1qg_1*(gamma1gq_2[ii2]*gamma1qg_2[ii2] + pow(gamma1qq_2[ii2],2))*LF*LQ - B2q*gamma1qg_1[ii1]*LF*LQ - B1q*eC1qq_2*gamma1qg_1[ii1]*LF*LQ + beta0*eC1qq_2*gamma1qg_1[ii1]*LF*LQ - 2.*(eC1qg_2*gamma1gq_2[ii2] + eC1qq_2*gamma1qq_2[ii2])*gamma1qg_1[ii1]*LF*LQ - B1q*(eC1qg_1*gamma1gg_1[ii1] + eC1qq_1*gamma1qg_1[ii1])*LF*LQ + beta0*(eC1qg_1*gamma1gg_1[ii1] + eC1qq_1*gamma1qg_1[ii1])*LF*LQ - 2.*gamma1qq_2[ii2]*(eC1qg_1*gamma1gg_1[ii1] + eC1qq_1*gamma1qg_1[ii1])*LF*LQ - eC1qq_2*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LF*LQ - 2.*gamma1qg_1[ii1]*gamma2qq_2[ii2]*LF*LQ - B1q*gamma2qg_1[ii1]*LF*LQ - 2.*gamma1qq_2[ii2]*gamma2qg_1[ii1]*LF*LQ - (gamma1qg_1[ii1]*gamma2gg_1[ii1] + gamma1qq_1[ii1]*gamma2qg_1[ii1])*LF*LQ - (B1q*beta0*gamma1qg_1[ii1]*LF2*LQ)/2. - B1q*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LF2*LQ - beta0*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LF2*LQ - (3.*(gamma1gq_2[ii2]*gamma1qg_2[ii2] + pow(gamma1qq_2[ii2],2))*gamma1qg_1[ii1]*LF2*LQ)/2. - (B1q*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LF2*LQ)/2. - (beta0*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LF2*LQ)/2. - (3.*gamma1qq_2[ii2]*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LF2*LQ)/2. + B1q*eC1qg_1*gamma1qq_2[ii2]*LQ2 - (3.*beta0*eC1qg_1*gamma1qq_2[ii2]*LQ2)/2. + (eC1qg_1*(gamma1gq_2[ii2]*gamma1qg_2[ii2] + pow(gamma1qq_2[ii2],2))*LQ2)/2. + B2q*gamma1qg_1[ii1]*LQ2 - (beta1*gamma1qg_1[ii1]*LQ2)/2. + B1q*eC1qq_2*gamma1qg_1[ii1]*LQ2 - (3.*beta0*eC1qq_2*gamma1qg_1[ii1]*LQ2)/2. + (eC1qg_2*gamma1gq_2[ii2] + eC1qq_2*gamma1qq_2[ii2])* gamma1qg_1[ii1]*LQ2 + B1q*(eC1qg_1*gamma1gg_1[ii1] + eC1qq_1*gamma1qg_1[ii1])*LQ2 - (3.*beta0*(eC1qg_1*gamma1gg_1[ii1] + eC1qq_1*gamma1qg_1[ii1])*LQ2)/2. + gamma1qq_2[ii2]*(eC1qg_1*gamma1gg_1[ii1] + eC1qq_1*gamma1qg_1[ii1])*LQ2 + (eC1qq_2*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LQ2)/2. + gamma1qg_1[ii1]*gamma2qq_2[ii2]*LQ2 + B1q*gamma2qg_1[ii1]*LQ2 - beta0*gamma2qg_1[ii1]*LQ2 + gamma1qq_2[ii2]*gamma2qg_1[ii1]*LQ2 + (gamma1qg_1[ii1]*gamma2gg_1[ii1] + gamma1qq_1[ii1]*gamma2qg_1[ii1])* LQ2 - (A1q*eC1qg_1*gamma1qq_2[ii2]*LF*LQ2)/2. - (A2q*gamma1qg_1[ii1]*LF*LQ2)/2. + (pow(B1q,2)*gamma1qg_1[ii1]*LF*LQ2)/2. - (B1q*beta0*gamma1qg_1[ii1]*LF*LQ2)/2. - (A1q*eC1qq_2*gamma1qg_1[ii1]*LF*LQ2)/2. + 2.*B1q*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LF*LQ2 - beta0*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LF*LQ2 + (3.*(gamma1gq_2[ii2]*gamma1qg_2[ii2] + pow(gamma1qq_2[ii2],2))*gamma1qg_1[ii1]*LF*LQ2)/2. - (A1q*(eC1qg_1*gamma1gg_1[ii1] + eC1qq_1*gamma1qg_1[ii1])*LF*LQ2)/2. + B1q*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LF*LQ2 - (beta0*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LF*LQ2)/2. + (3.*gamma1qq_2[ii2]*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LF*LQ2)/2. - (A1q*gamma2qg_1[ii1]*LF*LQ2)/2. - (A1q*beta0*gamma1qg_1[ii1]*LF2*LQ2)/4. - (A1q*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LF2*LQ2)/2. - (A1q*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LF2*LQ2)/4. + (A1q*eC1qg_1*gamma1qq_2[ii2]*LQ3)/2. + (A2q*gamma1qg_1[ii1]*LQ3)/2. - (pow(B1q,2)*gamma1qg_1[ii1]*LQ3)/2. + B1q*beta0*gamma1qg_1[ii1]*LQ3 - (pow(beta0,2)*gamma1qg_1[ii1]*LQ3)/3. + (A1q*eC1qq_2*gamma1qg_1[ii1]*LQ3)/2. - B1q*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LQ3 + beta0*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LQ3 - ((gamma1gq_2[ii2]*gamma1qg_2[ii2] + pow(gamma1qq_2[ii2],2))*gamma1qg_1[ii1]*LQ3)/2. + (A1q*(eC1qg_1*gamma1gg_1[ii1] + eC1qq_1*gamma1qg_1[ii1])*LQ3)/2. - (B1q*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LQ3)/2. + (beta0*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LQ3)/2. - (gamma1qq_2[ii2]*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LQ3)/2. + (A1q*gamma2qg_1[ii1]*LQ3)/2. + (A1q*B1q*gamma1qg_1[ii1]*LF*LQ3)/2. - (A1q*beta0*gamma1qg_1[ii1]*LF*LQ3)/3. + A1q*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LF*LQ3 + (A1q*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LF*LQ3)/2. - (A1q*B1q*gamma1qg_1[ii1]*LQ4)/2. + (7.*A1q*beta0*gamma1qg_1[ii1]*LQ4)/12. - (A1q*gamma1qq_2[ii2]*gamma1qg_1[ii1]*LQ4)/2. - (A1q*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1])*LQ4)/4. + (pow(A1q,2)*gamma1qg_1[ii1]*LF*LQ4)/8. - (pow(A1q,2)*gamma1qg_1[ii1]*LQ5)/8. + LF*(eC2qg_1*gamma1gg_1[ii1] + gamma1qg_1[ii1]*(eC2qq_1 + eC2qqb_1 + (eC2qqbp_1 + eC2qqp_1)*(-1. + nf))) - LQ*(eC2qg_1*gamma1gg_1[ii1] + gamma1qg_1[ii1]*(eC2qq_1 + eC2qqb_1 + (eC2qqbp_1 + eC2qqp_1)*(-1. + nf))) - LF*LQ*(gamma1gg_1[ii1]*gamma2qg_1[ii1] + gamma1qg_1[ii1]*(gamma2qq_1[ii1] + gamma2qqb_1[ii1] + (gamma2qqbp_1[ii1] + gamma2qqp_1[ii1])*(-1. + nf))) + (gamma1qg_1[ii1]*LF3*(pow(gamma1gg_1[ii1],2) + gamma1gg_1[ii1]*gamma1qq_1[ii1] + pow(gamma1qq_1[ii1],2) + 2.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*nf))/6. - (gamma1qg_1[ii1]*LF2*LQ*(pow(gamma1gg_1[ii1],2) + gamma1gg_1[ii1]*gamma1qq_1[ii1] + pow(gamma1qq_1[ii1],2) + 2.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*nf))/2. + (gamma1qg_1[ii1]*LF*LQ2*(pow(gamma1gg_1[ii1],2) + gamma1gg_1[ii1]*gamma1qq_1[ii1] + pow(gamma1qq_1[ii1],2) + 2.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*nf))/2. - (gamma1qg_1[ii1]*LQ3*(pow(gamma1gg_1[ii1],2) + gamma1gg_1[ii1]*gamma1qq_1[ii1] + pow(gamma1qq_1[ii1],2) + 2.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*nf))/6. + (LF2*(eC1qq_1*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1]) + eC1qg_1*(pow(gamma1gg_1[ii1],2) + 2.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*nf)))/2. - LF*LQ*(eC1qq_1*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1]) + eC1qg_1*(pow(gamma1gg_1[ii1],2) + 2.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*nf)) + (LQ2*(eC1qq_1*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + gamma1qq_1[ii1]) + eC1qg_1*(pow(gamma1gg_1[ii1],2) + 2.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*nf)))/2.;
	    complex <double> H3st_qg_2   = eC2qg_2*gamma1qq_1[ii1]*LF + eC1qg_2*(eC1qg_1*gamma1gq_1[ii1] + eC1qq_1*gamma1qq_1[ii1])*LF + eC2qq_1*gamma1qg_2[ii2]*LF + eC1qq_1*(eC1qg_2*gamma1gg_2[ii2] + eC1qq_2*gamma1qg_2[ii2])*LF + eC1qg_2*gamma2qq_1[ii1]*LF + eC1qq_1*gamma2qg_2[ii2]*LF + (eC1qg_2*gamma2gg_2[ii2] + eC1qq_2*gamma2qg_2[ii2])*LF + gamma3qg_2[ii2]*LF + (beta0*eC1qg_2*gamma1qq_1[ii1]*LF2)/2. + (eC1qg_2*(gamma1gq_1[ii1]*gamma1qg_1[ii1] + pow(gamma1qq_1[ii1],2))*LF2)/2. + (beta1*gamma1qg_2[ii2]*LF2)/2. + (beta0*eC1qq_1*gamma1qg_2[ii2]*LF2)/2. + (eC1qg_1*gamma1gq_1[ii1] + eC1qq_1*gamma1qq_1[ii1])*gamma1qg_2[ii2]*LF2 + (beta0*(eC1qg_2*gamma1gg_2[ii2] + eC1qq_2*gamma1qg_2[ii2])*LF2)/2. + gamma1qq_1[ii1]*(eC1qg_2*gamma1gg_2[ii2] + eC1qq_2*gamma1qg_2[ii2])*LF2 + (eC1qq_1*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LF2)/2. + gamma1qg_2[ii2]*gamma2qq_1[ii1]*LF2 + beta0*gamma2qg_2[ii2]*LF2 + gamma1qq_1[ii1]*gamma2qg_2[ii2]*LF2 + (gamma1qg_2[ii2]*gamma2gg_2[ii2] + gamma1qq_2[ii2]*gamma2qg_2[ii2])* LF2 + (pow(beta0,2)*gamma1qg_2[ii2]*LF3)/3. + beta0*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LF3 + ((gamma1gq_1[ii1]*gamma1qg_1[ii1] + pow(gamma1qq_1[ii1],2))*gamma1qg_2[ii2]*LF3)/2. + (beta0*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LF3)/2. + (gamma1qq_1[ii1]*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LF3)/2. - eC2qg_2*gamma1qq_1[ii1]*LQ - eC1qg_2*(eC1qg_1*gamma1gq_1[ii1] + eC1qq_1*gamma1qq_1[ii1])*LQ - eC2qq_1*gamma1qg_2[ii2]*LQ - eC1qq_1*(eC1qg_2*gamma1gg_2[ii2] + eC1qq_2*gamma1qg_2[ii2])*LQ - eC1qg_2*gamma2qq_1[ii1]*LQ - eC1qq_1*gamma2qg_2[ii2]*LQ - (eC1qg_2*gamma2gg_2[ii2] + eC1qq_2*gamma2qg_2[ii2])*LQ - gamma3qg_2[ii2]*LQ - B1q*eC1qg_2*gamma1qq_1[ii1]*LF*LQ + beta0*eC1qg_2*gamma1qq_1[ii1]*LF*LQ - eC1qg_2*(gamma1gq_1[ii1]*gamma1qg_1[ii1] + pow(gamma1qq_1[ii1],2))*LF*LQ - B2q*gamma1qg_2[ii2]*LF*LQ - B1q*eC1qq_1*gamma1qg_2[ii2]*LF*LQ + beta0*eC1qq_1*gamma1qg_2[ii2]*LF*LQ - 2.*(eC1qg_1*gamma1gq_1[ii1] + eC1qq_1*gamma1qq_1[ii1])*gamma1qg_2[ii2]*LF*LQ - B1q*(eC1qg_2*gamma1gg_2[ii2] + eC1qq_2*gamma1qg_2[ii2])*LF*LQ + beta0*(eC1qg_2*gamma1gg_2[ii2] + eC1qq_2*gamma1qg_2[ii2])*LF*LQ - 2.*gamma1qq_1[ii1]*(eC1qg_2*gamma1gg_2[ii2] + eC1qq_2*gamma1qg_2[ii2])*LF*LQ - eC1qq_1*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LF*LQ - 2.*gamma1qg_2[ii2]*gamma2qq_1[ii1]*LF*LQ - B1q*gamma2qg_2[ii2]*LF*LQ - 2.*gamma1qq_1[ii1]*gamma2qg_2[ii2]*LF*LQ - (gamma1qg_2[ii2]*gamma2gg_2[ii2] + gamma1qq_2[ii2]*gamma2qg_2[ii2])*LF*LQ - (B1q*beta0*gamma1qg_2[ii2]*LF2*LQ)/2. - B1q*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LF2*LQ - beta0*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LF2*LQ - (3.*(gamma1gq_1[ii1]*gamma1qg_1[ii1] + pow(gamma1qq_1[ii1],2))*gamma1qg_2[ii2]*LF2*LQ)/2. - (B1q*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LF2*LQ)/2. - (beta0*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LF2*LQ)/2. - (3.*gamma1qq_1[ii1]*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LF2*LQ)/2. + B1q*eC1qg_2*gamma1qq_1[ii1]*LQ2 - (3.*beta0*eC1qg_2*gamma1qq_1[ii1]*LQ2)/2. + (eC1qg_2*(gamma1gq_1[ii1]*gamma1qg_1[ii1] + pow(gamma1qq_1[ii1],2))*LQ2)/2. + B2q*gamma1qg_2[ii2]*LQ2 - (beta1*gamma1qg_2[ii2]*LQ2)/2. + B1q*eC1qq_1*gamma1qg_2[ii2]*LQ2 - (3.*beta0*eC1qq_1*gamma1qg_2[ii2]*LQ2)/2. + (eC1qg_1*gamma1gq_1[ii1] + eC1qq_1*gamma1qq_1[ii1])* gamma1qg_2[ii2]*LQ2 + B1q*(eC1qg_2*gamma1gg_2[ii2] + eC1qq_2*gamma1qg_2[ii2])*LQ2 - (3.*beta0*(eC1qg_2*gamma1gg_2[ii2] + eC1qq_2*gamma1qg_2[ii2])*LQ2)/2. + gamma1qq_1[ii1]*(eC1qg_2*gamma1gg_2[ii2] + eC1qq_2*gamma1qg_2[ii2])*LQ2 + (eC1qq_1*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LQ2)/2. + gamma1qg_2[ii2]*gamma2qq_1[ii1]*LQ2 + B1q*gamma2qg_2[ii2]*LQ2 - beta0*gamma2qg_2[ii2]*LQ2 + gamma1qq_1[ii1]*gamma2qg_2[ii2]*LQ2 + (gamma1qg_2[ii2]*gamma2gg_2[ii2] + gamma1qq_2[ii2]*gamma2qg_2[ii2])* LQ2 - (A1q*eC1qg_2*gamma1qq_1[ii1]*LF*LQ2)/2. - (A2q*gamma1qg_2[ii2]*LF*LQ2)/2. + (pow(B1q,2)*gamma1qg_2[ii2]*LF*LQ2)/2. - (B1q*beta0*gamma1qg_2[ii2]*LF*LQ2)/2. - (A1q*eC1qq_1*gamma1qg_2[ii2]*LF*LQ2)/2. + 2.*B1q*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LF*LQ2 - beta0*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LF*LQ2 + (3.*(gamma1gq_1[ii1]*gamma1qg_1[ii1] + pow(gamma1qq_1[ii1],2))*gamma1qg_2[ii2]*LF*LQ2)/2. - (A1q*(eC1qg_2*gamma1gg_2[ii2] + eC1qq_2*gamma1qg_2[ii2])*LF*LQ2)/2. + B1q*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LF*LQ2 - (beta0*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LF*LQ2)/2. + (3.*gamma1qq_1[ii1]*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LF*LQ2)/2. - (A1q*gamma2qg_2[ii2]*LF*LQ2)/2. - (A1q*beta0*gamma1qg_2[ii2]*LF2*LQ2)/4. - (A1q*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LF2*LQ2)/2. - (A1q*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LF2*LQ2)/4. + (A1q*eC1qg_2*gamma1qq_1[ii1]*LQ3)/2. + (A2q*gamma1qg_2[ii2]*LQ3)/2. - (pow(B1q,2)*gamma1qg_2[ii2]*LQ3)/2. + B1q*beta0*gamma1qg_2[ii2]*LQ3 - (pow(beta0,2)*gamma1qg_2[ii2]*LQ3)/3. + (A1q*eC1qq_1*gamma1qg_2[ii2]*LQ3)/2. - B1q*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LQ3 + beta0*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LQ3 - ((gamma1gq_1[ii1]*gamma1qg_1[ii1] + pow(gamma1qq_1[ii1],2))*gamma1qg_2[ii2]*LQ3)/2. + (A1q*(eC1qg_2*gamma1gg_2[ii2] + eC1qq_2*gamma1qg_2[ii2])*LQ3)/2. - (B1q*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LQ3)/2. + (beta0*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LQ3)/2. - (gamma1qq_1[ii1]*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LQ3)/2. + (A1q*gamma2qg_2[ii2]*LQ3)/2. + (A1q*B1q*gamma1qg_2[ii2]*LF*LQ3)/2. - (A1q*beta0*gamma1qg_2[ii2]*LF*LQ3)/3. + A1q*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LF*LQ3 + (A1q*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LF*LQ3)/2. - (A1q*B1q*gamma1qg_2[ii2]*LQ4)/2. + (7.*A1q*beta0*gamma1qg_2[ii2]*LQ4)/12. - (A1q*gamma1qq_1[ii1]*gamma1qg_2[ii2]*LQ4)/2. - (A1q*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2])*LQ4)/4. + (pow(A1q,2)*gamma1qg_2[ii2]*LF*LQ4)/8. - (pow(A1q,2)*gamma1qg_2[ii2]*LQ5)/8. + LF*(eC2qg_2*gamma1gg_2[ii2] + gamma1qg_2[ii2]*(eC2qq_2 + eC2qqb_2 + (eC2qqbp_2 + eC2qqp_2)*(-1. + nf))) - LQ*(eC2qg_2*gamma1gg_2[ii2] + gamma1qg_2[ii2]*(eC2qq_2 + eC2qqb_2 + (eC2qqbp_2 + eC2qqp_2)*(-1. + nf))) - LF*LQ*(gamma1gg_2[ii2]*gamma2qg_2[ii2] + gamma1qg_2[ii2]*(gamma2qq_2[ii2] + gamma2qqb_2[ii2] + (gamma2qqbp_2[ii2] + gamma2qqp_2[ii2])*(-1. + nf))) + (gamma1qg_2[ii2]*LF3*(pow(gamma1gg_2[ii2],2) + gamma1gg_2[ii2]*gamma1qq_2[ii2] + pow(gamma1qq_2[ii2],2) + 2.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*nf))/6. - (gamma1qg_2[ii2]*LF2*LQ*(pow(gamma1gg_2[ii2],2) + gamma1gg_2[ii2]*gamma1qq_2[ii2] + pow(gamma1qq_2[ii2],2) + 2.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*nf))/2. + (gamma1qg_2[ii2]*LF*LQ2*(pow(gamma1gg_2[ii2],2) + gamma1gg_2[ii2]*gamma1qq_2[ii2] + pow(gamma1qq_2[ii2],2) + 2.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*nf))/2. - (gamma1qg_2[ii2]*LQ3*(pow(gamma1gg_2[ii2],2) + gamma1gg_2[ii2]*gamma1qq_2[ii2] + pow(gamma1qq_2[ii2],2) + 2.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*nf))/6. + (LF2*(eC1qq_2*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2]) + eC1qg_2*(pow(gamma1gg_2[ii2],2) + 2.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*nf)))/2. - LF*LQ*(eC1qq_2*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2]) + eC1qg_2*(pow(gamma1gg_2[ii2],2) + 2.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*nf)) + (LQ2*(eC1qq_2*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + gamma1qq_2[ii2]) + eC1qg_2*(pow(gamma1gg_2[ii2],2) + 2.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*nf)))/2.;
	    complex <double> H3st_qq_1   = eC2qqb_1*gamma1qq_2[ii2]*LF + eC1qq_2*eC1qg_1*gamma1gq_1[ii1]*LF + (eC2qg_1*gamma1gq_1[ii1] + eC2qqb_1*gamma1qq_1[ii1])*LF + eC1qq_2*gamma2qqb_1[ii1]*LF + (eC1qg_1*gamma2gq_1[ii1] + eC1qq_1*gamma2qqb_1[ii1])*LF + gamma3qqb_1[ii1]*LF + (beta0*eC1qg_1*gamma1gq_1[ii1]*LF2)/2. + eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LF2 + (eC1qq_2*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2)/2. + (gamma1gq_1[ii1]*(eC1qq_1*gamma1qg_1[ii1] + eC1qg_1*(gamma1gg_1[ii1] + gamma1qq_1[ii1]))*LF2)/2. + beta0*gamma2qqb_1[ii1]*LF2 + gamma1qq_2[ii2]*gamma2qqb_1[ii1]*LF2 + (gamma1qg_1[ii1]*gamma2gq_1[ii1] + gamma1qq_1[ii1]*gamma2qqb_1[ii1])*LF2 + (beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF3)/2. + (gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF3)/2. + (gamma1gq_1[ii1]*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + 2.*gamma1qq_1[ii1])*LF3)/6. - eC2qqb_1*gamma1qq_2[ii2]*LQ - eC1qq_2*eC1qg_1*gamma1gq_1[ii1]*LQ - (eC2qg_1*gamma1gq_1[ii1] + eC2qqb_1*gamma1qq_1[ii1])*LQ - eC1qq_2*gamma2qqb_1[ii1]*LQ - (eC1qg_1*gamma2gq_1[ii1] + eC1qq_1*gamma2qqb_1[ii1])*LQ - gamma3qqb_1[ii1]*LQ - B1q*eC1qg_1*gamma1gq_1[ii1]*LF*LQ + beta0*eC1qg_1*gamma1gq_1[ii1]*LF*LQ - 2.*eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LF*LQ - eC1qq_2*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ - gamma1gq_1[ii1]*(eC1qq_1*gamma1qg_1[ii1] + eC1qg_1*(gamma1gg_1[ii1] + gamma1qq_1[ii1]))*LF*LQ - B1q*gamma2qqb_1[ii1]*LF*LQ - 2.*gamma1qq_2[ii2]*gamma2qqb_1[ii1]*LF*LQ - (gamma1qg_1[ii1]*gamma2gq_1[ii1] + gamma1qq_1[ii1]*gamma2qqb_1[ii1])*LF*LQ - (gamma1gq_1[ii1]*gamma2qg_1[ii1] + gamma1qq_1[ii1]*gamma2qqb_1[ii1])*LF*LQ - (B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2*LQ)/2. - (beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2*LQ)/ 2. - (3.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2*LQ)/2. - (gamma1gq_1[ii1]*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + 2.*gamma1qq_1[ii1])*LF2*LQ)/2. + B1q*eC1qg_1*gamma1gq_1[ii1]*LQ2 - (3.*beta0*eC1qg_1*gamma1gq_1[ii1]*LQ2)/2. + eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LQ2 + (eC1qq_2*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ2)/2. + (gamma1gq_1[ii1]*(eC1qq_1*gamma1qg_1[ii1] + eC1qg_1*(gamma1gg_1[ii1] + gamma1qq_1[ii1]))*LQ2)/2. + B1q*gamma2qqb_1[ii1]*LQ2 - beta0*gamma2qqb_1[ii1]*LQ2 + gamma1qq_2[ii2]*gamma2qqb_1[ii1]*LQ2 + (gamma1qg_1[ii1]*gamma2gq_1[ii1] + gamma1qq_1[ii1]*gamma2qqb_1[ii1])*LQ2 - (A1q*eC1qg_1*gamma1gq_1[ii1]*LF*LQ2)/2. + B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ2 - (beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ2)/2. + (3.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ2)/2. + (gamma1gq_1[ii1]*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + 2.*gamma1qq_1[ii1])*LF*LQ2)/2. - (A1q*gamma2qqb_1[ii1]*LF*LQ2)/2. - (A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2*LQ2)/4. + (A1q*eC1qg_1*gamma1gq_1[ii1]*LQ3)/2. - (B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ3)/2. + (beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ3)/2. - (gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ3)/2. - (gamma1gq_1[ii1]*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + 2.*gamma1qq_1[ii1])*LQ3)/6. + (A1q*gamma2qqb_1[ii1]*LQ3)/2. + (A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ3)/2. - (A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ4)/4.;
	    complex <double> H3st_qq_2   = eC2qqb_2*gamma1qq_1[ii1]*LF + eC1qq_1*eC1qg_2*gamma1gq_2[ii2]*LF + (eC2qg_2*gamma1gq_2[ii2] + eC2qqb_2*gamma1qq_2[ii2])*LF + eC1qq_1*gamma2qqb_2[ii2]*LF + (eC1qg_2*gamma2gq_2[ii2] + eC1qq_2*gamma2qqb_2[ii2])*LF + gamma3qqb_2[ii2]*LF + (beta0*eC1qg_2*gamma1gq_2[ii2]*LF2)/2. + eC1qg_2*gamma1qq_1[ii1]*gamma1gq_2[ii2]*LF2 + (eC1qq_1*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2)/2. + (gamma1gq_2[ii2]*(eC1qq_2*gamma1qg_2[ii2] + eC1qg_2*(gamma1gg_2[ii2] + gamma1qq_2[ii2]))*LF2)/2. + beta0*gamma2qqb_2[ii2]*LF2 + gamma1qq_1[ii1]*gamma2qqb_2[ii2]*LF2 + (gamma1qg_2[ii2]*gamma2gq_2[ii2] + gamma1qq_2[ii2]*gamma2qqb_2[ii2])*LF2 + (beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF3)/2. + (gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF3)/2. + (gamma1gq_2[ii2]*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + 2.*gamma1qq_2[ii2])*LF3)/6. - eC2qqb_2*gamma1qq_1[ii1]*LQ - eC1qq_1*eC1qg_2*gamma1gq_2[ii2]*LQ - (eC2qg_2*gamma1gq_2[ii2] + eC2qqb_2*gamma1qq_2[ii2])*LQ - eC1qq_1*gamma2qqb_2[ii2]*LQ - (eC1qg_2*gamma2gq_2[ii2] + eC1qq_2*gamma2qqb_2[ii2])*LQ - gamma3qqb_2[ii2]*LQ - B1q*eC1qg_2*gamma1gq_2[ii2]*LF*LQ + beta0*eC1qg_2*gamma1gq_2[ii2]*LF*LQ - 2.*eC1qg_2*gamma1qq_1[ii1]*gamma1gq_2[ii2]*LF*LQ - eC1qq_1*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ - gamma1gq_2[ii2]*(eC1qq_2*gamma1qg_2[ii2] + eC1qg_2*(gamma1gg_2[ii2] + gamma1qq_2[ii2]))*LF*LQ - B1q*gamma2qqb_2[ii2]*LF*LQ - 2.*gamma1qq_1[ii1]*gamma2qqb_2[ii2]*LF*LQ - (gamma1qg_2[ii2]*gamma2gq_2[ii2] + gamma1qq_2[ii2]*gamma2qqb_2[ii2])*LF*LQ - (gamma1gq_2[ii2]*gamma2qg_2[ii2] + gamma1qq_2[ii2]*gamma2qqb_2[ii2])*LF*LQ - (B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2*LQ)/2. - (beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2*LQ)/ 2. - (3.*gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2*LQ)/2. - (gamma1gq_2[ii2]*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + 2.*gamma1qq_2[ii2])*LF2*LQ)/2. + B1q*eC1qg_2*gamma1gq_2[ii2]*LQ2 - (3.*beta0*eC1qg_2*gamma1gq_2[ii2]*LQ2)/2. + eC1qg_2*gamma1qq_1[ii1]*gamma1gq_2[ii2]*LQ2 + (eC1qq_1*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ2)/2. + (gamma1gq_2[ii2]*(eC1qq_2*gamma1qg_2[ii2] + eC1qg_2*(gamma1gg_2[ii2] + gamma1qq_2[ii2]))*LQ2)/2. + B1q*gamma2qqb_2[ii2]*LQ2 - beta0*gamma2qqb_2[ii2]*LQ2 + gamma1qq_1[ii1]*gamma2qqb_2[ii2]*LQ2 + (gamma1qg_2[ii2]*gamma2gq_2[ii2] + gamma1qq_2[ii2]*gamma2qqb_2[ii2])*LQ2 - (A1q*eC1qg_2*gamma1gq_2[ii2]*LF*LQ2)/2. + B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ2 - (beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ2)/2. + (3.*gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ2)/2. + (gamma1gq_2[ii2]*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + 2.*gamma1qq_2[ii2])*LF*LQ2)/2. - (A1q*gamma2qqb_2[ii2]*LF*LQ2)/2. - (A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2*LQ2)/4. + (A1q*eC1qg_2*gamma1gq_2[ii2]*LQ3)/2. - (B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ3)/2. + (beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ3)/2. - (gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ3)/2. - (gamma1gq_2[ii2]*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + 2.*gamma1qq_2[ii2])*LQ3)/6. + (A1q*gamma2qqb_2[ii2]*LQ3)/2. + (A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ3)/2. - (A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ4)/4.;
	    complex <double> H3st_qqp_1  = (2.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + 2.*gamma1qq_1[ii1])*LQF3 + 6.*gamma1gq_1[ii1]*(eC1qq_1*gamma1qg_1[ii1] + eC1qg_1*(gamma1gg_1[ii1] + gamma1qq_1[ii1]))* (LF2 + LQ2) - 3.*(-4.*eC2qg_1*gamma1gq_1[ii1]*LF - 4.*eC1qg_1*gamma2gq_1[ii1]*LF - 4.*eC1qq_1*gamma2qqbp_1[ii1]*LF - 4.*gamma3qqbp_1[ii1]*LF - 2.*beta0*eC1qg_1*gamma1gq_1[ii1]*LF2 - 4.*eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LF2 - 4.*gamma1qg_1[ii1]*gamma2gq_1[ii1]*LF2 - 4.*beta0*gamma2qqbp_1[ii1]*LF2 - 4.*gamma1qq_2[ii2]*gamma2qqbp_1[ii1]*LF2 - 4.*gamma1qq_1[ii1]*gamma2qqbp_1[ii1]*LF2 - 2.*beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF3 - 2.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]* LF3 - 4.*eC2qqbp_1*(gamma1qq_2[ii2] + gamma1qq_1[ii1])*LQF - 2.*eC1qq_2*(2.*eC1qg_1*gamma1gq_1[ii1] + 2.*gamma2qqbp_1[ii1] + gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQF)*LQF + 4.*eC2qg_1*gamma1gq_1[ii1]*LQ + 4.*eC1qg_1*gamma2gq_1[ii1]*LQ + 4.*eC1qq_1*gamma2qqbp_1[ii1]*LQ + 4.*gamma3qqbp_1[ii1]*LQ + 4.*B1q*eC1qg_1*gamma1gq_1[ii1]*LF*LQ - 4.*beta0*eC1qg_1*gamma1gq_1[ii1]*LF*LQ + 8.*eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LF*LQ + 4.*eC1qg_1*gamma1gg_1[ii1]*gamma1gq_1[ii1]*LF* LQ + 4.*eC1qq_1*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ + 4.*eC1qg_1*gamma1gq_1[ii1]* gamma1qq_1[ii1]*LF*LQ + 4.*gamma1qg_1[ii1]*gamma2gq_1[ii1]*LF*LQ + 4.*gamma1gq_1[ii1]*gamma2qg_1[ii1]*LF*LQ + 4.*B1q*gamma2qqbp_1[ii1]*LF*LQ + 8.*gamma1qq_2[ii2]*gamma2qqbp_1[ii1]*LF*LQ + 8.*gamma1qq_1[ii1]*gamma2qqbp_1[ii1]*LF*LQ + 2.*B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2*LQ + 2.*beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2* LQ + 6.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2*LQ - 4.*B1q*eC1qg_1*gamma1gq_1[ii1]*LQ2 + 6.*beta0*eC1qg_1*gamma1gq_1[ii1]*LQ2 - 4.*eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LQ2 - 4.*gamma1qg_1[ii1]*gamma2gq_1[ii1]*LQ2 - 4.*B1q*gamma2qqbp_1[ii1]*LQ2 + 4.*beta0*gamma2qqbp_1[ii1]*LQ2 - 4.*gamma1qq_2[ii2]*gamma2qqbp_1[ii1]*LQ2 - 4.*gamma1qq_1[ii1]*gamma2qqbp_1[ii1]*LQ2 + 2.*A1q*eC1qg_1*gamma1gq_1[ii1]*LF*LQ2 - 4.*B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ2 + 2.*beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ2 - 6.*gamma1qq_2[ii2]*gamma1gq_1[ii1]* gamma1qg_1[ii1]*LF*LQ2 + 2.*A1q*gamma2qqbp_1[ii1]*LF*LQ2 + A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2*LQ2 - 2.*A1q*eC1qg_1*gamma1gq_1[ii1]*LQ3 + 2.*B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ3 - 2.*beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ3 + 2.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ3 - 2.*A1q*gamma2qqbp_1[ii1]*LQ3 - 2.*A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ3 + A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ4))/12.;
	    complex <double> H3st_qqp_2  = (2.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + 2.*gamma1qq_2[ii2])*LQF3 + 6.*gamma1gq_2[ii2]*(eC1qq_2*gamma1qg_2[ii2] + eC1qg_2*(gamma1gg_2[ii2] + gamma1qq_2[ii2]))* (LF2 + LQ2) - 3.*(-4.*eC2qg_2*gamma1gq_2[ii2]*LF - 4.*eC1qg_2*gamma2gq_2[ii2]*LF - 4.*eC1qq_2*gamma2qqbp_2[ii2]*LF - 4.*gamma3qqbp_2[ii2]*LF - 2.*beta0*eC1qg_2*gamma1gq_2[ii2]*LF2 - 4.*eC1qg_2*gamma1qq_1[ii1]*gamma1gq_2[ii2]*LF2 - 4.*gamma1qg_2[ii2]*gamma2gq_2[ii2]*LF2 - 4.*beta0*gamma2qqbp_2[ii2]*LF2 - 4.*gamma1qq_1[ii1]*gamma2qqbp_2[ii2]*LF2 - 4.*gamma1qq_2[ii2]*gamma2qqbp_2[ii2]*LF2 - 2.*beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF3 - 2.*gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]* LF3 - 4.*eC2qqbp_2*(gamma1qq_1[ii1] + gamma1qq_2[ii2])*LQF - 2.*eC1qq_1*(2.*eC1qg_2*gamma1gq_2[ii2] + 2.*gamma2qqbp_2[ii2] + gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQF)*LQF + 4.*eC2qg_2*gamma1gq_2[ii2]*LQ + 4.*eC1qg_2*gamma2gq_2[ii2]*LQ + 4.*eC1qq_2*gamma2qqbp_2[ii2]*LQ + 4.*gamma3qqbp_2[ii2]*LQ + 4.*B1q*eC1qg_2*gamma1gq_2[ii2]*LF*LQ - 4.*beta0*eC1qg_2*gamma1gq_2[ii2]*LF*LQ + 8.*eC1qg_2*gamma1qq_1[ii1]*gamma1gq_2[ii2]*LF*LQ + 4.*eC1qg_2*gamma1gg_2[ii2]*gamma1gq_2[ii2]*LF* LQ + 4.*eC1qq_2*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ + 4.*eC1qg_2*gamma1gq_2[ii2]* gamma1qq_2[ii2]*LF*LQ + 4.*gamma1qg_2[ii2]*gamma2gq_2[ii2]*LF*LQ + 4.*gamma1gq_2[ii2]*gamma2qg_2[ii2]*LF*LQ + 4.*B1q*gamma2qqbp_2[ii2]*LF*LQ + 8.*gamma1qq_1[ii1]*gamma2qqbp_2[ii2]*LF*LQ + 8.*gamma1qq_2[ii2]*gamma2qqbp_2[ii2]*LF*LQ + 2.*B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2*LQ + 2.*beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2* LQ + 6.*gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2*LQ - 4.*B1q*eC1qg_2*gamma1gq_2[ii2]*LQ2 + 6.*beta0*eC1qg_2*gamma1gq_2[ii2]*LQ2 - 4.*eC1qg_2*gamma1qq_1[ii1]*gamma1gq_2[ii2]*LQ2 - 4.*gamma1qg_2[ii2]*gamma2gq_2[ii2]*LQ2 - 4.*B1q*gamma2qqbp_2[ii2]*LQ2 + 4.*beta0*gamma2qqbp_2[ii2]*LQ2 - 4.*gamma1qq_1[ii1]*gamma2qqbp_2[ii2]*LQ2 - 4.*gamma1qq_2[ii2]*gamma2qqbp_2[ii2]*LQ2 + 2.*A1q*eC1qg_2*gamma1gq_2[ii2]*LF*LQ2 - 4.*B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ2 + 2.*beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ2 - 6.*gamma1qq_1[ii1]*gamma1gq_2[ii2]* gamma1qg_2[ii2]*LF*LQ2 + 2.*A1q*gamma2qqbp_2[ii2]*LF*LQ2 + A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2*LQ2 - 2.*A1q*eC1qg_2*gamma1gq_2[ii2]*LQ3 + 2.*B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ3 - 2.*beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ3 + 2.*gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ3 - 2.*A1q*gamma2qqbp_2[ii2]*LQ3 - 2.*A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ3 + A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ4))/12.;
	    complex <double> H3st_qqbp_1 = (2.*gamma1gq_1[ii1]*gamma1qg_1[ii1]*(gamma1gg_1[ii1] + 2.*gamma1qq_1[ii1])*LQF3 + 6.*gamma1gq_1[ii1]*(eC1qq_1*gamma1qg_1[ii1] + eC1qg_1*(gamma1gg_1[ii1] + gamma1qq_1[ii1]))* (LF2 + LQ2) - 3.*(-4.*eC2qg_1*gamma1gq_1[ii1]*LF - 4.*eC1qg_1*gamma2gq_1[ii1]*LF - 4.*eC1qq_1*gamma2qqp_1[ii1]*LF - 4.*gamma3qqp_1[ii1]*LF - 2.*beta0*eC1qg_1*gamma1gq_1[ii1]* LF2 - 4.*eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LF2 - 4.*gamma1qg_1[ii1]*gamma2gq_1[ii1]* LF2 - 4.*beta0*gamma2qqp_1[ii1]*LF2 - 4.*gamma1qq_2[ii2]*gamma2qqp_1[ii1]*LF2 - 4.*gamma1qq_1[ii1]*gamma2qqp_1[ii1]*LF2 - 2.*beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF3 - 2.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF3 - 4.*eC2qqp_1*(gamma1qq_2[ii2] + gamma1qq_1[ii1])*LQF - 2.*eC1qq_2*(2.*eC1qg_1*gamma1gq_1[ii1] + 2.*gamma2qqp_1[ii1] + gamma1gq_1[ii1]*gamma1qg_1[ii1]* LQF)*LQF + 4.*eC2qg_1*gamma1gq_1[ii1]*LQ + 4.*eC1qg_1*gamma2gq_1[ii1]*LQ + 4.*eC1qq_1*gamma2qqp_1[ii1]*LQ + 4.*gamma3qqp_1[ii1]*LQ + 4.*B1q*eC1qg_1*gamma1gq_1[ii1]*LF*LQ - 4.*beta0*eC1qg_1*gamma1gq_1[ii1]*LF*LQ + 8.*eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LF*LQ + 4.*eC1qg_1*gamma1gg_1[ii1]*gamma1gq_1[ii1]*LF* LQ + 4.*eC1qq_1*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ + 4.*eC1qg_1*gamma1gq_1[ii1]* gamma1qq_1[ii1]*LF*LQ + 4.*gamma1qg_1[ii1]*gamma2gq_1[ii1]*LF*LQ + 4.*gamma1gq_1[ii1]*gamma2qg_1[ii1]*LF*LQ + 4.*B1q*gamma2qqp_1[ii1]*LF*LQ + 8.*gamma1qq_2[ii2]*gamma2qqp_1[ii1]*LF*LQ + 8.*gamma1qq_1[ii1]*gamma2qqp_1[ii1]*LF*LQ + 2.*B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2*LQ + 2.*beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2* LQ + 6.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2*LQ - 4.*B1q*eC1qg_1*gamma1gq_1[ii1]*LQ2 + 6.*beta0*eC1qg_1*gamma1gq_1[ii1]*LQ2 - 4.*eC1qg_1*gamma1qq_2[ii2]*gamma1gq_1[ii1]*LQ2 - 4.*gamma1qg_1[ii1]*gamma2gq_1[ii1]*LQ2 - 4.*B1q*gamma2qqp_1[ii1]*LQ2 + 4.*beta0*gamma2qqp_1[ii1]*LQ2 - 4.*gamma1qq_2[ii2]*gamma2qqp_1[ii1]*LQ2 - 4.*gamma1qq_1[ii1]*gamma2qqp_1[ii1]*LQ2 + 2.*A1q*eC1qg_1*gamma1gq_1[ii1]*LF*LQ2 - 4.*B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ2 + 2.*beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ2 - 6.*gamma1qq_2[ii2]*gamma1gq_1[ii1]* gamma1qg_1[ii1]*LF*LQ2 + 2.*A1q*gamma2qqp_1[ii1]*LF*LQ2 + A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF2*LQ2 - 2.*A1q*eC1qg_1*gamma1gq_1[ii1]*LQ3 + 2.*B1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ3 - 2.*beta0*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ3 + 2.*gamma1qq_2[ii2]*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ3 - 2.*A1q*gamma2qqp_1[ii1]*LQ3 - 2.*A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LF*LQ3 + A1q*gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQ4))/12.;
	    complex <double> H3st_qqbp_2 = (2.*gamma1gq_2[ii2]*gamma1qg_2[ii2]*(gamma1gg_2[ii2] + 2.*gamma1qq_2[ii2])*LQF3 + 6.*gamma1gq_2[ii2]*(eC1qq_2*gamma1qg_2[ii2] + eC1qg_2*(gamma1gg_2[ii2] + gamma1qq_2[ii2]))* (LF2 + LQ2) - 3.*(-4.*eC2qg_2*gamma1gq_2[ii2]*LF - 4.*eC1qg_2*gamma2gq_2[ii2]*LF - 4.*eC1qq_2*gamma2qqp_2[ii2]*LF - 4.*gamma3qqp_2[ii2]*LF - 2.*beta0*eC1qg_2*gamma1gq_2[ii2]* LF2 - 4.*eC1qg_2*gamma1qq_1[ii1]*gamma1gq_2[ii2]*LF2 - 4.*gamma1qg_2[ii2]*gamma2gq_2[ii2]* LF2 - 4.*beta0*gamma2qqp_2[ii2]*LF2 - 4.*gamma1qq_1[ii1]*gamma2qqp_2[ii2]*LF2 - 4.*gamma1qq_2[ii2]*gamma2qqp_2[ii2]*LF2 - 2.*beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF3 - 2.*gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF3 - 4.*eC2qqp_2*(gamma1qq_1[ii1] + gamma1qq_2[ii2])*LQF - 2.*eC1qq_1*(2.*eC1qg_2*gamma1gq_2[ii2] + 2.*gamma2qqp_2[ii2] + gamma1gq_2[ii2]*gamma1qg_2[ii2]* LQF)*LQF + 4.*eC2qg_2*gamma1gq_2[ii2]*LQ + 4.*eC1qg_2*gamma2gq_2[ii2]*LQ + 4.*eC1qq_2*gamma2qqp_2[ii2]*LQ + 4.*gamma3qqp_2[ii2]*LQ + 4.*B1q*eC1qg_2*gamma1gq_2[ii2]*LF*LQ - 4.*beta0*eC1qg_2*gamma1gq_2[ii2]*LF*LQ + 8.*eC1qg_2*gamma1qq_1[ii1]*gamma1gq_2[ii2]*LF*LQ + 4.*eC1qg_2*gamma1gg_2[ii2]*gamma1gq_2[ii2]*LF* LQ + 4.*eC1qq_2*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ + 4.*eC1qg_2*gamma1gq_2[ii2]* gamma1qq_2[ii2]*LF*LQ + 4.*gamma1qg_2[ii2]*gamma2gq_2[ii2]*LF*LQ + 4.*gamma1gq_2[ii2]*gamma2qg_2[ii2]*LF*LQ + 4.*B1q*gamma2qqp_2[ii2]*LF*LQ + 8.*gamma1qq_1[ii1]*gamma2qqp_2[ii2]*LF*LQ + 8.*gamma1qq_2[ii2]*gamma2qqp_2[ii2]*LF*LQ + 2.*B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2*LQ + 2.*beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2* LQ + 6.*gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2*LQ - 4.*B1q*eC1qg_2*gamma1gq_2[ii2]*LQ2 + 6.*beta0*eC1qg_2*gamma1gq_2[ii2]*LQ2 - 4.*eC1qg_2*gamma1qq_1[ii1]*gamma1gq_2[ii2]*LQ2 - 4.*gamma1qg_2[ii2]*gamma2gq_2[ii2]*LQ2 - 4.*B1q*gamma2qqp_2[ii2]*LQ2 + 4.*beta0*gamma2qqp_2[ii2]*LQ2 - 4.*gamma1qq_1[ii1]*gamma2qqp_2[ii2]*LQ2 - 4.*gamma1qq_2[ii2]*gamma2qqp_2[ii2]*LQ2 + 2.*A1q*eC1qg_2*gamma1gq_2[ii2]*LF*LQ2 - 4.*B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ2 + 2.*beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ2 - 6.*gamma1qq_1[ii1]*gamma1gq_2[ii2]* gamma1qg_2[ii2]*LF*LQ2 + 2.*A1q*gamma2qqp_2[ii2]*LF*LQ2 + A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF2*LQ2 - 2.*A1q*eC1qg_2*gamma1gq_2[ii2]*LQ3 + 2.*B1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ3 - 2.*beta0*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ3 + 2.*gamma1qq_1[ii1]*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ3 - 2.*A1q*gamma2qqp_2[ii2]*LQ3 - 2.*A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LF*LQ3 + A1q*gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQ4))/12.;
	    complex <double> H3st_gg     = (LQF*(2.*eC2qg_2*gamma1qg_1[ii1] + 2.*eC2qg_1*gamma1qg_2[ii2] + 2.*eC1qg_2*gamma2qg_1[ii1] + beta0*eC1qg_2*gamma1qg_1[ii1]*LF + eC1qg_2*gamma1gg_1[ii1]*gamma1qg_1[ii1]*LF + eC1qg_2*gamma1qg_1[ii1]*gamma1qq_1[ii1]*LF + 2.*eC1qg_2*gamma1qg_1[ii1]*gamma1gg_2[ii2]*LF + 2.*eC1qq_2*gamma1qg_1[ii1]*gamma1qg_2[ii2]*LF + 2.*gamma1qg_2[ii2]*gamma2qg_1[ii1]*LF + 2.*gamma1qg_1[ii1]*gamma2qg_2[ii2]*LF + 2.*beta0*gamma1qg_1[ii1]*gamma1qg_2[ii2]*LF2 + gamma1gg_1[ii1]*gamma1qg_1[ii1]*gamma1qg_2[ii2]*LF2 + gamma1qg_1[ii1]*gamma1qq_1[ii1]*gamma1qg_2[ii2]* LF2 + gamma1qg_1[ii1]*gamma1gg_2[ii2]*gamma1qg_2[ii2]*LF2 + gamma1qg_1[ii1]*gamma1qg_2[ii2]*gamma1qq_2[ii2]*LF2 - 2.*B1q*eC1qg_2*gamma1qg_1[ii1]*LQ + 3.*beta0*eC1qg_2*gamma1qg_1[ii1]*LQ - eC1qg_2*gamma1gg_1[ii1]*gamma1qg_1[ii1]*LQ - eC1qg_2*gamma1qg_1[ii1]*gamma1qq_1[ii1]*LQ - 2.*eC1qg_2*gamma1qg_1[ii1]*gamma1gg_2[ii2]*LQ - 2.*eC1qq_2*gamma1qg_1[ii1]*gamma1qg_2[ii2]*LQ - 2.*gamma1qg_2[ii2]*gamma2qg_1[ii1]*LQ - 2.*gamma1qg_1[ii1]*gamma2qg_2[ii2]*LQ - 2.*B1q*gamma1qg_1[ii1]*gamma1qg_2[ii2]*LF*LQ - 2.*gamma1gg_1[ii1]*gamma1qg_1[ii1]*gamma1qg_2[ii2]*LF*LQ - 2.*gamma1qg_1[ii1]*gamma1qq_1[ii1]* gamma1qg_2[ii2]*LF*LQ - 2.*gamma1qg_1[ii1]*gamma1gg_2[ii2]*gamma1qg_2[ii2]*LF*LQ - 2.*gamma1qg_1[ii1]*gamma1qg_2[ii2]*gamma1qq_2[ii2]*LF*LQ - A1q*eC1qg_2*gamma1qg_1[ii1]*LQ2 + 2.*B1q*gamma1qg_1[ii1]*gamma1qg_2[ii2]*LQ2 - 2.*beta0*gamma1qg_1[ii1]*gamma1qg_2[ii2]*LQ2 + gamma1gg_1[ii1]*gamma1qg_1[ii1]*gamma1qg_2[ii2]*LQ2 + gamma1qg_1[ii1]*gamma1qq_1[ii1]*gamma1qg_2[ii2]* LQ2 + gamma1qg_1[ii1]*gamma1gg_2[ii2]*gamma1qg_2[ii2]*LQ2 + gamma1qg_1[ii1]*gamma1qg_2[ii2]*gamma1qq_2[ii2]*LQ2 - A1q*gamma1qg_1[ii1]*gamma1qg_2[ii2]*LF*LQ2 + A1q*gamma1qg_1[ii1]*gamma1qg_2[ii2]*LQ3 + 2.*eC1qq_1*gamma1qg_1[ii1]* (eC1qg_2 + gamma1qg_2[ii2]*LF - gamma1qg_2[ii2]*LQ) + eC1qg_1*(2.*eC1qg_2*(gamma1gg_1[ii1] + gamma1gg_2[ii2]) + 2.*eC1qq_2*gamma1qg_2[ii2] + 2.*gamma2qg_2[ii2] + beta0*gamma1qg_2[ii2]*LF + 2.*gamma1gg_1[ii1]*gamma1qg_2[ii2]*LF + gamma1gg_2[ii2]*gamma1qg_2[ii2]*LF + gamma1qg_2[ii2]*gamma1qq_2[ii2]*LF - 2.*B1q*gamma1qg_2[ii2]*LQ + 3.*beta0*gamma1qg_2[ii2]*LQ - 2.*gamma1gg_1[ii1]*gamma1qg_2[ii2]* LQ - gamma1gg_2[ii2]*gamma1qg_2[ii2]*LQ - gamma1qg_2[ii2]*gamma1qq_2[ii2]*LQ - A1q*gamma1qg_2[ii2]*LQ2)))/2.;
	    complex <double> H3st_qbg_1  = ((2.*eC2qqb_2*gamma1qg_1[ii1] + (2.*eC1qg_2*gamma1gq_2[ii2] + 2.*gamma2qqb_2[ii2] + gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQF)*(eC1qg_1 + gamma1qg_1[ii1]*LQF))* LQF)/2.;
	    complex <double> H3st_qbg_2  = ((2.*eC2qqb_1*gamma1qg_2[ii2] + (2.*eC1qg_1*gamma1gq_1[ii1] + 2.*gamma2qqb_1[ii1] + gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQF)*(eC1qg_2 + gamma1qg_2[ii2]*LQF))* LQF)/2.;
	    complex <double> H3st_qpg_1  = ((2.*eC2qqp_2*gamma1qg_1[ii1] + (2.*eC1qg_2*gamma1gq_2[ii2] + 2.*gamma2qqp_2[ii2] + gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQF)*(eC1qg_1 + gamma1qg_1[ii1]*LQF))* LQF)/2.;
	    complex <double> H3st_qpg_2  = ((2.*eC2qqp_1*gamma1qg_2[ii2] + (2.*eC1qg_1*gamma1gq_1[ii1] + 2.*gamma2qqp_1[ii1] + gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQF)*(eC1qg_2 + gamma1qg_2[ii2]*LQF))* LQF)/2.;
	    complex <double> H3st_qbpg_1 = ((2.*eC2qqbp_2*gamma1qg_1[ii1] + (2.*eC1qg_2*gamma1gq_2[ii2] + 2.*gamma2qqbp_2[ii2] + gamma1gq_2[ii2]*gamma1qg_2[ii2]*LQF)*(eC1qg_1 + gamma1qg_1[ii1]*LQF))* LQF)/2.;
	    complex <double> H3st_qbpg_2 = ((2.*eC2qqbp_1*gamma1qg_2[ii2] + (2.*eC1qg_1*gamma1gq_1[ii1] + 2.*gamma2qqbp_1[ii1] + gamma1gq_1[ii1]*gamma1qg_1[ii1]*LQF)*(eC1qg_2 + gamma1qg_2[ii2]*LQF))* LQF)/2.;

	    //Mixed LR*LQF variations
	    H3st_qqb    += -2*beta0*H2st_qqb*LR + (-beta1*LR - beta0*LR2)*H1st_qqb;
	    H3st_qq_1   += -2*beta0*H2st_qq_1*LR;
	    H3st_qq_2   += -2*beta0*H2st_qq_2*LR;
	    H3st_qqp_1  += -2*beta0*H2st_qqp_1*LR;
	    H3st_qqp_2  += -2*beta0*H2st_qqp_1*LR;
	    H3st_qqbp_1 += -2*beta0*H2st_qqbp_1*LR;
	    H3st_qqbp_2 += -2*beta0*H2st_qqbp_2*LR;
	    H3st_qg_1   += -2*beta0*H2st_qg_1*LR + (-beta1*LR - beta0*LR2)*H1st_qg_1;
	    H3st_qg_2   += -2*beta0*H2st_qg_2*LR + (-beta1*LR - beta0*LR2)*H1st_qg_2;
	    H3st_gg     += -2*beta0*H2st_gg*LR;
	    
	    qqb[idx]    += as3*H3st_qqb;
	    qg_1[idx]   += as3*H3st_qg_1;
	    qg_2[idx]   += as3*H3st_qg_2;
	    qq_1[idx]   += as3*H3st_qq_1;
	    qq_2[idx]   += as3*H3st_qq_2;
	    qqp_1[idx]  += as3*H3st_qqp_1;
	    qqp_2[idx]  += as3*H3st_qqp_2;
	    qqbp_1[idx] += as3*H3st_qqbp_1;
	    qqbp_2[idx] += as3*H3st_qqbp_2;
	    gg[idx]     += as3*H3st_gg;
	    qbg_1[idx]  += as3*H3st_qbg_1;
	    qbg_2[idx]  += as3*H3st_qbg_2;
	    qpg_1[idx]  += as3*H3st_qpg_1;
	    qpg_2[idx]  += as3*H3st_qpg_2;
	    qbpg_1[idx] += as3*H3st_qbpg_1;
	    qbpg_2[idx] += as3*H3st_qbpg_2;
	  }
}
