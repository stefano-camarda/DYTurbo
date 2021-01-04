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

using namespace ccoeff;
using namespace pmom;
using namespace resconst;
using namespace resint;

  //allocate memory
void muf::allocate()
{
  if (opts.order == 0)
    return;

  if (opts.mellin1d)
    {
      qqb  = new complex <double> [mellinint::mdim];
      qg   = new complex <double> [mellinint::mdim];
      qq   = new complex <double> [mellinint::mdim];
      qqp  = new complex <double> [mellinint::mdim];
      qqbp = new complex <double> [mellinint::mdim];
      gg   = new complex <double> [mellinint::mdim];
      qbg  = new complex <double> [mellinint::mdim];
      qpg  = new complex <double> [mellinint::mdim];
      qbpg = new complex <double> [mellinint::mdim];
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
  if (opts.order == 0)
    return;

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
  if (opts.order == 0)
    return;

  if (opts.mellin1d)
    {
      fill(qqb , qqb  + mellinint::mdim, 0.);
      fill(qg  , qg   + mellinint::mdim, 0.);
      fill(qq  , qq   + mellinint::mdim, 0.);
      fill(qqp , qqp  + mellinint::mdim, 0.);
      fill(qqbp, qqbp + mellinint::mdim, 0.);
      fill(gg  , gg   + mellinint::mdim, 0.);
      fill(qbg , qbg  + mellinint::mdim, 0.);
      fill(qpg , qpg  + mellinint::mdim, 0.);
      fill(qbpg, qbpg + mellinint::mdim, 0.);
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
  if (opts.order == 0)
    return;

  if (opts.evolmode == 1 && opts.mufevol)
    return;
  
  double LQ2 = pow(LQ,2);
  double LQ3 = pow(LQ,3);
  double LQ4 = pow(LQ,4);
  double LQ5 = pow(LQ,5);
  double LQ6 = pow(LQ,6);

  double LR2 = pow(LR,2);

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

  double as = aass;
  double as2 = pow(aass,2);
  double as3 = pow(aass,3);
  
  //cout << "mub " << real(pdfevol::mub) << "  " << LQb << "  " << LFb << "  " << LQFb << endl;

  if (opts.mellin1d)
    //rapidity integrated
    for (int i = 0; i < mellinint::mdim; i++)
      {
	int idx = pmom::index(i,mesq::positive);
      
	//NLO
	complex <double> H1st_qqb = 2.*gamma1qq[idx]*LQFb;
	complex <double> H1st_qg  = gamma1qg[idx]*LQFb;
	
	qqb[i] += as*H1st_qqb;
	qg[i]  += as*H1st_qg;
      
	if (opts.order == 1)
	  continue;

	//NNLO
	complex <double> H2st_qqb  = LQFb*(2.*C1qg[idx]*gamma1gq[idx] + 2.*gamma2qq[idx] + gamma1gq[idx]*gamma1qg[idx]*LQFb + gamma1qq[idx]*(4.*C1qq[idx] + 2.*gamma1qq[idx]*LQFb - 2.*B1q*LQ - A1q*LQ2 + beta0*(LFb + LQb)));
	complex <double> H2st_qq   = ((2.*C1qg[idx]*gamma1gq[idx] + 2.*gamma2qqb[idx] + gamma1gq[idx]*gamma1qg[idx]*LQFb)*LQFb)/2.;
	complex <double> H2st_qqp  = ((2.*C1qg[idx]*gamma1gq[idx] + 2.*gamma2qqp[idx] + gamma1gq[idx]*gamma1qg[idx]*LQFb)*LQFb)/2.;
	complex <double> H2st_qqbp = ((2.*C1qg[idx]*gamma1gq[idx] + 2.*gamma2qqp[idx] + gamma1gq[idx]*gamma1qg[idx]*LQFb)*LQFb)/2.;
	complex <double> H2st_qg   = (LQFb*(4.*C1qq[idx]*gamma1qg[idx] + 2.*C1qg[idx]*(gamma1gg[idx] + gamma1qq[idx]) + 2.*gamma2qg[idx] + beta0*gamma1qg[idx]*LFb + gamma1gg[idx]*gamma1qg[idx]*LFb + 3.*gamma1qg[idx]*gamma1qq[idx]*LFb - 2.*B1q*gamma1qg[idx]*LQ + beta0*gamma1qg[idx]*LQb - gamma1gg[idx]*gamma1qg[idx]*LQb - 3.*gamma1qg[idx]*gamma1qq[idx]*LQb - A1q*gamma1qg[idx]*LQ2))/2.;
	complex <double> H2st_gg   = gamma1qg[idx]*(2.*C1qg[idx] + gamma1qg[idx]*LQFb)*LQFb;

	qqb[i]  += as2*H2st_qqb;
	qg[i]   += as2*H2st_qg;
	qq[i]   += as2*H2st_qq;
	qqp[i]  += as2*H2st_qqp;
	qqbp[i] += as2*H2st_qqbp;
	gg[i]   += as2*H2st_gg;
      
	if (opts.order == 2)
	  continue;

	//NNNLO
	double nf = 5.;
	complex <double> H3st_qqb  = 2.*pow(C1qq[idx],2)*gamma1qq[idx]*LFb + 4.*C2qq[idx]*gamma1qq[idx]*LFb + 4.*C1qq[idx]*gamma2qq[idx]*LFb + 2.*gamma3qq[idx]*LFb + 2.*C1qq[idx]*gamma1gq[idx]*gamma1qg[idx]*LFb2 + beta1*gamma1qq[idx]*LFb2 + 2.*beta0*C1qq[idx]*gamma1qq[idx]*LFb2 + 4.*C1qq[idx]*pow(gamma1qq[idx],2)*LFb2 + 2.*gamma1qg[idx]*gamma2gq[idx]*LFb2 + 2.*beta0*gamma2qq[idx]*LFb2 + 4.*gamma1qq[idx]*gamma2qq[idx]*LFb2 + beta0*gamma1gq[idx]*gamma1qg[idx]*LFb3 + (gamma1gg[idx]*gamma1gq[idx]*gamma1qg[idx]*LFb3)/3. + (2.*pow(beta0,2)*gamma1qq[idx]*LFb3)/3. + (5.*gamma1gq[idx]*gamma1qg[idx]*gamma1qq[idx]*LFb3)/3. + 2.*beta0*pow(gamma1qq[idx],2)*LFb3 + (4.*pow(gamma1qq[idx],3)*LFb3)/3. + 2.*C2qg[idx]*gamma1gq[idx]*LQFb - 2.*pow(C1qq[idx],2)*gamma1qq[idx]*LQb - 4.*C2qq[idx]*gamma1qq[idx]*LQb - 4.*C1qq[idx]*gamma2qq[idx]*LQb - 2.*gamma3qq[idx]*LQb - 4.*C1qq[idx]*gamma1gq[idx]*gamma1qg[idx]*LFb*LQb - 2.*B2q*gamma1qq[idx]*LFb*LQb - 4.*B1q*C1qq[idx]*gamma1qq[idx]*LFb*LQb + 4.*beta0*C1qq[idx]*gamma1qq[idx]*LFb*LQb - 8.*C1qq[idx]*pow(gamma1qq[idx],2)*LFb*LQb - 2.*gamma1qg[idx]*gamma2gq[idx]*LFb*LQb - 2.*gamma1gq[idx]*gamma2qg[idx]*LFb*LQb - 2.*B1q*gamma2qq[idx]*LFb*LQb - 8.*gamma1qq[idx]*gamma2qq[idx]*LFb*LQb - B1q*gamma1gq[idx]*gamma1qg[idx]*LFb2*LQb - beta0*gamma1gq[idx]*gamma1qg[idx]*LFb2*LQb - gamma1gg[idx]*gamma1gq[idx]*gamma1qg[idx]*LFb2*LQb - B1q*beta0*gamma1qq[idx]*LFb2*LQb - 5.*gamma1gq[idx]*gamma1qg[idx]*gamma1qq[idx]*LFb2*LQb - 2.*B1q*pow(gamma1qq[idx],2)*LFb2*LQb - 2.*beta0*pow(gamma1qq[idx],2)*LFb2*LQb - 4.*pow(gamma1qq[idx],3)*LFb2*LQb + 2.*C1qq[idx]*gamma1gq[idx]*gamma1qg[idx]*LQb2 + 2.*B2q*gamma1qq[idx]*LQb2 - beta1*gamma1qq[idx]*LQb2 + 4.*B1q*C1qq[idx]*gamma1qq[idx]*LQb2 - 6.*beta0*C1qq[idx]*gamma1qq[idx]*LQb2 + 4.*C1qq[idx]*pow(gamma1qq[idx],2)*LQb2 + 2.*gamma1qg[idx]*gamma2gq[idx]*LQb2 + 2.*B1q*gamma2qq[idx]*LQb2 - 2.*beta0*gamma2qq[idx]*LQb2 + 4.*gamma1qq[idx]*gamma2qq[idx]*LQb2 + 2.*B1q*gamma1gq[idx]*gamma1qg[idx]*LFb*LQb2 - beta0*gamma1gq[idx]*gamma1qg[idx]*LFb*LQb2 + gamma1gg[idx]*gamma1gq[idx]*gamma1qg[idx]*LFb*LQb2 - A2q*gamma1qq[idx]*LFb*LQb2 + pow(B1q,2)*gamma1qq[idx]*LFb*LQb2 - B1q*beta0*gamma1qq[idx]*LFb*LQb2 - 2.*A1q*C1qq[idx]*gamma1qq[idx]*LFb*LQb2 + 5.*gamma1gq[idx]*gamma1qg[idx]*gamma1qq[idx]*LFb*LQb2 + 4.*B1q*pow(gamma1qq[idx],2)*LFb*LQb2 - 2.*beta0*pow(gamma1qq[idx],2)*LFb*LQb2 + 4.*pow(gamma1qq[idx],3)*LFb*LQb2 - A1q*gamma2qq[idx]*LFb*LQb2 - (A1q*gamma1gq[idx]*gamma1qg[idx]*LFb2*LQb2)/2. - (A1q*beta0*gamma1qq[idx]*LFb2*LQb2)/2. - A1q*pow(gamma1qq[idx],2)*LFb2*LQb2 - B1q*gamma1gq[idx]*gamma1qg[idx]*LQb3 + beta0*gamma1gq[idx]*gamma1qg[idx]*LQb3 - (gamma1gg[idx]*gamma1gq[idx]*gamma1qg[idx]*LQb3)/3. + A2q*gamma1qq[idx]*LQb3 - pow(B1q,2)*gamma1qq[idx]*LQb3 + 2.*B1q*beta0*gamma1qq[idx]*LQb3 - (2.*pow(beta0,2)*gamma1qq[idx]*LQb3)/3. + 2.*A1q*C1qq[idx]*gamma1qq[idx]*LQb3 - (5.*gamma1gq[idx]*gamma1qg[idx]*gamma1qq[idx]*LQb3)/3. - 2.*B1q*pow(gamma1qq[idx],2)*LQb3 + 2.*beta0*pow(gamma1qq[idx],2)*LQb3 - (4.*pow(gamma1qq[idx],3)*LQb3)/3. + A1q*gamma2qq[idx]*LQb3 + A1q*gamma1gq[idx]*gamma1qg[idx]*LFb*LQb3 + A1q*B1q*gamma1qq[idx]*LFb*LQb3 - (2.*A1q*beta0*gamma1qq[idx]*LFb*LQb3)/3. + 2.*A1q*pow(gamma1qq[idx],2)*LFb*LQb3 - (A1q*gamma1gq[idx]*gamma1qg[idx]*LQb4)/2. - A1q*B1q*gamma1qq[idx]*LQb4 + (7.*A1q*beta0*gamma1qq[idx]*LQb4)/6. - A1q*pow(gamma1qq[idx],2)*LQb4 + (pow(A1q,2)*gamma1qq[idx]*LFb*LQb4)/4. - (pow(A1q,2)*gamma1qq[idx]*LQb5)/4. + C1qg[idx]*LQFb*(2.*C1qq[idx]*gamma1gq[idx] + 2.*gamma2gq[idx] + gamma1gq[idx]*(beta0*LFb + gamma1gg[idx]*LFb + 3.*gamma1qq[idx]*LFb - 2.*B1q*LQb + 3.*beta0*LQb - gamma1gg[idx]*LQb - 3.*gamma1qq[idx]*LQb - A1q*LQb2));
	complex <double> H3st_qq   = (2.*LFb*(6.*(C2qg[idx]*gamma1gq[idx] + 2.*C2qqb[idx]*gamma1qq[idx] + 2.*C1qq[idx]*gamma2qqb[idx] + gamma3qqb[idx]) + 6.*(C1qq[idx]*gamma1gq[idx]*gamma1qg[idx] + gamma1qg[idx]*gamma2gq[idx] + (beta0 + 2.*gamma1qq[idx])* gamma2qqb[idx])*LFb + gamma1gq[idx]*gamma1qg[idx]*(3.*beta0 + gamma1gg[idx] + 5.*gamma1qq[idx])* LFb2) - 6.*(2.*C2qg[idx]*gamma1gq[idx] + 4.*C2qqb[idx]*gamma1qq[idx] + 4.*C1qq[idx]*gamma2qqb[idx] + 2.*gamma3qqb[idx] + 2.*(2.*C1qq[idx]*gamma1gq[idx]*gamma1qg[idx] + gamma1qg[idx]*gamma2gq[idx] + gamma1gq[idx]*gamma2qg[idx] + B1q*gamma2qqb[idx] + 4.*gamma1qq[idx]*gamma2qqb[idx])*LFb + gamma1gq[idx]*gamma1qg[idx]*(B1q + beta0 + gamma1gg[idx] + 5.*gamma1qq[idx])*LFb2)*LQb + 3.*(4.*C1qq[idx]*gamma1gq[idx]*gamma1qg[idx] + 4.*gamma1qg[idx]*gamma2gq[idx] + 4.*(B1q - beta0 + 2.*gamma1qq[idx])*gamma2qqb[idx] - 2.*A1q*gamma2qqb[idx]*LFb + gamma1gq[idx]*gamma1qg[idx]*LFb*(4.*B1q - 2.*beta0 + 2.*gamma1gg[idx] + 10.*gamma1qq[idx] - A1q*LFb))* LQb2 - 2.*(-3*A1q*gamma2qqb[idx] + gamma1gq[idx]*gamma1qg[idx]*(3.*B1q - 3.*beta0 + gamma1gg[idx] + 5.*gamma1qq[idx] - 3.*A1q*LFb))*LQb3 - 3.*A1q*gamma1gq[idx]*gamma1qg[idx]*LQb4 + 6.*C1qg[idx]*LQFb*(2.*gamma2gq[idx] + gamma1gq[idx]*(2.*C1qq[idx] + (beta0 + gamma1gg[idx] + 3.*gamma1qq[idx])*LFb - (2.*B1q - 3.*beta0 + gamma1gg[idx] + 3.*gamma1qq[idx])*LQb - A1q*LQb2)))/12.;
	complex <double> H3st_qqp  = (2.*LFb*(6.*(C2qg[idx]*gamma1gq[idx] + 2.*C2qqbp[idx]*gamma1qq[idx] + 2.*C1qq[idx]*gamma2qqbp[idx] + gamma3qqbp[idx]) + 6.*(C1qq[idx]*gamma1gq[idx]*gamma1qg[idx] + gamma1qg[idx]*gamma2gq[idx] + (beta0 + 2.*gamma1qq[idx])*gamma2qqbp[idx])*LFb + gamma1gq[idx]*gamma1qg[idx]* (3.*beta0 + gamma1gg[idx] + 5.*gamma1qq[idx])*LFb2) - 6.*(2.*C2qg[idx]*gamma1gq[idx] + 4.*C2qqbp[idx]*gamma1qq[idx] + 4.*C1qq[idx]*gamma2qqbp[idx] + 2.*gamma3qqbp[idx] + 2.*(2.*C1qq[idx]*gamma1gq[idx]*gamma1qg[idx] + gamma1qg[idx]*gamma2gq[idx] + gamma1gq[idx]*gamma2qg[idx] + B1q*gamma2qqbp[idx] + 4.*gamma1qq[idx]*gamma2qqbp[idx])*LFb + gamma1gq[idx]*gamma1qg[idx]*(B1q + beta0 + gamma1gg[idx] + 5.*gamma1qq[idx])*LFb2)*LQb + 3.*(4.*C1qq[idx]*gamma1gq[idx]*gamma1qg[idx] + 4.*gamma1qg[idx]*gamma2gq[idx] + 4.*(B1q - beta0 + 2.*gamma1qq[idx])*gamma2qqbp[idx] - 2.*A1q*gamma2qqbp[idx]*LFb + gamma1gq[idx]*gamma1qg[idx]*LFb*(4.*B1q - 2.*beta0 + 2.*gamma1gg[idx] + 10.*gamma1qq[idx] - A1q*LFb))* LQb2 - 2.*(-3*A1q*gamma2qqbp[idx] + gamma1gq[idx]*gamma1qg[idx]* (3.*B1q - 3.*beta0 + gamma1gg[idx] + 5.*gamma1qq[idx] - 3.*A1q*LFb))*LQb3 - 3.*A1q*gamma1gq[idx]*gamma1qg[idx]*LQb4 + 6.*C1qg[idx]*LQFb* (2.*gamma2gq[idx] + gamma1gq[idx]*(2.*C1qq[idx] + (beta0 + gamma1gg[idx] + 3.*gamma1qq[idx])*LFb - (2.*B1q - 3.*beta0 + gamma1gg[idx] + 3.*gamma1qq[idx])*LQb - A1q*LQb2)))/12.;
	complex <double> H3st_qqbp = (2.*LFb*(6.*(C2qg[idx]*gamma1gq[idx] + 2.*C2qqp[idx]*gamma1qq[idx] + 2.*C1qq[idx]*gamma2qqp[idx] + gamma3qqp[idx]) + 6.*(C1qq[idx]*gamma1gq[idx]*gamma1qg[idx] + gamma1qg[idx]*gamma2gq[idx] + (beta0 + 2.*gamma1qq[idx])* gamma2qqp[idx])*LFb + gamma1gq[idx]*gamma1qg[idx]*(3.*beta0 + gamma1gg[idx] + 5.*gamma1qq[idx])* LFb2) - 6.*(2.*C2qg[idx]*gamma1gq[idx] + 4.*C2qqp[idx]*gamma1qq[idx] + 4.*C1qq[idx]*gamma2qqp[idx] + 2.*gamma3qqp[idx] + 2.*(2.*C1qq[idx]*gamma1gq[idx]*gamma1qg[idx] + gamma1qg[idx]*gamma2gq[idx] + gamma1gq[idx]*gamma2qg[idx] + B1q*gamma2qqp[idx] + 4.*gamma1qq[idx]*gamma2qqp[idx])*LFb + gamma1gq[idx]*gamma1qg[idx]*(B1q + beta0 + gamma1gg[idx] + 5.*gamma1qq[idx])*LFb2)*LQb + 3.*(4.*C1qq[idx]*gamma1gq[idx]*gamma1qg[idx] + 4.*gamma1qg[idx]*gamma2gq[idx] + 4.*(B1q - beta0 + 2.*gamma1qq[idx])*gamma2qqp[idx] - 2.*A1q*gamma2qqp[idx]*LFb + gamma1gq[idx]*gamma1qg[idx]*LFb*(4.*B1q - 2.*beta0 + 2.*gamma1gg[idx] + 10.*gamma1qq[idx] - A1q*LFb))* LQb2 - 2.*(-3*A1q*gamma2qqp[idx] + gamma1gq[idx]*gamma1qg[idx]*(3.*B1q - 3.*beta0 + gamma1gg[idx] + 5.*gamma1qq[idx] - 3.*A1q*LFb))*LQb3 - 3.*A1q*gamma1gq[idx]*gamma1qg[idx]*LQb4 + 6.*C1qg[idx]*LQFb*(2.*gamma2gq[idx] + gamma1gq[idx]*(2.*C1qq[idx] + (beta0 + gamma1gg[idx] + 3.*gamma1qq[idx])*LFb - (2.*B1q - 3.*beta0 + gamma1gg[idx] + 3.*gamma1qq[idx])*LQb - A1q*LQb2)))/12.;
	complex <double> H3st_qg   = (12.*C1qq[idx]*(2.*C1qg[idx]*gamma1gg[idx] + 2.*C1qq[idx]*gamma1qg[idx] + 2.*gamma2qg[idx] + gamma1qg[idx]*(gamma1gg[idx] + gamma1qq[idx])*LQFb)*LQFb + 12.*C1qg[idx]*(2.*C1qg[idx]*gamma1gq[idx] + 2.*C1qq[idx]*gamma1qq[idx] + 2.*gamma2qq[idx] + (gamma1gq[idx]*gamma1qg[idx] + pow(gamma1qq[idx],2))*LQFb)*LQFb - 6.*A1q*gamma1qg[idx]*(gamma1gg[idx] + gamma1qq[idx])*LQb4 + 12.*beta0*LQFb*((2.*gamma2qg[idx] + gamma1qg[idx]*(gamma1gg[idx] + gamma1qq[idx])*LQFb)* (LFb + LQb) + C1qg[idx]*gamma1gg[idx]*(LFb + 3.*LQb) + C1qq[idx]*gamma1qg[idx]*(LFb + 3.*LQb)) + 12.*gamma1qq[idx]*LQFb*(2.*C2qg[idx] + (2.*C1qq[idx]*gamma1qg[idx] + 2.*gamma2qg[idx] + gamma1qg[idx]*(gamma1gg[idx] + gamma1qq[idx])*LQFb)*LQFb + beta0*C1qg[idx]*(LFb + 3.*LQb) + C1qg[idx]*(2.*gamma1gg[idx]*LFb - 2.*(B1q + gamma1gg[idx])*LQb - A1q*LQb2)) + 12.*gamma1qg[idx]*LQFb* (2.*C2qq[idx] + LFb*(2.*(C1qg[idx]*gamma1gq[idx] + C1qq[idx]*gamma1qq[idx] + gamma2qq[idx]) + (gamma1gq[idx]*gamma1qg[idx] + pow(gamma1qq[idx],2))*LFb) - 2.*(C1qg[idx]*gamma1gq[idx] + gamma2qq[idx] + gamma1gq[idx]*gamma1qg[idx]*LFb + (B1q + gamma1qq[idx])*(C1qq[idx] + gamma1qq[idx]*LFb))*LQb + (gamma1gq[idx]*gamma1qg[idx] + 2.*B1q*gamma1qq[idx] + pow(gamma1qq[idx],2) - A1q*(C1qq[idx] + gamma1qq[idx]*LFb))*LQb2 + A1q*gamma1qq[idx]*LQb3 + beta0*(LFb*(C1qq[idx] + 2.*gamma1qq[idx]*LFb) + 3.*C1qq[idx]*LQb - 2.*gamma1qq[idx]*LQb2)) + gamma1qg[idx]*LQFb*(12.*beta1*(LFb + LQb) + 8.*pow(beta0,2)*(LFb2 + LFb*LQb + LQb2) - 2.*beta0*LQb*(6.*B1q*(LFb + 2.*LQb) + A1q*LQb*(3.*LFb + 7.*LQb)) + 3.*LQb*(-8*B2q + LQb*(-4*A2q + pow(2.*B1q + A1q*LQb,2)))) - 24.*LQb*(C2qg[idx]*gamma1gg[idx] + C1qg[idx]*gamma2gg[idx] + C1qq[idx]*gamma2qg[idx] + gamma3qg[idx] + gamma1qg[idx]*(C2qq[idx] + C2qqb[idx] + (C2qqbp[idx] + C2qqp[idx])*(-1 + nf))) + 4.*gamma1qg[idx]*LFb3*(pow(gamma1gg[idx],2) + gamma1gg[idx]*gamma1qq[idx] + pow(gamma1qq[idx],2) + 2.*gamma1gq[idx]*gamma1qg[idx]*nf) + 24.*LFb*(C2qg[idx]*gamma1gg[idx] + C2qq[idx]*gamma1qg[idx] + C2qqb[idx]*gamma1qg[idx] - (C2qqbp[idx] + C2qqp[idx])*gamma1qg[idx] + C1qg[idx]*gamma2gg[idx] + C1qq[idx]*gamma2qg[idx] + gamma3qg[idx] - B1q*(C1qg[idx]*gamma1gg[idx] + C1qq[idx]*gamma1qg[idx] + gamma2qg[idx])* LQb + B1q*gamma1qg[idx]*(gamma1gg[idx] + gamma1qq[idx])*LQb2 + (LQb*(-2.*(C1qg[idx]*pow(gamma1gg[idx],2) + C1qq[idx]*gamma1qg[idx]*(gamma1gg[idx] + gamma1qq[idx]) + (gamma1gg[idx] + gamma1qq[idx])*gamma2qg[idx] + gamma1qg[idx]*(gamma2gg[idx] + gamma2qq[idx] + gamma2qqb[idx] - gamma2qqbp[idx] - gamma2qqp[idx])) + gamma1qg[idx]*(pow(gamma1gg[idx],2) + gamma1gg[idx]*gamma1qq[idx] + pow(gamma1qq[idx],2))*LQb - A1q*(C1qg[idx]*gamma1gg[idx] + C1qq[idx]*gamma1qg[idx] + gamma2qg[idx])*LQb + A1q*gamma1qg[idx]*(gamma1gg[idx] + gamma1qq[idx])*LQb2))/2. + gamma1qg[idx]*(C2qqbp[idx] + C2qqp[idx] - LQb*(2.*C1qg[idx]*gamma1gq[idx] + gamma2qqbp[idx] + gamma2qqp[idx] - gamma1gq[idx]*gamma1qg[idx]*LQb))*nf) + 12.*LQb2*(C1qq[idx]*gamma1gg[idx]*gamma1qg[idx] + C1qq[idx]*gamma1qg[idx]*gamma1qq[idx] + 2.*gamma1qg[idx]*gamma2gg[idx] + 2.*gamma1qq[idx]*gamma2qg[idx] + 2.*B1q*(C1qg[idx]*gamma1gg[idx] + C1qq[idx]*gamma1qg[idx] + gamma2qg[idx]) + C1qg[idx]*(pow(gamma1gg[idx],2) + 2.*gamma1gq[idx]*gamma1qg[idx]*nf)) + 4.*LQb3*(3.*A1q*(C1qg[idx]*gamma1gg[idx] + C1qq[idx]*gamma1qg[idx] + gamma2qg[idx]) - gamma1qg[idx]*(pow(gamma1gg[idx],2) + gamma1gg[idx]*gamma1qq[idx] + pow(gamma1qq[idx],2) + 3.*B1q*(gamma1gg[idx] + gamma1qq[idx]) + 2.*gamma1gq[idx]*gamma1qg[idx]*nf)) - 6.*LFb2*(-2.*C1qq[idx]*gamma1qg[idx]*(gamma1gg[idx] + gamma1qq[idx]) - 4.*gamma1qg[idx]*gamma2gg[idx] - 4.*gamma1qq[idx]*gamma2qg[idx] - 2.*C1qg[idx]*(pow(gamma1gg[idx],2) + 2.*gamma1gq[idx]*gamma1qg[idx]*nf) + gamma1qg[idx]*LQb*(2.*(gamma1gg[idx]*(B1q + gamma1gg[idx]) + (B1q + gamma1gg[idx])*gamma1qq[idx] + pow(gamma1qq[idx],2)) + A1q*(gamma1gg[idx] + gamma1qq[idx])*LQb + 4.*gamma1gq[idx]*gamma1qg[idx]*nf)))/24.;
	complex <double> H3st_qbg  = ((2.*C2qqb[idx]*gamma1qg[idx] + (C1qg[idx] + gamma1qg[idx]*LQFb)*(2.*C1qg[idx]*gamma1gq[idx] + 2.*gamma2qqb[idx] + gamma1gq[idx]*gamma1qg[idx]*LQFb))*LQFb)/2.;
	complex <double> H3st_qpg  = ((2.*C2qqp[idx]*gamma1qg[idx] + (C1qg[idx] + gamma1qg[idx]*LQFb)*(2.*C1qg[idx]*gamma1gq[idx] + 2.*gamma2qqp[idx] + gamma1gq[idx]*gamma1qg[idx]*LQFb))*LQFb)/2.;
	complex <double> H3st_qbpg = ((2.*C2qqbp[idx]*gamma1qg[idx] + (C1qg[idx] + gamma1qg[idx]*LQFb)*(2.*C1qg[idx]*gamma1gq[idx] + 2.*gamma2qqbp[idx] + gamma1gq[idx]*gamma1qg[idx]*LQFb))*LQFb)/2.;
	complex <double> H3st_gg   = (LQFb*(4.*pow(C1qg[idx],2)*gamma1gg[idx] + 4.*C2qg[idx]*gamma1qg[idx] + 2.*C1qg[idx]*(2.*C1qq[idx]*gamma1qg[idx] + 2.*gamma2qg[idx] + gamma1qg[idx]*(beta0*LFb + 3.*gamma1gg[idx]*LFb + gamma1qq[idx]*LFb - 2.*B1q*LQb + 3.*beta0*LQb - 3.*gamma1gg[idx]*LQb - gamma1qq[idx]*LQb - A1q*LQb2)) + gamma1qg[idx]*LQFb*(4.*C1qq[idx]*gamma1qg[idx] + 4.*gamma2qg[idx] + gamma1qg[idx]*(2.*gamma1gg[idx]*LFb + 2.*gamma1qq[idx]*LFb - 2.*B1q*LQb - 2.*gamma1gg[idx]*LQb - 2.*gamma1qq[idx]*LQb - A1q*LQb2 + 2.*beta0*(LFb + LQb)))))/2.;

	qqb[i]  += as3*H3st_qqb;
	qg[i]   += as3*H3st_qg;
	qq[i]   += as3*H3st_qq;
	qqp[i]  += as3*H3st_qqp;
	qqbp[i] += as3*H3st_qqbp;
	gg[i]   += as3*H3st_gg;
	qbg[i]  += as3*H3st_qbg;
	qpg[i]  += as3*H3st_qpg;
	qbpg[i] += as3*H3st_qbpg;
	//if (opts.order == 3)
	//continue;
    }
}
