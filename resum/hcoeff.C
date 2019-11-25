#include "hcoeff.h"
#include "resconst.h"
#include "anomalous.h"
#include "mesq.h"
#include "besselint.h"
#include "resint.h"
#include "settings.h"
#include "interface.h"
#include "phasespace.h"
#include <iostream>

complex <double> *hcoeff::Hqqb;
complex <double> *hcoeff::Hqg;
complex <double> *hcoeff::Hqq;
complex <double> *hcoeff::Hqqp;
complex <double> *hcoeff::Hgg;

complex <double> hcoeff::H1stgg;
complex <double> *hcoeff::H1stqg;
complex <double> *hcoeff::H1stqqb;

complex <double> *hcoeff::H2stqq;
complex <double> *hcoeff::H2stqqp;
complex <double> *hcoeff::H2stqqb;
complex <double> *hcoeff::H2stqg;
complex <double> *hcoeff::H2stgg;

complex <double> *hcoeff::aexpqq;
complex <double> *hcoeff::aexpqg;

complex <double> hcoeff::Hqqbz;
complex <double> hcoeff::Hqgz;
complex <double> hcoeff::Hqqz;
complex <double> hcoeff::Hqqpz;
complex <double> hcoeff::Hggz;

//using namespace hcoeff;
using namespace anomalous;
using namespace resconst;

void hcoeff::allocate()
{
  if (opts.order == 0)
    return;

  //allocate memory
  //LL
  Hqqb = new complex <double> [mellinint::mdim];

  //NLL
  Hqg = new complex <double> [mellinint::mdim];

  H1stqg = new complex <double> [mellinint::mdim];
  H1stqqb = new complex <double> [mellinint::mdim];

  //NNLL
  Hqq = new complex <double> [mellinint::mdim];
  Hqqp = new complex <double> [mellinint::mdim];
  Hgg = new complex <double> [mellinint::mdim];

  H2stqq = new complex <double> [mellinint::mdim];
  H2stqqp = new complex <double> [mellinint::mdim];
  H2stqqb = new complex <double> [mellinint::mdim];
  H2stqg = new complex <double> [mellinint::mdim];
  H2stgg = new complex <double> [mellinint::mdim];

  aexpqq = new complex <double> [mellinint::mdim];
  aexpqg = new complex <double> [mellinint::mdim];
}

void hcoeff::calc(double aass, complex <double> logmuf2q2, complex <double> logq2muf2, complex <double> logq2mur2, complex <double> loga)
{
  if (opts.order == 0)
    return;

  // logs of scales are computed in resint
  //  complex <double> logmuf2q2 = resint::logmuf2q2;
  //  complex <double> logq2muf2 = resint::logq2muf2;
  //  complex <double> logq2mur2 = resint::logq2mur2;
  //  complex <double> loga = resint::loga;
  //  double aass = resint::aass;

  
  //All the following coefficients need to be calculated at each resumm iteration only with fixed factorization and renormalization scale (variable logmuf2q2) or
  //with fixed resummation scale (variable loga: a = q2/mu_res). Otherwise they can be computed at initialization
  if (opts.order == 0)
      for (int i = 0; i < mellinint::mdim; i++)
	  {
	    Hqqb[i]=1.;
	    // channels which do not contribute at LL
	    Hqg[i] = 0.;
	    Hgg[i] = 0.;
	    Hqq[i] = 0.;
	    Hqqp[i] = 0.;
	  }
  
  if (opts.order == 1)
    for (int i = 0; i < mellinint::mdim; i++)
      {
	int idx = anomalous::index(i,mesq::positive);

	H1stqqb[i] = C1QQ[idx]
	  -gamma1qq[idx]*(logmuf2q2+2.*loga)
	  +(-2.*loga)*(B1q+A1q*loga);

	Hqqb[i]=1.+aass*H1stqqb[i];

	//	    H1stqg[i] = C1QG[idx]
	//	      -gamma1qg[idx]*(logmuf2q2+2.*loga);

	// channels which do not contribute at NLL
	Hgg[i] = 0.;
	Hqq[i] = 0.;
	Hqqp[i] = 0.;
      }

  if (opts.order >= 2)
    {
      complex <double> loga2 = pow(loga,2);
      complex <double> loga3 = pow(loga,3);
      complex <double> logq2muf22 = pow(logq2muf2,2);
      H1stgg=0;
      for (int i = 0; i < mellinint::mdim; i++)
	{
	  int idx = anomalous::index(i,mesq::positive);

	  H1stqqb[i] = C1QQ[idx]
	    -gamma1qq[idx]*(logmuf2q2+2.*loga)
	    +(-2.*loga)*(B1q+A1q*loga);

	  H1stqqb[i] *= 2.;
	    
	  H1stqg[i] = C1QG[idx]
	    -gamma1qg[idx]*(logmuf2q2+2.*loga);

	  //  qq  means   Qb Qb -> Qb Q =  Q Q -> Q Qb
	  H2stqq[i] = C2NSqqbM[idx] + C2SqqbM[idx]
	    +4.*(-2.*loga*1./4.*(gamma2qqbV[idx]+gamma2qqbS[idx])
		 +1./4.*(gamma2qqbV[idx]+gamma2qqbS[idx])*logq2muf2
		 +1./2.*(H1stqg[i]+C1QG[idx])/2.
		 *1./2.*gamma1gq[idx]*(logq2muf2-2.*loga));


	  //  qqp_2  means   Q Q' -> Q Qb  flavor in sigmaQQb determined by "first parton"
	  H2stqqp[i] = C2SqqbM[idx]
	    +4.*(-2.*loga*1./4.*gamma2qqbS[idx]
		 +1./4.*gamma2qqbS[idx]*logq2muf2
		 +1./2.*(H1stqg[i]+C1QG[idx])/2.
		 *1./2.*gamma1gq[idx]*(logq2muf2-2.*loga));
	    
	  //  All *4 because of normalization (as/2pi)^2 instead of as/pi
	  // normalization of as/2pi  implies  gamma^1->gamma^1/2
	  //                                   C^1, H^1 -> C^1/2, H^1/2 
	  //                                   gamma^2-> gamma^2/4
	  //                              beta,A,B -> beta,A,B are in as/pi already

	  H2stqqb[i] = C2NSqqM[idx] + C2NSqqM[idx] + C2SqqbM[idx] + C2SqqbM[idx] + C1QQ[idx]*C1QQ[idx]
	    + 4.*(+ 1./6.*A1q*beta0*8*loga3  
		  + 1./2.*4*loga2*(A2q-beta0*(B1q+2*A1q*loga+gamma1qq[idx]/2.+gamma1qq[idx]/2.))
		  - 2.*loga*(B2q+2*A2q*loga-beta0*(C1QQ[idx]+C1QQ[idx])/2.
			     + gamma2qqV[idx]/4. + gamma2qqS[idx]/4. 
			     + gamma2qqV[idx]/4. + gamma2qqS[idx]/4. )      
		  + beta0/2.*(gamma1qq[idx]+gamma1qq[idx])/2.*logq2muf22
		  + (gamma2qqV[idx]/4. + gamma2qqS[idx]/4. + gamma2qqV[idx]/4. + gamma2qqS[idx]/4.)*logq2muf2
		  - H1stqqb[i]/2.*beta0*logq2mur2
		  + 1./2.*(H1stqqb[i]+C1QQ[idx]+C1QQ[idx])/2.*((gamma1qq[idx]+gamma1qq[idx])/2.*(logq2muf2-2.*loga)-(B1q+A1q*loga)*2.*loga)
		  + 1./4.*(H1stqg[i]+C1QG[idx])
		  * gamma1gq[idx]/2.*(logq2muf2-2.*loga)
		  + 1./4.*(H1stqg[i]+C1QG[idx])
		  * gamma1gq[idx]/2.*(logq2muf2-2.*loga));
		
	  H2stqg[i] = C2qgM[idx] + C1QG[idx]*C1QQ[idx] 
	    + 4.*(+ 1./2.*beta0*4*loga2*(-gamma1qg[idx]/2.)
		  - 2.*loga*(-beta0 * C1QG[idx]/2. + gamma2qg[idx]/4.)
		  + 1./2.*beta0*logq2muf22*gamma1qg[idx]/2.
		  + gamma2qg[idx]/4.*logq2muf2
		  - beta0*logq2mur2*H1stqg[i]/2.
		  + 1./2.*(H1stqqb[i] + C1QQ[idx] + C1QQ[idx])/2.*(logq2muf2-2.*loga)*gamma1qg[idx]/2.
		  + 1./2.*(H1stqg[i] + C1QG[idx])/2.*((logq2muf2-2.*loga)*(gamma1qq[idx]+gamma1gg[idx])/2.-(B1q+A1q*loga)*2.*loga));
		
	  //  GG done
	  H2stgg[i] = C1QG[idx]*C1QG[idx]
	    -4.*1./2.*(logmuf2q2+2.*loga)*((H1stqg[i] + C1QG[idx])/2.*gamma1qg[idx]/2.+(H1stqg[i] + C1QG[idx])/2.*gamma1qg[idx]/2.);
	}
    }
}

//b-dependent coefficients
void hcoeff::calcb(double aass, complex <double> logmuf2q2, complex <double> loga, complex <double> alpq, complex <double> aexp, complex <double> aexpb)
{
  if (opts.order == 0)
    return;

  // logs of scales are computed in resint
  //  complex <double> logmuf2q2 = resint::logmuf2q2;
  //  complex <double> loga = resint::loga;
  //  double aass = resint::aass;

  // exponents are computes in alphasl from besselint
  // complex <double> aexpb = besselint:aexpb;
  // complex <double> aexp = besselint:aexp;

  //b-dependent quantities to be computed in invres(b)
  // complex <double> alpq = pdfevol::alpq; //(alpqf * alphasl(scale2))
  //--> alpq = alphas(b0^2/b^2)

  double aassh=aass/2.;
  double aasshsq = pow(aassh,2);
  double aass2=pow(aass,2);

  complex <double> aexp2 = aexp*aexp;

  // NLL
  if (opts.order == 1)
    {
      for (int i = 0; i < mellinint::mdim; i++)
	{         
	  int idx = anomalous::index(i,mesq::positive);
	  
	  //Hqg[i] = alpq*2.*C1QG[idx]
	  //  +(aass/2.)*(-gamma1qg[idx])*(logmuf2q2+2.*loga);

	  //Bug fix in DYRES (compare lines 659-660 of DYRes-v1.0/src/Res/main2.f and line 894 of DYqT-v1.0/enew.f)
	  Hqg[i] = aexp*(aass/2.)*C1QG[idx]
	    +(aass/2.)*(-gamma1qg[idx])*(logmuf2q2+2.*loga);
	}
    }
  
  // NNLL
  //complex <double> c1delta = pow(aexpb,aassh*(-2.*C1qqn));
  if (opts.order == 2)
    {
      for (int i = 0; i < mellinint::mdim; i++)
	{
	  int idx = anomalous::index(i,mesq::positive);

	  aexpqq[i]=pow(aexpb,aassh*(C1QQ[idx]-2.*C1qqn));
	  aexpqg[i]=pow(aexpb,aassh*(C2qgM[idx]/C1QG[idx]-2.*C1qqn));

	  //aexpqq[i]=pow(aexpb,aassh*C1QQ[idx]);
	  //aexpqg[i]=pow(aexpb,aassh*(C2qgM[idx]/C1QG[idx]));

	  //QQb
	  Hqqb[i] = (1.+aassh*H1stqqb[i]+aasshsq*H2stqqb[i])
	    *aexpqq[i]
	    *aexpqq[i];

	  //QG and GQ
	  Hqg[i] = (aassh*H1stqg[i]+aasshsq*H2stqg[i])
	    *aexp
	    *aexpqg[i]
	    *aexpqq[i];

	  //QQ
	  Hqq[i]=aasshsq*H2stqq[i]
	    *aexpqq[i]
	    //*c1delta
	    *aexp2; 

	  //QQ'
	  Hqqp[i]=aasshsq*H2stqqp[i]
	    *aexpqq[i]
	    //*c1delta
	    *aexp2;
	      
	  //GG
	  Hgg[i] = aasshsq*H2stgg[i]
	    *aexp*aexpqg[i]
	    *aexp*aexpqg[i];
	}
    }

  // NNNLL
  if (opts.order == 3)
    {

      for (int i = 0; i < mellinint::mdim; i++)
	{
	  int idx = anomalous::index(i,mesq::positive);

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
	   
	  aexpqq[i]=pow(aexpb,aassh*C1QQ[idx])
	    *pow(cx(exponent_.aexpc_),aass2*(0.5*pow(C1QQ[idx],2) - (C2NSqqM[idx] + C2SqqbM[idx])))
	    *pow(cx(exponent_.aexpd_),aass2*beta1/beta0*C1QQ[idx])
	    *pow(cx(exponent_.aexpc_),aass2*beta0*C1QQ[idx]*(resint::rlogq2mur2-2.*resint::rloga));
	    ;
	  
	  //aexpqq[i]=pow(aexpb,aassh*C1QQ[idx]);
	  aexpqg[i]=pow(aexpb,aassh*(C2qgM[idx]/C1QG[idx]));
	  
	  //QQb
	  Hqqb[i] = (1.+aassh*H1stqqb[i]+aasshsq*H2stqqb[i])
	    *aexpqq[i]
	    *aexpqq[i];

	  //QG and GQ
	  Hqg[i] = (aassh*H1stqg[i]+aasshsq*H2stqg[i])
	    *aexp
	    *aexpqg[i]
	    *aexpqq[i];

	  //QQ
	  Hqq[i]=aasshsq*H2stqq[i]
	    *aexpqq[i]
	    *aexp2; 

	  //QQ'
	  Hqqp[i]=aasshsq*H2stqqp[i]
	    *aexpqq[i]
	    *aexp2;
	      
	  //GG
	  Hgg[i] = aasshsq*H2stgg[i]
	    *aexp*aexpqg[i]
	    *aexp*aexpqg[i];
	  
	  //multiple parton case as in Eqs.(100-107) of https://arxiv.org/pdf/hep-ph/0508068.pdf
	}

    }
}

void hcoeff::truncate()
{
  //Calculate truncated moments
  //double x = phasespace::m/opts.sroot;

  double x = 1e-8;//pow(phasespace::m/opts.sroot,2);

  //double x = phasespace::m/opts.sroot*exp(phasespace::ymin);

  //double x = pow(phasespace::m/opts.sroot,2);

  double lx = log(x);
  
  //truncated moments (output)
  complex <double> Hqqb_t[mellinint::mdim] = {0.};
  complex <double> Hqg_t[mellinint::mdim]  = {0.};
  complex <double> Hqq_t[mellinint::mdim]  = {0.};
  complex <double> Hqqp_t[mellinint::mdim] = {0.};
  complex <double> Hgg_t[mellinint::mdim]  = {0.};

  //cache x^(N) values
  complex <double> xn[mellinint::mdim];
  for (int n = 0; n < mellinint::mdim; n++)
    xn[n] = pow(x,mellinint::Np[n]);

  //Normalisation times Jacobian
  complex <double> facp = mellinint::CCp/2./M_PI/complex <double>(0.,1);
  complex <double> facm = mellinint::CCm/2./M_PI/complex <double>(0.,1);
  
  //original moments times prefactor and weight
  complex <double> Hqqb_tp[mellinint::mdim];
  complex <double> Hqg_tp[mellinint::mdim];
  complex <double> Hqq_tp[mellinint::mdim];
  complex <double> Hqqp_tp[mellinint::mdim];
  complex <double> Hgg_tp[mellinint::mdim];
  complex <double> Hqqb_tm[mellinint::mdim];
  complex <double> Hqg_tm[mellinint::mdim];
  complex <double> Hqq_tm[mellinint::mdim];
  complex <double> Hqqp_tm[mellinint::mdim];
  complex <double> Hgg_tm[mellinint::mdim];
  for (int m = 0; m < mellinint::mdim; m++)
    {
      Hqqb_tp[m] = facp * Hqqb[m] * mellinint::wn[m];
      Hqg_tp[m]  = facp * Hqg[m]  * mellinint::wn[m];
      Hqq_tp[m]  = facp * Hqq[m]  * mellinint::wn[m];
      Hqqp_tp[m] = facp * Hqqp[m] * mellinint::wn[m];
      Hgg_tp[m]  = facp * Hgg[m]  * mellinint::wn[m];

      Hqqb_tm[m] = facm * conj(Hqqb[m]) * mellinint::wn[m];
      Hqg_tm[m]  = facm * conj(Hqg[m] ) * mellinint::wn[m];
      Hqq_tm[m]  = facm * conj(Hqq[m] ) * mellinint::wn[m];
      Hqqp_tm[m] = facm * conj(Hqqp[m]) * mellinint::wn[m];
      Hgg_tm[m]  = facm * conj(Hgg[m] ) * mellinint::wn[m];
    }

  //cache factor (1-x^(N-M))/(N-M) which limit is ln(x) when N-M -> 0
  complex <double> llxp[mellinint::mdim][mellinint::mdim];
  complex <double> llxm[mellinint::mdim][mellinint::mdim];
  for (int n = 0; n < mellinint::mdim; n++)
    for (int m = 0; m < mellinint::mdim; m++)
      {
	llxp[n][m] = (1.-xn[n]/xn[m])/(mellinint::Np[n]-mellinint::Np[m]);
	llxm[n][m] = (1.-xn[n]/conj(xn[m]))/(mellinint::Np[n]-mellinint::Nm[m]);
      }

  //overwrite divergent diagonal part
  for (int n = 0; n < mellinint::mdim; n++)
    llxp[n][n] = -lx;
  
  for (int n = 0; n < mellinint::mdim; n++)
    for (int m = 0; m < mellinint::mdim; m++)
      {
	Hqqb_t[m] += Hqqb_tp[m]*llxp[n][m] - Hqqb_tm[m]*llxm[n][m];
	Hqg_t[m]  += Hqg_tp[m] *llxp[n][m] - Hqg_tm[m] *llxm[n][m];
	Hqq_t[m]  += Hqq_tp[m] *llxp[n][m] - Hqq_tm[m] *llxm[n][m];
	Hqqp_t[m] += Hqqp_tp[m]*llxp[n][m] - Hqqp_tm[m]*llxm[n][m];
	Hgg_t[m]  += Hgg_tp[m] *llxp[n][m] - Hgg_tm[m] *llxm[n][m];
      }


  //replace moments
  for (int n = 0; n < mellinint::mdim; n++)
    {
      Hqqb[n] = Hqqb_t[n];
      Hqg[n]  = Hqg_t[n] ;
      Hqq[n]  = Hqq_t[n] ;
      Hqqp[n] = Hqqp_t[n];
      Hgg[n]  = Hgg_t[n] ;
    }
}

void hcoeff::invert(double q2)
{
  Hqqbz = 1;
  Hqgz = 0;
  Hggz = 0;
  Hqqz = 0;
  Hqqpz = 0;

  if (opts.order == 0)
    return;

  double bjx= q2/pow(opts.sroot,2);
  double ax = log(bjx);

  for (int i = 0; i < mellinint::mdim; i++)
    {
      complex <double> cexp = exp(-mellinint::Np[i] * ax)/M_PI * mellinint::CCp/complex <double>(0.,1);
      //complex <double> cexm = exp(-mellinint::Nm[i] * ax)/M_PI * mellinint::CCm/complex <double>(0.,1);
      
      //contribution starting at LL
      Hqqbz += real(Hqqb[i]*cexp*mellinint::wn[i]);
      //Hqqbz += real(cexp*mellinint::wn[i]);// * 2*sqrt(q2)/pow(opts.sroot,2);
      
      //contribution starting at NLL
      Hqgz += real(Hqg[i]*cexp*mellinint::wn[i]);
    
      //contributions starting at NNLL
      Hggz += real(Hgg[i]*cexp*mellinint::wn[i]);
      Hqqz += real(Hqq[i]*cexp*mellinint::wn[i]);
      Hqqpz += real(Hqqp[i]*cexp*mellinint::wn[i]);
    }

  //cout << "hcoeff::invert(q2) " << sqrt(q2) << " " << Hqqbz << endl;
}

void hcoeff::free()
{
  if (opts.order == 0)
    return;

  delete[] Hqqb;
  delete[] Hqg;
  delete[] Hqq;
  delete[] Hqqp;
  delete[] Hgg;
  delete[] H1stqg;
  delete[] H1stqqb;
  delete[] H2stqq;
  delete[] H2stqqp;
  delete[] H2stqqb;
  delete[] H2stqg;
  delete[] H2stgg;
  delete[] aexpqq;
  delete[] aexpqg;
}
