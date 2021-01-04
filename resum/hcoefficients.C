#include "hcoefficients.h"
#include "resconst.h"
#include "anomalous.h"
#include "resint.h"
#include "mesq.h"
#include "settings.h"
#include "interface.h"
#include <iostream>

const complex <double> H1q=(0.,0.);

complex <double> *hcoefficients::Hqqb;
complex <double> *hcoefficients::Hqg_1;
complex <double> *hcoefficients::Hqg_2;
complex <double> *hcoefficients::Hqq_1;
complex <double> *hcoefficients::Hqq_2;
complex <double> *hcoefficients::Hqqp_1;
complex <double> *hcoefficients::Hqqp_2;
complex <double> *hcoefficients::Hgg;

//complex <double> *hcoefficients::H1stqqb;
//complex <double> *hcoefficients::H1stqg;
//
//complex <double> *hcoefficients::H2stqqb;
//complex <double> *hcoefficients::H2stqg;
//complex <double> *hcoefficients::H2stqq;
//complex <double> *hcoefficients::H2stqqp;
//
//complex <double> *hcoefficients::aexpqq;
//complex <double> *hcoefficients::aexpqg;

using namespace anomalous;
using namespace resconst;

//using namespace hcoefficients;

//fortran interface
void hcoeff_calc_(double& aass, double& logmuf2q2, double& logq2muf2, double& logq2mur2, double& loga)
{
  hcoefficients::calc(aass, logmuf2q2, logq2muf2, logq2mur2, loga);
};
void hcoeff_calcb_(double& aass, double& logmuf2q2, double& loga, double& alpq, double &aexp, double &aexpb)
{
  hcoefficients::calcb(aass, logmuf2q2, loga, alpq, aexp, aexpb);
};

void hcoefficients::allocate()
{
  //allocate memory
  //LL
  Hqqb = new complex <double> [mellinint::mdim*mellinint::mdim*2];

  if (opts.order == 0)
    return;

  //NLL
  Hqg_1 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
  Hqg_2 = new complex <double> [mellinint::mdim*mellinint::mdim*2];

  //H1stqg = new complex <double> [mellinint::mdim*2];
  //H1stqqb = new complex <double> [mellinint::mdim*mellinint::mdim*2];

  //if (opts.order == 1)
  //return;

  //NNLL
  Hqq_1  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
  Hqq_2  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
  Hqqp_1 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
  Hqqp_2 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
  Hgg    = new complex <double> [mellinint::mdim*mellinint::mdim*2];

//  H2stqq = new complex <double> [mellinint::mdim*2];
//  H2stqqp = new complex <double> [mellinint::mdim*2];
//  H2stqqb = new complex <double> [mellinint::mdim*mellinint::mdim*2];
//  H2stqg_1 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
//  H2stqg_2 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
//  H2stgg = new complex <double> [mellinint::mdim*mellinint::mdim*2];
//
//  aexpqq = new complex <double> [mellinint::mdim*2];
//  aexpqg = new complex <double> [mellinint::mdim*2];
}

void hcoefficients::reset()
{
  fill(Hqqb,   Hqqb  +mellinint::mdim*mellinint::mdim*2, 1.);

  if (opts.order == 0)
    return;

  fill(Hqg_1,  Hqg_1 +mellinint::mdim*mellinint::mdim*2, 0.);
  fill(Hqg_2,  Hqg_2 +mellinint::mdim*mellinint::mdim*2, 0.);

  //if (opts.order == 1)
  //return;

  fill(Hqq_1,  Hqq_1 +mellinint::mdim*mellinint::mdim*2, 0.);
  fill(Hqq_2,  Hqq_2 +mellinint::mdim*mellinint::mdim*2, 0.);
  fill(Hqqp_1, Hqqp_1+mellinint::mdim*mellinint::mdim*2, 0.);
  fill(Hqqp_2, Hqqp_2+mellinint::mdim*mellinint::mdim*2, 0.);
  fill(Hgg ,   Hgg   +mellinint::mdim*mellinint::mdim*2, 0.);
}

void hcoefficients::calc(double aass, complex <double> logmuf2q2, complex <double> logq2muf2, complex <double> logq2mur2, complex <double> loga)
{
  if (opts.order == 0)
    return;

  // logs of scales are computed in resint
  //  complex <double> logmuf2q2 = resint::logmuf2q2;
  //  complex <double> logq2muf2 = resint::logq2muf2;
  //  complex <double> logq2mur2 = resint::logq2mur2;
  //  complex <double> loga = resint::loga;
  //  double aass = resint::aass;


  complex <double> LQF = logq2muf2-2.*loga;
  complex <double> LQ2F2 = pow(logq2muf2,2)-4.*pow(loga,2);
  if (opts.mufevol)
    {
      //LQF = 0.;
      //LQ2F2 = 0.;
      LQF = resint::LF;
      LQ2F2 = pow(resint::LF,2);
    }
  LQF = 0.;
  LQ2F2 = 0.;

  double aassh=aass/2.;
  double aasshsq = pow(aassh,2);
  
  
  //All the following coefficients need to be calculated at each resumm iteration only with fixed factorization and renormalization scale (variable logmuf2q2) or
  //with fixed resummation scale (variable loga: a = q2/mu_res). Otherwise they can be computed at initialization
  /*
  if (opts.order == 0)
    for (int sign = mesq::positive; sign <= mesq::negative; sign++)
      for (int i1 = 0; i1 < mellinint::mdim; i1++)
	for (int i2 = 0; i2 < mellinint::mdim; i2++)
	  {
	    Hqqb[index(i1,i2,sign)]=1.;
	    // channels which do not contribute at LL
	    Hqg_1[index(i1,i2,sign)] = 0.;
	    Hqg_2[index(i1,i2,sign)] = 0.;
	    Hgg[index(i1,i2,sign)] = 0.;//
	    Hqq[index(i1,i2,sign)] = 0.;//              !  Q Q -> Qb Q  = Qb Qb -> Q Qb
	    Hqq_1[index(i1,i2,sign)] = 0.;//
	    Hqq_2[index(i1,i2,sign)] = 0.;//
	    Hqqp_1[index(i1,i2,sign)] = 0.;//           !  qqp_1  means   Q' Q -> Qb Q  flavor in sigmaQQb determined by "second parton"
	    Hqqp_2[index(i1,i2,sign)] = 0.;//		!  qqp_2  means   Q Q' -> Q Qb  flavor in sigmaQQb determined by "first parton" 
	  }
  */

  complex <double> H1stqqb[mellinint::mdim*2];
  complex <double> H1stqg[mellinint::mdim*2];
  
  if (opts.order == 1)
    {
      for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	for (int i = 0; i < mellinint::mdim; i++)
	  {
	    int idx = anomalous::index(i,sign);
	    int ii = hcoefficients::index(i,sign);

	    H1stqqb[ii] = C1QQ[idx]
	      -gamma1qq[idx]*(-LQF)
	      +(-2.*loga)*(B1q+A1q*loga);

	    //Bug fix in DYRES (compare lines 659-660 of DYRes-v1.0/src/Res/main2.f and line 894 of DYqT-v1.0/enew.f)
	    H1stqg[ii] = C1QG[idx]
	      -gamma1qg[idx]*(-LQF);
	  }

      for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	for (int i1 = 0; i1 < mellinint::mdim; i1++)
	  for (int i2 = 0; i2 < mellinint::mdim; i2++)
	    {
	      //Hqqb[index(i1,i2,sign)]=1.+aass/2.*
	      //(H1q + anomalous::C1QQ[anomalous::index(i1,mesq::positive)] + anomalous::C1QQ[anomalous::index(i2,sign)])
	      //-aass/2.*(anomalous::gamma1qq[anomalous::index(i1,0)]+anomalous::gamma1qq[anomalous::index(i2,sign)])*(logmuf2q2+2.*loga)
	      //+aass/2.*(-4.*loga)*(resconst::B1q+resconst::A1q*loga);

	      int ii1 = index(i1,0);
	      int ii2 = index(i2,sign);
	      int idx = index(i1,i2,sign);
	      
	      Hqqb[idx] = 1. + aass/2.*(H1q+H1stqqb[ii1]+H1stqqb[ii2]);
  
	      Hqg_1[idx] = aass/2.*H1stqg[ii1];
	      Hqg_2[idx] = aass/2.*H1stqg[ii2];
	      
	      
	      /*
		H1stqqb[index(i1,i2,sign)]=
		(H1q+anomalous::C1QQ[anomalous::index(i1,mesq::positive)]+anomalous::C1QQ[anomalous::index(i2,sign)])
		-(anomalous::gamma1qq[anomalous::index(i1,mesq::positive)]+anomalous::gamma1qq[anomalous::index(i2,sign)])*(logmuf2q2+2.*loga)
		+(-4.*loga)*(resconst::B1q+resconst::A1q*loga);
		
		Hqqb[index(i1,i2,sign)]=1.+aass/2.*H1stqqb[index(i1,i2,sign)];
	      */
	    
	      //// channels which do not contribute at NLL
	      //Hgg[index(i1,i2,sign)] = 0.;//
	      //Hqq[index(i1,i2,sign)] = 0.;//              !  Q Q -> Qb Q  = Qb Qb -> Q Qb
	      //Hqq_1[index(i1,i2,sign)] = 0.;//
	      //Hqq_2[index(i1,i2,sign)] = 0.;//
	      //Hqqp_1[index(i1,i2,sign)] = 0.;//           !  qqp_1  means   Q' Q -> Qb Q  flavor in sigmaQQb determined by "second parton"
	      //Hqqp_2[index(i1,i2,sign)] = 0.;//		!  qqp_2  means   Q Q' -> Q Qb  flavor in sigmaQQb determined by "first parton" 
	    }
    }

  if (opts.order == 2)
    {
      complex <double> H1stqqb[mellinint::mdim*2];
      complex <double> H1stqg[mellinint::mdim*2];

      complex <double> H2stqqb[mellinint::mdim*2];
      complex <double> H2stqg[mellinint::mdim*2];
      complex <double> H2stqq[mellinint::mdim*2];
      complex <double> H2stqqp[mellinint::mdim*2];
      
      complex <double> loga2 = pow(loga,2);
      complex <double> loga3 = pow(loga,3);
      complex <double> logq2muf22 = pow(logq2muf2,2);
      for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	for (int i = 0; i < mellinint::mdim; i++)
	  {
	    int idx = anomalous::index(i,sign);
	    int ii = hcoefficients::index(i,sign);

	    H1stqqb[ii] = C1QQ[idx]
	      -gamma1qq[idx]*(-LQF)
	      +(-2.*loga)*(B1q+A1q*loga);

	    H1stqg[ii] = C1QG[idx]
	      -gamma1qg[idx]*(-LQF);

	    // All *4 because of normalization (as/2pi)^2 instead of as/pi
	    // normalization of as/2pi  implies  gamma^1->gamma^1/2
	    //                                   C^1, H^1 -> C^1/2, H^1/2 
	    //                                   gamma^2-> gamma^2/4
	    //                              beta,A,B -> beta,A,B are in as/pi already

	    //qqb
	    H2stqqb[ii] = C2NSqqM[idx]+C2SqqbM[idx] 
	      + 4.*(+ 0.5*1./6.*A1q*beta0*8*loga3  
		    + 0.5*1./2.*4*loga2*(A2q-beta0*(B1q+2*A1q*loga))
		    - 0.5*2.*loga*(B2q+2*A2q*loga-beta0*C1QQ[idx])
		    + 0.5*(gamma2qqV[idx]+gamma2qqS[idx])/2.*LQF
		    + 0.5*1./2.*beta0*gamma1qq[idx]*LQ2F2
		    - 0.5*H1stqqb[idx]*beta0*logq2mur2
		    + 0.5*1./2.*(H1stqqb[idx]+C1QQ[idx])*(-(B1q+A1q*loga)*2.*loga)
		    + 0.5*1./4.*(H1stqg[ii]+C1QG[idx])*gamma1gq[idx]*(LQF));

	    //qg
	    H2stqg[ii] = C2qgM[idx]
	      + 4.*(- 2.*loga*(-beta0 * C1QG[idx]/2.)
		    + 1./2.*beta0*gamma1qg[idx]/2.*LQ2F2
		    + gamma2qg[idx]/4.*LQF
		    - beta0*logq2mur2*H1stqg[ii]/2.
		    + 1./2.*(H1stqqb[ii]+C1QQ[idx])/2.*(LQF)*gamma1qg[idx]/2.
		    + 1./2.*(H1stqg[ii]+C1QG[idx])/2.*((LQF)*(gamma1gg[idx])/2.-(B1q+A1q*loga)*2.*loga)
		    );

	    //qq (Qb Qb -> Qb Q =  Q Q -> Q Qb)
	    H2stqq[ii] = C2NSqqbM[idx] + C2SqqbM[idx]
	      +4.*(+1./4.*(gamma2qqbV[idx]+gamma2qqbS[idx])*LQF
		   +1./2.*(H1stqg[ii]+C1QG[idx])/2.
		   *1./2.*gamma1gq[idx]*(LQF));


	    //qqp (Q Q' -> Q Qb)
	    H2stqqp[ii] = C2SqqbM[idx]
	      +4.*(+1./4.*gamma2qqbS[idx]*LQF
		   +1./2.*(H1stqg[ii]+C1QG[idx])/2.
		   *1./2.*gamma1gq[idx]*(LQF));
	  }

      for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	for (int i1 = 0; i1 < mellinint::mdim; i1++)
	  for (int i2 = 0; i2 < mellinint::mdim; i2++)
	    {
	      int idx1 = anomalous::index(i1,0);
	      int idx2 = anomalous::index(i2,sign);
	      int ii1 = hcoefficients::index(i1,0);
	      int ii2 = hcoefficients::index(i2,sign);
	      int ii12 = hcoefficients::index(i1,i2,sign);

	      complex <double> H2stqqb_mix = C1QQ[idx1]*C1QQ[idx2]
		+ 4.*(1./2.*(H1stqqb[ii1]+H1stqqb[ii2]+C1QQ[idx1]+C1QQ[idx2])/2.*(gamma1qq[idx1]+gamma1qq[idx2])/2.*(LQF));

	      complex <double> H2stqg_mix_1 = C1QG[idx1]*C1QQ[idx2]
		+4.*(+ 1./2.*(C1QQ[idx2])/2.*(LQF)*gamma1qg[idx1]/2.
		     + 1./2.*(H1stqg[ii1] + C1QG[idx1])/2.*((LQF)*(gamma1qq[idx2])/2.)
		     )
		;
	    
	      complex <double> H2stqg_mix_2 = C1QG[idx2]*C1QQ[idx1]
		+4.*(+ 1./2.*(C1QQ[idx1])/2.*(LQF)*gamma1qg[idx2]/2.
		     + 1./2.*(H1stqg[ii2] + C1QG[idx2])/2.*((LQF)*(gamma1qq[idx1])/2.)
		     )
		;
	      
	      complex <double> H2stgg_mix = C1QG[idx2]*C1QG[idx1]
		-4.*1./2.*(-LQF)*((H1stqg[ii1]+C1QG[idx1])/2.*gamma1qg[idx2]/2.
				  +(H1stqg[ii2]+C1QG[idx2])/2.*gamma1qg[idx1]/2.);
		
	      //QQb
	      Hqqb[ii12] = 1.+aassh*(H1stqqb[ii1]+H1stqqb[ii2]) + aasshsq*(H2stqqb[ii1]+H2stqqb[ii2]+H2stqqb_mix);

	      //qg_1 -> GQ initial state
	      Hqg_1[ii12] = aassh*H1stqg[ii1] + aasshsq*(H2stqg[ii1]+H2stqg_mix_1);

	      //qg_2 -> QG initial state
	      Hqg_2[ii12] = aassh*H1stqg[ii2] + aasshsq*(H2stqg[ii2]+H2stqg_mix_2);
	      
	      //QQ
	      Hqq_1[ii12] = aasshsq*H2stqq[ii1]; // Q Q -> Qb Q and Qb Qb -> Q Qb (sigmaQQb determined by "second parton")
	      Hqq_2[ii12] = aasshsq*H2stqq[ii2]; // Q Q -> Q Qb and Qb Qb -> Qb Q (sigmaQQb determined by "first parton")

	      //QQ'
	      Hqqp_1[ii12] = aasshsq*H2stqqp[ii1]; // qqp_1  means   Q' Q -> Qb Q  flavor in sigmaQQb determined by "second parton"
	      Hqqp_2[ii12] = aasshsq*H2stqqp[ii2]; // qqp_2  means   Q Q' -> Q Qb  flavor in sigmaQQb determined by "first parton"
	      
	      //GG
	      Hgg[ii12] = aasshsq*H2stgg_mix;
	    }
    }
}

//b-dependent coefficients
void hcoefficients::calcb(double aass, complex <double> logmuf2q2, complex <double> loga, complex <double> alpq, complex <double> aexp, complex <double> aexpb)
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
  /*
  double aassh=aass/2.;
  double aasshsq = pow(aassh,2);

  complex <double> aexp2 = aexp*aexp;

  double dC1qqn=2.*M_PI*M_PI/3.-16./3.; //= 2.*resconst::C1qqn;

  // NLL
  if (opts.order == 1)
    {
      for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	for (int i = 0; i < mellinint::mdim; i++)
	  {         
	    //Hqg[index(i,sign)] = alpq*2.*anomalous::C1QG[anomalous::index(i,sign)]
	    //  +(aass/2.)*(-anomalous::gamma1qg[anomalous::index(i,sign)])*(logmuf2q2+2.*loga);

	    //Bug fix in DYRES (compare lines 659-660 of DYRes-v1.0/src/Res/main2.f and line 894 of DYqT-v1.0/enew.f)
	    Hqg[index(i,sign)] = aexp*(aass/2.)*anomalous::C1QG[anomalous::index(i,sign)]
	      +(aass/2.)*(-anomalous::gamma1qg[anomalous::index(i,sign)])*(logmuf2q2+2.*loga);
	  }

      for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	for (int i1 = 0; i1 < mellinint::mdim; i1++)
	  for (int i2 = 0; i2 < mellinint::mdim; i2++)
	    {
	      Hqg_1[index(i1,i2,sign)]=Hqg[index(i1,0)];
	      Hqg_2[index(i1,i2,sign)]=Hqg[index(i2,sign)];
	    }
    }
  // NNLL
  if (opts.order == 2)
    {
      for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	for (int i = 0; i < mellinint::mdim; i++)
	  {
	    //precompute the CPU expensive powers
	    aexpqq[index(i,sign)]=pow(aexpb,aassh*(anomalous::C1QQ[anomalous::index(i,sign)]-dC1qqn));
	    aexpqg[index(i,sign)]=pow(aexpb,aassh*(anomalous::C2qgM[anomalous::index(i,sign)]/anomalous::C1QG[anomalous::index(i,sign)]-dC1qqn));
	    
	    Hqq_nnll[index(i,sign)]=aasshsq*H2stqq[index(i,sign)]*aexpqq[index(i,sign)]*aexp2;
	    Hqqp[index(i,sign)]=aasshsq*H2stqqp[index(i,sign)]*aexpqq[index(i,sign)]*aexp2;
	  }

      for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	for (int i1 = 0; i1 < mellinint::mdim; i1++)
	  for (int i2 = 0; i2 < mellinint::mdim; i2++)
	    {
	      int idx12 = index(i1,i2,sign);
	      int idx1 = index(i1,mesq::positive);
	      int idx2 = index(i2,sign);

	      //QQb
	      Hqqb[idx12] = (1.+aassh*H1stqqb[idx12]+aasshsq*H2stqqb[idx12])
		*aexpqq[idx1]           //   !leg q
		*aexpqq[idx2];         // !leg qb

		//  qg_1 means GQ initial state
	      Hqg_1[idx12] = (aassh*H1stqg[idx1]+aasshsq*H2stqg_1[idx12])
		*aexp
		*aexpqg[idx1]
		*aexpqq[idx1];

	      //  qg_2 means QG initial state
	      Hqg_2[idx12] = (aassh*H1stqg[idx2]+aasshsq*H2stqg_2[idx12])
		*aexp
		*aexpqg[idx2]
		*aexpqq[idx2];
	      
	      //  GG
	      Hgg[idx12] = aasshsq*H2stgg[idx12]
		*aexp*aexpqg[idx1]
		*aexp*aexpqg[idx2];

	      Hqq[idx12] = (Hqq_nnll[idx1]+Hqq_nnll[idx2])/2.;// Average QQ->QQb  and QbQb->QQb

	      Hqq_1[idx12] = Hqq_nnll[idx1];    // Q Q -> Qb Q and Qb Qb -> Q Qb (sigmaQQb determined by "second parton")
	      Hqq_2[idx12] = Hqq_nnll[idx2]; // Q Q -> Q Qb and Qb Qb -> Qb Q (sigmaQQb determined by "first parton")
	      
	      Hqqp_1[idx12] = Hqqp[idx1];    // qqp_1  means   Q' Q -> Qb Q  flavor in sigmaQQb determined by "second parton"
	      Hqqp_2[idx12] = Hqqp[idx2]; // qqp_2  means   Q Q' -> Q Qb  flavor in sigmaQQb determined by "first parton"
	    }
    }
*/
}
void hcoefficients::free()
{
  delete[] Hqqb;

  if (opts.order == 0)
    return;

  delete[] Hqg_1;
  delete[] Hqg_2;
  delete[] Hqq_1;
  delete[] Hqq_2;
  delete[] Hqqp_1;
  delete[] Hqqp_2;
  delete[] Hgg;
}
