#include "hcoefficients.h"
#include "resconst.h"
#include "anomalous.h"
#include "mesq.h"
#include "settings.h"
#include "interface.h"
#include <iostream>

const complex <double> H1q=(0.,0.);

complex <double> *hcoefficients::Hqqb;
complex <double> *hcoefficients::Hqg;
complex <double> *hcoefficients::Hqg_1;
complex <double> *hcoefficients::Hqg_2;
complex <double> *hcoefficients::Hqq_nnll;
complex <double> *hcoefficients::Hqq;
complex <double> *hcoefficients::Hqq_1;
complex <double> *hcoefficients::Hqq_2;
complex <double> *hcoefficients::Hqqp;
complex <double> *hcoefficients::Hqqp_1;
complex <double> *hcoefficients::Hqqp_2;
complex <double> *hcoefficients::Hgg;

complex <double> hcoefficients::H1stgg;
complex <double> *hcoefficients::H1stqg;
complex <double> *hcoefficients::H1stqqb;

complex <double> *hcoefficients::H2stqq;
complex <double> *hcoefficients::H2stqqp;
complex <double> *hcoefficients::H2stqqb;
complex <double> *hcoefficients::H2stqg_1;
complex <double> *hcoefficients::H2stqg_2;
complex <double> *hcoefficients::H2stgg;

complex <double> *hcoefficients::aexpqq;
complex <double> *hcoefficients::aexpqg;

//using namespace hcoefficients;
//#pragma omp threadprivate(Hqqb,Hqg,Hqg_1,Hqg_2,Hqq_nnll,Hqq,Hqq_1,Hqq_2,Hqqp,Hqqp_1,Hqqp_2,Hgg,H1stgg,H1stqg,H1stqqb,H2stqq,H2stqqp,H2stqqb,H2stqg_1,H2stqg_2,H2stgg,aexpqq,aexpqg)

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
  if (opts.order == 0)
    return;

  //allocate memory
  //LL
  Hqqb = new complex <double> [mellinint::mdim*mellinint::mdim*2];

  //NLL
  Hqg = new complex <double> [mellinint::mdim*2];
  Hqg_1 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
  Hqg_2 = new complex <double> [mellinint::mdim*mellinint::mdim*2];

  H1stqg = new complex <double> [mellinint::mdim*2];
  H1stqqb = new complex <double> [mellinint::mdim*mellinint::mdim*2];

  //NNLL
  Hqq_nnll = new complex <double> [mellinint::mdim*2];
  Hqq = new complex <double> [mellinint::mdim*mellinint::mdim*2];
  Hqq_1 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
  Hqq_2 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
  Hqqp = new complex <double> [mellinint::mdim*2];
  Hqqp_1 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
  Hqqp_2 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
  Hgg = new complex <double> [mellinint::mdim*mellinint::mdim*2];

  H2stqq = new complex <double> [mellinint::mdim*2];
  H2stqqp = new complex <double> [mellinint::mdim*2];
  H2stqqb = new complex <double> [mellinint::mdim*mellinint::mdim*2];
  H2stqg_1 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
  H2stqg_2 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
  H2stgg = new complex <double> [mellinint::mdim*mellinint::mdim*2];

  aexpqq = new complex <double> [mellinint::mdim*2];
  aexpqg = new complex <double> [mellinint::mdim*2];
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

  
  //All the following coefficients need to be calculated at each resumm iteration only with fixed factorization and renormalization scale (variable logmuf2q2) or
  //with fixed resummation scale (variable loga: a = q2/mu_res). Otherwise they can be computed at initialization
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
  
  if (opts.order == 1)
    for (int sign = mesq::positive; sign <= mesq::negative; sign++)
      for (int i1 = 0; i1 < mellinint::mdim; i1++)
	for (int i2 = 0; i2 < mellinint::mdim; i2++)
	  {
	    Hqqb[index(i1,i2,sign)]=1.+aass/2.*
	      (H1q + anomalous::C1QQ[anomalous::index(i1,mesq::positive)] + anomalous::C1QQ[anomalous::index(i2,sign)])
	      -aass/2.*(anomalous::gamma1qq[anomalous::index(i1,0)]+anomalous::gamma1qq[anomalous::index(i2,sign)])*(logmuf2q2+2.*loga)
	      +aass/2.*(-4.*loga)*(resconst::B1q+resconst::A1q*loga);
	    
	    /*
	    H1stqqb[index(i1,i2,sign)]=
	      (H1q+anomalous::C1QQ[anomalous::index(i1,mesq::positive)]+anomalous::C1QQ[anomalous::index(i2,sign)])
	      -(anomalous::gamma1qq[anomalous::index(i1,mesq::positive)]+anomalous::gamma1qq[anomalous::index(i2,sign)])*(logmuf2q2+2.*loga)
	      +(-4.*loga)*(resconst::B1q+resconst::A1q*loga);

	    Hqqb[index(i1,i2,sign)]=1.+aass/2.*H1stqqb[index(i1,i2,sign)];
	    */
	    
	    // channels which do not contribute at NLL
	    Hgg[index(i1,i2,sign)] = 0.;//
	    Hqq[index(i1,i2,sign)] = 0.;//              !  Q Q -> Qb Q  = Qb Qb -> Q Qb
	    Hqq_1[index(i1,i2,sign)] = 0.;//
	    Hqq_2[index(i1,i2,sign)] = 0.;//
	    Hqqp_1[index(i1,i2,sign)] = 0.;//           !  qqp_1  means   Q' Q -> Qb Q  flavor in sigmaQQb determined by "second parton"
	    Hqqp_2[index(i1,i2,sign)] = 0.;//		!  qqp_2  means   Q Q' -> Q Qb  flavor in sigmaQQb determined by "first parton" 
	  }

  if (opts.order == 2)
    {
      complex <double> loga2 = pow(loga,2);
      complex <double> loga3 = pow(loga,3);
      complex <double> logq2muf22 = pow(logq2muf2,2);
      H1stgg=0;
      for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	for (int i = 0; i < mellinint::mdim; i++)
	  {
	    H1stqg[index(i,sign)] = anomalous::C1QG[anomalous::index(i,sign)]
	      -anomalous::gamma1qg[anomalous::index(i,sign)]*(logmuf2q2+2.*loga);


	    //  qq  means   Qb Qb -> Qb Q =  Q Q -> Q Qb
	    H2stqq[index(i,sign)]= anomalous::C2NSqqbM[anomalous::index(i,sign)] + anomalous::C2SqqbM[anomalous::index(i,sign)];
	    H2stqq[index(i,sign)]=H2stqq[index(i,sign)]
	      +4.*(-2.*loga*1./4.*(anomalous::gamma2qqbV[anomalous::index(i,sign)]+anomalous::gamma2qqbS[anomalous::index(i,sign)])
		   +1./4.*(anomalous::gamma2qqbV[anomalous::index(i,sign)]+anomalous::gamma2qqbS[anomalous::index(i,sign)])*logq2muf2
		   +1./2.*(H1stqg[index(i,sign)]+anomalous::C1QG[anomalous::index(i,sign)])/2.
		   *1./2.*anomalous::gamma1gq[anomalous::index(i,sign)]*(logq2muf2-2.*loga));


	    //  qqp_2  means   Q Q' -> Q Qb  flavor in sigmaQQb determined by "first parton"
	    H2stqqp[index(i,sign)]= anomalous::C2SqqbM[anomalous::index(i,sign)];
	    H2stqqp[index(i,sign)]=H2stqqp[index(i,sign)]
	      +4.*(-2.*loga*1./4.*anomalous::gamma2qqbS[anomalous::index(i,sign)]
		   +1./4.*anomalous::gamma2qqbS[anomalous::index(i,sign)]*logq2muf2
		   +1./2.*(H1stqg[index(i,sign)]+anomalous::C1QG[anomalous::index(i,sign)])/2.
		   *1./2.*anomalous::gamma1gq[anomalous::index(i,sign)]*(logq2muf2-2.*loga));
	  }

      for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	for (int i1 = 0; i1 < mellinint::mdim; i1++)
	  for (int i2 = 0; i2 < mellinint::mdim; i2++)
	    {
	      H1stqqb[index(i1,i2,sign)]=(H1q+anomalous::C1QQ[anomalous::index(i1,0)]+anomalous::C1QQ[anomalous::index(i2,sign)])
		-(anomalous::gamma1qq[anomalous::index(i1,0)]+anomalous::gamma1qq[anomalous::index(i2,sign)])*(logmuf2q2+2.*loga)
		+(-4.*loga)*(resconst::B1q+resconst::A1q*loga);

	      //  All *4 because of normalization (as/2pi)^2 instead of as/pi
	      // normalization of as/2pi  implies  gamma^1->gamma^1/2
	      //                                   C^1, H^1 -> C^1/2, H^1/2 
	      //                                   gamma^2-> gamma^2/4
	      //                              beta,A,B -> beta,A,B are in as/pi already

	      H2stqqb[index(i1,i2,sign)] = anomalous::C2NSqqM[anomalous::index(i1,0)] + anomalous::C2NSqqM[anomalous::index(i2,sign)] + anomalous::C2SqqbM[anomalous::index(i1,0)] + anomalous::C2SqqbM[anomalous::index(i2,sign)] + anomalous::C1QQ[anomalous::index(i1,0)]*anomalous::C1QQ[anomalous::index(i2,sign)];
	      H2stqqb[index(i1,i2,sign)] = H2stqqb[index(i1,i2,sign)]
		+ 4.*(+ 1./6.*resconst::A1q*resconst::beta0*8*loga3  
		      + 1./2.*4*loga2*(resconst::A2q-resconst::beta0*(resconst::B1q+2*resconst::A1q*loga+anomalous::gamma1qq[anomalous::index(i1,0)]/2.+anomalous::gamma1qq[anomalous::index(i2,sign)]/2.))
		      - 2.*loga*(resconst::B2q+2*resconst::A2q*loga-resconst::beta0*(anomalous::C1QQ[anomalous::index(i1,0)]+anomalous::C1QQ[anomalous::index(i2,sign)])/2.
				+ anomalous::gamma2qqV[anomalous::index(i1,0)]/4. + anomalous::gamma2qqS[anomalous::index(i1,0)]/4. 
				+ anomalous::gamma2qqV[anomalous::index(i2,sign)]/4. + anomalous::gamma2qqS[anomalous::index(i2,sign)]/4. )      
		      + resconst::beta0/2.*(anomalous::gamma1qq[anomalous::index(i1,0)]+anomalous::gamma1qq[anomalous::index(i2,sign)])/2.*logq2muf22
		      + (anomalous::gamma2qqV[anomalous::index(i1,0)]/4. + anomalous::gamma2qqS[anomalous::index(i1,0)]/4. + anomalous::gamma2qqV[anomalous::index(i2,sign)]/4. + anomalous::gamma2qqS[anomalous::index(i2,sign)]/4.)*logq2muf2
		      - H1stqqb[index(i1,i2,sign)]/2.*resconst::beta0*logq2mur2
		      + 1./2.*(H1stqqb[index(i1,i2,sign)]+H1q+anomalous::C1QQ[anomalous::index(i1,0)]+anomalous::C1QQ[anomalous::index(i2,sign)])/2.*((anomalous::gamma1qq[anomalous::index(i1,0)]+anomalous::gamma1qq[anomalous::index(i2,sign)])/2.*(logq2muf2-2.*loga)-(resconst::B1q+resconst::A1q*loga)*2.*loga)
		      + 1./4.*(H1stqg[index(i1,0)]+anomalous::C1QG[anomalous::index(i1,0)])
		      * anomalous::gamma1gq[anomalous::index(i1,0)]/2.*(logq2muf2-2.*loga)
		      + 1./4.*(H1stqg[index(i2,sign)]+anomalous::C1QG[anomalous::index(i2,sign)])
		      * anomalous::gamma1gq[anomalous::index(i2,sign)]/2.*(logq2muf2-2.*loga));
		
	      H2stqg_1[index(i1,i2,sign)] = anomalous::C2qgM[anomalous::index(i1,0)] + anomalous::C1QG[anomalous::index(i1,0)]*anomalous::C1QQ[anomalous::index(i2,sign)] 
		+ 4.*(+ 1./2.*resconst::beta0*4*loga2*(-anomalous::gamma1qg[anomalous::index(i1,0)]/2.)
		     - 2.*loga*(-resconst::beta0 * anomalous::C1QG[anomalous::index(i1,0)]/2. + anomalous::gamma2qg[anomalous::index(i1,0)]/4.)
		     + 1./2.*resconst::beta0*logq2muf22*anomalous::gamma1qg[anomalous::index(i1,0)]/2.
		     + anomalous::gamma2qg[anomalous::index(i1,0)]/4.*logq2muf2
		     - resconst::beta0*logq2mur2*H1stqg[index(i1,0)]/2.
		     + 1./2.*(H1stqqb[index(i1,i2,sign)] + H1q + anomalous::C1QQ[anomalous::index(i1,0)] + anomalous::C1QQ[anomalous::index(i2,sign)])/2.*(logq2muf2-2.*loga)*anomalous::gamma1qg[anomalous::index(i1,0)]/2.
		     + 1./2.*(H1stqg[index(i1,0)] + anomalous::C1QG[anomalous::index(i1,0)])/2.*((logq2muf2-2.*loga)*(anomalous::gamma1qq[anomalous::index(i2,sign)]+anomalous::gamma1gg[anomalous::index(i1,0)])/2.-(resconst::B1q+resconst::A1q*loga)*2.*loga));
		
	      H2stqg_2[index(i1,i2,sign)] = anomalous::C2qgM[anomalous::index(i2,sign)] + anomalous::C1QG[anomalous::index(i2,sign)]*anomalous::C1QQ[anomalous::index(i1,0)]
		+ 4.*(+ 1./2.*resconst::beta0*4*loga2*(-anomalous::gamma1qg[anomalous::index(i2,sign)]/2.)
		     - 2.*loga*(-resconst::beta0 * anomalous::C1QG[anomalous::index(i2,sign)]/2. + anomalous::gamma2qg[anomalous::index(i2,sign)]/4.)
		     + 1./2.*resconst::beta0*logq2muf22*(anomalous::gamma1qg[anomalous::index(i2,sign)]/2.)
		     + anomalous::gamma2qg[anomalous::index(i2,sign)]/4.*logq2muf2
		     - resconst::beta0*logq2mur2*H1stqg[index(i2,sign)]/2.
		     + 1./2.*(H1stqqb[index(i1,i2,sign)] + H1q + anomalous::C1QQ[anomalous::index(i2,sign)] + anomalous::C1QQ[anomalous::index(i1,0)])/2.*(logq2muf2-2.*loga)*anomalous::gamma1qg[anomalous::index(i2,sign)]/2.
		     + 1./2.*(H1stqg[index(i2,sign)] + anomalous::C1QG[anomalous::index(i2,sign)])/2. * ((logq2muf2-2.*loga)*(anomalous::gamma1qq[anomalous::index(i1,0)] + anomalous::gamma1gg[anomalous::index(i2,sign)])/2.-(resconst::B1q+resconst::A1q*loga)*2.*loga));

	      //!!! simplification
	      //H2stqg_2[index(i1,i2,sign)] = H2stqg_1[index(i2,i1,0)]
	      //!!!
	      
	      //  GG done
	      H2stgg[index(i1,i2,sign)] = anomalous::C1QG[anomalous::index(i2,sign)]*anomalous::C1QG[anomalous::index(i1,0)]
		-4.*1./2.*(logmuf2q2+2.*loga)*((H1stqg[index(i1,0)] + anomalous::C1QG[anomalous::index(i1,0)])/2.*anomalous::gamma1qg[anomalous::index(i2,sign)]/2.+(H1stqg[index(i2,sign)] + anomalous::C1QG[anomalous::index(i2,sign)])/2.*anomalous::gamma1qg[anomalous::index(i1,0)]/2.);
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
}
void hcoefficients::free()
{
  if (opts.order == 0)
    return;

  delete[] Hqqb;
  delete[] Hqg;
  delete[] Hqg_1;
  delete[] Hqg_2;
  delete[] Hqq_nnll;
  delete[] Hqq;
  delete[] Hqq_1;
  delete[] Hqq_2;
  delete[] Hqqp;
  delete[] Hqqp_1;
  delete[] Hqqp_2;
  delete[] Hgg;
  delete[] H1stqg;
  delete[] H1stqqb;
  delete[] H2stqq;
  delete[] H2stqqp;
  delete[] H2stqqb;
  delete[] H2stqg_1;
  delete[] H2stqg_2;
  delete[] H2stgg;
  delete[] aexpqq;
  delete[] aexpqg;
}
