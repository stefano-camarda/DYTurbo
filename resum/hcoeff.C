#include "hcoeff.h"
#include "resconst.h"
#include "anomalous.h"
#include "ccoeff.h"
#include "mesq.h"
#include "besselint.h"
#include "resint.h"
#include "settings.h"
#include "interface.h"
#include "phasespace.h"
#include "parton.h"
#include "pmom.h"
#include "mellinint.h"
#include <ctime>

#include <iostream>
#include <map>

//LL
complex <double> *hcoeff::Hqqb;
//NLL
complex <double> *hcoeff::Hqg;
complex <double> *hcoeff::Hqg_1;
complex <double> *hcoeff::Hqg_2;
//NNLL
complex <double> *hcoeff::Hqq;
complex <double> *hcoeff::Hqq_1;
complex <double> *hcoeff::Hqq_2;
complex <double> *hcoeff::Hqqp;
complex <double> *hcoeff::Hqqp_1;
complex <double> *hcoeff::Hqqp_2;
complex <double> *hcoeff::Hqqbp;
complex <double> *hcoeff::Hqqbp_1;
complex <double> *hcoeff::Hqqbp_2;
complex <double> *hcoeff::Hgg;
//NNNLL
complex <double> *hcoeff::Hqbg;
complex <double> *hcoeff::Hqbg_1;
complex <double> *hcoeff::Hqbg_2;
complex <double> *hcoeff::Hqpg;
complex <double> *hcoeff::Hqpg_1;
complex <double> *hcoeff::Hqpg_2;
complex <double> *hcoeff::Hqbpg;
complex <double> *hcoeff::Hqbpg_1;
complex <double> *hcoeff::Hqbpg_2;


//complex <double> hcoeff::H1stgg;
//complex <double> *hcoeff::H1stqg;
//complex <double> *hcoeff::H1stqqb;
//
//complex <double> *hcoeff::H2stqq;
//complex <double> *hcoeff::H2stqqp;
//complex <double> *hcoeff::H2stqqb;
//complex <double> *hcoeff::H2stqg;
//complex <double> *hcoeff::H2stgg;
//
//complex <double> *hcoeff::aexpqq;
//complex <double> *hcoeff::aexpqg;

complex <double> hcoeff::Hqqbz;
complex <double> hcoeff::Hqgz;
complex <double> hcoeff::Hqqz;
complex <double> hcoeff::Hqqpz;
complex <double> hcoeff::Hggz;

//using namespace hcoeff;
using namespace anomalous;
using namespace resconst;
using namespace ccoeff;
using namespace resint;
using namespace parton;

//allocate memory
void hcoeff::allocate()
{
  if (opts.order == 0)
    return;

  if (opts.mellin1d)
    {
      Hqqb  = new complex <double> [mellinint::mdim];
      Hqg   = new complex <double> [mellinint::mdim];
      Hqq   = new complex <double> [mellinint::mdim];
      Hqqp  = new complex <double> [mellinint::mdim];
      Hqqbp = new complex <double> [mellinint::mdim];
      Hgg   = new complex <double> [mellinint::mdim];
      Hqbg  = new complex <double> [mellinint::mdim];
      Hqpg  = new complex <double> [mellinint::mdim];
      Hqbpg = new complex <double> [mellinint::mdim];
    }
  else
    {
      Hqqb    = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      Hqg_1   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      Hqg_2   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      Hqq_1   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      Hqq_2   = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      Hqqp_1  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      Hqqp_2  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      Hqqbp_1 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      Hqqbp_2 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      Hgg     = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      Hqbg_1  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      Hqbg_2  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      Hqpg_1  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      Hqpg_2  = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      Hqbpg_1 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
      Hqbpg_2 = new complex <double> [mellinint::mdim*mellinint::mdim*2];
    }

  //cleanup
  
  //H1stqg = new complex <double> [mellinint::mdim];
  //H1stqqb = new complex <double> [mellinint::mdim];
  //
  //H2stqq = new complex <double> [mellinint::mdim];
  //H2stqqp = new complex <double> [mellinint::mdim];
  //H2stqqb = new complex <double> [mellinint::mdim];
  //H2stqg = new complex <double> [mellinint::mdim];
  //H2stgg = new complex <double> [mellinint::mdim];
  //
  //aexpqq = new complex <double> [mellinint::mdim];
  //aexpqg = new complex <double> [mellinint::mdim];
}


void hcoeff::free()
{
  if (opts.order == 0)
    return;

  if (opts.mellin1d)
    {
      delete[] Hqqb;
      delete[] Hqg;
      delete[] Hqq;
      delete[] Hqqp;  
      delete[] Hqqbp;
      delete[] Hgg;
      delete[] Hqbg;
      delete[] Hqpg;
      delete[] Hqbpg;
    }
  else
    {
      delete[] Hqqb;
      delete[] Hqg_1;
      delete[] Hqg_2;
      delete[] Hqq_1;
      delete[] Hqq_2;
      delete[] Hqqp_1;
      delete[] Hqqp_2;
      delete[] Hqqbp_1;
      delete[] Hqqbp_2;
      delete[] Hgg;
      delete[] Hqbg_1;
      delete[] Hqbg_2;
      delete[] Hqpg_1;
      delete[] Hqpg_2;
      delete[] Hqbpg_1;
      delete[] Hqbpg_2;
    }
}

void hcoeff::reset()
{
  if (opts.order == 0)
    return;

  if (opts.mellin1d)
    {
      fill(Hqqb , Hqqb  + mellinint::mdim, 0.);
      fill(Hqg  , Hqg   + mellinint::mdim, 0.);
      fill(Hqq  , Hqq   + mellinint::mdim, 0.);
      fill(Hqqp , Hqqp  + mellinint::mdim, 0.);
      fill(Hqqbp, Hqqbp + mellinint::mdim, 0.);
      fill(Hgg  , Hgg   + mellinint::mdim, 0.);
      fill(Hqbg , Hqbg  + mellinint::mdim, 0.);
      fill(Hqpg , Hqpg  + mellinint::mdim, 0.);
      fill(Hqbpg, Hqbpg + mellinint::mdim, 0.);
    }
  else
    {
      fill(Hqqb   ,Hqqb    + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(Hqg_1  ,Hqg_1   + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(Hqg_2  ,Hqg_2   + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(Hqq_1  ,Hqq_1   + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(Hqq_2  ,Hqq_2   + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(Hqqp_1 ,Hqqp_1  + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(Hqqp_2 ,Hqqp_2  + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(Hqqbp_1,Hqqbp_1 + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(Hqqbp_2,Hqqbp_2 + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(Hgg    ,Hgg     + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(Hqbg_1 ,Hqbg_1  + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(Hqbg_2 ,Hqbg_2  + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(Hqpg_1 ,Hqpg_1  + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(Hqpg_2 ,Hqpg_2  + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(Hqbpg_1,Hqbpg_1 + mellinint::mdim*mellinint::mdim*2, 0.);
      fill(Hqbpg_2,Hqbpg_2 + mellinint::mdim*mellinint::mdim*2, 0.);
    }
}

void hcoeff::calc()
{
  if (opts.order == 0)
    return;

  /**************************************************/
  //LO
  if (opts.mellin1d)
    fill(Hqqb, Hqqb+mellinint::mdim, 1.);
  else 
    fill(Hqqb, Hqqb+mellinint::mdim*mellinint::mdim*2, 1.);
  
  if (opts.order == 0)
    return;
  /**************************************************/

  //complex <double> LQ = 2.*resint::loga;
  double LQ2 = pow(LQ,2);
  double LQ3 = pow(LQ,3);
  double LQ4 = pow(LQ,4);
  double LQ5 = pow(LQ,5);
  double LQ6 = pow(LQ,6);

  //complex <double> LR = resint::logq2mur2;
  double LR2 = pow(LR,2);

  double LF2 = pow(LF,2);
  double LF3 = pow(LF,3);
  //LF = 0.;
  //LF2 = 0.;
  //LF3 = 0.;
  //beta0 = 0.;
  //B1q = 0.;
  //A1q = 0.;
  //B2q = 0.;
  //A2q = 0.;

  double LQF  = LF - LQ;
  double LQF2 = pow(LQF,2);
  double LQF3 = pow(LQF,3);

  double as = aass;
  double as2 = pow(aass,2);
  double as3 = pow(aass,3);
  
  if (opts.mellin1d)
    //rapidity integrated
    for (int i = 0; i < mellinint::mdim; i++)
      {
	int idx = anomalous::index(i,mesq::positive);

	//if (opts.mufevol) //should be actually if (opts.mufevol || opts.evolmode == 2 || opts.evolmode == 3)
	//LQF = 0.;
	//LQF = resint::LF;
      
	//NLO
	complex <double> H1st_qqb = 2.*C1qq[idx];
	complex <double> H1st_qg = C1qg[idx];

	//Resummation scale variations
	H1st_qqb += -B1q*LQ -A1q/2.*pow(LQ,2);

	//Factorization scale variations
	if (!opts.mufevol)
	  {
	    H1st_qqb += 2.*pmom::gamma1qq[idx]*LQF;
	    H1st_qg  += pmom::gamma1qg[idx]*LQF;
	  }
      
	//cout << "c3new.txt and hst-channels.txt " << endl;
	//cout << "H1st_qqb " << H1st_qqb << endl;
	//cout << "H1st_qg  " << H1st_qg  << endl;

	Hqqb[i] += as*H1st_qqb;
	Hqg[i]  += as*H1st_qg;

	if (opts.order == 1)
	  continue;

	//NNLO
	complex <double> H2st_qqb  = 2.*C2qq[idx] + C1qq[idx]*C1qq[idx];
	complex <double> H2st_qg   = C2qg[idx] + C1qq[idx]*C1qg[idx];
	complex <double> H2st_qq   = C2qqb[idx];
	complex <double> H2st_qqp  = C2qqbp[idx];
	complex <double> H2st_qqbp = C2qqp[idx]; // (=H2st_qqp)
	complex <double> H2st_gg   = C1qg[idx]*C1qg[idx];
	  
	//Resummation scale variations
	H2st_qqb += (-B2q*LQ +(-A2q/2.+pow(B1q,2)/2.-B1q*beta0/2.)*LQ2 + A1q*(B1q/2. - beta0/3.)*LQ3 + pow(A1q,2)*LQ4/8.) + ((beta0 - B1q)*LQ - A1q/2.*LQ2)*2.*C1qq[idx];
	H2st_qg  += ((beta0 - B1q)*LQ -A1q/2.*LQ2)*C1qg[idx];

	//Renormalization scale variations
	H2st_qqb += -beta0*H1st_qqb*LR;
	H2st_qg  += -beta0*H1st_qg*LR;

	//Factorization scale variations
	if (!opts.mufevol)
	  {
	    H2st_qqb  += LQF*(2.*C1qg[idx]*pmom::gamma1gq[idx] + 2.*pmom::gamma2qq[idx] + pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LQF + pmom::gamma1qq[idx]*(4.*C1qq[idx] + 2.*pmom::gamma1qq[idx]*LF - 2.*(B1q + pmom::gamma1qq[idx])*LQ - A1q*LQ2 + beta0*(LF + LQ)));
	    H2st_qq   += ((2.*C1qg[idx]*pmom::gamma1gq[idx] + 2.*pmom::gamma2qqb[idx] + pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LQF)*LQF)/2.;
	    H2st_qqp  += ((2.*C1qg[idx]*pmom::gamma1gq[idx] + 2.*pmom::gamma2qqp[idx] + pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LQF)*LQF)/2.;
	    H2st_qqbp += ((2.*C1qg[idx]*pmom::gamma1gq[idx] + 2.*pmom::gamma2qqp[idx] + pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LQF)*LQF)/2.;
	    H2st_qg   += (LQF*(4.*C1qq[idx]*pmom::gamma1qg[idx] + 2.*C1qg[idx]*(pmom::gamma1gg[idx] + pmom::gamma1qq[idx]) + 2.*pmom::gamma2qg[idx] + beta0*pmom::gamma1qg[idx]*LF + pmom::gamma1gg[idx]*pmom::gamma1qg[idx]*LF + 3.*pmom::gamma1qg[idx]*pmom::gamma1qq[idx]*LF - 2.*B1q*pmom::gamma1qg[idx]*LQ + beta0*pmom::gamma1qg[idx]*LQ - pmom::gamma1gg[idx]*pmom::gamma1qg[idx]*LQ - 3.*pmom::gamma1qg[idx]*pmom::gamma1qq[idx]*LQ - A1q*pmom::gamma1qg[idx]*LQ2))/2.;
	    H2st_gg   += pmom::gamma1qg[idx]*(2.*C1qg[idx] + pmom::gamma1qg[idx]*LQF)*LQF;
	  }
      
	//cout << endl;
	//cout << "c3new.txt and hst-channels.txt " << endl;
	//cout << "H2st_qqb  " << H2st_qqb  << endl;
	//cout << "H2st_qg   " << H2st_qg   << endl;
	//cout << "H2st_qq   " << H2st_qq   << endl;
	//cout << "H2st_qqp  " << H2st_qqp  << endl;
	//cout << "H2st_qqbp " << H2st_qqbp << endl;
	//cout << "H2st_gg   " << H2st_gg   << endl;
      
	Hqqb[i]  += as2*H2st_qqb;
	Hqg[i]   += as2*H2st_qg;
	Hqq[i]   += as2*H2st_qq;
	Hqqp[i]  += as2*H2st_qqp;
	Hqqbp[i] += as2*H2st_qqbp;
	Hgg[i]   += as2*H2st_gg;

	if (opts.order == 2)
	  continue;

	//NNNLO
	complex <double> H3st_qqb  = 2.*C3qq[idx] + 2.*(C1qq[idx]*C2qq[idx]);
	complex <double> H3st_qg   = C3qg[idx] + (C1qq[idx]*C2qg[idx]) + (C2qq[idx] * C1qg[idx]);
	complex <double> H3st_qq   = C3qqb[idx] + (C1qq[idx]*C2qqb[idx]);
	complex <double> H3st_qqp  = C3qqbp[idx] + (C1qq[idx]*C2qqbp[idx]); // (=C3qqbp+ [C1qq * C2qqp]);
	complex <double> H3st_qqbp = C3qqp[idx] + (C1qq[idx]*C2qqp[idx]);
	complex <double> H3st_gg   = 2.*(C1qg[idx]*C2qg[idx]);
	complex <double> H3st_qbg  = (C2qqb[idx]*C1qg[idx]);
	complex <double> H3st_qpg  = (C2qqp[idx]*C1qg[idx]);
	complex <double> H3st_qbpg = (C2qqbp[idx]*C1qg[idx]);

	//Resummation scale variations
	H3st_qqb += (-B3q*LQ + (-A3q/2. + B1q*B2q - B2q*beta0 - B1q*beta1/2.)*LQ2
		     + (A2q*B1q/2. - pow(B1q,3)/6. + A1q*B2q/2. - 2.*A2q*beta0/3. + pow(B1q,2)*beta0/2. - B1q*pow(beta0,2)/3. - A1q*beta1/3.)*LQ3
		     + (A1q*A2q/4. - A1q*pow(B1q,2)/4. + 7*A1q*B1q*beta0/12. - A1q*pow(beta0,2)/4.)*LQ4 + (-pow(A1q,2)*B1q/8. + pow(A1q,2)*beta0/6.)*LQ5 - pow(A1q,3)*LQ6/48.)
	  + ((beta1 - B2q)*LQ + (-A2q/2. + pow(B1q,2)/2. - 3./2.*B1q*beta0 + pow(beta0,2))*LQ2 + A1q*(B1q/2. - 5*beta0/6.)*LQ3 + pow(A1q,2)*LQ4/8.)*2.*C1qq[idx]
	  + ((-B1q + 2.*beta0)*LQ - A1q/2.*LQ2)*(2.*C2qq[idx] + (C1qq[idx]*C1qq[idx]));

	H3st_qg += ((beta1 - B2q)*LQ + (-A2q/2. + pow(B1q,2)/2. - 3./2.*B1q*beta0 + pow(beta0,2))*LQ2 + A1q*(B1q/2. - 5*beta0/6.)*LQ3 + pow(A1q,2)*LQ4/8.)*C1qg[idx]
	  + ((-B1q + 2.*beta0)*LQ - A1q/2.*LQ2)*(C2qg[idx] + (C1qq[idx] * C1qg[idx]));
      
	H3st_qq += ((-B1q + 2.*beta0)*LQ - A1q/2.*LQ2)*C2qqb[idx];
	H3st_qqp +=  ((-B1q + 2.*beta0)*LQ - A1q/2.*LQ2)*C2qqbp[idx];
	H3st_qqbp +=  ((-B1q + 2*beta0)*LQ - A1q/2.*LQ2)*C2qqp[idx];

	H3st_gg += ((-B1q + 2.*beta0)*LQ - A1q/2.*LQ2)*(C1qg[idx] * C1qg[idx]);

	H3st_qbg  +=  ((-B1q + 2*beta0)*LQ - A1q/2.*LQ2)*(C1qqb[idx] * C1qg[idx]);
	H3st_qpg  +=  ((-B1q + 2*beta0)*LQ - A1q/2.*LQ2)*(C1qqp[idx] * C1qg[idx]);
	H3st_qbpg +=  ((-B1q + 2*beta0)*LQ - A1q/2.*LQ2)*(C1qqbp[idx] * C1qg[idx]);

	//Renormalization scale variations
	H3st_qqb  += -2*beta0*H2st_qqb*LR + (-beta1*LR - pow(beta0,2)*pow(LR,2))*H1st_qqb;
	H3st_qg   += -2*beta0*H2st_qg*LR  + (-beta1*LR - pow(beta0,2)*pow(LR,2))*H1st_qg;
	H3st_qq   += -2*beta0*H2st_qq*LR;
	H3st_qqp  += -2*beta0*H2st_qqp*LR;
	H3st_qqbp += -2*beta0*H2st_qqbp*LR;
	H3st_gg   += -2*beta0*H2st_gg*LR;

	//Factorization scale variations
	if (!opts.mufevol)
	  {
	    double nf = 5.;
	    H3st_qqb  += 2.*pow(C1qq[idx],2)*pmom::gamma1qq[idx]*LF + 4.*C2qq[idx]*pmom::gamma1qq[idx]*LF + 4.*C1qq[idx]*pmom::gamma2qq[idx]*LF + 2.*pmom::gamma3qq[idx]*LF + 2.*C1qq[idx]*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LF2 + beta1*pmom::gamma1qq[idx]*LF2 + 2.*beta0*C1qq[idx]*pmom::gamma1qq[idx]*LF2 + 4.*C1qq[idx]*pow(pmom::gamma1qq[idx],2)*LF2 + 2.*pmom::gamma1qg[idx]*pmom::gamma2gq[idx]*LF2 + 2.*beta0*pmom::gamma2qq[idx]*LF2 + 4.*pmom::gamma1qq[idx]*pmom::gamma2qq[idx]*LF2 + beta0*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LF3 + (pmom::gamma1gg[idx]*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LF3)/3. + (2.*pow(beta0,2)*pmom::gamma1qq[idx]*LF3)/3. + (5.*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*pmom::gamma1qq[idx]*LF3)/3. + 2.*beta0*pow(pmom::gamma1qq[idx],2)*LF3 + (4.*pow(pmom::gamma1qq[idx],3)*LF3)/3. + 2.*C2qg[idx]*pmom::gamma1gq[idx]*LQF - 2.*pow(C1qq[idx],2)*pmom::gamma1qq[idx]*LQ - 4.*C2qq[idx]*pmom::gamma1qq[idx]*LQ - 4.*C1qq[idx]*pmom::gamma2qq[idx]*LQ - 2.*pmom::gamma3qq[idx]*LQ - 4.*C1qq[idx]*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LF*LQ - 2.*B2q*pmom::gamma1qq[idx]*LF*LQ - 4.*B1q*C1qq[idx]*pmom::gamma1qq[idx]*LF*LQ + 4.*beta0*C1qq[idx]*pmom::gamma1qq[idx]*LF*LQ - 8.*C1qq[idx]*pow(pmom::gamma1qq[idx],2)*LF*LQ - 2.*pmom::gamma1qg[idx]*pmom::gamma2gq[idx]*LF*LQ - 2.*pmom::gamma1gq[idx]*pmom::gamma2qg[idx]*LF*LQ - 2.*B1q*pmom::gamma2qq[idx]*LF*LQ - 8.*pmom::gamma1qq[idx]*pmom::gamma2qq[idx]*LF*LQ - B1q*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LF2*LQ - beta0*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LF2*LQ - pmom::gamma1gg[idx]*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LF2*LQ - B1q*beta0*pmom::gamma1qq[idx]*LF2*LQ - 5.*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*pmom::gamma1qq[idx]*LF2*LQ - 2.*B1q*pow(pmom::gamma1qq[idx],2)*LF2*LQ - 2.*beta0*pow(pmom::gamma1qq[idx],2)*LF2*LQ - 4.*pow(pmom::gamma1qq[idx],3)*LF2*LQ + 2.*C1qq[idx]*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LQ2 + 2.*B2q*pmom::gamma1qq[idx]*LQ2 - beta1*pmom::gamma1qq[idx]*LQ2 + 4.*B1q*C1qq[idx]*pmom::gamma1qq[idx]*LQ2 - 6.*beta0*C1qq[idx]*pmom::gamma1qq[idx]*LQ2 + 4.*C1qq[idx]*pow(pmom::gamma1qq[idx],2)*LQ2 + 2.*pmom::gamma1qg[idx]*pmom::gamma2gq[idx]*LQ2 + 2.*B1q*pmom::gamma2qq[idx]*LQ2 - 2.*beta0*pmom::gamma2qq[idx]*LQ2 + 4.*pmom::gamma1qq[idx]*pmom::gamma2qq[idx]*LQ2 + 2.*B1q*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LF*LQ2 - beta0*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LF*LQ2 + pmom::gamma1gg[idx]*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LF*LQ2 - A2q*pmom::gamma1qq[idx]*LF*LQ2 + pow(B1q,2)*pmom::gamma1qq[idx]*LF*LQ2 - B1q*beta0*pmom::gamma1qq[idx]*LF*LQ2 - 2.*A1q*C1qq[idx]*pmom::gamma1qq[idx]*LF*LQ2 + 5.*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*pmom::gamma1qq[idx]*LF*LQ2 + 4.*B1q*pow(pmom::gamma1qq[idx],2)*LF*LQ2 - 2.*beta0*pow(pmom::gamma1qq[idx],2)*LF*LQ2 + 4.*pow(pmom::gamma1qq[idx],3)*LF*LQ2 - A1q*pmom::gamma2qq[idx]*LF*LQ2 - (A1q*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LF2*LQ2)/2. - (A1q*beta0*pmom::gamma1qq[idx]*LF2*LQ2)/2. - A1q*pow(pmom::gamma1qq[idx],2)*LF2*LQ2 - B1q*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LQ3 + beta0*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LQ3 - (pmom::gamma1gg[idx]*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LQ3)/3. + A2q*pmom::gamma1qq[idx]*LQ3 - pow(B1q,2)*pmom::gamma1qq[idx]*LQ3 + 2.*B1q*beta0*pmom::gamma1qq[idx]*LQ3 - (2.*pow(beta0,2)*pmom::gamma1qq[idx]*LQ3)/3. + 2.*A1q*C1qq[idx]*pmom::gamma1qq[idx]*LQ3 - (5.*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*pmom::gamma1qq[idx]*LQ3)/3. - 2.*B1q*pow(pmom::gamma1qq[idx],2)*LQ3 + 2.*beta0*pow(pmom::gamma1qq[idx],2)*LQ3 - (4.*pow(pmom::gamma1qq[idx],3)*LQ3)/3. + A1q*pmom::gamma2qq[idx]*LQ3 + A1q*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LF*LQ3 + A1q*B1q*pmom::gamma1qq[idx]*LF*LQ3 - (2.*A1q*beta0*pmom::gamma1qq[idx]*LF*LQ3)/3. + 2.*A1q*pow(pmom::gamma1qq[idx],2)*LF*LQ3 - (A1q*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LQ4)/2. - A1q*B1q*pmom::gamma1qq[idx]*LQ4 + (7.*A1q*beta0*pmom::gamma1qq[idx]*LQ4)/6. - A1q*pow(pmom::gamma1qq[idx],2)*LQ4 + (pow(A1q,2)*pmom::gamma1qq[idx]*LF*LQ4)/4. - (pow(A1q,2)*pmom::gamma1qq[idx]*LQ5)/4. + C1qg[idx]*LQF*(2.*C1qq[idx]*pmom::gamma1gq[idx] + 2.*pmom::gamma2gq[idx] + pmom::gamma1gq[idx]*(beta0*LF + pmom::gamma1gg[idx]*LF + 3.*pmom::gamma1qq[idx]*LF - 2.*B1q*LQ + 3.*beta0*LQ - pmom::gamma1gg[idx]*LQ - 3.*pmom::gamma1qq[idx]*LQ - A1q*LQ2));
	    H3st_qq   += (2.*LF*(6.*(C2qg[idx]*pmom::gamma1gq[idx] + 2.*C2qqb[idx]*pmom::gamma1qq[idx] + 2.*C1qq[idx]*pmom::gamma2qqb[idx] + pmom::gamma3qqb[idx]) + 6.*(C1qq[idx]*pmom::gamma1gq[idx]*pmom::gamma1qg[idx] + pmom::gamma1qg[idx]*pmom::gamma2gq[idx] + (beta0 + 2.*pmom::gamma1qq[idx])* pmom::gamma2qqb[idx])*LF + pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*(3.*beta0 + pmom::gamma1gg[idx] + 5.*pmom::gamma1qq[idx])* LF2) - 6.*(2.*C2qg[idx]*pmom::gamma1gq[idx] + 4.*C2qqb[idx]*pmom::gamma1qq[idx] + 4.*C1qq[idx]*pmom::gamma2qqb[idx] + 2.*pmom::gamma3qqb[idx] + 2.*(2.*C1qq[idx]*pmom::gamma1gq[idx]*pmom::gamma1qg[idx] + pmom::gamma1qg[idx]*pmom::gamma2gq[idx] + pmom::gamma1gq[idx]*pmom::gamma2qg[idx] + B1q*pmom::gamma2qqb[idx] + 4.*pmom::gamma1qq[idx]*pmom::gamma2qqb[idx])*LF + pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*(B1q + beta0 + pmom::gamma1gg[idx] + 5.*pmom::gamma1qq[idx])*LF2)*LQ + 3.*(4.*C1qq[idx]*pmom::gamma1gq[idx]*pmom::gamma1qg[idx] + 4.*pmom::gamma1qg[idx]*pmom::gamma2gq[idx] + 4.*(B1q - beta0 + 2.*pmom::gamma1qq[idx])*pmom::gamma2qqb[idx] - 2.*A1q*pmom::gamma2qqb[idx]*LF + pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LF*(4.*B1q - 2.*beta0 + 2.*pmom::gamma1gg[idx] + 10.*pmom::gamma1qq[idx] - A1q*LF))* LQ2 - 2.*(-3*A1q*pmom::gamma2qqb[idx] + pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*(3.*B1q - 3.*beta0 + pmom::gamma1gg[idx] + 5.*pmom::gamma1qq[idx] - 3.*A1q*LF))*LQ3 - 3.*A1q*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LQ4 + 6.*C1qg[idx]*LQF*(2.*pmom::gamma2gq[idx] + pmom::gamma1gq[idx]*(2.*C1qq[idx] + (beta0 + pmom::gamma1gg[idx] + 3.*pmom::gamma1qq[idx])*LF - (2.*B1q - 3.*beta0 + pmom::gamma1gg[idx] + 3.*pmom::gamma1qq[idx])*LQ - A1q*LQ2)))/12.;
	    H3st_qqp  += (2.*LF*(6.*(C2qg[idx]*pmom::gamma1gq[idx] + 2.*C2qqbp[idx]*pmom::gamma1qq[idx] + 2.*C1qq[idx]*pmom::gamma2qqbp[idx] + pmom::gamma3qqbp[idx]) + 6.*(C1qq[idx]*pmom::gamma1gq[idx]*pmom::gamma1qg[idx] + pmom::gamma1qg[idx]*pmom::gamma2gq[idx] + (beta0 + 2.*pmom::gamma1qq[idx])*pmom::gamma2qqbp[idx])*LF + pmom::gamma1gq[idx]*pmom::gamma1qg[idx]* (3.*beta0 + pmom::gamma1gg[idx] + 5.*pmom::gamma1qq[idx])*LF2) - 6.*(2.*C2qg[idx]*pmom::gamma1gq[idx] + 4.*C2qqbp[idx]*pmom::gamma1qq[idx] + 4.*C1qq[idx]*pmom::gamma2qqbp[idx] + 2.*pmom::gamma3qqbp[idx] + 2.*(2.*C1qq[idx]*pmom::gamma1gq[idx]*pmom::gamma1qg[idx] + pmom::gamma1qg[idx]*pmom::gamma2gq[idx] + pmom::gamma1gq[idx]*pmom::gamma2qg[idx] + B1q*pmom::gamma2qqbp[idx] + 4.*pmom::gamma1qq[idx]*pmom::gamma2qqbp[idx])*LF + pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*(B1q + beta0 + pmom::gamma1gg[idx] + 5.*pmom::gamma1qq[idx])*LF2)*LQ + 3.*(4.*C1qq[idx]*pmom::gamma1gq[idx]*pmom::gamma1qg[idx] + 4.*pmom::gamma1qg[idx]*pmom::gamma2gq[idx] + 4.*(B1q - beta0 + 2.*pmom::gamma1qq[idx])*pmom::gamma2qqbp[idx] - 2.*A1q*pmom::gamma2qqbp[idx]*LF + pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LF*(4.*B1q - 2.*beta0 + 2.*pmom::gamma1gg[idx] + 10.*pmom::gamma1qq[idx] - A1q*LF))* LQ2 - 2.*(-3*A1q*pmom::gamma2qqbp[idx] + pmom::gamma1gq[idx]*pmom::gamma1qg[idx]* (3.*B1q - 3.*beta0 + pmom::gamma1gg[idx] + 5.*pmom::gamma1qq[idx] - 3.*A1q*LF))*LQ3 - 3.*A1q*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LQ4 + 6.*C1qg[idx]*LQF* (2.*pmom::gamma2gq[idx] + pmom::gamma1gq[idx]*(2.*C1qq[idx] + (beta0 + pmom::gamma1gg[idx] + 3.*pmom::gamma1qq[idx])*LF - (2.*B1q - 3.*beta0 + pmom::gamma1gg[idx] + 3.*pmom::gamma1qq[idx])*LQ - A1q*LQ2)))/12.;
	    H3st_qqbp += (2.*LF*(6.*(C2qg[idx]*pmom::gamma1gq[idx] + 2.*C2qqp[idx]*pmom::gamma1qq[idx] + 2.*C1qq[idx]*pmom::gamma2qqp[idx] + pmom::gamma3qqp[idx]) + 6.*(C1qq[idx]*pmom::gamma1gq[idx]*pmom::gamma1qg[idx] + pmom::gamma1qg[idx]*pmom::gamma2gq[idx] + (beta0 + 2.*pmom::gamma1qq[idx])* pmom::gamma2qqp[idx])*LF + pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*(3.*beta0 + pmom::gamma1gg[idx] + 5.*pmom::gamma1qq[idx])* LF2) - 6.*(2.*C2qg[idx]*pmom::gamma1gq[idx] + 4.*C2qqp[idx]*pmom::gamma1qq[idx] + 4.*C1qq[idx]*pmom::gamma2qqp[idx] + 2.*pmom::gamma3qqp[idx] + 2.*(2.*C1qq[idx]*pmom::gamma1gq[idx]*pmom::gamma1qg[idx] + pmom::gamma1qg[idx]*pmom::gamma2gq[idx] + pmom::gamma1gq[idx]*pmom::gamma2qg[idx] + B1q*pmom::gamma2qqp[idx] + 4.*pmom::gamma1qq[idx]*pmom::gamma2qqp[idx])*LF + pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*(B1q + beta0 + pmom::gamma1gg[idx] + 5.*pmom::gamma1qq[idx])*LF2)*LQ + 3.*(4.*C1qq[idx]*pmom::gamma1gq[idx]*pmom::gamma1qg[idx] + 4.*pmom::gamma1qg[idx]*pmom::gamma2gq[idx] + 4.*(B1q - beta0 + 2.*pmom::gamma1qq[idx])*pmom::gamma2qqp[idx] - 2.*A1q*pmom::gamma2qqp[idx]*LF + pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LF*(4.*B1q - 2.*beta0 + 2.*pmom::gamma1gg[idx] + 10.*pmom::gamma1qq[idx] - A1q*LF))* LQ2 - 2.*(-3*A1q*pmom::gamma2qqp[idx] + pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*(3.*B1q - 3.*beta0 + pmom::gamma1gg[idx] + 5.*pmom::gamma1qq[idx] - 3.*A1q*LF))*LQ3 - 3.*A1q*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LQ4 + 6.*C1qg[idx]*LQF*(2.*pmom::gamma2gq[idx] + pmom::gamma1gq[idx]*(2.*C1qq[idx] + (beta0 + pmom::gamma1gg[idx] + 3.*pmom::gamma1qq[idx])*LF - (2.*B1q - 3.*beta0 + pmom::gamma1gg[idx] + 3.*pmom::gamma1qq[idx])*LQ - A1q*LQ2)))/12.;
	    H3st_qg   += (12.*C1qq[idx]*(2.*C1qg[idx]*pmom::gamma1gg[idx] + 2.*C1qq[idx]*pmom::gamma1qg[idx] + 2.*pmom::gamma2qg[idx] + pmom::gamma1qg[idx]*(pmom::gamma1gg[idx] + pmom::gamma1qq[idx])*LQF)*LQF + 12.*C1qg[idx]*(2.*C1qg[idx]*pmom::gamma1gq[idx] + 2.*C1qq[idx]*pmom::gamma1qq[idx] + 2.*pmom::gamma2qq[idx] + (pmom::gamma1gq[idx]*pmom::gamma1qg[idx] + pow(pmom::gamma1qq[idx],2))*LQF)*LQF - 6.*A1q*pmom::gamma1qg[idx]*(pmom::gamma1gg[idx] + pmom::gamma1qq[idx])*LQ4 + 12.*beta0*LQF*((2.*pmom::gamma2qg[idx] + pmom::gamma1qg[idx]*(pmom::gamma1gg[idx] + pmom::gamma1qq[idx])*LQF)* (LF + LQ) + C1qg[idx]*pmom::gamma1gg[idx]*(LF + 3.*LQ) + C1qq[idx]*pmom::gamma1qg[idx]*(LF + 3.*LQ)) + 12.*pmom::gamma1qq[idx]*LQF*(2.*C2qg[idx] + (2.*C1qq[idx]*pmom::gamma1qg[idx] + 2.*pmom::gamma2qg[idx] + pmom::gamma1qg[idx]*(pmom::gamma1gg[idx] + pmom::gamma1qq[idx])*LQF)*LQF + beta0*C1qg[idx]*(LF + 3.*LQ) + C1qg[idx]*(2.*pmom::gamma1gg[idx]*LF - 2.*(B1q + pmom::gamma1gg[idx])*LQ - A1q*LQ2)) + 12.*pmom::gamma1qg[idx]*LQF* (2.*C2qq[idx] + LF*(2.*(C1qg[idx]*pmom::gamma1gq[idx] + C1qq[idx]*pmom::gamma1qq[idx] + pmom::gamma2qq[idx]) + (pmom::gamma1gq[idx]*pmom::gamma1qg[idx] + pow(pmom::gamma1qq[idx],2))*LF) - 2.*(C1qg[idx]*pmom::gamma1gq[idx] + pmom::gamma2qq[idx] + pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LF + (B1q + pmom::gamma1qq[idx])*(C1qq[idx] + pmom::gamma1qq[idx]*LF))*LQ + (pmom::gamma1gq[idx]*pmom::gamma1qg[idx] + 2.*B1q*pmom::gamma1qq[idx] + pow(pmom::gamma1qq[idx],2) - A1q*(C1qq[idx] + pmom::gamma1qq[idx]*LF))*LQ2 + A1q*pmom::gamma1qq[idx]*LQ3 + beta0*(LF*(C1qq[idx] + 2.*pmom::gamma1qq[idx]*LF) + 3.*C1qq[idx]*LQ - 2.*pmom::gamma1qq[idx]*LQ2)) + pmom::gamma1qg[idx]*LQF*(12.*beta1*(LF + LQ) + 8.*pow(beta0,2)*(LF2 + LF*LQ + LQ2) - 2.*beta0*LQ*(6.*B1q*(LF + 2.*LQ) + A1q*LQ*(3.*LF + 7.*LQ)) + 3.*LQ*(-8*B2q + LQ*(-4*A2q + pow(2.*B1q + A1q*LQ,2)))) - 24.*LQ*(C2qg[idx]*pmom::gamma1gg[idx] + C1qg[idx]*pmom::gamma2gg[idx] + C1qq[idx]*pmom::gamma2qg[idx] + pmom::gamma3qg[idx] + pmom::gamma1qg[idx]*(C2qq[idx] + C2qqb[idx] + (C2qqbp[idx] + C2qqp[idx])*(-1 + nf))) + 4.*pmom::gamma1qg[idx]*LF3*(pow(pmom::gamma1gg[idx],2) + pmom::gamma1gg[idx]*pmom::gamma1qq[idx] + pow(pmom::gamma1qq[idx],2) + 2.*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*nf) + 24.*LF*(C2qg[idx]*pmom::gamma1gg[idx] + C2qq[idx]*pmom::gamma1qg[idx] + C2qqb[idx]*pmom::gamma1qg[idx] - (C2qqbp[idx] + C2qqp[idx])*pmom::gamma1qg[idx] + C1qg[idx]*pmom::gamma2gg[idx] + C1qq[idx]*pmom::gamma2qg[idx] + pmom::gamma3qg[idx] - B1q*(C1qg[idx]*pmom::gamma1gg[idx] + C1qq[idx]*pmom::gamma1qg[idx] + pmom::gamma2qg[idx])* LQ + B1q*pmom::gamma1qg[idx]*(pmom::gamma1gg[idx] + pmom::gamma1qq[idx])*LQ2 + (LQ*(-2.*(C1qg[idx]*pow(pmom::gamma1gg[idx],2) + C1qq[idx]*pmom::gamma1qg[idx]*(pmom::gamma1gg[idx] + pmom::gamma1qq[idx]) + (pmom::gamma1gg[idx] + pmom::gamma1qq[idx])*pmom::gamma2qg[idx] + pmom::gamma1qg[idx]*(pmom::gamma2gg[idx] + pmom::gamma2qq[idx] + pmom::gamma2qqb[idx] - pmom::gamma2qqbp[idx] - pmom::gamma2qqp[idx])) + pmom::gamma1qg[idx]*(pow(pmom::gamma1gg[idx],2) + pmom::gamma1gg[idx]*pmom::gamma1qq[idx] + pow(pmom::gamma1qq[idx],2))*LQ - A1q*(C1qg[idx]*pmom::gamma1gg[idx] + C1qq[idx]*pmom::gamma1qg[idx] + pmom::gamma2qg[idx])*LQ + A1q*pmom::gamma1qg[idx]*(pmom::gamma1gg[idx] + pmom::gamma1qq[idx])*LQ2))/2. + pmom::gamma1qg[idx]*(C2qqbp[idx] + C2qqp[idx] - LQ*(2.*C1qg[idx]*pmom::gamma1gq[idx] + pmom::gamma2qqbp[idx] + pmom::gamma2qqp[idx] - pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LQ))*nf) + 12.*LQ2*(C1qq[idx]*pmom::gamma1gg[idx]*pmom::gamma1qg[idx] + C1qq[idx]*pmom::gamma1qg[idx]*pmom::gamma1qq[idx] + 2.*pmom::gamma1qg[idx]*pmom::gamma2gg[idx] + 2.*pmom::gamma1qq[idx]*pmom::gamma2qg[idx] + 2.*B1q*(C1qg[idx]*pmom::gamma1gg[idx] + C1qq[idx]*pmom::gamma1qg[idx] + pmom::gamma2qg[idx]) + C1qg[idx]*(pow(pmom::gamma1gg[idx],2) + 2.*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*nf)) + 4.*LQ3*(3.*A1q*(C1qg[idx]*pmom::gamma1gg[idx] + C1qq[idx]*pmom::gamma1qg[idx] + pmom::gamma2qg[idx]) - pmom::gamma1qg[idx]*(pow(pmom::gamma1gg[idx],2) + pmom::gamma1gg[idx]*pmom::gamma1qq[idx] + pow(pmom::gamma1qq[idx],2) + 3.*B1q*(pmom::gamma1gg[idx] + pmom::gamma1qq[idx]) + 2.*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*nf)) - 6.*LF2*(-2.*C1qq[idx]*pmom::gamma1qg[idx]*(pmom::gamma1gg[idx] + pmom::gamma1qq[idx]) - 4.*pmom::gamma1qg[idx]*pmom::gamma2gg[idx] - 4.*pmom::gamma1qq[idx]*pmom::gamma2qg[idx] - 2.*C1qg[idx]*(pow(pmom::gamma1gg[idx],2) + 2.*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*nf) + pmom::gamma1qg[idx]*LQ*(2.*(pmom::gamma1gg[idx]*(B1q + pmom::gamma1gg[idx]) + (B1q + pmom::gamma1gg[idx])*pmom::gamma1qq[idx] + pow(pmom::gamma1qq[idx],2)) + A1q*(pmom::gamma1gg[idx] + pmom::gamma1qq[idx])*LQ + 4.*pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*nf)))/24.;
	    H3st_qbg  += ((2.*C2qqb[idx]*pmom::gamma1qg[idx] + (C1qg[idx] + pmom::gamma1qg[idx]*LQF)*(2.*C1qg[idx]*pmom::gamma1gq[idx] + 2.*pmom::gamma2qqb[idx] + pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LQF))*LQF)/2.;
	    H3st_qpg  += ((2.*C2qqp[idx]*pmom::gamma1qg[idx] + (C1qg[idx] + pmom::gamma1qg[idx]*LQF)*(2.*C1qg[idx]*pmom::gamma1gq[idx] + 2.*pmom::gamma2qqp[idx] + pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LQF))*LQF)/2.;
	    H3st_qbpg += ((2.*C2qqbp[idx]*pmom::gamma1qg[idx] + (C1qg[idx] + pmom::gamma1qg[idx]*LQF)*(2.*C1qg[idx]*pmom::gamma1gq[idx] + 2.*pmom::gamma2qqbp[idx] + pmom::gamma1gq[idx]*pmom::gamma1qg[idx]*LQF))*LQF)/2.;
	    H3st_gg   += (LQF*(4.*pow(C1qg[idx],2)*pmom::gamma1gg[idx] + 4.*C2qg[idx]*pmom::gamma1qg[idx] + 2.*C1qg[idx]*(2.*C1qq[idx]*pmom::gamma1qg[idx] + 2.*pmom::gamma2qg[idx] + pmom::gamma1qg[idx]*(beta0*LF + 3.*pmom::gamma1gg[idx]*LF + pmom::gamma1qq[idx]*LF - 2.*B1q*LQ + 3.*beta0*LQ - 3.*pmom::gamma1gg[idx]*LQ - pmom::gamma1qq[idx]*LQ - A1q*LQ2)) + pmom::gamma1qg[idx]*LQF*(4.*C1qq[idx]*pmom::gamma1qg[idx] + 4.*pmom::gamma2qg[idx] + pmom::gamma1qg[idx]*(2.*pmom::gamma1gg[idx]*LF + 2.*pmom::gamma1qq[idx]*LF - 2.*B1q*LQ - 2.*pmom::gamma1gg[idx]*LQ - 2.*pmom::gamma1qq[idx]*LQ - A1q*LQ2 + 2.*beta0*(LF + LQ)))))/2.;
	    
	    //cout << endl;
	    //cout << "c3new.txt and hst-channels.txt " << endl;
	    //cout << "H3st_qqb  " << H3st_qqb  << endl;
	    //cout << "H3st_qg   " << H3st_qg   << endl;
	    //cout << "H3st_qq   " << H3st_qq   << endl;
	    //cout << "H3st_qqp  " << H3st_qqp  << endl;
	    //cout << "H3st_qqbp " << H3st_qqbp << endl;
	    //cout << "H3st_gg   " << H3st_gg   << endl;
	    //cout << "H3st_qbg  " << H3st_qbg  << endl;
	    //cout << "H3st_qpg  " << H3st_qpg  << endl;
	    //cout << "H3st_qbpg " << H3st_qbpg << endl;
	  }
      
	Hqqb[i]  += as3*H3st_qqb;
	Hqg[i]   += as3*H3st_qg;
	Hqq[i]   += as3*H3st_qq;
	Hqqp[i]  += as3*H3st_qqp;
	Hqqbp[i] += as3*H3st_qqbp;
	Hgg[i]   += as3*H3st_gg;
	Hqbg[i]  += as3*H3st_qbg;
	Hqpg[i]  += as3*H3st_qpg;
	Hqbpg[i] += as3*H3st_qbpg;

      }
  else
    //rapidity dependent
    for (int sign = mesq::positive; sign <= mesq::negative; sign++)
      for (int i1 = 0; i1 < mellinint::mdim; i1++)
	for (int i2 = 0; i2 < mellinint::mdim; i2++)
	  {
	    int ii1 = anomalous::index(i1,mesq::positive);
	    int ii2 = anomalous::index(i2,sign);
	    int idx = index(i1,i2,sign);
	    
	    complex <double> H1st_qqb = C1qq_1[ii1] + C1qq_2[ii2];
	    complex <double> H1st_qg_1  = C1qg_1[ii1];
	    complex <double> H1st_qg_2  = C1qg_2[ii2];

	    //Resummation scale variations
	    H1st_qqb += (-B1q*LQ -A1q/2.*LQ2);

	    //Factorization scale variations
	    if (!opts.mufevol)
	      {
		H1st_qqb += (pmom::gamma1qq_1[ii1] + pmom::gamma1qq_2[ii2])*LQF ;
		H1st_qg_1  += pmom::gamma1qg_1[ii1]*LQF;
		H1st_qg_2  += pmom::gamma1qg_2[ii2]*LQF;
	      }
	    //	    cout << pmom::gamma1qq_1[ii1] << "  " << pmom::gamma1qq_2[ii2] << endl;
	    
	    Hqqb[idx]   += as*H1st_qqb;
	    Hqg_1[idx]  += as*H1st_qg_1;
	    Hqg_2[idx]  += as*H1st_qg_2;
	    
	    if (opts.order == 1)
	      continue;
	    
	    //NNLO
	    complex <double> H2st_qg_1  = C2qg_1[ii1] + C1qq_2[ii2]*C1qg_1[ii1];
	    complex <double> H2st_qg_2  = C2qg_2[ii2] + C1qq_1[ii1]*C1qg_2[ii2];
	    complex <double> H2st_qqb = C2qq_1[ii1] + C2qq_2[ii2] + C1qq_1[ii1]*C1qq_2[ii2];
	    complex <double> H2st_qq_1  = C2qqb_1[ii1];
	    complex <double> H2st_qq_2  = C2qqb_2[ii2];
	    complex <double> H2st_qqbp_1= C2qqp_1[ii1];
	    complex <double> H2st_qqbp_2= C2qqp_2[ii2];
	    complex <double> H2st_qqp_1 = C2qqbp_1[ii1];// (=H2st_qqbp_1)
	    complex <double> H2st_qqp_2 = C2qqbp_2[ii2];// (=H2st_qqbp_2)
	    complex <double> H2st_gg  = C1qg_1[ii1]*C1qg_2[ii2];

	    //Resummation scale variations
	    H2st_qqb += (- B2q*LQ + (-A2q/2. + pow(B1q,2)/2. - B1q*beta0/2.)*LQ2 + A1q*(B1q/2. - beta0/3.)*LQ3 + pow(A1q,2)*LQ4/8.)
	      + ((beta0 - B1q)*LQ - A1q/2.*LQ2)*(C1qq_1[ii1] + C1qq_2[ii2]);
	    H2st_qg_1 +=  ((beta0 - B1q)*LQ -A1q/2.*LQ2)*C1qg_1[ii1];
	    H2st_qg_2 +=  ((beta0 - B1q)*LQ -A1q/2.*LQ2)*C1qg_2[ii2];

	    //Renormalization scale variations
	    H2st_qqb  += - beta0*H1st_qqb*LR;
	    H2st_qg_1 += - beta0*H1st_qg_1*LR;
	    H2st_qg_2 += - beta0*H1st_qg_2*LR;

	    //Factorization scale variations
	    if (!opts.mufevol)
	      {
		H2st_qqb += LQF*(C1qg_1[ii1]*pmom::gamma1gq_1[ii1] + C1qq_2[ii2]*pmom::gamma1qq_1[ii1] + C1qg_2[ii2]*pmom::gamma1gq_2[ii2] + C1qq_2[ii2]*pmom::gamma1qq_2[ii2] + C1qq_1[ii1]*(pmom::gamma1qq_1[ii1] + pmom::gamma1qq_2[ii2]) + pmom::gamma2qq_1[ii1] + pmom::gamma2qq_2[ii2] + ((pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1] + pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2] + (pmom::gamma1qq_1[ii1] + pmom::gamma1qq_2[ii2])*(beta0 + pmom::gamma1qq_1[ii1] + pmom::gamma1qq_2[ii2]))*LQF)/ 2. - ((pmom::gamma1qq_1[ii1] + pmom::gamma1qq_2[ii2])*LQ*(2.*B1q - 2.*beta0 + A1q*LQ))/2.);
		H2st_qg_1 += (LQF*(2.*C1qg_1[ii1]*(pmom::gamma1qq_2[ii2] + pmom::gamma1gg_1[ii1]) + 2.*pmom::gamma2qg_1[ii1] + pmom::gamma1qg_1[ii1]*(2.*C1qq_2[ii2] + 2.*C1qq_1[ii1] + (beta0 + 2.*pmom::gamma1qq_2[ii2] + pmom::gamma1gg_1[ii1] + pmom::gamma1qq_1[ii1])*LF - (2.*B1q - beta0 + 2.*pmom::gamma1qq_2[ii2] + pmom::gamma1gg_1[ii1] + pmom::gamma1qq_1[ii1])*LQ - A1q*LQ2)))/2.;
		H2st_qg_2 += (LQF*(2.*C1qg_2[ii2]*(pmom::gamma1qq_1[ii1] + pmom::gamma1gg_2[ii2]) + 2.*pmom::gamma2qg_2[ii2] + pmom::gamma1qg_2[ii2]*(2.*C1qq_1[ii1] + 2.*C1qq_2[ii2] + (beta0 + 2.*pmom::gamma1qq_1[ii1] + pmom::gamma1gg_2[ii2] + pmom::gamma1qq_2[ii2])*LF - (2.*B1q - beta0 + 2.*pmom::gamma1qq_1[ii1] + pmom::gamma1gg_2[ii2] + pmom::gamma1qq_2[ii2])*LQ - A1q*LQ2)))/2.;
		H2st_qq_1 += ((2.*C1qg_1[ii1]*pmom::gamma1gq_1[ii1] + 2.*pmom::gamma2qqb_1[ii1] + pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LQF)* LQF)/2.;
		H2st_qq_2 += ((2.*C1qg_2[ii2]*pmom::gamma1gq_2[ii2] + 2.*pmom::gamma2qqb_2[ii2] + pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LQF)* LQF)/2.;
		H2st_qqp_1 += ((2.*C1qg_1[ii1]*pmom::gamma1gq_1[ii1] + 2.*pmom::gamma2qqbp_1[ii1] + pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LQF)* LQF)/2.;
		H2st_qqp_2 += ((2.*C1qg_2[ii2]*pmom::gamma1gq_2[ii2] + 2.*pmom::gamma2qqbp_2[ii2] + pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LQF)* LQF)/2.;
		H2st_qqbp_1 += ((2.*C1qg_1[ii1]*pmom::gamma1gq_1[ii1] + 2.*pmom::gamma2qqp_1[ii1] + pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LQF)* LQF)/2.;
		H2st_qqbp_2 += ((2.*C1qg_2[ii2]*pmom::gamma1gq_2[ii2] + 2.*pmom::gamma2qqp_2[ii2] + pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LQF)* LQF)/2.;
		H2st_gg += (C1qg_2[ii2]*pmom::gamma1qg_1[ii1] + pmom::gamma1qg_2[ii2]*(C1qg_1[ii1] + pmom::gamma1qg_1[ii1]*LQF))*LQF;
	      }
	    
	    Hqqb[idx]    += as2*H2st_qqb;
	    Hqg_1[idx]   += as2*H2st_qg_1;
	    Hqg_2[idx]   += as2*H2st_qg_2;
	    Hqq_1[idx]   += as2*H2st_qq_1;
	    Hqq_2[idx]   += as2*H2st_qq_2;
	    Hqqp_1[idx]  += as2*H2st_qqp_1;
	    Hqqp_2[idx]  += as2*H2st_qqp_2;
	    Hqqbp_1[idx] += as2*H2st_qqbp_1;
	    Hqqbp_2[idx] += as2*H2st_qqbp_2;
	    Hgg[idx]     += as2*H2st_gg;
	    
	    if (opts.order == 2)
	      continue;

	    //NNNLO
	    complex <double> H3st_qg_1  = C3qg_1[ii1] + C1qq_2[ii2]*C2qg_1[ii1] + C2qq_2[ii2]*C1qg_1[ii1];
	    complex <double> H3st_qg_2  = C3qg_2[ii2] + C1qq_1[ii1]*C2qg_2[ii2] + C2qq_1[ii1]*C1qg_2[ii2];
	    complex <double> H3st_qbg_1 = C2qqb_2[ii2]*C1qg_1[ii1];
	    complex <double> H3st_qbg_2 = C2qqb_1[ii1]*C1qg_2[ii2];
	    complex <double> H3st_qpg_1 = C2qqp_2[ii2]*C1qg_1[ii1];
	    complex <double> H3st_qpg_2 = C2qqp_1[ii1]*C1qg_2[ii2];
	    complex <double> H3st_qbpg_1= C2qqbp_2[ii2]*C1qg_1[ii1];
	    complex <double> H3st_qbpg_2= C2qqbp_1[ii1]*C1qg_2[ii2];
	    complex <double> H3st_qqb = C3qq_1[ii1] + C3qq_2[ii2] + C1qq_1[ii1]*C2qq_2[ii2] + C1qq_2[ii2]*C2qq_1[ii1];
	    complex <double> H3st_qq_1  = C3qqb_1[ii1] + C1qq_2[ii2]*C2qqb_1[ii1];
	    complex <double> H3st_qq_2  = C3qqb_2[ii2] + C1qq_1[ii1]*C2qqb_2[ii2];
	    complex <double> H3st_qqbp_1= C3qqp_1[ii1] + C1qq_2[ii2]*C2qqp_1[ii1];
	    complex <double> H3st_qqbp_2= C3qqp_2[ii2] + C1qq_1[ii1]*C2qqp_2[ii2];
	    complex <double> H3st_qqp_1 = C3qqbp_1[ii1]+ C1qq_2[ii2]*C2qqbp_1[ii1]; // (=C3qqbp_1+ C1qq_2[ii2]*C2qqp_1[ii1])
	    complex <double> H3st_qqp_2 = C3qqbp_2[ii2]+ C1qq_1[ii1]*C2qqbp_2[ii2]; // (=C3qqbp_2+ C1qq_1[ii1]*C2qqp_2[ii2])
	    complex <double> H3st_gg  = C1qg_1[ii1]*C2qg_2[ii2] + C1qg_2[ii2]*C2qg_1[ii1];

	    //Resummation scale variations
	    H3st_qqb += (-B3q*LQ + (-A3q/2. + B1q*B2q - B2q*beta0 - B1q*beta1/2.)*LQ2
			 + (A2q*B1q/2. - pow(B1q,3)/6. + A1q*B2q/2. - 2*A2q*beta0/3. + pow(B1q,2)*beta0/2. - B1q*pow(beta0,2)/3. - A1q*beta1/3.)*LQ3
			 + (A1q*A2q/4. - A1q*pow(B1q,2)/4. + 7.*A1q*B1q*beta0/12. - A1q*pow(beta0,2)/4.)*LQ4 + (-pow(A1q,2)*B1q/8. + pow(A1q,2)*beta0/6.)*LQ5 - pow(A1q,3)*LQ6/48.)
	      + ((beta1 - B2q)*LQ + (-A2q/2. + pow(B1q,2)/2. - 3./2.*B1q*beta0 + pow(beta0,2))*LQ2 + A1q*(B1q/2. - 5*beta0/6.)*LQ3 + pow(A1q,2)*LQ4/8.)*(C1qq_1[ii1] + C1qq_2[ii2])
	      + ((-B1q + 2*beta0)*LQ - A1q/2.*LQ2)*(C2qq_1[ii1] + C2qq_2[ii2] + C1qq_1[ii1]*C1qq_2[ii2]);

	    H3st_qq_1   +=  ((-B1q + 2*beta0)*LQ - A1q/2.*LQ2)*C2qqb_1[ii1];
	    H3st_qq_2   +=  ((-B1q + 2*beta0)*LQ - A1q/2.*LQ2)*C2qqb_2[ii2];

	    H3st_qqp_1  +=  ((-B1q + 2*beta0)*LQ - A1q/2.*LQ2)*C2qqbp_1[ii1];
	    H3st_qqp_2  +=  ((-B1q + 2*beta0)*LQ - A1q/2.*LQ2)*C2qqbp_2[ii2];

	    H3st_qqbp_1 +=  ((-B1q + 2*beta0)*LQ - A1q/2.*LQ2)*C2qqp_1[ii1];
	    H3st_qqbp_2 +=  ((-B1q + 2*beta0)*LQ - A1q/2.*LQ2)*C2qqp_2[ii2];

	    H3st_qg_1   += ((beta1 - B2q)*LQ + (-A2q/2. + pow(B1q,2)/2. - 3./2.*B1q*beta0 + pow(beta0,2))*LQ2 + A1q*(B1q/2. - 5*beta0/6.)*LQ3 + pow(A1q,2)*LQ4/8.)*C1qg_1[ii1]
	      + ((-B1q + 2*beta0)*LQ - A1q/2.*LQ2)*(C2qg_1[ii1] + C1qq_2[ii2]*C1qg_1[ii1]);

	    H3st_qg_2   += ((beta1 - B2q)*LQ + (-A2q/2. + pow(B1q,2)/2. - 3./2.*B1q*beta0 + pow(beta0,2))*LQ2 + A1q*(B1q/2. - 5*beta0/6.)*LQ3 + pow(A1q,2)*LQ4/8.)*C1qg_2[ii2]
	      + ((-B1q + 2*beta0)*LQ - A1q/2.*LQ2)*(C2qg_2[ii2] + C1qq_1[ii1]*C1qg_2[ii2]);

	    H3st_qbg_1 +=  ((-B1q + 2*beta0)*LQ - A1q/2.*LQ2)*C1qqb_2[ii2] * C1qg_1[ii1];
	    H3st_qbg_2 +=  ((-B1q + 2*beta0)*LQ - A1q/2.*LQ2)*C1qqb_1[ii1] * C1qg_2[ii2];

	    H3st_qpg_1 +=  ((-B1q + 2*beta0)*LQ - A1q/2.*LQ2)*C1qqp_2[ii2] * C1qg_1[ii1];
	    H3st_qpg_2 +=  ((-B1q + 2*beta0)*LQ - A1q/2.*LQ2)*C1qqp_1[ii1] * C1qg_2[ii2];

	    H3st_qbpg_1+=  ((-B1q + 2*beta0)*LQ - A1q/2.*LQ2)*C1qqbp_2[ii2] * C1qg_1[ii1];
	    H3st_qbpg_1+=  ((-B1q + 2*beta0)*LQ - A1q/2.*LQ2)*C1qqbp_1[ii1] * C1qg_2[ii2];

	    H3st_gg   += ((-B1q + 2*beta0)*LQ - A1q/2.*LQ2)*(C1qg_1[ii1]*C1qg_2[ii2]);

	    //Renormalization scale variations
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

	    //Factorization scale variations
	    if (!opts.mufevol)
	      {
		double nf = 5.;
		H3st_qqb += C1qg_2[ii2]*C1qq_1[ii1]*pmom::gamma1gq_2[ii2]*LF + C2qg_2[ii2]*pmom::gamma1gq_2[ii2]*LF + C1qq_2[ii2]*C1qq_1[ii1]*pmom::gamma1qq_2[ii2]*LF + C2qq_2[ii2]*pmom::gamma1qq_2[ii2]*LF + C2qq_1[ii1]*pmom::gamma1qq_2[ii2]*LF + C1qq_2[ii2]*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LF + C2qg_1[ii1]*pmom::gamma1gq_1[ii1]*LF + C1qq_2[ii2]*C1qq_1[ii1]*pmom::gamma1qq_1[ii1]*LF + C2qq_2[ii2]*pmom::gamma1qq_1[ii1]*LF + C2qq_1[ii1]*pmom::gamma1qq_1[ii1]*LF + C1qg_1[ii1]*pmom::gamma2gq_1[ii1]*LF + C1qg_2[ii2]*pmom::gamma2gq_2[ii2]*LF + C1qq_2[ii2]*pmom::gamma2qq_2[ii2]*LF + C1qq_1[ii1]*pmom::gamma2qq_2[ii2]*LF + C1qq_2[ii2]*pmom::gamma2qq_1[ii1]*LF + C1qq_1[ii1]*pmom::gamma2qq_1[ii1]*LF + pmom::gamma3qq_2[ii2]*LF + pmom::gamma3qq_1[ii1]*LF + (beta0*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LF2)/2. + (C1qg_2[ii2]*pmom::gamma1gg_2[ii2]*pmom::gamma1gq_2[ii2]*LF2)/2. + (C1qq_2[ii2]*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF2)/ 2. + (C1qq_1[ii1]*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF2)/2. + (beta1*pmom::gamma1qq_2[ii2]*LF2)/2. + (beta0*C1qq_2[ii2]*pmom::gamma1qq_2[ii2]*LF2)/2. + (beta0*C1qq_1[ii1]*pmom::gamma1qq_2[ii2]*LF2)/2. + (C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*pmom::gamma1qq_2[ii2]*LF2)/2. + (C1qq_2[ii2]*pow(pmom::gamma1qq_2[ii2],2)*LF2)/2. + (C1qq_1[ii1]*pow(pmom::gamma1qq_2[ii2],2)*LF2)/2. + (beta0*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LF2)/2. + C1qg_1[ii1]*pmom::gamma1qq_2[ii2]*pmom::gamma1gq_1[ii1]*LF2 + (C1qq_2[ii2]*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF2)/2. + (beta1*pmom::gamma1qq_1[ii1]*LF2)/2. + (beta0*C1qq_2[ii2]*pmom::gamma1qq_1[ii1]*LF2)/2. + (beta0*C1qq_1[ii1]*pmom::gamma1qq_1[ii1]*LF2)/2. + C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*pmom::gamma1qq_1[ii1]*LF2 + C1qq_2[ii2]*pmom::gamma1qq_2[ii2]*pmom::gamma1qq_1[ii1]*LF2 + C1qq_1[ii1]*pmom::gamma1qq_2[ii2]*pmom::gamma1qq_1[ii1]*LF2 + (C1qq_2[ii2]*pow(pmom::gamma1qq_1[ii1],2)*LF2)/2. + pmom::gamma1qg_2[ii2]*pmom::gamma2gq_2[ii2]*LF2 + beta0*pmom::gamma2qq_2[ii2]*LF2 + pmom::gamma1qq_2[ii2]*pmom::gamma2qq_2[ii2]*LF2 + pmom::gamma1qq_1[ii1]*pmom::gamma2qq_2[ii2]*LF2 + pmom::gamma1qg_1[ii1]*pmom::gamma2gq_1[ii1]*LF2 + beta0*pmom::gamma2qq_1[ii1]*LF2 + pmom::gamma1qq_2[ii2]*pmom::gamma2qq_1[ii1]*LF2 + pmom::gamma1qq_1[ii1]*pmom::gamma2qq_1[ii1]*LF2 + (beta0*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF3)/2. + (pow(beta0,2)*pmom::gamma1qq_2[ii2]*LF3)/3. + (beta0*pow(pmom::gamma1qq_2[ii2],2)*LF3)/2. + (beta0*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF3)/2. + (pmom::gamma1qq_2[ii2]*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF3)/2. + (pow(beta0,2)*pmom::gamma1qq_1[ii1]*LF3)/3. + (pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*pmom::gamma1qq_1[ii1]*LF3)/2. + beta0*pmom::gamma1qq_2[ii2]*pmom::gamma1qq_1[ii1]* LF3 + (pow(pmom::gamma1qq_2[ii2],2)*pmom::gamma1qq_1[ii1]*LF3)/2. + (beta0*pow(pmom::gamma1qq_1[ii1],2)*LF3)/2. + (pmom::gamma1qq_2[ii2]*pow(pmom::gamma1qq_1[ii1],2)*LF3)/2. + ((pmom::gamma1gg_1[ii1]*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1] + 2.*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*pmom::gamma1qq_1[ii1] + pow(pmom::gamma1qq_1[ii1],3))*LF3)/6. + ((pmom::gamma1gg_2[ii2]*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2] + 2.*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*pmom::gamma1qq_2[ii2] + pow(pmom::gamma1qq_2[ii2],3))*LQF3)/6. - C1qg_2[ii2]*C1qq_1[ii1]*pmom::gamma1gq_2[ii2]*LQ - C2qg_2[ii2]*pmom::gamma1gq_2[ii2]*LQ - C1qq_2[ii2]*C1qq_1[ii1]*pmom::gamma1qq_2[ii2]*LQ - C2qq_2[ii2]*pmom::gamma1qq_2[ii2]*LQ - C2qq_1[ii1]*pmom::gamma1qq_2[ii2]*LQ - C1qq_2[ii2]*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LQ - C2qg_1[ii1]*pmom::gamma1gq_1[ii1]*LQ - C1qq_2[ii2]*C1qq_1[ii1]*pmom::gamma1qq_1[ii1]*LQ - C2qq_2[ii2]*pmom::gamma1qq_1[ii1]*LQ - C2qq_1[ii1]*pmom::gamma1qq_1[ii1]*LQ - C1qg_1[ii1]*pmom::gamma2gq_1[ii1]*LQ - C1qg_2[ii2]*pmom::gamma2gq_2[ii2]*LQ - C1qq_2[ii2]*pmom::gamma2qq_2[ii2]*LQ - C1qq_1[ii1]*pmom::gamma2qq_2[ii2]*LQ - C1qq_2[ii2]*pmom::gamma2qq_1[ii1]*LQ - C1qq_1[ii1]*pmom::gamma2qq_1[ii1]*LQ - pmom::gamma3qq_2[ii2]*LQ - pmom::gamma3qq_1[ii1]*LQ - B1q*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LF*LQ + beta0*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LF*LQ - C1qg_2[ii2]*pmom::gamma1gg_2[ii2]*pmom::gamma1gq_2[ii2]*LF*LQ - C1qq_2[ii2]*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF*LQ - C1qq_1[ii1]*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF*LQ - B2q*pmom::gamma1qq_2[ii2]*LF*LQ - B1q*C1qq_2[ii2]*pmom::gamma1qq_2[ii2]*LF*LQ + beta0*C1qq_2[ii2]*pmom::gamma1qq_2[ii2]*LF*LQ - B1q*C1qq_1[ii1]*pmom::gamma1qq_2[ii2]*LF*LQ + beta0*C1qq_1[ii1]*pmom::gamma1qq_2[ii2]*LF*LQ - C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*pmom::gamma1qq_2[ii2]*LF*LQ - C1qq_2[ii2]*pow(pmom::gamma1qq_2[ii2],2)*LF*LQ - C1qq_1[ii1]*pow(pmom::gamma1qq_2[ii2],2)*LF*LQ - B1q*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LF*LQ + beta0*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LF*LQ - 2.*C1qg_1[ii1]*pmom::gamma1qq_2[ii2]*pmom::gamma1gq_1[ii1]*LF*LQ - C1qg_1[ii1]*pmom::gamma1gg_1[ii1]*pmom::gamma1gq_1[ii1]*LF*LQ - C1qq_2[ii2]*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF*LQ - C1qq_1[ii1]*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF*LQ - B2q*pmom::gamma1qq_1[ii1]*LF*LQ - B1q*C1qq_2[ii2]*pmom::gamma1qq_1[ii1]*LF*LQ + beta0*C1qq_2[ii2]*pmom::gamma1qq_1[ii1]*LF*LQ - B1q*C1qq_1[ii1]*pmom::gamma1qq_1[ii1]*LF*LQ + beta0*C1qq_1[ii1]*pmom::gamma1qq_1[ii1]*LF*LQ - 2.*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*pmom::gamma1qq_1[ii1]*LF*LQ - 2.*C1qq_2[ii2]*pmom::gamma1qq_2[ii2]*pmom::gamma1qq_1[ii1]*LF*LQ - 2.*C1qq_1[ii1]*pmom::gamma1qq_2[ii2]*pmom::gamma1qq_1[ii1]*LF* LQ - C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*pmom::gamma1qq_1[ii1]*LF*LQ - C1qq_2[ii2]*pow(pmom::gamma1qq_1[ii1],2)*LF*LQ - C1qq_1[ii1]*pow(pmom::gamma1qq_1[ii1],2)*LF*LQ - pmom::gamma1qg_2[ii2]*pmom::gamma2gq_2[ii2]*LF*LQ - pmom::gamma1gq_2[ii2]*pmom::gamma2qg_2[ii2]*LF*LQ - B1q*pmom::gamma2qq_2[ii2]*LF*LQ - 2.*pmom::gamma1qq_2[ii2]*pmom::gamma2qq_2[ii2]*LF*LQ - 2.*pmom::gamma1qq_1[ii1]*pmom::gamma2qq_2[ii2]*LF*LQ - pmom::gamma1qg_1[ii1]*pmom::gamma2gq_1[ii1]*LF*LQ - pmom::gamma1gq_1[ii1]*pmom::gamma2qg_1[ii1]*LF*LQ - B1q*pmom::gamma2qq_1[ii1]*LF*LQ - 2.*pmom::gamma1qq_2[ii2]*pmom::gamma2qq_1[ii1]*LF*LQ - 2.*pmom::gamma1qq_1[ii1]*pmom::gamma2qq_1[ii1]*LF*LQ - (B1q*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF2*LQ)/2. - (beta0*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF2*LQ)/2. - (B1q*beta0*pmom::gamma1qq_2[ii2]*LF2*LQ)/2. - (B1q*pow(pmom::gamma1qq_2[ii2],2)*LF2*LQ)/2. - (beta0*pow(pmom::gamma1qq_2[ii2],2)*LF2*LQ)/2. - (B1q*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF2*LQ)/2. - (beta0*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF2*LQ)/ 2. - (3.*pmom::gamma1qq_2[ii2]*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF2*LQ)/2. - (B1q*beta0*pmom::gamma1qq_1[ii1]*LF2*LQ)/2. - (3.*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*pmom::gamma1qq_1[ii1]*LF2*LQ)/2. - B1q*pmom::gamma1qq_2[ii2]*pmom::gamma1qq_1[ii1]*LF2*LQ - beta0*pmom::gamma1qq_2[ii2]*pmom::gamma1qq_1[ii1]*LF2*LQ - (3.*pow(pmom::gamma1qq_2[ii2],2)*pmom::gamma1qq_1[ii1]*LF2*LQ)/2. - (B1q*pow(pmom::gamma1qq_1[ii1],2)*LF2*LQ)/2. - (beta0*pow(pmom::gamma1qq_1[ii1],2)*LF2*LQ)/2. - (3.*pmom::gamma1qq_2[ii2]*pow(pmom::gamma1qq_1[ii1],2)*LF2*LQ)/2. - ((pmom::gamma1gg_1[ii1]*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1] + 2.*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*pmom::gamma1qq_1[ii1] + pow(pmom::gamma1qq_1[ii1],3))*LF2*LQ)/2. + B1q*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LQ2 - (3.*beta0*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LQ2)/2. + (C1qg_2[ii2]*pmom::gamma1gg_2[ii2]*pmom::gamma1gq_2[ii2]*LQ2)/2. + (C1qq_2[ii2]*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LQ2)/2. + (C1qq_1[ii1]*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LQ2)/ 2. + B2q*pmom::gamma1qq_2[ii2]*LQ2 - (beta1*pmom::gamma1qq_2[ii2]*LQ2)/2. + B1q*C1qq_2[ii2]*pmom::gamma1qq_2[ii2]*LQ2 - (3.*beta0*C1qq_2[ii2]*pmom::gamma1qq_2[ii2]*LQ2)/2. + B1q*C1qq_1[ii1]*pmom::gamma1qq_2[ii2]*LQ2 - (3.*beta0*C1qq_1[ii1]*pmom::gamma1qq_2[ii2]*LQ2)/2. + (C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*pmom::gamma1qq_2[ii2]*LQ2)/2. + (C1qq_2[ii2]*pow(pmom::gamma1qq_2[ii2],2)*LQ2)/2. + (C1qq_1[ii1]*pow(pmom::gamma1qq_2[ii2],2)*LQ2)/2. + B1q*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LQ2 - (3.*beta0*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LQ2)/2. + C1qg_1[ii1]*pmom::gamma1qq_2[ii2]*pmom::gamma1gq_1[ii1]*LQ2 + (C1qq_2[ii2]*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LQ2)/2. + B2q*pmom::gamma1qq_1[ii1]*LQ2 - (beta1*pmom::gamma1qq_1[ii1]*LQ2)/2. + B1q*C1qq_2[ii2]*pmom::gamma1qq_1[ii1]*LQ2 - (3.*beta0*C1qq_2[ii2]*pmom::gamma1qq_1[ii1]*LQ2)/2. + B1q*C1qq_1[ii1]*pmom::gamma1qq_1[ii1]*LQ2 - (3.*beta0*C1qq_1[ii1]*pmom::gamma1qq_1[ii1]*LQ2)/2. + C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*pmom::gamma1qq_1[ii1]*LQ2 + C1qq_2[ii2]*pmom::gamma1qq_2[ii2]*pmom::gamma1qq_1[ii1]*LQ2 + C1qq_1[ii1]*pmom::gamma1qq_2[ii2]*pmom::gamma1qq_1[ii1]*LQ2 + (C1qq_2[ii2]*pow(pmom::gamma1qq_1[ii1],2)*LQ2)/2. + pmom::gamma1qg_2[ii2]*pmom::gamma2gq_2[ii2]*LQ2 + B1q*pmom::gamma2qq_2[ii2]*LQ2 - beta0*pmom::gamma2qq_2[ii2]*LQ2 + pmom::gamma1qq_2[ii2]*pmom::gamma2qq_2[ii2]*LQ2 + pmom::gamma1qq_1[ii1]*pmom::gamma2qq_2[ii2]*LQ2 + pmom::gamma1qg_1[ii1]*pmom::gamma2gq_1[ii1]*LQ2 + B1q*pmom::gamma2qq_1[ii1]*LQ2 - beta0*pmom::gamma2qq_1[ii1]*LQ2 + pmom::gamma1qq_2[ii2]*pmom::gamma2qq_1[ii1]*LQ2 + pmom::gamma1qq_1[ii1]*pmom::gamma2qq_1[ii1]*LQ2 - (A1q*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LF*LQ2)/2. + B1q*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF*LQ2 - (beta0*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF*LQ2)/2. - (A2q*pmom::gamma1qq_2[ii2]*LF*LQ2)/2. + (pow(B1q,2)*pmom::gamma1qq_2[ii2]*LF*LQ2)/2. - (B1q*beta0*pmom::gamma1qq_2[ii2]*LF*LQ2)/2. - (A1q*C1qq_2[ii2]*pmom::gamma1qq_2[ii2]*LF*LQ2)/2. - (A1q*C1qq_1[ii1]*pmom::gamma1qq_2[ii2]*LF*LQ2)/2. + B1q*pow(pmom::gamma1qq_2[ii2],2)*LF*LQ2 - (beta0*pow(pmom::gamma1qq_2[ii2],2)*LF*LQ2)/2. - (A1q*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LF*LQ2)/2. + B1q*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF*LQ2 - (beta0*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF*LQ2)/2. + (3.*pmom::gamma1qq_2[ii2]*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF*LQ2)/2. - (A2q*pmom::gamma1qq_1[ii1]*LF*LQ2)/2. + (pow(B1q,2)*pmom::gamma1qq_1[ii1]*LF*LQ2)/2. - (B1q*beta0*pmom::gamma1qq_1[ii1]*LF*LQ2)/2. - (A1q*C1qq_2[ii2]*pmom::gamma1qq_1[ii1]*LF*LQ2)/2. - (A1q*C1qq_1[ii1]*pmom::gamma1qq_1[ii1]*LF*LQ2)/2. + (3.*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*pmom::gamma1qq_1[ii1]*LF*LQ2)/2. + 2.*B1q*pmom::gamma1qq_2[ii2]*pmom::gamma1qq_1[ii1]*LF*LQ2 - beta0*pmom::gamma1qq_2[ii2]*pmom::gamma1qq_1[ii1]*LF*LQ2 + (3.*pow(pmom::gamma1qq_2[ii2],2)*pmom::gamma1qq_1[ii1]*LF*LQ2)/2. + B1q*pow(pmom::gamma1qq_1[ii1],2)*LF*LQ2 - (beta0*pow(pmom::gamma1qq_1[ii1],2)*LF*LQ2)/2. + (3.*pmom::gamma1qq_2[ii2]*pow(pmom::gamma1qq_1[ii1],2)*LF*LQ2)/2. + ((pmom::gamma1gg_1[ii1]*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1] + 2.*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*pmom::gamma1qq_1[ii1] + pow(pmom::gamma1qq_1[ii1],3))*LF*LQ2)/2. - (A1q*pmom::gamma2qq_2[ii2]*LF*LQ2)/2. - (A1q*pmom::gamma2qq_1[ii1]*LF*LQ2)/2. - (A1q*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF2*LQ2)/4. - (A1q*beta0*pmom::gamma1qq_2[ii2]*LF2*LQ2)/4. - (A1q*pow(pmom::gamma1qq_2[ii2],2)*LF2*LQ2)/4. - (A1q*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF2*LQ2)/4. - (A1q*beta0*pmom::gamma1qq_1[ii1]*LF2*LQ2)/4. - (A1q*pmom::gamma1qq_2[ii2]*pmom::gamma1qq_1[ii1]*LF2*LQ2)/2. - (A1q*pow(pmom::gamma1qq_1[ii1],2)*LF2*LQ2)/4. + (A1q*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LQ3)/2. - (B1q*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LQ3)/2. + (beta0*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LQ3)/2. + (A2q*pmom::gamma1qq_2[ii2]*LQ3)/2. - (pow(B1q,2)*pmom::gamma1qq_2[ii2]*LQ3)/2. + B1q*beta0*pmom::gamma1qq_2[ii2]*LQ3 - (pow(beta0,2)*pmom::gamma1qq_2[ii2]*LQ3)/3. + (A1q*C1qq_2[ii2]*pmom::gamma1qq_2[ii2]*LQ3)/2. + (A1q*C1qq_1[ii1]*pmom::gamma1qq_2[ii2]*LQ3)/2. - (B1q*pow(pmom::gamma1qq_2[ii2],2)*LQ3)/2. + (beta0*pow(pmom::gamma1qq_2[ii2],2)*LQ3)/2. + (A1q*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LQ3)/2. - (B1q*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LQ3)/2. + (beta0*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LQ3)/2. - (pmom::gamma1qq_2[ii2]*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LQ3)/2. + (A2q*pmom::gamma1qq_1[ii1]*LQ3)/2. - (pow(B1q,2)*pmom::gamma1qq_1[ii1]*LQ3)/2. + B1q*beta0*pmom::gamma1qq_1[ii1]*LQ3 - (pow(beta0,2)*pmom::gamma1qq_1[ii1]*LQ3)/3. + (A1q*C1qq_2[ii2]*pmom::gamma1qq_1[ii1]*LQ3)/2. + (A1q*C1qq_1[ii1]*pmom::gamma1qq_1[ii1]*LQ3)/2. - (pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*pmom::gamma1qq_1[ii1]*LQ3)/2. - B1q*pmom::gamma1qq_2[ii2]*pmom::gamma1qq_1[ii1]*LQ3 + beta0*pmom::gamma1qq_2[ii2]*pmom::gamma1qq_1[ii1]*LQ3 - (pow(pmom::gamma1qq_2[ii2],2)*pmom::gamma1qq_1[ii1]*LQ3)/2. - (B1q*pow(pmom::gamma1qq_1[ii1],2)*LQ3)/2. + (beta0*pow(pmom::gamma1qq_1[ii1],2)*LQ3)/2. - (pmom::gamma1qq_2[ii2]*pow(pmom::gamma1qq_1[ii1],2)*LQ3)/2. - ((pmom::gamma1gg_1[ii1]*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1] + 2.*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*pmom::gamma1qq_1[ii1] + pow(pmom::gamma1qq_1[ii1],3))*LQ3)/6. + (A1q*pmom::gamma2qq_2[ii2]*LQ3)/2. + (A1q*pmom::gamma2qq_1[ii1]*LQ3)/2. + (A1q*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF*LQ3)/2. + (A1q*B1q*pmom::gamma1qq_2[ii2]*LF*LQ3)/2. - (A1q*beta0*pmom::gamma1qq_2[ii2]*LF*LQ3)/3. + (A1q*pow(pmom::gamma1qq_2[ii2],2)*LF*LQ3)/2. + (A1q*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF*LQ3)/2. + (A1q*B1q*pmom::gamma1qq_1[ii1]*LF*LQ3)/2. - (A1q*beta0*pmom::gamma1qq_1[ii1]*LF*LQ3)/3. + A1q*pmom::gamma1qq_2[ii2]*pmom::gamma1qq_1[ii1]*LF*LQ3 + (A1q*pow(pmom::gamma1qq_1[ii1],2)*LF*LQ3)/2. - (A1q*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LQ4)/4. - (A1q*B1q*pmom::gamma1qq_2[ii2]*LQ4)/2. + (7.*A1q*beta0*pmom::gamma1qq_2[ii2]*LQ4)/12. - (A1q*pow(pmom::gamma1qq_2[ii2],2)*LQ4)/4. - (A1q*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LQ4)/4. - (A1q*B1q*pmom::gamma1qq_1[ii1]*LQ4)/2. + (7.*A1q*beta0*pmom::gamma1qq_1[ii1]*LQ4)/12. - (A1q*pmom::gamma1qq_2[ii2]*pmom::gamma1qq_1[ii1]*LQ4)/2. - (A1q*pow(pmom::gamma1qq_1[ii1],2)*LQ4)/4. + (pow(A1q,2)*pmom::gamma1qq_2[ii2]*LF*LQ4)/8. + (pow(A1q,2)*pmom::gamma1qq_1[ii1]*LF*LQ4)/8. - (pow(A1q,2)*pmom::gamma1qq_2[ii2]*LQ5)/8. - (pow(A1q,2)*pmom::gamma1qq_1[ii1]*LQ5)/8. + ((C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*(pmom::gamma1gg_1[ii1] + pmom::gamma1qq_1[ii1]) + C1qq_1[ii1]*(pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1] + pow(pmom::gamma1qq_1[ii1],2)))*(LF2 + LQ2))/2.;
		H3st_qg_1 += C2qg_1[ii1]*pmom::gamma1qq_2[ii2]*LF + C1qg_1[ii1]*(C1qg_2[ii2]*pmom::gamma1gq_2[ii2] + C1qq_2[ii2]*pmom::gamma1qq_2[ii2])*LF + C2qq_2[ii2]*pmom::gamma1qg_1[ii1]*LF + C1qq_2[ii2]*(C1qg_1[ii1]*pmom::gamma1gg_1[ii1] + C1qq_1[ii1]*pmom::gamma1qg_1[ii1])*LF + C1qg_1[ii1]*pmom::gamma2qq_2[ii2]*LF + C1qq_2[ii2]*pmom::gamma2qg_1[ii1]*LF + (C1qg_1[ii1]*pmom::gamma2gg_1[ii1] + C1qq_1[ii1]*pmom::gamma2qg_1[ii1])*LF + pmom::gamma3qg_1[ii1]*LF + (beta0*C1qg_1[ii1]*pmom::gamma1qq_2[ii2]*LF2)/2. + (C1qg_1[ii1]*(pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2] + pow(pmom::gamma1qq_2[ii2],2))*LF2)/2. + (beta1*pmom::gamma1qg_1[ii1]*LF2)/2. + (beta0*C1qq_2[ii2]*pmom::gamma1qg_1[ii1]*LF2)/2. + (C1qg_2[ii2]*pmom::gamma1gq_2[ii2] + C1qq_2[ii2]*pmom::gamma1qq_2[ii2])*pmom::gamma1qg_1[ii1]*LF2 + (beta0*(C1qg_1[ii1]*pmom::gamma1gg_1[ii1] + C1qq_1[ii1]*pmom::gamma1qg_1[ii1])*LF2)/2. + pmom::gamma1qq_2[ii2]*(C1qg_1[ii1]*pmom::gamma1gg_1[ii1] + C1qq_1[ii1]*pmom::gamma1qg_1[ii1])*LF2 + (C1qq_2[ii2]*pmom::gamma1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + pmom::gamma1qq_1[ii1])*LF2)/2. + pmom::gamma1qg_1[ii1]*pmom::gamma2qq_2[ii2]*LF2 + beta0*pmom::gamma2qg_1[ii1]*LF2 + pmom::gamma1qq_2[ii2]*pmom::gamma2qg_1[ii1]*LF2 + (pmom::gamma1qg_1[ii1]*pmom::gamma2gg_1[ii1] + pmom::gamma1qq_1[ii1]*pmom::gamma2qg_1[ii1])* LF2 + (pow(beta0,2)*pmom::gamma1qg_1[ii1]*LF3)/3. + beta0*pmom::gamma1qq_2[ii2]*pmom::gamma1qg_1[ii1]*LF3 + ((pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2] + pow(pmom::gamma1qq_2[ii2],2))*pmom::gamma1qg_1[ii1]*LF3)/2. + (beta0*pmom::gamma1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + pmom::gamma1qq_1[ii1])*LF3)/2. + (pmom::gamma1qq_2[ii2]*pmom::gamma1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + pmom::gamma1qq_1[ii1])*LF3)/2. - C2qg_1[ii1]*pmom::gamma1qq_2[ii2]*LQ - C1qg_1[ii1]*(C1qg_2[ii2]*pmom::gamma1gq_2[ii2] + C1qq_2[ii2]*pmom::gamma1qq_2[ii2])*LQ - C2qq_2[ii2]*pmom::gamma1qg_1[ii1]*LQ - C1qq_2[ii2]*(C1qg_1[ii1]*pmom::gamma1gg_1[ii1] + C1qq_1[ii1]*pmom::gamma1qg_1[ii1])*LQ - C1qg_1[ii1]*pmom::gamma2qq_2[ii2]*LQ - C1qq_2[ii2]*pmom::gamma2qg_1[ii1]*LQ - (C1qg_1[ii1]*pmom::gamma2gg_1[ii1] + C1qq_1[ii1]*pmom::gamma2qg_1[ii1])*LQ - pmom::gamma3qg_1[ii1]*LQ - B1q*C1qg_1[ii1]*pmom::gamma1qq_2[ii2]*LF*LQ + beta0*C1qg_1[ii1]*pmom::gamma1qq_2[ii2]*LF*LQ - C1qg_1[ii1]*(pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2] + pow(pmom::gamma1qq_2[ii2],2))*LF*LQ - B2q*pmom::gamma1qg_1[ii1]*LF*LQ - B1q*C1qq_2[ii2]*pmom::gamma1qg_1[ii1]*LF*LQ + beta0*C1qq_2[ii2]*pmom::gamma1qg_1[ii1]*LF*LQ - 2.*(C1qg_2[ii2]*pmom::gamma1gq_2[ii2] + C1qq_2[ii2]*pmom::gamma1qq_2[ii2])*pmom::gamma1qg_1[ii1]*LF*LQ - B1q*(C1qg_1[ii1]*pmom::gamma1gg_1[ii1] + C1qq_1[ii1]*pmom::gamma1qg_1[ii1])*LF*LQ + beta0*(C1qg_1[ii1]*pmom::gamma1gg_1[ii1] + C1qq_1[ii1]*pmom::gamma1qg_1[ii1])*LF*LQ - 2.*pmom::gamma1qq_2[ii2]*(C1qg_1[ii1]*pmom::gamma1gg_1[ii1] + C1qq_1[ii1]*pmom::gamma1qg_1[ii1])*LF*LQ - C1qq_2[ii2]*pmom::gamma1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + pmom::gamma1qq_1[ii1])*LF*LQ - 2.*pmom::gamma1qg_1[ii1]*pmom::gamma2qq_2[ii2]*LF*LQ - B1q*pmom::gamma2qg_1[ii1]*LF*LQ - 2.*pmom::gamma1qq_2[ii2]*pmom::gamma2qg_1[ii1]*LF*LQ - (pmom::gamma1qg_1[ii1]*pmom::gamma2gg_1[ii1] + pmom::gamma1qq_1[ii1]*pmom::gamma2qg_1[ii1])*LF*LQ - (B1q*beta0*pmom::gamma1qg_1[ii1]*LF2*LQ)/2. - B1q*pmom::gamma1qq_2[ii2]*pmom::gamma1qg_1[ii1]*LF2*LQ - beta0*pmom::gamma1qq_2[ii2]*pmom::gamma1qg_1[ii1]*LF2*LQ - (3.*(pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2] + pow(pmom::gamma1qq_2[ii2],2))*pmom::gamma1qg_1[ii1]*LF2*LQ)/2. - (B1q*pmom::gamma1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + pmom::gamma1qq_1[ii1])*LF2*LQ)/2. - (beta0*pmom::gamma1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + pmom::gamma1qq_1[ii1])*LF2*LQ)/2. - (3.*pmom::gamma1qq_2[ii2]*pmom::gamma1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + pmom::gamma1qq_1[ii1])*LF2*LQ)/2. + B1q*C1qg_1[ii1]*pmom::gamma1qq_2[ii2]*LQ2 - (3.*beta0*C1qg_1[ii1]*pmom::gamma1qq_2[ii2]*LQ2)/2. + (C1qg_1[ii1]*(pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2] + pow(pmom::gamma1qq_2[ii2],2))*LQ2)/2. + B2q*pmom::gamma1qg_1[ii1]*LQ2 - (beta1*pmom::gamma1qg_1[ii1]*LQ2)/2. + B1q*C1qq_2[ii2]*pmom::gamma1qg_1[ii1]*LQ2 - (3.*beta0*C1qq_2[ii2]*pmom::gamma1qg_1[ii1]*LQ2)/2. + (C1qg_2[ii2]*pmom::gamma1gq_2[ii2] + C1qq_2[ii2]*pmom::gamma1qq_2[ii2])* pmom::gamma1qg_1[ii1]*LQ2 + B1q*(C1qg_1[ii1]*pmom::gamma1gg_1[ii1] + C1qq_1[ii1]*pmom::gamma1qg_1[ii1])*LQ2 - (3.*beta0*(C1qg_1[ii1]*pmom::gamma1gg_1[ii1] + C1qq_1[ii1]*pmom::gamma1qg_1[ii1])*LQ2)/2. + pmom::gamma1qq_2[ii2]*(C1qg_1[ii1]*pmom::gamma1gg_1[ii1] + C1qq_1[ii1]*pmom::gamma1qg_1[ii1])*LQ2 + (C1qq_2[ii2]*pmom::gamma1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + pmom::gamma1qq_1[ii1])*LQ2)/2. + pmom::gamma1qg_1[ii1]*pmom::gamma2qq_2[ii2]*LQ2 + B1q*pmom::gamma2qg_1[ii1]*LQ2 - beta0*pmom::gamma2qg_1[ii1]*LQ2 + pmom::gamma1qq_2[ii2]*pmom::gamma2qg_1[ii1]*LQ2 + (pmom::gamma1qg_1[ii1]*pmom::gamma2gg_1[ii1] + pmom::gamma1qq_1[ii1]*pmom::gamma2qg_1[ii1])* LQ2 - (A1q*C1qg_1[ii1]*pmom::gamma1qq_2[ii2]*LF*LQ2)/2. - (A2q*pmom::gamma1qg_1[ii1]*LF*LQ2)/2. + (pow(B1q,2)*pmom::gamma1qg_1[ii1]*LF*LQ2)/2. - (B1q*beta0*pmom::gamma1qg_1[ii1]*LF*LQ2)/2. - (A1q*C1qq_2[ii2]*pmom::gamma1qg_1[ii1]*LF*LQ2)/2. + 2.*B1q*pmom::gamma1qq_2[ii2]*pmom::gamma1qg_1[ii1]*LF*LQ2 - beta0*pmom::gamma1qq_2[ii2]*pmom::gamma1qg_1[ii1]*LF*LQ2 + (3.*(pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2] + pow(pmom::gamma1qq_2[ii2],2))*pmom::gamma1qg_1[ii1]*LF*LQ2)/2. - (A1q*(C1qg_1[ii1]*pmom::gamma1gg_1[ii1] + C1qq_1[ii1]*pmom::gamma1qg_1[ii1])*LF*LQ2)/2. + B1q*pmom::gamma1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + pmom::gamma1qq_1[ii1])*LF*LQ2 - (beta0*pmom::gamma1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + pmom::gamma1qq_1[ii1])*LF*LQ2)/2. + (3.*pmom::gamma1qq_2[ii2]*pmom::gamma1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + pmom::gamma1qq_1[ii1])*LF*LQ2)/2. - (A1q*pmom::gamma2qg_1[ii1]*LF*LQ2)/2. - (A1q*beta0*pmom::gamma1qg_1[ii1]*LF2*LQ2)/4. - (A1q*pmom::gamma1qq_2[ii2]*pmom::gamma1qg_1[ii1]*LF2*LQ2)/2. - (A1q*pmom::gamma1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + pmom::gamma1qq_1[ii1])*LF2*LQ2)/4. + (A1q*C1qg_1[ii1]*pmom::gamma1qq_2[ii2]*LQ3)/2. + (A2q*pmom::gamma1qg_1[ii1]*LQ3)/2. - (pow(B1q,2)*pmom::gamma1qg_1[ii1]*LQ3)/2. + B1q*beta0*pmom::gamma1qg_1[ii1]*LQ3 - (pow(beta0,2)*pmom::gamma1qg_1[ii1]*LQ3)/3. + (A1q*C1qq_2[ii2]*pmom::gamma1qg_1[ii1]*LQ3)/2. - B1q*pmom::gamma1qq_2[ii2]*pmom::gamma1qg_1[ii1]*LQ3 + beta0*pmom::gamma1qq_2[ii2]*pmom::gamma1qg_1[ii1]*LQ3 - ((pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2] + pow(pmom::gamma1qq_2[ii2],2))*pmom::gamma1qg_1[ii1]*LQ3)/2. + (A1q*(C1qg_1[ii1]*pmom::gamma1gg_1[ii1] + C1qq_1[ii1]*pmom::gamma1qg_1[ii1])*LQ3)/2. - (B1q*pmom::gamma1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + pmom::gamma1qq_1[ii1])*LQ3)/2. + (beta0*pmom::gamma1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + pmom::gamma1qq_1[ii1])*LQ3)/2. - (pmom::gamma1qq_2[ii2]*pmom::gamma1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + pmom::gamma1qq_1[ii1])*LQ3)/2. + (A1q*pmom::gamma2qg_1[ii1]*LQ3)/2. + (A1q*B1q*pmom::gamma1qg_1[ii1]*LF*LQ3)/2. - (A1q*beta0*pmom::gamma1qg_1[ii1]*LF*LQ3)/3. + A1q*pmom::gamma1qq_2[ii2]*pmom::gamma1qg_1[ii1]*LF*LQ3 + (A1q*pmom::gamma1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + pmom::gamma1qq_1[ii1])*LF*LQ3)/2. - (A1q*B1q*pmom::gamma1qg_1[ii1]*LQ4)/2. + (7.*A1q*beta0*pmom::gamma1qg_1[ii1]*LQ4)/12. - (A1q*pmom::gamma1qq_2[ii2]*pmom::gamma1qg_1[ii1]*LQ4)/2. - (A1q*pmom::gamma1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + pmom::gamma1qq_1[ii1])*LQ4)/4. + (pow(A1q,2)*pmom::gamma1qg_1[ii1]*LF*LQ4)/8. - (pow(A1q,2)*pmom::gamma1qg_1[ii1]*LQ5)/8. + LF*(C2qg_1[ii1]*pmom::gamma1gg_1[ii1] + pmom::gamma1qg_1[ii1]*(C2qq_1[ii1] + C2qqb_1[ii1] + (C2qqbp_1[ii1] + C2qqp_1[ii1])*(-1. + nf))) - LQ*(C2qg_1[ii1]*pmom::gamma1gg_1[ii1] + pmom::gamma1qg_1[ii1]*(C2qq_1[ii1] + C2qqb_1[ii1] + (C2qqbp_1[ii1] + C2qqp_1[ii1])*(-1. + nf))) - LF*LQ*(pmom::gamma1gg_1[ii1]*pmom::gamma2qg_1[ii1] + pmom::gamma1qg_1[ii1]*(pmom::gamma2qq_1[ii1] + pmom::gamma2qqb_1[ii1] + (pmom::gamma2qqbp_1[ii1] + pmom::gamma2qqp_1[ii1])*(-1. + nf))) + (pmom::gamma1qg_1[ii1]*LF3*(pow(pmom::gamma1gg_1[ii1],2) + pmom::gamma1gg_1[ii1]*pmom::gamma1qq_1[ii1] + pow(pmom::gamma1qq_1[ii1],2) + 2.*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*nf))/6. - (pmom::gamma1qg_1[ii1]*LF2*LQ*(pow(pmom::gamma1gg_1[ii1],2) + pmom::gamma1gg_1[ii1]*pmom::gamma1qq_1[ii1] + pow(pmom::gamma1qq_1[ii1],2) + 2.*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*nf))/2. + (pmom::gamma1qg_1[ii1]*LF*LQ2*(pow(pmom::gamma1gg_1[ii1],2) + pmom::gamma1gg_1[ii1]*pmom::gamma1qq_1[ii1] + pow(pmom::gamma1qq_1[ii1],2) + 2.*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*nf))/2. - (pmom::gamma1qg_1[ii1]*LQ3*(pow(pmom::gamma1gg_1[ii1],2) + pmom::gamma1gg_1[ii1]*pmom::gamma1qq_1[ii1] + pow(pmom::gamma1qq_1[ii1],2) + 2.*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*nf))/6. + (LF2*(C1qq_1[ii1]*pmom::gamma1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + pmom::gamma1qq_1[ii1]) + C1qg_1[ii1]*(pow(pmom::gamma1gg_1[ii1],2) + 2.*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*nf)))/2. - LF*LQ*(C1qq_1[ii1]*pmom::gamma1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + pmom::gamma1qq_1[ii1]) + C1qg_1[ii1]*(pow(pmom::gamma1gg_1[ii1],2) + 2.*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*nf)) + (LQ2*(C1qq_1[ii1]*pmom::gamma1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + pmom::gamma1qq_1[ii1]) + C1qg_1[ii1]*(pow(pmom::gamma1gg_1[ii1],2) + 2.*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*nf)))/2.;
		H3st_qg_2 += C2qg_2[ii2]*pmom::gamma1qq_1[ii1]*LF + C1qg_2[ii2]*(C1qg_1[ii1]*pmom::gamma1gq_1[ii1] + C1qq_1[ii1]*pmom::gamma1qq_1[ii1])*LF + C2qq_1[ii1]*pmom::gamma1qg_2[ii2]*LF + C1qq_1[ii1]*(C1qg_2[ii2]*pmom::gamma1gg_2[ii2] + C1qq_2[ii2]*pmom::gamma1qg_2[ii2])*LF + C1qg_2[ii2]*pmom::gamma2qq_1[ii1]*LF + C1qq_1[ii1]*pmom::gamma2qg_2[ii2]*LF + (C1qg_2[ii2]*pmom::gamma2gg_2[ii2] + C1qq_2[ii2]*pmom::gamma2qg_2[ii2])*LF + pmom::gamma3qg_2[ii2]*LF + (beta0*C1qg_2[ii2]*pmom::gamma1qq_1[ii1]*LF2)/2. + (C1qg_2[ii2]*(pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1] + pow(pmom::gamma1qq_1[ii1],2))*LF2)/2. + (beta1*pmom::gamma1qg_2[ii2]*LF2)/2. + (beta0*C1qq_1[ii1]*pmom::gamma1qg_2[ii2]*LF2)/2. + (C1qg_1[ii1]*pmom::gamma1gq_1[ii1] + C1qq_1[ii1]*pmom::gamma1qq_1[ii1])*pmom::gamma1qg_2[ii2]*LF2 + (beta0*(C1qg_2[ii2]*pmom::gamma1gg_2[ii2] + C1qq_2[ii2]*pmom::gamma1qg_2[ii2])*LF2)/2. + pmom::gamma1qq_1[ii1]*(C1qg_2[ii2]*pmom::gamma1gg_2[ii2] + C1qq_2[ii2]*pmom::gamma1qg_2[ii2])*LF2 + (C1qq_1[ii1]*pmom::gamma1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + pmom::gamma1qq_2[ii2])*LF2)/2. + pmom::gamma1qg_2[ii2]*pmom::gamma2qq_1[ii1]*LF2 + beta0*pmom::gamma2qg_2[ii2]*LF2 + pmom::gamma1qq_1[ii1]*pmom::gamma2qg_2[ii2]*LF2 + (pmom::gamma1qg_2[ii2]*pmom::gamma2gg_2[ii2] + pmom::gamma1qq_2[ii2]*pmom::gamma2qg_2[ii2])* LF2 + (pow(beta0,2)*pmom::gamma1qg_2[ii2]*LF3)/3. + beta0*pmom::gamma1qq_1[ii1]*pmom::gamma1qg_2[ii2]*LF3 + ((pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1] + pow(pmom::gamma1qq_1[ii1],2))*pmom::gamma1qg_2[ii2]*LF3)/2. + (beta0*pmom::gamma1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + pmom::gamma1qq_2[ii2])*LF3)/2. + (pmom::gamma1qq_1[ii1]*pmom::gamma1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + pmom::gamma1qq_2[ii2])*LF3)/2. - C2qg_2[ii2]*pmom::gamma1qq_1[ii1]*LQ - C1qg_2[ii2]*(C1qg_1[ii1]*pmom::gamma1gq_1[ii1] + C1qq_1[ii1]*pmom::gamma1qq_1[ii1])*LQ - C2qq_1[ii1]*pmom::gamma1qg_2[ii2]*LQ - C1qq_1[ii1]*(C1qg_2[ii2]*pmom::gamma1gg_2[ii2] + C1qq_2[ii2]*pmom::gamma1qg_2[ii2])*LQ - C1qg_2[ii2]*pmom::gamma2qq_1[ii1]*LQ - C1qq_1[ii1]*pmom::gamma2qg_2[ii2]*LQ - (C1qg_2[ii2]*pmom::gamma2gg_2[ii2] + C1qq_2[ii2]*pmom::gamma2qg_2[ii2])*LQ - pmom::gamma3qg_2[ii2]*LQ - B1q*C1qg_2[ii2]*pmom::gamma1qq_1[ii1]*LF*LQ + beta0*C1qg_2[ii2]*pmom::gamma1qq_1[ii1]*LF*LQ - C1qg_2[ii2]*(pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1] + pow(pmom::gamma1qq_1[ii1],2))*LF*LQ - B2q*pmom::gamma1qg_2[ii2]*LF*LQ - B1q*C1qq_1[ii1]*pmom::gamma1qg_2[ii2]*LF*LQ + beta0*C1qq_1[ii1]*pmom::gamma1qg_2[ii2]*LF*LQ - 2.*(C1qg_1[ii1]*pmom::gamma1gq_1[ii1] + C1qq_1[ii1]*pmom::gamma1qq_1[ii1])*pmom::gamma1qg_2[ii2]*LF*LQ - B1q*(C1qg_2[ii2]*pmom::gamma1gg_2[ii2] + C1qq_2[ii2]*pmom::gamma1qg_2[ii2])*LF*LQ + beta0*(C1qg_2[ii2]*pmom::gamma1gg_2[ii2] + C1qq_2[ii2]*pmom::gamma1qg_2[ii2])*LF*LQ - 2.*pmom::gamma1qq_1[ii1]*(C1qg_2[ii2]*pmom::gamma1gg_2[ii2] + C1qq_2[ii2]*pmom::gamma1qg_2[ii2])*LF*LQ - C1qq_1[ii1]*pmom::gamma1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + pmom::gamma1qq_2[ii2])*LF*LQ - 2.*pmom::gamma1qg_2[ii2]*pmom::gamma2qq_1[ii1]*LF*LQ - B1q*pmom::gamma2qg_2[ii2]*LF*LQ - 2.*pmom::gamma1qq_1[ii1]*pmom::gamma2qg_2[ii2]*LF*LQ - (pmom::gamma1qg_2[ii2]*pmom::gamma2gg_2[ii2] + pmom::gamma1qq_2[ii2]*pmom::gamma2qg_2[ii2])*LF*LQ - (B1q*beta0*pmom::gamma1qg_2[ii2]*LF2*LQ)/2. - B1q*pmom::gamma1qq_1[ii1]*pmom::gamma1qg_2[ii2]*LF2*LQ - beta0*pmom::gamma1qq_1[ii1]*pmom::gamma1qg_2[ii2]*LF2*LQ - (3.*(pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1] + pow(pmom::gamma1qq_1[ii1],2))*pmom::gamma1qg_2[ii2]*LF2*LQ)/2. - (B1q*pmom::gamma1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + pmom::gamma1qq_2[ii2])*LF2*LQ)/2. - (beta0*pmom::gamma1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + pmom::gamma1qq_2[ii2])*LF2*LQ)/2. - (3.*pmom::gamma1qq_1[ii1]*pmom::gamma1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + pmom::gamma1qq_2[ii2])*LF2*LQ)/2. + B1q*C1qg_2[ii2]*pmom::gamma1qq_1[ii1]*LQ2 - (3.*beta0*C1qg_2[ii2]*pmom::gamma1qq_1[ii1]*LQ2)/2. + (C1qg_2[ii2]*(pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1] + pow(pmom::gamma1qq_1[ii1],2))*LQ2)/2. + B2q*pmom::gamma1qg_2[ii2]*LQ2 - (beta1*pmom::gamma1qg_2[ii2]*LQ2)/2. + B1q*C1qq_1[ii1]*pmom::gamma1qg_2[ii2]*LQ2 - (3.*beta0*C1qq_1[ii1]*pmom::gamma1qg_2[ii2]*LQ2)/2. + (C1qg_1[ii1]*pmom::gamma1gq_1[ii1] + C1qq_1[ii1]*pmom::gamma1qq_1[ii1])* pmom::gamma1qg_2[ii2]*LQ2 + B1q*(C1qg_2[ii2]*pmom::gamma1gg_2[ii2] + C1qq_2[ii2]*pmom::gamma1qg_2[ii2])*LQ2 - (3.*beta0*(C1qg_2[ii2]*pmom::gamma1gg_2[ii2] + C1qq_2[ii2]*pmom::gamma1qg_2[ii2])*LQ2)/2. + pmom::gamma1qq_1[ii1]*(C1qg_2[ii2]*pmom::gamma1gg_2[ii2] + C1qq_2[ii2]*pmom::gamma1qg_2[ii2])*LQ2 + (C1qq_1[ii1]*pmom::gamma1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + pmom::gamma1qq_2[ii2])*LQ2)/2. + pmom::gamma1qg_2[ii2]*pmom::gamma2qq_1[ii1]*LQ2 + B1q*pmom::gamma2qg_2[ii2]*LQ2 - beta0*pmom::gamma2qg_2[ii2]*LQ2 + pmom::gamma1qq_1[ii1]*pmom::gamma2qg_2[ii2]*LQ2 + (pmom::gamma1qg_2[ii2]*pmom::gamma2gg_2[ii2] + pmom::gamma1qq_2[ii2]*pmom::gamma2qg_2[ii2])* LQ2 - (A1q*C1qg_2[ii2]*pmom::gamma1qq_1[ii1]*LF*LQ2)/2. - (A2q*pmom::gamma1qg_2[ii2]*LF*LQ2)/2. + (pow(B1q,2)*pmom::gamma1qg_2[ii2]*LF*LQ2)/2. - (B1q*beta0*pmom::gamma1qg_2[ii2]*LF*LQ2)/2. - (A1q*C1qq_1[ii1]*pmom::gamma1qg_2[ii2]*LF*LQ2)/2. + 2.*B1q*pmom::gamma1qq_1[ii1]*pmom::gamma1qg_2[ii2]*LF*LQ2 - beta0*pmom::gamma1qq_1[ii1]*pmom::gamma1qg_2[ii2]*LF*LQ2 + (3.*(pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1] + pow(pmom::gamma1qq_1[ii1],2))*pmom::gamma1qg_2[ii2]*LF*LQ2)/2. - (A1q*(C1qg_2[ii2]*pmom::gamma1gg_2[ii2] + C1qq_2[ii2]*pmom::gamma1qg_2[ii2])*LF*LQ2)/2. + B1q*pmom::gamma1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + pmom::gamma1qq_2[ii2])*LF*LQ2 - (beta0*pmom::gamma1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + pmom::gamma1qq_2[ii2])*LF*LQ2)/2. + (3.*pmom::gamma1qq_1[ii1]*pmom::gamma1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + pmom::gamma1qq_2[ii2])*LF*LQ2)/2. - (A1q*pmom::gamma2qg_2[ii2]*LF*LQ2)/2. - (A1q*beta0*pmom::gamma1qg_2[ii2]*LF2*LQ2)/4. - (A1q*pmom::gamma1qq_1[ii1]*pmom::gamma1qg_2[ii2]*LF2*LQ2)/2. - (A1q*pmom::gamma1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + pmom::gamma1qq_2[ii2])*LF2*LQ2)/4. + (A1q*C1qg_2[ii2]*pmom::gamma1qq_1[ii1]*LQ3)/2. + (A2q*pmom::gamma1qg_2[ii2]*LQ3)/2. - (pow(B1q,2)*pmom::gamma1qg_2[ii2]*LQ3)/2. + B1q*beta0*pmom::gamma1qg_2[ii2]*LQ3 - (pow(beta0,2)*pmom::gamma1qg_2[ii2]*LQ3)/3. + (A1q*C1qq_1[ii1]*pmom::gamma1qg_2[ii2]*LQ3)/2. - B1q*pmom::gamma1qq_1[ii1]*pmom::gamma1qg_2[ii2]*LQ3 + beta0*pmom::gamma1qq_1[ii1]*pmom::gamma1qg_2[ii2]*LQ3 - ((pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1] + pow(pmom::gamma1qq_1[ii1],2))*pmom::gamma1qg_2[ii2]*LQ3)/2. + (A1q*(C1qg_2[ii2]*pmom::gamma1gg_2[ii2] + C1qq_2[ii2]*pmom::gamma1qg_2[ii2])*LQ3)/2. - (B1q*pmom::gamma1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + pmom::gamma1qq_2[ii2])*LQ3)/2. + (beta0*pmom::gamma1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + pmom::gamma1qq_2[ii2])*LQ3)/2. - (pmom::gamma1qq_1[ii1]*pmom::gamma1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + pmom::gamma1qq_2[ii2])*LQ3)/2. + (A1q*pmom::gamma2qg_2[ii2]*LQ3)/2. + (A1q*B1q*pmom::gamma1qg_2[ii2]*LF*LQ3)/2. - (A1q*beta0*pmom::gamma1qg_2[ii2]*LF*LQ3)/3. + A1q*pmom::gamma1qq_1[ii1]*pmom::gamma1qg_2[ii2]*LF*LQ3 + (A1q*pmom::gamma1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + pmom::gamma1qq_2[ii2])*LF*LQ3)/2. - (A1q*B1q*pmom::gamma1qg_2[ii2]*LQ4)/2. + (7.*A1q*beta0*pmom::gamma1qg_2[ii2]*LQ4)/12. - (A1q*pmom::gamma1qq_1[ii1]*pmom::gamma1qg_2[ii2]*LQ4)/2. - (A1q*pmom::gamma1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + pmom::gamma1qq_2[ii2])*LQ4)/4. + (pow(A1q,2)*pmom::gamma1qg_2[ii2]*LF*LQ4)/8. - (pow(A1q,2)*pmom::gamma1qg_2[ii2]*LQ5)/8. + LF*(C2qg_2[ii2]*pmom::gamma1gg_2[ii2] + pmom::gamma1qg_2[ii2]*(C2qq_2[ii2] + C2qqb_2[ii2] + (C2qqbp_2[ii2] + C2qqp_2[ii2])*(-1. + nf))) - LQ*(C2qg_2[ii2]*pmom::gamma1gg_2[ii2] + pmom::gamma1qg_2[ii2]*(C2qq_2[ii2] + C2qqb_2[ii2] + (C2qqbp_2[ii2] + C2qqp_2[ii2])*(-1. + nf))) - LF*LQ*(pmom::gamma1gg_2[ii2]*pmom::gamma2qg_2[ii2] + pmom::gamma1qg_2[ii2]*(pmom::gamma2qq_2[ii2] + pmom::gamma2qqb_2[ii2] + (pmom::gamma2qqbp_2[ii2] + pmom::gamma2qqp_2[ii2])*(-1. + nf))) + (pmom::gamma1qg_2[ii2]*LF3*(pow(pmom::gamma1gg_2[ii2],2) + pmom::gamma1gg_2[ii2]*pmom::gamma1qq_2[ii2] + pow(pmom::gamma1qq_2[ii2],2) + 2.*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*nf))/6. - (pmom::gamma1qg_2[ii2]*LF2*LQ*(pow(pmom::gamma1gg_2[ii2],2) + pmom::gamma1gg_2[ii2]*pmom::gamma1qq_2[ii2] + pow(pmom::gamma1qq_2[ii2],2) + 2.*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*nf))/2. + (pmom::gamma1qg_2[ii2]*LF*LQ2*(pow(pmom::gamma1gg_2[ii2],2) + pmom::gamma1gg_2[ii2]*pmom::gamma1qq_2[ii2] + pow(pmom::gamma1qq_2[ii2],2) + 2.*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*nf))/2. - (pmom::gamma1qg_2[ii2]*LQ3*(pow(pmom::gamma1gg_2[ii2],2) + pmom::gamma1gg_2[ii2]*pmom::gamma1qq_2[ii2] + pow(pmom::gamma1qq_2[ii2],2) + 2.*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*nf))/6. + (LF2*(C1qq_2[ii2]*pmom::gamma1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + pmom::gamma1qq_2[ii2]) + C1qg_2[ii2]*(pow(pmom::gamma1gg_2[ii2],2) + 2.*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*nf)))/2. - LF*LQ*(C1qq_2[ii2]*pmom::gamma1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + pmom::gamma1qq_2[ii2]) + C1qg_2[ii2]*(pow(pmom::gamma1gg_2[ii2],2) + 2.*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*nf)) + (LQ2*(C1qq_2[ii2]*pmom::gamma1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + pmom::gamma1qq_2[ii2]) + C1qg_2[ii2]*(pow(pmom::gamma1gg_2[ii2],2) + 2.*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*nf)))/2.;
		H3st_qq_1 += C2qqb_1[ii1]*pmom::gamma1qq_2[ii2]*LF + C1qq_2[ii2]*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LF + (C2qg_1[ii1]*pmom::gamma1gq_1[ii1] + C2qqb_1[ii1]*pmom::gamma1qq_1[ii1])*LF + C1qq_2[ii2]*pmom::gamma2qqb_1[ii1]*LF + (C1qg_1[ii1]*pmom::gamma2gq_1[ii1] + C1qq_1[ii1]*pmom::gamma2qqb_1[ii1])*LF + pmom::gamma3qqb_1[ii1]*LF + (beta0*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LF2)/2. + C1qg_1[ii1]*pmom::gamma1qq_2[ii2]*pmom::gamma1gq_1[ii1]*LF2 + (C1qq_2[ii2]*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF2)/2. + (pmom::gamma1gq_1[ii1]*(C1qq_1[ii1]*pmom::gamma1qg_1[ii1] + C1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + pmom::gamma1qq_1[ii1]))*LF2)/2. + beta0*pmom::gamma2qqb_1[ii1]*LF2 + pmom::gamma1qq_2[ii2]*pmom::gamma2qqb_1[ii1]*LF2 + (pmom::gamma1qg_1[ii1]*pmom::gamma2gq_1[ii1] + pmom::gamma1qq_1[ii1]*pmom::gamma2qqb_1[ii1])*LF2 + (beta0*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF3)/2. + (pmom::gamma1qq_2[ii2]*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF3)/2. + (pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + 2.*pmom::gamma1qq_1[ii1])*LF3)/6. - C2qqb_1[ii1]*pmom::gamma1qq_2[ii2]*LQ - C1qq_2[ii2]*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LQ - (C2qg_1[ii1]*pmom::gamma1gq_1[ii1] + C2qqb_1[ii1]*pmom::gamma1qq_1[ii1])*LQ - C1qq_2[ii2]*pmom::gamma2qqb_1[ii1]*LQ - (C1qg_1[ii1]*pmom::gamma2gq_1[ii1] + C1qq_1[ii1]*pmom::gamma2qqb_1[ii1])*LQ - pmom::gamma3qqb_1[ii1]*LQ - B1q*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LF*LQ + beta0*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LF*LQ - 2.*C1qg_1[ii1]*pmom::gamma1qq_2[ii2]*pmom::gamma1gq_1[ii1]*LF*LQ - C1qq_2[ii2]*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF*LQ - pmom::gamma1gq_1[ii1]*(C1qq_1[ii1]*pmom::gamma1qg_1[ii1] + C1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + pmom::gamma1qq_1[ii1]))*LF*LQ - B1q*pmom::gamma2qqb_1[ii1]*LF*LQ - 2.*pmom::gamma1qq_2[ii2]*pmom::gamma2qqb_1[ii1]*LF*LQ - (pmom::gamma1qg_1[ii1]*pmom::gamma2gq_1[ii1] + pmom::gamma1qq_1[ii1]*pmom::gamma2qqb_1[ii1])*LF*LQ - (pmom::gamma1gq_1[ii1]*pmom::gamma2qg_1[ii1] + pmom::gamma1qq_1[ii1]*pmom::gamma2qqb_1[ii1])*LF*LQ - (B1q*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF2*LQ)/2. - (beta0*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF2*LQ)/ 2. - (3.*pmom::gamma1qq_2[ii2]*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF2*LQ)/2. - (pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + 2.*pmom::gamma1qq_1[ii1])*LF2*LQ)/2. + B1q*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LQ2 - (3.*beta0*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LQ2)/2. + C1qg_1[ii1]*pmom::gamma1qq_2[ii2]*pmom::gamma1gq_1[ii1]*LQ2 + (C1qq_2[ii2]*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LQ2)/2. + (pmom::gamma1gq_1[ii1]*(C1qq_1[ii1]*pmom::gamma1qg_1[ii1] + C1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + pmom::gamma1qq_1[ii1]))*LQ2)/2. + B1q*pmom::gamma2qqb_1[ii1]*LQ2 - beta0*pmom::gamma2qqb_1[ii1]*LQ2 + pmom::gamma1qq_2[ii2]*pmom::gamma2qqb_1[ii1]*LQ2 + (pmom::gamma1qg_1[ii1]*pmom::gamma2gq_1[ii1] + pmom::gamma1qq_1[ii1]*pmom::gamma2qqb_1[ii1])*LQ2 - (A1q*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LF*LQ2)/2. + B1q*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF*LQ2 - (beta0*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF*LQ2)/2. + (3.*pmom::gamma1qq_2[ii2]*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF*LQ2)/2. + (pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + 2.*pmom::gamma1qq_1[ii1])*LF*LQ2)/2. - (A1q*pmom::gamma2qqb_1[ii1]*LF*LQ2)/2. - (A1q*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF2*LQ2)/4. + (A1q*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LQ3)/2. - (B1q*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LQ3)/2. + (beta0*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LQ3)/2. - (pmom::gamma1qq_2[ii2]*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LQ3)/2. - (pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + 2.*pmom::gamma1qq_1[ii1])*LQ3)/6. + (A1q*pmom::gamma2qqb_1[ii1]*LQ3)/2. + (A1q*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF*LQ3)/2. - (A1q*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LQ4)/4.;
		H3st_qq_2 += C2qqb_2[ii2]*pmom::gamma1qq_1[ii1]*LF + C1qq_1[ii1]*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LF + (C2qg_2[ii2]*pmom::gamma1gq_2[ii2] + C2qqb_2[ii2]*pmom::gamma1qq_2[ii2])*LF + C1qq_1[ii1]*pmom::gamma2qqb_2[ii2]*LF + (C1qg_2[ii2]*pmom::gamma2gq_2[ii2] + C1qq_2[ii2]*pmom::gamma2qqb_2[ii2])*LF + pmom::gamma3qqb_2[ii2]*LF + (beta0*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LF2)/2. + C1qg_2[ii2]*pmom::gamma1qq_1[ii1]*pmom::gamma1gq_2[ii2]*LF2 + (C1qq_1[ii1]*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF2)/2. + (pmom::gamma1gq_2[ii2]*(C1qq_2[ii2]*pmom::gamma1qg_2[ii2] + C1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + pmom::gamma1qq_2[ii2]))*LF2)/2. + beta0*pmom::gamma2qqb_2[ii2]*LF2 + pmom::gamma1qq_1[ii1]*pmom::gamma2qqb_2[ii2]*LF2 + (pmom::gamma1qg_2[ii2]*pmom::gamma2gq_2[ii2] + pmom::gamma1qq_2[ii2]*pmom::gamma2qqb_2[ii2])*LF2 + (beta0*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF3)/2. + (pmom::gamma1qq_1[ii1]*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF3)/2. + (pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + 2.*pmom::gamma1qq_2[ii2])*LF3)/6. - C2qqb_2[ii2]*pmom::gamma1qq_1[ii1]*LQ - C1qq_1[ii1]*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LQ - (C2qg_2[ii2]*pmom::gamma1gq_2[ii2] + C2qqb_2[ii2]*pmom::gamma1qq_2[ii2])*LQ - C1qq_1[ii1]*pmom::gamma2qqb_2[ii2]*LQ - (C1qg_2[ii2]*pmom::gamma2gq_2[ii2] + C1qq_2[ii2]*pmom::gamma2qqb_2[ii2])*LQ - pmom::gamma3qqb_2[ii2]*LQ - B1q*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LF*LQ + beta0*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LF*LQ - 2.*C1qg_2[ii2]*pmom::gamma1qq_1[ii1]*pmom::gamma1gq_2[ii2]*LF*LQ - C1qq_1[ii1]*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF*LQ - pmom::gamma1gq_2[ii2]*(C1qq_2[ii2]*pmom::gamma1qg_2[ii2] + C1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + pmom::gamma1qq_2[ii2]))*LF*LQ - B1q*pmom::gamma2qqb_2[ii2]*LF*LQ - 2.*pmom::gamma1qq_1[ii1]*pmom::gamma2qqb_2[ii2]*LF*LQ - (pmom::gamma1qg_2[ii2]*pmom::gamma2gq_2[ii2] + pmom::gamma1qq_2[ii2]*pmom::gamma2qqb_2[ii2])*LF*LQ - (pmom::gamma1gq_2[ii2]*pmom::gamma2qg_2[ii2] + pmom::gamma1qq_2[ii2]*pmom::gamma2qqb_2[ii2])*LF*LQ - (B1q*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF2*LQ)/2. - (beta0*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF2*LQ)/ 2. - (3.*pmom::gamma1qq_1[ii1]*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF2*LQ)/2. - (pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + 2.*pmom::gamma1qq_2[ii2])*LF2*LQ)/2. + B1q*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LQ2 - (3.*beta0*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LQ2)/2. + C1qg_2[ii2]*pmom::gamma1qq_1[ii1]*pmom::gamma1gq_2[ii2]*LQ2 + (C1qq_1[ii1]*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LQ2)/2. + (pmom::gamma1gq_2[ii2]*(C1qq_2[ii2]*pmom::gamma1qg_2[ii2] + C1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + pmom::gamma1qq_2[ii2]))*LQ2)/2. + B1q*pmom::gamma2qqb_2[ii2]*LQ2 - beta0*pmom::gamma2qqb_2[ii2]*LQ2 + pmom::gamma1qq_1[ii1]*pmom::gamma2qqb_2[ii2]*LQ2 + (pmom::gamma1qg_2[ii2]*pmom::gamma2gq_2[ii2] + pmom::gamma1qq_2[ii2]*pmom::gamma2qqb_2[ii2])*LQ2 - (A1q*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LF*LQ2)/2. + B1q*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF*LQ2 - (beta0*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF*LQ2)/2. + (3.*pmom::gamma1qq_1[ii1]*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF*LQ2)/2. + (pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + 2.*pmom::gamma1qq_2[ii2])*LF*LQ2)/2. - (A1q*pmom::gamma2qqb_2[ii2]*LF*LQ2)/2. - (A1q*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF2*LQ2)/4. + (A1q*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LQ3)/2. - (B1q*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LQ3)/2. + (beta0*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LQ3)/2. - (pmom::gamma1qq_1[ii1]*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LQ3)/2. - (pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + 2.*pmom::gamma1qq_2[ii2])*LQ3)/6. + (A1q*pmom::gamma2qqb_2[ii2]*LQ3)/2. + (A1q*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF*LQ3)/2. - (A1q*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LQ4)/4.;
		H3st_qqp_1 += (2.*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + 2.*pmom::gamma1qq_1[ii1])*LQF3 + 6.*pmom::gamma1gq_1[ii1]*(C1qq_1[ii1]*pmom::gamma1qg_1[ii1] + C1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + pmom::gamma1qq_1[ii1]))* (LF2 + LQ2) - 3.*(-4.*C2qg_1[ii1]*pmom::gamma1gq_1[ii1]*LF - 4.*C1qg_1[ii1]*pmom::gamma2gq_1[ii1]*LF - 4.*C1qq_1[ii1]*pmom::gamma2qqbp_1[ii1]*LF - 4.*pmom::gamma3qqbp_1[ii1]*LF - 2.*beta0*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LF2 - 4.*C1qg_1[ii1]*pmom::gamma1qq_2[ii2]*pmom::gamma1gq_1[ii1]*LF2 - 4.*pmom::gamma1qg_1[ii1]*pmom::gamma2gq_1[ii1]*LF2 - 4.*beta0*pmom::gamma2qqbp_1[ii1]*LF2 - 4.*pmom::gamma1qq_2[ii2]*pmom::gamma2qqbp_1[ii1]*LF2 - 4.*pmom::gamma1qq_1[ii1]*pmom::gamma2qqbp_1[ii1]*LF2 - 2.*beta0*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF3 - 2.*pmom::gamma1qq_2[ii2]*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]* LF3 - 4.*C2qqbp_1[ii1]*(pmom::gamma1qq_2[ii2] + pmom::gamma1qq_1[ii1])*LQF - 2.*C1qq_2[ii2]*(2.*C1qg_1[ii1]*pmom::gamma1gq_1[ii1] + 2.*pmom::gamma2qqbp_1[ii1] + pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LQF)*LQF + 4.*C2qg_1[ii1]*pmom::gamma1gq_1[ii1]*LQ + 4.*C1qg_1[ii1]*pmom::gamma2gq_1[ii1]*LQ + 4.*C1qq_1[ii1]*pmom::gamma2qqbp_1[ii1]*LQ + 4.*pmom::gamma3qqbp_1[ii1]*LQ + 4.*B1q*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LF*LQ - 4.*beta0*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LF*LQ + 8.*C1qg_1[ii1]*pmom::gamma1qq_2[ii2]*pmom::gamma1gq_1[ii1]*LF*LQ + 4.*C1qg_1[ii1]*pmom::gamma1gg_1[ii1]*pmom::gamma1gq_1[ii1]*LF* LQ + 4.*C1qq_1[ii1]*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF*LQ + 4.*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]* pmom::gamma1qq_1[ii1]*LF*LQ + 4.*pmom::gamma1qg_1[ii1]*pmom::gamma2gq_1[ii1]*LF*LQ + 4.*pmom::gamma1gq_1[ii1]*pmom::gamma2qg_1[ii1]*LF*LQ + 4.*B1q*pmom::gamma2qqbp_1[ii1]*LF*LQ + 8.*pmom::gamma1qq_2[ii2]*pmom::gamma2qqbp_1[ii1]*LF*LQ + 8.*pmom::gamma1qq_1[ii1]*pmom::gamma2qqbp_1[ii1]*LF*LQ + 2.*B1q*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF2*LQ + 2.*beta0*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF2* LQ + 6.*pmom::gamma1qq_2[ii2]*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF2*LQ - 4.*B1q*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LQ2 + 6.*beta0*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LQ2 - 4.*C1qg_1[ii1]*pmom::gamma1qq_2[ii2]*pmom::gamma1gq_1[ii1]*LQ2 - 4.*pmom::gamma1qg_1[ii1]*pmom::gamma2gq_1[ii1]*LQ2 - 4.*B1q*pmom::gamma2qqbp_1[ii1]*LQ2 + 4.*beta0*pmom::gamma2qqbp_1[ii1]*LQ2 - 4.*pmom::gamma1qq_2[ii2]*pmom::gamma2qqbp_1[ii1]*LQ2 - 4.*pmom::gamma1qq_1[ii1]*pmom::gamma2qqbp_1[ii1]*LQ2 + 2.*A1q*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LF*LQ2 - 4.*B1q*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF*LQ2 + 2.*beta0*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF*LQ2 - 6.*pmom::gamma1qq_2[ii2]*pmom::gamma1gq_1[ii1]* pmom::gamma1qg_1[ii1]*LF*LQ2 + 2.*A1q*pmom::gamma2qqbp_1[ii1]*LF*LQ2 + A1q*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF2*LQ2 - 2.*A1q*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LQ3 + 2.*B1q*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LQ3 - 2.*beta0*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LQ3 + 2.*pmom::gamma1qq_2[ii2]*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LQ3 - 2.*A1q*pmom::gamma2qqbp_1[ii1]*LQ3 - 2.*A1q*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF*LQ3 + A1q*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LQ4))/12.;
		H3st_qqp_2 += (2.*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + 2.*pmom::gamma1qq_2[ii2])*LQF3 + 6.*pmom::gamma1gq_2[ii2]*(C1qq_2[ii2]*pmom::gamma1qg_2[ii2] + C1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + pmom::gamma1qq_2[ii2]))* (LF2 + LQ2) - 3.*(-4.*C2qg_2[ii2]*pmom::gamma1gq_2[ii2]*LF - 4.*C1qg_2[ii2]*pmom::gamma2gq_2[ii2]*LF - 4.*C1qq_2[ii2]*pmom::gamma2qqbp_2[ii2]*LF - 4.*pmom::gamma3qqbp_2[ii2]*LF - 2.*beta0*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LF2 - 4.*C1qg_2[ii2]*pmom::gamma1qq_1[ii1]*pmom::gamma1gq_2[ii2]*LF2 - 4.*pmom::gamma1qg_2[ii2]*pmom::gamma2gq_2[ii2]*LF2 - 4.*beta0*pmom::gamma2qqbp_2[ii2]*LF2 - 4.*pmom::gamma1qq_1[ii1]*pmom::gamma2qqbp_2[ii2]*LF2 - 4.*pmom::gamma1qq_2[ii2]*pmom::gamma2qqbp_2[ii2]*LF2 - 2.*beta0*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF3 - 2.*pmom::gamma1qq_1[ii1]*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]* LF3 - 4.*C2qqbp_2[ii2]*(pmom::gamma1qq_1[ii1] + pmom::gamma1qq_2[ii2])*LQF - 2.*C1qq_1[ii1]*(2.*C1qg_2[ii2]*pmom::gamma1gq_2[ii2] + 2.*pmom::gamma2qqbp_2[ii2] + pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LQF)*LQF + 4.*C2qg_2[ii2]*pmom::gamma1gq_2[ii2]*LQ + 4.*C1qg_2[ii2]*pmom::gamma2gq_2[ii2]*LQ + 4.*C1qq_2[ii2]*pmom::gamma2qqbp_2[ii2]*LQ + 4.*pmom::gamma3qqbp_2[ii2]*LQ + 4.*B1q*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LF*LQ - 4.*beta0*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LF*LQ + 8.*C1qg_2[ii2]*pmom::gamma1qq_1[ii1]*pmom::gamma1gq_2[ii2]*LF*LQ + 4.*C1qg_2[ii2]*pmom::gamma1gg_2[ii2]*pmom::gamma1gq_2[ii2]*LF* LQ + 4.*C1qq_2[ii2]*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF*LQ + 4.*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]* pmom::gamma1qq_2[ii2]*LF*LQ + 4.*pmom::gamma1qg_2[ii2]*pmom::gamma2gq_2[ii2]*LF*LQ + 4.*pmom::gamma1gq_2[ii2]*pmom::gamma2qg_2[ii2]*LF*LQ + 4.*B1q*pmom::gamma2qqbp_2[ii2]*LF*LQ + 8.*pmom::gamma1qq_1[ii1]*pmom::gamma2qqbp_2[ii2]*LF*LQ + 8.*pmom::gamma1qq_2[ii2]*pmom::gamma2qqbp_2[ii2]*LF*LQ + 2.*B1q*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF2*LQ + 2.*beta0*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF2* LQ + 6.*pmom::gamma1qq_1[ii1]*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF2*LQ - 4.*B1q*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LQ2 + 6.*beta0*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LQ2 - 4.*C1qg_2[ii2]*pmom::gamma1qq_1[ii1]*pmom::gamma1gq_2[ii2]*LQ2 - 4.*pmom::gamma1qg_2[ii2]*pmom::gamma2gq_2[ii2]*LQ2 - 4.*B1q*pmom::gamma2qqbp_2[ii2]*LQ2 + 4.*beta0*pmom::gamma2qqbp_2[ii2]*LQ2 - 4.*pmom::gamma1qq_1[ii1]*pmom::gamma2qqbp_2[ii2]*LQ2 - 4.*pmom::gamma1qq_2[ii2]*pmom::gamma2qqbp_2[ii2]*LQ2 + 2.*A1q*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LF*LQ2 - 4.*B1q*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF*LQ2 + 2.*beta0*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF*LQ2 - 6.*pmom::gamma1qq_1[ii1]*pmom::gamma1gq_2[ii2]* pmom::gamma1qg_2[ii2]*LF*LQ2 + 2.*A1q*pmom::gamma2qqbp_2[ii2]*LF*LQ2 + A1q*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF2*LQ2 - 2.*A1q*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LQ3 + 2.*B1q*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LQ3 - 2.*beta0*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LQ3 + 2.*pmom::gamma1qq_1[ii1]*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LQ3 - 2.*A1q*pmom::gamma2qqbp_2[ii2]*LQ3 - 2.*A1q*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF*LQ3 + A1q*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LQ4))/12.;
		H3st_qqbp_1 += (2.*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + 2.*pmom::gamma1qq_1[ii1])*LQF3 + 6.*pmom::gamma1gq_1[ii1]*(C1qq_1[ii1]*pmom::gamma1qg_1[ii1] + C1qg_1[ii1]*(pmom::gamma1gg_1[ii1] + pmom::gamma1qq_1[ii1]))* (LF2 + LQ2) - 3.*(-4.*C2qg_1[ii1]*pmom::gamma1gq_1[ii1]*LF - 4.*C1qg_1[ii1]*pmom::gamma2gq_1[ii1]*LF - 4.*C1qq_1[ii1]*pmom::gamma2qqp_1[ii1]*LF - 4.*pmom::gamma3qqp_1[ii1]*LF - 2.*beta0*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]* LF2 - 4.*C1qg_1[ii1]*pmom::gamma1qq_2[ii2]*pmom::gamma1gq_1[ii1]*LF2 - 4.*pmom::gamma1qg_1[ii1]*pmom::gamma2gq_1[ii1]* LF2 - 4.*beta0*pmom::gamma2qqp_1[ii1]*LF2 - 4.*pmom::gamma1qq_2[ii2]*pmom::gamma2qqp_1[ii1]*LF2 - 4.*pmom::gamma1qq_1[ii1]*pmom::gamma2qqp_1[ii1]*LF2 - 2.*beta0*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF3 - 2.*pmom::gamma1qq_2[ii2]*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF3 - 4.*C2qqp_1[ii1]*(pmom::gamma1qq_2[ii2] + pmom::gamma1qq_1[ii1])*LQF - 2.*C1qq_2[ii2]*(2.*C1qg_1[ii1]*pmom::gamma1gq_1[ii1] + 2.*pmom::gamma2qqp_1[ii1] + pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]* LQF)*LQF + 4.*C2qg_1[ii1]*pmom::gamma1gq_1[ii1]*LQ + 4.*C1qg_1[ii1]*pmom::gamma2gq_1[ii1]*LQ + 4.*C1qq_1[ii1]*pmom::gamma2qqp_1[ii1]*LQ + 4.*pmom::gamma3qqp_1[ii1]*LQ + 4.*B1q*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LF*LQ - 4.*beta0*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LF*LQ + 8.*C1qg_1[ii1]*pmom::gamma1qq_2[ii2]*pmom::gamma1gq_1[ii1]*LF*LQ + 4.*C1qg_1[ii1]*pmom::gamma1gg_1[ii1]*pmom::gamma1gq_1[ii1]*LF* LQ + 4.*C1qq_1[ii1]*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF*LQ + 4.*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]* pmom::gamma1qq_1[ii1]*LF*LQ + 4.*pmom::gamma1qg_1[ii1]*pmom::gamma2gq_1[ii1]*LF*LQ + 4.*pmom::gamma1gq_1[ii1]*pmom::gamma2qg_1[ii1]*LF*LQ + 4.*B1q*pmom::gamma2qqp_1[ii1]*LF*LQ + 8.*pmom::gamma1qq_2[ii2]*pmom::gamma2qqp_1[ii1]*LF*LQ + 8.*pmom::gamma1qq_1[ii1]*pmom::gamma2qqp_1[ii1]*LF*LQ + 2.*B1q*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF2*LQ + 2.*beta0*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF2* LQ + 6.*pmom::gamma1qq_2[ii2]*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF2*LQ - 4.*B1q*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LQ2 + 6.*beta0*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LQ2 - 4.*C1qg_1[ii1]*pmom::gamma1qq_2[ii2]*pmom::gamma1gq_1[ii1]*LQ2 - 4.*pmom::gamma1qg_1[ii1]*pmom::gamma2gq_1[ii1]*LQ2 - 4.*B1q*pmom::gamma2qqp_1[ii1]*LQ2 + 4.*beta0*pmom::gamma2qqp_1[ii1]*LQ2 - 4.*pmom::gamma1qq_2[ii2]*pmom::gamma2qqp_1[ii1]*LQ2 - 4.*pmom::gamma1qq_1[ii1]*pmom::gamma2qqp_1[ii1]*LQ2 + 2.*A1q*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LF*LQ2 - 4.*B1q*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF*LQ2 + 2.*beta0*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF*LQ2 - 6.*pmom::gamma1qq_2[ii2]*pmom::gamma1gq_1[ii1]* pmom::gamma1qg_1[ii1]*LF*LQ2 + 2.*A1q*pmom::gamma2qqp_1[ii1]*LF*LQ2 + A1q*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF2*LQ2 - 2.*A1q*C1qg_1[ii1]*pmom::gamma1gq_1[ii1]*LQ3 + 2.*B1q*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LQ3 - 2.*beta0*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LQ3 + 2.*pmom::gamma1qq_2[ii2]*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LQ3 - 2.*A1q*pmom::gamma2qqp_1[ii1]*LQ3 - 2.*A1q*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LF*LQ3 + A1q*pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LQ4))/12.;
		H3st_qqbp_2 += (2.*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + 2.*pmom::gamma1qq_2[ii2])*LQF3 + 6.*pmom::gamma1gq_2[ii2]*(C1qq_2[ii2]*pmom::gamma1qg_2[ii2] + C1qg_2[ii2]*(pmom::gamma1gg_2[ii2] + pmom::gamma1qq_2[ii2]))* (LF2 + LQ2) - 3.*(-4.*C2qg_2[ii2]*pmom::gamma1gq_2[ii2]*LF - 4.*C1qg_2[ii2]*pmom::gamma2gq_2[ii2]*LF - 4.*C1qq_2[ii2]*pmom::gamma2qqp_2[ii2]*LF - 4.*pmom::gamma3qqp_2[ii2]*LF - 2.*beta0*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]* LF2 - 4.*C1qg_2[ii2]*pmom::gamma1qq_1[ii1]*pmom::gamma1gq_2[ii2]*LF2 - 4.*pmom::gamma1qg_2[ii2]*pmom::gamma2gq_2[ii2]* LF2 - 4.*beta0*pmom::gamma2qqp_2[ii2]*LF2 - 4.*pmom::gamma1qq_1[ii1]*pmom::gamma2qqp_2[ii2]*LF2 - 4.*pmom::gamma1qq_2[ii2]*pmom::gamma2qqp_2[ii2]*LF2 - 2.*beta0*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF3 - 2.*pmom::gamma1qq_1[ii1]*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF3 - 4.*C2qqp_2[ii2]*(pmom::gamma1qq_1[ii1] + pmom::gamma1qq_2[ii2])*LQF - 2.*C1qq_1[ii1]*(2.*C1qg_2[ii2]*pmom::gamma1gq_2[ii2] + 2.*pmom::gamma2qqp_2[ii2] + pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]* LQF)*LQF + 4.*C2qg_2[ii2]*pmom::gamma1gq_2[ii2]*LQ + 4.*C1qg_2[ii2]*pmom::gamma2gq_2[ii2]*LQ + 4.*C1qq_2[ii2]*pmom::gamma2qqp_2[ii2]*LQ + 4.*pmom::gamma3qqp_2[ii2]*LQ + 4.*B1q*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LF*LQ - 4.*beta0*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LF*LQ + 8.*C1qg_2[ii2]*pmom::gamma1qq_1[ii1]*pmom::gamma1gq_2[ii2]*LF*LQ + 4.*C1qg_2[ii2]*pmom::gamma1gg_2[ii2]*pmom::gamma1gq_2[ii2]*LF* LQ + 4.*C1qq_2[ii2]*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF*LQ + 4.*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]* pmom::gamma1qq_2[ii2]*LF*LQ + 4.*pmom::gamma1qg_2[ii2]*pmom::gamma2gq_2[ii2]*LF*LQ + 4.*pmom::gamma1gq_2[ii2]*pmom::gamma2qg_2[ii2]*LF*LQ + 4.*B1q*pmom::gamma2qqp_2[ii2]*LF*LQ + 8.*pmom::gamma1qq_1[ii1]*pmom::gamma2qqp_2[ii2]*LF*LQ + 8.*pmom::gamma1qq_2[ii2]*pmom::gamma2qqp_2[ii2]*LF*LQ + 2.*B1q*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF2*LQ + 2.*beta0*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF2* LQ + 6.*pmom::gamma1qq_1[ii1]*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF2*LQ - 4.*B1q*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LQ2 + 6.*beta0*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LQ2 - 4.*C1qg_2[ii2]*pmom::gamma1qq_1[ii1]*pmom::gamma1gq_2[ii2]*LQ2 - 4.*pmom::gamma1qg_2[ii2]*pmom::gamma2gq_2[ii2]*LQ2 - 4.*B1q*pmom::gamma2qqp_2[ii2]*LQ2 + 4.*beta0*pmom::gamma2qqp_2[ii2]*LQ2 - 4.*pmom::gamma1qq_1[ii1]*pmom::gamma2qqp_2[ii2]*LQ2 - 4.*pmom::gamma1qq_2[ii2]*pmom::gamma2qqp_2[ii2]*LQ2 + 2.*A1q*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LF*LQ2 - 4.*B1q*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF*LQ2 + 2.*beta0*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF*LQ2 - 6.*pmom::gamma1qq_1[ii1]*pmom::gamma1gq_2[ii2]* pmom::gamma1qg_2[ii2]*LF*LQ2 + 2.*A1q*pmom::gamma2qqp_2[ii2]*LF*LQ2 + A1q*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF2*LQ2 - 2.*A1q*C1qg_2[ii2]*pmom::gamma1gq_2[ii2]*LQ3 + 2.*B1q*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LQ3 - 2.*beta0*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LQ3 + 2.*pmom::gamma1qq_1[ii1]*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LQ3 - 2.*A1q*pmom::gamma2qqp_2[ii2]*LQ3 - 2.*A1q*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LF*LQ3 + A1q*pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LQ4))/12.;
		H3st_gg += (LQF*(2.*C2qg_2[ii2]*pmom::gamma1qg_1[ii1] + 2.*C2qg_1[ii1]*pmom::gamma1qg_2[ii2] + 2.*C1qg_2[ii2]*pmom::gamma2qg_1[ii1] + beta0*C1qg_2[ii2]*pmom::gamma1qg_1[ii1]*LF + C1qg_2[ii2]*pmom::gamma1gg_1[ii1]*pmom::gamma1qg_1[ii1]*LF + C1qg_2[ii2]*pmom::gamma1qg_1[ii1]*pmom::gamma1qq_1[ii1]*LF + 2.*C1qg_2[ii2]*pmom::gamma1qg_1[ii1]*pmom::gamma1gg_2[ii2]*LF + 2.*C1qq_2[ii2]*pmom::gamma1qg_1[ii1]*pmom::gamma1qg_2[ii2]*LF + 2.*pmom::gamma1qg_2[ii2]*pmom::gamma2qg_1[ii1]*LF + 2.*pmom::gamma1qg_1[ii1]*pmom::gamma2qg_2[ii2]*LF + 2.*beta0*pmom::gamma1qg_1[ii1]*pmom::gamma1qg_2[ii2]*LF2 + pmom::gamma1gg_1[ii1]*pmom::gamma1qg_1[ii1]*pmom::gamma1qg_2[ii2]*LF2 + pmom::gamma1qg_1[ii1]*pmom::gamma1qq_1[ii1]*pmom::gamma1qg_2[ii2]* LF2 + pmom::gamma1qg_1[ii1]*pmom::gamma1gg_2[ii2]*pmom::gamma1qg_2[ii2]*LF2 + pmom::gamma1qg_1[ii1]*pmom::gamma1qg_2[ii2]*pmom::gamma1qq_2[ii2]*LF2 - 2.*B1q*C1qg_2[ii2]*pmom::gamma1qg_1[ii1]*LQ + 3.*beta0*C1qg_2[ii2]*pmom::gamma1qg_1[ii1]*LQ - C1qg_2[ii2]*pmom::gamma1gg_1[ii1]*pmom::gamma1qg_1[ii1]*LQ - C1qg_2[ii2]*pmom::gamma1qg_1[ii1]*pmom::gamma1qq_1[ii1]*LQ - 2.*C1qg_2[ii2]*pmom::gamma1qg_1[ii1]*pmom::gamma1gg_2[ii2]*LQ - 2.*C1qq_2[ii2]*pmom::gamma1qg_1[ii1]*pmom::gamma1qg_2[ii2]*LQ - 2.*pmom::gamma1qg_2[ii2]*pmom::gamma2qg_1[ii1]*LQ - 2.*pmom::gamma1qg_1[ii1]*pmom::gamma2qg_2[ii2]*LQ - 2.*B1q*pmom::gamma1qg_1[ii1]*pmom::gamma1qg_2[ii2]*LF*LQ - 2.*pmom::gamma1gg_1[ii1]*pmom::gamma1qg_1[ii1]*pmom::gamma1qg_2[ii2]*LF*LQ - 2.*pmom::gamma1qg_1[ii1]*pmom::gamma1qq_1[ii1]* pmom::gamma1qg_2[ii2]*LF*LQ - 2.*pmom::gamma1qg_1[ii1]*pmom::gamma1gg_2[ii2]*pmom::gamma1qg_2[ii2]*LF*LQ - 2.*pmom::gamma1qg_1[ii1]*pmom::gamma1qg_2[ii2]*pmom::gamma1qq_2[ii2]*LF*LQ - A1q*C1qg_2[ii2]*pmom::gamma1qg_1[ii1]*LQ2 + 2.*B1q*pmom::gamma1qg_1[ii1]*pmom::gamma1qg_2[ii2]*LQ2 - 2.*beta0*pmom::gamma1qg_1[ii1]*pmom::gamma1qg_2[ii2]*LQ2 + pmom::gamma1gg_1[ii1]*pmom::gamma1qg_1[ii1]*pmom::gamma1qg_2[ii2]*LQ2 + pmom::gamma1qg_1[ii1]*pmom::gamma1qq_1[ii1]*pmom::gamma1qg_2[ii2]* LQ2 + pmom::gamma1qg_1[ii1]*pmom::gamma1gg_2[ii2]*pmom::gamma1qg_2[ii2]*LQ2 + pmom::gamma1qg_1[ii1]*pmom::gamma1qg_2[ii2]*pmom::gamma1qq_2[ii2]*LQ2 - A1q*pmom::gamma1qg_1[ii1]*pmom::gamma1qg_2[ii2]*LF*LQ2 + A1q*pmom::gamma1qg_1[ii1]*pmom::gamma1qg_2[ii2]*LQ3 + 2.*C1qq_1[ii1]*pmom::gamma1qg_1[ii1]* (C1qg_2[ii2] + pmom::gamma1qg_2[ii2]*LF - pmom::gamma1qg_2[ii2]*LQ) + C1qg_1[ii1]*(2.*C1qg_2[ii2]*(pmom::gamma1gg_1[ii1] + pmom::gamma1gg_2[ii2]) + 2.*C1qq_2[ii2]*pmom::gamma1qg_2[ii2] + 2.*pmom::gamma2qg_2[ii2] + beta0*pmom::gamma1qg_2[ii2]*LF + 2.*pmom::gamma1gg_1[ii1]*pmom::gamma1qg_2[ii2]*LF + pmom::gamma1gg_2[ii2]*pmom::gamma1qg_2[ii2]*LF + pmom::gamma1qg_2[ii2]*pmom::gamma1qq_2[ii2]*LF - 2.*B1q*pmom::gamma1qg_2[ii2]*LQ + 3.*beta0*pmom::gamma1qg_2[ii2]*LQ - 2.*pmom::gamma1gg_1[ii1]*pmom::gamma1qg_2[ii2]* LQ - pmom::gamma1gg_2[ii2]*pmom::gamma1qg_2[ii2]*LQ - pmom::gamma1qg_2[ii2]*pmom::gamma1qq_2[ii2]*LQ - A1q*pmom::gamma1qg_2[ii2]*LQ2)))/2.;
		H3st_qbg_1 += ((2.*C2qqb_2[ii2]*pmom::gamma1qg_1[ii1] + (2.*C1qg_2[ii2]*pmom::gamma1gq_2[ii2] + 2.*pmom::gamma2qqb_2[ii2] + pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LQF)*(C1qg_1[ii1] + pmom::gamma1qg_1[ii1]*LQF))* LQF)/2.;
		H3st_qbg_2 += ((2.*C2qqb_1[ii1]*pmom::gamma1qg_2[ii2] + (2.*C1qg_1[ii1]*pmom::gamma1gq_1[ii1] + 2.*pmom::gamma2qqb_1[ii1] + pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LQF)*(C1qg_2[ii2] + pmom::gamma1qg_2[ii2]*LQF))* LQF)/2.;
		H3st_qpg_1 += ((2.*C2qqp_2[ii2]*pmom::gamma1qg_1[ii1] + (2.*C1qg_2[ii2]*pmom::gamma1gq_2[ii2] + 2.*pmom::gamma2qqp_2[ii2] + pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LQF)*(C1qg_1[ii1] + pmom::gamma1qg_1[ii1]*LQF))* LQF)/2.;
		H3st_qpg_2 += ((2.*C2qqp_1[ii1]*pmom::gamma1qg_2[ii2] + (2.*C1qg_1[ii1]*pmom::gamma1gq_1[ii1] + 2.*pmom::gamma2qqp_1[ii1] + pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LQF)*(C1qg_2[ii2] + pmom::gamma1qg_2[ii2]*LQF))* LQF)/2.;
		H3st_qbpg_1 += ((2.*C2qqbp_2[ii2]*pmom::gamma1qg_1[ii1] + (2.*C1qg_2[ii2]*pmom::gamma1gq_2[ii2] + 2.*pmom::gamma2qqbp_2[ii2] + pmom::gamma1gq_2[ii2]*pmom::gamma1qg_2[ii2]*LQF)*(C1qg_1[ii1] + pmom::gamma1qg_1[ii1]*LQF))* LQF)/2.;
		H3st_qbpg_2 += ((2.*C2qqbp_1[ii1]*pmom::gamma1qg_2[ii2] + (2.*C1qg_1[ii1]*pmom::gamma1gq_1[ii1] + 2.*pmom::gamma2qqbp_1[ii1] + pmom::gamma1gq_1[ii1]*pmom::gamma1qg_1[ii1]*LQF)*(C1qg_2[ii2] + pmom::gamma1qg_2[ii2]*LQF))* LQF)/2.;
	      }
	    /*
	    cout << endl;
	    cout << "hcoeff " << endl;
	    cout << "H3st_qqb  " << H3st_qqb  << endl;
	    cout << "H3st_qg_1   " << H3st_qg_1   << endl;
	    cout << "H3st_qq_1   " << H3st_qq_1   << endl;
	    cout << "H3st_qqp_1  " << H3st_qqp_1  << endl;
	    cout << "H3st_qqbp_1 " << H3st_qqbp_1 << endl;
	    cout << "H3st_gg   " << H3st_gg   << endl;
	    cout << "H3st_qbg_1  " << H3st_qbg_1  << endl;
	    cout << "H3st_qpg_1  " << H3st_qpg_1  << endl;
	    cout << "H3st_qbpg_1 " << H3st_qbpg_1 << endl;
	    */
	    
	    Hqqb[idx]    += as3*H3st_qqb;
	    Hqg_1[idx]   += as3*H3st_qg_1;
	    Hqg_2[idx]   += as3*H3st_qg_2;
	    Hqq_1[idx]   += as3*H3st_qq_1;
	    Hqq_2[idx]   += as3*H3st_qq_2;
	    Hqqp_1[idx]  += as3*H3st_qqp_1;
	    Hqqp_2[idx]  += as3*H3st_qqp_2;
	    Hqqbp_1[idx] += as3*H3st_qqbp_1;
	    Hqqbp_2[idx] += as3*H3st_qqbp_2;
	    Hgg[idx]     += as3*H3st_gg;
	    Hqbg_1[idx]  += as3*H3st_qbg_1;
	    Hqbg_2[idx]  += as3*H3st_qbg_2;
	    Hqpg_1[idx]  += as3*H3st_qpg_1;
	    Hqpg_2[idx]  += as3*H3st_qpg_2;
	    Hqbpg_1[idx] += as3*H3st_qbpg_1;
	    Hqbpg_2[idx] += as3*H3st_qbpg_2;
	  }
  
  return;

  /******  ALL THE FOLLOWING CAN BE CLEANED UP *****/
  // logs of scales are computed in resint
  //  complex <double> logmuf2q2 = resint::logmuf2q2;
  //  complex <double> logq2muf2 = resint::logq2muf2;
  //  complex <double> logq2mur2 = resint::logq2mur2;
  //  complex <double> loga = resint::loga;
  //  double aass = resint::aass;

  
  //All the following coefficients need to be calculated at each resumm iteration only with fixed factorization and renormalization scale (variable logmuf2q2) or
  //with fixed resummation scale (variable loga: a = q2/mu_res). Otherwise they can be computed at initialization
  /*
  if (opts.order == 0)
    fill(Hqqb, Hqqb+mellinint::mdim, 1.);
  */
  //for (int i = 0; i < mellinint::mdim; i++)
  //Hqqb[i] = 1.;
  
  //complex <double> LQF = logq2muf2-2.*loga;
  //complex <double> LQ2F2 = pow(logq2muf2,2)-4.*pow(loga,2);
  //double LQF = LF-LQ;
  //double LQ2F2 = LF2-LQ2;
  //if (opts.mufevol) //No! togliere le gamma invece! --> should be actually if (opts.evolmode == 2 || opts.evolmode == 3)    
  //{
      //LQF = 0.;
      //LQ2F2 = 0.;
      //LQF = resint::LF;
      //LQ2F2 = pow(resint::LF,2);
  //}
  /*
  double aassh=aass/2.;
  double aasshsq = pow(aassh,2);
  double aass2=pow(aass,2);
      
  complex <double> loga2 = pow(loga,2);
  complex <double> loga3 = pow(loga,3);
  complex <double> logq2muf22 = pow(logq2muf2,2);
  H1stgg=0;
  */
//  if (opts.order == 1)
//    for (int i = 0; i < mellinint::mdim; i++)
//      {
//	int idx = anomalous::index(i,mesq::positive);
//	
//	complex <double> H1st_qqb = C1QQ[idx]
//	  -gamma1qq[idx]*(-LQF)
//	  +(-2.*loga)*(B1q+A1q*loga);
//	
//	complex <double> H1st_qg = C1QG[idx]
//	  -gamma1qg[idx]*(-LQF);
//	
//	Hqqb[i] = 1.+aass*H1st_qqb;
//	Hqg[i] = aass/2.*H1st_qg;
//      }
//
//  if (opts.order >= 2)
//    {
//      for (int i = 0; i < mellinint::mdim; i++)
//	{
//	  int idx = anomalous::index(i,mesq::positive);
//
//	  complex <double> H1st_qqb = C1QQ[idx]
//	    -gamma1qq[idx]*(-LQF)
//	    +(-2.*loga)*(B1q+A1q*loga);
//	  
//	  H1st_qqb *= 2.;
//	    
//	  complex <double> H1st_qg = C1QG[idx]
//	    -gamma1qg[idx]*(-LQF);
//
//
////	  H2stqq[i] = C2NSqqbM[idx] + C2SqqbM[idx]
////	    +4.*(-2.*loga*1./4.*(gamma2qqbV[idx]+gamma2qqbS[idx])
////		 +1./4.*(gamma2qqbV[idx]+gamma2qqbS[idx])*logq2muf2
////		 +1./2.*(H1stqg[i]+C1QG[idx])/2.
////		 *1./2.*gamma1gq[idx]*(LQF));
////
////	  H2stqqp[i] = C2SqqbM[idx]
////	    +4.*(-2.*loga*1./4.*gamma2qqbS[idx]
////		 +1./4.*gamma2qqbS[idx]*logq2muf2
////		 +1./2.*(H1stqg[i]+C1QG[idx])/2.
////		 *1./2.*gamma1gq[idx]*(LQF));
////
////	  H2stqqb[i] = C2NSqqM[idx] + C2NSqqM[idx] + C2SqqbM[idx] + C2SqqbM[idx] + C1QQ[idx]*C1QQ[idx]
////	    + 4.*(+ 1./6.*A1q*beta0*8*loga3  
////		  + 1./2.*4*loga2*(A2q-beta0*(B1q+2*A1q*loga+gamma1qq[idx]/2.+gamma1qq[idx]/2.))
////		  - 2.*loga*(B2q+2*A2q*loga-beta0*(C1QQ[idx]+C1QQ[idx])/2.
////			     + gamma2qqV[idx]/4. + gamma2qqS[idx]/4. 
////			     + gamma2qqV[idx]/4. + gamma2qqS[idx]/4. )      
////		  + beta0/2.*(gamma1qq[idx]+gamma1qq[idx])/2.*logq2muf22
////		  + (gamma2qqV[idx]/4. + gamma2qqS[idx]/4. + gamma2qqV[idx]/4. + gamma2qqS[idx]/4.)*logq2muf2
////		  - H1stqqb[i]/2.*beta0*logq2mur2
////		  + 1./2.*(H1stqqb[i]+C1QQ[idx]+C1QQ[idx])/2.*((gamma1qq[idx]+gamma1qq[idx])/2.*(LQF)-(B1q+A1q*loga)*2.*loga)
////		  + 1./4.*(H1stqg[i]+C1QG[idx])
////		  * gamma1gq[idx]/2.*(LQF)
////		  + 1./4.*(H1stqg[i]+C1QG[idx])
////		  * gamma1gq[idx]/2.*(LQF));
////
////	  H2stqg[i] = C2qgM[idx] + C1QG[idx]*C1QQ[idx] 
////	    + 4.*(+ 1./2.*beta0*4*loga2*(-gamma1qg[idx]/2.)
////		  - 2.*loga*(-beta0 * C1QG[idx]/2. + gamma2qg[idx]/4.)
////		  + 1./2.*beta0*logq2muf22*gamma1qg[idx]/2.
////		  + gamma2qg[idx]/4.*logq2muf2
////		  - beta0*logq2mur2*H1stqg[i]/2.
////		  + 1./2.*(H1stqqb[i] + C1QQ[idx] + C1QQ[idx])/2.*(LQF)*gamma1qg[idx]/2.
////		  + 1./2.*(H1stqg[i] + C1QG[idx])/2.*((LQF)*(gamma1qq[idx]+gamma1gg[idx])/2.-(B1q+A1q*loga)*2.*loga));
////
////	  H2stgg[i] = C1QG[idx]*C1QG[idx]
////	    -4.*1./2.*(logmuf2q2+2.*loga)*((H1stqg[i] + C1QG[idx])/2.*gamma1qg[idx]/2.+(H1stqg[i] + C1QG[idx])/2.*gamma1qg[idx]/2.);
//
//	  
//	    
//	  
//	  //  All *4 because of normalization (as/2pi)^2 instead of as/pi
//	  // normalization of as/2pi  implies  gamma^1->gamma^1/2
//	  //                                   C^1, H^1 -> C^1/2, H^1/2 
//	  //                                   gamma^2-> gamma^2/4
//	  //                              beta,A,B -> beta,A,B are in as/pi already
//
//	  H2stqqb[i] = C2NSqqM[idx] + C2NSqqM[idx] + C2SqqbM[idx] + C2SqqbM[idx] + C1QQ[idx]*C1QQ[idx]
//	    + 4.*(+ 1./6.*A1q*beta0*8*loga3  
//		  + 1./2.*4*loga2*(A2q-beta0*(B1q+2*A1q*loga))
//		  - 2.*loga*(B2q+2*A2q*loga-beta0*(C1QQ[idx]+C1QQ[idx])/2.)
//		  + (gamma2qqV[idx]/4. + gamma2qqS[idx]/4. + gamma2qqV[idx]/4. + gamma2qqS[idx]/4.)*LQF
//		  + beta0/2.*(gamma1qq[idx]+gamma1qq[idx])/2.*LQ2F2
//		  - H1st_qqb/2.*beta0*logq2mur2
//		  + 1./2.*(H1st_qqb+C1QQ[idx]+C1QQ[idx])/2.*((gamma1qq[idx]+gamma1qq[idx])/2.*(LQF)-(B1q+A1q*loga)*2.*loga)
//		  + 1./4.*(H1st_qg+C1QG[idx])*gamma1gq[idx]/2.*(LQF)
//		  + 1./4.*(H1st_qg+C1QG[idx])*gamma1gq[idx]/2.*(LQF)
//		  );
//
//
//	  H2stqg[i] = C2qgM[idx] + C1QG[idx]*C1QQ[idx]
//	    + 4.*(- 2.*loga*(-beta0 * C1QG[idx]/2.)
//		  + 1./2.*beta0*gamma1qg[idx]/2.*LQ2F2
//		  + gamma2qg[idx]/4.*LQF
//		  - beta0*logq2mur2*H1st_qg/2.
//		  + 1./2.*(H1st_qqb+C1QQ[idx]+C1QQ[idx])/2.*(LQF)*gamma1qg[idx]/2.
//		  + 1./2.*(H1st_qg+C1QG[idx])/2.*((LQF)*(gamma1qq[idx]+gamma1gg[idx])/2.-(B1q+A1q*loga)*2.*loga)
//		  );
//
//	  //qq  (Qb Qb -> Qb Q =  Q Q -> Q Qb)
//	  H2stqq[i] = C2NSqqbM[idx] + C2SqqbM[idx]
//	    +4.*(1./4.*(gamma2qqbV[idx]+gamma2qqbS[idx])*LQF
//		 +1./2.*(H1st_qg+C1QG[idx])/2.
//		 *1./2.*gamma1gq[idx]*(LQF));
//
//	  
//	  //qqp (Q Q' -> Q Qb)
//	  H2stqqp[i] = C2SqqbM[idx]
//	    +4.*(1./4.*gamma2qqbS[idx]*LQF
//		 +1./2.*(H1st_qg+C1QG[idx])/2.
//		 *1./2.*gamma1gq[idx]*(LQF));
//	  
//	  
//	  
//	  //gg
//	  H2stgg[i] = C1QG[idx]*C1QG[idx]
//	    +4.*1./2.*(LQF)*((H1st_qg + C1QG[idx])/2.*gamma1qg[idx]/2.+(H1st_qg + C1QG[idx])/2.*gamma1qg[idx]/2.);
//
//
////	  //QQb
////	  Hqqb[i] = 1.+aassh*H1st_qqb+aasshsq*H2stqqb[i];
////	    
////	  //QG and GQ
////	  Hqg[i] = aassh*H1st_qg+aasshsq*H2stqg[i];
////
////	  //QQ
////	  Hqq[i] = aasshsq*H2stqq[i];
////
////	  //QQ'
////	  Hqqp[i] = aasshsq*H2stqqp[i];
////
////	  //QQ'
////	  Hqqbp[i] = aasshsq*H2stqqp[i];
////	  
////	  //GG
////	  Hgg[i] = aasshsq*H2stgg[i];
//
//	}
//    }
}

//b-dependent coefficients --> clean up, this is now done in expc.C
void hcoeff::calcb(double aass, complex <double> logmuf2q2, complex <double> loga, complex <double> alpq, complex <double> aexp, complex <double> aexpb)
{
  if (opts.order == 0)
    return;
  /*
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
	  //Hqg[i] = aexp*(aass/2.)*C1QG[idx]
	  //+(aass/2.)*(-gamma1qg[idx])*(logmuf2q2+2.*loga);

	  //I would write instead // -> No, should be kept as original
	  Hqg[i] = (aassh*H1stqg[i])
	    *aexp;
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
      complex <double> aexpC = cx(exponent_.aexpc_);
      complex <double> aexpD = cx(exponent_.aexpd_);
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

	  //aexpdod[i]=pow(aexpb,1./2.*aass*C3/C2)
	  //  *pow(aexpC,1./2.*aass2*(pow(C3/C2,2)/2. - C4/C2))
	  //  *pow(aexpD,1./2.*aass2*beta1/beta0*C3/C2)
	  //  *pow(aexpC,1./2.*aass2*beta0*C3/C2*(resint::rlogq2mur2-2.*resint::rloga))
	  //  ;
	  
	  aexpqq[i]=pow(aexpb,1./2.*aass*C1QQ[idx])
	    *pow(aexpC,1./2.*aass2*(0.5*pow(C1QQ[idx],2) - (C2NSqqM[idx] + C2SqqbM[idx])))
	    *pow(aexpD,-1./2.*aass2*beta1/beta0*C1QQ[idx])
	    *pow(aexpC,1./2.*aass2*beta0*C1QQ[idx]*(resint::rlogq2mur2-2.*resint::rloga))
	  ;
	  
	  aexpqg[i]=pow(aexpb,1./2.*aass*(C2qgM[idx]/C1QG[idx]))
	    //*pow(aexpC,1./2.*aass2*(pow(C2qgM[idx]/C1QG[idx],2)/2. - C3/C1QG[idx]))
	    *pow(aexpD,-1./2.*aass2*beta1/beta0*C2qgM[idx]/C1QG[idx])
	    *pow(aexpC,1./2.*aass2*(beta0*C2qgM[idx]/C1QG[idx])*(resint::rlogq2mur2-2.*resint::rloga))
	    ;
	  
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
  */
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

