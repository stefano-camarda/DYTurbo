#include "anomalous.h"
#include "mesq.h"
#include "interface.h"
#include "resconst.h"

#include <iostream>

complex <double> *anomalous::ans;
complex <double> *anomalous::am;
complex <double> *anomalous::ap;
complex <double> *anomalous::al;
complex <double> *anomalous::be;
complex <double> *anomalous::ab;
complex <double> *anomalous::rmin;
complex <double> *anomalous::rplus;
complex <double> *anomalous::rqq;
complex <double> *anomalous::rqg;
complex <double> *anomalous::rgq;
complex <double> *anomalous::rgg;

complex <double> *anomalous::RMMQQ;
complex <double> *anomalous::RMMQG;
complex <double> *anomalous::RMMGQ;
complex <double> *anomalous::RMMGG;
complex <double> *anomalous::RMPQQ;
complex <double> *anomalous::RMPQG;
complex <double> *anomalous::RMPGQ;
complex <double> *anomalous::RMPGG;
complex <double> *anomalous::RPMQQ;
complex <double> *anomalous::RPMQG;
complex <double> *anomalous::RPMGQ;
complex <double> *anomalous::RPMGG;
complex <double> *anomalous::RPPQQ;
complex <double> *anomalous::RPPQG;
complex <double> *anomalous::RPPGQ;
complex <double> *anomalous::RPPGG;

complex <double> *anomalous::C1QQ;
complex <double> *anomalous::C1QG;
complex <double> *anomalous::C1GQ;
complex <double> anomalous::C1GG;

complex <double> *anomalous::gamma1qq;
complex <double> *anomalous::gamma1qg;
complex <double> *anomalous::gamma1gq;
complex <double> *anomalous::gamma1gg;
complex <double> *anomalous::gamma2qq;
complex <double> *anomalous::gamma2qqV;
complex <double> *anomalous::gamma2qqbV;
complex <double> *anomalous::gamma2qqS;
complex <double> *anomalous::gamma2qqbS;
complex <double> *anomalous::gamma2qg;
complex <double> *anomalous::gamma2gq;
complex <double> *anomalous::gamma2gg;

complex <double> *anomalous::C2qgM;
complex <double> *anomalous::C2NSqqM;
complex <double> *anomalous::C2SqqbM;
complex <double> *anomalous::C2NSqqbM;

void anomalous::init()
{
  // ANOMALOUS DIMENSIONS FOR LEADING AND NEXT TO LEADING ORDER
  // EVOLUTION OF PARTON DENSITIES AND WILSON COEFFICIENTS FOR
  // NLO STRUCTURE FUNCTIONS :

  // all values are precalculated at the values of the z points used for the complex contour of the Mellin inverse transform

  //settings: nf = 5 (flavours)
  int nf = resconst::NF;

  ans = new complex <double>[mellinint::mdim*2];
  am = new complex <double> [mellinint::mdim*2];
  ap = new complex <double> [mellinint::mdim*2];
  al = new complex <double> [mellinint::mdim*2];
  be = new complex <double> [mellinint::mdim*2];
  ab = new complex <double> [mellinint::mdim*2];
  rmin = new complex <double> [mellinint::mdim*2];
  rplus = new complex <double> [mellinint::mdim*2];
  rqq = new complex <double> [mellinint::mdim*2];
  rqg = new complex <double> [mellinint::mdim*2];
  rgq = new complex <double> [mellinint::mdim*2];
  rgg = new complex <double> [mellinint::mdim*2];

  RMMQQ = new complex <double> [mellinint::mdim*2];
  RMMQG = new complex <double> [mellinint::mdim*2];
  RMMGQ = new complex <double> [mellinint::mdim*2];
  RMMGG = new complex <double> [mellinint::mdim*2];
  RMPQQ = new complex <double> [mellinint::mdim*2];
  RMPQG = new complex <double> [mellinint::mdim*2];
  RMPGQ = new complex <double> [mellinint::mdim*2];
  RMPGG = new complex <double> [mellinint::mdim*2];
  RPMQQ = new complex <double> [mellinint::mdim*2];
  RPMQG = new complex <double> [mellinint::mdim*2];
  RPMGQ = new complex <double> [mellinint::mdim*2];
  RPMGG = new complex <double> [mellinint::mdim*2];
  RPPQQ = new complex <double> [mellinint::mdim*2];
  RPPQG = new complex <double> [mellinint::mdim*2];
  RPPGQ = new complex <double> [mellinint::mdim*2];
  RPPGG = new complex <double> [mellinint::mdim*2];

  C1QQ = new complex <double> [mellinint::mdim*2];
  C1QG = new complex <double> [mellinint::mdim*2];
  C1GQ = new complex <double> [mellinint::mdim*2];

  gamma1qq = new complex <double> [mellinint::mdim*2];
  gamma1qg = new complex <double> [mellinint::mdim*2];
  gamma1gq = new complex <double> [mellinint::mdim*2];
  gamma1gg = new complex <double> [mellinint::mdim*2];
  gamma2qq = new complex <double> [mellinint::mdim*2];
  gamma2qqV = new complex <double> [mellinint::mdim*2];
  gamma2qqbV = new complex <double> [mellinint::mdim*2];
  gamma2qqS = new complex <double> [mellinint::mdim*2];
  gamma2qqbS = new complex <double> [mellinint::mdim*2];
  gamma2qg = new complex <double> [mellinint::mdim*2];
  gamma2gq = new complex <double> [mellinint::mdim*2];
  gamma2gg = new complex <double> [mellinint::mdim*2];
  
  C2qgM = new complex <double> [mellinint::mdim*2];
  C2NSqqM = new complex <double> [mellinint::mdim*2];
  C2SqqbM = new complex <double> [mellinint::mdim*2];
  C2NSqqbM = new complex <double> [mellinint::mdim*2];

  for (int sign = mesq::positive; sign <= mesq::negative; sign++)
    for (int i = 0; i < mellinint::mdim; i++)
      {
	fcomplex fxn; //input of ancalc
	complex <double> cxn;
	if (sign == mesq::positive)
	  {
	    fxn.real = real(mellinint::Np[i]);
	    fxn.imag = imag(mellinint::Np[i]);
	    cxn = mellinint::Np[i];
	  }
	else
	  {
	    fxn.real = real(mellinint::Nm[i]);
	    fxn.imag = imag(mellinint::Nm[i]);
	    cxn = mellinint::Nm[i];
	  }
	fcomplex QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F, QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F;
	//input: fxn
	//output: QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F, QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F;
	ancalc_(QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F, QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, fxn);

	fcomplex ANS, AM, AP, AL, BE, AB, RMIN, RPLUS, RQQ, RQG, RGQ, RGG, C2Q, C2G, CDYQ, CDYG; //output of anom
	fcomplex CDYQI;// = 0.; //unused dummy input
	fcomplex CDYGI;// = 0.; //unused dummy input
	fcomplex C2QI;// = 0.; //unused dummy input --> check
	fcomplex C2GF;// = 0.; //unused dummy input --> check

	//inputs: QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F, QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, C2QI, C2GF, CDYQI, CDYGI (only function of fxn and nf)
	//outputs ans, am, ap, al, be, ab, rmin, rplus, rqq, rgg, rgq, rgg (are also only function of fxn)
	anom_(ANS, AM, AP, AL, BE, AB, RMIN, RPLUS, RQQ, RQG,RGQ, RGG, C2Q, C2G, CDYQ, CDYG,
	      fxn, nf, QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F, QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, C2QI, C2GF, CDYQI, CDYGI);

	ans[index(i,sign)] = cx(ANS);
	am[index(i,sign)] = cx(AM);
	ap[index(i,sign)] = cx(AP); 
	al[index(i,sign)] = cx(AL);
	be[index(i,sign)] = cx(BE);
	ab[index(i,sign)] =  cx(AB);
	rmin[index(i,sign)] = cx(RMIN);
	rplus[index(i,sign)] = cx(RPLUS);
	rqq[index(i,sign)] = cx(RQQ);
	rqg[index(i,sign)] = cx(RQG);
	rgq[index(i,sign)] = cx(RGQ);
	rgg[index(i,sign)] = cx(RGG);

	// Compute R coefficients
	complex <double> AC   = 1.- cx(AL);
	complex <double> NMP  = 1.- cx(AM) + cx(AP);
	complex <double> NPM  = 1.- cx(AP) + cx(AM);
	complex <double> DMQQ =  cx(AL) *  cx(RQQ) + cx(BE) * cx(RGQ);
	complex <double> DMQG =  cx(AL) *  cx(RQG) + cx(BE) * cx(RGG);
	complex <double> DMGQ =  cx(AB) *  cx(RQQ) +    AC  * cx(RGQ);
	complex <double> DMGG =  cx(AB) *  cx(RQG) +    AC  * cx(RGG);
	complex <double> DPQQ =     AC  *  cx(RQQ) - cx(BE) * cx(RGQ);
	complex <double> DPQG =     AC  *  cx(RQG) - cx(BE) * cx(RGG);
	complex <double> DPGQ = -cx(AB) *  cx(RQQ) + cx(AL) * cx(RGQ);
	complex <double> DPGG = -cx(AB) *  cx(RQG) + cx(AL) * cx(RGG);

	RMMQQ[index(i,sign)]=  cx(AL) * DMQQ+cx(AB)*DMQG;
	RMMQG[index(i,sign)]=  cx(BE) * DMQQ+  AC  *DMQG;
	RMMGQ[index(i,sign)]=  cx(AL) * DMGQ+cx(AB)*DMGG;
	RMMGG[index(i,sign)]=  cx(BE) * DMGQ+  AC  *DMGG;
	RMPQQ[index(i,sign)]= (  AC   * DMQQ-cx(AB)*DMQG)/NMP;
	RMPQG[index(i,sign)]=(-cx(BE) * DMQQ+cx(AL)*DMQG)/NMP;
	RMPGQ[index(i,sign)]= (  AC   * DMGQ-cx(AB)*DMGG)/NMP;
	RMPGG[index(i,sign)]=(-cx(BE) * DMGQ+cx(AL)*DMGG)/NMP;
	RPMQQ[index(i,sign)]= (cx(AL) * DPQQ+cx(AB)*DPQG)/NPM;
	RPMQG[index(i,sign)]= (cx(BE) * DPQQ+  AC  *DPQG)/NPM;
	RPMGQ[index(i,sign)]= (cx(AL) * DPGQ+cx(AB)*DPGG)/NPM;
	RPMGG[index(i,sign)]= (cx(BE) * DPGQ+  AC  *DPGG)/NPM;
	RPPQQ[index(i,sign)]=    AC   * DPQQ-cx(AB)*DPQG;
	RPPQG[index(i,sign)]= -cx(BE) * DPQQ+cx(AL)*DPQG;
	RPPGQ[index(i,sign)]=    AC   * DPGQ-cx(AB)*DPGG;
	RPPGG[index(i,sign)]= -cx(BE) * DPGQ+cx(AL)*DPGG;

      
	// C1 coefficients      
	double pi2 = M_PI*M_PI;
	C1QG[index(i,sign)]=1./((cxn+1.)*(cxn+2.));
	C1GQ[index(i,sign)]=4./3./(cxn+1.);
	C1QQ[index(i,sign)]=2.*pi2/3.-16./3.+4./3./(cxn*(cxn+1.)); // = 2. * resconst::C1qqn + 4./3./(cxn*(cxn+1.));
	C1GG=pi2/2.+11./2.+pi2; // = 2. * resconst::C1ggn;

	// gamma1 gamma2: NORMALIZED ANOMALOUS DIMENSIONS AND COEFFICIENTS
	gamma1qq[index(i,sign)]=-1.*(cx(QQI)/4.);
	gamma1qg[index(i,sign)]=-1.*(cx(QGF)/8.);
	gamma1gq[index(i,sign)]=-1.*(cx(GQI)/4.);
	gamma1gg[index(i,sign)]=-1.*((cx(GGI)+(double)nf*cx(GGF))/4.);
	gamma2qq[index(i,sign)]=-1.*(((cx(NS1PI)+(double)nf*cx(NS1F))+(double)nf*cx(QQ1F))/8.);
	gamma2qqV[index(i,sign)]=-1.*(cx(NS1PI)+2*(double)nf*cx(NS1F)+cx(NS1MI))/16.;
	gamma2qqbV[index(i,sign)]=-1.*(cx(NS1PI)-cx(NS1MI))/16.;
	gamma2qqS[index(i,sign)]=-1.*(cx(QQ1F)/16.);
	gamma2qqbS[index(i,sign)]=gamma2qqS[index(i,sign)];
	gamma2qg[index(i,sign)]=-1.*(cx(QG1F)/16.);
	gamma2gq[index(i,sign)]=-1.*((cx(GQ1I)+(double)nf*cx(GQ1F))/8.);
	gamma2gg[index(i,sign)]=-1.*((cx(GG1I)+(double)nf*cx(GG1F))/8.);

	// C2 coefficients      
	//input: fxn and nf from the common block
	fcomplex C2qg,C2NSqqb,C2NSqq,C2Sqqb;
	h2calc_(C2qg,C2NSqqb,C2NSqq,C2Sqqb,fxn);
	C2qgM[index(i,sign)]    = cx(C2qg);
	C2NSqqM[index(i,sign)]  = cx(C2NSqq);
	C2SqqbM[index(i,sign)]  = cx(C2Sqqb);
	C2NSqqbM[index(i,sign)] = cx(C2NSqqb);
      }
}
