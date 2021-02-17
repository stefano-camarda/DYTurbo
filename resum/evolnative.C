#include "evolnative.h"

#include "pdfevol.h"
#include "mellinint.h"
#include "mellinpdf.h"
#include "settings.h"
#include "scales.h"
#include "mesq.h"
#include "resconst.h"
#include "anomalous.h"
#include "resint.h"
#include "phasespace.h"
#include "psi.h"
#include "constants.h"
#include "melfun.h"
//#include "clock_real.h"

#include <iostream>

int evolnative::ndim;

//PDFs mellin moments at the factorisation scale
complex <double> *evolnative::UVP;
complex <double> *evolnative::DVP;
complex <double> *evolnative::SVP;
complex <double> *evolnative::CVP;
complex <double> *evolnative::BVP;
complex <double> *evolnative::USP;
complex <double> *evolnative::DSP;
complex <double> *evolnative::SSP;
complex <double> *evolnative::CHP;
complex <double> *evolnative::BOP;
complex <double> *evolnative::GLP;

complex <double> *evolnative::SIP;
complex <double> *evolnative::NS3P;
complex <double> *evolnative::NS8P;
complex <double> *evolnative::NS15P;
complex <double> *evolnative::NS24P;
complex <double> *evolnative::NS35P;

complex <double> evolnative::XL;
complex <double> evolnative::XL1;
complex <double> evolnative::SALP;

complex <double> evolnative::alpr;
complex <double> evolnative::alpq;
complex <double> evolnative::alpqf;

complex <double> *evolnative::ans;
complex <double> *evolnative::am;
complex <double> *evolnative::ap;
complex <double> *evolnative::al;
complex <double> *evolnative::be;
complex <double> *evolnative::ab;
complex <double> *evolnative::rmin;
complex <double> *evolnative::rplus;
complex <double> *evolnative::rqq;
complex <double> *evolnative::rqg;
complex <double> *evolnative::rgq;
complex <double> *evolnative::rgg;

complex <double> *evolnative::rmmqq;
complex <double> *evolnative::rmmqg;
complex <double> *evolnative::rmmgq;
complex <double> *evolnative::rmmgg;
complex <double> *evolnative::rmpqq;
complex <double> *evolnative::rmpqg;
complex <double> *evolnative::rmpgq;
complex <double> *evolnative::rmpgg;
complex <double> *evolnative::rpmqq;
complex <double> *evolnative::rpmqg;
complex <double> *evolnative::rpmgq;
complex <double> *evolnative::rpmgg;
complex <double> *evolnative::rppqq;
complex <double> *evolnative::rppqg;
complex <double> *evolnative::rppgq;
complex <double> *evolnative::rppgg;

//separate allocation and update of PDF moments at the starting scale and evolution engine

//allocate PDFs at the starting scale
void evolnative::allocate_pdfs()
{
  UVP = new complex <double>[ndim];
  DVP = new complex <double>[ndim];
  SVP = new complex <double>[ndim];
  CVP = new complex <double>[ndim];
  BVP = new complex <double>[ndim];
  USP = new complex <double>[ndim];
  DSP = new complex <double>[ndim];
  SSP = new complex <double>[ndim];
  GLP = new complex <double>[ndim];
  CHP = new complex <double>[ndim];
  BOP = new complex <double>[ndim];
  
  SIP   = new complex <double>[ndim];
  NS3P  = new complex <double>[ndim];
  NS8P  = new complex <double>[ndim];
  NS15P = new complex <double>[ndim];
  NS24P = new complex <double>[ndim];
  NS35P = new complex <double>[ndim];
}

//allocate evolution engine
void evolnative::allocate_engine()
{
  ans   = new complex <double> [ndim*2];
  am    = new complex <double> [ndim*2];
  ap    = new complex <double> [ndim*2];
  al    = new complex <double> [ndim*2];
  be    = new complex <double> [ndim*2];
  ab    = new complex <double> [ndim*2];
  rmin  = new complex <double> [ndim*2];
  rplus = new complex <double> [ndim*2];
  rqq   = new complex <double> [ndim*2];
  rqg   = new complex <double> [ndim*2];
  rgq   = new complex <double> [ndim*2];
  rgg   = new complex <double> [ndim*2];

  rmmqq = new complex <double> [ndim*2];
  rmmqg = new complex <double> [ndim*2];
  rmmgq = new complex <double> [ndim*2];
  rmmgg = new complex <double> [ndim*2];
  rmpqq = new complex <double> [ndim*2];
  rmpqg = new complex <double> [ndim*2];
  rmpgq = new complex <double> [ndim*2];
  rmpgg = new complex <double> [ndim*2];
  rpmqq = new complex <double> [ndim*2];
  rpmqg = new complex <double> [ndim*2];
  rpmgq = new complex <double> [ndim*2];
  rpmgg = new complex <double> [ndim*2];
  rppqq = new complex <double> [ndim*2];
  rppqg = new complex <double> [ndim*2];
  rppgq = new complex <double> [ndim*2];
  rppgg = new complex <double> [ndim*2];
}
  
void evolnative::allocate()
{
  allocate_pdfs();
  if (opts.melup == 2)
    allocate_engine();
}
void evolnative::free_pdfs()
{
  delete[] UVP;
  delete[] DVP;
  delete[] SVP;
  delete[] CVP;
  delete[] BVP;
  delete[] USP;
  delete[] DSP;
  delete[] SSP;
  delete[] GLP;
  delete[] CHP;
  delete[] BOP;
  
  delete[] SIP;
  delete[] NS3P;
  delete[] NS8P;
  delete[] NS15P;
  delete[] NS24P;
  delete[] NS35P;
}
void evolnative::free_engine()
{
  delete[] ans;
  delete[] am;
  delete[] ap;
  delete[] al;
  delete[] be;
  delete[] ab;
  delete[] rmin;
  delete[] rplus;
  delete[] rqq;
  delete[] rqg;
  delete[] rgq;
  delete[] rgg;

  delete[] rmmqq;
  delete[] rmmqg;
  delete[] rmmgq;
  delete[] rmmgg;
  delete[] rmpqq;
  delete[] rmpqg;
  delete[] rmpgq;
  delete[] rmpgg;
  delete[] rpmqq;
  delete[] rpmqg;
  delete[] rpmgq;
  delete[] rpmgg;
  delete[] rppqq;
  delete[] rppqg;
  delete[] rppgq;
  delete[] rppgg;
}

//free PDFs at the starting scale
void evolnative::free()
{
  free_pdfs();
  if (opts.melup == 2)
    free_engine();
}

void evolnative::init()
{
  //if (opts.evolmode != 0) //--> Always init for ctmellin
  //return;

  if (opts.mellin1d)
    ndim = mellinint::mdim;
  else
    ndim = 2*mellinint::mdim;
  
  //  if (opts.fmufac != 0) //probably better to drop distinction between fixed and mll scale... ? It creates a lot of trouble, thread private needs copyin, etc...
  //    return;

  /*
  //calculate Mellin moments of PDFs
  if (!opts.silent) cout << "Initialise PDF moments with numerical integration... " << flush;
  clock_t begin_time, end_time;
  begin_time = clock();  

  allocate_pdfs();
  scales::set(opts.rmass);
  update_pdfs();
  */

  if (opts.melup <= 1)
    allocate_engine();
  if (opts.melup == 0)
    update_engine();
    
  //end_time = clock();  
  //if (!opts.silent) cout << "Done in "  << float(end_time - begin_time)/CLOCKS_PER_SEC*1000. << "ms" << endl;
}

void evolnative::release()
{
  //if (opts.evolmode != 0)  //--> should always free evolnative for the counterterm in Mellin space?
  //return;

  //  if (opts.fmufac != 0)
  //    return;
  
  //free(); //Do not free at the end of DYTurbo as all variables are local thread
  if (opts.melup <= 1)
    free_engine();
}

//reimplementation of ancalc and anom
void anom(
	  complex <double> XN,
	  complex <double> &ans,
	  complex <double> &am,
	  complex <double> &ap,
	  complex <double> &al,
	  complex <double> &be,
	  complex <double> &ab,
	  complex <double> &rmin,
	  complex <double> &rplus,
	  complex <double> &rqq,
	  complex <double> &rqg,
	  complex <double> &rgq,
	  complex <double> &rgg
	  )
{
  //Calculation of anomalous dimensions and wilson coefficients
  //up to their dependence of the number of active flavours f :
  complex <double> XNS = XN * XN;
  complex <double> XN1 = XN + 1.;
  complex <double> XN2 = XN + 2.;
  complex <double> XNM = XN - 1.;
  //Leading order
  complex <double> CPSI = cpsi0(XN1) + constants::euler;
  complex <double> QQI = (8./3.) * (-3.- 2./(XN * XN1) + 4.* CPSI);
  complex <double> QGF = -4.* (XNS + XN +2.) / (XN * XN1 * XN2);
  complex <double> GQI = -(16./3.) * (XNS + XN + 2.) / (XN * XN1 * XNM);
  complex <double> GGI = -22.- 24./(XN * XNM) - 24./(XN1 * XN2) + 24.* CPSI;
  complex <double> GGF = 4./3.;
  //Next to leading order
  complex <double> XNT = XNS * XN;
  complex <double> XNFO = XNT * XN;
  complex <double> XN1S = XN1 * XN1;
  complex <double> XN1T = XN1S * XN1;
  //Analytic continuations of n-sums as given in glueck et al. (1990)
  double ZETA2 = constants::zeta2;
  double ZETA3 = constants::zeta3;
  complex <double> CPSI1 = ZETA2 - cpsi(1,XN1);
//       SPMOM = 1.01/XN1 - 0.846/XN2 + 1.155/(XN+3.) - 1.074/(XN+4.) +
//     1         0.55/(XN+5.)
//       SPMOM = 1.004D0 / XN1 - 0.846D0 / XN2 + 1.342D0 / (XN+3.) -
//     1         1.532D0 / (XN+4.) + 0.839D0 / (XN+5.)
  complex <double> SPMOM = fun6(XN-1.); // --> Improved approximation for the Mellin transform of: Li2(x)/(1+x)

  complex <double> SLC = -5./8.* ZETA3;
  complex <double> SLV = - ZETA2/2.* (cpsi0 (XN1/2.) - cpsi0 (XN/2.))
    + CPSI/XNS + SPMOM;
  complex <double> SSCHLM = SLC - SLV;
  complex <double> SSTR2M = ZETA2 - cpsi(1,XN1/2.);
  complex <double> SSTR3M = 0.5 * cpsi(2,XN1/2.) + ZETA3;
  complex <double> SSCHLP = SLC + SLV;
  complex <double> SSTR2P = ZETA2 - cpsi(1,XN2/2.);
  complex <double> SSTR3P = 0.5 * cpsi(2,XN2/2.) + ZETA3;
  //Non-singlet pieces as given in curci et al. (1980)
  complex <double> NS1MA = 16.* CPSI * (2.* XN + 1.) / (XNS * XN1S) +
    16.* (2.* CPSI - 1./(XN * XN1)) * ( CPSI1 - SSTR2M ) +
    64.* SSCHLM + 24.* CPSI1 - 3. - 8.* SSTR3M -
    8.* (3.* XNT + XNS -1.) / (XNT * XN1T) +
    16.* (2.* XNS + 2.* XN +1.) / (XNT * XN1T);
  complex <double> NS1PA = 16.* CPSI * (2.* XN + 1.) / (XNS * XN1S) +
    16.* (2.* CPSI - 1./(XN * XN1)) * ( CPSI1 - SSTR2P ) +
    64.* SSCHLP + 24.* CPSI1 - 3. - 8.* SSTR3P -
    8.* (3.* XNT + XNS -1.) / (XNT * XN1T) -
    16.* (2.* XNS + 2.* XN +1.) / (XNT * XN1T);
  complex <double> NS1B = CPSI * (536./9. + 8.* (2.* XN + 1.) / (XNS * XN1S)) -
    (16.* CPSI + 52./3.- 8./(XN * XN1)) * CPSI1 - 43./6. -
    (151.* XNFO + 263.* XNT + 97.* XNS + 3.* XN + 9.) *
    4./ (9.* XNT * XN1T);
  complex <double> NS1C = -160./9.* CPSI + 32./3.* CPSI1 + 4./3. +
    16.* (11.* XNS + 5.* XN - 3.) / (9.* XNS * XN1S);
  complex <double> NS1MI = -2./9.* NS1MA + 4.* NS1B;
  complex <double> NS1PI = -2./9.* NS1PA + 4.* NS1B;
  complex <double> NS1F = 2./3. * NS1C;
  //Singlet pieces as given in floratos et al. (1981)
  complex <double> XNFI = XNFO * XN;
  complex <double> XNSI = XNFI * XN;
  complex <double> XNSE = XNSI * XN;
  complex <double> XNE = XNSE * XN;
  complex <double> XNN = XNE * XN;
  complex <double> XNMS = XNM * XNM;
  complex <double> XN2S = XN2 * XN2;
  complex <double> XN2T = XN2S * XN2;
  complex <double> QQ1F = (5.* XNFI + 32.* XNFO + 49.* XNT + 38.* XNS + 28.* XN
			   + 8.) / (XNM * XNT * XN1T * XN2S) * (-32./3.);
  complex <double> QG1A = (-2.* CPSI * CPSI + 2.* CPSI1 - 2.* SSTR2P)
    * (XNS + XN + 2.) / (XN * XN1 * XN2)
    + (8.* CPSI * (2.* XN + 3.)) / (XN1S * XN2S)
    + 2.* (XNN + 6.* XNE + 15. * XNSE + 25.* XNSI + 36.* XNFI
	   + 85.* XNFO + 128.* XNT + 104.* XNS + 64.* XN + 16.)
    / (XNM * XNT * XN1T * XN2T);
  complex <double> QG1B = (2.* CPSI * CPSI - 2.* CPSI1 + 5.) * (XNS + XN + 2.)
    / (XN * XN1 * XN2)   -   4.* CPSI / XNS
    + (11.* XNFO + 26.* XNT + 15.* XNS + 8.* XN + 4.)
    / (XNT * XN1T * XN2);
  complex <double> QG1F = - 12.* QG1A - 16./3.* QG1B;
  complex <double> GQ1A = (-2.* CPSI * CPSI + 10.* CPSI - 2.* CPSI1)
    * (XNS + XN + 2.) / (XNM * XN * XN1)  -  4.* CPSI / XN1S
    - (12.* XNSI + 30.* XNFI + 43.* XNFO + 28.* XNT - XNS
       - 12.* XN - 4.) / (XNM * XNT * XN1T);
  complex <double> GQ1B = (CPSI * CPSI + CPSI1 - SSTR2P) * (XNS + XN + 2.)
    / (XNM * XN * XN1)
    - CPSI * (17.* XNFO + 41.* XNS - 22.* XN - 12.)
    / (3.* XNMS * XNS * XN1)
    + (109.* XNN + 621.* XNE + 1400.* XNSE + 1678.* XNSI
       + 695.* XNFI - 1031.* XNFO - 1304.* XNT - 152.* XNS
       + 432.* XN + 144.) / (9.* XNMS * XNT * XN1T * XN2S);
  complex <double> GQ1C = (CPSI - 8./3.) * (XNS + XN + 2.) / (XNM * XN * XN1)
    + 1./ XN1S;
  complex <double> GQ1I = - 64./9.* GQ1A - 32.* GQ1B;
  complex <double> GQ1F = - 64./9.* GQ1C;
  complex <double> GG1A = 16./9.* (38.* XNFO + 76.* XNT + 94.* XNS + 56.* XN + 12.)
    / (XNM * XNS * XN1S * XN2)   -   160./9.* CPSI + 32./3.;
  complex <double> GG1B = (2.* XNSI + 4.* XNFI + XNFO - 10.* XNT - 5.* XNS - 4.* XN
			   - 4.) * 16. / (XNM * XNT * XN1T * XN2)   +   8.;
  complex <double> GG1C = (2.* XNFI + 5.* XNFO + 8.* XNT + 7.* XNS - 2.* XN - 2.)
    * 64.* CPSI / (XNMS * XNS * XN1S * XN2S)
    + 536./9.* CPSI - 64./3.
    + 32.* SSTR2P * (XNS + XN + 1.) / (XNM * XN * XN1 * XN2)
    - 16.* CPSI * SSTR2P + 32.* SSCHLP - 4.* SSTR3P
    - 4.* (457.* XNN + 2742.* XNE + 6040.* XNSE + 6098.* XNSI
	   + 1567.* XNFI - 2344.* XNFO - 1632.* XNT + 560.* XNS
	   + 1488.* XN + 576.) / (9.* XNMS * XNT * XN1T * XN2T);
  complex <double> GG1I = 9.* GG1C;
  complex <double> GG1F = 3./2.* GG1A + 2./3.* GG1B;

  //Anomalous dimensions for leading and next to leading order
  //evolution of parton densities and wilson coefficients for
  //NLO structure functions

  //Anomalous dimensions and related quantities in leading order
  int nf = resconst::NF;
  int N = 1;
  double B0 = 11.- 2./3.* double(nf);
  double B02 = 2.* B0;
  complex <double> QQ = QQI;
  complex <double> QG = double(nf) * QGF;
  complex <double> GQ = GQI;
  complex <double> GG = GGI + double(nf) * GGF;
  complex <double> SQ = sqrt((GG - QQ) * (GG - QQ) + 4.* QG * GQ);
  complex <double> GP = 0.5 * (QQ + GG + SQ);
  complex <double> GM = 0.5 * (QQ + GG - SQ);
  ans = QQ / B02;
  am = GM / B02;
  ap = GP / B02;
  al = (QQ - GP) / (GM - GP);
  be = QG / (GM - GP);
  ab = GQ / (GM - GP);
  //Next to leading order : anomalous dimensions and wilson coefficients
  //in the ms-bar factorization scheme of bardeen et al. (1981)
  complex <double> NS1M = NS1MI + double(nf) * NS1F;
  complex <double> NS1P = NS1PI + double(nf) * NS1F;
  complex <double> QQ1 = NS1P + double(nf) * QQ1F;
  complex <double> QG1 = double(nf) * QG1F;
  complex <double> GQ1 = GQ1I + double(nf) * GQ1F;
  complex <double> GG1 = GG1I + double(nf) * GG1F;
  //Combinations of anomalous dimensions for NLO - singlet evolution
  double B1 = 102 - 38./3.* double(nf);
  double B10 = B1 / B0;
  rmin = (NS1M - QQ * B10) / B02;
  rplus = (NS1P - QQ * B10) / B02;
  rqq = (QQ1 - QQ * B10) / B02;
  rqg = (QG1 - QG * B10) / B02;
  rgq = (GQ1 - GQ * B10) / B02;
  rgg = (GG1 - GG * B10) / B02;
}

//Update the PDFs at the factorisation scale by performing the Mellin transform from x to N
void evolnative::update_pdfs()
{
  clock_t begin_time, end_time;

  //Restrict the integration of moments to xmin = m/sqrt(s)*exp(-ymax) = (m/sqrt(s))^2
  //double xmin = pow(phasepace::m/opts.sroot,2);

  //scales::set(phasespace::m); //assume scales were already set
  mellinpdf::allocate();
  //mellinpdf::evalpdfs(scales::fac);
  mellinpdf::evalpdfs(scales::fac,phasespace::m,phasespace::y);
  mellinpdf::update_mellin();
  mellinpdf::transform();

  //Here it would be better to use a common PDF basis between evolnative and mellinpdf (and all other PDF modules in N space)
  if (opts.mellin1d)
    for (int k = 0; k < mellinint::mdim; k++)
      {
	UVP[k] = mellinpdf::UP[k] - mellinpdf::UB[k];
	DVP[k] = mellinpdf::DO[k] - mellinpdf::DB[k];
	SVP[k] = mellinpdf::ST[k] - mellinpdf::SB[k];
	CVP[k] = mellinpdf::CH[k] - mellinpdf::CB[k];
	BVP[k] = mellinpdf::BO[k] - mellinpdf::BB[k];
	USP[k] = mellinpdf::UB[k];
	DSP[k] = mellinpdf::DB[k];
	SSP[k] = mellinpdf::SB[k];
	CHP[k] = mellinpdf::CB[k];
	BOP[k] = mellinpdf::BB[k];
	GLP[k] = mellinpdf::GL[k];
	//cout << "moment " << k << "  " << mellinint::Np[k] << "  ";
	//cout << "uval  " << mellinpdf::UV[k] << endl;
	//cout << "dval  " << mellinpdf::DV[k] << endl;
	//cout << "usea  " << mellinpdf::US[k] << endl;
	//cout << "dsea  " << mellinpdf::DS[k] << endl;
	//cout << "ssea  " << mellinpdf::SM[k] << endl;
	//cout << "gluon " << mellinpdf::GL[k] << endl;
	//cout << "charm " << mellinpdf::CP[k] << endl;
	//cout << "bottom" << mellinpdf::BP[k] << endl;
    }
  else
    for (int k = 0; k < mellinint::mdim; k++)
      {
	UVP[k] = mellinpdf::UP_1[k] - mellinpdf::UB_1[k];
	DVP[k] = mellinpdf::DO_1[k] - mellinpdf::DB_1[k];
	SVP[k] = mellinpdf::ST_1[k] - mellinpdf::SB_1[k];
	CVP[k] = mellinpdf::CH_1[k] - mellinpdf::CB_1[k];
	BVP[k] = mellinpdf::BO_1[k] - mellinpdf::BB_1[k];
	USP[k] = mellinpdf::UB_1[k];
	DSP[k] = mellinpdf::DB_1[k];
	SSP[k] = mellinpdf::SB_1[k];
	CHP[k] = mellinpdf::CB_1[k];
	BOP[k] = mellinpdf::BB_1[k];
	GLP[k] = mellinpdf::GL_1[k];

	UVP[k+mellinint::mdim] = mellinpdf::UP_2[k] - mellinpdf::UB_2[k];
	DVP[k+mellinint::mdim] = mellinpdf::DO_2[k] - mellinpdf::DB_2[k];
	SVP[k+mellinint::mdim] = mellinpdf::ST_2[k] - mellinpdf::SB_2[k];
	CVP[k+mellinint::mdim] = mellinpdf::CH_2[k] - mellinpdf::CB_2[k];
	BVP[k+mellinint::mdim] = mellinpdf::BO_2[k] - mellinpdf::BB_2[k];
	USP[k+mellinint::mdim] = mellinpdf::UB_2[k];
	DSP[k+mellinint::mdim] = mellinpdf::DB_2[k];
	SSP[k+mellinint::mdim] = mellinpdf::SB_2[k];
	CHP[k+mellinint::mdim] = mellinpdf::CB_2[k];
	BOP[k+mellinint::mdim] = mellinpdf::BB_2[k];
	GLP[k+mellinint::mdim] = mellinpdf::GL_2[k];

	//cout << "moment " << k << "  " << mellinint::Np_1[k] << "  ";
	//cout << "uval  " << mellinpdf::UV_1[k] << endl;
	//cout << "dval  " << mellinpdf::DV_1[k] << endl;
	//cout << "usea  " << mellinpdf::US_1[k] << endl;
	//cout << "dsea  " << mellinpdf::DS_1[k] << endl;
	//cout << "ssea  " << mellinpdf::SM_1[k] << endl;
	//cout << "gluon " << mellinpdf::GL_1[k] << endl;
	//cout << "charm " << mellinpdf::CP_1[k] << endl;
	//cout << "bottom" << mellinpdf::BP_1[k] << endl;
	
      }
    

  mellinpdf::free();
  
  if (opts.order_evol == 0)
    return;
  //Precompute singlet/non-singlet decomposition
  
  //N flavour dependence
  int nf = resconst::NF;
  for (int k = 0; k < ndim; k++)
    {
      /*
      NS3P[k] = UVP[k] + 2.*USP[k] - DVP[k] - 2.*DSP[k];               //(u+ub-d-db)               u-d
      NS8P[k] = UVP[k] + 2.*USP[k] + DVP[k] + 2.*DSP[k] - 4.*SSP[k];   //(u+ub+d+db-2(s+sb)        u+d-2s

      if (nf == 5)
	{
	  SIP[k]   = UVP[k] + DVP[k] + 2.*USP[k] + 2.*DSP[k] + 2.*SSP[k] + 2.*CHP[k] + 2.*BOP[k]; //(u+ub+d+db+s+sb+c+cb+b+bb) -> all quarks      u+d+s+c+b
	  NS15P[k] = UVP[k] + DVP[k] + 2.*USP[k] + 2.*DSP[k] + 2.*SSP[k] - 6.*CHP[k];             //(u+ub+d+db+s+sb-3(c+cb))                      u+d+s-3c
	  NS24P[k] = UVP[k] + DVP[k] + 2.*USP[k] + 2.*DSP[k] + 2.*SSP[k] + 2.*CHP[k] - 8.*BOP[k]; //(u+ub+d+db+s+sb+c+cb-4(b+bb))                 u+d+s+c-4b
	  NS35P[k] = SIP[k];
	}
      if (nf == 4)
	{
	  SIP[k]   = UVP[k] + DVP[k] + 2.*USP[k] + 2.*DSP[k] + 2.*SSP[k] + 2.*CHP[k];          //(u+ub+d+db+s+sb+c+cb) -> all quarks
	  NS15P[k] = UVP[k] + DVP[k] + 2.*USP[k] + 2.*DSP[k] + 2.*SSP[k] - 6.*CHP[k];
	  NS24P[k] = UVP[k] + DVP[k] + 2.*USP[k] + 2.*DSP[k] + 2.*SSP[k] + 2.*CHP[k];
	  NS35P[k] = SIP[k];
	}
      if (nf == 3)
	{
	  SIP[k]   = UVP[k] + DVP[k] + 2.*USP[k] + 2.*DSP[k] + 2.*SSP[k];     //(u+ub+d+db+s+sb) -> all quarks
	  NS15P[k] = SIP[k];
	  NS24P[k] = SIP[k];
	  NS35P[k] = SIP[k];
	}
      */

      //--> Version which acounts for s,c,b asymmetries
      NS3P[k] = UVP[k] + 2.*USP[k] - DVP[k] - 2.*DSP[k];                             //(u+ub-d-db)               u-d
      NS8P[k] = UVP[k] + 2.*USP[k] + DVP[k] + 2.*DSP[k] - 2.*(SVP[k] + 2.*SSP[k]);   //(u+ub+d+db-2s-2sb)        u+d-2s
      if (nf == 5)
	{
	  SIP[k]   = UVP[k] + DVP[k] + SVP[k] + CVP[k] + BVP[k] + 2.*USP[k] + 2.*DSP[k] + 2.*SSP[k] + 2.*CHP[k] + 2.*BOP[k];      //(u+ub+d+db+s+sb+c+cb+b+bb) -> all quarks      u+d+s+c+b
	  NS15P[k] = UVP[k] + DVP[k] + SVP[k] + 2.*USP[k] + 2.*DSP[k] + 2.*SSP[k] - 3.*(CVP[k] + 2.*CHP[k]);                      //(u+ub+d+db+s+sb-3(c+cb))                      u+d+s-3c
	  NS24P[k] = UVP[k] + DVP[k] + SVP[k] + CVP[k] + 2.*USP[k] + 2.*DSP[k] + 2.*SSP[k] + 2.*CHP[k] - 4.*(BVP[k] + 2.*BOP[k]); //(u+ub+d+db+s+sb+c+cb-4(b+bb))                 u+d+s+c-4b
	  NS35P[k] = SIP[k];
	}
      if (nf == 4)
	{
	  SIP[k]   = UVP[k] + DVP[k] + SVP[k] + CVP[k] + 2.*USP[k] + 2.*DSP[k] + 2.*SSP[k] + 2.*CHP[k];     //(u+ub+d+db+s+sb+c+cb) -> all quarks
	  NS15P[k] = UVP[k] + DVP[k] + SVP[k] + 2.*USP[k] + 2.*DSP[k] + 2.*SSP[k] - 3.*(CVP[k] + 2.*CHP[k]); //(u+ub+d+db+s+sb-3(c+cb))                      u+d+s-3c
	  NS24P[k] = SIP[k];
	  NS35P[k] = SIP[k];
	}
      if (nf == 3)
	{
	  SIP[k]   = UVP[k] + DVP[k] + SVP[k] + 2.*USP[k] + 2.*DSP[k] + 2.*SSP[k];     //(u+ub+d+db+s+sb) -> all quarks
	  NS15P[k] = SIP[k];
	  NS24P[k] = SIP[k];
	  NS35P[k] = SIP[k];
	}
      
    }
}
void evolnative::update_engine()
{
  if (opts.order_evol == 0)
    return;

  int nf = resconst::NF;
  fcomplex QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F, QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F;
  fcomplex ANS, AM, AP, AL, BE, AB, RMIN, RPLUS, RQQ, RQG, RGQ, RGG, C2Q, C2G, CDYQ, CDYG; //output of anom
  fcomplex CDYQI;// = 0.; //unused dummy input
  fcomplex CDYGI;// = 0.; //unused dummy input
  fcomplex C2QI;// = 0.; //unused dummy input --> check
  fcomplex C2GF;// = 0.; //unused dummy input --> check
  complex <double> AC,NMP,NPM,DMQQ,DMQG,DMGQ,DMGG,DPQQ,DPQG,DPGQ,DPGG ;

  complex <double> cxn[ndim*2] = {0.};
  if (opts.mellin1d)
    for (int sign = mesq::positive; sign <= mesq::negative; sign++)
      for (int i = 0; i < mellinint::mdim; i++)
	{
	  int idx = index(i,sign);
	  if (sign == mesq::positive)
	    cxn[idx] = mellinint::Np[i];
	  else
	    cxn[idx] = mellinint::Nm[i];
	  //cout << idx << "  " << cxn[idx] << endl;
	}
  else
    for (int sign = mesq::positive; sign <= mesq::negative; sign++)
      for (int beam = 1; beam <= 2; beam++)
	for (int i = 0; i < mellinint::mdim; i++)
	  {
	    int idx = index(i,beam,sign);
	    if (beam == 1)
	      if (sign == mesq::positive)
		cxn[idx] = mellinint::Np_1[i];
	      else
		cxn[idx] = mellinint::Nm_1[i];
	    else
	      if (sign == mesq::positive)
		cxn[idx] = mellinint::Np_2[i];
	      else
		cxn[idx] = mellinint::Nm_2[i];
	  }

  for (int i = 0; i < ndim*2; i++)
    {
      fcomplex fxn = fcx(cxn[i]);
      
      //input: fxn
      //output: QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F, QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F;
      ancalc_(QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F, QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, fxn);

      //inputs: QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F, QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, C2QI, C2GF, CDYQI, CDYGI (only function of fxn and nf)
      //outputs ans, am, ap, al, be, ab, rmin, rplus, rqq, rgg, rgq, rgg (are also only function of fxn)
      anom_(ANS, AM, AP, AL, BE, AB, RMIN, RPLUS, RQQ, RQG,RGQ, RGG, C2Q, C2G, CDYQ, CDYG,
	    fxn, nf, QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F, QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, C2QI, C2GF, CDYQI, CDYGI);

      complex <double> Ans, Am, Ap, Al, Be, Ab, Rmin, Rplus, Rqq, Rqg, Rgq, Rgg; //output of anom
      anom(cxn[i], Ans, Am, Ap, Al, Be, Ab, Rmin, Rplus, Rqq, Rqg, Rgq, Rgg);

      //cout << endl;
      //cout << " N " << cxn[i] << endl;
      //
      //cout << "  " << Ans   - cx(ANS) << endl;
      //cout << "  " << Am    - cx(AM) << endl;
      //cout << "  " << Ap    - cx(AP) << endl; 
      //cout << "  " << Al    - cx(AL) << endl;
      //cout << "  " << Be    - cx(BE) << endl;
      //cout << "  " << Ab    - cx(AB) << endl;
      //cout << "  " << Rmin  - cx(RMIN) << endl;
      //cout << "  " << Rplus - cx(RPLUS) << endl;
      //cout << "  " << Rqq   - cx(RQQ) << endl;
      //cout << "  " << Rqg   - cx(RQG) << endl;
      //cout << "  " << Rgq   - cx(RGQ) << endl;
      //cout << "  " << Rgg   - cx(RGG) << endl;
      
      ans[i] = cx(ANS);
      am[i] = cx(AM);
      ap[i] = cx(AP); 
      al[i] = cx(AL);
      be[i] = cx(BE);
      ab[i] =  cx(AB);
      rmin[i] = cx(RMIN);
      rplus[i] = cx(RPLUS);
      rqq[i] = cx(RQQ);
      rqg[i] = cx(RQG);
      rgq[i] = cx(RGQ);
      rgg[i] = cx(RGG);
      
      // Compute R coefficients
      AC   = 1.- cx(AL);
      NMP  = 1.- cx(AM) + cx(AP);
      NPM  = 1.- cx(AP) + cx(AM);
      DMQQ =  cx(AL) *  cx(RQQ) + cx(BE) * cx(RGQ);
      DMQG =  cx(AL) *  cx(RQG) + cx(BE) * cx(RGG);
      DMGQ =  cx(AB) *  cx(RQQ) +    AC  * cx(RGQ);
      DMGG =  cx(AB) *  cx(RQG) +    AC  * cx(RGG);
      DPQQ =     AC  *  cx(RQQ) - cx(BE) * cx(RGQ);
      DPQG =     AC  *  cx(RQG) - cx(BE) * cx(RGG);
      DPGQ = -cx(AB) *  cx(RQQ) + cx(AL) * cx(RGQ);
      DPGG = -cx(AB) *  cx(RQG) + cx(AL) * cx(RGG);

      rmmqq[i]=  cx(AL) * DMQQ+cx(AB)*DMQG;
      rmmqg[i]=  cx(BE) * DMQQ+  AC  *DMQG;
      rmmgq[i]=  cx(AL) * DMGQ+cx(AB)*DMGG;
      rmmgg[i]=  cx(BE) * DMGQ+  AC  *DMGG;
      rmpqq[i]= (  AC   * DMQQ-cx(AB)*DMQG)/NMP;
      rmpqg[i]=(-cx(BE) * DMQQ+cx(AL)*DMQG)/NMP;
      rmpgq[i]= (  AC   * DMGQ-cx(AB)*DMGG)/NMP;
      rmpgg[i]=(-cx(BE) * DMGQ+cx(AL)*DMGG)/NMP;
      rpmqq[i]= (cx(AL) * DPQQ+cx(AB)*DPQG)/NPM;
      rpmqg[i]= (cx(BE) * DPQQ+  AC  *DPQG)/NPM;
      rpmgq[i]= (cx(AL) * DPGQ+cx(AB)*DPGG)/NPM;
      rpmgg[i]= (cx(BE) * DPGQ+  AC  *DPGG)/NPM;
      rppqq[i]=    AC   * DPQQ-cx(AB)*DPQG;
      rppqg[i]= -cx(BE) * DPQQ+cx(AL)*DPQG;
      rppgq[i]=    AC   * DPGQ-cx(AB)*DPGG;
      rppgg[i]= -cx(BE) * DPGQ+cx(AL)*DPGG;

//      if (opts.mellin1d)
//	continue;
//
//      fxn.real = real(mellinint::Np_2[i]);
//      fxn.imag = imag(mellinint::Np_2[i]);
//
//      ancalc_(QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F, QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, fxn);
//      anom_(ANS, AM, AP, AL, BE, AB, RMIN, RPLUS, RQQ, RQG,RGQ, RGG, C2Q, C2G, CDYQ, CDYG,
//	    fxn, nf, QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F, QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, C2QI, C2GF, CDYQI, CDYGI);
//
//      ans[i+mellinint::mdim] = cx(ANS);
//      am[i+mellinint::mdim] = cx(AM);
//      ap[i+mellinint::mdim] = cx(AP); 
//      al[i+mellinint::mdim] = cx(AL);
//      be[i+mellinint::mdim] = cx(BE);
//      ab[i+mellinint::mdim] =  cx(AB);
//      rmin[i+mellinint::mdim] = cx(RMIN);
//      rplus[i+mellinint::mdim] = cx(RPLUS);
//      rqq[i+mellinint::mdim] = cx(RQQ);
//      rqg[i+mellinint::mdim] = cx(RQG);
//      rgq[i+mellinint::mdim] = cx(RGQ);
//      rgg[i+mellinint::mdim] = cx(RGG);
//      
//      AC   = 1.- cx(AL);
//      NMP  = 1.- cx(AM) + cx(AP);
//      NPM  = 1.- cx(AP) + cx(AM);
//      DMQQ =  cx(AL) *  cx(RQQ) + cx(BE) * cx(RGQ);
//      DMQG =  cx(AL) *  cx(RQG) + cx(BE) * cx(RGG);
//      DMGQ =  cx(AB) *  cx(RQQ) +    AC  * cx(RGQ);
//      DMGG =  cx(AB) *  cx(RQG) +    AC  * cx(RGG);
//      DPQQ =     AC  *  cx(RQQ) - cx(BE) * cx(RGQ);
//      DPQG =     AC  *  cx(RQG) - cx(BE) * cx(RGG);
//      DPGQ = -cx(AB) *  cx(RQQ) + cx(AL) * cx(RGQ);
//      DPGG = -cx(AB) *  cx(RQG) + cx(AL) * cx(RGG);
//
//      rmmqq[i+mellinint::mdim]=  cx(AL) * DMQQ+cx(AB)*DMQG;
//      rmmqg[i+mellinint::mdim]=  cx(BE) * DMQQ+  AC  *DMQG;
//      rmmgq[i+mellinint::mdim]=  cx(AL) * DMGQ+cx(AB)*DMGG;
//      rmmgg[i+mellinint::mdim]=  cx(BE) * DMGQ+  AC  *DMGG;
//      rmpqq[i+mellinint::mdim]= (  AC   * DMQQ-cx(AB)*DMQG)/NMP;
//      rmpqg[i+mellinint::mdim]=(-cx(BE) * DMQQ+cx(AL)*DMQG)/NMP;
//      rmpgq[i+mellinint::mdim]= (  AC   * DMGQ-cx(AB)*DMGG)/NMP;
//      rmpgg[i+mellinint::mdim]=(-cx(BE) * DMGQ+cx(AL)*DMGG)/NMP;
//      rpmqq[i+mellinint::mdim]= (cx(AL) * DPQQ+cx(AB)*DPQG)/NPM;
//      rpmqg[i+mellinint::mdim]= (cx(BE) * DPQQ+  AC  *DPQG)/NPM;
//      rpmgq[i+mellinint::mdim]= (cx(AL) * DPGQ+cx(AB)*DPGG)/NPM;
//      rpmgg[i+mellinint::mdim]= (cx(BE) * DPGQ+  AC  *DPGG)/NPM;
//      rppqq[i+mellinint::mdim]=    AC   * DPQQ-cx(AB)*DPQG;
//      rppqg[i+mellinint::mdim]= -cx(BE) * DPQQ+cx(AL)*DPQG;
//      rppgq[i+mellinint::mdim]=    AC   * DPGQ-cx(AB)*DPGG;
//      rppgg[i+mellinint::mdim]= -cx(BE) * DPGQ+cx(AL)*DPGG;
    } 
}
void evolnative::update()
{
  //  if (opts.evolmode != 0) //--> always update for the counterterm in Mellin space
  //    return;
  //  if (opts.fmufac == 0) //This will not work when the Mellin inversion points are updated
  //    return;
  
  update_pdfs();
  update_engine();
}

void evolnative::scales()
{
  //Set scales for evolution in pdfevol
  //alphasl gives the LL/NLL/NNLL evolution of alpha from Qres=Q/a_param to  q2=b0^2/b^2

  /*******************  alpqf *************************/
  //evolnative::alpqf is used as starting scale of the PDF evolution in evolnative::evolution
  //according to Eq. 42 of arXiv:hep-ph/0508068 it should be the factorisation scale,
  //but in Eq. 98 in the resummation scale is used
  //if (opts.mufevol)
  //alpqf = resint::alpqfac;            //alphas at factorisation scale (as in Eq. 42 of arXiv:hep-ph/0508068, but see also Eq. 98)
  //else
  alpqf = resint::alpqres;              //alphas at resummation scale   (as in dyres)
  /*******************  end alpqf *************************/

  
  /**********************  alpq  **********************
  //alpq is not actually used anymore, as it is more correct to use aexp instead
  //alpq was used in hcoefficients::calcb, it is alpq = alphas(b0^2/b^2)
  //it was used only at NLL, at NNLL instead aexp*alphas(mur) is used (aexp is the same as alphasl, but with a different blim)
  //complex <double> alpq;
  alpq = resint::alpqres * cx(alphasl_(fscale2_mub)); //--> Attention! in DYRES it is alphas(qres), in DyQT it is alphas(mur)

  if (opts.evolmode == 2 || opts.evolmode == 3)
    //In order to have variable flavour evolution, use here a VFN definition of alpq
    //There is possibly an issue here when the renormalisation scale is (very) different from mll, since aass=alphas(muren) is used in xlambda = beta0*aass*blog
    {
      //      double M2 = pow(fabs(evolnative::bstartilde),2);
      double M2 = pow(fabs(evolnative::qbstar),2);                     //qbstar = b0/bstar (without a_param)
      double M2tilde;
      if (opts.modlog)
	M2tilde = M2 * resint::mures2 / (M2 + resint::mures2); //modified logs L -> L~
      else
	M2tilde = M2;
      double M2prime = M2tilde * resint::muren2/resint::mures2;     //scale for differences between muren and mures
      M2 = M2prime;
      double R2 = M2;
      double ASI;
      int NF;
      alpq = pegasus::alphas(M2, R2, ASI, NF);

      //alpq = LHAPDF::alphasPDF(sqrt(M2)) / (4.* M_PI);


      //      //tweak this to look as in dyres, where alphas(muren) is used:
      //      double R20 = resint::_m;
      //      alpq = as_(M2, R20, resint::alpqren, NF);
      //      alpq = alpq/resint::alpqren*resint::alpqres;

      
      //Based on the definition of aexp, may be need to rescale alpq for as(qres)/as(muren)?
      //alpq = alpq/resint::alpqren*resint::alpqres;
      //cout << "alpq " << alpq << " M2 " << M2 << endl;
    }
  **********************  end alpq  **********************/ 

//  /********* Scales for the evolution *****************/
//  //convert to fortran complex number
//  fcomplex fscale2_mub = fcx(pow(mub_a,2));
//  //evolnative::XL = evolnative::alpqf / alpq; // = 1./cx(alphasl_(fscale2));
//
//  XL = 1./cx(alphasl_(fscale2_mub)); //XL = alphas(mures2)/alphas(b0^2/b^2)
  XL = 1./pdfevol::asl; //XL = alphas(mures2)/alphas(b0^2/b^2)
  XL1 = 1.- XL;
  SALP = log(XL);
//  //cout << pdfevol::bscale << "  " << evolnative::XL << endl;
//  // SELECT ORDER FOR EVOLUTION LO/NLO
  if (opts.order_evol <= 1)
    alpr = 0.;
  else
    //alpr = alpqf * cx(alphasl_(fscale2_mub));
    alpr = alpqf * pdfevol::asl;
//  //force LO evolution
//  //alpr = 0.
//  //cout << b << "  " << scale2 << "  " << evolnative::SALP << "  " << log(1./cx(alphasl_(fscale2))) << "  " << evolnative::alpr << "  " << alpq <<  endl;
//  //  cout << "new " << b << "  " << evolnative::SALP << "  " << log(1./asl) << "  " << evolnative::alpr << "  " << evolnative::alpq <<  endl;
//  //**************************************

}

//Evolve the Mellin moment of PDFs at muF from Q to mub (from reno2 in DYRes)
//The output is stored in pdfevol::fx1 and pdfevol::fx2
void evolnative::evolve()
{
  if (opts.evolmode != 0)
    return;

  scales();
  
  //N flavour dependence
  int nf = resconst::NF;

  //complex <double> alpq, ALPr;

  //evolve only the moments of the positive branch, the negative branch is obtained by complex conjugation
  //int sign = mesq::positive;

  complex <double> fx[11] = {0.};

  //At LL there is no PDF evolution, PDFs are evaluated at the factorisation scale
//  if (opts.evolmode == 1 && opts.order == 0)
//    {
//      //N flavour dependence
//      int nf = resconst::NF;
//
//      complex <double> fx[11];
//      //At LL there is no PDF evolution, PDFs are evaluated at the factorisation scale
//      for (int i = 0; i < mellinint::mdim; i++)
//	pdfevol::retrievemuf(i);
//    }
  
  if (opts.order_evol == 0)
    {
      if (opts.mellin1d)
	for (int i = 0; i < mellinint::mdim; i++)
	  {
	    //XP[i] are moments of PDFs at the starting scale (factorisation scale)
	    fx[1+MAXNF] = UVP[i] + USP[i];
	    fx[2+MAXNF] = DVP[i] + DSP[i];
	    fx[3+MAXNF] = SVP[i] + SSP[i];
	    fx[0+MAXNF] = GLP[i];
	    fx[-1+MAXNF] = USP[i];
	    fx[-2+MAXNF] = DSP[i];
	    fx[-3+MAXNF] = SSP[i];
	    if (nf >= 4)
	      {
		fx[4+MAXNF] = CVP[i] + CHP[i];
		fx[-4+MAXNF] = CHP[i];
	      }
	    if (nf >= 5)
	      {
		fx[5+MAXNF] = BVP[i] + BOP[i];
		fx[-5+MAXNF] = BOP[i];
	      }
	  
	    pdfevol::storemoments(i, fx);
	  }
      else
	for (int i = 0; i < mellinint::mdim; i++)
	  {
	    //XP[i] are moments of PDFs at the starting scale (factorisation scale)
	    fx[1+MAXNF] = UVP[i] + USP[i];
	    fx[2+MAXNF] = DVP[i] + DSP[i];
	    fx[3+MAXNF] = SVP[i] + SSP[i];
	    fx[0+MAXNF] = GLP[i];
	    fx[-1+MAXNF] = USP[i];
	    fx[-2+MAXNF] = DSP[i];
	    fx[-3+MAXNF] = SSP[i];
	    if (nf >= 4)
	      {
		fx[4+MAXNF] = CVP[i] + CHP[i];
		fx[-4+MAXNF] = CHP[i];
	      }
	    if (nf >= 5)
	      {
		fx[5+MAXNF] = BVP[i] + BOP[i];
		fx[-5+MAXNF] = BOP[i];
	      }
	    pdfevol::storemoments_1(i, fx);
	    
	    fx[1+MAXNF] = UVP[i+mellinint::mdim] + USP[i+mellinint::mdim];
	    fx[2+MAXNF] = DVP[i+mellinint::mdim] + DSP[i+mellinint::mdim];
	    fx[3+MAXNF] = SVP[i+mellinint::mdim] + SSP[i+mellinint::mdim];
	    fx[0+MAXNF] = GLP[i+mellinint::mdim];
	    fx[-1+MAXNF] = USP[i+mellinint::mdim];
	    fx[-2+MAXNF] = DSP[i+mellinint::mdim];
	    fx[-3+MAXNF] = SSP[i+mellinint::mdim];
	    if (nf >= 4)
	      {
		fx[4+MAXNF] = CVP[i+mellinint::mdim] + CHP[i+mellinint::mdim];
		fx[-4+MAXNF] = CHP[i+mellinint::mdim];
	      }
	    if (nf >= 5)
	      {
		fx[5+MAXNF] = BVP[i+mellinint::mdim] + BOP[i+mellinint::mdim];
		fx[-5+MAXNF] = BOP[i+mellinint::mdim];
	      }
	    pdfevol::storemoments_2(i, fx);
	  }
	
      return;
    }
  
  //i is the index of the complex mellin moment in the z-space for the gaussian quadrature used for the mellin inversion
  for (int sign = mesq::positive; sign <= mesq::negative; sign++)
    for (int n = 0; n < ndim; n++)
      {
	// Singlet/non-singlet decomposition at the starting scale (factorisation scale)
	complex <double> UVN,DVN,GLN,SIN,NS3N,NS8N,NS15N,NS24N,NS35N;
	if (sign == mesq::positive)
	  {
	    UVN   = UVP[n];
	    DVN   = DVP[n];
	    GLN   = GLP[n];
	    SIN   = SIP[n];
	    NS3N  = NS3P[n];
	    NS8N  = NS8P[n];
	    NS15N = NS15P[n];
	    NS24N = NS24P[n];
	    NS35N = NS35P[n];
	  }
	else
	  {
	    UVN   = conj(UVP[n]);
	    DVN   = conj(DVP[n]);
	    GLN   = conj(GLP[n]);
	    SIN   = conj(SIP[n]);
	    NS3N  = conj(NS3P[n]);
	    NS8N  = conj(NS8P[n]);
	    NS15N = conj(NS15P[n]);
	    NS24N = conj(NS24P[n]);
	    NS35N = conj(NS35P[n]);
	  }
	int i = index(n,sign);
	complex <double> SG = SIN;
	complex <double> GL = GLN;

	// retrieved values cached in anomalous.C --> can move here to this module, and drop anomalous
	complex <double> ANS = ans[i];
	complex <double> AM = am[i];
	complex <double> AP = ap[i];
	complex <double> AL = al[i];
	complex <double> BE = be[i];
	complex <double> AB = ab[i];
	complex <double> AC  = 1. -AL;
	complex <double> RMIN = rmin[i];
	complex <double> RPLUS = rplus[i];
	complex <double> RMMQQ = rmmqq[i];
	complex <double> RMMQG = rmmqg[i];
	complex <double> RMMGQ = rmmgq[i];
	complex <double> RMMGG = rmmgg[i];
	complex <double> RMPQQ = rmpqq[i];
	complex <double> RMPQG = rmpqg[i];
	complex <double> RMPGQ = rmpgq[i];
	complex <double> RMPGG = rmpgg[i];
	complex <double> RPMQQ = rpmqq[i];
	complex <double> RPMQG = rpmqg[i];
	complex <double> RPMGQ = rpmgq[i];
	complex <double> RPMGG = rpmgg[i];
	complex <double> RPPQQ = rppqq[i];
	complex <double> RPPQG = rppqg[i];
	complex <double> RPPGQ = rppgq[i];
	complex <double> RPPGG = rppgg[i];
	// **************************************

	// **************************************
	//     b-dependence
	//resummation scale
	//  complex <double> XL = 1./cx(alphasl_(fcx(scale2)));
	//  complex <double> XL1 = 1.- XL;
	//  complex <double> SALP = log(XL);
	//--> SALP ~ log[alphas(Q)/alphas(b0/b)]
  
	complex <double> S = SALP;
	//  cout << S << "  " << <<alpr <<  endl;

	complex <double> ENS = exp(-ANS*S);
	complex <double> EM  = exp(-AM*S);
	complex <double> EP  = exp(-AP*S);
	complex <double> EMP = EM/EP;
	complex <double> EPM = EP/EM;

	//...EVOLUTION OF LIGHT PARTON DENSITIES
	//double q2s = q2/pow(resint::a,2);                //resummation scale
	//double alpqf = dyalphas_lhapdf_(sqrt(q2s))/4./M_PI; //alphas at the resummation scale
	//complex <double> alpq = alpqf * alphasl(scale2);              //alphas at the resummation scale times alphas at 1/b
	//complex <double> alpr= alpq * 1 *(opts.order-1);
	//--> alpr = 0 at NLL; alpr = alphas(Q) * alphasl ~ alphas(b0/b) at NNLL

	if (opts.order_evol == 1)
	  {
	    UVN  = UVN  * ENS;
	    DVN  = DVN  * ENS;
	    NS3N = NS3N * ENS;
	    NS8N = NS8N * ENS;
	    
	    SIN = EM * (AL*SG + BE*GL) + EP * (AC*SG - BE*GL);
	    GLN = EM * (AB*SG + AC*GL) + EP * (-AB*SG + AL*GL);
	    
	    NS15N = NS15N * ENS;
	    NS24N = NS24N * ENS;
	  }
	else if (opts.order_evol == 2)
	  {
	    UVN  = UVN  * ENS * (1.+  alpr * XL1 * RMIN);
	    DVN  = DVN  * ENS * (1.+  alpr * XL1 * RMIN);
	    NS3N = NS3N * ENS * (1.+  alpr * XL1 * RPLUS);
	    NS8N = NS8N * ENS * (1.+  alpr * XL1 * RPLUS);
	    
	    SIN = EM * ((AL + alpr * (RMMQQ*XL1 + RMPQQ*(EPM-XL)))* SG + (BE + alpr * (RMMQG*XL1 + RMPQG*(EPM-XL))) * GL)
	      + EP * ((AC + alpr * (RPPQQ*XL1 + RPMQQ*(EMP-XL)))* SG + (-BE + alpr * (RPPQG*XL1 + RPMQG*(EMP-XL))) * GL);
	    GLN = EM * ((AB + alpr * (RMMGQ*XL1 + RMPGQ*(EPM-XL)))* SG + (AC + alpr * (RMMGG*XL1 + RMPGG*(EPM-XL))) * GL)
	      + EP *((-AB + alpr * (RPPGQ*XL1 + RPMGQ*(EMP-XL)))* SG + (AL + alpr * (RPPGG*XL1 + RPMGG*(EMP-XL))) * GL);
  
	    NS15N = NS15N * ENS * (1.+  alpr * XL1 * RPLUS);
	    NS24N = NS24N * ENS * (1.+  alpr * XL1 * RPLUS);
	  }

      
	NS35N = SIN;

	//...  FLAVOUR DECOMPOSITION OF THE QUARK SEA :
	complex <double> SSN, DSN, USN, CHN, BON;
	SSN = (10.* SIN + 2.* NS35N + 3.* NS24N + 5.* NS15N - 20.* NS8N) / 120.;
	DSN = (10.* SIN + 2.* NS35N + 3.* NS24N + 5.* NS15N + 10.* NS8N - 30.* NS3N - 60.* DVN) / 120.;
	USN = (10.* SIN + 2.* NS35N + 3.* NS24N + 5.* NS15N + 10.* NS8N + 30.* NS3N - 60.* UVN) / 120.;
      
	CHN = (UVN + DVN + 2.*USN + 2.*DSN + 2.*SSN - NS15N)/6.;
	BON = (UVN + DVN + 2.*USN + 2.*DSN + 2.*SSN + 2.*CHN - NS24N)/8.;
	//equivalent to:
	//CHN = (10.* SIN + 2. *NS35N + 3.* NS24N - 15.* NS15N) / 120.;
	//BON = (10.* SIN + 2. *NS35N - 12.* NS24N) / 120.;

	if (nf == 3) //GRV
	  {
	    SSN= (20.* SIN - 20.* NS8N)/120.;
	    DSN = (20.* SIN + 10.* NS8N - 30.* NS3N - 60.* DVN) / 120.;
	    USN = (20.* SIN + 10.* NS8N + 30.* NS3N - 60.* UVN) / 120.;
	    CHN=0.;
	    BON=0.;
	  }

	//if (fabs(bstarscale) < LHAPDF::getThreshold(4))
	//    CHN *= exp(-pow((LHAPDF::getThreshold(4)-fabs(bstarscale)),2)/pow((LHAPDF::getThreshold(4)/10.) ,2)); //=0
	//  double delta = 1./5.;
	//  if (fabs(bstarscale) < LHAPDF::getThreshold(5)*(1+delta/2.))
	//    BON *= exp(-pow((LHAPDF::getThreshold(5)*(1+delta/2.)-fabs(bstarscale)),2)/pow((LHAPDF::getThreshold(5)*delta),2)); //=0

	// **************************************

	// output:
	// bbar cbar sbar dbar ubar gluon  u   d   s   c   b
	// -5    -4   -3   -2   -1    0    1   2   3   4   5

	fx[0+MAXNF] = GLN;
	fx[1+MAXNF] = UVN + USN;
	fx[-1+MAXNF] = USN;
	fx[2+MAXNF] = DVN + DSN;
	fx[-2+MAXNF] = DSN;
	fx[3+MAXNF] = SSN;
	fx[-3+MAXNF] = SSN;
	if (nf >= 4)
	  {
	    fx[4+MAXNF] = CHN;
	    fx[-4+MAXNF] = CHN;
	  }
	else
	  {
	    fx[4+MAXNF] = 0.;
	    fx[-4+MAXNF] = 0.;
	  }
  
	if (nf >= 5)
	  {
	    fx[5+MAXNF] = BON;
	    fx[-5+MAXNF] = BON;
	  }
	else
	  {
	    fx[5+MAXNF] = 0.;
	    fx[-5+MAXNF] = 0.;
	  }

	if (opts.mellin1d)
	  pdfevol::storemoments(n, sign, fx);
	else
	  {
	    if (n < mellinint::mdim)
	      pdfevol::storemoments_1(n, sign, fx);
	    else
	      pdfevol::storemoments_2(n-mellinint::mdim, sign, fx);
	  }
      }
}


void evolnative::init_fortran()
{
  //Old Fortran code for the moments
  fcomplex uval,dval,usea,dsea,splus,ssea,glu,charm,bot;
  clock_t begin_time, end_time;
  begin_time = clock();
  double xmin = pow(bins.mbins.front()/opts.sroot,2); //Restrict the integration of moments to xmin = m/sqrt(s)*exp(-ymax) = (m/sqrt(s))^2
  for (int k = 0; k < mellinint::mdim; k++)
    {
      int hadron = 1; //opts.ih1;
      fcomplex XN = fcx(mellinint::Np[k]); //compute positive branch only, the negative branch is obtained by complex conjugation
      double facscale = opts.kmufac*opts.rmass;
      pdfmoments_(hadron,facscale,XN,uval,dval,usea,dsea,splus,ssea,glu,charm,bot,xmin);
      // cout << "moment " << k << "  " << cx(XN) << "  ";
      // cout << "uval  " << cx(uval) << endl;
      // cout << "dval  " << dval << endl;
      // cout << "usea  " << usea << endl;
      // cout << "dsea  " << dsea << endl;
      // cout << "gluon " << cx(glu) << endl;
      // cout << "charm " << charm << endl;
      // cout << "bottom" << bot << endl;
      UVP[k] = cx(uval);
      DVP[k] = cx(dval);
      USP[k] = cx(usea);
      DSP[k] = cx(dsea);
      SSP[k] = cx(ssea);//SSP[k] = cx(splus); !! issue for PDFs with s-sbar asymmetry !!
      GLP[k] = cx(glu);
      CHP[k] = cx(charm);
      BOP[k] = cx(bot);
    }
  end_time = clock();
  cout << "Done " << float(end_time - begin_time)/CLOCKS_PER_SEC << endl;

  /*
  //polynomial interpolation for analytical continuation
  int order = 500;
  complex <double> a = opts.cpoint+1-2*opts.zmax*1i;
  complex <double> b = opts.cpoint+1+2*opts.zmax*1i;
  complex <double> c = 0.5*(a+b);
  complex <double> m = 0.5*(b-a);
  complex <double> f[order];
  for (int i = 1; i <= order; i++)
    {
      complex <double> x = c+m*cheb::xxx[order-1][i-1];
      int hadron = 1; //opts.ih1;
      fcomplex XN = fcx(x);
      double facscale = 1.0;//opts.kmufac*opts.rmass;
      pdfmoments_(hadron,facscale,XN,uval,dval,usea,dsea,splus,ssea,glu,charm,bot);
      f[i-1] = cx(uval);
    }
  
  for (int i = 0; i < 100; i++)
    {
      complex <double> x = a+(b-a)*double(i)/100.;
      int hadron = 1; //opts.ih1;
      fcomplex XN = fcx(x);
      double facscale = 1.0;//opts.kmufac*opts.rmass;
      pdfmoments_(hadron,facscale,XN,uval,dval,usea,dsea,splus,ssea,glu,charm,bot);
      cout << setw(30) << x
	   << setw(30) << cx(uval)
	   << setw(30) << cheb::ipol(a,b,order, f, x)
	   << setw(30) << cx(uval)-cheb::ipol(a,b,order, f, x) <<  endl;
    }
  for (int i = 1; i <= order; i++)
    {
      complex <double> x = c+m*cheb::xxx[order-1][i-1];
      cout << x << "  " << f[i-1] << "  " << cheb::ipol(a,b,order, f, x) << endl;
    }

  //verify
  for (int k = 0; k < mellinint::mdim; k++)
    {
      complex <double> z = mellinint::Np[k];
      cout << setw(30) << z
	   << setw(30) << UVP[k]
	   << setw(30) << cheb::ipol(a,b,order, f, z)
	   << endl;
    }
  */
}
