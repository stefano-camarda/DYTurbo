#include "pdfevol.h"
#include "interface.h"
#include "settings.h"
#include "mesq.h"
#include "anomalous.h"
#include "resconst.h"
#include "chebyshev.h"
#include "phasespace.h"

#include <LHAPDF/LHAPDF.h>

#include <iostream>
#include <complex>

//PDFs mellin moments at the factorisation scale
complex <double> *pdfevol::UVP;
complex <double> *pdfevol::DVP;
complex <double> *pdfevol::USP;
complex <double> *pdfevol::DSP;
complex <double> *pdfevol::SSP;
complex <double> *pdfevol::GLP;
complex <double> *pdfevol::CHP;
complex <double> *pdfevol::BOP;

complex <double> pdfevol::fn1[2*MAXNF+1];
complex <double> pdfevol::fn2[2*MAXNF+1];

//scales
complex <double> pdfevol::bscale;
complex <double> pdfevol::bstarscale;
complex <double> pdfevol::bstartilde;
complex <double> pdfevol::qbstar;
complex <double> pdfevol::bcomplex;

complex <double> pdfevol::XL;
complex <double> pdfevol::XL1;
complex <double> pdfevol::SALP;

complex <double> pdfevol::alpr;

//fortran interface
void pdfevol_(int& i1, int& i2, int& sign)
{
  pdfevol::retrieve(i1-1, i2-1, sign-1);
};

void pdfevol::init()
{
  UVP = new complex <double>[mellinint::mdim];
  DVP = new complex <double>[mellinint::mdim];
  USP = new complex <double>[mellinint::mdim];
  DSP = new complex <double>[mellinint::mdim];
  SSP = new complex <double>[mellinint::mdim];
  GLP = new complex <double>[mellinint::mdim];
  CHP = new complex <double>[mellinint::mdim];
  BOP = new complex <double>[mellinint::mdim];

  //calculate Mellin moments of PDFs
  cout << "Initialise PDF moments with numerical integration (C++)... " << flush;
  
  //double xmin = 1e-8;
  double xmin = pow(bins.mbins.front()/opts.sroot,2); //Restrict the integration of moments to xmin = m/sqrt(s)*exp(-ymax) = (m/sqrt(s))^2

  fcomplex uval,dval,usea,dsea,splus,ssea,glu,charm,bot;
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
      SSP[k] = cx(ssea);
      GLP[k] = cx(glu);
      CHP[k] = cx(charm);
      BOP[k] = cx(bot);
    }

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
  
  cout << "Done" << endl;
}

// PROVIDES MOMENTS OF DENSITIES AT A GIVEN SCALE (INCLUDES EVOLUTION)
// AND ANOMALOUS DIMENSIONS
// EVERYTHING IN MELLIN SPACE
//
//     input: I point in z grid, ISIGN (+ or -), IBEAM (1 or 2), SCALE2 (b), alps = alpqf (mass dependent)
//     output: FN, alpq = alps * alphasl(scale2)
//     dependence: FN(b, mass, I)
//
// ***************************************
// Can completely cache the output FN in FN(I,IFIT,ISIG), and calculate in the init? No because of b
// reno2 is cached and looped only on I, before entering the I1 I2 doube loop into cfx1 cfx2p cfx2m

//Evolve the Mellin moment of PDF corresponding to the index i, sign sign, beam beam, from the scale q2 to the scale scale2(b)
//Output in fx
void pdfevol::evolution(int i) //from reno2
{
  // i is the index of the complex mellin moment in the z-space for the gaussian quadrature used for the mellin inversion

  //N flavour dependence
  int nf = resconst::NF;

  //At LL there is no PDF evolution, PDFs are evaluated at the factorisation scale
  if (opts.order == 0)
    {
      //XP[i] are moments of PDFs at the starting scale (factorisation scale)
      complex <double> fx[11];
      fx[0+MAXNF] = GLP[i];
      fx[1+MAXNF] = UVP[i] + USP[i];
      fx[-1+MAXNF] = USP[i];
      fx[2+MAXNF] = DVP[i] + DSP[i];
      fx[-2+MAXNF] = DSP[i];
      fx[3+MAXNF] = SSP[i];
      fx[-3+MAXNF] = SSP[i];
      if (nf >= 4)
	{
	  fx[4+MAXNF] = CHP[i];
	  fx[-4+MAXNF] = CHP[i];
	}
      else
	{
	  fx[4+MAXNF] = 0.;
	  fx[-4+MAXNF] = 0.;
	}
      if (nf >= 5)
	{
	  fx[5+MAXNF] = BOP[i];
	  fx[-5+MAXNF] = BOP[i];
	}
      else
	{
	  fx[5+MAXNF] = 0.;
	  fx[-5+MAXNF] = 0.;
	}

      storemoments(i, fx);
      return;
    }

  complex <double> alpq, ALPr;

  //Moments at the factorisation scale
  complex <double> UVI,	DVI, USI, DSI, SSI, GLI, CHI, BOI;

  //evolve only the moments of the positive branch, the negative branch is obtained by complex conjugation
  int sign = mesq::positive;

  //Moments of PDFs at the starting scale (factorisation scale)
  UVI = UVP[i];
  DVI = DVP[i];
  USI = USP[i];
  DSI = DSP[i];
  SSI = SSP[i];
  GLI = GLP[i];
  CHI = CHP[i];
  BOI = BOP[i];

  // ***********************
  // this part can be precomputed
  complex <double> UVN = UVI;
  complex <double> DVN = DVI;
  complex <double> NS3N = UVI + 2.*USI - DVI - 2.*DSI;            //(u+ub-d-db)               u-d
  complex <double> NS8N = UVI + 2.*USI + DVI + 2.*DSI - 4.*SSI;   //(u+ub+d+db-s-sb)          u+d-s
  complex <double> GLN = GLI;

  complex <double> SIN, NS15N, NS24N, NS35N;
  if (nf == 5)
    {
      SIN = UVI + DVI + 2.*USI + 2.*DSI + 2.*SSI + 2.*CHI + 2.*BOI;   //(u+ub+d+db+s+sb+c+cb+b+bb) -> all quarks      u+d+s+c+b
      NS15N = UVI + DVI + 2.*USI + 2.*DSI + 2.*SSI - 6.*CHI;          //(u+ub+d+db+s+sb-3(c+cb))                      u+d+s-3c
      NS24N = UVI + DVI + 2.*USI + 2.*DSI + 2.*SSI + 2.*CHI - 8.*BOI; //(u+ub+d+db+s+sb+c+cb-4(b+bb))                 u+d+s+c-4b
      NS35N = SIN;
    }
  if (nf == 4)
    {
      SIN = UVI + DVI + 2.*USI + 2.*DSI + 2.*SSI + 2.*CHI;          //(u+ub+d+db+s+sb+c+cb) -> all quarks
      NS15N = UVI + DVI + 2.*USI + 2.*DSI + 2.*SSI - 6.*CHI;
      NS24N = UVI + DVI + 2.*USI + 2.*DSI + 2.*SSI + 2.*CHI;
      NS35N = SIN;
    }
  if (nf == 3)
    {
      SIN = UVI + DVI + 2.*USI + 2.*DSI + 2.*SSI;     //(u+ub+d+db+s+sb) -> all quarks
      NS15N = SIN;
      NS24N = SIN;
      NS35N = SIN;
    }
  complex <double> SG = SIN;
  complex <double> GL = GLN;

  // **************************************
  // retrieved values cached in anomalous.C
  complex <double> ANS = anomalous::ans[anomalous::index(i,sign)];
  complex <double> AM = anomalous::am[anomalous::index(i,sign)];
  complex <double> AP = anomalous::ap[anomalous::index(i,sign)];
  complex <double> AL = anomalous::al[anomalous::index(i,sign)];
  complex <double> BE = anomalous::be[anomalous::index(i,sign)];
  complex <double> AB = anomalous::ab[anomalous::index(i,sign)];
  complex <double> AC  = 1. -AL;
  complex <double> RMIN = anomalous::rmin[anomalous::index(i,sign)];
  complex <double> RPLUS = anomalous::rplus[anomalous::index(i,sign)];
  complex <double> RMMQQ = anomalous::RMMQQ[anomalous::index(i,sign)];
  complex <double> RMMQG = anomalous::RMMQG[anomalous::index(i,sign)];
  complex <double> RMMGQ = anomalous::RMMGQ[anomalous::index(i,sign)];
  complex <double> RMMGG = anomalous::RMMGG[anomalous::index(i,sign)];
  complex <double> RMPQQ = anomalous::RMPQQ[anomalous::index(i,sign)];
  complex <double> RMPQG = anomalous::RMPQG[anomalous::index(i,sign)];
  complex <double> RMPGQ = anomalous::RMPGQ[anomalous::index(i,sign)];
  complex <double> RMPGG = anomalous::RMPGG[anomalous::index(i,sign)];
  complex <double> RPMQQ = anomalous::RPMQQ[anomalous::index(i,sign)];
  complex <double> RPMQG = anomalous::RPMQG[anomalous::index(i,sign)];
  complex <double> RPMGQ = anomalous::RPMGQ[anomalous::index(i,sign)];
  complex <double> RPMGG = anomalous::RPMGG[anomalous::index(i,sign)];
  complex <double> RPPQQ = anomalous::RPPQQ[anomalous::index(i,sign)];
  complex <double> RPPQG = anomalous::RPPQG[anomalous::index(i,sign)];
  complex <double> RPPGQ = anomalous::RPPGQ[anomalous::index(i,sign)];
  complex <double> RPPGG = anomalous::RPPGG[anomalous::index(i,sign)];
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

  complex <double> fx[11];
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

  storemoments(i, fx);
}


void pdfevol::calculate(int i)
{
  //N flavour dependence
  int nf = resconst::NF;

  fcomplex uval,dval,usea,dsea,s,sbar,glu,charm,bot;

  int hadron = 1;
  //double facscale = fabs(bscale);
  double facscale = fabs(bstarscale);//better qbstar here? actually better bstartilde = qbstar * resint::mures / sqrt(pow(qbstar,2) + resint::mures2);
  //double facscale = fabs(opts.muf);
  facscale = fabs(pdfevol::bstartilde);
  fcomplex XN = fcx(mellinint::Np[i]);
  double xmin = 1e-8;
  pdfmoments_(hadron,facscale,XN,uval,dval,usea,dsea,s,sbar,glu,charm,bot,xmin);

  complex <double> fx[11];
  fx[0+MAXNF] = cx(glu);
  fx[1+MAXNF] = cx(uval) + cx(usea);
  fx[-1+MAXNF] = cx(usea);
  fx[2+MAXNF] = cx(dval) + cx(dsea);
  fx[-2+MAXNF] = cx(dsea);
  fx[3+MAXNF] = cx(s);
  fx[-3+MAXNF] = cx(sbar);
  if (nf >= 4)
    {
      fx[4+MAXNF] = cx(charm);
      fx[-4+MAXNF] = cx(charm);
    }
  else
    {
      fx[4+MAXNF] = 0.;
      fx[-4+MAXNF] = 0.;
    }
  
  if (nf >= 5)
    {
      fx[5+MAXNF] = cx(bot);
      fx[-5+MAXNF] = cx(bot);
    }
  else
    {
      fx[5+MAXNF] = 0.;
      fx[-5+MAXNF] = 0.;
    }

  storemoments(i, fx);
}


void pdfevol::storemoments(int i, complex <double> fx[11])
{
  //Save the evolved PDFs into the fortran common block (can actually use a C++ data format since it is only accessed in C++)

  //beam 1
  creno_.cfx1_[i][-5+MAXNF] = fcx(fx[-5+MAXNF]);
  creno_.cfx1_[i][-4+MAXNF] = fcx(fx[-4+MAXNF]);
  creno_.cfx1_[i][-3+MAXNF] = fcx(fx[-3+MAXNF]);
  creno_.cfx1_[i][ 0+MAXNF] = fcx(fx[ 0+MAXNF]);
  creno_.cfx1_[i][ 3+MAXNF] = fcx(fx[ 3+MAXNF]);
  creno_.cfx1_[i][ 4+MAXNF] = fcx(fx[ 4+MAXNF]);
  creno_.cfx1_[i][ 5+MAXNF] = fcx(fx[ 5+MAXNF]);
  if (opts.ih1 == 1)
    {
      creno_.cfx1_[i][-2+MAXNF] = fcx(fx[-2+MAXNF]);
      creno_.cfx1_[i][-1+MAXNF] = fcx(fx[-1+MAXNF]);
      creno_.cfx1_[i][ 1+MAXNF] = fcx(fx[ 1+MAXNF]);
      creno_.cfx1_[i][ 2+MAXNF] = fcx(fx[ 2+MAXNF]);
    }
  else if (opts.ih1 == -1)
    {
      creno_.cfx1_[i][-2+MAXNF] = fcx(fx[ 2+MAXNF]);
      creno_.cfx1_[i][-1+MAXNF] = fcx(fx[ 1+MAXNF]);
      creno_.cfx1_[i][ 1+MAXNF] = fcx(fx[-1+MAXNF]);
      creno_.cfx1_[i][ 2+MAXNF] = fcx(fx[-2+MAXNF]);
    }

  //beam 2 positive
  creno_.cfx2p_[i][-5+MAXNF] = fcx(fx[-5+MAXNF]);
  creno_.cfx2p_[i][-4+MAXNF] = fcx(fx[-4+MAXNF]);
  creno_.cfx2p_[i][-3+MAXNF] = fcx(fx[-3+MAXNF]);
  creno_.cfx2p_[i][ 0+MAXNF] = fcx(fx[ 0+MAXNF]);
  creno_.cfx2p_[i][ 3+MAXNF] = fcx(fx[ 3+MAXNF]);
  creno_.cfx2p_[i][ 4+MAXNF] = fcx(fx[ 4+MAXNF]);
  creno_.cfx2p_[i][ 5+MAXNF] = fcx(fx[ 5+MAXNF]);
  if (opts.ih2 == 1)
    {
      creno_.cfx2p_[i][-2+MAXNF] = fcx(fx[-2+MAXNF]);
      creno_.cfx2p_[i][-1+MAXNF] = fcx(fx[-1+MAXNF]);
      creno_.cfx2p_[i][ 1+MAXNF] = fcx(fx[ 1+MAXNF]);
      creno_.cfx2p_[i][ 2+MAXNF] = fcx(fx[ 2+MAXNF]);
    }
  else if (opts.ih2 == -1)
    {
      creno_.cfx2p_[i][-2+MAXNF] = fcx(fx[ 2+MAXNF]);
      creno_.cfx2p_[i][-1+MAXNF] = fcx(fx[ 1+MAXNF]);
      creno_.cfx2p_[i][ 1+MAXNF] = fcx(fx[-1+MAXNF]);
      creno_.cfx2p_[i][ 2+MAXNF] = fcx(fx[-2+MAXNF]);
    }

  //beam 2 negative
  creno_.cfx2m_[i][-5+MAXNF] = fcx(conj(fx[-5+MAXNF]));
  creno_.cfx2m_[i][-4+MAXNF] = fcx(conj(fx[-4+MAXNF]));
  creno_.cfx2m_[i][-3+MAXNF] = fcx(conj(fx[-3+MAXNF]));
  creno_.cfx2m_[i][ 0+MAXNF] = fcx(conj(fx[ 0+MAXNF]));
  creno_.cfx2m_[i][ 3+MAXNF] = fcx(conj(fx[ 3+MAXNF]));
  creno_.cfx2m_[i][ 4+MAXNF] = fcx(conj(fx[ 4+MAXNF]));
  creno_.cfx2m_[i][ 5+MAXNF] = fcx(conj(fx[ 5+MAXNF]));
  if (opts.ih2 == 1)
    {
      creno_.cfx2m_[i][-2+MAXNF] = fcx(conj(fx[-2+MAXNF]));
      creno_.cfx2m_[i][-1+MAXNF] = fcx(conj(fx[-1+MAXNF]));
      creno_.cfx2m_[i][ 1+MAXNF] = fcx(conj(fx[ 1+MAXNF]));
      creno_.cfx2m_[i][ 2+MAXNF] = fcx(conj(fx[ 2+MAXNF]));
    }
  else if (opts.ih2 == -1)
    {
      creno_.cfx2m_[i][-2+MAXNF] = fcx(conj(fx[ 2+MAXNF]));
      creno_.cfx2m_[i][-1+MAXNF] = fcx(conj(fx[ 1+MAXNF]));
      creno_.cfx2m_[i][ 1+MAXNF] = fcx(conj(fx[-1+MAXNF]));
      creno_.cfx2m_[i][ 2+MAXNF] = fcx(conj(fx[-2+MAXNF]));
    }
}

void pdfevol::retrieve(int i1, int i2, int sign)
{
  //  cout << i1 << endl;
  //  cout << creno_.cfx1_[i1][5].real << "  " << creno_.cfx1_[i1][5].imag << endl;
  fn1[-5+MAXNF] = cx(creno_.cfx1_[i1][-5+MAXNF]);
  fn1[-4+MAXNF] = cx(creno_.cfx1_[i1][-4+MAXNF]);
  fn1[-3+MAXNF] = cx(creno_.cfx1_[i1][-3+MAXNF]);
  fn1[-2+MAXNF] = cx(creno_.cfx1_[i1][-2+MAXNF]);
  fn1[-1+MAXNF] = cx(creno_.cfx1_[i1][-1+MAXNF]);
  fn1[ 0+MAXNF] = cx(creno_.cfx1_[i1][ 0+MAXNF]);
  fn1[ 1+MAXNF] = cx(creno_.cfx1_[i1][ 1+MAXNF]);
  fn1[ 2+MAXNF] = cx(creno_.cfx1_[i1][ 2+MAXNF]);
  fn1[ 3+MAXNF] = cx(creno_.cfx1_[i1][ 3+MAXNF]);
  fn1[ 4+MAXNF] = cx(creno_.cfx1_[i1][ 4+MAXNF]);
  fn1[ 5+MAXNF] = cx(creno_.cfx1_[i1][ 5+MAXNF]);
  if (sign == mesq::positive)
    {
      fn2[-5+MAXNF] = cx(creno_.cfx2p_[i2][-5+MAXNF]);
      fn2[-4+MAXNF] = cx(creno_.cfx2p_[i2][-4+MAXNF]);
      fn2[-3+MAXNF] = cx(creno_.cfx2p_[i2][-3+MAXNF]);
      fn2[-2+MAXNF] = cx(creno_.cfx2p_[i2][-2+MAXNF]);
      fn2[-1+MAXNF] = cx(creno_.cfx2p_[i2][-1+MAXNF]);
      fn2[ 0+MAXNF] = cx(creno_.cfx2p_[i2][ 0+MAXNF]);
      fn2[ 1+MAXNF] = cx(creno_.cfx2p_[i2][ 1+MAXNF]);
      fn2[ 2+MAXNF] = cx(creno_.cfx2p_[i2][ 2+MAXNF]);
      fn2[ 3+MAXNF] = cx(creno_.cfx2p_[i2][ 3+MAXNF]);
      fn2[ 4+MAXNF] = cx(creno_.cfx2p_[i2][ 4+MAXNF]);
      fn2[ 5+MAXNF] = cx(creno_.cfx2p_[i2][ 5+MAXNF]);
    }
  else if (sign == mesq::negative)
    {
      fn2[-5+MAXNF] = cx(creno_.cfx2m_[i2][-5+MAXNF]);
      fn2[-4+MAXNF] = cx(creno_.cfx2m_[i2][-4+MAXNF]);
      fn2[-3+MAXNF] = cx(creno_.cfx2m_[i2][-3+MAXNF]);
      fn2[-2+MAXNF] = cx(creno_.cfx2m_[i2][-2+MAXNF]);
      fn2[-1+MAXNF] = cx(creno_.cfx2m_[i2][-1+MAXNF]);
      fn2[ 0+MAXNF] = cx(creno_.cfx2m_[i2][ 0+MAXNF]);
      fn2[ 1+MAXNF] = cx(creno_.cfx2m_[i2][ 1+MAXNF]);
      fn2[ 2+MAXNF] = cx(creno_.cfx2m_[i2][ 2+MAXNF]);
      fn2[ 3+MAXNF] = cx(creno_.cfx2m_[i2][ 3+MAXNF]);
      fn2[ 4+MAXNF] = cx(creno_.cfx2m_[i2][ 4+MAXNF]);
      fn2[ 5+MAXNF] = cx(creno_.cfx2m_[i2][ 5+MAXNF]);
    }
  //set b to 0
  //  fn2[-5+MAXNF] = 0;  fn1[-5+MAXNF] = 0;
  //  fn2[5+MAXNF]  = 0;  fn1[5+MAXNF]  = 0;
  
  //set s and c to 0
  //  fn2[-4+MAXNF] = 0;  fn1[-4+MAXNF] = 0;
  //  fn2[-3+MAXNF] = 0;  fn1[-3+MAXNF] = 0;
  //  fn2[3+MAXNF]  = 0;  fn1[3+MAXNF]  = 0;
  //  fn2[4+MAXNF]  = 0;  fn1[4+MAXNF]  = 0;

  //set u and d to 0
  //  fn2[-2+MAXNF] = 0; fn1[-2+MAXNF] = 0;
  //  fn2[-1+MAXNF] = 0; fn1[-1+MAXNF] = 0;
  //  fn2[1+MAXNF]  = 0; fn1[1+MAXNF]  = 0;
  //  fn2[2+MAXNF]  = 0; fn1[2+MAXNF]  = 0;
}


void pdfevol::retrieve1d(int i, int sign)
{
  //  cout << i1 << endl;
  //  cout << creno_.cfx1_[i1][5].real << "  " << creno_.cfx1_[i1][5].imag << endl;
  if (sign == mesq::positive)
    {
      fn1[-5+MAXNF] = cx(creno_.cfx1_[i][-5+MAXNF]);
      fn1[-4+MAXNF] = cx(creno_.cfx1_[i][-4+MAXNF]);
      fn1[-3+MAXNF] = cx(creno_.cfx1_[i][-3+MAXNF]);
      fn1[-2+MAXNF] = cx(creno_.cfx1_[i][-2+MAXNF]);
      fn1[-1+MAXNF] = cx(creno_.cfx1_[i][-1+MAXNF]);
      fn1[ 0+MAXNF] = cx(creno_.cfx1_[i][ 0+MAXNF]);
      fn1[ 1+MAXNF] = cx(creno_.cfx1_[i][ 1+MAXNF]);
      fn1[ 2+MAXNF] = cx(creno_.cfx1_[i][ 2+MAXNF]);
      fn1[ 3+MAXNF] = cx(creno_.cfx1_[i][ 3+MAXNF]);
      fn1[ 4+MAXNF] = cx(creno_.cfx1_[i][ 4+MAXNF]);
      fn1[ 5+MAXNF] = cx(creno_.cfx1_[i][ 5+MAXNF]);
    }
  else if (sign == mesq::negative)
    {
      fn1[-5+MAXNF] = conj(cx(creno_.cfx1_[i][-5+MAXNF]));
      fn1[-4+MAXNF] = conj(cx(creno_.cfx1_[i][-4+MAXNF]));
      fn1[-3+MAXNF] = conj(cx(creno_.cfx1_[i][-3+MAXNF]));
      fn1[-2+MAXNF] = conj(cx(creno_.cfx1_[i][-2+MAXNF]));
      fn1[-1+MAXNF] = conj(cx(creno_.cfx1_[i][-1+MAXNF]));
      fn1[ 0+MAXNF] = conj(cx(creno_.cfx1_[i][ 0+MAXNF]));
      fn1[ 1+MAXNF] = conj(cx(creno_.cfx1_[i][ 1+MAXNF]));
      fn1[ 2+MAXNF] = conj(cx(creno_.cfx1_[i][ 2+MAXNF]));
      fn1[ 3+MAXNF] = conj(cx(creno_.cfx1_[i][ 3+MAXNF]));
      fn1[ 4+MAXNF] = conj(cx(creno_.cfx1_[i][ 4+MAXNF]));
      fn1[ 5+MAXNF] = conj(cx(creno_.cfx1_[i][ 5+MAXNF]));
    }
  if (sign == mesq::positive)
    {
      fn2[-5+MAXNF] = cx(creno_.cfx2p_[i][-5+MAXNF]);
      fn2[-4+MAXNF] = cx(creno_.cfx2p_[i][-4+MAXNF]);
      fn2[-3+MAXNF] = cx(creno_.cfx2p_[i][-3+MAXNF]);
      fn2[-2+MAXNF] = cx(creno_.cfx2p_[i][-2+MAXNF]);
      fn2[-1+MAXNF] = cx(creno_.cfx2p_[i][-1+MAXNF]);
      fn2[ 0+MAXNF] = cx(creno_.cfx2p_[i][ 0+MAXNF]);
      fn2[ 1+MAXNF] = cx(creno_.cfx2p_[i][ 1+MAXNF]);
      fn2[ 2+MAXNF] = cx(creno_.cfx2p_[i][ 2+MAXNF]);
      fn2[ 3+MAXNF] = cx(creno_.cfx2p_[i][ 3+MAXNF]);
      fn2[ 4+MAXNF] = cx(creno_.cfx2p_[i][ 4+MAXNF]);
      fn2[ 5+MAXNF] = cx(creno_.cfx2p_[i][ 5+MAXNF]);
    }
  else if (sign == mesq::negative)
    {
      fn2[-5+MAXNF] = cx(creno_.cfx2m_[i][-5+MAXNF]);
      fn2[-4+MAXNF] = cx(creno_.cfx2m_[i][-4+MAXNF]);
      fn2[-3+MAXNF] = cx(creno_.cfx2m_[i][-3+MAXNF]);
      fn2[-2+MAXNF] = cx(creno_.cfx2m_[i][-2+MAXNF]);
      fn2[-1+MAXNF] = cx(creno_.cfx2m_[i][-1+MAXNF]);
      fn2[ 0+MAXNF] = cx(creno_.cfx2m_[i][ 0+MAXNF]);
      fn2[ 1+MAXNF] = cx(creno_.cfx2m_[i][ 1+MAXNF]);
      fn2[ 2+MAXNF] = cx(creno_.cfx2m_[i][ 2+MAXNF]);
      fn2[ 3+MAXNF] = cx(creno_.cfx2m_[i][ 3+MAXNF]);
      fn2[ 4+MAXNF] = cx(creno_.cfx2m_[i][ 4+MAXNF]);
      fn2[ 5+MAXNF] = cx(creno_.cfx2m_[i][ 5+MAXNF]);
    }
  //set b to 0
  //  fn2[-5+MAXNF] = 0;  fn1[-5+MAXNF] = 0;
  //  fn2[5+MAXNF]  = 0;  fn1[5+MAXNF]  = 0;
  
  //set s and c to 0
  //  fn2[-4+MAXNF] = 0;  fn1[-4+MAXNF] = 0;
  //  fn2[-3+MAXNF] = 0;  fn1[-3+MAXNF] = 0;
  //  fn2[3+MAXNF]  = 0;  fn1[3+MAXNF]  = 0;
  //  fn2[4+MAXNF]  = 0;  fn1[4+MAXNF]  = 0;

  //set u and d to 0
  //  fn2[-2+MAXNF] = 0; fn1[-2+MAXNF] = 0;
  //  fn2[-1+MAXNF] = 0; fn1[-1+MAXNF] = 0;
  //  fn2[1+MAXNF]  = 0; fn1[1+MAXNF]  = 0;
  //  fn2[2+MAXNF]  = 0; fn1[2+MAXNF]  = 0;
}

//Retrieve PDFs at the starting scale (muf)
void pdfevol::retrievemuf(int i, int sign)
{
  // i is the index of the complex mellin moment in the z-space for the gaussian quadrature used for the mellin inversion

  //N flavour dependence
  int nf = resconst::NF;

  //XP[i] are moments of PDFs at the starting scale (factorisation scale)
  complex <double> fx[11];
  fx[0+MAXNF] = GLP[i];
  fx[1+MAXNF] = UVP[i] + USP[i];
  fx[-1+MAXNF] = USP[i];
  fx[2+MAXNF] = DVP[i] + DSP[i];
  fx[-2+MAXNF] = DSP[i];
  fx[3+MAXNF] = SSP[i];
  fx[-3+MAXNF] = SSP[i];
  if (nf >= 4)
    {
      fx[4+MAXNF] = CHP[i];
      fx[-4+MAXNF] = CHP[i];
    }
  else
    {
      fx[4+MAXNF] = 0.;
      fx[-4+MAXNF] = 0.;
    }
  if (nf >= 5)
    {
      fx[5+MAXNF] = BOP[i];
      fx[-5+MAXNF] = BOP[i];
    }
  else
    {
      fx[5+MAXNF] = 0.;
      fx[-5+MAXNF] = 0.;
    }
  
  storemoments(i, fx);
  retrieve1d(i, sign);
  //  cout << i << "  " << GLP[i] << "  " << fx[0+MAXNF] << "  " << fn1[MAXNF] << "  " << fn2[MAXNF] << endl;
  return;
}
void pdfevol::truncate()
{
  //Calculate truncated moments
  double x1 = phasespace::m/opts.sroot*exp(phasespace::ymin);
  double x2 = phasespace::m/opts.sroot*exp(-phasespace::ymax);

  //double x1 = 1e-8;//pow(phasespace::m/opts.sroot,2);
  //double x2 = 1e-8;//pow(phasespace::m/opts.sroot,2);

  double lx1 = log(x1);
  double lx2 = log(x2);
  
  //truncated moments (output)
  complex <double> fx1[mellinint::mdim][2*MAXNF+1] = {0.};
  complex <double> fx2[mellinint::mdim][2*MAXNF+1] = {0.};

  //cache x^(N) values
  complex <double> x1n[mellinint::mdim];
  complex <double> x2n[mellinint::mdim];
  for (int n = 0; n < mellinint::mdim; n++)
    {
      x1n[n] = pow(x1,mellinint::Np[n]);
      x2n[n] = pow(x2,mellinint::Np[n]);
    }

  //Normalisation times Jacobian
  complex <double> facp = mellinint::CCp/2./M_PI/complex <double>(0.,1);
  complex <double> facm = mellinint::CCm/2./M_PI/complex <double>(0.,1);
  
  //original moments times prefactor and weight
  complex <double> fm1p[mellinint::mdim][2*MAXNF+1];
  complex <double> fm2p[mellinint::mdim][2*MAXNF+1];
  complex <double> fm1m[mellinint::mdim][2*MAXNF+1];
  complex <double> fm2m[mellinint::mdim][2*MAXNF+1];
  for (int m = 0; m < mellinint::mdim; m++)
    for (int f = 0; f < 2*MAXNF+1; f++)
      {
	fm1p[m][f] = facp * cx(creno_.cfx1_[m][f] ) * mellinint::wn[m];
	fm2p[m][f] = facp * cx(creno_.cfx2p_[m][f]) * mellinint::wn[m];
	fm1m[m][f] = facm * conj(cx(creno_.cfx1_[m][f]) ) * mellinint::wn[m];
	fm2m[m][f] = facm * conj(cx(creno_.cfx2p_[m][f])) * mellinint::wn[m];
      }

  //cache factor (1-x^(N-M))/(N-M) which limit is ln(x) when N-M -> 0
  complex <double> llx1p[mellinint::mdim][mellinint::mdim];
  complex <double> llx2p[mellinint::mdim][mellinint::mdim];
  complex <double> llx1m[mellinint::mdim][mellinint::mdim];
  complex <double> llx2m[mellinint::mdim][mellinint::mdim];
  for (int n = 0; n < mellinint::mdim; n++)
    for (int m = 0; m < mellinint::mdim; m++)
      {
	llx1p[n][m] = (1.-x1n[n]/x1n[m])/(mellinint::Np[n]-mellinint::Np[m]);
	llx2p[n][m] = (1.-x2n[n]/x2n[m])/(mellinint::Np[n]-mellinint::Np[m]);
	llx1m[n][m] = (1.-x1n[n]/conj(x1n[m]))/(mellinint::Np[n]-mellinint::Nm[m]);
	llx2m[n][m] = (1.-x2n[n]/conj(x2n[m]))/(mellinint::Np[n]-mellinint::Nm[m]);
      }

  //overwrite divergent diagonal part
  for (int n = 0; n < mellinint::mdim; n++)
    {
      llx1p[n][n] = -lx1;
      llx2p[n][n] = -lx2;
    }
  
  for (int n = 0; n < mellinint::mdim; n++)
    {
      //positive branch
      for (int m = 0; m < mellinint::mdim; m++)
	for (int f = 0; f < 2*MAXNF+1; f++)
	  {
	    /*
	    if (m == n)
	      {
		fx1[n][f] += fm1p[m][f] * (-lx1);
		fx2[n][f] += fm2p[m][f] * (-lx2);
	      }
	    else
	      {
		fx1[n][f] += fm1p[m][f] * (1.-x1n[n]/x1n[m])/(mellinint::Np[n]-mellinint::Np[m]);
		fx2[n][f] += fm2p[m][f] * (1.-x2n[n]/x2n[m])/(mellinint::Np[n]-mellinint::Np[m]);
	      }
	    */
	    //	    fx1[n][f] += fm1p[m][f]*llx1p[n][m]; 
	    //	    fx2[n][f] += fm2p[m][f]*llx2p[n][m];
	    fx1[n][f] += fm1p[m][f]*llx1p[n][m] - fm1m[m][f]*llx1m[n][m]; 
	    fx2[n][f] += fm2p[m][f]*llx2p[n][m] - fm2m[m][f]*llx2m[n][m];
	  }

      /*
      //negative branch
      for (int m = 0; m < mellinint::mdim; m++)
	for (int f = 0; f < 2*MAXNF+1; f++)
	  {
	    //	    fx1[n][f] -= fm1m[m][f]*(1.-x1n[n]/conj(x1n[m]))/(mellinint::Np[n]-mellinint::Nm[m]);
	    //	    fx2[n][f] -= fm2m[m][f]*(1.-x2n[n]/conj(x2n[m]))/(mellinint::Np[n]-mellinint::Nm[m]);
	    fx1[n][f] -= fm1m[m][f]*llx1m[n][m];
	    fx2[n][f] -= fm2m[m][f]*llx2m[n][m];
	  }
      */

  
      //cout << "truncated " << n << fx1[5] << endl;
    }

  //replace moments
  for (int n = 0; n < mellinint::mdim; n++)
    storemoments(n, fx1[n]);
  //storemoments(n, fx1[n], fx2[n]);
}


void pdfevol::uppertruncate()
{
  //Calculate truncated moments
  double xmin1 = pow(phasespace::m/opts.sroot,2);
  double xmin2 = pow(phasespace::m/opts.sroot,2);

  double xmax1 = phasespace::m/opts.sroot*exp(phasespace::ymax);
  double xmax2 = phasespace::m/opts.sroot*exp(-phasespace::ymin);

  double lx1 = log(xmin1/xmax1);
  double lx2 = log(xmin2/xmax2);
  
  //truncated moments (output)
  complex <double> fx1[mellinint::mdim][2*MAXNF+1] = {0.};
  complex <double> fx2[mellinint::mdim][2*MAXNF+1] = {0.};

  //cache x^(N) values
  complex <double> xmin1n[mellinint::mdim];
  complex <double> xmin2n[mellinint::mdim];
  complex <double> xmax1n[mellinint::mdim];
  complex <double> xmax2n[mellinint::mdim];
  for (int n = 0; n < mellinint::mdim; n++)
    {
      xmin1n[n] = pow(xmin1,mellinint::Np[n]);
      xmin2n[n] = pow(xmin2,mellinint::Np[n]);
      xmax1n[n] = pow(xmax1,mellinint::Np[n]);
      xmax2n[n] = pow(xmax2,mellinint::Np[n]);
    }

  //Normalisation times Jacobian
  complex <double> facp = mellinint::CCp/2./M_PI/complex <double>(0.,1);
  complex <double> facm = mellinint::CCm/2./M_PI/complex <double>(0.,1);
  
  //original moments times prefactor and weight
  complex <double> fm1p[mellinint::mdim][2*MAXNF+1];
  complex <double> fm2p[mellinint::mdim][2*MAXNF+1];
  complex <double> fm1m[mellinint::mdim][2*MAXNF+1];
  complex <double> fm2m[mellinint::mdim][2*MAXNF+1];
  for (int m = 0; m < mellinint::mdim; m++)
    for (int f = 0; f < 2*MAXNF+1; f++)
      {
	fm1p[m][f] = facp * cx(creno_.cfx1_[m][f] ) * mellinint::wn[m];
	fm2p[m][f] = facp * cx(creno_.cfx2p_[m][f]) * mellinint::wn[m];
	fm1m[m][f] = facm * conj(cx(creno_.cfx1_[m][f]) ) * mellinint::wn[m];
	fm2m[m][f] = facm * conj(cx(creno_.cfx2p_[m][f])) * mellinint::wn[m];
      }

  //cache factor (xmax^(N-M)-xmin^(N-M))/(N-M) which limit is ln(xmin/xmax) when N-M -> 0
  complex <double> llx1p[mellinint::mdim][mellinint::mdim];
  complex <double> llx2p[mellinint::mdim][mellinint::mdim];
  complex <double> llx1m[mellinint::mdim][mellinint::mdim];
  complex <double> llx2m[mellinint::mdim][mellinint::mdim];
  for (int n = 0; n < mellinint::mdim; n++)
    for (int m = 0; m < mellinint::mdim; m++)
      {
	llx1p[n][m] = (xmax1n[n]/xmax1n[m]-xmin1n[n]/xmin1n[m])/(mellinint::Np[n]-mellinint::Np[m]);
	llx2p[n][m] = (xmax2n[n]/xmax2n[m]-xmin2n[n]/xmin2n[m])/(mellinint::Np[n]-mellinint::Np[m]);
	llx1m[n][m] = (xmax1n[n]/conj(xmax1n[m])-xmin1n[n]/conj(xmin1n[m]))/(mellinint::Np[n]-mellinint::Nm[m]);
	llx2m[n][m] = (xmax2n[n]/conj(xmax2n[m])-xmin2n[n]/conj(xmin2n[m]))/(mellinint::Np[n]-mellinint::Nm[m]);
      }

  //overwrite divergent diagonal part
  for (int n = 0; n < mellinint::mdim; n++)
    {
      llx1p[n][n] = -lx1;
      llx2p[n][n] = -lx2;
    }
  
  for (int n = 0; n < mellinint::mdim; n++)
    {
      //positive branch
      for (int m = 0; m < mellinint::mdim; m++)
	for (int f = 0; f < 2*MAXNF+1; f++)
	  {
	    /*
	    if (m == n)
	      {
		fx1[n][f] += fm1p[m][f] * (-lx1);
		fx2[n][f] += fm2p[m][f] * (-lx2);
	      }
	    else
	      {
		fx1[n][f] += fm1p[m][f] * (1.-x1n[n]/x1n[m])/(mellinint::Np[n]-mellinint::Np[m]);
		fx2[n][f] += fm2p[m][f] * (1.-x2n[n]/x2n[m])/(mellinint::Np[n]-mellinint::Np[m]);
	      }
	    */
	    //	    fx1[n][f] += fm1p[m][f]*llx1p[n][m]; 
	    //	    fx2[n][f] += fm2p[m][f]*llx2p[n][m];
	    fx1[n][f] += fm1p[m][f]*llx1p[n][m] - fm1m[m][f]*llx1m[n][m]; 
	    fx2[n][f] += fm2p[m][f]*llx2p[n][m] - fm2m[m][f]*llx2m[n][m];
	  }

      /*
      //negative branch
      for (int m = 0; m < mellinint::mdim; m++)
	for (int f = 0; f < 2*MAXNF+1; f++)
	  {
	    //	    fx1[n][f] -= fm1m[m][f]*(1.-x1n[n]/conj(x1n[m]))/(mellinint::Np[n]-mellinint::Nm[m]);
	    //	    fx2[n][f] -= fm2m[m][f]*(1.-x2n[n]/conj(x2n[m]))/(mellinint::Np[n]-mellinint::Nm[m]);
	    fx1[n][f] -= fm1m[m][f]*llx1m[n][m];
	    fx2[n][f] -= fm2m[m][f]*llx2m[n][m];
	  }
      */

  
      //cout << "truncated " << n << fx1[5] << endl;
    }

  //replace moments
  for (int n = 0; n < mellinint::mdim; n++)
    storemoments(n, fx1[n]);
    //storemoments(n, fx1[n], fx2[n]);
}

/*
void pdfevol::invert(double x, double fx[2*MAXNF+1])
{
  fx = {0.};

  complex <double> fm[2*MAXNF+1];
  complex <double> fm2[2*MAXNF+1];

  double lx = log(x);
  
  //positive branch
  for (int m = 0; m < mellinint::mdim; m++)
    {
      complex <double> cexp = mellinint::CCp/2./M_PI/complex <double>(0.,1) * exp(-mellinint::Np[m] * lx);
      for (int f = 0; f < 2*MAXNF+1; f++)
	{
	  fm[f] = cx(creno_.cfx1_[m][f]);
	  fx[f] += cexp * fm[f] * mellinint::wn[m];
	}

  //negative branch
  for (int m = 0; m < mellinint::mdim; m++)
    {
      complex <double> cexm = mellinint::CCm/2./M_PI/complex <double>(0.,1) * exp(-mellinint::Nm[m] * lx);
      for (int f = 0; f < 2*MAXNF+1; f++)
	{
	  fm[f] = cx(creno_.cfx2m_[m][f]);
	fx1[f] -= mellinint::CCm/2./M_PI/complex <double>(0.,1) * conj(fm1[f])*(1.-pow(x1,mellinint::Np[n]-mellinint::Nm[m]))/(mellinint::Np[n]-mellinint::Nm[m]) * mellinint::wn[m];
	fx2[f] -= mellinint::CCm/2./M_PI/complex <double>(0.,1) * conj(fm2[f])*(1.-pow(x2,mellinint::Np[n]-mellinint::Nm[m]))/(mellinint::Np[n]-mellinint::Nm[m]) * mellinint::wn[m];
      }
  
  cout << "truncated " << n << fx1[5] << endl;


}
*/
