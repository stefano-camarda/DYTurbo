#include "pdfevol.h"
#include "interface.h"
#include "mesq.h"
#include "anomalous.h"
#include "resconst.h"

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
  cout << "Initialise PDF moments with numerical integration (C++) " << endl;

  fcomplex uval,dval,usea,dsea,splus,ssea,glu,charm,bot;
  for (int k = 0; k < mellinint::mdim; k++)
    {
      int hadron = 1; //opts.ih1;
      fcomplex XN = fcx(mellinint::Np[k]); //compute positive branch only, the negative branch is obtained by complex conjugation
      double facscale = opts.kmufac*opts.rmass;
      pdfmoments_(hadron,facscale,XN,uval,dval,usea,dsea,splus,ssea,glu,charm,bot);
      // cout << "moment " << k << "  " << cx(XN) << "  ";
      // cout << "uval  " << uval << endl;
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

  cout << "End PDF moments initialization" << endl;
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
  complex <double> alpq, ALPr;

  //N flavour dependence
  int nf = resconst::NF;

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
  double facscale = fabs(bstarscale);//better qbstar here? actually better qbstar * resint::mures / sqrt(pow(qbstar,2) + resint::mures2);
  //double facscale = fabs(opts.muf);
  fcomplex XN = fcx(mellinint::Np[i]);
  pdfmoments_(hadron,facscale,XN,uval,dval,usea,dsea,s,sbar,glu,charm,bot);

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
