#include "pdfevol.h"
#include "interface.h"
#include "mesq.h"
#include "anomalous.h"
#include <iostream>

//PDFs mellin moments at the factorisation scale
//Beam 1
complex <double> *pdfevol::UVP1;
complex <double> *pdfevol::DVP1;
complex <double> *pdfevol::USP1;
complex <double> *pdfevol::DSP1;
complex <double> *pdfevol::SSP1;
complex <double> *pdfevol::GLP1;
complex <double> *pdfevol::CHP1;
complex <double> *pdfevol::BOP1;
complex <double> *pdfevol::UVM1;
complex <double> *pdfevol::DVM1;
complex <double> *pdfevol::USM1;
complex <double> *pdfevol::DSM1;
complex <double> *pdfevol::SSM1;
complex <double> *pdfevol::GLM1;
complex <double> *pdfevol::CHM1;
complex <double> *pdfevol::BOM1;

//Beam 2
complex <double> *pdfevol::UVP2;
complex <double> *pdfevol::DVP2;
complex <double> *pdfevol::USP2;
complex <double> *pdfevol::DSP2;
complex <double> *pdfevol::SSP2;
complex <double> *pdfevol::GLP2;
complex <double> *pdfevol::CHP2;
complex <double> *pdfevol::BOP2;
complex <double> *pdfevol::UVM2;
complex <double> *pdfevol::DVM2;
complex <double> *pdfevol::USM2;
complex <double> *pdfevol::DSM2;
complex <double> *pdfevol::SSM2;
complex <double> *pdfevol::GLM2;
complex <double> *pdfevol::CHM2;
complex <double> *pdfevol::BOM2;

complex <double> *pdfevol::fn1;
complex <double> *pdfevol::fn2;

//scales
complex <double> pdfevol::XL;
complex <double> pdfevol::XL1;
complex <double> pdfevol::SALP;

complex <double> pdfevol::alpr;


//fortran interface
void pdfevol_(int& i1, int& i2, int& sign)
{
  pdfevol::evolve(i1-1, i2-1, sign-1);
};

void pdfevol::init()
{
  fn1 = new complex <double>[11];
  fn2 = new complex <double>[11];

  UVP1 = new complex <double>[mellinint::mdim];
  DVP1 = new complex <double>[mellinint::mdim];
  USP1 = new complex <double>[mellinint::mdim];
  DSP1 = new complex <double>[mellinint::mdim];
  SSP1 = new complex <double>[mellinint::mdim];
  GLP1 = new complex <double>[mellinint::mdim];
  CHP1 = new complex <double>[mellinint::mdim];
  BOP1 = new complex <double>[mellinint::mdim];
  UVM1 = new complex <double>[mellinint::mdim];
  DVM1 = new complex <double>[mellinint::mdim];
  USM1 = new complex <double>[mellinint::mdim];
  DSM1 = new complex <double>[mellinint::mdim];
  SSM1 = new complex <double>[mellinint::mdim];
  GLM1 = new complex <double>[mellinint::mdim];
  CHM1 = new complex <double>[mellinint::mdim];
  BOM1 = new complex <double>[mellinint::mdim];

  UVP2 = new complex <double>[mellinint::mdim];
  DVP2 = new complex <double>[mellinint::mdim];
  USP2 = new complex <double>[mellinint::mdim];
  DSP2 = new complex <double>[mellinint::mdim];
  SSP2 = new complex <double>[mellinint::mdim];
  GLP2 = new complex <double>[mellinint::mdim];
  CHP2 = new complex <double>[mellinint::mdim];
  BOP2 = new complex <double>[mellinint::mdim];
  UVM2 = new complex <double>[mellinint::mdim];
  DVM2 = new complex <double>[mellinint::mdim];
  USM2 = new complex <double>[mellinint::mdim];
  DSM2 = new complex <double>[mellinint::mdim];
  SSM2 = new complex <double>[mellinint::mdim];
  GLM2 = new complex <double>[mellinint::mdim];
  CHM2 = new complex <double>[mellinint::mdim];
  BOM2 = new complex <double>[mellinint::mdim];

  //     calculate Mellin moments of PDFs
  cout << "Initialise PDF moments with numerical integration (C++) " << endl;

  fcomplex uval,dval,usea,dsea,ssea,glu,charm,bot;
  //  cout << "beam 1 positive" << endl;
  //beam 1
  for (int k = 0; k < mellinint::mdim; k++)
    {
      int hadron = opts.ih1;
      // positive branch
      fcomplex XN = fcx(mellinint::Np[k]);
      pdfmoments_(hadron,XN,uval,dval,usea,dsea,ssea,glu,charm,bot);
      // cout << "moment " << k << "  " << cx(XN) << endl;
      // cout << "uval  " << uval << endl;
      // cout << "dval  " << dval << endl;
      // cout << "usea  " << usea << endl;
      // cout << "dsea  " << dsea << endl;
      // cout << "gluon " << cx(glu) << endl;
      // cout << "charm " << charm << endl;
      // cout << "bottom" << bot << endl;
      UVP1[k] = cx(uval);
      DVP1[k] = cx(dval);
      USP1[k] = cx(usea);
      DSP1[k] = cx(dsea);
      SSP1[k] = cx(ssea);
      GLP1[k] = cx(glu);
      CHP1[k] = cx(charm);
      BOP1[k] = cx(bot);
      // negative branch
      XN = fcx(mellinint::Nm[k]);
      pdfmoments_(hadron,XN,uval,dval,usea,dsea,ssea,glu,charm,bot);
      UVM1[k] = cx(uval);
      DVM1[k] = cx(dval);
      USM1[k] = cx(usea);
      DSM1[k] = cx(dsea);
      SSM1[k] = cx(ssea);
      GLM1[k] = cx(glu);
      CHM1[k] = cx(charm);
      BOM1[k] = cx(bot);
    }

  //beam 2
  for (int k = 0; k < mellinint::mdim; k++)
    {
      int hadron = opts.ih2;
      // positive branch
      fcomplex XN = fcx(mellinint::Np[k]);
      pdfmoments_(hadron,XN,uval,dval,usea,dsea,ssea,glu,charm,bot);
      UVP2[k] = cx(uval);
      DVP2[k] = cx(dval);
      USP2[k] = cx(usea);
      DSP2[k] = cx(dsea);
      SSP2[k] = cx(ssea);
      GLP2[k] = cx(glu);
      CHP2[k] = cx(charm);
      BOP2[k] = cx(bot);
      // negative branch
      XN = fcx(mellinint::Nm[k]);
      pdfmoments_(hadron,XN,uval,dval,usea,dsea,ssea,glu,charm,bot);
      UVM2[k] = cx(uval);
      DVM2[k] = cx(dval);
      USM2[k] = cx(usea);
      DSM2[k] = cx(dsea);
      SSM2[k] = cx(ssea);
      GLM2[k] = cx(glu);
      CHM2[k] = cx(charm);
      BOM2[k] = cx(bot);
    }
  cout << "End PDF moments initialization" << endl;
}

void pdfevol::evolve(int i1, int i2, int sign)
{
  //  cout << i1 << endl;
  //  cout << creno_.cfx1_[i1][5].real << "  " << creno_.cfx1_[i1][5].imag << endl;
  fn1[0] = cx(creno_.cfx1_[i1][0]);
  fn1[1] = cx(creno_.cfx1_[i1][1]);
  fn1[2] = cx(creno_.cfx1_[i1][2]);
  fn1[3] = cx(creno_.cfx1_[i1][3]);
  fn1[4] = cx(creno_.cfx1_[i1][4]);
  fn1[5] = cx(creno_.cfx1_[i1][5]);
  fn1[6] = cx(creno_.cfx1_[i1][6]);
  fn1[7] = cx(creno_.cfx1_[i1][7]);
  fn1[8] = cx(creno_.cfx1_[i1][8]);
  fn1[9] = cx(creno_.cfx1_[i1][9]);
  fn1[10] = cx(creno_.cfx1_[i1][10]);
  if (sign == mesq::positive)
    {
      fn2[0] = cx(creno_.cfx2p_[i2][0]);
      fn2[1] = cx(creno_.cfx2p_[i2][1]);
      fn2[2] = cx(creno_.cfx2p_[i2][2]);
      fn2[3] = cx(creno_.cfx2p_[i2][3]);
      fn2[4] = cx(creno_.cfx2p_[i2][4]);
      fn2[5] = cx(creno_.cfx2p_[i2][5]);
      fn2[6] = cx(creno_.cfx2p_[i2][6]);
      fn2[7] = cx(creno_.cfx2p_[i2][7]);
      fn2[8] = cx(creno_.cfx2p_[i2][8]);
      fn2[9] = cx(creno_.cfx2p_[i2][9]);
      fn2[10] = cx(creno_.cfx2p_[i2][10]);
    }
  else if (sign == mesq::negative)
    {
      fn2[0] = cx(creno_.cfx2m_[i2][0]);
      fn2[1] = cx(creno_.cfx2m_[i2][1]);
      fn2[2] = cx(creno_.cfx2m_[i2][2]);
      fn2[3] = cx(creno_.cfx2m_[i2][3]);
      fn2[4] = cx(creno_.cfx2m_[i2][4]);
      fn2[5] = cx(creno_.cfx2m_[i2][5]);
      fn2[6] = cx(creno_.cfx2m_[i2][6]);
      fn2[7] = cx(creno_.cfx2m_[i2][7]);
      fn2[8] = cx(creno_.cfx2m_[i2][8]);
      fn2[9] = cx(creno_.cfx2m_[i2][9]);
      fn2[10] = cx(creno_.cfx2m_[i2][10]);
    }
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
void pdfevol::evolution (int i, int sign, int beam) //from reno2
{
  complex <double> alpq,ALPr;
  //       COMMON/alphasldata/XL,XL1,SALP,alpq,ALPr

  //N flavour dependence
  int nnf = 5;
  int f = 5;

  //Moments at the factorisation scale
  complex <double> UVI,	DVI, USI, DSI, SSI, GLI, CHI, BOI;

  // **************************************
  //     IFIT is only rapidity dependent, I is the point in the z grid of the gaussian quadrature loop
  //int IFIT
  if (beam == 1)
    if (sign == mesq::positive)
      {
	UVI = UVP1[i];
	DVI = DVP1[i];
	USI = USP1[i];
	DSI = DSP1[i];
	SSI = SSP1[i];
	GLI = GLP1[i];
	CHI = CHP1[i];
	BOI = BOP1[i];
      }
    else
      {
        UVI = UVM1[i];
        DVI = DVM1[i];
        USI = USM1[i];
        DSI = DSM1[i];
        SSI = SSM1[i];
        GLI = GLM1[i];
        CHI = CHM1[i];
        BOI = BOM1[i];
      }
  else if (beam == 2)
    if (sign == mesq::positive)
      {
        UVI = UVP2[i];
        DVI = DVP2[i];
        USI = USP2[i];
        DSI = DSP2[i];
        SSI = SSP2[i];
        GLI = GLP2[i];
        CHI = CHP2[i];
        BOI = BOP2[i];
      }
    else
      {
        UVI = UVM2[i];
        DVI = DVM2[i];
        USI = USM2[i];
        DSI = DSM2[i];
        SSI = SSM2[i];
        GLI = GLM2[i];
        CHI = CHM2[i];
        BOI = BOM2[i];
      }
  else
    {
      cout << "Wrong beam in evolution:  " << beam << endl;
      exit(-1);
    }
             
  // ***********************
  // this part can be precomputed
  complex <double> UVN = UVI;
  complex <double> DVN = DVI;
  complex <double> NS3N = UVI + 2.*USI - DVI - 2.*DSI;
  complex <double> NS8N = UVI + 2.*USI + DVI + 2.*DSI - 4.*SSI;
  complex <double> GLN = GLI;
  complex <double> SIN = UVI + DVI + 2.*USI + 2.*DSI + 2.*SSI + 2.*CHI + 2.*BOI;
  complex <double> NS15N = UVI + DVI + 2.*USI + 2.*DSI + 2.*SSI - 6.*CHI;
  complex <double> NS24N = UVI + DVI + 2.*USI + 2.*DSI + 2.*SSI + 2.*CHI - 8.*BOI;
  complex <double> NS35N = SIN;
 
  if (nnf == 3)
    {
      SIN = UVI + DVI + 2.*USI + 2.*DSI + 2.*SSI;
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
  complex <double> RMIN = anomalous::rmin[anomalous::index(i,sign)];
  complex <double> RPLUS = anomalous::rplus[anomalous::index(i,sign)];
  complex <double> AC  = 1. -AL;
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
  //double q2s = q2/pow(opts.a_param,2);                //resummation scale
  //double alpqf = dyalphas_lhapdf_(sqrt(q2s))/4./M_PI; //alphas at the resummation scale
  //complex <double> alpq = alpqf * alphasl(scale2);              //alphas at the resummation scale times alphas at 1/b
  //complex <double> alpr= alpq * 1 *(opts.order-1);
  
  UVN  = UVN  * ENS * (1.+  alpr * XL1 * RMIN);
  DVN  = DVN  * ENS * (1.+  alpr * XL1 * RMIN);
  NS3N = NS3N * ENS * (1.+  alpr * XL1 * RPLUS);
  NS8N = NS8N * ENS * (1.+  alpr * XL1 * RPLUS);

  SIN = EM * ((AL + alpr * (RMMQQ * XL1 + RMPQQ * (EPM-XL)))* SG
	      + (BE + alpr * (RMMQG * XL1 + RMPQG * (EPM-XL))) * GL)
    + EP * ((AC + alpr * (RPPQQ * XL1 + RPMQQ * (EMP-XL))) * SG
	    +(-BE + alpr * (RPPQG * XL1 + RPMQG * (EMP-XL))) * GL);
  GLN = EM * ((AB + alpr * (RMMGQ * XL1 + RMPGQ * (EPM-XL))) * SG
	      + (AC + alpr * (RMMGG * XL1 + RMPGG * (EPM-XL))) * GL)
    + EP *((-AB + alpr * (RPPGQ * XL1 + RPMGQ * (EMP-XL))) * SG
	   + (AL + alpr * (RPPGG * XL1 + RPMGG * (EMP-XL))) * GL);
  
  NS15N = NS15N * ENS * (1.+  alpr * XL1 * RPLUS);
  NS24N = NS24N * ENS * (1.+  alpr * XL1 * RPLUS);
  NS35N = SIN;

  //...  FLAVOUR DECOMPOSITION OF THE QUARK SEA :
  complex <double> SSN = (10.* SIN + 2.* NS35N + 3.* NS24N + 5.* NS15N - 20.* NS8N) / 120.;
  complex <double> DSN = (10.* SIN + 2.* NS35N + 3.* NS24N + 5.* NS15N + 10.* NS8N - 30.* NS3N - 60.* DVN) / 120.;
  complex <double> USN = (10.* SIN + 2.* NS35N + 3.* NS24N + 5.* NS15N + 10.* NS8N + 30.* NS3N - 60.* UVN) / 120.;
  complex <double> CHN = (UVN + DVN + 2.*USN + 2.*DSN + 2.*SSN - NS15N)/6.;
  complex <double> BON = (UVN + DVN + 2.*USN + 2.*DSN + 2.*SSN + 2.*CHN - NS24N)/8.;
      
  if (nnf == 3) //GRV
    {
      SSN= (20.* SIN - 20.* NS8N)/120.;
      DSN = (20.* SIN + 10.* NS8N - 30.* NS3N - 60.* DVN) / 120.;
      USN = (20.* SIN + 10.* NS8N + 30.* NS3N - 60.* UVN) / 120.;
      CHN=0.;
      BON=0.;
    }

  // **************************************

  // ...  OUTPUT  
  //   U   UB   D   DB   S   C    B    G 
  //   1    2   3   4    5   6    7    8 
  // ************************************

  complex <double> fx[11];
  fx[0+5] = GLN;
  fx[1+5] = UVN + USN;
  fx[-1+5] = USN;
  fx[2+5] = DVN + DSN;
  fx[-2+5] = DSN;
  fx[3+5] = SSN;
  fx[-3+5] = SSN;
  if (nnf >= 4)
    {
      fx[4+5] = CHN;
      fx[-4+5] = CHN;
    }
  else
    {
      fx[4+5] = 0.;
      fx[-4+5] = 0.;
    }
  
  if (nnf >= 5)
    {
      fx[5+5] = BON;
      fx[-5+5] = BON;
    }
  else
    {
      fx[5+5] = 0.;
      fx[-5+5] = 0.;
    }

  //convert proton into antiproton (need to switch only u and d, since the sea is assumed symmetric)
  if ((beam == 1 && opts.ih1 == -1) || (beam == 2 && opts.ih2 == -1))
    {
      complex <double> utemp = fx[1+5];
      fx[1+5] = fx[-1+5];
      fx[-1+5] = utemp;
      complex <double> dtemp = fx[2+5];
      fx[2+5] = fx[-2+5];
      fx[-2+5] = dtemp;
    }

  //Save the evolved PDFs into the fortran common block (can actually use a C++ data format since it is only accessed in C++)
  if (beam == 1)
    {
      creno_.cfx1_[i][0]  = fcx(fx[0]);
      creno_.cfx1_[i][1]  = fcx(fx[1]);
      creno_.cfx1_[i][2]  = fcx(fx[2]);
      creno_.cfx1_[i][3]  = fcx(fx[3]);
      creno_.cfx1_[i][4]  = fcx(fx[4]);
      creno_.cfx1_[i][5]  = fcx(fx[5]);
      creno_.cfx1_[i][6]  = fcx(fx[6]);
      creno_.cfx1_[i][7]  = fcx(fx[7]);
      creno_.cfx1_[i][8]  = fcx(fx[8]);
      creno_.cfx1_[i][9]  = fcx(fx[9]);
      creno_.cfx1_[i][10] = fcx(fx[10]);
    }
  else if (beam == 2)
    if (sign == mesq::positive)
      {
      creno_.cfx2p_[i][0]  = fcx(fx[0] );
      creno_.cfx2p_[i][1]  = fcx(fx[1] );
      creno_.cfx2p_[i][2]  = fcx(fx[2] );
      creno_.cfx2p_[i][3]  = fcx(fx[3] );
      creno_.cfx2p_[i][4]  = fcx(fx[4] );
      creno_.cfx2p_[i][5]  = fcx(fx[5] );
      creno_.cfx2p_[i][6]  = fcx(fx[6] );
      creno_.cfx2p_[i][7]  = fcx(fx[7] );
      creno_.cfx2p_[i][8]  = fcx(fx[8] );
      creno_.cfx2p_[i][9]  = fcx(fx[9] );
      creno_.cfx2p_[i][10] = fcx(fx[10]);
      }
    else if (sign == mesq::negative)
      {
	creno_.cfx2m_[i][0]  = fcx(fx[0] );
	creno_.cfx2m_[i][1]  = fcx(fx[1] );
	creno_.cfx2m_[i][2]  = fcx(fx[2] );
	creno_.cfx2m_[i][3]  = fcx(fx[3] );
	creno_.cfx2m_[i][4]  = fcx(fx[4] );
	creno_.cfx2m_[i][5]  = fcx(fx[5] );
	creno_.cfx2m_[i][6]  = fcx(fx[6] );
	creno_.cfx2m_[i][7]  = fcx(fx[7] );
	creno_.cfx2m_[i][8]  = fcx(fx[8] );
	creno_.cfx2m_[i][9]  = fcx(fx[9] );
	creno_.cfx2m_[i][10] = fcx(fx[10]);
      }
  // **************************************
}
