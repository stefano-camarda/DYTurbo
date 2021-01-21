#include "ccoeff.h"
#include "dynnlo_interface.h"
#include "dyres_interface.h"
#include "anomalous.h"
#include "mellinint.h"
#include "mesq.h"
#include "settings.h"
#include "gaussrules.h"
#include "resconst.h"
#include "constants.h"
#include "icoeff.h"
#include "phasespace.h"
//#include "hs.h"
#include "psi.h"
#include "cmom.h"
#include "melfun.h"
#include <iostream>
#include <iomanip>


//fotran interface for dymellinh2.f
extern "C"
{
  fcomplex fun1_(fcomplex &xn)  {return fcx(fun1(cx(xn)));};
  fcomplex fun2_(fcomplex &xn)  {return fcx(fun2(cx(xn)));};
  fcomplex fun3_(fcomplex &xn)  {return fcx(fun3(cx(xn)));};
  fcomplex fun4_(fcomplex &xn)  {return fcx(fun4(cx(xn)));};
  fcomplex fun5_(fcomplex &xn)  {return fcx(fun5(cx(xn)));};
  fcomplex fun6_(fcomplex &xn)  {return fcx(fun6(cx(xn)));};
  fcomplex fun7_(fcomplex &xn)  {return fcx(fun7(cx(xn)));};
  fcomplex fun8_(fcomplex &xn)  {return fcx(fun8(cx(xn)));};
  fcomplex fun9_(fcomplex &xn)  {return fcx(fun9(cx(xn)));};
  fcomplex fun10_(fcomplex &xn) {return fcx(fun10(cx(xn)));};
  fcomplex fun11_(fcomplex &xn) {return fcx(fun11(cx(xn)));};
  fcomplex fun12_(fcomplex &xn) {return fcx(fun12(cx(xn)));};
  fcomplex fun13_(fcomplex &xn) {return fcx(fun13(cx(xn)));};
  fcomplex fun14_(fcomplex &xn) {return fcx(fun14(cx(xn)));};
  fcomplex fun15_(fcomplex &xn) {return fcx(fun15(cx(xn)));};
  fcomplex fun16_(fcomplex &xn) {return fcx(fun16(cx(xn)));};

  fcomplex dypsi0_(fcomplex &xn)  {return fcx(psi0(cx(xn)));};
  fcomplex dypsi1_(fcomplex &xn)  {return fcx(psi1(cx(xn)));};
  fcomplex dypsi2_(fcomplex &xn)  {return fcx(psi2(cx(xn)));};
  fcomplex dypsi3_(fcomplex &xn)  {return fcx(psi3(cx(xn)));};
  fcomplex dybe0_(fcomplex &xn)  {return fcx(be0(cx(xn)));};
  fcomplex dybe1_(fcomplex &xn)  {return fcx(be1(cx(xn)));};
  fcomplex dybe2_(fcomplex &xn)  {return fcx(be2(cx(xn)));};
  fcomplex dybe3_(fcomplex &xn)  {return fcx(be3(cx(xn)));};
}

using namespace std;
using namespace resconst;
using namespace constants;

// overload int - complex operators to allow arithmetics
// std::complex<double> operator+( const int& lhs, const std::complex<double>& rhs)
// {
//   return double (lhs) + rhs;
// }
// std::complex<double> operator-( const int& lhs, const std::complex<double>& rhs)
// {
//   return double (lhs) - rhs;
// }
// std::complex<double> operator*( const int& lhs, const std::complex<double>& rhs)
// {
//   return double (lhs) * rhs;
// }
// std::complex<double> operator+( const std::complex<double>& lhs, const int& rhs)
// {
//   return lhs + double(rhs);
// }
// std::complex<double> operator-( const std::complex<double>& lhs, const int& rhs)
// {
//   return lhs - double(rhs);
// }
// std::complex<double> operator*( const std::complex<double>& lhs, const int& rhs )
// {
//   return lhs * double(rhs);
// }


double ccoeff::C1qq_delta;
double ccoeff::C2qq_delta;
double ccoeff::C3qq_delta;
double ccoeff::C3qq_delta_NFV;

complex <double> *ccoeff::C1qg;
complex <double> *ccoeff::C1qq;
complex <double> *ccoeff::C1qqb;
complex <double> *ccoeff::C1qqp;
complex <double> *ccoeff::C1qqbp;

complex <double> *ccoeff::C2qg;
complex <double> *ccoeff::C2qq;
complex <double> *ccoeff::C2qqb;
complex <double> *ccoeff::C2qqp;
complex <double> *ccoeff::C2qqbp;

complex <double> *ccoeff::C3qg;
complex <double> *ccoeff::C3qq;
complex <double> *ccoeff::C3qqb;
complex <double> *ccoeff::C3qqp;
complex <double> *ccoeff::C3qqbp;

complex <double> *ccoeff::C1qg_1;
complex <double> *ccoeff::C1qq_1;
complex <double> *ccoeff::C1qqb_1;
complex <double> *ccoeff::C1qqp_1;
complex <double> *ccoeff::C1qqbp_1;
complex <double> *ccoeff::C2qg_1;
complex <double> *ccoeff::C2qq_1;
complex <double> *ccoeff::C2qqb_1;
complex <double> *ccoeff::C2qqp_1;
complex <double> *ccoeff::C2qqbp_1;
complex <double> *ccoeff::C3qg_1;
complex <double> *ccoeff::C3qq_1;
complex <double> *ccoeff::C3qqb_1;
complex <double> *ccoeff::C3qqp_1;
complex <double> *ccoeff::C3qqbp_1;

complex <double> *ccoeff::C1qg_2;
complex <double> *ccoeff::C1qq_2;
complex <double> *ccoeff::C1qqb_2;
complex <double> *ccoeff::C1qqp_2;
complex <double> *ccoeff::C1qqbp_2;
complex <double> *ccoeff::C2qg_2;
complex <double> *ccoeff::C2qq_2;
complex <double> *ccoeff::C2qqb_2;
complex <double> *ccoeff::C2qqp_2;
complex <double> *ccoeff::C2qqbp_2;
complex <double> *ccoeff::C3qg_2;
complex <double> *ccoeff::C3qq_2;
complex <double> *ccoeff::C3qqb_2;
complex <double> *ccoeff::C3qqp_2;
complex <double> *ccoeff::C3qqbp_2;

//const int rule = 500;
//double t[rule];              //gauss nodes
//double fac[rule];            //overall factor
//complex <double> kern[rule]; //kernel of the Mellin transform

void ccoeff::allocate()
{
  if (opts.mellin1d)
    {
      C1qg = new complex <double>[mellinint::mdim*2];
      C1qq = new complex <double>[mellinint::mdim*2];
      C1qqb = new complex <double>[mellinint::mdim*2];
      C1qqp = new complex <double>[mellinint::mdim*2];
      C1qqbp = new complex <double>[mellinint::mdim*2];

      C2qg = new complex <double>[mellinint::mdim*2];
      C2qq = new complex <double>[mellinint::mdim*2];
      C2qqb = new complex <double>[mellinint::mdim*2];
      C2qqp = new complex <double>[mellinint::mdim*2];
      C2qqbp = new complex <double>[mellinint::mdim*2];

      C3qg = new complex <double>[mellinint::mdim*2];
      C3qq = new complex <double>[mellinint::mdim*2];
      C3qqb = new complex <double>[mellinint::mdim*2];
      C3qqp = new complex <double>[mellinint::mdim*2];
      C3qqbp = new complex <double>[mellinint::mdim*2];

      //fill(C1qg,   C1qg+2*mellinint::mdim, 0.);
      //fill(C1qq,   C1qq+2*mellinint::mdim, 0.);
      //fill(C1qqb,  C1qqb+2*mellinint::mdim, 0.);
      //fill(C1qqp,  C1qqp+2*mellinint::mdim, 0.);
      //fill(C1qqbp, C1qqbp+2*mellinint::mdim, 0.);
      //
      //fill(C2qg,   C2qg+2*mellinint::mdim, 0.);
      //fill(C2qq,   C2qq+2*mellinint::mdim, 0.);
      //fill(C2qqb,  C2qqb+2*mellinint::mdim, 0.);
      //fill(C2qqp,  C2qqp+2*mellinint::mdim, 0.);
      //fill(C2qqbp, C2qqbp+2*mellinint::mdim, 0.);
      //
      //fill(C3qg,   C2qg+2*mellinint::mdim, 0.);
      //fill(C3qq,   C2qq+2*mellinint::mdim, 0.);
      //fill(C3qqb,  C2qqb+2*mellinint::mdim, 0.);
      //fill(C3qqp,  C2qqp+2*mellinint::mdim, 0.);
      //fill(C3qqbp, C2qqbp+2*mellinint::mdim, 0.);
    }
  else
    {
      C1qg_1 = new complex <double>[mellinint::mdim*2];
      C1qq_1 = new complex <double>[mellinint::mdim*2];
      C1qqb_1 = new complex <double>[mellinint::mdim*2];
      C1qqp_1 = new complex <double>[mellinint::mdim*2];
      C1qqbp_1 = new complex <double>[mellinint::mdim*2];
      C2qg_1 = new complex <double>[mellinint::mdim*2];
      C2qq_1 = new complex <double>[mellinint::mdim*2];
      C2qqb_1 = new complex <double>[mellinint::mdim*2];
      C2qqp_1 = new complex <double>[mellinint::mdim*2];
      C2qqbp_1 = new complex <double>[mellinint::mdim*2];
      C3qg_1 = new complex <double>[mellinint::mdim*2];
      C3qq_1 = new complex <double>[mellinint::mdim*2];
      C3qqb_1 = new complex <double>[mellinint::mdim*2];
      C3qqp_1 = new complex <double>[mellinint::mdim*2];
      C3qqbp_1 = new complex <double>[mellinint::mdim*2];
      C1qg_2 = new complex <double>[mellinint::mdim*2];
      C1qq_2 = new complex <double>[mellinint::mdim*2];
      C1qqb_2 = new complex <double>[mellinint::mdim*2];
      C1qqp_2 = new complex <double>[mellinint::mdim*2];
      C1qqbp_2 = new complex <double>[mellinint::mdim*2];
      C2qg_2 = new complex <double>[mellinint::mdim*2];
      C2qq_2 = new complex <double>[mellinint::mdim*2];
      C2qqb_2 = new complex <double>[mellinint::mdim*2];
      C2qqp_2 = new complex <double>[mellinint::mdim*2];
      C2qqbp_2 = new complex <double>[mellinint::mdim*2];
      C3qg_2 = new complex <double>[mellinint::mdim*2];
      C3qq_2 = new complex <double>[mellinint::mdim*2];
      C3qqb_2 = new complex <double>[mellinint::mdim*2];
      C3qqp_2 = new complex <double>[mellinint::mdim*2];
      C3qqbp_2 = new complex <double>[mellinint::mdim*2];
    }
}

void ccoeff::init()
{
  //Evaluate N-independent part
  delta();
  icoeff::allocate();
  icoeff::init();

  if (opts.melup <= 1)
    allocate();
  if (opts.melup == 0)
    if (opts.mellin1d)
      calc1d();
    else
      calc2d();
}

void ccoeff::delta()
{
  //delta pieces -> N-independent part (proportional to delta(1-z))
  int NF2 = NF*NF; //Number of flavours from resconst
  double NFV;
  //NFV is the charge weighted sum of the quark flavours (in the case of purely electromagnetic interactions we have NFgamma=1/eq*Sum_{i=,nf} ei).
  // --> ??? Check this, I do not understand the formula
  //NFV = Sum_{i=1,nf} eq = +1/3 -2/3 +1/3 -2/3 -2/3
  //double NFV = -4./3.;
  //NFV = Sum_{i=1,nf} |eq| = +1/3 +2/3 +1/3 +2/3 +2/3
  NFV = 0.;//8./3.;

  double NFV_gamma[2*MAXNF+1];
  double NFV_Z[2*MAXNF+1];

  //no quark box contributions for W (--> to be checked)
  if (opts.nproc == 1 || opts.nproc == 2)
    NFV = 0;
  
  double sum_gamma = 0;
  double sum_Z = 0;
  for (int f = 1; f <= MAXNF; f++)
    {
      sum_gamma += ewcharge_.Q_[MAXNF+f]   + ewcharge_.Q_[MAXNF-f];
      sum_Z     += ewcharge_.tau_[MAXNF+f] + ewcharge_.tau_[MAXNF-f];
    }
  for (int f = 1; f <= MAXNF; f++)
    {
      NFV_gamma[f] = ewcharge_.Q_[MAXNF+f]/sum_gamma;
      NFV_Z[f]     = ewcharge_.tau_[MAXNF+f]/sum_Z;
    }  
  double S1 = -CF*pi2/12.;
  double S2 = CF*(9.*CF*pi4 + CA*(4856. - 603.*pi2 + 18.*pi4 - 2772.*zeta3) + NF*(-656. + 90.*pi2 + 504.*zeta3))/2592.;
  double S3 = CF*(pow(CA,2)*(689661*pi4 - 55548*pi6 + 105*pi2*(-297481 + 89100*zeta3) + 35*(5211949 - 8161128*zeta3 + 1353024*pow(zeta3,2) + 2630232*zeta5))
		  - 7*(405*pow(CF,2)*pi6 + 4*NF2*(640 + 13770*pi2 + 297*pi4 + 37800*zeta3) - 27*CF*NF*(6*pi4 - 5*pi2*(-3131 + 2664*zeta3) + 5*(-42727 + 20928*zeta3 + 12096*zeta5)))
		  - 7*CA*(135*CF*pi2*(4856 - 603*pi2 + 18*pi4 - 2772*zeta3) + 2*NF*(5616*pi4 - 75*pi2*(7453 + 324*zeta3) + 5*(412765 - 660312*zeta3 + 402408*zeta5))))/29393280.;

  double K1 = CF*(-4 + 7*pi2/12.);
  double K2 = CF*(27.*CF*(7665. - 1660.*pi2 + 134.*pi4 - 3600.*zeta3) + 10.*NF*(4085. - 546.*pi2 + 72.*zeta3) + CA*(-255785. + 31830.*pi2 - 288.*pi4 + 112680.*zeta3))/25920.;
  double K3 = CF*(NF2*(-190931./419904. + (403*pi2)/3888. + (43*pi4)/38880. - (13*zeta3)/243.)
		  + CA*CF*(824281./20736. - (406507*pi2)/62208. + (92237*pi4)/155520. - (739*pi6)/54432. - (13141*zeta3)/432. + (845*pi2*zeta3)/288. + (37*pow(zeta3,2))/12. - (689*zeta5)/72.)
		  + pow(CA,2)*(-51082685./1679616. + (596513*pi2)/139968. - (4303*pi4)/311040. + (299*pi6)/102060. + (505087*zeta3)/15552. - (73*pi2*zeta3)/36. - (71*pow(zeta3,2))/18. - (217*zeta5)/144.)
		  + CF*NF*(-56963./31104. + (13705*pi2)/15552. - (1463*pi4)/15552. +  (815*zeta3)/162. - (37*pi2*zeta3)/144. - (13*zeta5)/9.)
		  + ((-4 + pow(NC,2))*NFV*(1./8. + (5*pi2)/96. - pi4/2880. + (7*zeta3)/48. - (5*zeta5)/6.))/NC
		  + CA*NF*(1700171./209952. - (201749*pi2)/139968. -  (35*pi4)/15552. - (134*zeta3)/27. + (37*pi2*zeta3)/144. -  zeta5/24.)
		  + pow(CF,2)*(-5599/384. + (4339*pi2)/2304. - (173*pi4)/480. +  (27403*pi6)/1088640. - (115*zeta3)/16. - (35*pi2*zeta3)/48. + pow(zeta3,2)/2. + (83*zeta5)/4.));
  double K3_NFV = ( (-4 + pow(NC,2)) *(1./8. + (5*pi2)/96. - pi4/2880. + (7*zeta3)/48. - (5*zeta5)/6.) ) /NC; 

  C1qq_delta = (K1+S1)/2.;
  C2qq_delta = (K2+S2+K1*S1)/2. - pow(C1qq_delta,2)/2.;
  C3qq_delta = (K3 + S3 + K2*S1 + K1*S2)/2. - C1qq_delta*C2qq_delta;

  //C3qq_delta_NFV = K3_NFV/2.;
  C3qq_delta_NFV = 0.;
    
  icoeff::delta();
}

void ccoeff::calc1d()
{
  //analytical C1 and C2 from h2calc
  if (false)
    for (int m = 0; m < mellinint::mdim; m++)
      {
	fcomplex fxn; //input of ancalc
	fxn.real = real(mellinint::Np[m]);
	fxn.imag = imag(mellinint::Np[m]);
	complex <double> cxn = mellinint::Np[m];
	fcomplex fC2qg,fC2NSqqb,fC2NSqq,fC2Sqqb;
	dyh2calc_(fC2qg,fC2NSqqb,fC2NSqq,fC2Sqqb,fxn);

	int idx = anomalous::index(m,mesq::positive);
	
	C1qg[idx] = 1./((cxn+1.)*(cxn+2.))                             /2.;
	C1qq[idx] = (2.*constants::pi2/3.-16./3.+4./3./(cxn*(cxn+1.))) /2.;

	C2qg[idx]   =  cx(fC2qg)                 /4.;
	C2qq[idx]   = (cx(fC2NSqq)+cx(fC2Sqqb))  /4.;
	C2qqb[idx]  = (cx(fC2NSqqb)+cx(fC2Sqqb)) /4.;
	C2qqp[idx]  =  cx(fC2Sqqb)               /4.;
	C2qqbp[idx] =  cx(fC2Sqqb)               /4.;
      }

  //Evaluate numerically C1, C2 and C3 with I functions
  if (false)
    {
      icoeff::calc1d();
      for (int m = 0; m < mellinint::mdim; m++)
	{
	  int idx = anomalous::index(m,mesq::positive);
	  C1qg[idx]   = icoeff::c1qg[idx];
	  C1qq[idx]   = icoeff::c1qq[idx];
	  
	  C2qg[idx]   = icoeff::c2qg[idx];
	  C2qq[idx]   = icoeff::c2qq[idx];
	  C2qqb[idx]  = icoeff::c2qqb[idx];
	  C2qqp[idx]  = icoeff::c2qqp[idx];
	  C2qqbp[idx] = icoeff::c2qqbp[idx];
	  
	  C3qg[idx]   = icoeff::c3qg[idx];
	  C3qq[idx]   = icoeff::c3qq[idx];
	  C3qqb[idx]  = icoeff::c3qqb[idx];
	  C3qqp[idx]  = icoeff::c3qqp[idx];
	  C3qqbp[idx] = icoeff::c3qqbp[idx];
	}
    }

  //Analytic C1, C2 and C3 from systematic expansion of harmonic sums
  for (int m = 0; m < mellinint::mdim; m++)
    {
      int idx = anomalous::index(m,mesq::positive);
      complex <double> N = mellinint::Np[m]-1.;
      C1qg[idx]   = C1qgN(N);
      C1qq[idx]   = C1qqN(N);

      C2qg[idx]   = C2qgN(N);
      C2qq[idx]   = C2qqN(N);
      C2qqb[idx]  = C2qqbN(N);
      C2qqp[idx]  = C2qqpN(N);
      C2qqbp[idx] = C2qqp[idx];

      C3qg[idx]   = C3qgN(N);
      C3qq[idx]   = C3qqN(N);
      C3qqb[idx]  = C3qqbN(N);
      C3qqp[idx]  = C3qqpN(N);
      C3qqbp[idx] = C3qqbpN(N);
    }
  
  //Compute negative branch
  for (int m = 0; m < mellinint::mdim; m++)
   {
     int idxp = anomalous::index(m,mesq::positive);
     int idxm = anomalous::index(m,mesq::negative);

     C1qg[idxm]   = conj(C1qg[idxp]);
     C1qq[idxm]   = conj(C1qq[idxp]);

     C2qg[idxm]   = conj(C2qg[idxp]);
     C2qq[idxm]   = conj(C2qq[idxp]);
     C2qqb[idxm]  = conj(C2qqb[idxp]);
     C2qqp[idxm]  = conj(C2qqp[idxp]);
     C2qqbp[idxm] = conj(C2qqbp[idxp]);

     C3qg[idxm]   = conj(C3qg[idxp]);
     C3qq[idxm]   = conj(C3qq[idxp]);
     C3qqb[idxm]  = conj(C3qqb[idxp]);
     C3qqp[idxm]  = conj(C3qqp[idxp]);
     C3qqbp[idxm] = conj(C3qqbp[idxp]);
   }

  //Checks
  /*
  for (int m = 0; m < mellinint::mdim; m++)
   {
     complex <double> cxn = mellinint::Np[m];// - 1.;
     fcomplex fxn; //input of ancalc
     fxn.real = real(cxn);
     fxn.imag = imag(cxn);
     fcomplex fC2qg,fC2NSqqb,fC2NSqq,fC2Sqqb;
     h2calc_(fC2qg,fC2NSqqb,fC2NSqq,fC2Sqqb,fxn);
     fcomplex dyfC2qg,dyfC2NSqqb,dyfC2NSqq,dyfC2Sqqb;
     dyh2calc_(dyfC2qg,dyfC2NSqqb,dyfC2NSqq,dyfC2Sqqb,fxn);
     
      int idx = anomalous::index(m,mesq::positive);
      complex <double> N = cxn-1.;
      cout << endl;
      cout << " anomalous    expanded " << endl;
      cout << cxn << " C1qg   " << 1./((cxn+1.)*(cxn+2.))                             /2.  - C1qgN(N)  << endl;
      cout << cxn << " C1qq   " << (2.*constants::pi2/3.-16./3.+4./3./(cxn*(cxn+1.))) /2.  - C1qqN(N)  << endl;


      cout << cxn << " C2qg   " << cx(fC2qg)/4.                 - C2qgN(N)   << "  " << cx(dyfC2qg)/4.                    - C2qgN(N)  << endl;
      cout << cxn << " C2qq   " << (cx(fC2NSqq)+cx(fC2Sqqb))/4. - C2qqN(N)   << "  " << (cx(dyfC2NSqq)+cx(dyfC2Sqqb))/4.  - C2qqN(N)  << endl;
      cout << cxn << " C2qqp  " << cx(fC2Sqqb)/4.               - C2qqpN(N)  << "  " << cx(dyfC2Sqqb)/4.                  - C2qqpN(N) << endl;
      cout << cxn << " C2qqb  " << (cx(fC2NSqqb)+cx(fC2Sqqb))/4.- C2qqbN(N)  << "  " << (cx(dyfC2NSqqb)+cx(dyfC2Sqqb))/4. - C2qqbN(N) << endl;
    }
  */
  
  /*
  for (int m = 0; m < mellinint::mdim; m++)
   {
      int idx = anomalous::index(m,mesq::positive);
      complex <double> N = mellinint::Np[m]-1.;
      cout << endl;
      cout << " anomalous    expanded    numerical" << endl;
      cout << mellinint::Np[m] << " C2qg   " << anomalous::C2qgM[idx]/4.                              << "  " << C2qgN(N)  << "  " << icoeff::c2qg[idx]  << endl;
      cout << mellinint::Np[m] << " C2qq   " << (anomalous::C2NSqqM[idx]+anomalous::C2SqqbM[idx])/4.  << "  " << C2qqN(N)  << "  " << icoeff::c2qq[idx]  << endl;
      cout << mellinint::Np[m] << " C2qqp  " << anomalous::C2SqqbM[idx]/4.                            << "  " << C2qqpN(N) << "  " << icoeff::c2qqp[idx] << endl;
      cout << mellinint::Np[m] << " C2qqb  " << (anomalous::C2NSqqbM[idx]+anomalous::C2SqqbM[idx])/4. << "  " << C2qqbN(N) << "  " << icoeff::c2qqb[idx] << endl;

      cout << endl;      
      cout << mellinint::Np[m] << " C3qg   " << icoeff::c3qg[idx] - C3qgN(N)	   << endl;
      cout << mellinint::Np[m] << " C3qq   " << icoeff::c3qq[idx] - C3qqN(N)	   << endl;
      cout << mellinint::Np[m] << " C3qqp  " << icoeff::c3qqp[idx] - C3qqpN(N)   << endl;
      cout << mellinint::Np[m] << " C3qqb  " << icoeff::c3qqb[idx] - C3qqbN(N)   << endl;
      cout << mellinint::Np[m] << " C3qqbp " << icoeff::c3qqbp[idx] - C3qqbpN(N) << endl;
    }
  */
  
/*
  int prec = 8;
  int dist = 2*prec+20;
  cout << setprecision(8) << endl;
  cout << setw(10) << "C2qqn = " << resconst::C2qqn  << endl;
  cout << setw(10) << "H2qqD0 "    << resconst::H2qqD0         << endl;
  cout << setw(10) << "H2qqdelta " << resconst::H2qqdelta      << endl;
  cout << setw(dist) << " N "            << setw(10) << "coeff "     << setw(dist) << " analytical " << " function " << setw(dist) << " numerical "       << setw(dist/2) <<" Rel Error " << endl;
  for (int m = 0; m < mellinint::mdim; m++)
   {
     fcomplex xn = fcx(mellinint::Np[m]);
     fcomplex ps0n;
     psi0_(xn,ps0n);
     complex <double> s1nm1 = cx(ps0n)+resconst::Euler;

     int idx = anomalous::index(m,mesq::positive);
     cout << setw(dist) << mellinint::Np[m] << setw(10) << "C2QG "      << setw(dist) << anomalous::C2qgM[idx]    << " 4*c2qg          " << setw(dist) << 4.*c2qg[idx]        << setw(dist/2) << fabs(anomalous::C2qgM[idx]/4./c2qg[idx]-1.) << endl;
     cout << setw(dist) << mellinint::Np[m] << setw(10) << "C2NSQQ "    << setw(dist) << anomalous::C2NSqqM[idx]  << " 4*(c2qqreg+C2qqn+H2qqD0/2-c2qqp)       " << setw(dist) << 4.*(c2NSqq[idx]+resconst::C2qqn+(-s1nm1*resconst::H2qqD0/2.))      << setw(dist/2) << fabs(anomalous::C2NSqqM[idx]/(4.*(c2NSqq[idx]+resconst::C2qqn+(-s1nm1*resconst::H2qqD0/2.)))-1.) << endl;
     cout << setw(dist) << mellinint::Np[m] << setw(10) << "C2SQQB "    << setw(dist) << anomalous::C2SqqbM[idx]  << " 4*c2qqp         " << setw(dist) << 4.*c2Sqqb[idx]      << setw(dist/2) << fabs(anomalous::C2SqqbM[idx]/4./c2Sqqb[idx]-1.) << endl;
     cout << setw(dist) << mellinint::Np[m] << setw(10) << "C2NSQQB "   << setw(dist) << anomalous::C2NSqqbM[idx] << " 4*(c2qqb-c2qqp)         " << setw(dist) << 4.*c2NSqqb[idx]     << setw(dist/2) << fabs(anomalous::C2NSqqbM[idx]/4./c2NSqqb[idx]-1.) << endl;
     cout << setw(dist) << mellinint::Np[m] << setw(10) << "ps0n "      << setw(dist) << cx(ps0n) << endl;
     cout << setw(dist) << mellinint::Np[m] << setw(10) << "s1nm1 "     << setw(dist) << s1nm1 << endl;
     cout << endl;
   }
*/

  /*
  int idx;
  for (int m = 0; m < mellinint::mdim; m++)
   {
     fcomplex xn = fcx(mellinint::Np[m]);
     fcomplex ps0n;
     psi0_(xn,ps0n);
     complex <double> s1nm1 = cx(ps0n)+resconst::Euler;

     int idxp = anomalous::index(m,mesq::positive);
     anomalous::C2qgM[idxp]    = 4.*c2qg[idxp];
     anomalous::C2NSqqM[idxp]  = 4.*(c2NSqq[idxp]+resconst::C2qqn+(-s1nm1*resconst::H2qqD0/2.));
     anomalous::C2SqqbM[idxp]  = 4.*c2Sqqb[idxp];
     anomalous::C2NSqqbM[idxp] = 4.*c2NSqqb[idxp];                                                
     int idxm = anomalous::index(m,mesq::negative);
     anomalous::C2qgM[idxm]    = conj(anomalous::C2qgM[idxp]);
     anomalous::C2NSqqM[idxm]  = conj(anomalous::C2NSqqM[idxp]);
     anomalous::C2SqqbM[idxm]  = conj(anomalous::C2SqqbM[idxp]);
     anomalous::C2NSqqbM[idxm] = conj(anomalous::C2NSqqbM[idxp]);

     //cout << conj(anomalous::C2qgM[idxp])    <<	 "  " << 4.*conj(c2qg[idxp]) << endl;                                                         
     //cout << conj(anomalous::C2NSqqM[idxp])  <<	 "  " << 4.*(conj(c2NSqq[idxp])+resconst::C2qqn+(-conj(s1nm1)*resconst::H2qqD0/2.)) << endl;  
     //cout << conj(anomalous::C2SqqbM[idxp])  <<	 "  " << 4.*conj(c2Sqqb[idxp]) << endl;						       
     //cout << conj(anomalous::C2NSqqbM[idxp]) <<	 "  " << 4.*conj(c2NSqqb[idxp]) << endl;						       
     //cout << endl;
   }
  */
}

void ccoeff::calc2d()
{
  //analytical C1 and C2
  if (false)
    {
      for (int m = 0; m < mellinint::mdim; m++)
	{
	  int idx = anomalous::index(m,mesq::positive);
      
	  C1qg_1[idx] = anomalous::C1QG_1[idx]/2.;
	  C1qq_1[idx] = anomalous::C1QQ_1[idx]/2.;
	  C2qg_1[idx]   = anomalous::C2qgM_1[idx]/4.;
	  C2qq_1[idx]   = (anomalous::C2NSqqM_1[idx]+anomalous::C2SqqbM_1[idx])/4.;
	  C2qqb_1[idx]  = (anomalous::C2NSqqbM_1[idx]+anomalous::C2SqqbM_1[idx])/4.;
	  C2qqp_1[idx]  = anomalous::C2SqqbM_1[idx]/4.;
	  C2qqbp_1[idx] = anomalous::C2SqqbM_1[idx]/4.;

	  C2qg_2[idx]   = anomalous::C2qgM_2[idx]/4.;
	  C2qq_2[idx]   = (anomalous::C2NSqqM_2[idx]+anomalous::C2SqqbM_2[idx])/4.;
	  C2qqb_2[idx]  = (anomalous::C2NSqqbM_2[idx]+anomalous::C2SqqbM_2[idx])/4.;
	  C2qqp_2[idx]  = anomalous::C2SqqbM_2[idx]/4.;
	  C2qqbp_2[idx] = anomalous::C2SqqbM_2[idx]/4.;
	  C1qg_2[idx] = anomalous::C1QG_2[idx]/2.;
	  C1qq_2[idx] = anomalous::C1QQ_2[idx]/2.;
	}

      //Evaluate C3 coefficients with I functions
      icoeff::calc2d();

      for (int m = 0; m < mellinint::mdim; m++)
	{
	  int idx = anomalous::index(m,mesq::positive);

	  C3qg_1[idx]   = icoeff::c3qg_1[idx];
	  C3qq_1[idx]   = icoeff::c3qq_1[idx];
	  C3qqb_1[idx]  = icoeff::c3qqb_1[idx];
	  C3qqp_1[idx]  = icoeff::c3qqp_1[idx];
	  C3qqbp_1[idx] = icoeff::c3qqbp_1[idx];

	  C3qg_2[idx]   = icoeff::c3qg_2[idx];
	  C3qq_2[idx]   = icoeff::c3qq_2[idx];
	  C3qqb_2[idx]  = icoeff::c3qqb_2[idx];
	  C3qqp_2[idx]  = icoeff::c3qqp_2[idx];
	  C3qqbp_2[idx] = icoeff::c3qqbp_2[idx];
	}
    }

  //Analytic C1, C2 and C3 from systematic expansion of harmonic sums
  //if (false)
    for (int m = 0; m < mellinint::mdim; m++)
      {
	int idx = anomalous::index(m,mesq::positive);

	complex <double> N_1 = mellinint::Np_1[m]-1.;
	C1qg_1[idx]   = C1qgN(N_1);
	C1qq_1[idx]   = C1qqN(N_1);
	C2qg_1[idx]   = C2qgN(N_1);
	C2qq_1[idx]   = C2qqN(N_1);
	C2qqb_1[idx]  = C2qqbN(N_1);
	C2qqp_1[idx]  = C2qqpN(N_1);
	C2qqbp_1[idx] = C2qqp_1[idx];
	C3qg_1[idx]   = C3qgN(N_1);
	C3qq_1[idx]   = C3qqN(N_1);
	C3qqb_1[idx]  = C3qqbN(N_1);
	C3qqp_1[idx]  = C3qqpN(N_1);
	C3qqbp_1[idx] = C3qqbpN(N_1);

	complex <double> N_2 = mellinint::Np_2[m]-1.;
	C1qg_2[idx]   = C1qgN(N_2);
	C1qq_2[idx]   = C1qqN(N_2);
	C2qg_2[idx]   = C2qgN(N_2);
	C2qq_2[idx]   = C2qqN(N_2);
	C2qqb_2[idx]  = C2qqbN(N_2);
	C2qqp_2[idx]  = C2qqpN(N_2);
	C2qqbp_2[idx] = C2qqp_2[idx];
	C3qg_2[idx]   = C3qgN(N_2);
	C3qq_2[idx]   = C3qqN(N_2);
	C3qqb_2[idx]  = C3qqbN(N_2);
	C3qqp_2[idx]  = C3qqpN(N_2);
	C3qqbp_2[idx] = C3qqbpN(N_2);
      }
  
  //Compute negative branch
  for (int m = 0; m < mellinint::mdim; m++)
   {
     int idxp = anomalous::index(m,mesq::positive);
     int idxm = anomalous::index(m,mesq::negative);

     C1qg_1[idxm]   = conj(C1qg_1[idxp]);
     C1qq_1[idxm]   = conj(C1qq_1[idxp]);
     C2qg_1[idxm]   = conj(C2qg_1[idxp]);
     C2qq_1[idxm]   = conj(C2qq_1[idxp]);
     C2qqb_1[idxm]  = conj(C2qqb_1[idxp]);
     C2qqp_1[idxm]  = conj(C2qqp_1[idxp]);
     C2qqbp_1[idxm] = conj(C2qqbp_1[idxp]);
     C3qg_1[idxm]   = conj(C3qg_1[idxp]);
     C3qq_1[idxm]   = conj(C3qq_1[idxp]);
     C3qqb_1[idxm]  = conj(C3qqb_1[idxp]);
     C3qqp_1[idxm]  = conj(C3qqp_1[idxp]);
     C3qqbp_1[idxm] = conj(C3qqbp_1[idxp]);

     C1qg_2[idxm]   = conj(C1qg_2[idxp]);
     C1qq_2[idxm]   = conj(C1qq_2[idxp]);
     C2qg_2[idxm]   = conj(C2qg_2[idxp]);
     C2qq_2[idxm]   = conj(C2qq_2[idxp]);
     C2qqb_2[idxm]  = conj(C2qqb_2[idxp]);
     C2qqp_2[idxm]  = conj(C2qqp_2[idxp]);
     C2qqbp_2[idxm] = conj(C2qqbp_2[idxp]);
     C3qg_2[idxm]   = conj(C3qg_2[idxp]);
     C3qq_2[idxm]   = conj(C3qq_2[idxp]);
     C3qqb_2[idxm]  = conj(C3qqb_2[idxp]);
     C3qqp_2[idxm]  = conj(C3qqp_2[idxp]);
     C3qqbp_2[idxm] = conj(C3qqbp_2[idxp]);
   }
}

void ccoeff::num_calc()
{
  //Evaluate C1 and C2 coefficients numerically
  int rule = 500;
  double *tlog,*tlin;              //gauss nodes
  double *faclog,*faclin;            //overall factor
  complex <double> *kernlog,*kernlin; //kernel of the Mellin transform
  //memory leak here!!!
  tlog = new double [rule];
  faclog = new double [rule];
  kernlog = new complex <double> [rule*mellinint::mdim*2];
  tlin = new double [rule];
  faclin = new double [rule];
  kernlin = new complex <double> [rule*mellinint::mdim*2];

  // boundaries of integration
  double xmax = 1;
  //double xminlog = pow(bins.mbins.front()/opts.sroot,2);
  //double xminlin = pow(bins.mbins.front()/opts.sroot,2);
  double xminlog = 1e-14;
  double xminlin = 0;
  double ll = log(xmax/xminlog);
  double cc = 0.5;
  double mm = 0.5;
  for (int i = 0; i < rule; i++)
    {
      double x = cc+mm*gr::xxx[rule-1][i];
      tlog[i] = xminlog*exp(ll*x);
      tlin[i] = xminlin+(xmax-xminlin)*x;
      double jaclog = mm * tlog[i] * ll;
      double jaclin = mm*(xmax-xminlin);
      faclog[i] = jaclog * gr::www[rule-1][i];
      faclin[i] = jaclin * gr::www[rule-1][i];
      for (int m = 0; m < mellinint::mdim; m++)
	{
	  kernlog[i*mellinint::mdim+m] = pow(tlog[i], mellinint::Np[m]-1.);
	  kernlin[i*mellinint::mdim+m] = pow(tlin[i], mellinint::Np[m]-1.);
	}
    }

  fill(C1qg,   C1qg+2*mellinint::mdim, 0.);
  fill(C1qq,   C1qq+2*mellinint::mdim, 0.);
  fill(C1qqb,  C1qqb+2*mellinint::mdim, 0.);
  fill(C1qqp,  C1qqp+2*mellinint::mdim, 0.);
  fill(C1qqbp, C1qqbp+2*mellinint::mdim, 0.);

  fill(C2qg,   C2qg+2*mellinint::mdim, 0.);
  fill(C2qq,   C2qq+2*mellinint::mdim, 0.);
  fill(C2qqb,  C2qqb+2*mellinint::mdim, 0.);
  fill(C2qqp,  C2qqp+2*mellinint::mdim, 0.);
  fill(C2qqbp, C2qqbp+2*mellinint::mdim, 0.);

  //Calculate Mellin moments as:
  //integral_0^1{ x^(N-1) fx dx}
  for (int i = 0; i < rule; i++)
    {
      double zlog = tlog[i];
      double zlin = tlin[i];

      double c1qgzlin = cqg_(zlin);
      double c1qqzlin = cqq_(zlin);
      
      double c2qgzlin = c2qg_(zlin);
      double c2qqzlin = c2qqreg_(zlin);
      double c2qqbzlog = c2qqb_(zlog);
      double c2qqpzlog = c2qqp_(zlog);
      double c2qqbpzlog = c2qqpzlog; //=c2qqp_(zlog);
      for (int m = 0; m < mellinint::mdim; m++)
	{
	  int idx = anomalous::index(m,mesq::positive);
	  
	  C1qg[idx] += faclin[i]*c1qgzlin * kernlin[i*mellinint::mdim+m];
	  C1qq[idx] += faclin[i]*c1qqzlin * kernlin[i*mellinint::mdim+m];

	  C2qg[idx]   += faclin[i]*kernlin[i*mellinint::mdim+m] * c2qgzlin;
	  C2qq[idx]   += faclin[i]*kernlin[i*mellinint::mdim+m] * c2qqzlin;
	  C2qqb[idx]  += faclog[i]*kernlog[i*mellinint::mdim+m] * c2qqbzlog;
	  C2qqp[idx]  += faclog[i]*kernlog[i*mellinint::mdim+m] * c2qqpzlog;
	  C2qqbp[idx] += faclog[i]*kernlog[i*mellinint::mdim+m] * c2qqbpzlog;

	  //cout << i << "  " << kern[i*mellinint::mdim+m] << endl;
	}
    }

  //Add delta pieces
  for (int m = 0; m < mellinint::mdim; m++)
    {
      int idx = anomalous::index(m,mesq::positive);
      C1qq[idx] += resconst::C1qqn;
      C2qq[idx] += resconst::C2qqn;
    }

  //Add plus-distribution pieces: the Mellin transform of 1/(1-z)+ is -S1(N-1)
  for (int m = 0; m < mellinint::mdim; m++)
    {
      complex <double> s1nm1 = cpsi0(mellinint::Np[m])+constants::euler;
      int idx = anomalous::index(m,mesq::positive);
      C2qq[idx] += (-s1nm1)*resconst::H2qqD0/2.;
    }
}

void ccoeff::truncate()
{
  //start from original untruncated moments
  for (int m = 0; m < mellinint::mdim; m++)
    {
      int idxp = anomalous::index(m,mesq::positive);
	  
      C1qg[idxp] = anomalous::C1QG[idxp]/2.;
      C1qq[idxp] = anomalous::C1QQ[idxp]/2.;

      int idxm = anomalous::index(m,mesq::negative);
      C1qg[idxm]   = conj(C1qg[idxp]);
      C1qq[idxm]   = conj(C1qq[idxp]);
    }
  
  //Calculate truncated moments
  complex <double> C1qg_tr[mellinint::mdim];
  complex <double> C1qq_tr[mellinint::mdim];

  double x1 = phasespace::m/opts.sroot*exp(phasespace::ymin);
  double x2 = phasespace::m/opts.sroot*exp(-phasespace::ymax);

  //double x1 = 1e-8;//pow(phasespace::m/opts.sroot,2);
  //double x2 = 1e-8;//pow(phasespace::m/opts.sroot,2);

  double lx1 = log(x1);
  double lx2 = log(x2);
  
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
  complex <double> c1qq_p[mellinint::mdim];
  complex <double> c1qg_p[mellinint::mdim];
  complex <double> c1qq_m[mellinint::mdim];
  complex <double> c1qg_m[mellinint::mdim];
  for (int m = 0; m < mellinint::mdim; m++)
    {
      c1qq_p[m] = facp * C1qq[m] * mellinint::wn[m];
      c1qg_p[m] = facp * C1qg[m] * mellinint::wn[m];
      c1qq_m[m] = facm * conj(C1qq[m]) * mellinint::wn[m];
      c1qg_m[m] = facm * conj(C1qg[m]) * mellinint::wn[m];
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
	{
	    /*
	    if (m == n)
	      {
		ffx1[n][f] += fm1p[m][f] * (-lx1);
		ffx2[n][f] += fm2p[m][f] * (-lx2);
	      }
	    else
	      {
		ffx1[n][f] += fm1p[m][f] * (1.-x1n[n]/x1n[m])/(mellinint::Np[n]-mellinint::Np[m]);
		ffx2[n][f] += fm2p[m][f] * (1.-x2n[n]/x2n[m])/(mellinint::Np[n]-mellinint::Np[m]);
	      }
	    */
	    //	    ffx1[n][f] += fm1p[m][f]*llx1p[n][m]; 
	    //	    ffx2[n][f] += fm2p[m][f]*llx2p[n][m];
	  C1qq_tr[n] += c1qq_p[m]*llx1p[n][m] - c1qq_m[m]*llx1m[n][m]; 
	  C1qg_tr[n] += c1qg_p[m]*llx2p[n][m] - c1qg_m[m]*llx2m[n][m];
	}

      /*
      //negative branch
      for (int m = 0; m < mellinint::mdim; m++)
	for (int f = 0; f < 2*MAXNF+1; f++)
	  {
	    //	    ffx1[n][f] -= fm1m[m][f]*(1.-x1n[n]/conj(x1n[m]))/(mellinint::Np[n]-mellinint::Nm[m]);
	    //	    ffx2[n][f] -= fm2m[m][f]*(1.-x2n[n]/conj(x2n[m]))/(mellinint::Np[n]-mellinint::Nm[m]);
	    ffx1[n][f] -= fm1m[m][f]*llx1m[n][m];
	    ffx2[n][f] -= fm2m[m][f]*llx2m[n][m];
	  }
      */

  
      //cout << "truncated " << n << ffx1[5] << endl;
    }

  //replace moments
  for (int m = 0; m < mellinint::mdim; m++)
    {
      int idxp = anomalous::index(m,mesq::positive);
	  
      C1qg[idxp] = C1qg_tr[m];
      C1qq[idxp] = C1qq_tr[m];

      int idxm = anomalous::index(m,mesq::negative);
      C1qg[idxm]   = conj(C1qg[idxp]);
      C1qq[idxm]   = conj(C1qq[idxp]);
    }

}

void ccoeff::free()
{
  delete[] C1qg;
  delete[] C1qq;
  delete[] C1qqb;
  delete[] C1qqp;
  delete[] C1qqbp;
  delete[] C2qg;
  delete[] C2qq;
  delete[] C2qqb;
  delete[] C2qqp;
  delete[] C2qqbp;
  delete[] C3qg;
  delete[] C3qq;
  delete[] C3qqb;
  delete[] C3qqp;
  delete[] C3qqbp;

}
void ccoeff::release()
{
  if (opts.melup <= 1)
    free();
  
  icoeff::release();
}
