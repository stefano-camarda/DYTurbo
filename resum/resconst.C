#include "resconst.h"
#include "constants.h"
#include <math.h>
#include <iostream>

using namespace std;
using namespace constants;

// constants.f
const int resconst::NF = 5;

//DYRes precision
//const double resconst::Euler = 0.57721566;
//const double resconst::Z2 = 1.644934;
//const double resconst::Z3 = 1.202057;

//higher precision
const double resconst::Z2 = zeta2; //M_PI*M_PI/6.;
const double resconst::Z3 = zeta3; //1.20205690315959428540;
const double resconst::Z4 = zeta4; //1.08232323371113819152;
const double resconst::Z5 = zeta5; //1.03692775514336992633;
const double resconst::b0 = 2.*exp(-euler);

double resconst::beta0, resconst::beta1, resconst::beta2, resconst::beta3, resconst::beta4, resconst::Kappa;
double resconst::A1g, resconst::A2g, resconst::A3g, resconst::B1g, resconst::B2g, resconst::C1ggn;
double resconst::A1q, resconst::A2q, resconst::A3q, resconst:: A4q, resconst:: A5q;
double resconst::B1q, resconst::B2q, resconst::B3q, resconst::B4q;
double resconst::C1qqn, resconst::C2qqn;
double resconst::C1qqdelta, resconst::Delta2qq;
double resconst::D0qqqq, resconst::D1qqqq;
double resconst::Deltaqqqq;
double resconst::H2qqdelta;
double resconst::H2qqD0;

void resconst::init()
{
  double nf2 = NF*NF;
  double nf3 = nf2*NF;

  //Resummation coefficients
  beta0=(33.-2.*NF)/12.;
  beta1=(153.-19.*NF)/24.;
  beta2=2857./128.-5033.*NF/1152.+325.*nf2/3456.;
  beta3=(149753./6. + 3564.*Z3 + NF*(-1078361./162.-6508./27.*Z3)+nf2*(50065./162.+6472./81.*Z3)+1093./729.*nf3)/256.; //from https://arxiv.org/pdf/1701.01404.pdf Eq.(3.6), with a 4^4=256 normalisation factor
  beta4=(8157455./16. + 621885./2.*Z3 - 88209./2.*Z4 - 288090.*Z5
	 +NF* (-336460813./1944. - 4811164./81.*Z3 + 33935./6.*Z4 + 1358995./27.*Z5)
	 +NF*NF*(25960913./1944. + 698531./81.*Z3 - 10526./9.*Z4 - 381760./81.*Z5)
	 +NF*NF*NF*(-630559./5832. - 48722./243.*Z3 + 1618./27.*Z4 + 460./9.*Z5)
	 +NF*NF*NF*NF*(1205./2916. - 152./81.*Z3))/1024.;  //from https://arxiv.org/pdf/1701.01404.pdf Eq.(3.7), with a 4^5=1024 normalisation factor

  //cout << beta0/M_PI << "  " << beta1/pow(M_PI,2) << "  " << beta2/pow(M_PI,3) << "  " << beta3/pow(M_PI,4) << "  " << beta4/pow(M_PI,5) << endl;

  Kappa=67./6.-pi2/2.-5./9.*NF;

  //gluon coefficients
  A1g=CA;
  A2g=CA/2.*(67./6.-pi2/2.-5./9.*NF); //CA/2.*Kappa;
  A3g=CA*(13.81-2.15*NF-nf2/108.);
  B1g=-(11.*CA-2.*NF)/6.;
  B2g=CA*CA*(23./24.+(11.*pi2)/18.-3.*Z3/2.)+CF*NF/2.-CA*NF*(1./12.+pi2/9.)-11./8.*CF*CA;
  C1ggn=(pi2/2.+11./2.+pi2)/2.;

  //quark coefficients
  A1q=CF; //4./3;
  A2q=CF/2.*Kappa; //(67./6.-pi2/2.-5./9.*NF);
  //A3q=CF*(13.81-2.15*NF-nf2/108.)+CF*(CA*(29.9259-28.*Z3)-8.2963*NF/2.)*2.*(beta0*4.)/64.; //A3 from Becher & Neubert
  //cout << "A3q orig " << A3q << endl;
  //Better approximation for A3q:
  //A3q=CF*(13.768339123-2.15*NF-nf2/108.)+CF*(CA*(29.9259-28.*Z3)-8.2963*NF/2.)*2.*(beta0*4.)/64.;
  A3q=CF*( (2*CA*CA*(245./24.-67./9.*Z2+11./6.*Z3+11./5.*Z2*Z2)/8.)
	  +(2*(CF*(-55./24.+2*Z3)+CA*(-209./108.+10./9.*Z2-7./3.*Z3))/8.)*NF-nf2/108.)
    +CF*(CA*((808./27.)-28.*Z3)-(224./54.*2.)*NF/2.)*2.*(beta0*4.)/64.; //A3 from Becher & Neubert (https://arxiv.org/pdf/1007.4005.pdf Eq.(74))
  //cout << "A3q fix " << A3q << endl;

  //A4q from https://arxiv.org/pdf/1705.09127.pdf
  //A4q=(CA*CA*CA*CF*(121./3.*Z3*Z2 - 8789./162.*Z2 - 19093./54.*Z3 - 847./24.*Z4 + 132.*Z5 + 3761815./11664.)
  //     + CA*CA*CF*NF*(-22./3.*Z3*Z2 + 2731./162.*Z2 + 4955./54.*Z3 + 11./6.*Z4 -24*Z5 - 31186./243.)
  //     + CA*CF*CF*NF*(272./9.*Z3 + 11.*Z4 - 7351./144.)
  //     + CA*CF*nf2*(-103./81.*Z2 - 47./27.*Z3 +5./6.*Z4 + 13819./972.)
  //     + CF*CF*nf2*(-38./9.*Z3 - 2*Z4 + 215./24.)
  //     + CF*NF*NF*NF*(-4./9.*Z3-232./729.))/16.;
  //A4q=(-1189.19)/16.; //simplied formula for NF=5
    
  
  //from https://arxiv.org/pdf/1604.01404.pdf
  double gamma1r = CF*CA*(28*Z3-808./27.) + 112.*CF*NF/27.;
  double gamma2r = CF*CA*CA*(-176./3*Z3*Z2 + 6392./81.*Z2 + 12328./27.*Z3 +154./3.*Z4 - 192.*Z5 - 297029./729.)
    + CF*CA*NF* ( -824./81.*Z2 - 904./27.*Z3 +20./3.*Z4 + 62626./729.)
    + CF*nf2* ( -32./9.*Z3 - 1856./729.)
    + CF*CF*NF*(-304./9.*Z3-16.*Z4 + 1711./27.);

  double D2q = -gamma1r;
  double D3q = -gamma2r;

  double A4q_th;
  //A4q_th = 1553./256. //Pade approximant from https://arxiv.org/pdf/hep-ph/0506288.pdf
  //A4q_th = 0.42441*0.0133*pow(M_PI,4); //from Eq. (4.3) of https://arxiv.org/pdf/1805.09638.pdf

  //Exact result for A4q_th from Eq.(3.6) of https://arxiv.org/pdf/1912.12920.pdf
  double DfRAnc = 5./2.;
  double DfRRnc = 5./36.;
  A4q_th =  (CF*CA*CA*CA*(84278./81. - 88400./81.*zeta2 + 20944./27.*zeta3 + 1804.*zeta4 - 352./3.*zeta2*zeta3 - 3608./9.*zeta5 - 16.*pow(zeta3,2) - 2504./3.*zeta6)
	     + DfRAnc*(- 128.*zeta2 + 128./3.*zeta3 + 3520./3.*zeta5 - 384.*pow(zeta3,2) - 992.*zeta6)
	     + CF*CF*CF*NF*(572./9. + 592./3.*zeta3 - 320.*zeta5)
	     + CF*CF*CA*NF*(- 34066./81. + 440./3.*zeta2 + 3712./9.*zeta3 - 176.*zeta4 - 128.*zeta2*zeta3 + 160.*zeta5)
	     + CF*CA*CA*NF*(- 24137./81. + 20320./81.*zeta2 - 23104./27.*zeta3 - 176./3.*zeta4 + 448./3.*zeta2*zeta3 + 2096./9.*zeta5)
	     + NF*DfRRnc*(256.*zeta2 - 256./3.*zeta3 - 1280./3.*zeta5)
	     + CF*CF*NF*NF*(2392./81. - 640./9.*zeta3 + 32.*zeta4)
	     + CF*CA*NF*NF*(923./81. - 608./81.*zeta2 + 2240./27.*zeta3 - 112./3.*zeta4)
	     - CF*NF*NF*NF*(32./81. - 64./27.*zeta3))/pow(4.,4);
  
  A4q = A4q_th
    + (3*beta0*4*D3q + 2*beta1*16*D2q)/256.; //formula from Eq. (71) of https://arxiv.org/pdf/1007.4005.pdf with 4^4 normalisation.

  double A5q_th = 0.42441*0.5*pow(M_PI,5); //from Eq. (3.8) of https://arxiv.org/pdf/1912.12920.pdf
  A5q = 0.;//A5q_th;
  
  B1q=-(3.*CF)/2.;                                                                                //-2;
  B2q=CF*CF*(pi2/4.-3./16.-3.*Z3)+CA*CF*(11.*pi2/36.-193./48.+3.*Z3/2.)+CF*NF*(17./24.-pi2/18.);  //B2q=4./9.*(pi2-3./4.-12.*Z3)+(11./9.*pi2-193./12.+6.*Z3)+NF/6.*(17./3.-4./9.*pi2);
  //B2q=CF*CF*(3./2.*Z2-3./16.-3.*Z3)+CA*CF*(11.*Z2/6.-193./48.+3.*Z3/2.)+CF*NF*(17./24.-Z2/3.);
  //B2q = (13.3447 + 3.4138*NF)/16.; from https://arxiv.org/pdf/1604.01404.pdf, Eq.(13)
  //B3q=-116.685/8. // from https://arxiv.org/pdf/1705.09127.pdf, but does not look OK
  B3q=(7358.86 -721.516*NF + 20.5951* nf2)/64.; // from https://arxiv.org/pdf/1604.01404.pdf, Eq.(13)

  //In Eqs.(10-12) of https://arxiv.org/pdf/1912.12920.pdf there is a B4q coefficient, but it does not seem to be in the same scheme
  B4q = 0.;
  
  C1qqn=CF/2.*(pi2/2.-4.);          // Only delta(1-z) part, i.e. N independent part

  // Only delta(1-z) part, i.e. N independent part
//  C2qqn = 1./2.* (CA*CF*(59.*Z3/18. - 1535./192. + 215.*pi2/216. - pi4/240.)
//		  + CF*NF*(192.*Z3 + 1143. - 152.*pi2)/864.
//		  + CF*CF*(-15.*Z3 + 511./16. - 67.*pi2/12. + 17.*pi4/45.)/4. - CF*CF* pow(pi2-8.,2)/16.);

  C2qqn = 1./2.* (CA*CF*(59.*Z3/18. - 1535./192. + 215*pi2/216. - pi4/240.)
		  + CF*NF*(192.*Z3 + 1143. - 152.*pi2)/864.
		  + CF*CF*(-15.*Z3 + 511./16. - 67.*pi2/12. + 17.*pi4/45.)/4. - CF*CF*pow(pi2-8.,2)/16.);
  
  // Delta term in c1qq coefficient
  C1qqdelta=(pi2-8.)/3.;

  // Delta term in P2qq splitting function (as/M_PI normalization)
  Delta2qq=16./9.*(3./8.-pi2/2.+6.*Z3)+4.*(17./24.+11.*pi2/18.-3.*Z3)-2./3.*NF*(1./6.+2.*pi2/9.);
  Delta2qq=Delta2qq/4.;

  // Coefficients of D0 and D1 in P*P (as/M_PI normalization)
  D0qqqq=8./3.;
  D1qqqq=32./9.;

  // Coefficients of delta(1-z) in P*P
  Deltaqqqq=4./9.*(9./4.-2.*pi2/3.);

  // H2qq contribution: coefficient of delta(1-z)
  H2qqdelta=-2561./144.+127.*NF/72.+3.*pi2/2.-19.*NF*pi2/81.+49.*pi4/324.+58.*Z3/9.+8.*NF*Z3/27.;

  // H2qq contribution: coefficient of D0(z)
  H2qqD0=-404./27.+(56.*NF)/81.+14.*Z3;
  //(23) of https://arxiv.org/pdf/1209.0158.pdf (1/1-z)+)

  //A1q=0; //<--
  //A2q=0; //<--
  //A3q=0;
  //A4q=0;
  //B1q=0; //<--
  //B2q=0;
  //B3q=0;
  //C1qqn=0;
}
