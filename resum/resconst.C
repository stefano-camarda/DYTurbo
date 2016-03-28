#include "resconst.h"
#include "interface.h"
#include <math.h>

// constants.f
const int resconst::NF = 5;

const double resconst::CA = 3.;
const double resconst::Cf = 4./3.;

const double resconst::Euler = 0.57721566;
const double resconst::Z2 = 1.64493406685; //1.644934;
const double resconst::Z3 = 1.20205690316;// 1.202057;
//higher precision
/*
const double resconst::Euler = 0.57721566490153286;
const double resconst::Z2 = M_PI*M_PI/6.;
const double resconst::Z3 = 1.2020569031595942853;
*/

double resconst::b0;
double resconst::beta0, resconst::beta1, resconst::beta2, resconst::Kappa;
double resconst::A1g, resconst::A2g, resconst::A3g, resconst::B1g, resconst::B2g, resconst::C1ggn;
double resconst::A1q, resconst::A2q, resconst::A3q, resconst::B1q, resconst::B2q, resconst::C1qqn;
double resconst::C1qqdelta, resconst::Delta2qq;
double resconst::D0qqqq, resconst::D1qqqq;
double resconst::Deltaqqqq;
double resconst::H2qqdelta;
double resconst::H2qqD0;

void resconst::init()
{
  double pi2 = M_PI*M_PI;
  double pi4 = pi2*pi2;
  double nf2 = NF*NF;

  b0=2*exp(-Euler);

  //Resummation coefficients
  beta0=(33.-2.*NF)/12.;
  beta1=(153.-19.*NF)/24.;
  beta2=2857./128.-5033.*NF/1152.+325.*nf2/3456.;
  Kappa=67./6.-pi2/2.-5./9.*NF;

  //gluon coefficients
  A1g=CA;
  A2g=CA/2.*(67./6.-pi2/2.-5./9.*NF); //CA/2.*Kappa;
  A3g=CA*(13.81-2.15*NF-nf2/108.);
  B1g=-(11.*CA-2.*NF)/6.;
  B2g=CA*CA*(23./24.+(11.*pi2)/18.-3.*Z3/2.)+Cf*NF/2.-CA*NF*(1./12.+pi2/9.)-11./8.*Cf*CA;
  C1ggn=(pi2/2.+11./2.+pi2)/2.;

  //quark coefficients
  A1q=Cf; //4./3;
  A2q=Cf/2.*Kappa; //(67./6.-pi2/2.-5./9.*NF);
  A3q=Cf*(13.81-2.15*NF-nf2/108.)+Cf*(CA*(29.9259-28.*Z3)-8.2963*NF/2.)*2.*(beta0*4.)/64.; //A3 from Becher & Neubert
  B1q=-(3.*Cf)/2.; //-2;
  B2q=Cf*Cf*(pi2/4.-3./16.-3.*Z3)+CA*Cf*(11*pi2/36.-193./48.+3.*Z3/2.)+Cf*NF*(17./24.-pi2/18.);
  //B2q=4./9.*(pi2-3./4.-12.*Z3)+(11./9.*pi2-193./12.+6.*Z3)+NF/6.*(17./3.-4./9.*pi2);
  C1qqn=Cf/2.*(pi2/2.-4.);          // Only delta(1-z) part, i.e. N independent part

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
}
