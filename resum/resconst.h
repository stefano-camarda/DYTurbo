#ifndef resconst_h
#define resconst_h

extern const double CA;
extern const double Cf;
extern const double Euler;
extern const double Z2;
extern const double Z3;

namespace resconst
{
  extern void init();

  //Resummation coefficients
  extern double b0;
  extern double beta0, beta1, beta2, Kappa;
  extern double A1g, A2g, A3g, B1g, B2g, C1ggn;
  extern double A1q, A2q, A3q, B1q, B2q, C1qqn;

  // Delta terms
  extern double C1qqdelta, Delta2qq;

  // Coefficients of D0 and D1 in P*P (as/M_PI normalization)
  extern double D0qqqq, D1qqqq;

  // Coefficients of delta(1-z) in P*P
  extern double Deltaqqqq;

  // H2qq contribution: coefficient of delta(1-z)
  extern double H2qqdelta;

  // H2qq contribution: coefficient of D0(z)
  extern double H2qqD0;
}

#endif
