#ifndef itilde_h
#define itilde_h

#include "resint.h"

//workspace for intdeo DEQUAD integration
namespace itilde
{
  extern double aw[lenaw];

  extern double iqt,iq2,ib02;
  extern int in;
#pragma omp threadprivate(iqt,iq2,ib02,in)
  
  void init();
  double besselk(double x, int n);
  double calc(double qt, double q, int n);
  double integrand(double b);
  double logn(double b);
  double integral(double qt, double q, int n);
  double integrand_int(double b);
  double lognob(double b);
}

#endif
