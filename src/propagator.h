#ifndef propagator_h
#define propagator_h

#include "fcomplex.h"

namespace prop
{
  extern double m;
  extern double w;
  extern double m2;
  extern double w2;

  extern double W;
  extern double Z;
  extern double G;
  extern double ZG;
#pragma omp threadprivate(W,Z,G,ZG)
  
  extern void init();

  extern complex <double> bos(double q2);
  extern double gamma(double q2);
  extern void set(double q2);
  extern void set_real(double q2);
  extern void set_offset(double q2);
}

//fortran interface
extern "C"
{
  fcomplex bosprop_(double &q2);
}

#endif
