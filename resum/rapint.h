#ifndef rapint_h
#define rapint_h
#include <complex>

#include "mellinint.h"
#include "settings.h"

using namespace std;

namespace rapint
{
  extern int ydim;
  void init();
  void release();

  extern complex <double> * cfpy;
  extern complex <double> * cfmy;
  void cache(double ymin, double ymax);

  extern complex <double> *Ith0p;
  extern complex <double> *Ith1p;
  extern complex <double> *Ith2p;
  extern complex <double> *Ith0m;
  extern complex <double> *Ith1m;
  extern complex <double> *Ith2m;
#pragma omp threadprivate(Ith0p,Ith1p,Ith2p,Ith0m,Ith1m,Ith2m)

  void allocate();
  void integrate(double m, double ymin, double ymax, int helicity = -1);
  void numint(double m, double ymin, double ymax, int helicity = -1);
  void free();

  inline int index(int i, int j, int i1, int i2)
  {return i2 + mellinint::mdim*(i1 + mellinint::mdim*(j + opts.yrule*i));}
}

extern "C"
{
  inline void rapint_cache_(double& ymin, double& ymax) {rapint::cache(ymin, ymax);};
  inline void rapint_integrate_(double& ymin, double& ymax, double& m) {rapint::integrate(ymin, ymax, m);};
}


#endif
