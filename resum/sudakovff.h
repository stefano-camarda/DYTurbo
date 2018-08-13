#ifndef sudakovff_h
#define sudakovff_h

#include "fcomplex.h"
#include "dyres_interface.h"
#include "settings.h"

extern "C"
{
  fcomplex f0_(fcomplex y);
  fcomplex f1_(fcomplex y);
  fcomplex f2_(fcomplex y);

  extern struct {
    fcomplex log1y_;
  } csud_;

  
}
#pragma omp threadprivate(csud_)


namespace sudakov
{
  extern int nf;
  extern double beta0,beta1,beta2;
#pragma omp threadprivate(nf,beta0,beta1,beta2)

  void setnf(int nff);

  complex <double> f0(complex <double> y);
  complex <double> f1(complex <double> y);
  complex <double> sff(complex <double> b);
};

#endif
