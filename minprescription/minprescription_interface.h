#ifndef minprescription_h
#define minprescription_h

#include "fcomplex.h"

extern "C"
{

  double adzint_ (double (*F)(double &x), double &A, double &B, double &AERR, double &RERR, double &ERREST, double &IER, double &IACTA, double &IACTB);

  fcomplex h1_(fcomplex &z);
  fcomplex h2_(fcomplex &z);

  //  void ch12n_ (int &n, fcomplex &z, int &nm, fcomplex *chf1, fcomplex *chd1, fcomplex *chf2, fcomplex *chd2 );

  struct {
    double v_;
  } v_;
  //#pragma omp threadprivate(v_)  
}

#endif
