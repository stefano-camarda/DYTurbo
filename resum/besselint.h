#ifndef besselint_h
#define besselint_h

#include "interface.h"
#include "fcomplex.h"

namespace besselint
{
  extern complex <double> bint(complex <double> b);
}

//fortran interface
extern "C"
{
  void besselint_(double &b, double &qt, double &q2);
  double fort_besj0_(double &x);
  double fort_besy0_(double &x);
  double fort_besj1_(double &x);
  double fort_besy1_(double &x);
  double fort_besjn_(int &n, double &x);
  double fort_besyn_(int &n, double &x);
  fcomplex s_(fcomplex &b);
  fcomplex alphasl_(fcomplex &q2);

  //exponents computed in alphasl
  extern struct {
    fcomplex aexp_;
    fcomplex aexpb_;
  } exponent_;
}
#pragma omp threadprivate(exponent_)

#endif
