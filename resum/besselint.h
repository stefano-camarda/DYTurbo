#ifndef besselint_h
#define besselint_h

#include "interface.h"
#include "fcomplex.h"

namespace besselint
{
  extern double bint(double b);
}

//fortran interface
extern "C"
{
  void besselint_(double &b, double &qt, double &q2);
  double fort_besj0_(double &x);
  fcomplex s_(fcomplex &b);
  fcomplex alphasl_(fcomplex &q2);

  //exponents computed in alphasl
  extern struct {
    fcomplex aexp_;
    fcomplex aexpb_;
  } exponent_;
}

#endif
