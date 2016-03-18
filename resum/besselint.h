#ifndef besselint_h
#define besselint_h

#include "interface.h"

namespace besselint
{
  extern double bint(double b, double qt, double q2);
}

//fortran interface
extern "C"
{
  void besselint_(double &b, double &qt, double &q2);
  double fort_besj0_(double &x);
  fcomplex s_(fcomplex &b);
  fcomplex alphasl_(fcomplex &q2);
}

#endif
