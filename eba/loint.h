#ifndef loint_h
#define loint_h

#include "mcfm_interface.h"

namespace loint
{
  void init();
  //  double lint(double costh, double m, double y, int mode, double f[]);
  void lint(double costh, double m, double y, int mode, double f[2]);
  double convolute(double fx1[2*MAXNF+1], double fx2[2*MAXNF+1], complex <double> mesq[12]);
}

#endif
