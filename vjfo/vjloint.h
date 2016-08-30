#ifndef vjloint_h
#define vjloint_h
#include "mcfm_interface.h"

namespace vjloint
{
  extern void init();
  extern double calc(const double x[5]);
  extern double calcvegas(const double x[7]);
  extern void fillp(double p[4][12]);
  extern double convolute(double fx1[2*MAXNF+1], double fx2[2*MAXNF+1], double msq[2*MAXNF+1][2*MAXNF+1]);

  extern double muf, mur;
}

#endif
