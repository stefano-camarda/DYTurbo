#ifndef vjloint_h
#define vjloint_h

namespace vjloint
{
  extern void init();
  extern double calc(const double x[5]);
  extern void fillp(double p[4][12]);
  extern int x2int(unsigned ndim, const double x[], void *data, unsigned ncomp, double f[]);

  extern double muf, mur;
}

#endif
