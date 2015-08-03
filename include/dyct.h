#ifndef dyct_h
#define dyct_h
#include "cuba.h"

extern double dyct(double m, double y, double qt, double phicm, double phiZ, double cos_th, double alpha, double beta);

integrand_t ctintegrand(const int &ndim, const double x[], const int &ncomp, double f[]);


#endif
