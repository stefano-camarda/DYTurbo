#ifndef finintegr_h
#define finintegr_h
#include "cuba.h"

integrand_t lowintegrand(const int &ndim, const double x[], const int &ncomp, double f[]);
integrand_t realintegrand(const int &ndim, const double x[], const int &ncomp, double f[]);
integrand_t virtintegrand(const int &ndim, const double x[], const int &ncomp, double f[]);
integrand_t ctintegrand(const int &ndim, const double x[], const int &ncomp, double f[]);
integrand_t ctintegrand3d(const int &ndim, const double x[], const int &ncomp, double f[]);

integrand_t ctintegrandab(const int &ndim, const double x[], const int &ncomp, double f[]);

#endif
