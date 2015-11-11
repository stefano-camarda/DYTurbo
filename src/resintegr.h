#ifndef resintegr_h
#define resintegr_h
#include "cuba.h"

integrand_t resintegrand2d(const int &ndim, const double x[], const int &ncomp, double f[]);
integrand_t resintegrand3d(const int &ndim, const double x[], const int &ncomp, double f[]);
integrand_t resintegrand4d(const int &ndim, const double x[], const int &ncomp, double f[],
			   void* userdata, const int &nvec, const int &core,
			   double &weight, const int &iter);
#endif