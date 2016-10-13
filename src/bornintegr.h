#ifndef bornintegr_h
#define bornintegr_h
#include "cuba.h"

//LO V integrands
integrand_t lointegrandMC(const int &ndim, const double x[], const int &ncomp, double f[], void* userdata, const int &nvec, const int &core, double &weight, const int &iter);
integrand_t lointegrand2d(const int &ndim, const double x[], const int &ncomp, double f[]);
int lointegrand2d_cubature_v(unsigned ndim, long unsigned npts, const double x[], void *data, unsigned ncomp, double f[]);
int lointegrand2d_cubature(unsigned ndim, const double x[], void *data, unsigned ncomp, double f[]);

extern "C" {
  double lowinthst_dynnlo_(double r[22], double &wgt, double f[]);
}

#endif
