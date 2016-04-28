#ifndef ctintegr_h
#define ctintegr_h
#include "cuba.h"

integrand_t ctintegrand   (const int &ndim, const double x[], const int &ncomp, double f[], void* userdata, const int &nvec, const int &core, double &weight, const int &iter);
integrand_t ctintegrandMC (const int &ndim, const double x[], const int &ncomp, double f[], void* userdata, const int &nvec, const int &core, double &weight, const int &iter);
integrand_t ctintegrand3d (const int &ndim, const double x[], const int &ncomp, double f[]);
integrand_t ctintegrand2d (const int &ndim, const double x[], const int &ncomp, double f[]);

int ctintegrand2d_cubature_v(unsigned ndim, long unsigned npts, const double x[], void *data, unsigned ncomp, double f[]);
int ctintegrand2d_cubature(unsigned ndim, const double x[], void *data, unsigned ncomp, double f[]);

extern "C" {
  double countint_(double r[22], double &wgt);
  double ctint_(double &costh, double &mm, double &qtt, double &yy, int &mode, double f[]);
}

#endif
