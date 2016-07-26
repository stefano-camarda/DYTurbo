#ifndef finintegr_h
#define finintegr_h
#include "cuba.h"

integrand_t vjintegrand(const int &ndim, const double x[], const int &ncomp, double f[]);
integrand_t lowintegrand  (const int &ndim, const double x[], const int &ncomp, double f[], void* userdata, const int &nvec, const int &core, double &weight, const int &iter);
integrand_t realintegrand (const int &ndim, const double x[], const int &ncomp, double f[], void* userdata, const int &nvec, const int &core, double &weight, const int &iter);
integrand_t virtintegrand (const int &ndim, const double x[], const int &ncomp, double f[], void* userdata, const int &nvec, const int &core, double &weight, const int &iter);
integrand_t doublevirtintegrand(const int &ndim, const double x[], const int &ncomp, double f[], void* userdata, const int &nvec, const int &core, double &weight, const int &iter);
integrand_t v2jintegrand(const int &ndim, const double x[], const int &ncomp, double f[], void* userdata, const int &nvec, const int &core, double &weight, const int &iter);
integrand_t lointegrand(const int &ndim, const double x[], const int &ncomp, double f[], void* userdata, const int &nvec, const int &core, double &weight, const int &iter);
integrand_t lointegrand2d(const int &ndim, const double x[], const int &ncomp, double f[]);

int vjintegrand_cubature_v(unsigned ndim, long unsigned npts, const double x[], void *data, unsigned ncomp, double f[]);
int vjintegrand_cubature(unsigned ndim, const double x[], void *data, unsigned ncomp, double f[]);

int lointegrand2d_cubature_v(unsigned ndim, long unsigned npts, const double x[], void *data, unsigned ncomp, double f[]);
int lointegrand2d_cubature(unsigned ndim, const double x[], void *data, unsigned ncomp, double f[]);

extern "C" {
  double lowint_(double r[22], double &wgt, double f[]);
  double realint_(double r[22], double &wgt, double f[]);
  double virtint_(double r[22], double &wgt, double f[]);
  double lowinthst_dynnlo_(double r[22], double &wgt, double f[]);
  double vjfo_(double &m, double &pt, double &y);
  double v2jint_(double r[22], double &wgt, double f[]);
}

#endif
