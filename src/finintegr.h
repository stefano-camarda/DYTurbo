#ifndef finintegr_h
#define finintegr_h
#include "cuba.h"

integrand_t vjintegrand(const int &ndim, const double x[], const int &ncomp, double f[]);
integrand_t lowintegrand  (const int &ndim, const double x[], const int &ncomp, double f[], void* userdata, const int &nvec, const int &core, double &weight, const int &iter);
integrand_t realintegrand (const int &ndim, const double x[], const int &ncomp, double f[], void* userdata, const int &nvec, const int &core, double &weight, const int &iter);
integrand_t virtintegrand (const int &ndim, const double x[], const int &ncomp, double f[], void* userdata, const int &nvec, const int &core, double &weight, const int &iter);
integrand_t doublevirtintegrand(const int &ndim, const double x[], const int &ncomp, double f[], void* userdata, const int &nvec, const int &core, double &weight, const int &iter);

extern "C" {
  double lowint_(double r[22], double &wgt, double f[]);
  double realint_(double r[22], double &wgt, double f[]);
  double virtint_(double r[22], double &wgt, double f[]);
  double lowinthst_(double r[22], double &wgt);
  double vjfo_(double &m, double &pt, double &y);
}

#endif
