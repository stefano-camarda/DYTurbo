#ifndef finintegr_h
#define finintegr_h
#include "cuba.h"

integrand_t vjintegrand(const int &ndim, const double x[], const int &ncomp, double f[]);
integrand_t lowintegrand  (const int &ndim, const double x[], const int &ncomp, double f[], void* userdata, const int &nvec, const int &core, double &weight, const int &iter);
integrand_t realintegrand (const int &ndim, const double x[], const int &ncomp, double f[], void* userdata, const int &nvec, const int &core, double &weight, const int &iter);
integrand_t virtintegrand (const int &ndim, const double x[], const int &ncomp, double f[], void* userdata, const int &nvec, const int &core, double &weight, const int &iter);
integrand_t ctintegrand   (const int &ndim, const double x[], const int &ncomp, double f[], void* userdata, const int &nvec, const int &core, double &weight, const int &iter);
integrand_t ctintegrandMC(const int &ndim, const double x[], const int &ncomp, double f[], void* userdata, const int &nvec, const int &core, double &weight, const int &iter);
integrand_t ctintegrand3d (const int &ndim, const double x[], const int &ncomp, double f[]);
integrand_t ctintegrand2d (const int &ndim, const double x[], const int &ncomp, double f[]);
integrand_t doublevirtintegrand(const int &ndim, const double x[], const int &ncomp, double f[], void* userdata, const int &nvec, const int &core, double &weight, const int &iter);

#endif
