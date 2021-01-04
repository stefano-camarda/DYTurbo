#ifndef resintegr_h
#define resintegr_h
#include "cuba.h"

//extern double const qtcutoff;

integrand_t resintegrand1d(const int &ndim, const double x[], const int &ncomp, double f[]);
integrand_t resintegrand2d(const int &ndim, const double x[], const int &ncomp, double f[]);
integrand_t resintegrand2d_my(const int &ndim, const double x[], const int &ncomp, double f[]);
integrand_t resintegrand3d(const int &ndim, const double x[], const int &ncomp, double f[]);
integrand_t resintegrandMC(const int &ndim, const double x[], const int &ncomp, double f[],
			   void* userdata, const int &nvec, const int &core,
			   double &weight, const int &iter);

int resintegrand1d_cubature(unsigned ndim, const double x[], void *data, unsigned ncomp, double f[]);
int resintegrand1d_cubature_v(unsigned ndim, long unsigned npts, const double x[], void *data, unsigned ncomp, double f[]);

int resintegrand2d_cubature(unsigned ndim, const double x[], void *data, unsigned ncomp, double f[]);
int resintegrand2d_cubature_v(unsigned ndim, long unsigned npts, const double x[], void *data, unsigned ncomp, double f[]);

int resintegrand2d_my_cubature(unsigned ndim, const double x[], void *data, unsigned ncomp, double f[]);
int resintegrand2d_my_cubature_v(unsigned ndim, long unsigned npts, const double x[], void *data, unsigned ncomp, double f[]);

int resintegrand3d_cubature(unsigned ndim, const double x[], void *data, unsigned ncomp, double f[]);
int resintegrand3d_cubature_v(unsigned ndim, long unsigned npts, const double x[], void *data, unsigned ncomp, double f[]);

#endif
