#include <cuba.h>

#include "cubacall.h"
#include "settings.h"
#include "integr.h"
#include "finintegr.h"

void integr2d(double &res, double &err)
{
  const int ndim = 2;     //dimensions of the integral
  const int ncomp = 1;  //components of the integrand
  void *userdata;
  const int nvec = 1;
  const double epsrel = 0.;
  const double epsabs = 0.;
  const char *statefile = "";
  void *spin=NULL;
  int neval;
  int fail;
  double integral[1];
  double error[1];
  double prob[1];
  const int flags = 0+opts.cubaverbosity;
  const int mineval = 65+2*65*opts.niter;
  const int maxeval = 65+2*65*opts.niter;
  const int key = 13;
  int nregions;
  Cuhre(ndim, ncomp,
	(integrand_t) resintegrand2d, userdata, nvec,
	epsrel, epsabs,
	flags,
	mineval, maxeval,
	key, statefile, NULL,
	&nregions, &neval, &fail,
  	integral, error, prob);

  res = integral[0];
  err = error[0];
  return;
}

void integr3d(double &res, double &err)
{
  const int ndim = 3;     //dimensions of the integral
  const int ncomp = 1;  //components of the integrand
  void *userdata;
  const int nvec = 1;
  const double epsrel = 0.;
  const double epsabs = 0.;
  const char *statefile = "";
  void *spin=NULL;
  int neval;
  int fail;
  double integral[1];
  double error[1];
  double prob[1];
  const int flags = 0+opts.cubaverbosity;
  const int mineval = 127+2*127*opts.niter;
  const int maxeval = 127+2*127*opts.niter;
  const int key = 13;
  int nregions;
  Cuhre(ndim, ncomp,
	(integrand_t) resintegrand3d, userdata, nvec,
	epsrel, epsabs,
	flags,
	mineval, maxeval,
	key, statefile, NULL,
	&nregions, &neval, &fail,
  	integral, error, prob);

  res = integral[0];
  err = error[0];
  return;
}

//original integration, can use this to sample phase space and fill histograms
void integr4d(double &res, double &err)
{
  const int ndim = 4;   //dimensions of the integral
  const int ncomp = 1;  //components of the integrand
  void *userdata;
  const int nvec = 1;
  const double epsrel = 0.;
  const double epsabs = 0.;
  const char *statefile = "";
  void *spin=NULL;
  int neval;
  int fail;
  double integral[1];
  double error[1];
  double prob[1];
  const int flags = 4+opts.cubaverbosity;
  const int seed = 1;
  const int mineval = opts.vegasncalls;
  const int maxeval = opts.vegasncalls;
  const int nstart = 1000;
  const int nincrease = 1000;
  const int nbatch = 10000;
  const int gridno = 1;
  Vegas(ndim, ncomp, (integrand_t)resintegrand4d, userdata, nvec,
	epsrel, epsabs,
	flags, seed,
	mineval, maxeval,
	nstart, nincrease, nbatch,
	gridno, statefile, NULL,
	&neval, &fail,
	integral, error, prob);
  res = integral[0];
  err = error[0];
  return;
}

//Cuba integration of Z+j LO
void lowintegr(double &res, double &err)
{
  const int ndim = 7;   //dimensions of the integral
  const int ncomp = 1;  //components of the integrand
  void *userdata;
  const int nvec = 1;
  const double epsrel = 0.;
  const double epsabs = 0.;
  const char *statefile = "";
  void *spin=NULL;
  int neval;
  int fail;
  double integral[1];
  double error[1];
  double prob[1];
  const int flags = 8+4+opts.cubaverbosity;
  const int seed = 1;
  const int mineval = 10000000;
  const int maxeval = 10000000;
  const int nstart = 100000;
  const int nincrease = 100000;
  const int nbatch = 1000;
  const int gridno = 1;
  Vegas(ndim, ncomp, (integrand_t)lowintegrand, userdata, nvec,
	epsrel, epsabs,
	flags, seed,
	mineval, maxeval,
	nstart, nincrease, nbatch,
	gridno, statefile, spin,
	&neval, &fail,
	integral, error, prob);
  res = integral[0];
  err = error[0];

  return;
}

//Cuba integration of the real part
void realintegr(double &res, double &err)
{
  const int ndim = 10;   //dimensions of the integral
  const int ncomp = 1;  //components of the integrand
  void *userdata;
  const int nvec = 1;
  const double epsrel = 0.;
  const double epsabs = 0.;
  const char *statefile = "";
  void *spin=NULL;
  int neval;
  int fail;
  double integral[1];
  double error[1];
  double prob[1];
  const int flags = 8+4+opts.cubaverbosity;
  const int seed = 1;
  const int mineval = 1000000;
  const int maxeval = 1000000;
  const int nstart = 10000;
  const int nincrease = 10000;
  const int nbatch = 1000;
  const int gridno = 1;
  Vegas(ndim, ncomp, (integrand_t)realintegrand, userdata, nvec,
	epsrel, epsabs,
	flags, seed,
	mineval, maxeval,
	nstart, nincrease, nbatch,
	gridno, statefile, spin,
	&neval, &fail,
	integral, error, prob);
  res = integral[0];
  err = error[0];

  return;
}

//Cuba integration of the virtual part
void virtintegr(double &res, double &err)
{
  const int ndim = 8;   //dimensions of the integral
  const int ncomp = 1;  //components of the integrand
  void *userdata;
  const int nvec = 1;
  const double epsrel = 0.;
  const double epsabs = 0.;
  const char *statefile = "";
  void *spin=NULL;
  int neval;
  int fail;
  double integral[1];
  double error[1];
  double prob[1];
  const int flags = 8+4+opts.cubaverbosity;
  const int seed = 1;
  const int mineval = 10000000;
  const int maxeval = 10000000;
  const int nstart = 100000;
  const int nincrease = 100000;
  const int nbatch = 1000;
  const int gridno = 1;
  Vegas(ndim, ncomp, (integrand_t)virtintegrand, userdata, nvec,
	epsrel, epsabs,
	flags, seed,
	mineval, maxeval,
	nstart, nincrease, nbatch,
	gridno, statefile, spin,
	&neval, &fail,
	integral, error, prob);
  res = integral[0];
  err = error[0];

  return;
}

//Cuba integration of the counterterm
void ctintegr(double &res, double &err)
{
  const int ndim = 8;   //dimensions of the integral
  const int ncomp = 1;  //components of the integrand
  void *userdata;
  const int nvec = 1;
  const double epsrel = 0.;
  const double epsabs = 0.;
  const char *statefile = "";
  void *spin=NULL;
  int neval;
  int fail;
  double integral[1];
  double error[1];
  double prob[1];
  const int flags = 8+4+opts.cubaverbosity;
  const int seed = 1;
  const int mineval = 1000000;
  const int maxeval = 1000000;
  const int nstart = 100000;
  const int nincrease = 100000;
  const int nbatch = 1000;
  const int gridno = 1;
  Vegas(ndim, ncomp, (integrand_t)ctintegrand, userdata, nvec,
	epsrel, epsabs,
	flags, seed,
	mineval, maxeval,
	nstart, nincrease, nbatch,
	gridno, statefile, spin,
	&neval, &fail,
	integral, error, prob);
  res = integral[0];
  err = error[0];

  return;
}

