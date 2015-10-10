#include <cuba.h>

#include "cubacall.h"
#include "settings.h"
#include "integr.h"
#include "finintegr.h"
#include "resintegr.h"
#include "plotter.h"

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
  const int mineval = 65+2*65*opts.niterRES;
  const int maxeval = 65+2*65*opts.niterRES;
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
  const int mineval = 127+2*127*opts.niterRES;
  const int maxeval = 127+2*127*opts.niterRES;
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
  const int seed = opts.rseed;
  const int mineval = opts.vegasncallsRES;
  const int maxeval = opts.vegasncallsRES;
  const int nstart = max(10, int(opts.vegasncallsRES/10));
  const int nincrease = max(10, int(opts.vegasncallsRES/10));
  const int nbatch = 1000;
  const int gridno = 0;
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
  const int seed = opts.rseed;
  const int mineval   = opts.vegasncallsLO;
  const int maxeval   = opts.vegasncallsLO;
  const int nstart    = max(10, int(opts.vegasncallsLO/10));
  const int nincrease = max(10, int(opts.vegasncallsLO/10));
  const int nbatch    = 1000;
  const int gridno = 0;
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
  const int seed = opts.rseed;
  const int mineval = opts.vegasncallsREAL;
  const int maxeval = opts.vegasncallsREAL;
  const int nstart = max(10, int(opts.vegasncallsREAL/10));
  const int nincrease = max(10, int(opts.vegasncallsREAL/10));
  const int nbatch = 1000;
  const int gridno = 0;
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
  const int seed = opts.rseed;
  const int mineval = opts.vegasncallsVIRT;
  const int maxeval = opts.vegasncallsVIRT;
  const int nstart = max(10, int(opts.vegasncallsVIRT/10));
  const int nincrease = max(10, int(opts.vegasncallsVIRT/10));
  const int nbatch = 1000;
  const int gridno = 0;
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
  const int seed = opts.rseed;
  const int mineval = opts.vegasncallsCT;
  const int maxeval = opts.vegasncallsCT;
  const int nstart = max(10, int(opts.vegasncallsCT/10));
  const int nincrease = max(10, int(opts.vegasncallsCT/10));
  const int nbatch = 1000;
  const int gridno = 0;
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

void ctintegr3d(double &res, double &err)
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
  const int mineval = 127+2*127*opts.niterCT;
  const int maxeval = 127+2*127*opts.niterCT;
  const int key = 13;
  int nregions;
  Cuhre(ndim, ncomp,
	(integrand_t) ctintegrand3d, userdata, nvec,
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

void ctintegr2d(double &res, double &err)
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
  const int mineval = 65+2*65*opts.niterCT;
  const int maxeval = 65+2*65*opts.niterCT;
  const int key = 13;
  int nregions;
  Cuhre(ndim, ncomp,
	(integrand_t) ctintegrand2d, userdata, nvec,
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

//Cuba integration of double virtual
void doublevirtintegr(double &res, double &err)
{
  const int ndim = 6;   //dimensions of the integral
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
  const int seed = opts.rseed;
  const int mineval   = opts.vegasncallsVV;
  const int maxeval   = opts.vegasncallsVV;
  const int nstart    = max(10, int(opts.vegasncallsVV/10));
  const int nincrease = max(10, int(opts.vegasncallsVV/10));
  const int nbatch    = 1000;
  const int gridno = 0;
  Vegas(ndim, ncomp, (integrand_t)doublevirtintegrand, userdata, nvec,
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

void exitfun(void * input, const int &core){
    hists.Finalise(core);
    //hists.Merge();
    //printf(" WORKER exit process %d\n", core);
}
