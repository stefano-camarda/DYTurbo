#include "cubacall.h"
#include "settings.h"
#include "resintegr.h"
#include "ctintegr.h"
#include "finintegr.h"
#include "plotter.h"
#include "cubature.h"

#include <cuba.h>
#include <iostream>

void resintegr2d(double &res, double &err)
{
  const int ndim = 2;     //dimensions of the integral
  const int ncomp = 1;  //components of the integrand
  void *userdata = NULL;
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

  if (!opts.pcubature)
    Cuhre(ndim, ncomp,
	  (integrand_t) resintegrand2d, userdata, nvec,
	  epsrel, epsabs,
	  flags,
	  mineval, maxeval,
	  key, statefile, NULL,
	  &nregions, &neval, &fail,
	  integral, error, prob);
  else
    {
      const int eval = 0;
      const double epsrel = opts.pcubaccuracy;
      const double epsabs = 0.;
      double xmin[2] = {0, 0};
      double xmax[2] = {1, 1};
      if (opts.cubacores == 0)
	pcubature(ncomp, resintegrand2d_cubature, userdata, 
		  ndim, xmin, xmax, 
		  eval, epsabs, epsrel, ERROR_INDIVIDUAL, integral, error);
      else
	pcubature_v(ncomp, resintegrand2d_cubature_v, userdata, 
		    ndim, xmin, xmax, 
		    eval, epsabs, epsrel, ERROR_INDIVIDUAL, integral, error);
    }
  
  res = integral[0];
  err = error[0];
  hists.FillQuadrature(res,err);
  return;
}

void resintegr3d(double &res, double &err)
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
  hists.FillQuadrature(res,err);
  return;
}

//original integration, can use this to sample phase space and fill histograms
void resintegrMC(double &res, double &err)
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
  const int flags = 4+opts.cubaverbosity;
  const int seed = opts.rseed;
  const int mineval = opts.vegasncallsRES;
  const int maxeval = opts.vegasncallsRES;
  const int nstart = max(10, int(opts.vegasncallsRES/10));
  const int nincrease = max(10, int(opts.vegasncallsRES/10));
  const int nbatch = 1000;
  const int gridno = 0;
  Vegas(ndim, ncomp, (integrand_t)resintegrandMC, userdata, nvec,
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

//Cuba integration of analytical Z+j
void vjintegr3d(double &res, double &err)
{
  const int ndim = 3;     //dimensions of the integral
  const int ncomp = 1;    //components of the integrand
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
  const int mineval = 127+2*127*opts.niterVJ;
  const int maxeval = 127+2*127*opts.niterVJ;
  const int key = 13;
  int nregions;
  Cuhre(ndim, ncomp,
	(integrand_t) vjintegrand, userdata, nvec,
	epsrel, epsabs,
	flags,
	mineval, maxeval,
	key, statefile, NULL,
	&nregions, &neval, &fail,
  	integral, error, prob);

  res = integral[0];
  err = error[0];
  hists.FillQuadrature(res,err);
  return;
}


//Cuba integration of Z+j LO
void lowintegr(vector <double> &res, double &err)
{
  const int ndim = 7;   //dimensions of the integral
  const int ncomp = opts.totpdf;  //components of the integrand
  void *userdata;
  const int nvec = 1;
  const double epsrel = 0.;
  const double epsabs = 0.;
  const char *statefile = "";
  void *spin=NULL;
  int neval;
  int fail;
  double integral[ncomp];
  double error[ncomp];
  double prob[ncomp];
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
  res.clear();
  for (int i = 0; i < opts.totpdf; i++)
    res.push_back(integral[i]);
  err = error[0];

  return;
}

//Cuba integration of the real part
void realintegr(vector <double> &res, double &err)
{
  const int ndim = 10;   //dimensions of the integral
  const int ncomp = opts.totpdf;  //components of the integrand
  void *userdata;
  const int nvec = 1;
  const double epsrel = 0.;
  const double epsabs = 0.;
  const char *statefile = "";
  void *spin=NULL;
  int neval;
  int fail;
  double integral[ncomp];
  double error[ncomp];
  double prob[ncomp];
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
  res.clear();
  for (int i = 0; i < opts.totpdf; i++)
    res.push_back(integral[i]);
  err = error[0];

  return;
}

//Cuba integration of the virtual part
void virtintegr(vector <double> &res, double &err)
{
  const int ndim = 8;   //dimensions of the integral
  const int ncomp = opts.totpdf;  //components of the integrand
  void *userdata;
  const int nvec = 1;
  const double epsrel = 0.;
  const double epsabs = 0.;
  const char *statefile = "";
  void *spin=NULL;
  int neval;
  int fail;
  double integral[ncomp];
  double error[ncomp];
  double prob[ncomp];
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
  res.clear();
  for (int i = 0; i < opts.totpdf; i++)
    res.push_back(integral[i]);
  err = error[0];

  return;
}

//Cuba integration of double virtual
void doublevirtintegr(vector <double> &res, double &err)
{
  const int ndim = 6;   //dimensions of the integral
  const int ncomp = opts.totpdf;  //components of the integrand
  void *userdata;
  const int nvec = 1;
  const double epsrel = 0.;
  const double epsabs = 0.;
  const char *statefile = "";
  void *spin=NULL;
  int neval;
  int fail;
  double integral[ncomp];
  double error[ncomp];
  double prob[ncomp];
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
  res.clear();
  for (int i = 0; i < opts.totpdf; i++)
    res.push_back(integral[i]);
  err = error[0];
  hists.FillQuadrature(res[0],err);

  return;
}

//Cuba integration of the counterterm
void ctintegr(vector <double> &res, double &err)
{
  const int ndim = 8;   //dimensions of the integral
  const int ncomp = opts.totpdf;  //components of the integrand
  void *userdata;
  const int nvec = 1;
  const double epsrel = 0.;
  const double epsabs = 0.;
  const char *statefile = "";
  void *spin=NULL;
  int neval;
  int fail;
  double integral[ncomp];
  double error[ncomp];
  double prob[ncomp];
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
  res.clear();
  for (int i = 0; i < opts.totpdf; i++)
    res.push_back(integral[i]);
  err = error[0];

  return;
}

//Cuba integration of the counterterm
void ctintegrMC(vector <double> &res, double &err)
{
  const int ndim = 6;   //dimensions of the integral
  const int ncomp = opts.totpdf;  //components of the integrand
  void *userdata;
  const int nvec = 1;
  const double epsrel = 0.;
  const double epsabs = 0.;
  const char *statefile = "";
  void *spin=NULL;
  int neval;
  int fail;
  double integral[ncomp];
  double error[ncomp];
  double prob[ncomp];
  const int flags = 8+4+opts.cubaverbosity;
  const int seed = opts.rseed;
  const int mineval = opts.vegasncallsCT;
  const int maxeval = opts.vegasncallsCT;
  const int nstart = max(10, int(opts.vegasncallsCT/10));
  const int nincrease = max(10, int(opts.vegasncallsCT/10));
  const int nbatch = 1000;
  const int gridno = 0;
  Vegas(ndim, ncomp, (integrand_t)ctintegrandMC, userdata, nvec,
	epsrel, epsabs,
	flags, seed,
	mineval, maxeval,
	nstart, nincrease, nbatch,
	gridno, statefile, spin,
	&neval, &fail,
	integral, error, prob);
  res.clear();
  for (int i = 0; i < opts.totpdf; i++)
    res.push_back(integral[i]);
  err = error[0];

  return;
}

void ctintegr3d(vector <double> &res, double &err)
{
  const int ndim = 3;     //dimensions of the integral
  const int ncomp = opts.totpdf;  //components of the integrand
  void *userdata;
  const int nvec = 1;
  const double epsrel = 0.;
  const double epsabs = 0.;
  const char *statefile = "";
  void *spin=NULL;
  int neval;
  int fail;
  double integral[ncomp];
  double error[ncomp];
  double prob[ncomp];
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
  res.clear();
  for (int i = 0; i < opts.totpdf; i++)
    res.push_back(integral[i]);
  err = error[0];

  hists.FillQuadrature(res[0],err);
  return;
}

void ctintegr2d(vector <double> &res, double &err)
{
  const int ndim = 2;     //dimensions of the integral
  const int ncomp = opts.totpdf;  //components of the integrand
  void *userdata;
  const int nvec = 1;
  const double epsrel = 0.;
  const double epsabs = 0.;
  const char *statefile = "";
  void *spin=NULL;
  int neval;
  int fail;
  double integral[ncomp];
  double error[ncomp];
  double prob[ncomp];
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
  res.clear();
  for (int i = 0; i < opts.totpdf; i++)
    res.push_back(integral[i]);
  err = error[0];

  hists.FillQuadrature(res[0],err);
  return;
}

void exitfun(void * input, const int &core){
    if (opts.cubacores!=0) hists.Finalise(core);
}
