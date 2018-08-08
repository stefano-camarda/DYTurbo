#include "cubacall.h"
#include "settings.h"
#include "resintegr.h"
#include "ctintegr.h"
#include "finintegr.h"
#include "bornintegr.h"
#include "cubature.h"
#include "smolpack.h"
#include "HistoHandler.h"

#include <cuba.h>
#include <iostream>

//flags for cuba Vegas integration:
//flags += 0 or 4; //collect only weights from final iteration (4) or from all iterations (0)
//flags = 8; //smoothing of importance sampling (0) or not (8)


/***************************************************************/
//resummation
void resintegr1d(vector <double> &res, double &err)
{
  const int ndim = 1;     //dimensions of the integral
  const int ncomp = 1;  //components of the integrand
  void *userdata = NULL;
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
  const int mineval = 65+2*65*opts.niterBORN;
  const int maxeval = 65+2*65*opts.niterBORN;
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
      double xmin[1] = {0};
      double xmax[1] = {1};
      if (opts.cubacores == 0)
	pcubature(ncomp, resintegrand1d_cubature, userdata,
		  ndim, xmin, xmax,
		  eval, epsabs, epsrel, ERROR_INDIVIDUAL, integral, error);
      else
	pcubature_v(ncomp, resintegrand1d_cubature_v, userdata,
		    ndim, xmin, xmax,
		    eval, epsabs, epsrel, ERROR_INDIVIDUAL, integral, error);
    }
  res.clear();
  res.push_back(integral[0]);
  for (int i = 1; i < opts.totpdf; i++)
    res.push_back(0);
  err = error[0];
  return;
}

void resintegr2d(vector <double> &res, double &err)
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
  double integral[ncomp];
  double error[ncomp];
  double prob[ncomp];
  const int flags = 0+opts.cubaverbosity;
  const int mineval = 65+2*65*opts.niterBORN;
  const int maxeval = 65+2*65*opts.niterBORN;
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
  res.clear();
  res.push_back(integral[0]);
  for (int i = 1; i < opts.totpdf; i++)
    res.push_back(0);
  err = error[0];
  return;
}

void resintegr3d(vector <double> &res, double &err)
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
  double integral[ncomp];
  double error[ncomp];
  double prob[ncomp];
  const int flags = 0+opts.cubaverbosity;
  const int mineval = 127+2*127*opts.niterBORN;
  const int maxeval = 127+2*127*opts.niterBORN;
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
  res.clear();
  res.push_back(integral[0]);
  for (int i = 1; i < opts.totpdf; i++)
    res.push_back(0);
  err = error[0];
  return;
}

//original integration, can use this to sample phase space and fill histograms
void resintegrMC(vector <double> &res, double &err)
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
  double integral[ncomp];
  double error[ncomp];
  double prob[ncomp];
  const int flags = 4+opts.cubaverbosity;
  const int seed = opts.rseed;
  const int mineval = opts.vegasncallsBORN;
  const int maxeval = opts.vegasncallsBORN;
  const int nstart = max(10, int(opts.vegasncallsBORN/10));
  const int nincrease = max(10, int(opts.vegasncallsBORN/10));
  const int nbatch = opts.cubanbatch;
  const int gridno = 0;
  Vegas(ndim, ncomp, (integrand_t)resintegrandMC, userdata, nvec,
	epsrel, epsabs,
	flags, seed,
	mineval, maxeval,
	nstart, nincrease, nbatch,
	gridno, statefile, spin,
	&neval, &fail,
	integral, error, prob);
  res.clear();
  res.push_back(integral[0]);
  for (int i = 1; i < opts.totpdf; i++)
    res.push_back(0);
  err = error[0];
  return;
}
/***************************************************************/

/***************************************************************/
// fixed order
//Cuba integration of analytical V+j
void vjintegr3d(vector <double> &res, double &err)
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
  double integral[ncomp];
  double error[ncomp];
  double prob[ncomp];
  const int flags = 0+opts.cubaverbosity;
  const int mineval = 127+2*127*opts.niterVJ;
  const int maxeval = 127+2*127*opts.niterVJ;
  const int key = 13;
  int nregions;
  if (!opts.pcubature)
    Cuhre(ndim, ncomp,
	  (integrand_t) vjintegrand, userdata, nvec,
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
      double xmin[3] = {0, 0., 0.};
      double xmax[3] = {1, 1., 1.};
      if (opts.cubacores == 0)
	pcubature(ncomp, vjintegrand_cubature, userdata,
		  //hcubature(ncomp, vjintegrand_cubature, userdata, 
		  ndim, xmin, xmax, 
		  eval, epsabs, epsrel, ERROR_INDIVIDUAL, integral, error);
      else
	pcubature_v(ncomp, vjintegrand_cubature_v, userdata, 
		    ndim, xmin, xmax, 
		    eval, epsabs, epsrel, ERROR_INDIVIDUAL, integral, error);
      /*
      //smolyak
      int print_stats = 0;
      int dim = ndim;
      int l = 10;
      integral[0] = int_smolyak (ndim, ndim+l, vjintegrand_smolyak, print_stats );
      error[0] = 0.000000001;
      */
    }
  res.clear();
  res.push_back(integral[0]);
  for (int i = 1; i < opts.totpdf; i++)
    res.push_back(0);
  err = error[0];
  return;
}

//Cuba integration of V+j LO
void vjlointegr(vector <double> &res, double &err)
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
  const int mineval   = opts.vegasncallsVJLO;
  const int maxeval   = opts.vegasncallsVJLO;
  const int nstart    = max(10, int(opts.vegasncallsVJLO/10));
  const int nincrease = max(10, int(opts.vegasncallsVJLO/10));
  const int nbatch    = opts.cubanbatch;
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

void vjlointegr7d(vector <double> &res, double &err)
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
  const int mineval   = opts.vegasncallsVJLO;
  const int maxeval   = opts.vegasncallsVJLO;
  const int nstart    = max(10, int(opts.vegasncallsVJLO/10));
  const int nincrease = max(10, int(opts.vegasncallsVJLO/10));
  const int nbatch    = opts.cubanbatch;
  const int gridno = 0;
  Vegas(ndim, ncomp, (integrand_t)vjlointegrandMC, userdata, nvec,
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

void vjlointegr5d(vector <double> &res, double &err)
{
  const int ndim = 4;//5;     //dimensions of the integral
  const int ncomp = (opts.helicity >= 0 ? 2 : 1);    //components of the integrand
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
  const int mineval = 153+2*153*opts.niterVJ;
  const int maxeval = 153+2*153*opts.niterVJ;
  const int key = 13;
  int nregions;
  if (!opts.pcubature)
    Cuhre(ndim, ncomp,
	  (integrand_t) vjlointegrand, userdata, nvec,
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
      double tiny = 0.;//1e-6;
      /*
      double xmin[5] = {0., 0., tiny,    0., 0.};
      double xmax[5] = {1., 1., 1.-tiny, 1., 1.};
      //if (opts.cubacores == 0)
      //pcubature(ncomp, vjlointegrand_cubature, userdata, //--> pcubature has an issue when phi is symmetric, it is resonant for nested quadrature rules, and pcubature misses the substructure of the phi distribution
      hcubature(ncomp, vjlointegrand_cubature, userdata, 
		ndim, xmin, xmax, 
		eval, epsabs, epsrel, ERROR_INDIVIDUAL, integral, error);
      */

      //      /*
      //4d integration (phi_lep integrated inside) works better in full phase space, and to calculate moments
      double xmin[4] = {0., 0., tiny,    0.};
      double xmax[4] = {1., 1., 1.-tiny, 1.};
      if (opts.cubacores == 0)
	//pcubature(ncomp, vjlointegrand_cubature, userdata, //--> pcubature has an issue when phi is symmetric, it is resonant for nested quadrature rules, and pcubature misses the substructure of the phi distribution
	hcubature(ncomp, vjlointegrand_cubature, userdata, 
		  4, xmin, xmax, 
		  eval, epsabs, epsrel, ERROR_LINF, integral, error);
      else
	//pcubature_v(ncomp, vjlointegrand_cubature_v, userdata, 
	hcubature_v(ncomp, vjlointegrand_cubature_v, userdata, 
		    4, xmin, xmax, 
		    eval, epsabs, epsrel, ERROR_LINF, integral, error);
      //      */

      /*
      //smolyak
      int print_stats = 0;
      int dim = 4;
      int l = 20;
      integral[0] = int_smolyak (ndim, ndim+l, vjlointegrand_smolyak, print_stats );
      error[0] = 0.000001;
      */
    }
  res.clear();
  if (opts.helicity >= 0)
    res.push_back(integral[0] != 0? integral[1]/integral[0] : 0);
  else
    res.push_back(integral[0]);
  for (int i = 1; i < opts.totpdf; i++)
    res.push_back(0);
  
  if (opts.helicity >= 0)
    err = integral[0] != 0? error[0]/integral[0]*integral[1]/integral[0] : 0; // error[1]/integral[1]*integral[1]/integral[0]
  else
    err = error[0];
  return;
}


//Cuba integration of the V+J NLO real part
void vjrealintegr(vector <double> &res, double &err)
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
  const int mineval = opts.vegasncallsVJREAL;
  const int maxeval = opts.vegasncallsVJREAL;
  const int nstart = max(10, int(opts.vegasncallsVJREAL/10));
  const int nincrease = max(10, int(opts.vegasncallsVJREAL/10));
  const int nbatch = opts.cubanbatch;
  const int gridno = 0;

  Vegas(ndim, ncomp, (integrand_t)realintegrand, userdata, nvec,
	epsrel, epsabs,
	flags, seed,
	mineval, maxeval,
	nstart, nincrease, nbatch,
	gridno, statefile, spin,
	&neval, &fail,
	integral, error, prob);

  /*
  //smolyak
  int print_stats = 0;
  int dim = 10;
  int l = 12;
  integral[0] = int_smolyak (ndim, ndim+l, realintegrand_smolyak, print_stats );
  error[0] = 0.000001;
  */
  
  res.clear();
  for (int i = 0; i < opts.totpdf; i++)
    res.push_back(integral[i]);
  err = error[0];
  return;
}

//Cuba integration of the V+J NLO virtual part
void vjvirtintegr(vector <double> &res, double &err)
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
  const int mineval = opts.vegasncallsVJVIRT;
  const int maxeval = opts.vegasncallsVJVIRT;
  const int nstart = max(10, int(opts.vegasncallsVJVIRT/10));
  const int nincrease = max(10, int(opts.vegasncallsVJVIRT/10));
  const int nbatch = opts.cubanbatch;
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

//Cuba integration of double virtual born configuration
void bornintegrMC6d(vector <double> &res, double &err)
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
  const int mineval   = opts.vegasncallsBORN;
  const int maxeval   = opts.vegasncallsBORN;
  const int nstart    = max(10, int(opts.vegasncallsBORN/10));
  const int nincrease = max(10, int(opts.vegasncallsBORN/10));
  const int nbatch    = opts.cubanbatch;
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
  return;
}

//Cuba integration of LO
void bornintegrMC4d(vector <double> &res, double &err)
{
  const int ndim = 4;     //dimensions of the integral
  const int ncomp = 1;    //components of the integrand
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
  const int mineval   = opts.vegasncallsBORN;
  const int maxeval   = opts.vegasncallsBORN;
  const int nstart    = max(10, int(opts.vegasncallsBORN/10));
  const int nincrease = max(10, int(opts.vegasncallsBORN/10));
  const int nbatch    = opts.cubanbatch;
  const int gridno = 0;
  Vegas(ndim, ncomp, (integrand_t)lointegrandMC, userdata, nvec,
	epsrel, epsabs,
	flags, seed,
	mineval, maxeval,
	nstart, nincrease, nbatch,
	gridno, statefile, spin,
	&neval, &fail,
	integral, error, prob);
  res.clear();
  res.push_back(integral[0]);
  for (int i = 1; i < opts.totpdf; i++)
    res.push_back(0);
  err = error[0];
  return;
}

//Cuba integration of LO
void bornintegr2d(vector <double> &res, double &err)
{
  const int ndim = 2;     //dimensions of the integral
  const int ncomp = (opts.helicity >= 0 ? 2 : 1);    //components of the integrand
  void *userdata = NULL;
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
  const int mineval = 65+2*65*opts.niterBORN;
  const int maxeval = 65+2*65*opts.niterBORN;
  const int key = 13;
  int nregions;

  if (!opts.pcubature)
    Cuhre(ndim, ncomp,
	  (integrand_t) lointegrand2d, userdata, nvec,
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
	pcubature(ncomp, lointegrand2d_cubature, userdata, 
		  ndim, xmin, xmax, 
		  eval, epsabs, epsrel, ERROR_INDIVIDUAL, integral, error);
      else
	pcubature_v(ncomp, lointegrand2d_cubature_v, userdata, 
		    ndim, xmin, xmax, 
		    eval, epsabs, epsrel, ERROR_INDIVIDUAL, integral, error);

      /*
      //smolyak
      int print_stats = 0;
      int dim = 2;
      int l = 10;
      integral[0] = int_smolyak (ndim, ndim+l, lointegrand2d_smolyak, print_stats );
      error[0] = 0.000001;
      */
    }
  res.clear();

  /*
  res.push_back(integral[0]);
  for (int i = 1; i < opts.totpdf; i++)
    res.push_back(0);
  err = error[0];
  */

  if (opts.helicity >= 0)
    res.push_back(integral[0] != 0? integral[1]/integral[0] : 0);
  else
    res.push_back(integral[0]);
  for (int i = 1; i < opts.totpdf; i++)
    res.push_back(0);
  
  if (opts.helicity >= 0)
    err = integral[0] != 0? error[0]/integral[0]*integral[1]/integral[0] : 0; // error[1]/integral[1]*integral[1]/integral[0]
  else
    err = error[0];
  return;

  
  return;
}

//Cuba integration of V + 2 jets
void v2jintegr(vector <double> &res, double &err)
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
  const int mineval = opts.vegasncallsVJREAL;
  const int maxeval = opts.vegasncallsVJREAL;
  const int nstart = max(10, int(opts.vegasncallsVJREAL/10));
  const int nincrease = max(10, int(opts.vegasncallsVJREAL/10));
  const int nbatch = opts.cubanbatch;
  const int gridno = 0;
  Vegas(ndim, ncomp, (integrand_t)v2jintegrand, userdata, nvec,
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
/***************************************************************/


/***************************************************************/
// Counter term
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
  const int nbatch = opts.cubanbatch;
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
  const int nbatch = opts.cubanbatch;
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

  //hists.FillQuadrature(res[0],err);
  return;
}

void ctintegr2d(vector <double> &res, double &err)
{
  const int ndim = 2;     //dimensions of the integral
  const int ncomp = opts.totpdf;  //components of the integrand
  void *userdata = NULL;
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
  if (!opts.pcubature)
    Cuhre(ndim, ncomp,
	  (integrand_t) ctintegrand2d, userdata, nvec,
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
      double xmin[2] = {0, 0.};
      double xmax[2] = {1, 1.};
      if (opts.cubacores == 0)
	pcubature(ncomp, ctintegrand2d_cubature, userdata, 
		  ndim, xmin, xmax, 
		  eval, epsabs, epsrel, ERROR_INDIVIDUAL, integral, error);
      else
	pcubature_v(ncomp, ctintegrand2d_cubature_v, userdata, 
		    ndim, xmin, xmax, 
		    eval, epsabs, epsrel, ERROR_INDIVIDUAL, integral, error);
    }

  res.clear();
  for (int i = 0; i < opts.totpdf; i++)
    res.push_back(integral[i]);
  err = error[0];

  //hists.FillQuadrature(res[0],err);
  return;
}

void ctintegr1d(vector <double> &res, double &err)
{
  const int ndim = 1;     //dimensions of the integral
  const int ncomp = opts.totpdf;  //components of the integrand
  void *userdata = NULL;
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
  if (!opts.pcubature)
    Cuhre(ndim, ncomp,
	  (integrand_t) ctintegrand1d, userdata, nvec,
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
      double xmin[1] = {0.};
      double xmax[1] = {1.};
      if (opts.cubacores == 0)
	pcubature(ncomp, ctintegrand1d_cubature, userdata,
		  ndim, xmin, xmax,
		  eval, epsabs, epsrel, ERROR_INDIVIDUAL, integral, error);
      else
	pcubature_v(ncomp, ctintegrand1d_cubature_v, userdata,
		    ndim, xmin, xmax,
		    eval, epsabs, epsrel, ERROR_INDIVIDUAL, integral, error);
    }

  res.clear();
  for (int i = 0; i < opts.totpdf; i++)
    res.push_back(integral[i]);
  err = error[0];

  //hists.FillQuadrature(res[0],err);
  return;
}
/***************************************************************/

void initfun(void * input, const int &core){
    if(core != PARENT_PROC) HistoHandler::Reset();
}


void exitfun(void * input, const int &core){
    HistoHandler::Save(core);
}

void tell_to_grid_we_are_alive(){
  if(opts.verbose && ICALL % 100000==0) 
      printf (" Hi Grid, we are sitll alive! Look, our event is %d\n",ICALL);
  ICALL++;
}
