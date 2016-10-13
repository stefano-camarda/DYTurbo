#include "bornintegr.h"
#include "interface.h"
#include "phasespace.h"
#include "settings.h"
#include "loint.h"
#include "cubacall.h"
//#include "isnan.h"

//#include <iostream>
//#include <iomanip>
//#include <ctime>
//#include <math.h>


#include "histo/KinematicCuts.h"

//using namespace std;

integrand_t doublevirtintegrand(const int &ndim, const double x[], const int &ncomp, double f[],
                        void* userdata, const int &nvec, const int &core,
                        double &weight, const int &iter)
{
  
  double rvv[22];
  for (int i = 0; i < ndim; i++)
    rvv[i]=x[i];

  dofill_.doFill_ = int(iter==last_iter);
  f[0] = lowinthst_dynnlo_(rvv,weight,f);

  tell_to_grid_we_are_alive();
  return 0;
}

integrand_t lointegrandMC(const int &ndim, const double x[], const int &ncomp, double f[],
			 void* userdata, const int &nvec, const int &core,
			 double &weight, const int &iter)
{
  clock_t begin_time, end_time;

  begin_time = clock();

  if (opts.fixedorder && phasespace::qtmin > 0)
    {
      f[0]=0.;
      return 0;
    }
  
  //Jacobian of the change of variables from the unitary hypercube x[3] to the m, y boundaries
  double jac = 1.;
  bool status = true;

  //generate phase space
  double r1[2] = {x[0], x[1]};
  status = phasespace::gen_my(r1, jac);
  if (!status)
    {
      f[0] = 0.;
      return 0;
    }
  phasespace::set_qt(0.);
  phasespace::set_phiV(0.);
  
  double r2[2] = {x[2], x[3]};
  phasespace::gen_costhphi(r2, jac);

  //generate boson and leptons
  phasespace::calcexpy();
  phasespace::calcmt();
  phasespace::genV4p();
  phasespace::genl4p();

  //apply cuts
  if (!Kinematics::Cuts::KeepThisEvent(phasespace::p3,phasespace::p4))
    {
      f[0] = 0.;
      return 0;
    }
  
  //calculate integrand
  f[0] = loint::lint(phasespace::costh,phasespace::m,phasespace::y,0)*jac;

  //fill histograms
  if (iter==last_iter)
    {
      double wt = weight*f[0];
      hists_fill_(phasespace::p3,phasespace::p4,&wt);
    }
  
  tell_to_grid_we_are_alive();
  return 0;
}

int lointegrand2d_cubature_v(unsigned ndim, long unsigned npts, const double x[], void *data, unsigned ncomp, double f[])
{
#pragma omp parallel for num_threads(opts.cubacores) copyin(scale_,facscale_,qcdcouple_)
  for (unsigned i = 0; i < npts; i++)
    {
      // evaluate the integrand for npts points
      double xi[ndim];
      double fi[ncomp];
      for (unsigned j = 0; j < ndim; j++)
	xi[j] = x[i*ndim + j];

      lointegrand2d(ndim, xi, ncomp, fi);
      
      for (unsigned k = 0; k < ncomp; ++k)
	f[i*ncomp + k] = fi[k];
    }
  tell_to_grid_we_are_alive();
  return 0;
}

int lointegrand2d_cubature(unsigned ndim, const double x[], void *data, unsigned ncomp, double f[])
{
  lointegrand2d(ndim, x, ncomp, f);
  tell_to_grid_we_are_alive();
  return 0;
}

integrand_t lointegrand2d(const int &ndim, const double x[], const int &ncomp, double f[])
{
  clock_t begin_time, end_time;

  begin_time = clock();

  if (opts.fixedorder && phasespace::qtmin > 0)
    {
      f[0]=0.;
      return 0;
    }

  //Jacobian of the change of variables from the unitary hypercube x[3] to the m, y boundaries
  double jac = 1.;
  bool status = true;

  double r1[2] = {x[0], x[1]};
  status = phasespace::gen_my(r1, jac);
  if (!status)
    {
      f[0] = 0.;
      return 0;
    }
  phasespace::set_qt(0.);

  double m = phasespace::m;
  double y = phasespace::y;

  f[0] = loint::lint(0.,phasespace::m,phasespace::y,1)*jac;

  //perform phi integration
  f[0] = f[0] * 2.*M_PI;
  
  tell_to_grid_we_are_alive();
  return 0;
}
