#include "interface.h"
#include "finintegr.h"
#include "integr.h"
#include "settings.h"
#include "vjint.h"

#include <iostream>
#include <iomanip>
#include <ctime>
#include <math.h>


using namespace std;

int vjintegrand_cubature_v(unsigned ndim, long unsigned npts, const double x[], void *data, unsigned ncomp, double f[])
{
  //  cout << "parallel " << npts << endl;
#pragma omp parallel for num_threads(opts.cubacores) copyin(scale_,facscale_,qcdcouple_)
  for (unsigned i = 0; i < npts; i++)
    {
      // evaluate the integrand for npts points
      double xi[ndim];
      double fi[ncomp];
      for (unsigned j = 0; j < ndim; j++)
	xi[j] = x[i*ndim + j];

      vjintegrand(ndim, xi, ncomp, fi);
      
      for (unsigned k = 0; k < ncomp; ++k)
	f[i*ncomp + k] = fi[k];
    }
  tell_to_grid_we_are_alive();
  return 0;
}

int vjintegrand_cubature(unsigned ndim, const double x[], void *data, unsigned ncomp, double f[])
{
  vjintegrand(ndim, x, ncomp, f);
  tell_to_grid_we_are_alive();
  return 0;
}


integrand_t vjintegrand(const int &ndim, const double x[], const int &ncomp, double f[])
//Generates the phase space 4 vectors
//Calculates the V+j LO/NLO integrand as a function of m, qt, and y
//dOmega integration is already performed in the integrand
//works only without cuts on the leptons
{
  //  cout << x[0] << "  " << x[1] << "  " << x[2] << endl;
  clock_t begin_time, end_time;

  begin_time = clock();

  //Jacobian of the change of variables from the unitary square x[2] to the m, qt boundaries
  double jac = 1.;

  // Generate the boson invariant mass between the integration boundaries
  double mcut = qtmax/qtcut_.xqtcut_;
  double wsqmin = pow(mmin,2);
  double wsqmax = pow(min(mmax,mcut),2);
  if (wsqmin >= wsqmax)
    {
      f[0] = 0.;
      return 0;
    }
  double x1=x[0];
  double m2,wt;
  breitw_(x1,wsqmin,wsqmax,opts.rmass,opts.rwidth,m2,wt);
  double m=sqrt(m2);
  jac=jac*wt;

  //Dynamic scale (not implemented yet)
  //if(dynamicscale) call scaleset(m2)

  //Limit y boundaries to the kinematic limit in y

  //fodyqt limits
  /*  //.....kinematical limits on qt
  double z=q2/pow(energy_.sroot_,2);
  double xr=pow((1-z),2)-4*z*pow(qt/q,2);
  if (xr < 0)
    {
      f[0] = 0.;
      return 0;
    }
  
  //.....kinematical limits on y
  double tm=sqrt(m2+pow(qt,2));
  double tmpx=(q2+pow(energy_.sroot_,2))/energy_.sroot_/tm;
  ymax=log((tmpx+sqrt(pow(tmpx,2)-4.))/2.);*/

  double ylim = 0.5*log(pow(energy_.sroot_,2)/m2);
  double ymn = min(max(-ylim, ymin),ylim);
  double ymx = max(min(ylim, ymax),-ylim);
  if (ymn >= ymx)
    {
      f[0] = 0.;
      return 0;
    }

  //integrate between ymin and ymax
  double y=ymn+(ymx-ymn)*x[1];
  jac=jac*(ymx-ymn);

  
  //integrate between qtmin and qtmax
  //use qt to get the correct jacobian
  double qtcut = qtcut_.xqtcut_*m;
  double qtmn = max(qtcut, qtmin);
  double cosh2y34=pow((exp(y)+exp(-y))*0.5,2);
  double kinqtlim = sqrt(max(0.,pow(pow(energy_.sroot_,2)+m*m,2)/(4*pow(energy_.sroot_,2)*cosh2y34)-m*m)); //introduced max to avoid neqative argument of sqrt when y=ymax
  double qtmx = min(kinqtlim, qtmax);
  if (qtmn >= qtmx)
    {
      f[0] = 0.;
      return 0;
    }
  double qtmn2 = pow(qtmn,2);
  double qtmx2 = pow(qtmx,2);

  double tiny = 1E-5;
  double a = 1./(1+log(qtmx2/tiny));
  double b = 1./(1+log(qtmn2/tiny));
  double x2 = a + (b-a) * x[2];
  jac = jac * (b-a);
  double qt2=tiny*exp(1./x2 - 1);
  jac=jac*qt2/pow(x2,2);
  double qt=sqrt(qt2);
  jac=jac/qt/2.;

  //  double qt=qtmn+(qtmx-qtmn)*x[2];
  //  jac=jac*(qtmax-qtmn);

  //evaluate the Vj (N)LO cross section
  if (opts.pcubature)
    f[0] = vjint::vint(m,qt,y);
  else
    f[0] = vjfo_(m,qt,y);

  if (f[0] != f[0])
    f[0] = 0.;  //avoid nans
	   
  f[0] = f[0]*jac;

  end_time = clock();
  if (opts.timeprofile)
    cout << setw (3) << "m" << setw(10) << m << setw(4) << "qt" << setw(10) <<  qt << setw(4) << "y" << setw(10) <<  y
	 << setw(8) << "result" << setw(12) << f[0]
	 << setw(10) << "tot time" << setw(10) << float( end_time - begin_time ) /  CLOCKS_PER_SEC
	 << endl;
  
  tell_to_grid_we_are_alive();
  return 0;
}

integrand_t lowintegrand(const int &ndim, const double x[], const int &ncomp, double f[],
                        void* userdata, const int &nvec, const int &core,
                        double &weight, const int &iter)
{
  
  double rlo[22];
  for (int i = 0; i < ndim; i++)
    rlo[i]=x[i];


  dofill_.doFill_ = int(iter==last_iter);
  f[0] = lowint_(rlo,weight,f);
  tell_to_grid_we_are_alive();
  return 0;
}

integrand_t realintegrand(const int &ndim, const double x[], const int &ncomp, double f[],
                        void* userdata, const int &nvec, const int &core,
                        double &weight, const int &iter)
{
  double rre[22];
  for (int i = 0; i < ndim; i++)
    rre[i]=x[i];

  dofill_.doFill_ = int(iter==last_iter);
  f[0] = realint_(rre,weight,f);

  tell_to_grid_we_are_alive();
  return 0;
}

integrand_t virtintegrand(const int &ndim, const double x[], const int &ncomp, double f[],
                        void* userdata, const int &nvec, const int &core,
                        double &weight, const int &iter)
{
  
  double rvi[22];
  for (int i = 0; i < ndim; i++)
    rvi[i]=x[i];
  rvi[9] = rvi[7];
  
  dofill_.doFill_ = int(iter==last_iter);
  f[0] = virtint_(rvi,weight,f);
  
  tell_to_grid_we_are_alive();
  return 0;
}

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
