#include "interface.h"
#include "finintegr.h"
#include "omegaintegr.h"
#include "phasespace.h"
#include "settings.h"
#include "vjint.h"
#include "loint.h"

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
  double mcut = (opts.xqtcut != 0.) ? phasespace::qtmax/opts.xqtcut : phasespace::mmax;
  double wsqmin = pow(phasespace::mmin,2);
  double wsqmax = pow(min(phasespace::mmax,mcut),2);
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
  double ymn = min(max(-ylim, phasespace::ymin),ylim);
  double ymx = max(min(ylim, phasespace::ymax),-ylim);
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
  double qtcut = max(opts.qtcut,opts.xqtcut*m);
  double qtmn = max(qtcut, phasespace::qtmin);
  double cosh2y34=pow((exp(y)+exp(-y))*0.5,2);
  double kinqtlim = sqrt(max(0.,pow(pow(energy_.sroot_,2)+m*m,2)/(4*pow(energy_.sroot_,2)*cosh2y34)-m*m)); //introduced max to avoid neqative argument of sqrt when y=ymax
  double qtmx = min(kinqtlim, phasespace::qtmax);
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
  //  jac=jac*(phasespace::qtmax-qtmn);

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

integrand_t v2jintegrand(const int &ndim, const double x[], const int &ncomp, double f[],
			 void* userdata, const int &nvec, const int &core,
			 double &weight, const int &iter)
{
  /*
  clock_t begin_time, end_time;

  begin_time = clock();

  //Jacobian of the change of variables from the unitary hypercube x[3] to the m, y, qt boundaries
  double jac = 1.;

  // Generate the boson invariant mass between the integration boundaries
  double wsqmin = pow(phasespace::mmin,2);
  double wsqmax = pow(phasespace::mmax,2);
  double x1=x[0];
  double m2,wt;
  breitw_(x1,wsqmin,wsqmax,opts.rmass,opts.rwidth,m2,wt);
  double m=sqrt(m2);
  jac=jac*wt;

  //Limit y boundaries to the kinematic limit in y
  double ylim = 0.5*log(pow(energy_.sroot_,2)/m2);
  double ymn = min(max(-ylim, phasespace::ymin),ylim);
  double ymx = max(min(ylim, phasespace::ymax),-ylim);
  if (ymn >= ymx)
    {
      f[0]=0.;
      return 0;
    }

  //integrate between ymin and ymax
  double y=ymn+(ymx-ymn)*x[1];
  jac=jac*(ymx-ymn);

  //integrate between qtmin and qtmax
  double qtmn = phasespace::qtmin;
  double cosh2y34=pow((exp(y)+exp(-y))*0.5,2);
  double qtlim = sqrt(pow(pow(energy_.sroot_,2)+m*m,2)/(4*pow(energy_.sroot_,2)*cosh2y34)-m*m);
  double qtmx = min(qtlim, phasespace::qtmax);
  if (qtmn >= qtmx)
    {
      f[0]=0.;
      return 0;
    }
  double qt=qtmn+(qtmx-qtmn)*x[2];
  jac=jac*(qtmx-qtmn);

  //angular variables
  double phicm = x[3]*2.*M_PI;
  double cos_th = x[4]*2.-1.;
  double phiZ = x[5]*2.*M_PI;

  //variables to be integrated (4 dimensions in total)
  //1 collinear PDF dimension
  double zcth;
  zcth = x[6]*2.-1.;   //polar angle of Z in the CM frame
  //3 dimensions of real radiation
  double mjjmin = 0.;
  double mjjmax = opts.sroot;
  double mjj = mjjmin + x[7]*(mjjmax - mjjmin);    //invariant mass of the dijet (set to 0 to have the correct virtual phase space mapping)
  jac=jac*(mjjmax-mjjmin);
  double phijj = x[8];  //phi of the dijet in dijet rest frame
  double costhjj = x[9];//costh of the dijet in dijet rest frame

  //Here should create the l1, l2, j1, j2, 4-momenta from the above parametrisation
  //Instead, for now,
  //convert to 10-dim hypercube used in gen4()
  //need to calculate tau and yboost
  
  //evaluate four momentum of Z in lab frame, given y, m, qt
  double zql, zqt, zp, ze;
  zql = 0.5*sqrt(m*m+qt*qt)*(exp(y)-exp(-y));
  zqt = qt;
  zp = sqrt(zqt*zqt + zql*zql);
  ze = sqrt(m*m+zqt*zqt+zql*zql);

  //Compute Z four momentum in CM, given y, m, qt and zcth
  double zqtcm, zpcm, zqlcm, zecm, zycm;
  zqtcm = qt;
  zpcm = zqtcm / sin(acos(zcth)); //absolute momentum of Z in CM
  zqlcm = zpcm * zcth;         //longitudinal momentum of Z in CM
  zecm = sqrt(m*m + zpcm*zpcm);//Energy of Z in CM
  zycm = 0.5 * log((zecm + zqlcm)/(zecm - zqlcm)); //rapidity of Z in CM

  //boost from CM to lab
  double yboost;
  yboost = y - zycm;

  //Compute recoil four momentum in CM (massless recoil with mjj = 0)
  double rqtcm, rqlcm, recm, rycm;
  rqtcm = zqtcm;
  rqlcm = zqlcm;
  recm = sqrt(rqtcm*rqtcm+rqlcm*rqlcm+mjj*mjj);
  rycm = 0.5*log((recm-rqlcm)/(recm+rqlcm));

  //now apply yboost to the recoil to have the longitudinal momentum of the recoil in the lab frame
  double rqt, ry, rql, re;
  rqt = rqtcm;
  ry = rycm+yboost;
  rql = 0.5*sqrt(rqtcm*rqtcm+mjj*mjj)*(exp(ry)-exp(-ry));
  re = sqrt(rql*rql+rqt*rqt+mjj*mjj);
  
  //compute tau = mtot^2/s
  double s = pow(opts.sroot,2);
  double mtot, tau;
  mtot = sqrt(pow(ze+re,2)-pow(zql+rql,2));
  tau = mtot*mtot/s;

  if (tau > 1.)
    {
      f[0] = 0.;
      return 0;
    }
    

  //cross check that rapidity of the system is the rapidity of the boost
  //  double ytot;
  //  ytot = 0.5* log((ze+re + zql+rql)/(ze+re - (zql+rql)));
  //  cout << ytot << "  " << yboost << endl;

  //convert invariant mass to 0-1 for Breit-Wigner weighting
  double rmass = opts.rmass;
  double rwidth = opts.rwidth;
  double M2, M1, s3max, almin, almax, tanal, al, X1;
  M2 = mjj;
  M1 = mtot;
  s3max=(M2-M1)*(M2-M1);
  almin=atan((0.-rmass*rmass)/rmass/rwidth);
  almax=atan((s3max-rmass*rmass)/rmass/rwidth);
  tanal = (m*m - rmass*rmass)/rmass/rwidth;
  al = atan(tanal);
  X1 = (al-almin)/(almax-almin);

  //************** REAL RADIATION ***************
  //phase space mapping
  double minmass = sqrt(taumin_.taumin_)*energy_.sroot_;
  double rre[22];
  rre[1] = X1;                                     //mll
  rre[3] = x[3];//phicm/2./M_PI;                                  //phi of Z in CM frame
  rre[6] = x[4];//(cos_th + 1.)/2.;                       //cos_th of dilepton in Z rest frame
  rre[7] = x[5];//phiZ/2./M_PI;                                   //phi of dilepton in Z rest frame
  rre[8] = log(tau)/log(minmass*minmass/s);        //tau of the system
  rre[9] = 0.5 - yboost/log(tau);                  //y of the system
  //dimension to be integrated
  rre[2] = (zcth + 1.)/2.;                         //cos th of Z in CM frame
  //dimensions of real radiation to be integrated
  rre[0] = mjj*mjj / (mtot*mtot);                  //mjj
  rre[4] = phijj;                                  //phi jj
  rre[5] = costhjj;                                //costh jj
  */

  double rre[22];
  for (int i = 0; i < ndim; i++)
    rre[i]=x[i];

  dofill_.doFill_ = int(iter==last_iter);
  f[0] = v2jint_(rre,weight,f);

  tell_to_grid_we_are_alive();
  return 0;
}

integrand_t lointegrand(const int &ndim, const double x[], const int &ncomp, double f[],
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

  // Generate the boson invariant mass between the integration boundaries
  double wsqmin = pow(phasespace::mmin,2);
  double wsqmax = pow(phasespace::mmax,2);
  double x1=x[0];
  double m2,wt;
  breitw_(x1,wsqmin,wsqmax,opts.rmass,opts.rwidth,m2,wt);
  double m=sqrt(m2);
  jac=jac*wt;

  //Limit y boundaries to the kinematic limit in y
  double ylim = 0.5*log(pow(energy_.sroot_,2)/m2);
  double ymn = min(max(-ylim, phasespace::ymin),ylim);
  double ymx = max(min(ylim, phasespace::ymax),-ylim);
  if (ymn >= ymx)
    {
      f[0]=0.;
      return 0;
    }

  //integrate between ymin and ymax
  double y=ymn+(ymx-ymn)*x[1];
  jac=jac*(ymx-ymn);

  //costh
  double costh = x[2]*2.-1.;
  jac=jac*2.;
  
  dofill_.doFill_ = int(iter==last_iter);
  f[0] = loint::lint(costh,m,y,0)*jac;

  //phi integration is not needed because an overall 2*pi factor is already in the integrand
  //f[0] = f[0] * 2.*M_PI;
  
  tell_to_grid_we_are_alive();
  return 0;
}

int lointegrand2d_cubature_v(unsigned ndim, long unsigned npts, const double x[], void *data, unsigned ncomp, double f[])
{
#pragma omp parallel for num_threads(opts.cubacores)
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

  // Generate the boson invariant mass between the integration boundaries
  double wsqmin = pow(phasespace::mmin,2);
  double wsqmax = pow(phasespace::mmax,2);
  double x1=x[0];
  double m2,wt;
  breitw_(x1,wsqmin,wsqmax,opts.rmass,opts.rwidth,m2,wt);
  double m=sqrt(m2);
  jac=jac*wt;

  //Limit y boundaries to the kinematic limit in y
  double ylim = 0.5*log(pow(energy_.sroot_,2)/m2);
  double ymn = min(max(-ylim, phasespace::ymin),ylim);
  double ymx = max(min(ylim, phasespace::ymax),-ylim);
  if (ymn >= ymx)
    {
      f[0]=0.;
      return 0;
    }

  //integrate between ymin and ymax
  double y=ymn+(ymx-ymn)*x[1];
  jac=jac*(ymx-ymn);

  //dofill_.doFill_ = int(iter==last_iter);
  f[0] = loint::lint(0.,m,y,1)*jac;

  //phi integration is not needed because an overall 2*pi factor is already in the integrand
  //f[0] = f[0] * 2.*M_PI;
  
  tell_to_grid_we_are_alive();
  return 0;
}
