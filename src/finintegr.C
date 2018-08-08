#include "interface.h"
#include "finintegr.h"
#include "omegaintegr.h"
#include "phasespace.h"
#include "settings.h"
#include "vjint.h"
#include "loint.h"
#include "vjloint.h"
#include "cubacall.h"
#include "isnan.h"

#include <iostream>
#include <iomanip>
#include <ctime>
#include <math.h>


#include "KinematicCuts.h"

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

double vjintegrand_smolyak(int ndim, double x[])
{
  int ncomp = 1;
  double f[ncomp];
  vjintegrand(ndim, x, ncomp, f);
  tell_to_grid_we_are_alive();
  return f[0];
}


integrand_t vjintegrand(const int &ndim, const double x[], const int &ncomp, double f[])
//Generates the phase space 4 vectors
//Calculates the V+j LO/NLO integrand as a function of m, qt, and y
//dOmega integration is already performed in the integrand
//works only without cuts on the leptons
{
  clock_t begin_time, end_time;

  begin_time = clock();

  //Jacobian of the change of variables from the unitary square x[2] to the m, qt boundaries
  double jac = 1.;
  bool status = true;
  
  double r3[3] = {x[0], x[1], x[2]};
  status = phasespace::gen_mqty(r3, jac, true);
  //status = phasespace::gen_myqt(r3, jac, true);
  if (!status)
    {
      f[0] = 0.;
      return 0;
    }
  
  double m = phasespace::m;
  double qt = phasespace::qt;
  double y = phasespace::y;
  phasespace::calcexpy();
  
  //evaluate the Vj (N)LO cross section
  f[0] = vjint::vint(m,qt,y); //C++ interface
  
  if (isnan_ofast(f[0]))
    {
      cout << m << " " << qt << " " << y << " nan in vjint:vint() " << endl;
      f[0] = 0.;  //avoid nans
    }
	   
  if (isnan_ofast(jac))
    {
      cout << m << " " << qt << " " << y << ": nan in gen_myqt() " << endl;
      f[0] = 0.;  //avoid nans
    }
  else
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

double realintegrand_smolyak(int ndim, double x[])
{
  int ncomp = 1;
  double f[ncomp];
  void *userdata;
  const int nvec = 1;
  int core;
  double weight;
  int iter;
  realintegrand(ndim, x, ncomp, f,
		userdata, nvec, core, weight, iter);
  return f[0];
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

integrand_t vjlointegrandMC(const int &ndim, const double x[], const int &ncomp, double f[],
			    void* userdata, const int &nvec, const int &core,
			    double &weight, const int &iter)
{
  f[0] = vjloint::calcvegas(x);

  if (isnan_ofast(f[0]))
    f[0] = 0.;  //avoid nans

  if (iter==last_iter)
    {
      double wt = weight*f[0];
      hists_fill_(phasespace::p3,phasespace::p4,&wt);
    }
  
  tell_to_grid_we_are_alive();
  return 0;
}


integrand_t vjlointegrand(const int &ndim, const double x[], const int &ncomp, double f[])
{
  clock_t begin_time, end_time;

  begin_time = clock();

  double ff[2];
  vjloint::calc(x, ff);

  f[0] = ff[0];
  if (opts.helicity >= 0)
    f[1] = ff[1];

  if (isnan_ofast(f[0]))
    {
      f[0] = 0.;  //avoid nans
      if (opts.helicity >= 0)
	f[1]=0.;
    }
  
  end_time = clock();
  if (opts.timeprofile)
    cout << setw (3) << "m" << setw(10) << phasespace::m << setw(4) << "qt" << setw(10) <<  phasespace::qt << setw(4) << "y" << setw(10) <<  phasespace::y
	 << setw (6) << "costh" << setw(10) << phasespace::costh << setw(4) << "phi" << setw(10) <<  phasespace::phi_lep << setw(3) << "x2" << setw(10) <<  phasespace::x2
	 << setw(8) << "result" << setw(12) << f[0]
	 << setw(10) << "tot time" << setw(10) << float( end_time - begin_time ) /  CLOCKS_PER_SEC
	 << endl;
  
  tell_to_grid_we_are_alive();
  return 0;
}

int vjlointegrand_cubature(unsigned ndim, const double x[], void *data, unsigned ncomp, double f[])
{
  vjlointegrand(ndim, x, ncomp, f);
  tell_to_grid_we_are_alive();
  return 0;
}

int vjlointegrand_cubature_v(unsigned ndim, long unsigned npts, const double x[], void *data, unsigned ncomp, double f[])
{
#pragma omp parallel for num_threads(opts.cubacores) copyin(scale_,facscale_,qcdcouple_,vjloint::mur,vjloint::muf)
  for (unsigned i = 0; i < npts; i++)
    {
      // evaluate the integrand for npts points
      double xi[ndim];
      double fi[ncomp];
      for (unsigned j = 0; j < ndim; j++)
	xi[j] = x[i*ndim + j];

      vjlointegrand(ndim, xi, ncomp, fi);
      
      for (unsigned k = 0; k < ncomp; ++k)
	f[i*ncomp + k] = fi[k];
    }
  tell_to_grid_we_are_alive();
  return 0;
}

double vjlointegrand_smolyak(int ndim, double x[])
{
  int ncomp = 1;
  double f[ncomp];
  vjlointegrand(ndim, x, ncomp, f);
  tell_to_grid_we_are_alive();
  return f[0];
}
