#include <iostream>
#include <LHAPDF/LHAPDF.h>
#include <dyct.h>
#include <dyfin.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <ctime>

using namespace std;

double dyct(double m, double y, double qt, double phicm, double phiZ, double cos_th, double alpha, double beta)
{
  //set up constants
  double rmass = 91.1876;
  double rwidth = 2.495;
  double minmass = 40.;
  double sqrts = 7000.;
  double s = sqrts*sqrts;
  double mmin = 66.;
  double mmax = 116.;
  double xqtcut = 0.008;

  //  std::cout << std::setprecision(15);
  clock_t begin_time, end_time;
  double value;
  double rre[22];
  double rvi[22];
  double rct[22];
  double wgt = 1;

  /*
  //variables to be integrated
  //2 dimensions for the counterterm
  double alpha,beta;
  beta = 0.1;
  alpha = 0.1;
  */
  
  //phase space for the counterterm
  double qtcut,xth,tau0,phict;
  qtcut = xqtcut*m;
  xth = 1/(log((qt*qt)/(qtcut*qtcut))+1);
  tau0 = m*m/s;
  phict = -phicm + 0.25;
  if (phict < 0.)
    phict = phict+1.;

  //************** COUNTERTERM ***************
  //phase space mapping
  rct[5] = (m*m - mmin*mmin)/(mmax*mmax-mmin*mmin);     //q2 (why not breit wigner?)
  rct[2] = xth;                                         //costh of the system (qt)
  rct[6] = 0.5 - y/log(tau0);                           //y of the Z boson
  rct[0] = phict;                                       //phi of the system
  rct[3] = (cos_th + 1.)/2.;                            //costh of dilepton in Z rest frame
  rct[4] = phiZ;                                        //phi of dilepton in Z rest frame
  //dummy dimension
  rct[1] = 0.;
  //dimensions to be integrated (alpha and beta of the counterterm)
  rct[7] = beta;
  rct[8] = alpha;
  begin_time = clock();
  value = countint_(rct,wgt);
  end_time = clock();
  //  cout << "Subtraction: " << value << "  " << "time " << float( end_time - begin_time ) /  CLOCKS_PER_SEC << endl;
  //******************************************

  return value;
}

integrand_t ctintegrand(const int &ndim, const double x[], const int &ncomp, double f[])
{
  double wgt = 1;

  double rct[ndim];
  for (int i = 0; i < ndim; i++)
    rct[i]=x[i];

  double value = countint_(rct,wgt);
  
  f[0] = value;
  return 0;
}


//generate the phase space 4 vectors
//write a fortran function which calculates the ct as a function of m, qt, y, and costh moments, (and alpha beta)
//perform the integration in alpha and beta
double dyctint(double m, double y, double qt, double phicm, double phiZ, double costh)
{

  double cthmom0 = 0;
  double cthmom1 = 0;
  double cthmom2 = 0;
  double alpha = 0.5;
  double beta = 0.5;
  countterm_(costh,m,qt,y,alpha,beta,cthmom0,cthmom1,cthmom2);

}
