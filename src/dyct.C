#include <iostream>
#include <LHAPDF/LHAPDF.h>
#include <dyct.h>
#include <dyfin.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <ctime>
#include "interface.h"

using namespace std;

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
