#include <iostream>
#include <LHAPDF/LHAPDF.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <ctime>

#include "interface.h"
#include "finintegr.h"
#include "integr.h"

using namespace std;

integrand_t ctintegrand(const int &ndim, const double x[], const int &ncomp, double f[])
{
  double wgt = 1;
  //here generate the phase space according to x[], and pass the p vector to countint_
  
  double rct[22];
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

int binner_(double p3[4], double p4[4])
{
  double qt = sqrt((float)pow(p3[0]+p4[0],2)+pow(p3[1]+p4[1],2));
  if (qt < qtmin || qt > qtmax)
    return false;

  double y = 0.5 *log((p3[3] + p4[3] + p3[2] + p4[2]) / (p3[3] + p4[3] - p3[2] - p4[2]));
  if (y < ymin || y > ymax)
    return false;

  //  cout << "qt " << qt << " y " << y << endl;
  return true;
}
