#include "abint.h"

#include "gaussrules.h"
#include "settings.h"
#include <math.h>

double *abint::abx;
double *abint::abw;
int abint::abdim;

void abint::init()
{
  //initialize the points of the gaussian quadrature for the alfa and beta integration
  int abintervals = 1; //--> make a setting
  int abrule = 64;     //--> make a setting
  abdim = abintervals*abrule;

  abx = new double [abdim];
  abw = new double [abdim];
    
  double x,t,jac;
  int i,j;

  double min = 1e-7;
  double max = 1.;
  double lmm = log(max/min);
    
  for (int i = 0; i < abintervals; i++)
    {
      double a = min+(max-min)*i/abintervals;
      double b = min+(max-min)*(i+1)/abintervals;
      double c = 0.5*(a+b);
      double m = 0.5*(b-a);
      for (int j = 0; j < abrule; j++)
	{
	  double x = c+m*gr::xxx[abrule-1][j];
	  double t = min*pow(max/min,x);
	  double jac=t*lmm;
	  abx[j+i*abrule]=t;
	  abw[j+i*abrule]=gr::www[abrule-1][j]*m*jac;
	}
    }
}