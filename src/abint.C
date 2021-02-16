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
  abdim = opts.abintervals*opts.abrule;

  abx = new double [abdim];
  abw = new double [abdim];
    
  double min = 1e-7;
  double max = 1.;
  double lmm = log(max/min);
    
  for (int i = 0; i < opts.abintervals; i++)
    {
      double a = 0.+(1.-0.)*i/opts.abintervals;     
      double b = 0.+(1.-0.)*(i+1)/opts.abintervals;
      double c = 0.5*(a+b);
      double m = 0.5*(b-a);
      for (int j = 0; j < opts.abrule; j++)
	{
	  //with change of variable t = eps^(1-x), and with lower cutoff eps on z1 z2
	  double x = c+m*gr::xxx[opts.abrule-1][j];
	  double t = min*pow(max/min,x);
	  double jac=t*lmm;
	  abx[j+i*opts.abrule]=t;
	  abw[j+i*opts.abrule]=gr::www[opts.abrule-1][j]*m*jac;

	  //without change of variable, and without lower cutoff on z1 z2
	  //abx[j+i*opts.abrule]=c+m*gr::xxx[opts.abrule-1][j];
	  //abw[j+i*opts.abrule]=gr::www[opts.abrule-1][j]*m;
	}
    }
}
void abint::release()
{
  delete[] abx;
  delete[] abw;

}
