#include "switch.h"
#include "settings.h"

#include <math.h>

//default values
double switching::k = 3./4.;
double switching::delta = 1./2.;
int switching::mode = 1;

const double switching::cutoff = 0.01;

void switching::init()
{
  switching::k = opts.dampk;
  switching::delta = opts.dampdelta;
  switching::mode = opts.dampmode;
}

double switching::swtch(double qt, double m)
{
  double swtch=1.;

  if (opts.fixedorder)
    return swtch;

  if (qt < m*k)
    return swtch;
  
  //Gaussian damping
  if (mode == 1)
    swtch = exp(-pow((k*m-qt),2)/pow((delta*m),2));

  //Exponential damping
  if (mode == 2)
    swtch = exp((pow(k*m,2)-(qt*qt))/pow((delta*m),2));
  
  //Cos damping
  if (mode == 3)
    {
      if (qt > m*k+delta*m)
	return 0;
      swtch = (cos(M_PI/(delta*m)*(qt-k*m))+1.)/2.;
    }
  
  if (swtch <= cutoff)
    swtch = 0;

  return swtch;
}

double switching::qtlimit(double m)
{
  double limit=opts.sroot/2.;

  if (opts.fixedorder)
    return limit;
  
  if (mode == 1)
    limit = m*(k+delta*sqrt(log(1./cutoff)));

  if (mode == 2)
    limit = sqrt(pow(k*m,2) + log(1./cutoff)*pow((delta*m),2));
  
  if (mode == 3)
    limit = k*m+delta*m;

  return limit;
}

double switching_(double &qt, double &m)
{
  switching::swtch(qt, m);
}
