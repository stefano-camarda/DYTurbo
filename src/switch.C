#include "switch.h"
#include "settings.h"
#include <iostream>
#include <iomanip>

#include <math.h>

//default values
double switching::k = 3./4.;
double switching::delta = 1./2.;
int switching::mode = 1;

const double switching::cutoff = 0.01;

//this tolerance on the cutoff is introduced to avoid 0 values at the upper pt limit of integration,
//it is needed because pcubature evaluate the integrals at the borders of the phase space
const double switching::tolerance = 1.-1e-10; //--> probably not needed anymore

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

  if (!opts.damp)
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
  
  if (swtch < cutoff*tolerance)  //--> probably not needed anymore, since the phase space is generated up to the qt and m limits
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

double switching::mlimit(double qt)
{
  //This is a lower limit
  double limit=0.;

  if (opts.fixedorder)
    return limit;
  
  if (mode == 1)
    limit = qt/(k+delta*sqrt(log(1./cutoff)));

  if (mode == 2)
    limit = sqrt(pow(qt,2)/(pow(k,2)+log(1./cutoff)*pow(delta,2)));
  
  if (mode == 3)
    limit = qt/(k+delta);

  return limit;
}

double switching_(double &qt, double &m)
{
  return switching::swtch(qt, m);
}
