#include "vjint.h"
#include "codata.h"
#include "phasespace.h"
#include "settings.h"
#include "coupling.h"
#include "luminosity.h"
#include "mesq.h"
#include "constants.h"

#include <iostream>

double vjint::logx2min;

double vjint::delta(double x)
{
  double fac = gevfb/1000.*coupling::aemmz*asp_.asp_*constants::CF;
  double q2 = phasespace::m2;

  //change of variable to allow integration over 0-1
  //"esp" introduced for better behavior at small qt
  int esp = 8; //16;
  double x2 = exp((1.-pow(x,esp))*logx2min);
  double xjac = x2*logx2min*esp*pow(x,(esp-1));

  //double eps = 1e-7;
  //double x2 = exp(pow(eps,1.-x)*logx2min);
  //double xjac = x2*logx2min*(-log(eps))*pow(eps,1.-x);

  //x2=x2min*exp(1d0/x2min*t)
  //xjac=xjac*x2(1d0-x2min)

  //compute x1 as a function of x2
  double x1 = (opts.sroot*phasespace::mt*phasespace::exppy*x2-q2)/(x2*pow(opts.sroot,2)-opts.sroot*phasespace::mt*phasespace::expmy);

  //check x1,x2<1
  double tiny = 0.;// !1d-8
  if (x1 > 1.-tiny || x2 > 1.-tiny)
    return 0.;

  //Partonic mandelstam invariants (bug fix in DYqT: th <-> uh)
  double sh = x1*x2*pow(opts.sroot,2);
  double th = q2-opts.sroot*x1*phasespace::mt*phasespace::expmy;
  double uh = q2-opts.sroot*x2*phasespace::mt*phasespace::exppy;

  //1/xjj is the jacobian of the integration in dx1 of delta(q2-sh-uh-th) dx1
  //following from the change of variable t = q2-sh-uh-th
  //dx1 = dt * dx1(t)/dt  -> 1/xjj = dx1(t)/dt
  double xjj = abs(x2*pow(opts.sroot,2)-opts.sroot*phasespace::mt*phasespace::expmy);

  //common factor for all the contributions
  double factor = fac/sh;
  
  //compute parton luminosity
  luminosity::pdf1(x1);
  luminosity::pdf2(x2);
  luminosity::calc();
  
  //leading order delta(s2) contributions
  double xloqg=(factor/(pow(coupling::NC,2)-1.))*(aqg0_(sh,th,uh,q2)+agq0_(sh,th,uh,q2));
  double xloqqb=(factor/coupling::NC)*aqqb0_(sh,th,uh,q2);
  double xlo=xloqg+xloqqb;
            
  //next to leading order delta(s2) contributions
  double xnlo = 0.;
  if (opts.order >= 2)
    {
      double s2 = 0.;
      utils_fu_(uh,q2);
      utils_(sh,th,uh,q2,s2);
      utils_dilog_(sh,th,uh,q2);
      double xnloqg = factor/(pow(coupling::NC,2)-1.)
	*(bqg1_(sh,th,uh,q2)+bgq1_(sh,th,uh,q2)+
          bqg2_(sh,th,uh,q2)+bgq2_(sh,th,uh,q2)+
          cqg1_(sh,th,uh,q2)+cgq1_(sh,th,uh,q2)+
          cqg2_(sh,th,uh,q2)+cgq2_(sh,th,uh,q2)+
          bqg3_(sh,th,uh,q2)+bgq3_(sh,th,uh,q2));
      double xnloqqb = factor/coupling::NC
	*(bqqb1_(sh,th,uh,q2)+bqqb2_(sh,th,uh,q2)+
	  cqqb1_(sh,th,uh,q2)+d0aa_(sh,th,uh,q2)+
	  bqqb3_(sh,th,uh,q2));
      
      xnlo = xnloqg+xnloqqb;
    }

  return -M_PI*xjac/xjj * (xlo+(asp_.asp_/2./M_PI)*xnlo);
}
