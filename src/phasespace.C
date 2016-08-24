#include "phasespace.h"

#include <iostream>
#include <math.h>

//integration boundaries
double phasespace::mmin;
double phasespace::mmax;
double phasespace::qtmin;
double phasespace::qtmax;
double phasespace::ymin;
double phasespace::ymax;

double phasespace::cthmin = -1;
double phasespace::cthmax = +1;

//global variables for the generation of the phase space
double phasespace::m;
double phasespace::qt;
double phasespace::y;
double phasespace::phiV;

//square variables
double phasespace::m2;
double phasespace::qt2;
double phasespace::mt2;

double phasespace::costh;

using namespace std;

void phasespace::setbounds(double m1, double m2, double qt1, double qt2, double y1, double y2)
{
  mmin = m1;
  mmax = m2;
  qtmin = qt1;
  qtmax = qt2;
  ymin = y1;
  ymax = y2;
    
  //Check ordering
  if (mmin > mmax || qtmin > qtmax || ymin > ymax)
    {
      cout << "Error on integration boundaries" << endl;
      exit(-1);
    }
}

void phasespace::setcthbounds(double cth1, double cth2)
{
  cthmin = cth1;
  cthmax = cth2;
    
  //Check ordering
  if (cthmin > cthmax)
    {
      cout << "Error on costh integration boundaries" << endl;
      exit(-1);
    }
}

void phasespace::set_mqtyphi(double M, double Qt, double Y, double PhiV)
{
  m = M;
  qt = Qt;
  y = Y;
  phiV = PhiV;
}
void phasespace::set_m(double M) {m = M;}
void phasespace::set_qt(double Qt) {qt = Qt;}
void phasespace::set_y(double Y) {y = Y;}
void phasespace::set_phiV(double PhiV) {phiV = PhiV;}

void phasespace::set_cth(double Costh) {costh = Costh;}

//fortran functions
void setqt_(double &qtt) {phasespace::qt = qtt;}
void sety_(double &yy) {phasespace::y = yy;}

void phasespace::gen_mqty(const double x[3], double& jac)
{
  //generate phase space as m, qt, y
  //jac is the Jacobian of the change of variables from the unitary hypercube x[3] to the m, pt, y boundaries
  jac = 1.;

  // Generate the boson invariant mass between the integration boundaries
  double mcut = (opts.xqtcut != 0.) ? phasespace::qtmax/opts.xqtcut : phasespace::mmax;
  double wsqmin = pow(phasespace::mmin,2);
  double wsqmax = pow(min(phasespace::mmax,mcut),2);
  if (wsqmin >= wsqmax)
    {
      jac = 0.;
      return;
    }
  double xm=x[0];
  double m2,wtm;
  breitw_(xm,wsqmin,wsqmax,opts.rmass,opts.rwidth,m2,wtm);
  m=sqrt(m2);
  jac=jac*wtm;

  //integrate between qtmin and qtmax
  double qtcut = max(opts.qtcut,opts.xqtcut*m);
  double qtmn = max(qtcut, phasespace::qtmin);
  double kinqtlim = sqrt(max(0.,pow(pow(energy_.sroot_,2)+m2,2)/(4*pow(energy_.sroot_,2))-m2)); //introduced max to avoid neqative argument of sqrt when y=ymax
  double qtmx = min(kinqtlim, phasespace::qtmax);
  if (qtmn >= qtmx)
    {
      jac = 0.;
      return;
    }

  double xqt = x[2];
  qtweight_(xqt,qtmn,qtmx,qt,jac);
  //qtweight_flat_(xqt,qtmn,qtmx,qt,jac)
  qt2 = qt*qt;
  
  //Limit y boundaries to the kinematic limit in y
  mt2 = m2+qt2;
  double mt = sqrt(mt2);
  double tmpx=(m2+pow(energy_.sroot_,2))/energy_.sroot_/mt;
  double ylim=log((tmpx+sqrt(pow(tmpx,2)-4.))/2.);

  //  double ylim = 0.5*log(pow(energy_.sroot_,2)/m2);
  double ymn = min(max(-ylim, phasespace::ymin),ylim);
  double ymx = max(min(ylim, phasespace::ymax),-ylim);
  if (ymn >= ymx)
    {
      jac = 0.;
      return;
    }

  //integrate between ymin and ymax
  double y=ymn+(ymx-ymn)*x[1];
  jac=jac*(ymx-ymn);
}


//This is a boson variables binner: add mass here (not yet needed, because the mass is always in the phase space boundaries, will be needed when doing also mass bins)
//the boson binner function is used only with the MCFM integrands, it is needed because the generation of the phase space is not done with respect to m, pt, y
int binner_(double p3[4], double p4[4])
{
  double qt = sqrt((float)pow(p3[0]+p4[0],2)+pow(p3[1]+p4[1],2));
  if (qt < phasespace::qtmin || qt > phasespace::qtmax)
    return false;

  double y = 0.5 *log((p3[3] + p4[3] + p3[2] + p4[2]) / (p3[3] + p4[3] - p3[2] - p4[2]));
  if (y < phasespace::ymin || y > phasespace::ymax)
    return false;

  double m = sqrt(pow(p3[3]+p4[3],2) - pow(p3[2]+p4[2],2) - qt*qt);
  if (m < phasespace::mmin || m > phasespace::mmax)
    return false;

  //move the costh cut to the lepton variables binner
  if (phasespace::cthmin != -1 && phasespace::cthmax != 1)
    {
      double costh = ((p3[3]+p3[2])*(p4[3]-p4[2])-(p4[3]+p4[2])*(p3[3]-p3[2]))/sqrt(m*m*(m*m+qt*qt));
      costh *= (y < 0. ? -1 : 1); //sign flip according to boson rapidity
      if (costh < phasespace::cthmin || costh > phasespace::cthmax)
	return false;
    }
  
  //cout << "qt " << qt << " y " << y << endl;
  return true;
}

//Make also a lepton variables binner, which is needed for producing distributions as a function of lepton variables with quadrature integration
