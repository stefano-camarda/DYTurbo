#include "phasespace.h"
#include "settings.h"
#include "switch.h"
#include "isnan.h"

#include <iostream>
#include <math.h>

//additional useful variables
double phasespace::m2;
double phasespace::qt2;
double phasespace::mt2;
double phasespace::mt;
double phasespace::exppy;
double phasespace::expmy;

//Bjorken-x
double phasespace::x1;
double phasespace::x2;

//Vector boson 4-momentum and boost
double phasespace::pV[4];
double phasespace::gam;
double phasespace::beta[3];

//leptons 4-momenta
double phasespace::p4[4];
double phasespace::p3[4];
double phasespace::p1[4];
double phasespace::p2[4];
double phasespace::p5[4];

//exponent for the change of variable in x2
const double esp=8.;

bool phasespace::gen_m(double x, double& jac, double mlim, bool qtcut, bool qtswitching)
{
  //Generate the boson invariant mass between the integration boundaries,
  //jac gets multiplied by the Jacobian of the change of variable from the unitary interval to the m boundaries
  double qtcutlim = (qtcut && opts.xqtcut != 0.) ? phasespace::qtmax/opts.xqtcut : phasespace::mmax;
  double switchlim = (qtswitching) ? switching::mlimit(phasespace::qtmin) : phasespace::mmin;
  double wsqmin = pow(max(phasespace::mmin,switchlim),2);
  double wsqmax = pow(min(mlim,min(phasespace::mmax,qtcutlim)),2);
  if (wsqmin >= wsqmax)
    return false;
  mweight_breitw_(x,wsqmin,wsqmax,opts.rmass,opts.rwidth,m2,jac);
  m = sqrt(m2);
  return true;
}

bool phasespace::gen_y(double x, double& jac, double ylim)
{
  //Generate the boson rapidity between the integration boundaries,
  //jac gets multiplied by the Jacobian of the change of variable from the unitary interval to the y boundaries
  double ymn = min(max(-ylim, phasespace::ymin),ylim);
  double ymx = max(min(ylim, phasespace::ymax),-ylim);
  if (ymn >= ymx)
    return false;
  y = ymn+(ymx-ymn)*x;
  jac = jac*(ymx-ymn);
  return true;
}

bool phasespace::gen_qt(double x, double& jac, double qtlim, bool qtcut)
{
  //Generate the boson rapidity between the integration boundaries,
  //jac gets multiplied by the Jacobian of the change of variable from the unitary interval to the qt boundaries
  double qtmn = qtcut ?
    max(max(opts.qtcut,opts.xqtcut*m), phasespace::qtmin)
    : max(opts.qtcutoff,phasespace::qtmin);
  double qtmx = min(qtlim, phasespace::qtmax);
  if (qtmn >= qtmx)
    return false;
  if (qtcut) //phase space generation for counterterm and V+j fixed order
    qtweight_(x,qtmn,qtmx,qt,jac);
  else       //phase space generation for resummed cross section
    qtweight_res_(x,qtmn,qtmx,qt,jac);
  //qtweight_flat_(x,qtmn,qtmx,qt,jac);
  qt2 = qt*qt;
  return true;
}

bool phasespace::gen_qt_ctfo(double x, double& jac)
{
  //Generate the boson rapidity between the integration boundaries,
  //jac gets multiplied by the Jacobian of the change of variable from the unitary interval to the qt boundaries
  double qtmn = max(opts.qtcut,opts.xqtcut*m);
  double qtmx = 1e10;
  qtweight_(x,qtmn,qtmx,qt,jac);
  //qtweight_flat_(xqt,qtmn,qtmx,qt,jac);
  qt2 = qt*qt;
  return true;
}

bool phasespace::gen_mqty(const double x[3], double& jac, bool qtcut, bool qtswitching)
{
  //generate phase space as m, qt, y
  //jac gets multiplied by the Jacobian of the change of variables from the unitary cube x[3] to the m, qt, y boundaries
  bool status = true;

  //Generate the boson invariant mass between the integration boundaries
  //kinematic limit from the relation x(1,2) < 1 ==> exp(fabs(y)) < sqrt(s)/m
  double ymn;
  if (ymin*ymax <= 0.)
    ymn = 0.;
  else
    ymn = min(fabs(ymin),fabs(ymax));
  double mlim = opts.sroot/exp(ymn);
  status = gen_m(x[0], jac, mlim, qtcut, qtswitching);
  if (!status)
    return false;

  //Generate the boson transverse momentum between the integration boundaries
  //double qtlim = sqrt(max(0.,pow(pow(opts.sroot,2)-m2,2))/(4*pow(opts.sroot,2))); //introduced max to avoid negative argument of sqrt when  y=ymax
  //Account for ymin-ymax
  double exppylim = exp(ymn);
  double expmylim = 1./exppylim;
  double cosh2y=pow((exppylim+expmylim)*0.5,2);

  //kinematic limit from the relation E < sqrt(s) ==> (mt + qt)*cosh(y) < sqrt(s)
  double qtlim;
  if (qtswitching)
    {
      //in switching mode (ct and res) do not account for y in the qt limit (actually the qtlimit should be 1e10, i.e. unphysical)
      double kinqtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+m2,2)/(4*pow(opts.sroot,2))-m2)); //introduced max to avoid neqative argument of sqrt
      double switchqtlim = switching::qtlimit(m);
      qtlim = min(kinqtlim, switchqtlim);
    }
  else
    //qtlim = sqrt(max(0.,pow(pow(opts.sroot,2)-m2*cosh2y,2)/(4*pow(opts.sroot,2)*cosh2y))); //introduced max to avoid negative argument of sqrt when y=ymax --> wrong formula calculated by me
    qtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+m2,2)/(4*pow(opts.sroot,2)*cosh2y)-m2)); //introduced max to avoid neqative argument of sqrt when y=ymax --> apparently correct formula, cannot check
  
  
  status = gen_qt(x[1], jac, qtlim, qtcut);
  if (!status)
    return false;
  
  //Generate the boson rapidity between the integration boundaries
  calcmt();
  double tmpx=(m2+pow(opts.sroot,2))/opts.sroot/mt;
  double ylim=log((tmpx+sqrt(max(0.,pow(tmpx,2)-4.)))/2.); //Limit y boundaries to the kinematic limit in y --> introduced max to avoid neqative argument of sqrt
  status = gen_y(x[2], jac, ylim);
  if (!status)
    return false;

  return status;
}

bool phasespace::gen_myqt(const double x[3], double& jac, bool qtcut, bool qtswitching)
{
  //With respect to the previous routine, the order in which y and qt are generated is reversed
  //generate phase space as m, y, qt
  //jac gets multiplied by the Jacobian of the change of variables from the unitary cube x[3] to the m, y, qt boundaries
  bool status = true;

  //Generate the boson invariant mass between the integration boundaries
  //kinematic limit from the relation x(1,2) < 1 ==> exp(fabs(y)) < sqrt(s)/m
  double ymn;
  if (ymin*ymax <= 0.)
    ymn = 0.;
  else
    ymn = min(fabs(ymin),fabs(ymax));
  double mlim = opts.sroot/exp(ymn);
  status = gen_m(x[0], jac, mlim, qtcut, qtswitching);
  if (!status)
    return false;

  //Generate the boson rapidity between the integration boundaries
  double ylim;
  if (qtswitching)
    ylim = 0.5*log(pow(opts.sroot,2)/m2); //Limit y boundaries to the kinematic limit in y
  else
    {
      //Limit y boundaries to the kinematic limit in y, accounting for phasespace::qtmin  
      double qtmn = qtcut ?
	max(max(opts.qtcut,opts.xqtcut*m), phasespace::qtmin)
	: phasespace::qtmin;
      double mtmin = sqrt(m2+pow(qtmn,2));
      double tmpx=(m2+pow(opts.sroot,2))/opts.sroot/mtmin;
      ylim=log((tmpx+sqrt(max(0.,pow(tmpx,2)-4.)))/2.); //introduced max to avoid neqative argument of sqrt
    }
  status = gen_y(x[1], jac, ylim);
  if (!status)
    return false;
//  //fodyqt limits
//  //.....kinematical limits on qt
//  double z=q2/pow(energy_.sroot_,2);
//  double xr=pow((1-z),2)-4*z*pow(qt/q,2);
//  if (xr < 0)
//    {
//      f[0] = 0.;
//      return 0;
//    }
//  

  //Generate the boson transverse momentum between the integration boundaries
  double qtlim;
  if (qtswitching)
    {
      //in switching mode (ct and res) do not account for y in the qt limit (actually the qtlimit should be 1e10, i.e. unphysical)
      double kinqtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+m2,2)/(4*pow(opts.sroot,2))-m2)); //introduced max to avoid neqative argument of sqrt
      double switchqtlim = switching::qtlimit(m);
      qtlim = min(kinqtlim, switchqtlim);
    }
  else
    {
      calcexpy();
      double cosh2y=pow((exppy+expmy)*0.5,2);
      qtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+m2,2)/(4*pow(opts.sroot,2)*cosh2y)-m2)); //introduced max to avoid neqative argument of sqrt when y=ymax
    }
  status = gen_qt(x[2], jac, qtlim, qtcut);
  if (!status)
    return false;

  return status;
}

bool phasespace::gen_mqt(const double x[2], double& jac, bool qtcut, bool qtswitching)
{
  //generate phase space as m, qt (for resummed part)
  //jac gets multiplied by the Jacobian of the change of variables from the unitary square x[2] to the m, qt boundaries
  //Current version works for born kinematics, i.e. qt=0 (which is, probably, the correct implementation for born, resummed and counterterm)
  bool status = true;

  //Generate the boson invariant mass between the integration boundaries
  //kinematic limit from the relation x(1,2) < 1 ==> exp(fabs(y)) < sqrt(s)/m
  double ymn;
  if (ymin*ymax <= 0.)
    ymn = 0.;
  else
    ymn = min(fabs(ymin),fabs(ymax));
  double mlim = opts.sroot/exp(ymn);
  status = gen_m(x[0], jac, mlim, qtcut, qtswitching);
  if (!status)
    return false;

  //Generate the boson pt between the integration boundaries
  //Here there is a delicate point, should consider the value of y != 0 when setting the kinematical limit for qt?
  //For counterterm and resummed piece, the cross section is calculated at born level with pt=0, and the qt is generated independently
  //However, this will determine some inconsistencies at large rapidities, because of nonzero cross sections above the kinematical limit for qt
  
  double kinqtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+m2,2)/(4*pow(opts.sroot,2))-m2)); //introduced max to avoid neqative argument of sqrt
  double switchqtlim = switching::qtlimit(m);
  double qtlim = min(kinqtlim, switchqtlim);
  status = gen_qt(x[1], jac, qtlim, qtcut);
  if (!status)
    return false;

  return status;
}  

bool phasespace::gen_my(const double x[2], double& jac, bool qtcut, bool qtswitching)
{
  //generate phase space as m, y
  //jac gets multiplied by the Jacobian of the change of variables from the unitary square x[2] to the m, y boundaries
  //Current version works for born kinematics, i.e. qt=0 (which is, probably, the correct implementation for born, resummed and counterterm)
  bool status = true;

  //Generate the boson invariant mass between the integration boundaries
  //kinematic limit from the relation x(1,2) < 1 ==> exp(fabs(y)) < sqrt(s)/m
  double ymn;
  if (ymin*ymax <= 0.)
    ymn = 0.;
  else
    ymn = min(fabs(ymin),fabs(ymax));
  double mlim = opts.sroot/exp(ymn);
  status = gen_m(x[0], jac, mlim, qtcut, qtswitching);
  if (!status)
    return false;

  //Generate the boson rapidity between the integration boundaries
  //Here there is a delicate point, should consider the value of qt != 0 when setting the kinematical limit for y?
  //For counterterm and resummed piece, the cross section is calculated at born level with pt=0, and the qt is generated independently
  //However, this will determine some inconsistencies at large rapidities, because of nonzero cross sections above the kinematical limit for y
  
  //double mtmin = sqrt(m2+pow(phasespace::qtmin,2));
  //double tmpx=(m2+pow(opts.sroot,2))/opts.sroot/mtmin;
  //double ylim=log((tmpx+sqrt(pow(tmpx,2)-4.))/2.); //Limit y boundaries to the kinematic limit in y

  double ylim = 0.5*log(pow(opts.sroot,2)/m2); //Limit y boundaries to the kinematic limit in y

  status = gen_y(x[1], jac, ylim);
  if (!status)
    return false;

  return status;
}  

void phasespace::gen_costhphi(const double x[2], double& jac)
{
  //generate phase space as costh, phi
  //jac gets multiplied by the Jacobian of the change of variables from the unitary square x[2] to the costh, phi_lep boundaries
  gen_costh(x[0], jac);
  gen_phi(x[1], jac);  
}

void phasespace::gen_costh(const double x, double& jac)
{
  //generate phase space as costh
  //cosine of the polar angle of the lepton in the boson rest frame
  costh=-1.+2.*x;
  jac=jac*2.;
}

void phasespace::gen_phi(const double x, double& jac)
{
  //generate phase space as phi
  //lepton azimuthal angle in the boson rest frame
  phi_lep = -M_PI+2.*M_PI*x; // --> important to generate phi in [-pi,pi] to use fast trigonometric relations as sin(phi) = sqrt(1 - cos2(phi)) * sign(phi)
  jac=jac*2.*M_PI;
}

void phasespace::gen_x2(const double x, double& jac)
{
  //Calculate Bjorken x1 and x2
  //Phase space is generated for integration performed in dx2
  double ss = opts.sroot;
  double s = ss*ss;
  double x2min=(m2-ss*mt*expmy)/(ss*mt*exppy-s);
  double lx2min = log(x2min);
    
  if (x2min > 1. || x2min < 0.)
    {
      cout << "error in x2min " << x2min << " m " << m << " pt " << qt  << " y " << y << endl;
      jac=0.;
      return;
    }
  
  x2=exp((1.-pow(x,esp))*lx2min);
  jac=jac*fabs(x2*lx2min*esp*pow(x,(esp-1.)));
      
  //define x1 as a function of x2 from energy conservation
  x1=(ss*mt*exppy*x2-m2)/(x2*s-ss*mt*expmy);

  if (x1 > 1. || x2 > 1.)
    {
      jac=0.;
      return;
    }

  //Jacobian from the change of variable in delta(p1+p2-p3-p4-p5)
  double xjj=fabs(x2*s-ss*mt*expmy);
  jac = jac/xjj;
}
