#include "phasespace.h"
#include "settings.h"

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

void phasespace::gen_mqty(const double x[3], double& jac)
{
  //This routine employs the qtcut, should do also a version without qtcut
  
  //generate phase space as m, qt, y
  //jac gets multiplied by the Jacobian of the change of variables from the unitary cube x[3] to the m, pt, y boundaries

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
  mweight_breitw_(xm,wsqmin,wsqmax,opts.rmass,opts.rwidth,m2,jac);
  m = sqrt(m2);

  //integrate between qtmin and qtmax
  double qtcut = max(opts.qtcut,opts.xqtcut*m);
  double qtmn = max(qtcut, phasespace::qtmin);
  double kinqtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+m2,2)/(4*pow(opts.sroot,2))-m2)); //introduced max to avoid negative argument of sqrt when y=ymax
  double qtmx = min(kinqtlim, phasespace::qtmax);
  if (qtmn >= qtmx)
    {
      jac = 0.;
      return;
    }

  double xqt = x[1];
  qtweight_(xqt,qtmn,qtmx,qt,jac);
  //qtweight_flat_(xqt,qtmn,qtmx,qt,jac)
  qt2 = qt*qt;
  
  //Limit y boundaries to the kinematic limit in y
  calcmt();
  double tmpx=(m2+pow(opts.sroot,2))/opts.sroot/mt;
  double ylim=log((tmpx+sqrt(pow(tmpx,2)-4.))/2.);

  //  double ylim = 0.5*log(pow(opts.sroot,2)/m2);
  double ymn = min(max(-ylim, phasespace::ymin),ylim);
  double ymx = max(min(ylim, phasespace::ymax),-ylim);
  if (ymn >= ymx)
    {
      jac = 0.;
      return;
    }

  //integrate between ymin and ymax
  y=ymn+(ymx-ymn)*x[2];
  jac=jac*(ymx-ymn);

  calcexpy();
}

void phasespace::gen_costhphi(const double x[2], double& jac)
{
  //generate phase space as costh, phi
  //jac gets multiplied by the Jacobian of the change of variables from the unitary square x[2] to the costh, phi_lep boundaries

  //cos of the polar angle of the lepton in the boson rest frame
  costh=-1.+2.*x[0];
  jac=jac*2.;
  
  //lepton azimuthal angle in the boson rest frame
  phi_lep = 2.*M_PI*x[1];
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
      cout << "error in x2min " << x2min
	   << "  " << "m " << m
	   << "pt " << qt
	   << "y " << y << endl;
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


//the m, qt, y, and phi variables are taken from the phasespace:: namelist
void phasespace::genV4p()
{
  //Generate the boson 4-momentum

  //vector boson momentum: pV[3]^2-pV[0]^2-pV[1]^2-pV[2]^2=m2
  double cosphiV, sinphiV;
  if (phiV == 0.)
    {
      cosphiV = 1.;
      sinphiV = 0.;
    }
  else
    {
      cosphiV = cos(phiV);
      sinphiV = sqrt(max(0.,1.-pow(cosphiV,2)))*(phiV>0 ? 1 : -1);
    }
  pV[0]=qt*cosphiV;        //px
  pV[1]=qt*sinphiV;        //py
  pV[2]=0.5*mt*(exppy-expmy);   //pz
  pV[3]=0.5*mt*(exppy+expmy);   //E

  //Calculate the boost 4-vector from the boson rest frame to the laboratory frame (-pV)
  boostv_(m,pV,gam,beta);
}

void phasespace::genl4p()
{
  //simple phase space generation in the naive dilepton rest frame
  //with axes x = (1,0,0); y = (0,1,0); z=(0,0,1)
  //p4 is the lepton and p3 is the antilepton

  //Generate the 4-momentum of the first lepton in the boson rest frame, using costh and phi_lep
  double p4cm[4];
  genp_(costh, phi_lep, m, p4cm);

  //Boost to go in the lab frame
  dyboost_(gam, beta, p4cm, p4);
  
  //momentum of the second lepton
  p3[3]=pV[3]-p4[3];      //E
  p3[0]=pV[0]-p4[0];      //px
  p3[1]=pV[1]-p4[1];      //py
  p3[2]=pV[2]-p4[2];      //pz
}

void phasespace::genp12()
{
  //Generate the incoming partons
  p1[3]=-x1*opts.sroot*0.5;
  p1[0]=0.;
  p1[1]=0.;
  p1[2]=-x1*opts.sroot*0.5;

  p2[3]=-x2*opts.sroot*0.5;
  p2[0]=0.;
  p2[1]=0.;
  p2[2]=+x2*opts.sroot*0.5;
}

void phasespace::genp5()
{
  //Generate jet four momentum
  p5[0]=-p1[0]-p2[0]-p3[0]-p4[0];
  p5[1]=-p1[1]-p2[1]-p3[1]-p4[1];
  p5[2]=-p1[2]-p2[2]-p3[2]-p4[2];
  p5[3]=-p1[3]-p2[3]-p3[3]-p4[3];
}
