#include "vjint.h"
#include "settings.h"
#include "interface.h"
#include "resconst.h"
#include "cubature.h"
#include "gaussrules.h"
#include <iomanip>

#include "LHAPDF/LHAPDF.h"

#include <math.h>

const double zmin = 1e-13;
const double zmax = 1.-1e-10;
const double lz = log(zmax/zmin);

const int z1rule = 64;
const int z2rule = 64;
double t1[z1rule];
double t2[z2rule];

double vjint::brz;
double vjint::brw;

void vjint::init()
{
  double az = zmin;
  double bz = zmax;
  double cz = 0.5*(az+bz);
  double mz = 0.5*(bz-az);
  for (int j = 0; j < z1rule; j++)
    t1[j] = zmin*pow(zmax/zmin, cz+mz*gr::xxx[z1rule-1][j]);

  for (int j = 0; j < z2rule; j++)
    t2[j] = zmin*pow(zmax/zmin, cz+mz*gr::xxx[z2rule-1][j]);

  //gevpb_.gevpb_ = 3.8937966e8; //MCFM 6.8 value
  gevpb_.gevpb_ = 0.389379e9; //dyres value

  flagch_.flagch_ = 0;

  vegint_.iter_ = 20;
  vegint_.locall_ = 10000;
  vegint_.nlocall_ = 20000;

  dypara_.ss_ =  opts.sroot;
  dypara_.s_ = pow(dypara_.ss_,2);

  dypdf_.ih1_ = opts.ih1;
  dypdf_.ih2_ = opts.ih2;

  vjorder_.iord_ = opts.order - 1;

  if (opts.nproc == 3)
    {
      if (opts.useGamma && !opts.zerowidth)
	prodflag_.prodflag_ = 5; 
      else
	prodflag_.prodflag_ = 3;
    }
  else if (opts.nproc == 1)
    prodflag_.prodflag_ = 21;
  else if (opts.nproc == 2)
    prodflag_.prodflag_ = 22;

  ckm_.vud_ = cabib_.Vud_;
  ckm_.vus_ = cabib_.Vus_;
  ckm_.vub_ = cabib_.Vub_;
  ckm_.vcd_ = cabib_.Vcd_;
  ckm_.vcs_ = cabib_.Vcs_;
  ckm_.vcb_ = cabib_.Vcb_;
  ckm_.vtd_ = 0.;
  ckm_.vts_ = 0.;
  ckm_.vtb_ = 0.;

  double cw2 = 1.- ewcouple_.xw_;
  em_.aemmz_ = sqrt(2.)* ewcouple_.Gf_ *pow(dymasses_.wmass_,2)* ewcouple_.xw_ /M_PI;

  brz = 1./dymasses_.zwidth_*em_.aemmz_/24.*dymasses_.zmass_ *(pow(-1.+2.*ewcouple_.xw_,2)+pow(2.*ewcouple_.xw_,2))/(ewcouple_.xw_*cw2); //=0.033638
  brw = 1./dymasses_.wwidth_*em_.aemmz_*dymasses_.wmass_/(12.*ewcouple_.xw_); //=0.10906 

  if (opts.zerowidth)
    {
      sigs_.sigz_ = brz; // /(16.*cw2)
      sigs_.sigw_ = brw; // /4.
    }
  
  //quarks are ordered according to mass:
  //1,2,3,4,5,6
  //u,d,s,c,b,t
  
  //double equ = 2./3.;  //up-quarks electric charge
  //double eqd = -1./3.; //down-quarks electric charge
  //quarks_.eq_[0]=equ;                 
  //quarks_.eq_[1]=eqd;
  //quarks_.eq_[2]=eqd;
  //quarks_.eq_[3]=equ;
  //quarks_.eq_[4]=eqd;
  quarks_.eq_[0] = ewcharge_.Q_[MAXNF+2];
  quarks_.eq_[1] = ewcharge_.Q_[MAXNF+1];
  quarks_.eq_[2] = ewcharge_.Q_[MAXNF+3];
  quarks_.eq_[3] = ewcharge_.Q_[MAXNF+4];
  quarks_.eq_[4] = ewcharge_.Q_[MAXNF+5];

  //definition of 'generalized' ckm matrix:
  //    (uu ud us uc ub ut)
  //    (du dd ds dc db dt)
  //ckm=(su sd ss sc sb st)
  //    (cu cd cs cc cb ct)
  //    (bu bd bs bc bb bt)
  //    (tu td ts tc tb tt)
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
      quarks_.ckm_[j][i] = 0.;
  
  quarks_.ckm_[1][0]=ckm_.vud_;
  quarks_.ckm_[2][0]=ckm_.vus_;
  quarks_.ckm_[4][0]=ckm_.vub_;
  quarks_.ckm_[0][1]=ckm_.vud_;
  quarks_.ckm_[3][1]=ckm_.vcd_;
  quarks_.ckm_[5][1]=ckm_.vtd_;
  quarks_.ckm_[0][2]=ckm_.vus_;
  quarks_.ckm_[3][2]=ckm_.vcs_;
  quarks_.ckm_[5][2]=ckm_.vts_;
  quarks_.ckm_[1][3]=ckm_.vcd_;
  quarks_.ckm_[2][3]=ckm_.vcs_;
  quarks_.ckm_[4][3]=ckm_.vcb_;
  quarks_.ckm_[0][4]=ckm_.vub_;
  quarks_.ckm_[3][4]=ckm_.vcb_;
  quarks_.ckm_[5][4]=ckm_.vtb_;
  quarks_.ckm_[1][5]=ckm_.vtd_;
  quarks_.ckm_[2][5]=ckm_.vts_;
  quarks_.ckm_[4][5]=ckm_.vtb_;
      
  //definition of 'delta' matrix
  for (int i = 0; i < MAXNF; i++)
    for (int j = 0; j < MAXNF; j++)
      quarks_.delta_[j][i] = 0.;

  for (int i = 0; i < MAXNF; i++)
    quarks_.delta_[i][i] = 1.;

  //definition of tau3's Pauli matrix
  for (int i = 0; i < MAXNF; i++)
    for (int j = 0; j < MAXNF; j++)
      quarks_.tau3_[j][i] = 0.;

  //quarks_.tau3_[0][0]=1.;
  //quarks_.tau3_[1][1]=-1.;
  //quarks_.tau3_[2][2]=-1.;
  //quarks_.tau3_[3][3]=1.;
  //quarks_.tau3_[4][4]=-1.;
  quarks_.tau3_[0][0] = ewcharge_.tau_[MAXNF+2];
  quarks_.tau3_[1][1] = ewcharge_.tau_[MAXNF+1];
  quarks_.tau3_[2][2] = ewcharge_.tau_[MAXNF+3];
  quarks_.tau3_[3][3] = ewcharge_.tau_[MAXNF+4];
  quarks_.tau3_[4][4] = ewcharge_.tau_[MAXNF+5];

  //definition of constants and couplings
  const2_.ca_ = 3.;
  const2_.xnc_ = 3.;
  const2_.cf_ = 4./3.;
  const2_.tr_ = resconst::NF/2.;
  const2_.pi_ = M_PI;
  dycouplings_.alpha0_ = em_.aemmz_;
  dycouplings_.xw_ = ewcouple_.xw_;
  dycouplings_.sw_ = sqrt(ewcouple_.xw_); // sin_w
  dycouplings_.cw_ = sqrt(1.-ewcouple_.xw_); // cos_w
}


double vjint::vint(double m, double pt, double y)
{
  //set scales and alpha strong
  if (opts.dynamicscale)
    {
      scales2_.xmur_ = m*opts.kmuren;
      scales2_.xmuf_ = m*opts.kmufac;
    }
  else
    {
      scales2_.xmur_ = opts.rmass*opts.kmuren;
      scales2_.xmuf_ = opts.rmass*opts.kmufac;
    }
  
  scales2_.xmur2_ = pow(scales2_.xmur_,2);
  scales2_.xmuf2_ = pow(scales2_.xmuf_,2);
  asnew_.as_ = LHAPDF::alphasPDF(scales2_.xmur_)/M_PI;

  //calculate propagators
  double q2 = m*m;
  double cw2 = 1.- ewcouple_.xw_;
  sigs_.sigz_ = brz*dymasses_.zwidth_*q2/(M_PI*dymasses_.zmass_) / (pow(q2-pow(dymasses_.zmass_,2),2)+pow(dymasses_.zmass_,2)*pow(dymasses_.zwidth_,2)); // /(16.*cw2)
  sigs_.sigw_ = brw*dymasses_.wwidth_*q2/(M_PI*dymasses_.wmass_) / (pow(q2-pow(dymasses_.wmass_,2),2)+pow(dymasses_.wmass_,2)*pow(dymasses_.wwidth_,2)); // !/4.
  sigs_.siggamma_ = em_.aemmz_/(3.*M_PI*q2);
  sigs_.sigint_ = -em_.aemmz_/(6.*M_PI)*(q2-pow(dymasses_.zmass_,2)) /(pow(q2-pow(dymasses_.zmass_,2),2)+pow(dymasses_.zmass_,2)*pow(dymasses_.zwidth_,2)) * (-1.+4.*ewcouple_.xw_)/(2.*sqrt(ewcouple_.xw_*cw2));  // !sqrt(sw2/cw2)/2.

  //set phase space variables
  internal_.q_ = m;
  internal_.q2_ = pow(m,2);
  internal_.qt_ = pt;
  yv_.yv_ = y;
  yv_.expyp_ = exp(y);
  yv_.expym_ = exp(-y);//1./yv_.expyp_;

  if (opts.zerowidth)
    internal_.q_ = opts.rmass;
  
  //call cross section calculation
  //  int ord = opts.order - 1;
  //  double res, err, chi2a;
  //  qtdy_(res,err,chi2a,y,y,ord);
  //  cout << setprecision(16) << "fortran " << m << "  " << pt << "  " << "  " << y << "  " << res << endl;
    
  double res = calc(m, pt, y);
  //cout << setprecision(16) << "C++ " << m << "  " << pt << "  " << "  " << y << "  " << res << endl;
  
  //Apply conversion factors:
  //ds/dqt = ds/dqt2 * dqt2/dqt
  //fb = 1000 * pb
  return 2*pt*res*1000.;
}

int xdelta_cubature(unsigned ndim, const double x[], void *data, unsigned ncomp, double f[])
{
  f[0] = xdelta_(x);
  return 0;
}
int sing_cubature(unsigned ndim, const double x[], void *data, unsigned ncomp, double f[])
{
  double zz[2];
  zz[0] = zmin*pow(zmax/zmin,x[0]);
  double jacz1 = zz[0]*lz;

  zz[1] = zmin*pow(zmax/zmin,x[1]);
  double jacz2 = zz[1]*lz;
  
  f[0] = sing_(zz)*jacz1*jacz2;
  return 0;
}


double vjint::calc(double m, double pt, double y)
{
  double q2 = m*m;
  tm_.tm_ = sqrt(q2+pt*pt);
      
  //kinematical limits on qt (at y = 0)
  double z = q2/dypara_.s_;
  double xr = pow(1-z,2)-4*z*pow(pt/m,2);
  if (xr < 0)
    return 0.;

  //kinematical limits on y including qt
  double tmpx = (q2+dypara_.s_)/dypara_.ss_/tm_.tm_;
  double ymax = log((tmpx+sqrt(max(0.,pow(tmpx,2)-4.)))/2.);

  double ay = fabs(y);
  if (fabs(y) > ymax ||
      fabs(*(short*)& ay - *(short*)& ymax) < 2 ) //check also equality
    return 0.;
  
  //...compute as at order=nloop
  asp_.asp_ = asnew_.as_*M_PI;

  //pcubature integration
  const int ndimx = 1;     //dimensions of the integral
  const int ncomp = 1;  //components of the integrand
  void *userdata = NULL;
  double integral[1];
  double error[1];
  const int eval = 0;
  const double epsrel = min(1e-4,opts.pcubaccuracy);
  const double epsabs = 0.;
  //     boundaries of integration      
  double xmin[1] = {0.};
  double xmax[1] = {1.};

  pcubature(ncomp, xdelta_cubature, userdata, 
	    ndimx, xmin, xmax, 
	    eval, epsabs, epsrel, ERROR_INDIVIDUAL, integral, error);
  double rdelta = integral[0];
  double err = error[0];

  //integration of non-delta(s2) terms (only at NLO)

  //initialise
  double rsing = 0.;

  if (opts.order == 2)
    {
      //pcubature integration
      /*
      const int ndims = 2;     //dimensions of the integral
      //     boundaries of integration      
      double zmn[2] = {zmin,zmin};
      double zmx[2] = {zmax,zmax};
      pcubature(ncomp, sing_cubature, userdata, 
		ndims, zmn, zmx, 
		eval, epsabs, 1e-4, ERROR_INDIVIDUAL, integral, error);
      rsing = integral[0];
      err = sqrt(err*err + error[0]*error[0]);
      */
      //      cout << m << "  " << pt << "  " << y << "  " << integral[0] << "  " << error[0] << endl;
      
      //gaussian quadrature
      double zz[2];
      double az1 = zmin;
      double bz1 = zmax;
      double cz1=0.5*(az1+bz1);
      double mz1=0.5*(bz1-az1);
      for (int j = 0; j < z1rule; j++)
	{
	  double z1 = cz1 + mz1*gr::xxx[z1rule-1][j];
	  double t = t1[j];//zmin*pow(zmax/zmin,z1);
	  double jacz1=t*lz;
	  zz[0]=t;
	  
	  double az2 = zmin;
	  double bz2 = zmax;
	  double cz2 = 0.5*(az2+bz2);
	  double mz2 = 0.5*(bz2-az2);
	  for (int jj = 0; jj < z2rule; jj++)
	    {
	      double z2 = cz2 + mz2*gr::xxx[z2rule-1][jj];
	      double t = t2[jj];//zmin*pow(zmax/zmin,z2);
	      double jacz2 = t*lz;
	      zz[1] = t;
	      rsing=rsing+sing_(zz)*gr::www[z1rule-1][j]*gr::www[z2rule-1][jj]
		*jacz1*jacz2*mz1*mz2;
	    }
	}
      //      cout << m << "  " << pt << "  " << y << "  " << rsing << endl;
      
      rsing = rsing*(asp_.asp_/2./M_PI);
    }
  else if (opts.order == 1)
    rsing = 0.;

  //final result
  double res = rdelta+rsing;
  //cout << m << "  " << pt << "  " << y << "  " << res << "  " << err << endl;
  //      print *,q,qt,yv,rdelta,rsing
  return res;
}
