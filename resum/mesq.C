#include "mesq.h"
#include "settings.h"
#include "interface.h"
#include "coupling.h"
#include "integr.h"
#include "rapint.h"
#include <iostream>
#include <cmath>

//constants
const double eequ =2./3.;
const double eeqd = -1./3.;
const double gevfb = 3.8937966e11;

double mesq::fac;
//Z couplings
double mesq::mZ2;
double mesq::wZ2;
double mesq::gLZu;
double mesq::gLZd;
double mesq::gRZu;
double mesq::gRZd;
double mesq::fLZ;
double mesq::fRZ;
double mesq::fLpfR;
double mesq::fLmfR;
double mesq::ugLpgR;
double mesq::ugLmgR;
double mesq::dgLpgR;
double mesq::dgLmgR;
//W coupling
double mesq::mW2;
double mesq::wW2;
double mesq::gLWfLW;
//gamma* coupling
double mesq::aem2pi;
double mesq::aem2pi2;

//mass dependent variables
double mesq::q2;
double mesq::propZ;
double mesq::propW;
double mesq::propG;
double mesq::propZG;

complex <double> mesq::mesqij[12];
complex <double>* mesq::mesqij_expy;

//fortran interface
void setmesq_expy_(int& mode, double& m, double& costh, double& y)
{
  mesq::setmesq_expy(mode, m, costh, y);
};

void mesq::init()
{
  mZ2 = pow(coupling::zmass, 2);
  wZ2 = pow(dymasses_.zwidth_,2);
  mW2 = pow(coupling::wmass, 2);
  wW2 = pow(dymasses_.wwidth_,2);

  double gZ=sqrt(sqrt(2.)*coupling::Gf*mZ2);
  double gW=sqrt(4.*sqrt(2.)*coupling::Gf*mW2);

  gLZu=gZ*(1./2.-eequ*coupling::xw);
  gLZd=gZ*(-1./2.-eeqd*coupling::xw);
  gRZu=-gZ*eequ*coupling::xw;
  gRZd=-gZ*eeqd*coupling::xw;
  fLZ=gZ*(-1./2.+coupling::xw);
  fRZ=gZ*coupling::xw;

  double gLW=gW/sqrt(2.);
  double fLW=gW/sqrt(2.);

  //Initialize constants
  fac=1./9./M_PI*gevfb;

  //Z couplings
  fLpfR = pow(fLZ,2)+pow(fRZ,2);
  fLmfR = pow(fLZ,2)-pow(fRZ,2);
  ugLpgR = pow(gLZu,2)+pow(gRZu,2);
  ugLmgR = pow(gLZu,2)-pow(gRZu,2);
  dgLpgR = pow(gLZd,2)+pow(gRZd,2);
  dgLmgR = pow(gLZd,2)-pow(gRZd,2);

  //W coupling
  gLWfLW = pow(gLW,2)*pow(fLW,2)/16.;

  //gamma* coupling
  aem2pi = 2.*M_PI*coupling::aemmz;
  aem2pi2 = pow(aem2pi,2);

  //allocate memory
  mesqij_expy = new complex <double> [mellinint::mdim*mellinint::mdim*2*12];

  q2 = 0;
  propZ = 0;
  propW = 0;
  propG = 0;
  propZG = 0;
}

void mesq::setmesq_expy(int mode, double m, double costh, double y)
{
  //Number of partonic channels
  int totpch;
  if (opts.nproc == 3)
    totpch = 10; //only 4 partonic channels are actually needed
  else
    totpch = 12;

  //mass dependent part
  q2 = pow(m,2);
  if (opts.nproc == 3)
    propZ = q2/(pow(q2-mZ2,2)+mZ2*wZ2);
  else
    propW = q2/(pow(q2-mW2,2)+mW2*wW2);
  if (opts.useGamma)
    {
      propG = 1./q2;
      propZG = (q2-mZ2)/(pow(q2-mZ2,2)+mZ2*wZ2);
    }

  if (mode < 2)
    {
      //costh dependent part
      double one, costh1, costh2;
      if (mode == 0)
	{
	  one = 1.;
	  costh1 = costh;
	  costh2 = pow(costh,2);
	}
      else if (mode == 1)
	{
	  //retrieve cos theta moments
	  double cthmom0, cthmom1, cthmom2;
	  cthmoments_(cthmom0, cthmom1, cthmom2);
	  one = cthmom0;
	  costh1 = cthmom1;
	  costh2 = cthmom2;
	}

      setmesq(one, costh1, costh2);
      //mass and rapidity dependent part
      double bjx= q2/pow(opts.sroot,2);
      double ax = log(bjx);
      double ax1 = (ax+2*y)/2.; 
      double ax2 = (ax-2*y)/2.; 
      complex <double> cex1[mellinint::mdim];
      complex <double> cex2p[mellinint::mdim];
      complex <double> cex2m[mellinint::mdim];
      for (int i = 0; i < mellinint::mdim; i++)
	{
	  cex1[i] = exp(-mellinint::Np[i] * ax1) / M_PI * mellinint::CCp;
	  cex2p[i] = exp(-mellinint::Np[i] * ax2) / M_PI * mellinint::CCp;
	  cex2m[i] = exp(-mellinint::Nm[i] * ax2) / M_PI * mellinint::CCm; //Is this the complex conjugate of the above?
	}

      //product of amplitudes and x1^-z1 * x2^-z2 piece of the Mellin inverse transform
      for (int i1 = 0; i1 < mellinint::mdim; i1++)
	for (int i2 = 0; i2 < mellinint::mdim; i2++)
	  for (int pch = 0; pch < totpch; pch++)
	    {
	      mesqij_expy[mesq::index(pch, i1, i2, mesq::positive)] = mesqij[pch] * cex1[i1] * cex2p[i2] * mellinint::wn[i1] * mellinint::wn[i2];
	      mesqij_expy[mesq::index(pch, i1, i2, mesq::negative)] = mesqij[pch] * cex1[i1] * cex2m[i2] * mellinint::wn[i1] * mellinint::wn[i2];
	    }
    }
  else //if mode == 2 -> rapidity integrated mode
    {
      for (int i1 = 0; i1 < mellinint::mdim; i1++)
	for (int i2 = 0; i2 < mellinint::mdim; i2++)
	  {
	    complex <double> one, costh1, costh2;

	    //retrieve Ithxp Ithxm positive branch integrals of the x1^-z1 * x2^-z2 piece 
	    one = rapint::Ith0p[mellinint::index(i1,i2)];
	    costh1 = rapint::Ith1p[mellinint::index(i1,i2)];
	    costh2 = rapint::Ith2p[mellinint::index(i1,i2)];
	    setmesq(one, costh1, costh2);
	    for (int pch = 0; pch < totpch; pch++)
	      mesqij_expy[mesq::index(pch, i1, i2, mesq::positive)] = mesqij[pch];

	    //retrieve Ithxp Ithxm positive branch integrals of the x1^-z1 * x2^-z2 piece 
	    one = rapint::Ith0m[mellinint::index(i1,i2)];
	    costh1 = rapint::Ith1m[mellinint::index(i1,i2)];
	    costh2 = rapint::Ith2m[mellinint::index(i1,i2)];
	    setmesq(one, costh1, costh2);
	    for (int pch = 0; pch < totpch; pch++)
	      mesqij_expy[mesq::index(pch, i1, i2, mesq::negative)] = mesqij[pch];
	  }
    }
}

template <class T>
void mesq::setmesq(T one, T costh1, T costh2)
{
  T omcosth2 = one - 2.*costh1 + costh2;
  T opcosth2 = one + 2.*costh1 + costh2;
  if (opts.nproc == 1) //W+ 
    {
      mesqij[0]=fac*propW*gLWfLW*pow(cabib_.Vud_,2)*omcosth2;
      mesqij[1]=fac*propW*gLWfLW*pow(cabib_.Vud_,2)*opcosth2;
      mesqij[2]=fac*propW*gLWfLW*pow(cabib_.Vus_,2)*omcosth2;
      mesqij[3]=fac*propW*gLWfLW*pow(cabib_.Vus_,2)*opcosth2;
      mesqij[4]=fac*propW*gLWfLW*pow(cabib_.Vub_,2)*omcosth2;
      mesqij[5]=fac*propW*gLWfLW*pow(cabib_.Vub_,2)*opcosth2;
      mesqij[6]=fac*propW*gLWfLW*pow(cabib_.Vcs_,2)*omcosth2;
      mesqij[7]=fac*propW*gLWfLW*pow(cabib_.Vcs_,2)*opcosth2;
      mesqij[8]=fac*propW*gLWfLW*pow(cabib_.Vcd_,2)*omcosth2;
      mesqij[9]=fac*propW*gLWfLW*pow(cabib_.Vcd_,2)*opcosth2;
      mesqij[10]=fac*propW*gLWfLW*pow(cabib_.Vcb_,2)*omcosth2;
      mesqij[11]=fac*propW*gLWfLW*pow(cabib_.Vcb_,2)*opcosth2;
    }
  if (opts.nproc == 2) //W-
    {
      mesqij[0]=fac*propW*gLWfLW*pow(cabib_.Vud_,2)*omcosth2;
      mesqij[1]=fac*propW*gLWfLW*pow(cabib_.Vud_,2)*opcosth2;
      mesqij[2]=fac*propW*gLWfLW*pow(cabib_.Vus_,2)*omcosth2;
      mesqij[3]=fac*propW*gLWfLW*pow(cabib_.Vus_,2)*opcosth2;
      mesqij[4]=fac*propW*gLWfLW*pow(cabib_.Vub_,2)*omcosth2;
      mesqij[5]=fac*propW*gLWfLW*pow(cabib_.Vub_,2)*opcosth2;
      mesqij[6]=fac*propW*gLWfLW*pow(cabib_.Vcs_,2)*omcosth2;
      mesqij[7]=fac*propW*gLWfLW*pow(cabib_.Vcs_,2)*opcosth2;
      mesqij[8]=fac*propW*gLWfLW*pow(cabib_.Vcd_,2)*omcosth2;
      mesqij[9]=fac*propW*gLWfLW*pow(cabib_.Vcd_,2)*opcosth2;
      mesqij[10]=fac*propW*gLWfLW*pow(cabib_.Vcb_,2)*omcosth2;
      mesqij[11]=fac*propW*gLWfLW*pow(cabib_.Vcb_,2)*opcosth2;
    }
  if (opts.nproc == 3 && !opts.useGamma) //Z 
    {
      mesqij[0]=fac*propZ*(ugLpgR*fLpfR*(one + costh2)-ugLmgR*fLmfR*(2.*costh1)); //u-ubar
      mesqij[6]=mesqij[0]; //c-cbar
	    
      mesqij[1]=fac*propZ*(ugLpgR*fLpfR*(one + costh2)+ugLmgR*fLmfR*(2.*costh1)); //ubar-u
      mesqij[7]=mesqij[1]; //cbar-c

      mesqij[2]=fac*propZ*(dgLpgR*fLpfR*(one + costh2)-dgLmgR*fLmfR*(2.*costh1)); //d-dbar
      mesqij[4]=mesqij[2]; //s-sbar
      mesqij[8]=mesqij[2]; //b-bbar

      mesqij[3]=fac*propZ*(dgLpgR*fLpfR*(one + costh2)+dgLmgR*fLmfR*(2.*costh1)); //dbar-d
      mesqij[5]=mesqij[3]; //sbar-s
      mesqij[9]=mesqij[3]; //bbar-b
    }
  if (opts.nproc == 3 && opts.useGamma) // Z/gamma*
    {
      //                           Z                 gamma*                  interference
      mesqij[0]=fac*((propZ*ugLpgR*fLpfR + propG*aem2pi2*pow(eequ,2) - propZG*aem2pi*eequ*(gLZu+gRZu)*(fLZ+fRZ)) * (one + costh2)
		     -(propZ*ugLmgR*fLmfR                            - propZG*aem2pi*eequ*(gLZu-gRZu)*(fLZ-fRZ)) * (2.*costh1));
      mesqij[6]=mesqij[0];

      mesqij[1]=fac*((propZ*ugLpgR*fLpfR + propG*aem2pi2*pow(eequ,2) - propZG*aem2pi*eequ*(gLZu+gRZu)*(fLZ+fRZ)) * (one + costh2)
		     +(propZ*ugLmgR*fLmfR                            - propZG*aem2pi*eequ*(gLZu-gRZu)*(fLZ-fRZ)) * (2.*costh1));
      mesqij[7]=mesqij[1];

      mesqij[2]=fac*((propZ*dgLpgR*fLpfR + propG*aem2pi2*pow(eeqd,2) - propZG*aem2pi*eeqd*(gLZd+gRZd)*(fLZ+fRZ)) * (one + costh2)
		     -(propZ*dgLmgR*fLmfR                            - propZG*aem2pi*eeqd*(gLZd-gRZd)*(fLZ-fRZ)) * (2.*costh1));
      mesqij[4]=mesqij[2];
      mesqij[8]=mesqij[2];

      mesqij[3]=fac*((propZ*dgLpgR*fLpfR + propG*aem2pi2*pow(eeqd,2) - propZG*aem2pi*eeqd*(gLZd+gRZd)*(fLZ+fRZ)) * (one + costh2)
		     +(propZ*dgLmgR*fLmfR                            - propZG*aem2pi*eeqd*(gLZd-gRZd)*(fLZ-fRZ)) * (2.*costh1));
      mesqij[5]=mesqij[3];
      mesqij[9]=mesqij[3];
    }
}
