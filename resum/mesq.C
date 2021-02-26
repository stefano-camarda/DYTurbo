#include "mesq.h"
#include "settings.h"
#include "dyres_interface.h"
#include "coupling.h"
#include "omegaintegr.h"
#include "rapint.h"
#include "parton.h"
#include "codata.h"
#include "clenshawcurtisrules.h"
#include "gaussrules.h"
#include "propagator.h"
#include "pdf.h"
#include <iostream>
#include <cmath>

using namespace parton;

//constants
const double eequ =2./3.;
const double eeqd = -1./3.;
//const double gevfb = 3.8937966e11;

double mesq::fac;
//Z couplings
double mesq::mZ2;
double mesq::wZ2;
double mesq::gLZu;
double mesq::gLZd;
double mesq::gRZu;
double mesq::gRZd;

double mesq::gLZ[MAXNF];
double mesq::gRZ[MAXNF];

double mesq::fLZ;
double mesq::fRZ;
double mesq::fLpfR;
double mesq::fLmfR;
double mesq::ugLpgR;
double mesq::ugLmgR;
double mesq::dgLpgR;
double mesq::dgLmgR;

double mesq::gLpgR[MAXNF];
double mesq::gLmgR[MAXNF];

//W coupling
double mesq::mW2;
double mesq::wW2;
double mesq::gLWfLW;

//gamma* coupling
double mesq::aem2pi;
double mesq::aem2pi2;
double mesq::Q[MAXNF];

//box at 3-loops
double mesq::gLpgR_sum[MAXNF];
double mesq::gLmgR_sum[MAXNF];
double mesq::Q2_sum[MAXNF];
double mesq::QgZp_sum[MAXNF];
double mesq::QgZm_sum[MAXNF];

//mass dependent variables
double mesq::q2;
double mesq::propZ;
double mesq::propW;
double mesq::propG;
double mesq::propZG;

complex <double> mesq::mesqij[12];
complex <double>* mesq::mesqij_expy;

complex <double> mesq::mesqij_NFZ[10];
complex <double>* mesq::mesqij_expy_NFZ;

//parton mapping for subprocesses
int mesq::totpch;

//parton id in x space (PDG convention)
pdgid *mesq::pid1;
pdgid *mesq::pid2;
pdgid p1Z[10] = {U,Ub,D,Db,S,Sb,C,Cb,B,Bb};
pdgid p2Z[10] = {Ub,U,Db,D,Sb,S,Cb,C,Bb,B};
pdgid p1Wp[12] = {U,Db,U,Sb,U,Bb,C,Sb,C,Db,C,Bb};
pdgid p2Wp[12] = {Db,U,Sb,U,Bb,U,Sb,C,Db,C,Bb,C};
pdgid p1Wm[12] = {D,Ub,S,Ub,B,Ub,S,Cb,D,Cb,B,Cb};
pdgid p2Wm[12] = {Ub,D,Ub,S,Ub,B,Cb,S,Cb,D,Cb,B};

//parton id in N space (DYRES convention)
partid *mesq::pidn1;
partid *mesq::pidn2;
partid p1Zn[10] = {u,ub,d,db,s,sb,c,cb,b,bb};
partid p2Zn[10] = {ub,u,db,d,sb,s,cb,c,bb,b};
partid p1Wpn[12] = {u,db,u,sb,u,bb,c,sb,c,db,c,bb};
partid p2Wpn[12] = {db,u,sb,u,bb,u,sb,c,db,c,bb,c};
partid p1Wmn[12] = {d,ub,s,ub,b,ub,s,cb,d,cb,b,cb};
partid p2Wmn[12] = {ub,d,ub,s,ub,b,cb,s,cb,d,cb,b};

//fortran interface
void setmesq_expy_(int& mode, double& m, double& costh, double& y)
{
  mesq::setmesq_expy(mode, m, costh, y);
};

void mesq::init()
{
  //Number of partonic channels
  if (opts.nproc == 3)
    totpch = 10; //only 4 partonic channels are actually needed
  else
    totpch = 12;

  if (opts.nproc == 3)
    {
      pid1 = p1Z;
      pid2 = p2Z;
      pidn1 = p1Zn;
      pidn2 = p2Zn;
    }
  else if (opts.nproc == 1)
    {
      pid1 = p1Wp;
      pid2 = p2Wp;
      pidn1 = p1Wpn;
      pidn2 = p2Wpn;
    }
  else if (opts.nproc == 2)
    {
      pid1 = p1Wm;
      pid2 = p2Wm;
      pidn1 = p1Wmn;
      pidn2 = p2Wmn;
    }

  mZ2 = pow(coupling::zmass, 2);
  mW2 = pow(coupling::wmass, 2);
  wZ2 = pow(coupling::zwidth,2);
  wW2 = pow(coupling::wwidth,2);

  double gZ=sqrt(sqrt(2.)*coupling::Gf*mZ2);
  double gW=sqrt(4.*sqrt(2.)*coupling::Gf*mW2);

  gLZu=gZ*(1./2.-eequ*coupling::xw);
  gLZd=gZ*(-1./2.-eeqd*coupling::xw);
  gRZu=-gZ*eequ*coupling::xw;
  gRZd=-gZ*eeqd*coupling::xw;

  //put here ka(q)
  for (int j = 0; j < MAXNF; j++)
    {
      gLZ[j]=gZ*(ewcharge_.tau_[j+MAXNF+1]*1./2.-ewcharge_.Q_[j+MAXNF+1]*coupling::xw);
      gRZ[j]=-gZ*ewcharge_.Q_[j+MAXNF+1]*coupling::xw;
    }

  //put here ka(e)
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

  for (int j = 0; j < MAXNF; j++)
    {
      gLpgR[j] = pow(gLZ[j],2)+pow(gRZ[j],2);
      gLmgR[j] = pow(gLZ[j],2)-pow(gRZ[j],2);
    }

  //W coupling
  gLWfLW = pow(gLW,2)*pow(fLW,2)/16.;

  //gamma* coupling
  aem2pi = 2.*M_PI*coupling::aemmz;
  aem2pi2 = pow(aem2pi,2);
  for (int j = 0; j < MAXNF; j++)
    Q[j] = ewcharge_.Q_[j+MAXNF+1];

  //Sum couplings for the NFZ component (virtual box at 3-loops)
  //Pure Z
  for (int j = 0; j < MAXNF; j++)
    {
      gLpgR_sum[j] = 0.;
      for (int f = 0; f < MAXNF; f++)
	gLpgR_sum[j] += gLZ[j]*gLZ[f]+gRZ[j]*gRZ[f];
      gLmgR_sum[j] = 0.;
      for (int f = 0; f < MAXNF; f++)
	gLmgR_sum[j] += gLZ[j]*gLZ[f]-gRZ[j]*gRZ[f];
    }
  //gamma*
  for (int j = 0; j < MAXNF; j++)
    {
      Q2_sum[j] = 0.;
      for (int f = 0; f < MAXNF; f++)
	Q2_sum[j] += Q[j]*Q[f];
    }
  //interference (mixed terms in Q and gZ)
  for (int j = 0; j < MAXNF; j++)
    {
      QgZp_sum[j] = 0.;
      for (int f = 0; f < MAXNF; f++)
	QgZp_sum[j] += 1./2.*(Q[j]*(gLZ[f]+gRZ[f])+Q[f]*(gLZ[j]+gRZ[j]));
      QgZm_sum[j] = 0.;
      for (int f = 0; f < MAXNF; f++)
	QgZm_sum[j] += 1./2.*(Q[j]*(gLZ[f]-gRZ[f])+Q[f]*(gLZ[j]-gRZ[j]));
    }
  
  q2 = 0;
  propZ = 0;
  propW = 0;
  propG = 0;
  propZG = 0;
}

void mesq::setpropagators(double m)
{
  q2 = pow(m,2);

  prop::set(q2);
  //prop::set_real(q2);
  //prop::set_offset(q2);
  propW = prop::W;
  propZ = prop::Z;
  propG = prop::G;
  propZG = prop::ZG;
  
  /*
  //fixed width propagators
  if (opts.nproc == 3)
    propZ = q2/(pow(q2-mZ2,2)+mZ2*wZ2);
  else
    propW = q2/(pow(q2-mW2,2)+mW2*wW2);
  if (opts.useGamma)
    {
      propG = 1./q2;
      propZG = (q2-mZ2)/(pow(q2-mZ2,2)+mZ2*wZ2);
    }
  */

  //flat propagators, for tests
  //propZ = 1./m*(8./3.)*pow(opts.sroot,2)/2/(1./9./M_PI*gevfb);
  //propW = 1./m*(8./3.)*pow(opts.sroot,2)/2/(1./9./M_PI*gevfb);
  //propG = 1./m*(8./3.)*pow(opts.sroot,2)/2/(1./9./M_PI*gevfb);
}

void mesq::allocate()
{
  //allocate memory
  mesqij_expy     = new complex <double> [mellinint::mdim*mellinint::mdim*2*12];
  if (opts.order >= 3)
    mesqij_expy_NFZ = new complex <double> [mellinint::mdim*mellinint::mdim*2*10];
}

void mesq::setmesq_expy(int mode, double m, double costh, double y, int helicity)
{
  //mode 0: differential
  //mode 1: integrated in costh
  //mode 2: integrated in costh and y
  //mode 3: integrated in costh, y, and pt
  //mode 4: integrated in costh and pt

  //mass dependent part
  setpropagators(m);
  
  if (mode < 2 || mode == 4)
    {
      //costh dependent part
      double one, costh1, costh2;
      if (mode == 0)
	{
	  one = 1.;
	  costh1 = costh;
	  costh2 = pow(costh,2);
	}
      else if (mode == 1 || mode == 4)
	{
	  //retrieve cos theta moments
	  double cthmom0, cthmom1, cthmom2;
	  cthmoments_(cthmom0, cthmom1, cthmom2);
	  //omegaintegr::cthmoments(cthmom0, cthmom1, cthmom2);
	  if (helicity >= 0)
	    {
	      if (helicity == 0)
		{
		  cthmom0 = 4./3.; //20./3. * (0.5*2. - 1.5*2./3.  ) +2./3.*2.;
		  cthmom1 = 0.;
		  cthmom2 = -4./3.; //20./3. * (0.5*2./3.-1.5*2./5.  ) +2./3.*2./3.;
		}
	      else if (helicity == 4)
		{
		  cthmom0=0.;
		  //cthmom1=4*cthmom2*(y < 0 ? -1. : 1.);
		  cthmom1=4*cthmom2;
		  cthmom2=0.;
		}
	      else
		{
		  cthmom0=0.;
		  cthmom1=0.;
		  cthmom2=0.;
		}
	    }
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
//	  cex1[i] = exp(-mellinint::Np[i] * ax1) / M_PI * mellinint::CCp;
//	  cex2p[i] = exp(-mellinint::Np[i] * ax2) / M_PI * mellinint::CCp;
//	  cex2m[i] = exp(-mellinint::Nm[i] * ax2) / M_PI * mellinint::CCm; //Is this the complex conjugate of the above?
	  cex1[i] = exp(-mellinint::Np_1[i] * ax1) / M_PI * mellinint::CCp;
	  cex2p[i] = exp(-mellinint::Np_2[i] * ax2) / M_PI * mellinint::CCp;
	  cex2m[i] = exp(-mellinint::Nm_2[i] * ax2) / M_PI * mellinint::CCm; //Is this the complex conjugate of the above?
	}

      //product of amplitudes and x1^-z1 * x2^-z2 piece of the Mellin inverse transform
      for (int i1 = 0; i1 < mellinint::mdim; i1++)
	for (int i2 = 0; i2 < mellinint::mdim; i2++)
	  for (int pch = 0; pch < totpch; pch++)
	    {
	      //mesqij_expy[mesq::index(pch, i1, i2, mesq::positive)] = mesqij[pch] * cex1[i1] * cex2p[i2] * mellinint::wn[i1] * mellinint::wn[i2];
	      //mesqij_expy[mesq::index(pch, i1, i2, mesq::negative)] = mesqij[pch] * cex1[i1] * cex2m[i2] * mellinint::wn[i1] * mellinint::wn[i2];
	      mesqij_expy[mesq::index(pch, i1, i2, mesq::positive)] = mesqij[pch] * cex1[i1] * cex2p[i2] * mellinint::wn_1[i1] * mellinint::wn_2[i2];
	      mesqij_expy[mesq::index(pch, i1, i2, mesq::negative)] = mesqij[pch] * cex1[i1] * cex2m[i2] * mellinint::wn_1[i1] * mellinint::wn_2[i2];
	      //cout << i1 << "  " << i2 << "  " << mesqij[pch] << "  " <<  cex1[i1] << "  " <<  cex2m[i2] << "  " <<  mellinint::wn[i1] << "  " <<  mellinint::wn[i2] << endl;
	      if (opts.nproc == 3 && opts.order >= 3)
		{
		  mesqij_expy_NFZ[mesq::index(pch, i1, i2, mesq::positive)] = mesqij_NFZ[pch] * cex1[i1] * cex2p[i2] * mellinint::wn_1[i1] * mellinint::wn_2[i2];
		  mesqij_expy_NFZ[mesq::index(pch, i1, i2, mesq::negative)] = mesqij_NFZ[pch] * cex1[i1] * cex2m[i2] * mellinint::wn_1[i1] * mellinint::wn_2[i2];
		}
	    }
    }
  else //if mode == 2 or mode == 3 -> rapidity integrated mode
    {
      //1D mellin
      if (opts.mellin1d) //helicity option is not valid in mellin 1d if collisions are pp, because cannot apply the y <> 0 switch!
	{
	  double cthmom0, cthmom1, cthmom2;
	  omegaintegr::cthmoments(cthmom0,cthmom1,cthmom2);
	  //cthmom0 = 2.;
	  //cthmom1 = 0.;
	  //cthmom2 = 2./3.;
	  if (helicity >= 0)
	    {
	      if (helicity == 0)
		{
		  cthmom0 = 4./3.; //20./3. * (0.5*2. - 1.5*2./3.  ) +2./3.*2.;
		  cthmom1 = 0.;
		  cthmom2 = -4./3.; //20./3. * (0.5*2./3.-1.5*2./5.  ) +2./3.*2./3.;
		}
	      else if (helicity == 4)
		{
		  cthmom0=0.;
		  //cthmom1=4*cthmom2*(y < 0 ? -1. : 1.);
		  cthmom1=4*cthmom2;
		  cthmom2=0.;
		}
	      else
		{
		  cthmom0=0.;
		  cthmom1=0.;
		  cthmom2=0.;
		}
	    }
	  setmesq(cthmom0, cthmom1, cthmom2);
	  double bjx= q2/pow(opts.sroot,2);
	  double ax = log(bjx);

	  for (int i = 0; i < mellinint::mdim; i++)
	    for (int pch = 0; pch < totpch; pch++)
	      {
		complex <double> cexp = exp(-mellinint::Np[i] * ax)/M_PI * mellinint::CCp/complex <double>(0.,1);
		complex <double> cexm = exp(-mellinint::Nm[i] * ax)/M_PI * mellinint::CCm/complex <double>(0.,1);
		mesqij_expy[mesq::index(pch, i, i, mesq::positive)] = mesqij[pch] * cexp * mellinint::wn[i];
		mesqij_expy[mesq::index(pch, i, i, mesq::negative)] = mesqij[pch] * cexm * conj(mellinint::wn[i]);
		if (opts.nproc == 3 && opts.order >= 3)
		  {
		    mesqij_expy_NFZ[mesq::index(pch, i, i, mesq::positive)] = mesqij_NFZ[pch] * cexp * mellinint::wn[i];
		    mesqij_expy_NFZ[mesq::index(pch, i, i, mesq::negative)] = mesqij_NFZ[pch] * cexm * conj(mellinint::wn[i]);
		  }
	      }

	  /*
	  //clenshaw-curtis with weight functions
	  double m = 0.5;
	  double jac = opts.zmax;
	  cc::setw(ax*jac*m);
	  for (int i = 0; i < mellinint::mdim; i++)
	    for (int pch = 0; pch < totpch; pch++)
	      {
		complex <double> cexp = exp(-ax*(opts.cpoint+1.+complex <double>(0.,1)*jac*m))/M_PI * mellinint::CCp/complex <double>(0.,1);
		complex <double> cexm = exp(-ax*(opts.cpoint+1.+complex <double>(0.,1)*jac*m))/M_PI * mellinint::CCm/complex <double>(0.,1);
		mesqij_expy[mesq::index(pch, i, i, mesq::positive)] = mesqij[pch] * cexp * (cc::cosw[opts.mellinrule-1][i] - complex<double>(0.,1)*cc::sinw[opts.mellinrule-1][i])*m*jac;
		mesqij_expy[mesq::index(pch, i, i, mesq::negative)] = mesqij[pch] * cexm * (cc::cosw[opts.mellinrule-1][i] + complex<double>(0.,1)*cc::sinw[opts.mellinrule-1][i])*m*jac;
	      }
	  */

	  
	}
      else
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
		  {
		    mesqij_expy[mesq::index(pch, i1, i2, mesq::positive)] = mesqij[pch];
		    if (opts.nproc == 3 && opts.order >= 3)
		      mesqij_expy_NFZ[mesq::index(pch, i1, i2, mesq::positive)] = mesqij_NFZ[pch];
		  }
		//retrieve Ithxp Ithxm negative branch integrals of the x1^-z1 * x2^-z2 piece 
		one = rapint::Ith0m[mellinint::index(i1,i2)];
		costh1 = rapint::Ith1m[mellinint::index(i1,i2)];
		costh2 = rapint::Ith2m[mellinint::index(i1,i2)];
		setmesq(one, costh1, costh2);
		for (int pch = 0; pch < totpch; pch++)
		  {
		    mesqij_expy[mesq::index(pch, i1, i2, mesq::negative)] = mesqij[pch];
		    if (opts.nproc == 3 && opts.order >= 3)
		      mesqij_expy_NFZ[mesq::index(pch, i1, i2, mesq::negative)] = mesqij_NFZ[pch];
		  }
	      }
	}
    }
}
void mesq::free()
{
  delete[] mesqij_expy;
  if (opts.order >= 3)
    delete[] mesqij_expy_NFZ;
}

template <class T>
void mesq::setmesq(T one, T costh1, T costh2)
{
  //Important! in this matrix elements definition costh is the angle between the antilepton
  //and the parton from the incoming beam 1, that is opposite sign with respect to the MCFM convention (angle between lepton and parton 1).
  //To restore the usual convention, need to flip the sign of costh.
  //This would be better done in the formulas below.
  costh1 = -costh1;
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
      /*
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
      */
      //0 1 2 3 4
      //d u s c b
      mesqij[0]=fac*propZ*(gLpgR[1]*fLpfR*(one + costh2)-gLmgR[1]*fLmfR*(2.*costh1)); //u-ubar //+(ka(eq)-ka(e)*ka(q))*sin2thetaw*
      mesqij[6]=fac*propZ*(gLpgR[3]*fLpfR*(one + costh2)-gLmgR[3]*fLmfR*(2.*costh1)); //c-cbar
	    
      mesqij[1]=fac*propZ*(gLpgR[1]*fLpfR*(one + costh2)+gLmgR[1]*fLmfR*(2.*costh1)); //ubar-u
      mesqij[7]=fac*propZ*(gLpgR[3]*fLpfR*(one + costh2)+gLmgR[3]*fLmfR*(2.*costh1)); //cbar-c

      mesqij[2]=fac*propZ*(gLpgR[0]*fLpfR*(one + costh2)-gLmgR[0]*fLmfR*(2.*costh1)); //d-dbar
      mesqij[4]=fac*propZ*(gLpgR[2]*fLpfR*(one + costh2)-gLmgR[2]*fLmfR*(2.*costh1)); //s-sbar
      mesqij[8]=fac*propZ*(gLpgR[4]*fLpfR*(one + costh2)-gLmgR[4]*fLmfR*(2.*costh1)); //b-bbar

      mesqij[3]=fac*propZ*(gLpgR[0]*fLpfR*(one + costh2)+gLmgR[0]*fLmfR*(2.*costh1)); //dbar-d
      mesqij[5]=fac*propZ*(gLpgR[2]*fLpfR*(one + costh2)+gLmgR[2]*fLmfR*(2.*costh1)); //sbar-s
      mesqij[9]=fac*propZ*(gLpgR[4]*fLpfR*(one + costh2)+gLmgR[4]*fLmfR*(2.*costh1)); //bbar-b
    }
  if (opts.nproc == 3 && opts.useGamma) // Z/gamma*
    {
      //                           Z                 gamma*                  interference
      mesqij[0]=fac*((propZ*gLpgR[1]*fLpfR + propG*aem2pi2*pow(Q[1],2) - propZG*aem2pi*Q[1]*(gLZ[1]+gRZ[1])*(fLZ+fRZ)) * (one + costh2)
		     -(propZ*gLmgR[1]*fLmfR                            - propZG*aem2pi*Q[1]*(gLZ[1]-gRZ[1])*(fLZ-fRZ)) * (2.*costh1));
      mesqij[6]=fac*((propZ*gLpgR[3]*fLpfR + propG*aem2pi2*pow(Q[3],2) - propZG*aem2pi*Q[3]*(gLZ[3]+gRZ[3])*(fLZ+fRZ)) * (one + costh2)
		     -(propZ*gLmgR[3]*fLmfR                            - propZG*aem2pi*Q[3]*(gLZ[3]-gRZ[3])*(fLZ-fRZ)) * (2.*costh1));

      mesqij[1]=fac*((propZ*gLpgR[1]*fLpfR + propG*aem2pi2*pow(Q[1],2) - propZG*aem2pi*Q[1]*(gLZ[1]+gRZ[1])*(fLZ+fRZ)) * (one + costh2)
		     +(propZ*gLmgR[1]*fLmfR                            - propZG*aem2pi*Q[1]*(gLZ[1]-gRZ[1])*(fLZ-fRZ)) * (2.*costh1));
      mesqij[7]=fac*((propZ*gLpgR[3]*fLpfR + propG*aem2pi2*pow(Q[3],2) - propZG*aem2pi*Q[3]*(gLZ[3]+gRZ[3])*(fLZ+fRZ)) * (one + costh2)
		     +(propZ*gLmgR[3]*fLmfR                            - propZG*aem2pi*Q[3]*(gLZ[3]-gRZ[3])*(fLZ-fRZ)) * (2.*costh1));

      mesqij[2]=fac*((propZ*gLpgR[0]*fLpfR + propG*aem2pi2*pow(Q[0],2) - propZG*aem2pi*Q[0]*(gLZ[0]+gRZ[0])*(fLZ+fRZ)) * (one + costh2)
		     -(propZ*gLmgR[0]*fLmfR                            - propZG*aem2pi*Q[0]*(gLZ[0]-gRZ[0])*(fLZ-fRZ)) * (2.*costh1));
      mesqij[4]=fac*((propZ*gLpgR[2]*fLpfR + propG*aem2pi2*pow(Q[2],2) - propZG*aem2pi*Q[2]*(gLZ[2]+gRZ[2])*(fLZ+fRZ)) * (one + costh2)
		     -(propZ*gLmgR[2]*fLmfR                            - propZG*aem2pi*Q[2]*(gLZ[2]-gRZ[2])*(fLZ-fRZ)) * (2.*costh1));
      mesqij[8]=fac*((propZ*gLpgR[4]*fLpfR + propG*aem2pi2*pow(Q[4],2) - propZG*aem2pi*Q[4]*(gLZ[4]+gRZ[4])*(fLZ+fRZ)) * (one + costh2)
		     -(propZ*gLmgR[4]*fLmfR                            - propZG*aem2pi*Q[4]*(gLZ[4]-gRZ[4])*(fLZ-fRZ)) * (2.*costh1));

      //!!! bug for flavour decomposition of the Z/gamma* interference !!!
      mesqij[3]=fac*((propZ*gLpgR[0]*fLpfR + propG*aem2pi2*pow(Q[0],2) - propZG*aem2pi*Q[0]*(gLZd+gRZd)*(fLZ+fRZ)) * (one + costh2)
		     +(propZ*gLmgR[0]*fLmfR                            - propZG*aem2pi*Q[0]*(gLZd-gRZd)*(fLZ-fRZ)) * (2.*costh1));
      mesqij[5]=fac*((propZ*gLpgR[2]*fLpfR + propG*aem2pi2*pow(Q[2],2) - propZG*aem2pi*Q[2]*(gLZd+gRZd)*(fLZ+fRZ)) * (one + costh2)
		     +(propZ*gLmgR[2]*fLmfR                            - propZG*aem2pi*Q[2]*(gLZd-gRZd)*(fLZ-fRZ)) * (2.*costh1));
      mesqij[9]=fac*((propZ*gLpgR[4]*fLpfR + propG*aem2pi2*pow(Q[4],2) - propZG*aem2pi*Q[4]*(gLZd+gRZd)*(fLZ+fRZ)) * (one + costh2)
		     +(propZ*gLmgR[4]*fLmfR                            - propZG*aem2pi*Q[4]*(gLZd-gRZd)*(fLZ-fRZ)) * (2.*costh1));

      //this should be correct --> to be tested
      /*
      mesqij[3]=fac*((propZ*gLpgR[0]*fLpfR + propG*aem2pi2*pow(Q[0],2) - propZG*aem2pi*Q[0]*(gLZ[0]+gRZ[0])*(fLZ+fRZ)) * (one + costh2)
		     +(propZ*gLmgR[0]*fLmfR                            - propZG*aem2pi*Q[0]*(gLZ[0]-gRZ[0])*(fLZ-fRZ)) * (2.*costh1));
      mesqij[5]=fac*((propZ*gLpgR[2]*fLpfR + propG*aem2pi2*pow(Q[2],2) - propZG*aem2pi*Q[2]*(gLZ[2]+gRZ[2])*(fLZ+fRZ)) * (one + costh2)
		     +(propZ*gLmgR[2]*fLmfR                            - propZG*aem2pi*Q[2]*(gLZ[2]-gRZ[2])*(fLZ-fRZ)) * (2.*costh1));
      mesqij[9]=fac*((propZ*gLpgR[4]*fLpfR + propG*aem2pi2*pow(Q[4],2) - propZG*aem2pi*Q[4]*(gLZ[4]+gRZ[4])*(fLZ+fRZ)) * (one + costh2)
		     +(propZ*gLmgR[4]*fLmfR                            - propZG*aem2pi*Q[4]*(gLZ[4]-gRZ[4])*(fLZ-fRZ)) * (2.*costh1));
      */
      
      /*
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
      */
    }

  //NFZ terms: box contributions at 3-loops
  if (opts.nproc == 3 && !opts.useGamma && opts.order >= 3)
    {
      mesqij_NFZ[0]=fac*propZ*(gLpgR_sum[1]*fLpfR*(one + costh2)-gLmgR_sum[1]*fLmfR*(2.*costh1)); //u-ubar
      mesqij_NFZ[6]=fac*propZ*(gLpgR_sum[3]*fLpfR*(one + costh2)-gLmgR_sum[3]*fLmfR*(2.*costh1)); //c-cbar
	    
      mesqij_NFZ[1]=fac*propZ*(gLpgR_sum[1]*fLpfR*(one + costh2)+gLmgR_sum[1]*fLmfR*(2.*costh1)); //ubar-u
      mesqij_NFZ[7]=fac*propZ*(gLpgR_sum[3]*fLpfR*(one + costh2)+gLmgR_sum[3]*fLmfR*(2.*costh1)); //cbar-c

      mesqij_NFZ[2]=fac*propZ*(gLpgR_sum[0]*fLpfR*(one + costh2)-gLmgR_sum[0]*fLmfR*(2.*costh1)); //d-dbar
      mesqij_NFZ[4]=fac*propZ*(gLpgR_sum[2]*fLpfR*(one + costh2)-gLmgR_sum[2]*fLmfR*(2.*costh1)); //s-sbar
      mesqij_NFZ[8]=fac*propZ*(gLpgR_sum[4]*fLpfR*(one + costh2)-gLmgR_sum[4]*fLmfR*(2.*costh1)); //b-bbar

      mesqij_NFZ[3]=fac*propZ*(gLpgR_sum[0]*fLpfR*(one + costh2)+gLmgR_sum[0]*fLmfR*(2.*costh1)); //dbar-d
      mesqij_NFZ[5]=fac*propZ*(gLpgR_sum[2]*fLpfR*(one + costh2)+gLmgR_sum[2]*fLmfR*(2.*costh1)); //sbar-s
      mesqij_NFZ[9]=fac*propZ*(gLpgR_sum[4]*fLpfR*(one + costh2)+gLmgR_sum[4]*fLmfR*(2.*costh1)); //bbar-b
    }
  if (opts.nproc == 3 && opts.useGamma && opts.order >= 3)
    {
      //                           Z                 gamma*                  interference
      mesqij_NFZ[0]=fac*(+(propZ*gLpgR_sum[1]*fLpfR + propG*aem2pi2*Q2_sum[1] - propZG*aem2pi*QgZp_sum[1]*(fLZ+fRZ)) * (one + costh2)
			 -(propZ*gLmgR_sum[1]*fLmfR                           - propZG*aem2pi*QgZm_sum[1]*(fLZ-fRZ)) * (2.*costh1));
      mesqij_NFZ[6]=fac*(+(propZ*gLpgR_sum[3]*fLpfR + propG*aem2pi2*Q2_sum[3] - propZG*aem2pi*QgZp_sum[3]*(fLZ+fRZ)) * (one + costh2)
			 -(propZ*gLmgR_sum[3]*fLmfR                           - propZG*aem2pi*QgZm_sum[3]*(fLZ-fRZ)) * (2.*costh1));

      mesqij_NFZ[1]=fac*(+(propZ*gLpgR_sum[1]*fLpfR + propG*aem2pi2*Q2_sum[1] - propZG*aem2pi*QgZp_sum[1]*(fLZ+fRZ)) * (one + costh2)
			 +(propZ*gLmgR_sum[1]*fLmfR                           - propZG*aem2pi*QgZm_sum[1]*(fLZ-fRZ)) * (2.*costh1));
      mesqij_NFZ[7]=fac*(+(propZ*gLpgR_sum[3]*fLpfR + propG*aem2pi2*Q2_sum[3] - propZG*aem2pi*QgZp_sum[3]*(fLZ+fRZ)) * (one + costh2)
			 +(propZ*gLmgR_sum[3]*fLmfR                           - propZG*aem2pi*QgZm_sum[3]*(fLZ-fRZ)) * (2.*costh1));

      mesqij_NFZ[2]=fac*(+(propZ*gLpgR_sum[0]*fLpfR + propG*aem2pi2*Q2_sum[0] - propZG*aem2pi*QgZp_sum[0]*(fLZ+fRZ)) * (one + costh2)
			 -(propZ*gLmgR_sum[0]*fLmfR                           - propZG*aem2pi*QgZm_sum[0]*(fLZ-fRZ)) * (2.*costh1));
      mesqij_NFZ[4]=fac*(+(propZ*gLpgR_sum[2]*fLpfR + propG*aem2pi2*Q2_sum[2] - propZG*aem2pi*QgZp_sum[2]*(fLZ+fRZ)) * (one + costh2)
			 -(propZ*gLmgR_sum[2]*fLmfR                           - propZG*aem2pi*QgZm_sum[2]*(fLZ-fRZ)) * (2.*costh1));
      mesqij_NFZ[8]=fac*(+(propZ*gLpgR_sum[4]*fLpfR + propG*aem2pi2*Q2_sum[4] - propZG*aem2pi*QgZp_sum[4]*(fLZ+fRZ)) * (one + costh2)
			 -(propZ*gLmgR_sum[4]*fLmfR                           - propZG*aem2pi*QgZm_sum[4]*(fLZ-fRZ)) * (2.*costh1));

      mesqij_NFZ[3]=fac*(+(propZ*gLpgR_sum[0]*fLpfR + propG*aem2pi2*Q2_sum[0] - propZG*aem2pi*QgZp_sum[0]*(fLZ+fRZ)) * (one + costh2)
			 +(propZ*gLmgR_sum[0]*fLmfR                           - propZG*aem2pi*QgZm_sum[0]*(fLZ-fRZ)) * (2.*costh1));
      mesqij_NFZ[5]=fac*(+(propZ*gLpgR_sum[2]*fLpfR + propG*aem2pi2*Q2_sum[2] - propZG*aem2pi*QgZp_sum[2]*(fLZ+fRZ)) * (one + costh2)
			 +(propZ*gLmgR_sum[2]*fLmfR                           - propZG*aem2pi*QgZm_sum[2]*(fLZ-fRZ)) * (2.*costh1));
      mesqij_NFZ[9]=fac*(+(propZ*gLpgR_sum[4]*fLpfR + propG*aem2pi2*Q2_sum[4] - propZG*aem2pi*QgZp_sum[4]*(fLZ+fRZ)) * (one + costh2)
			 +(propZ*gLmgR_sum[4]*fLmfR                           - propZG*aem2pi*QgZm_sum[4]*(fLZ-fRZ)) * (2.*costh1));
    }
}

//fortran interfaces
void initsigma_cpp_(double &m, double &cthmom0, double &cthmom1, double &cthmom2)
{
  //mass dependent part
  mesq::setpropagators(m);

  //  double one, costh1, costh2;
  //  one = 1.;
  //  costh1 = costh;
  //  costh2 = pow(costh,2);
  mesq::setmesq(cthmom0, cthmom1, cthmom2);

  //fortan array indices are inverted: [ub][u] means (u,ub)
  if (opts.nproc == 3)
    {
      sigmaij_.sigmaij_[ub][u ] = real(mesq::mesqij[0]);
      sigmaij_.sigmaij_[cb][c ] = real(mesq::mesqij[6]);
      sigmaij_.sigmaij_[u ][ub] = real(mesq::mesqij[1]);
      sigmaij_.sigmaij_[c ][cb] = real(mesq::mesqij[7]);
      sigmaij_.sigmaij_[db][d ] = real(mesq::mesqij[2]);
      sigmaij_.sigmaij_[sb][s ] = real(mesq::mesqij[4]);
      sigmaij_.sigmaij_[bb][b ] = real(mesq::mesqij[8]);
      sigmaij_.sigmaij_[d ][db] = real(mesq::mesqij[3]);
      sigmaij_.sigmaij_[s ][sb] = real(mesq::mesqij[5]);
      sigmaij_.sigmaij_[b ][bb] = real(mesq::mesqij[9]);
    }
  else if (opts.nproc == 1)
    {
      sigmaij_.sigmaij_[db][u ] = real(mesq::mesqij[0]);
      sigmaij_.sigmaij_[u ][db] = real(mesq::mesqij[1]);
      sigmaij_.sigmaij_[sb][u ] = real(mesq::mesqij[2]);
      sigmaij_.sigmaij_[u ][sb] = real(mesq::mesqij[3]);
      sigmaij_.sigmaij_[bb][u ] = real(mesq::mesqij[4]);
      sigmaij_.sigmaij_[u ][bb] = real(mesq::mesqij[5]);
      sigmaij_.sigmaij_[sb][c ] = real(mesq::mesqij[6]);
      sigmaij_.sigmaij_[c ][sb] = real(mesq::mesqij[7]);
      sigmaij_.sigmaij_[db][c ] = real(mesq::mesqij[8]);
      sigmaij_.sigmaij_[c ][db] = real(mesq::mesqij[9]);
      sigmaij_.sigmaij_[bb][c ] = real(mesq::mesqij[10]);
      sigmaij_.sigmaij_[c ][bb] = real(mesq::mesqij[11]);
    }
  else if (opts.nproc == 2)
    {
      sigmaij_.sigmaij_[ub][d ] = real(mesq::mesqij[0]);
      sigmaij_.sigmaij_[d ][ub] = real(mesq::mesqij[1]);
      sigmaij_.sigmaij_[ub][s ] = real(mesq::mesqij[2]);
      sigmaij_.sigmaij_[s ][ub] = real(mesq::mesqij[3]);
      sigmaij_.sigmaij_[ub][b ] = real(mesq::mesqij[4]);
      sigmaij_.sigmaij_[b ][ub] = real(mesq::mesqij[5]);
      sigmaij_.sigmaij_[cb][s ] = real(mesq::mesqij[6]);
      sigmaij_.sigmaij_[s ][cb] = real(mesq::mesqij[7]);
      sigmaij_.sigmaij_[cb][d ] = real(mesq::mesqij[8]);
      sigmaij_.sigmaij_[d ][cb] = real(mesq::mesqij[9]);
      sigmaij_.sigmaij_[cb][b ] = real(mesq::mesqij[10]);
      sigmaij_.sigmaij_[b ][cb] = real(mesq::mesqij[11]);
    }
}

//rapidity differential
double mesq::loxs(double x1, double x2, double muf)
{
  double xs = 0.;

  //PDFs
  double fx1[2*MAXNF+1],fx2[2*MAXNF+1];
  fdist_(opts.ih1,x1,muf,fx1);
  fdist_(opts.ih2,x2,muf,fx2);
	  
  for (int sp = 0; sp < totpch; sp++)
    xs += fx1[pid1[sp]]*fx2[pid2[sp]]*real(mesqij[sp]);

  return xs;
}
  
 //rapidity integrated
double mesq::loxs(double tau, double muf)
{
  double xs = 0.;

  double fx1[2*MAXNF+1],fx2[2*MAXNF+1];
  for (int i=0; i < opts.yintervals; i++)
    {
      double ya = phasespace::ymin+(phasespace::ymax-phasespace::ymin)*i/opts.yintervals;
      double yb = phasespace::ymin+(phasespace::ymax-phasespace::ymin)*(i+1)/opts.yintervals;
      double xc = 0.5*(ya+yb);
      double xm = 0.5*(yb-ya);
      for (int j=0; j < opts.yrule; j++)
	{
	  double y = xc+xm*gr::xxx[opts.yrule-1][j];
	  double exppy = exp(y);
	  double expmy = 1./exppy;
	  double x1 = tau*exppy;
	  double x2 = tau*expmy;
		      
	  //PDFs
	  fdist_(opts.ih1,x1,muf,fx1);
	  fdist_(opts.ih2,x2,muf,fx2);
		      
	  for (int sp = 0; sp < totpch; sp++)
	    xs += fx1[pid1[sp]]*fx2[pid2[sp]]*real(mesqij[sp])*gr::www[opts.yrule-1][j]*xm;;
	}
    }
  return xs;
}
