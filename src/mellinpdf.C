#include "mellinpdf.h"
#include "gaussrules.h"
#include "clenshawcurtisrules.h"
#include "settings.h"
#include "mellinint.h"
#include "pdf.h"
#include "parton.h"
#ifdef USECERES
#include "cerespdf.h"
#endif
#include "scales.h"
#include "phasespace.h"
#include "clock_real.h"

#include <complex>
#include <cmath>
#include <iostream>
#include <string.h>

using namespace parton;

//Cache PDFs, t values, jacobian, kernels
double *mellinpdf::t;              //gauss nodes
double *mellinpdf::fac;            //overall factor
complex <double> *mellinpdf::kern; //kernel of the Mellin transform

double *mellinpdf::t_1;
double *mellinpdf::fac_1;
complex <double> *mellinpdf::kern_1;
double *mellinpdf::t_2;
double *mellinpdf::fac_2;
complex <double> *mellinpdf::kern_2;

//Mellin moments of PDFs
complex <double> *mellinpdf::UP;
complex <double> *mellinpdf::DO;
complex <double> *mellinpdf::ST;
complex <double> *mellinpdf::CH;
complex <double> *mellinpdf::BO;
complex <double> *mellinpdf::GL;
complex <double> *mellinpdf::UB;
complex <double> *mellinpdf::DB;
complex <double> *mellinpdf::SB;
complex <double> *mellinpdf::CB;
complex <double> *mellinpdf::BB;

complex <double> *mellinpdf::UP_1;
complex <double> *mellinpdf::DO_1;
complex <double> *mellinpdf::ST_1;
complex <double> *mellinpdf::CH_1;
complex <double> *mellinpdf::BO_1;
complex <double> *mellinpdf::GL_1;
complex <double> *mellinpdf::UB_1;
complex <double> *mellinpdf::DB_1;
complex <double> *mellinpdf::SB_1;
complex <double> *mellinpdf::CB_1;
complex <double> *mellinpdf::BB_1;
complex <double> *mellinpdf::UP_2;
complex <double> *mellinpdf::DO_2;
complex <double> *mellinpdf::ST_2;
complex <double> *mellinpdf::CH_2;
complex <double> *mellinpdf::BO_2;
complex <double> *mellinpdf::GL_2;
complex <double> *mellinpdf::UB_2;
complex <double> *mellinpdf::DB_2;
complex <double> *mellinpdf::SB_2;
complex <double> *mellinpdf::CB_2;
complex <double> *mellinpdf::BB_2;

//Values of the PDFs at the nodes of the Gauss rule
double *mellinpdf::fup;
double *mellinpdf::fdo;
double *mellinpdf::fst;
double *mellinpdf::fch;
double *mellinpdf::fbo;
double *mellinpdf::fgl;
double *mellinpdf::fub;
double *mellinpdf::fdb;
double *mellinpdf::fsb;
double *mellinpdf::fcb;
double *mellinpdf::fbb;

double *mellinpdf::fup_1;
double *mellinpdf::fdo_1;
double *mellinpdf::fst_1;
double *mellinpdf::fch_1;
double *mellinpdf::fbo_1;
double *mellinpdf::fgl_1;
double *mellinpdf::fub_1;
double *mellinpdf::fdb_1;
double *mellinpdf::fsb_1;
double *mellinpdf::fcb_1;
double *mellinpdf::fbb_1;
double *mellinpdf::fup_2;
double *mellinpdf::fdo_2;
double *mellinpdf::fst_2;
double *mellinpdf::fch_2;
double *mellinpdf::fbo_2;
double *mellinpdf::fgl_2;
double *mellinpdf::fub_2;
double *mellinpdf::fdb_2;
double *mellinpdf::fsb_2;
double *mellinpdf::fcb_2;
double *mellinpdf::fbb_2;

//Allocate thread private memory
void mellinpdf::allocate()
{
  if (opts.mellin1d)
    {
      t = new double [opts.pdfrule];
      fac = new double [opts.pdfrule];
      kern = new complex <double> [opts.pdfrule*mellinint::mdim];
    }
  else
    {
      t_1 = new double [opts.pdfrule];
      fac_1 = new double [opts.pdfrule];
      kern_1 = new complex <double> [opts.pdfrule*mellinint::mdim];
      t_2 = new double [opts.pdfrule];
      fac_2 = new double [opts.pdfrule];
      kern_2 = new complex <double> [opts.pdfrule*mellinint::mdim];
    }
  
  if (opts.mellin1d)
    {
      UP = new complex <double> [mellinint::mdim];
      DO = new complex <double> [mellinint::mdim];
      ST = new complex <double> [mellinint::mdim];
      CH = new complex <double> [mellinint::mdim];
      BO = new complex <double> [mellinint::mdim];
      GL = new complex <double> [mellinint::mdim];
      UB = new complex <double> [mellinint::mdim];
      DB = new complex <double> [mellinint::mdim];
      SB = new complex <double> [mellinint::mdim];
      CB = new complex <double> [mellinint::mdim];
      BB = new complex <double> [mellinint::mdim];

      fup = new double [opts.pdfrule];
      fdo = new double [opts.pdfrule];
      fst = new double [opts.pdfrule];
      fch = new double [opts.pdfrule];
      fbo = new double [opts.pdfrule];
      fgl = new double [opts.pdfrule];
      fub = new double [opts.pdfrule];
      fdb = new double [opts.pdfrule];
      fsb = new double [opts.pdfrule];
      fcb = new double [opts.pdfrule];
      fbb = new double [opts.pdfrule];
    }
  else
    {
      UP_1 = new complex <double> [mellinint::mdim];
      DO_1 = new complex <double> [mellinint::mdim];
      ST_1 = new complex <double> [mellinint::mdim];
      CH_1 = new complex <double> [mellinint::mdim];
      BO_1 = new complex <double> [mellinint::mdim];
      GL_1 = new complex <double> [mellinint::mdim];
      UB_1 = new complex <double> [mellinint::mdim];
      DB_1 = new complex <double> [mellinint::mdim];
      SB_1 = new complex <double> [mellinint::mdim];
      CB_1 = new complex <double> [mellinint::mdim];
      BB_1 = new complex <double> [mellinint::mdim];
      UP_2 = new complex <double> [mellinint::mdim];
      DO_2 = new complex <double> [mellinint::mdim];
      ST_2 = new complex <double> [mellinint::mdim];
      CH_2 = new complex <double> [mellinint::mdim];
      BO_2 = new complex <double> [mellinint::mdim];
      GL_2 = new complex <double> [mellinint::mdim];
      UB_2 = new complex <double> [mellinint::mdim];
      DB_2 = new complex <double> [mellinint::mdim];
      SB_2 = new complex <double> [mellinint::mdim];
      CB_2 = new complex <double> [mellinint::mdim];
      BB_2 = new complex <double> [mellinint::mdim];

      fup_1 = new double [opts.pdfrule];
      fdo_1 = new double [opts.pdfrule];
      fst_1 = new double [opts.pdfrule];
      fch_1 = new double [opts.pdfrule];
      fbo_1 = new double [opts.pdfrule];
      fgl_1 = new double [opts.pdfrule];
      fub_1 = new double [opts.pdfrule];
      fdb_1 = new double [opts.pdfrule];
      fsb_1 = new double [opts.pdfrule];
      fcb_1 = new double [opts.pdfrule];
      fbb_1 = new double [opts.pdfrule];
      fup_2 = new double [opts.pdfrule];
      fdo_2 = new double [opts.pdfrule];
      fst_2 = new double [opts.pdfrule];
      fch_2 = new double [opts.pdfrule];
      fbo_2 = new double [opts.pdfrule];
      fgl_2 = new double [opts.pdfrule];
      fub_2 = new double [opts.pdfrule];
      fdb_2 = new double [opts.pdfrule];
      fsb_2 = new double [opts.pdfrule];
      fcb_2 = new double [opts.pdfrule];
      fbb_2 = new double [opts.pdfrule];
    }
}

//Release thread private memory
void mellinpdf::free()
{
  if (opts.mellin1d)
    {
      delete [] t;
      delete [] fac;
      delete [] kern;
    }
  else
    {
      delete [] t_1;
      delete [] fac_1;
      delete [] kern_1;
      delete [] t_2;
      delete [] fac_2;
      delete [] kern_2;
    }
  if (opts.mellin1d)
    {
      delete [] UP;
      delete [] DO;
      delete [] ST;
      delete [] CH;
      delete [] BO;
      delete [] GL;
      delete [] UB;
      delete [] DB;
      delete [] SB;
      delete [] CB;
      delete [] BB;
      delete [] fup;
      delete [] fdo;
      delete [] fst;
      delete [] fch;
      delete [] fbo;
      delete [] fgl;
      delete [] fub;
      delete [] fdb;
      delete [] fsb;
      delete [] fcb;
      delete [] fbb;
    }
  else
    {
      delete [] UP_1;
      delete [] DO_1;
      delete [] ST_1;
      delete [] CH_1;
      delete [] BO_1;
      delete [] GL_1;
      delete [] UB_1;
      delete [] DB_1;
      delete [] SB_1;
      delete [] CB_1;
      delete [] BB_1;
      delete [] UP_2;
      delete [] DO_2;
      delete [] ST_2;
      delete [] CH_2;
      delete [] BO_2;
      delete [] GL_2;
      delete [] UB_2;
      delete [] DB_2;
      delete [] SB_2;
      delete [] CB_2;
      delete [] BB_2;
      delete [] fup_1;
      delete [] fdo_1;
      delete [] fst_1;
      delete [] fch_1;
      delete [] fbo_1;
      delete [] fgl_1;
      delete [] fub_1;
      delete [] fdb_1;
      delete [] fsb_1;
      delete [] fcb_1;
      delete [] fbb_1;
      delete [] fup_2;
      delete [] fdo_2;
      delete [] fst_2;
      delete [] fch_2;
      delete [] fbo_2;
      delete [] fgl_2;
      delete [] fub_2;
      delete [] fdb_2;
      delete [] fsb_2;
      delete [] fcb_2;
      delete [] fbb_2;
    }
}

//Release non-thread-private memory
void mellinpdf::release()
{
#ifdef USECERES
  cerespdf::free();
#endif
}

//set up the x-spacing
void mellinpdf::init()
{
  // boundaries of integration
  //double xmin = 1e-8;
  //double xmin = pdf::xmin;
  double xmin = pow(bins.mbins.front()/opts.sroot,2); //Restrict the integration of moments to xmin = m/sqrt(s)*exp(-ymax) = (m/sqrt(s))^2
  /*
  double xmax = 1;
  double ll = log(xmax/xmin);
  double cx = 0.5;
  double mx = 0.5;
  for (int i = 0; i < opts.pdfrule; i++)
    {
      double x = cx+mx*gr::xxx[opts.pdfrule-1][i];
      if (opts.mellin1d)
	{
	  t[i] = xmin*exp(ll*x);
	  fac[i] = mx * t[i] * ll * gr::www[opts.pdfrule-1][i];
	}
      else
	{
	  t_1[i] = xmin*exp(ll*x);
	  fac_1[i] = mx * t_1[i] * ll * gr::www[opts.pdfrule-1][i];
	  t_2[i] = xmin*exp(ll*x);
	  fac_2[i] = mx * t_2[i] * ll * gr::www[opts.pdfrule-1][i];
	}
    }
  */

#ifdef USECERES
  //Approximate with ceres
  if (opts.mellintr == 2)
    cerespdf::init(xmin);
#endif
}

//Evaluate PDFs in x-space at the required factorisation scale
void mellinpdf::evalpdfs(double scale, double m, double y)
{
  //There is a problem here when evalpdfs is called without m and y from pegasus, because t,fac,kern are not set!!!
  //--> Should always call evalpdfs with some values of m and y
  clock_t begin_time, end_time;
  begin_time = clock();  

  //factorization scale
  double muf = scale;

  if (muf < 0)
    muf = scales::fac;

  if (m < 0)
    {
      m = bins.mbins.front();
      y = 0.;
    }

  if (opts.melup == 0)
    {
      m = bins.mbins.front();
      //y = 0.;
      double ylim = log(opts.sroot/bins.mbins.back());  //here take bins.mbins.back() so that this is the smallest ylim --> avoid z < 0
      double ymn = min(max(-ylim, bins.ybins.front()),ylim);
      double ymx = max(min(ylim, bins.ybins.back()),-ylim);
      y = (ymn+ymx)/2.;
    }
  if (opts.melup == 1)
    {
      //m = sqrt(phasespace::mmin*phasespace::mmax);
      //double ylim = log(opts.sroot/m);
      m = phasespace::mmin;
      double ylim = log(opts.sroot/phasespace::mmax);
      double ymn = min(max(-ylim, phasespace::ymin),ylim);
      double ymx = max(min(ylim, phasespace::ymax),-ylim);
      y = (ymn+ymx)/2.;
    }
  if (opts.melup == 2)
    {
      m = phasespace::m;
      //y = phasespace::y;
      double ylim = log(opts.sroot/phasespace::mmax);
      double ymn = min(max(-ylim, phasespace::ymin),ylim);
      double ymx = max(min(ylim, phasespace::ymax),-ylim);
      y = (ymn+ymx)/2.;
      if (phasespace::ymax >= ylim && phasespace::ymin > -ylim)
	y = ymn;
      if (phasespace::ymin <= -ylim && phasespace::ymax < ylim)
	y = ymx;
    }

  //y = 0.;
  if (opts.mellin1d)
    {
      int p = 2; //2;
      //double xmin = 1e-8;
      //double xmin = pdf::xmin;
      double xmin = pow(m/opts.sroot,p);
      double ll = log(1./xmin);
      for (int i = 0; i < opts.pdfrule; i++)
	{
	  double x = 0.5+0.5*gr::xxx[opts.pdfrule-1][i];
	  t[i] = xmin*exp(ll*x);
	  fac[i] = 0.5 * t[i] * ll * gr::www[opts.pdfrule-1][i];
	}
    }
  else
    {
      int p = 4;
      double xmin_1 = pow(m/opts.sroot*exp(y),p);
      double xmin_2 = pow(m/opts.sroot*exp(-y),p);
      //double xmin_1 = pow(m/opts.sroot*exp(phasespace::ymin),p);
      //double xmin_2 = pow(m/opts.sroot*exp(-phasespace::ymax),p);
      
      double ll_1 = log(1./xmin_1);
      double ll_2 = log(1./xmin_2);
      for (int i = 0; i < opts.pdfrule; i++)
	{
	  double x = 0.5+0.5*gr::xxx[opts.pdfrule-1][i];
	  t_1[i] = xmin_1*exp(ll_1*x);
	  t_2[i] = xmin_2*exp(ll_2*x);
	  fac_1[i] = 0.5 * t_1[i] * ll_1 * gr::www[opts.pdfrule-1][i];
	  fac_2[i] = 0.5 * t_2[i] * ll_2 * gr::www[opts.pdfrule-1][i];
	}
    }
  
#ifdef USECERES
  if (opts.mellintr == 2)
    {
      cerespdf::update(muf);
      return;
    }
#endif
  
  //else if (opts.mellintr == 0 || opts.mellintr == 1)
  
  //Call the PDFs and store values at the knots
  if (opts.mellin1d)
    {
      double fx[11];
      for (int i = 0; i < opts.pdfrule; i++)
	{
	  fdist_(opts.ih1,t[i],muf,fx);
	  fup[i] = fac[i] * fx[U];
	  fdo[i] = fac[i] * fx[D];
	  fst[i] = fac[i] * fx[S];
	  fch[i] = fac[i] * fx[C];
	  fbo[i] = fac[i] * fx[B];
	  fgl[i] = fac[i] * fx[G];
	  fub[i] = fac[i] * fx[Ub];
	  fdb[i] = fac[i] * fx[Db];
	  fsb[i] = fac[i] * fx[Sb];
	  fcb[i] = fac[i] * fx[Cb];
	  fbb[i] = fac[i] * fx[Bb];
	  //cout << i << "  " << t[i] << "  " << fac[i] << "  " << fx[G] << endl;
	}
    }
  else
    {
      double fx_1[11];
      double fx_2[11];
      for (int i = 0; i < opts.pdfrule; i++)
	{
	  fdist_(opts.ih1,t_1[i],muf,fx_1);
	  fup_1[i] = fac_1[i] * fx_1[U];
	  fdo_1[i] = fac_1[i] * fx_1[D];
	  fst_1[i] = fac_1[i] * fx_1[S];
	  fch_1[i] = fac_1[i] * fx_1[C];
	  fbo_1[i] = fac_1[i] * fx_1[B];
	  fgl_1[i] = fac_1[i] * fx_1[G];
	  fub_1[i] = fac_1[i] * fx_1[Ub];
	  fdb_1[i] = fac_1[i] * fx_1[Db];
	  fsb_1[i] = fac_1[i] * fx_1[Sb];
	  fcb_1[i] = fac_1[i] * fx_1[Cb];
	  fbb_1[i] = fac_1[i] * fx_1[Bb];

	  fdist_(opts.ih1,t_2[i],muf,fx_2);
	  fup_2[i] = fac_2[i] * fx_2[U];
	  fdo_2[i] = fac_2[i] * fx_2[D];
	  fst_2[i] = fac_2[i] * fx_2[S];
	  fch_2[i] = fac_2[i] * fx_2[C];
	  fbo_2[i] = fac_2[i] * fx_2[B];
	  fgl_2[i] = fac_2[i] * fx_2[G];
	  fub_2[i] = fac_2[i] * fx_2[Ub];
	  fdb_2[i] = fac_2[i] * fx_2[Db];
	  fsb_2[i] = fac_2[i] * fx_2[Sb];
	  fcb_2[i] = fac_2[i] * fx_2[Cb];
	  fbb_2[i] = fac_2[i] * fx_2[Bb];
	}
    }

  end_time = clock();  
  //cout << "muf is " << muf << " PDFs evaluation in x-space done in "  << float(end_time - begin_time)/CLOCKS_PER_SEC*1000. << "ms" << endl;
}

//Update the the N points for the Mellin inversion
void mellinpdf::update_mellin()
{
  if (opts.mellin1d)
    for (int i = 0; i < opts.pdfrule; i++)
      for (int m = 0; m < mellinint::mdim; m++)
	kern[i*mellinint::mdim+m] = pow(t[i], mellinint::Np[m]-1.);
  else
    for (int i = 0; i < opts.pdfrule; i++)
      for (int m = 0; m < mellinint::mdim; m++)
	{
	  kern_1[i*mellinint::mdim+m] = pow(t_1[i], mellinint::Np_1[m]-1.);
	  kern_2[i*mellinint::mdim+m] = pow(t_2[i], mellinint::Np_2[m]-1.);
	}
}

//Update PDFs in N space //void mellinpdf::update(double muf)
//Transform PDFs from x to N
void mellinpdf::transform()
{
  clock_t begin_time, end_time;
  begin_time = clock();  

  //evalpdfs();
  switch (opts.mellintr)
    {
    case 0: gauss_quad(); break;
    case 1: laguerre_ipol(); break;
#ifdef USECERES
    case 2: ceres_pdf(); break;
#endif
    }
  
  end_time = clock();  
  //cout << "muf is " << muf << " x to N done in "  << float(end_time - begin_time)/CLOCKS_PER_SEC*1000. << "ms" << endl;
}

//Computes the complex N Mellin moments of pdfs by means of Gauss-Legendre quadrature
void mellinpdf::gauss_quad()
{
  if (opts.mellin1d)
    {
      fill(UP,UP+mellinint::mdim, 0.);
      fill(DO,DO+mellinint::mdim, 0.);
      fill(ST,ST+mellinint::mdim, 0.);
      fill(CH,CH+mellinint::mdim, 0.);
      fill(BO,BO+mellinint::mdim, 0.);
      fill(GL,GL+mellinint::mdim, 0.);
      fill(UB,UB+mellinint::mdim, 0.);
      fill(DB,DB+mellinint::mdim, 0.);
      fill(SB,SB+mellinint::mdim, 0.);
      fill(CB,CB+mellinint::mdim, 0.);
      fill(BB,BB+mellinint::mdim, 0.);
  
      //Calculate Mellin moments as:
      //integral_0^1{ x^(N-1) fx dx}
      for (int i = 0; i < opts.pdfrule; i++)
	for (int m = 0; m < mellinint::mdim; m++)
	  {
	    UP[m] += fup[i] * kern[i*mellinint::mdim+m];
	    DO[m] += fdo[i] * kern[i*mellinint::mdim+m];
	    ST[m] += fst[i] * kern[i*mellinint::mdim+m];
	    CH[m] += fch[i] * kern[i*mellinint::mdim+m];
	    BO[m] += fbo[i] * kern[i*mellinint::mdim+m];
	    GL[m] += fgl[i] * kern[i*mellinint::mdim+m];
	    UB[m] += fub[i] * kern[i*mellinint::mdim+m];
	    DB[m] += fdb[i] * kern[i*mellinint::mdim+m];
	    SB[m] += fsb[i] * kern[i*mellinint::mdim+m];
	    CB[m] += fcb[i] * kern[i*mellinint::mdim+m];
	    BB[m] += fbb[i] * kern[i*mellinint::mdim+m];
	    //cout << i << "  " << kern[i*mellinint::mdim+m] << endl;
	    //if (m == 0) cout << m << "  " << i << "  " << GL[m] << endl;
	  }
    }
  else
    {
      fill(UP_1,UP_1+mellinint::mdim, 0.);
      fill(DO_1,DO_1+mellinint::mdim, 0.);
      fill(ST_1,ST_1+mellinint::mdim, 0.);
      fill(CH_1,CH_1+mellinint::mdim, 0.);
      fill(BO_1,BO_1+mellinint::mdim, 0.);
      fill(GL_1,GL_1+mellinint::mdim, 0.);
      fill(UB_1,UB_1+mellinint::mdim, 0.);
      fill(DB_1,DB_1+mellinint::mdim, 0.);
      fill(SB_1,SB_1+mellinint::mdim, 0.);
      fill(CB_1,CB_1+mellinint::mdim, 0.);
      fill(BB_1,BB_1+mellinint::mdim, 0.);
      fill(UP_2,UP_2+mellinint::mdim, 0.);
      fill(DO_2,DO_2+mellinint::mdim, 0.);
      fill(ST_2,ST_2+mellinint::mdim, 0.);
      fill(CH_2,CH_2+mellinint::mdim, 0.);
      fill(BO_2,BO_2+mellinint::mdim, 0.);
      fill(GL_2,GL_2+mellinint::mdim, 0.);
      fill(UB_2,UB_2+mellinint::mdim, 0.);
      fill(DB_2,DB_2+mellinint::mdim, 0.);
      fill(SB_2,SB_2+mellinint::mdim, 0.);
      fill(CB_2,CB_2+mellinint::mdim, 0.);
      fill(BB_2,BB_2+mellinint::mdim, 0.);
  
      //Calculate Mellin moments as:
      //integral_0^1{ x^(N-1) fx dx}
      for (int i = 0; i < opts.pdfrule; i++)
	for (int m = 0; m < mellinint::mdim; m++)
	  {
	    UP_1[m] += fup_1[i] * kern_1[i*mellinint::mdim+m];
	    DO_1[m] += fdo_1[i] * kern_1[i*mellinint::mdim+m];
	    ST_1[m] += fst_1[i] * kern_1[i*mellinint::mdim+m];
	    CH_1[m] += fch_1[i] * kern_1[i*mellinint::mdim+m];
	    BO_1[m] += fbo_1[i] * kern_1[i*mellinint::mdim+m];
	    GL_1[m] += fgl_1[i] * kern_1[i*mellinint::mdim+m];
	    UB_1[m] += fub_1[i] * kern_1[i*mellinint::mdim+m];
	    DB_1[m] += fdb_1[i] * kern_1[i*mellinint::mdim+m];
	    SB_1[m] += fsb_1[i] * kern_1[i*mellinint::mdim+m];
	    CB_1[m] += fcb_1[i] * kern_1[i*mellinint::mdim+m];
	    BB_1[m] += fbb_1[i] * kern_1[i*mellinint::mdim+m];
	    UP_2[m] += fup_2[i] * kern_2[i*mellinint::mdim+m];
	    DO_2[m] += fdo_2[i] * kern_2[i*mellinint::mdim+m];
	    ST_2[m] += fst_2[i] * kern_2[i*mellinint::mdim+m];
	    CH_2[m] += fch_2[i] * kern_2[i*mellinint::mdim+m];
	    BO_2[m] += fbo_2[i] * kern_2[i*mellinint::mdim+m];
	    GL_2[m] += fgl_2[i] * kern_2[i*mellinint::mdim+m];
	    UB_2[m] += fub_2[i] * kern_2[i*mellinint::mdim+m];
	    DB_2[m] += fdb_2[i] * kern_2[i*mellinint::mdim+m];
	    SB_2[m] += fsb_2[i] * kern_2[i*mellinint::mdim+m];
	    CB_2[m] += fcb_2[i] * kern_2[i*mellinint::mdim+m];
	    BB_2[m] += fbb_2[i] * kern_2[i*mellinint::mdim+m];
	  }
    }
}

void laguerre(int n, double x, double *L)
{
  L[0] = 1.;
  double lkm2 = 0.; //--> L[-1]
  double lkm1 = L[0];// = 1.;
  for (int k = 1; k <= n; k++)
    {
      L[k] = ((2*k-1.-x)*lkm1-(k-1.)*lkm2)/k;
      lkm2=lkm1;
      lkm1=L[k];
    }
  return;
}

//Computes the complex N Mellin moments of pdfs by means of Laguerre interpolation and Laplace transform
void mellinpdf::laguerre_ipol()
{
  //order of the Laguerre approximation --> make this an option
  int N = 150;

  //scaling x -> x^p
  double p = 2.;

  //Calculate Laguerre projection
  double lup[N+1] = {0.};
  double ldo[N+1] = {0.};
  double lst[N+1] = {0.};
  double lch[N+1] = {0.};
  double lbo[N+1] = {0.};
  double lgl[N+1] = {0.};
  double lub[N+1] = {0.};
  double ldb[N+1] = {0.};
  double lsb[N+1] = {0.};
  double lcb[N+1] = {0.};
  double lbb[N+1] = {0.};
  
  double cx = 0.5;
  double mx = 0.5;

  //calculate the x^alpha*(1-x)^beta approximation
  //calculate the 0, 1, 2 Mellin moments
  double fnup[3] = {0.};
  double fndo[3] = {0.};
  double fnst[3] = {0.};
  double fnch[3] = {0.};
  double fnbo[3] = {0.};
  double fngl[3] = {0.};
  double fnub[3] = {0.};
  double fndb[3] = {0.};
  double fnsb[3] = {0.};
  double fncb[3] = {0.};
  double fnbb[3] = {0.};


  for (int i = 0; i < opts.pdfrule; i++)
    {
      //double x = cx+mx*xxx[rule-1][i];
      //double t = a*pow(1./a,x);
      //double jac = 0.5* t * log(1./a);
      //f0 += fknots[i]      * fac[i];
      //f1 += fknots[i] *t   * fac[i];
      //f2 += fknots[i] *t*t * fac[i];
      for (int n = 0; n < 3; n++)
	{
	  fnup[n] += fup[i] * pow(t[i],n) *fac[i];
	  fndo[n] += fdo[i] * pow(t[i],n) *fac[i];
	  fnst[n] += fst[i] * pow(t[i],n) *fac[i];
	  fnch[n] += fch[i] * pow(t[i],n) *fac[i];
	  fnbo[n] += fbo[i] * pow(t[i],n) *fac[i];
	  fngl[n] += fgl[i] * pow(t[i],n) *fac[i];
	  fnub[n] += fub[i] * pow(t[i],n) *fac[i];
	  fndb[n] += fdb[i] * pow(t[i],n) *fac[i];
	  fnsb[n] += fsb[i] * pow(t[i],n) *fac[i];
	  fncb[n] += fcb[i] * pow(t[i],n) *fac[i];
	  fnbb[n] += fbb[i] * pow(t[i],n) *fac[i];
	}
    }

  //This is used to regularise the PDFs with the same x^alpha behaviour at low x
  double alpha = 2;
  double powup = alpha - fnup[1]*(fnup[2]-fnup[1])/(fnup[1]*fnup[1]-fnup[0]*fnup[2]);
  double powdo = alpha - fndo[1]*(fndo[2]-fndo[1])/(fndo[1]*fndo[1]-fndo[0]*fndo[2]);
  double powst = alpha - fnst[1]*(fnst[2]-fnst[1])/(fnst[1]*fnst[1]-fnst[0]*fnst[2]);
  double powch = alpha - fnch[1]*(fnch[2]-fnch[1])/(fnch[1]*fnch[1]-fnch[0]*fnch[2]);
  double powbo = alpha - fnbo[1]*(fnbo[2]-fnbo[1])/(fnbo[1]*fnbo[1]-fnbo[0]*fnbo[2]);
  double powgl = alpha - fngl[1]*(fngl[2]-fngl[1])/(fngl[1]*fngl[1]-fngl[0]*fngl[2]);
  double powub = alpha - fnub[1]*(fnub[2]-fnub[1])/(fnub[1]*fnub[1]-fnub[0]*fnub[2]);
  double powdb = alpha - fndb[1]*(fndb[2]-fndb[1])/(fndb[1]*fndb[1]-fndb[0]*fndb[2]);
  double powsb = alpha - fnsb[1]*(fnsb[2]-fnsb[1])/(fnsb[1]*fnsb[1]-fnsb[0]*fnsb[2]);
  double powcb = alpha - fncb[1]*(fncb[2]-fncb[1])/(fncb[1]*fncb[1]-fncb[0]*fncb[2]);
  double powbb = alpha - fnbb[1]*(fnbb[2]-fnbb[1])/(fnbb[1]*fnbb[1]-fnbb[0]*fnbb[2]);

  //The largest poles of the Mellin transform are at these values of z:
  //cout << powup << "  "
  //     << powdo << "  "
  //     << powub << "  "
  //     << powdb << "  "
  //     << powst << "  "
  //     << powsb << "  "
  //     << powgl << "  "
  //     << powch << "  "
  //     << powbo << endl;

  
  //double alpha = 2;
  for (int i = 0; i < opts.pdfrule; i++)
    {
      double L[N+1];
      double u = -log(t[i])/p;
      laguerre(N, u, L);// --> cache it as L[n][i]
      //      cout << endl << i << "  " << u << endl;
      //double po = pow(t,alpha);
      for (int n = 0; n <= N; n++)
	{
	  lup[n] += fup[i] * L[n] * pow(t[i],1./p-1.+powup);
	  ldo[n] += fdo[i] * L[n] * pow(t[i],1./p-1.+powdo);
	  lst[n] += fst[i] * L[n] * pow(t[i],1./p-1.+powst);
	  lch[n] += fch[i] * L[n] * pow(t[i],1./p-1.+powch);
	  lbo[n] += fbo[i] * L[n] * pow(t[i],1./p-1.+powbo);
	  lgl[n] += fgl[i] * L[n] * pow(t[i],1./p-1.+powgl);
	  lub[n] += fub[i] * L[n] * pow(t[i],1./p-1.+powub);
	  ldb[n] += fdb[i] * L[n] * pow(t[i],1./p-1.+powdb);
	  lsb[n] += fsb[i] * L[n] * pow(t[i],1./p-1.+powsb);
	  lcb[n] += fcb[i] * L[n] * pow(t[i],1./p-1.+powcb);
	  lbb[n] += fbb[i] * L[n] * pow(t[i],1./p-1.+powbb);

//	  cout << i << "  " << n << "  " << L[n] << endl;
	}
    }
  //for (int n = 0; n <= N; n++)
  //printf("Laguerre projection %d \t (%e) \n", n, lgl[n]);
  
  //initialise the moments
  for (int m = 0; m < mellinint::mdim; m++)
    {
      UP[m] = 0.;
      DO[m] = 0.;
      ST[m] = 0.;
      CH[m] = 0.;
      BO[m] = 0.;
      GL[m] = 0.;
      UB[m] = 0.;
      DB[m] = 0.;
      SB[m] = 0.;
      CB[m] = 0.;
      BB[m] = 0.;
    }
  //Calculate Mellin moments from the
  //Mellin transform of Laguerre polynomials
  for (int m = 0; m < mellinint::mdim; m++)
    {
      //complex <double> z = (mellinint::Np[m] - alpha)*p;
      complex <double> zup = (mellinint::Np[m] - powup)*p;
      complex <double> zdo = (mellinint::Np[m] - powdo)*p;
      complex <double> zst = (mellinint::Np[m] - powst)*p;
      complex <double> zch = (mellinint::Np[m] - powch)*p;
      complex <double> zbo = (mellinint::Np[m] - powbo)*p;
      complex <double> zgl = (mellinint::Np[m] - powgl)*p;
      complex <double> zub = (mellinint::Np[m] - powub)*p;
      complex <double> zdb = (mellinint::Np[m] - powdb)*p;
      complex <double> zsb = (mellinint::Np[m] - powsb)*p;
      complex <double> zcb = (mellinint::Np[m] - powcb)*p;
      complex <double> zbb = (mellinint::Np[m] - powbb)*p;
      for (int n = 0; n <= N; n++)
	{
	  //complex <double> lapl = 1./z * pow(((z-1.)/z),n);
	  //complex <double> lapl = 1./z * pow(1.-1./z,n); //cache this as lapl[m][n]

	  complex <double> laplup = 1./zup * pow(1.-1./zup,n); //cache this as laplup[m][n]
	  complex <double> lapldo = 1./zdo * pow(1.-1./zdo,n); //cache this as lapldo[m][n]
	  complex <double> laplst = 1./zst * pow(1.-1./zst,n); //cache this as laplst[m][n]
	  complex <double> laplch = 1./zch * pow(1.-1./zch,n); //cache this as laplch[m][n]
	  complex <double> laplbo = 1./zbo * pow(1.-1./zbo,n); //cache this as laplbo[m][n]
	  complex <double> laplgl = 1./zgl * pow(1.-1./zgl,n); //cache this as laplgl[m][n]
	  complex <double> laplub = 1./zub * pow(1.-1./zub,n); //cache this as laplub[m][n]
	  complex <double> lapldb = 1./zdb * pow(1.-1./zdb,n); //cache this as lapldb[m][n]
	  complex <double> laplsb = 1./zsb * pow(1.-1./zsb,n); //cache this as laplsb[m][n]
	  complex <double> laplcb = 1./zcb * pow(1.-1./zcb,n); //cache this as laplcb[m][n]
	  complex <double> laplbb = 1./zbb * pow(1.-1./zbb,n); //cache this as laplbb[m][n]

	  UP[m] += lup[n] * laplup;
	  DO[m] += ldo[n] * lapldo;
	  ST[m] += lst[n] * laplst;
	  CH[m] += lch[n] * laplch;
	  BO[m] += lbo[n] * laplbo;
	  GL[m] += lgl[n] * laplgl;
	  UB[m] += lub[n] * laplub;
	  DB[m] += ldb[n] * lapldb;
	  SB[m] += lsb[n] * laplsb;
	  CB[m] += lcb[n] * laplcb;
	  BB[m] += lbb[n] * laplbb;
	  //cout << m << "  " << n << "  " << real(fli[n] * lapl) << "  " << z << "  " << lapl << endl;
	}
    }
  
}

void mellinpdf::ceres_pdf()
{
#ifdef USECERES
  for (int k = 0; k < mellinint::mdim; k++)
    {
      UP[k] = cerespdf::U.mellin(mellinint::Np[k]);
      DO[k] = cerespdf::D.mellin(mellinint::Np[k]);
      ST[k] = cerespdf::S.mellin(mellinint::Np[k]);
      CH[k] = cerespdf::C.mellin(mellinint::Np[k]);
      BO[k] = cerespdf::B.mellin(mellinint::Np[k]);
      GL[k] = cerespdf::G.mellin(mellinint::Np[k]);
      UB[k] = cerespdf::UB.mellin(mellinint::Np[k]);
      DB[k] = cerespdf::DB.mellin(mellinint::Np[k]);
      SB[k] = cerespdf::SB.mellin(mellinint::Np[k]);
      CB[k] = cerespdf::CB.mellin(mellinint::Np[k]);
      BB[k] = cerespdf::BB.mellin(mellinint::Np[k]);
    }
#endif
}
