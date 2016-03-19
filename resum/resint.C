#include "resint.h"
#include "settings.h"
#include "interface.h"
#include "resconst.h"
#include "mesq.h"
#include "hcoefficients.h"
#include "besselint.h"
#include <iostream>

#include "intde2.c"

#include <LHAPDF/LHAPDF.h>

double tiny = 1.0e-307;
double aw[lenaw];

double resint::_qt;
double resint::_m;
double resint::_y;
double resint::_costh;

int resint::_mode;

double resint::mur;
double resint::muf;
double resint::mur2;
double resint::muf2;

complex <double> resint::loga;
double resint::rloga;
complex <double> resint::logmuf2q2;
complex <double> resint::logq2muf2;
complex <double> resint::logq2mur2;
double resint::rlogq2mur2;

double resint::alpqr;
double resint::alpqf;
double resint::aass;

//void intdeoini(int lenaw, double tiny, double eps, double *aw);
//void intdeo(double (*f)(double), double a, double omega, double *aw, double *i, double *err);

void resint::init()
{
  // precompute log of scales
  loga = log(opts.a_param);
  rloga = log(opts.a_param);
      
  //Set values in the common block for the sudakov
  flag1_.flag1_ = opts.order;              //order for the sudakov
  flagrealcomplex_.flagrealcomplex_ = 0;   // choose real axis (complex plane) integration of bstar (b) (always 0)
  modified_.imod_ = 1;                     // normal (imod=0) or modified (imod=1) sudakov      
  a_param_.a_param_ = opts.a_param;        //dynamic resummation scale
  a_param_.b0p_ = resconst::b0*opts.a_param;
  np_.g_ = opts.g_param;                   // non perturbative parameter
  rlogs_.rloga_ = rloga;
  
  // define points for quadratures integration
  intdeoini(lenaw, tiny, 1.0e-2, aw);
}

double resint::rint(double costh, double m, double qt, double y, int mode)
{
  //point in phase space
  _qt = qt;
  _y = y;
  _m = m;
  _costh = costh;

  //mode 0: differential mode 1: integrated in costh mode 2: integrated in costh and y  
  _mode = mode;

  // Kinematical limit
  double q2=m*m;
  double x = q2/pow(opts.sroot,2);
  double ax = log(x);
  double ylim=-0.5*ax;
  if (abs(y) > ylim)
    return 0;
  
  //Set factorization and renormalization scales (set dynamic scale here)
  //if (opts.dynamicscale) else
  mur = opts.mur;
  muf = opts.muf;

  mur2=pow(mur,2);
  muf2=pow(muf,2);

  //alpqr is alphas at the renormalization scale
  if (opts_.approxpdf_ == 1)
    {
      int nloop = 3;
      double scale = sqrt(mur2);
      alpqr=dyalphas_mcfm_(scale,couple_.amz_,nloop)/4./M_PI;
    }
  else
    alpqr=LHAPDF::alphasPDF(mur)/4./M_PI;

  //squared of resummation scale (used for the evolution)       
  double q2s=q2/pow(opts.a_param,2);
  //...  COUPLING CONSTANTS AT INPUT SCALE = ALPHAS/4/PI
  //     ALPQF = ALPHA AT RESUMMATION SCALE (used to start evolution)
  if (opts_.approxpdf_ == 1)
    {
      int nloop = 3;
      double scale = sqrt(q2s);
      alpqf=dyalphas_mcfm_(scale,couple_.amz_,nloop)/4./M_PI;
    }
  else
    alpqf=LHAPDF::alphasPDF(sqrt(q2s))/4./M_PI;
         
  //aaas is alphas/pi
  aass = alpqr*4.;

  /***************************************************/
  //set values in the common block for the sudakov
  aass_.aass_ = aass; 
  double q = m;
  blimit_.rblim_=a_param_.b0p_*(1./q)*exp(1./(2.*aass*resconst::beta0)); // avoid Landau pole     
  if (q2 > 2.2*mur2) //AVOID LARGE VALE OF b IN f2(y) WHEN q2>mur2 
    blimit_.rblim_ = a_param_.b0p_*(1./q)*exp(1./(2.*aass*resconst::beta0))*sqrt(mur2/q2);

  blimit_.cblim_= fcx(a_param_.b0p_*(1./q)*exp(1./(2.*aass*resconst::beta0)));
  if(q2 > 2.2*mur2)
    blimit_.cblim_ = fcx(a_param_.b0p_*(1./q)*exp(1./(2.*aass*resconst::beta0))*sqrt(mur2/q2));

  //common block for scales (in the sudakov)
  scaleh_.qt_ = _qt;
  scaleh_.qt2_ = _qt*_qt;
  scaleh_.xtau_ = 0.;
  scaleh_.q_ = _m;
  scaleh_.q2_ = _m*_m;
  scaleh_.shad_ = pow(opts.sroot,2);
  scaleh_.sroot_ = opts.sroot;
  scaleh_.mur_ = mur;
  scaleh_.mur2_ = mur2;
  scaleh_.muf_ = muf;
  scaleh_.muf2_ = muf2;
  /***************************************************/

  //******************************************
  //     Initialise the born-level matrix elements mesqij
  //     ---> mesqij depends on mass and costh
  //     when costh is integrated and costh moments are provided, mesqij_expy
  //     are convoluted with the rapidity depent expression
  //      and they depend also on rapidity
  mesq::setmesq_expy(mode, m, costh, y);
  //******************************************

  //*****************************************
  //precompute scales: there is mass dependence in logq2muf2 logq2mur2, only with fixed factorization and renormalization scales
  logmuf2q2=log(muf2/q2);
  logq2muf2=log(q2/muf2);
  logq2mur2=log(q2/mur2);
  rlogq2mur2=log(q2/mur2);

  rlogs_.rlogq2mur2_ = real(logq2mur2);
  
  //no need to recompute coefficients with variable scales
  //if (!dynamicscale || fixedresummationscale)
  hcoefficients::calc(aass,logmuf2q2,logq2muf2,logq2mur2,loga);
  //*****************************************

  //*****************************************
  //dependence on qt, m, also y, and costh unless integrated
  //perform b integration (int_0^inf db)
  double res = bintegral(qt);
  // cout << "phase space point in resumm" << "  " << _qt << "  " <<  _y << "  " <<  _m << "  " << _costh << "  " << res << endl;
  //*****************************************

  // Normalization
  res *= qt/2./pow(opts.sroot,2);

  return res;
}

double resint::bintegral(double qt)
{
  // Compute integral using double exponential quadratures for oscillatory functions (intde2.c)
  double res, err;
  intdeo(besselint::bint, 0.0, qt, aw, &res, &err);
  //  cout << "dequad result of inverse bessel transform  " << res << "  " << err << endl;
  return res;
}
