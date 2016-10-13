#include "resint.h"
#include "settings.h"
//#include "interface.h"
#include "dyres_interface.h"
#include "resconst.h"
#include "mesq.h"
#include "hcoefficients.h"
#include "pdfevol.h"
#include "besselint.h"
#include "cubature.h"
#include <iostream>

#include "intde2.c"

#include <LHAPDF/LHAPDF.h>

double tiny = 1.0e-307;
double awinf[lenaw];
//double awfin[lenaw];

double resint::_qt;
double resint::_m;
double resint::_y;
double resint::_costh;

int resint::_mode;

double resint::muren;
double resint::mufac;
double resint::mures;
double resint::muren2;
double resint::mufac2;
double resint::mures2;

double resint::a;
complex <double> resint::loga;
double resint::rloga;
complex <double> resint::logmuf2q2;
complex <double> resint::logq2muf2;
complex <double> resint::logq2mur2;
double resint::rlogq2mur2;

double resint::alpqr;
double resint::alpqf;
double resint::aass;

double resint::alpqfac;
double resint::alpqren;
double resint::alpqres;

//using namespace resint;
//#pragma omp threadprivate(_qt,_m,_y,_costh,_mode,muren,mufac,mures,muren2,mufac2,mures2,a,loga,rloga,logmuf2q2,logq2muf2,logq2mur2,rlogq2mur2,alpqr,alpqf,aass,alpqfac,alpqren,alpqres)


int besselint_cubature(unsigned ndim, const double x[], void *data, unsigned ncomp, double f[])
{
  f[0] = besselint::bint(x[0]);
  return 0;
}
void resint::init()
{
  //Set values in the common block for the sudakov
  flag1_.flag1_ = opts.order;              //order for the sudakov
  iorder_.iord_ = opts.order - 1;          //order for LL/NLL running of alphas
  flagrealcomplex_.flagrealcomplex_ = 0;   // choose real axis (complex plane) integration of bstar (b) (always 0)
  modified_.imod_ = 1;                     // normal (imod=0) or modified (imod=1) sudakov      

  // g-parameter of the non perturbative form factor
  np_.g_ = opts.g_param;

  //a-parameter of the resummation scale, set it for the dynamic case
  if (opts.dynamicresscale)
    {
      a = 1./opts.kmures;
      a_param_.a_param_ = a;
      a_param_.b0p_ = resconst::b0*a;

      // precompute log of scales
      loga = log(a);
      rloga = log(a);
      rlogs_.rloga_ = rloga;
      //  clogs_.loga_ = loga; //complex loga is not used
    }
  
  // define points for quadratures integration
  intdeoini(lenaw, tiny, opts.bintaccuracy, awinf);
  //intdeini(lenaw, tiny, 1.0e-2, awfin);
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
  
  //move the scale setting to an independent namelist, so it can be called by all the calculations
  //Set factorization and renormalization scales (set dynamic scale here)
  if (opts.dynamicscale)
    {
      muren = m*opts.kmuren;
      mufac = m*opts.kmufac;
    }
  else
    {
      muren = opts.rmass*opts.kmuren;
      mufac = opts.rmass*opts.kmufac;
    }

  if (opts.dynamicresscale)
    mures = m*opts.kmures;
  else
    //for fixed resummation scale need to recompute a_param
    {
      mures = opts.rmass*opts.kmures;
      a = m/mures;
      loga = log(a);
      rloga = log(a);
      a_param_.a_param_ = a;
      a_param_.b0p_ = resconst::b0*a;
      rlogs_.rloga_ = rloga;
    }

  muren2=pow(muren,2);
  mufac2=pow(mufac,2);
  mures2=pow(mures,2);

  //alphas at various scales (alphas convention is alphas/4/pi)
  alpqren=LHAPDF::alphasPDF(muren)/4./M_PI;
  alpqfac=LHAPDF::alphasPDF(mufac)/4./M_PI;
  alpqres=LHAPDF::alphasPDF(mures)/4./M_PI;

  //alpqr is alphas at the renormalization scale
  if (opts_.approxpdf_ == 1)
    {
      int nloop = 3;
      double scale = muren;
      alpqr=dyalphas_mcfm_(scale,couple_.amz_,nloop)/4./M_PI;
    }
  else
    alpqr=alpqren;

  //alpqf is alphas at the resummation scale. It is used as starting scale for the PDF evolution
  if (opts_.approxpdf_ == 1)
    {
      int nloop = 3;
      double scale = mures;
      alpqf=dyalphas_mcfm_(scale,couple_.amz_,nloop)/4./M_PI;
    }
  else
    alpqf=alpqres;

  /***************************************************/
  //set values in the common blocks for the sudakov

  //aaas is alphas/pi at renormalization scale
  aass = alpqr*4.;
  aass_.aass_ = aass; 
  double q = m;
  double blim;
  if (q2 > 2.2*muren2) //avoid large values of b in f2(y) when q2>mur2 
    blim = a_param_.b0p_*(1./q)*exp(1./(2.*aass*resconst::beta0))*sqrt(muren2/q2);
  else
    blim = a_param_.b0p_*(1./q)*exp(1./(2.*aass*resconst::beta0)); // avoid Landau pole

  if (opts.blim > 0)
    {
      double lambdaqcd         = muren/(exp(1./(2.*aass*resconst::beta0))); //--> corresponds to a divergence in alphas
      double lambdasudakov     = mures/(exp(1./(2.*aass*resconst::beta0))); //--> correspond to a divergence in the Sudakov

      //should blim depend or not on a_param? In principle yes because PDFs are evolved to b0/bstar*a
      //      blim = resconst::b0/lambdasudakov
      blim = min(resconst::b0/lambdaqcd,resconst::b0/lambdasudakov);
      //cout << _m << "  " << resconst::b0/lambdaqcd << "  " << resconst::b0/lambdasudakov << "  " << a_param_.b0p_*(1./q)*exp(1./(2.*aass*resconst::beta0)) << endl;

      blim = blim/opts.blim;
    }

  //  blim = 100;
  //  blim = 4.0;

  blimit_.rblim_ = blim;
  blimit_.cblim_ = fcx(blim);

  //  cout << q << "  " << mur << "  " << blimit_.rblim_ << "  " << a_param_.b0p_/lambdaqcd << "  " << lambdaqcd <<"  " <<  a_param_.b0p_/blim << endl;

  //fill in common block for scales (used in the sudakov)
  scaleh_.qt_ = _qt;
  scaleh_.qt2_ = _qt*_qt;
  scaleh_.xtau_ = 0.;
  scaleh_.q_ = _m;
  scaleh_.q2_ = _m*_m;
  scaleh_.shad_ = pow(opts.sroot,2);
  scaleh_.sroot_ = opts.sroot;
  scaleh_.mur_ = muren;
  scaleh_.mur2_ = muren2;
  scaleh_.muf_ = mufac;
  scaleh_.muf2_ = mufac2;
  /***************************************************/

  //******************************************
  //     Initialise the born-level matrix elements mesqij
  //     ---> mesqij depends on mass and costh
  //     when costh is integrated and costh moments are provided, mesqij_expy
  //     are convoluted with the rapidity depent expression
  //      and they depend also on rapidity

  mesq::allocate();
  mesq::setmesq_expy(mode, m, costh, y);
  //******************************************

  //*****************************************
  //precompute scales: there is mass dependence in logq2muf2 logq2mur2, only with fixed factorization and renormalization scales
  logmuf2q2=log(mufac2/q2);
  logq2muf2=log(q2/mufac2);
  logq2mur2=log(q2/muren2);
  rlogq2mur2=log(q2/muren2);

  rlogs_.rlogq2mur2_ = real(logq2mur2);
  
  //no need to recompute coefficients with variable scales
  //if (!dynamicscale || fixedresummationscale)
  hcoefficients::allocate();
  hcoefficients::calc(aass,logmuf2q2,logq2muf2,logq2mur2,loga);
  //*****************************************


  //*****************************************
  //dependence on qt, m, also y, and costh unless integrated
  //perform b integration (int_0^inf db)
  double res = bintegral(qt);
  //  cout << "phase space point in resumm" << "  " << _qt << "  " <<  _y << "  " <<  _m << "  " << _costh << "  " << res << endl;
  //*****************************************

  // Normalization
  res *= qt/2./pow(opts.sroot,2);

  //free allocated Local Thread Storage (LTS) memory
  hcoefficients::free();
  mesq::free();

  return res;
}

double resint::bintegral(double qt)
{
  // Compute integral using double exponential quadratures for oscillatory functions (intde2.c)

  double res, err;
  intdeo(besselint::bint, 0.0, qt, awinf, &res, &err);
  //  cout << "dequad result of inverse bessel transform, pt=" << _qt << " m=" << _m << " y=" << _y << " : " << setprecision(16) << res << " +- " << err/res*100 << "%" << endl;
  //  cout.precision(6); cout.unsetf(ios_base::floatfield);

  /*
  //split the integral above and below the b mass
  //scale b0/b without a_param
  double bstar_mb = resconst::b0/LHAPDF::getThreshold(5);
  
  //convert bstar to b
  double b_mb = bstar_mb / sqrt(1-(bstar_mb*bstar_mb)/(blimit_.rblim_*blimit_.rblim_));


  double res1, err1;
  intdeo(besselint::bint, b_mb, qt, awinf, &res1, &err1);


  double res2, err2;
  //dequad seems not to be appropriate for this integral 
  //  intde(besselint::bint, 0.0, b_mb, awfin, &res2, &err2);


  //try simple pcubature
  const int ndim = 1;     //dimensions of the integral
  const int ncomp = 1;  //components of the integrand
  void *userdata = NULL;
  double integral[1];
  double error[1];
  const int eval = 0;
  const double epsrel = 1e-2;//opts.bintaccuracy;
  const double epsabs = 0.;
  double xmin[1] = {0};
  double xmax[2] = {b_mb};

  pcubature(ncomp, besselint_cubature, userdata, 
	    ndim, xmin, xmax, 
	    eval, epsabs, epsrel, ERROR_INDIVIDUAL, integral, error);
  res2 = integral[0];
  err2 = error[0];

  
  cout << "split: dequad result of inverse bessel transform  " << setprecision(16) << res1+res2 << " +- " << res1 << "  " << err1/res1*100 << "% "  << res2 << "  " << err2/res2*100 << "%" << endl;
  cout.precision(6); cout.unsetf(ios_base::floatfield);
  res = res1+res2;
  cout << endl;

  //plot the b-integrand
  cout << "{" << endl;
  cout << "TGraph *g = new TGraph();" << endl;
  for (int i = 0; i < 1000; i++)
    {
      //      double b = b_mb*0.95+i*(b_mb*1.05 -b_mb*0.95)/100.;
      double b = 0.+i*(b_mb*50 -0)/1000.;
  //      cout << b << "  " << b_mb << "  " << besselint::bint(b) << endl;;
      cout << "g->SetPoint(g->GetN(), " << b << ", " << besselint::bint(b) << ");" << endl;
    }
  cout << "g->Draw();" << endl;
  cout << "}" << endl;
  */

  return res;
}
