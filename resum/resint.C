#include "resint.h"
#include "phasespace.h"
#include "settings.h"
#include "dyres_interface.h"
#include "resconst.h"
#include "mesq.h"
#include "omegaintegr.h"
#include "rapint.h"
#include "gaussrules.h"
#include "hcoefficients.h"
#include "hcoeff.h"
#include "hcoeff_check.h"
#include "expc.h"
#include "muf.h"
#include "gint.h"
#include "blim.h"
#include "pdfevol.h"
#include "mellinpdf.h"
#include "besselint.h"
#include "cubature.h"
#include "pdf.h"
#include "scales.h"
#include "sudakovff.h"
#include "besche_interface.h"
#include "hankel.h"
#include "minprescription_interface.h"
#include "specialfunctions_interface.h"
#include "588_interface.h"
#include "bequad.h"
#include "ccoeff.h"
#include "pegasus.h"
#include "pmom.h"
#include <iostream>

#include "intde2_c.h"

#include <LHAPDF/LHAPDF.h>

double tiny = 1.0e-307;
double awinf[lenaw];
double awdei[lenaw];
//double awfin[lenaw];

double resint::_qt;
double resint::_m;
double resint::_y;
double resint::_costh;

double resint::tau;
double resint::x1;
double resint::x2;

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

double resint::LQ;
double resint::LF;
double resint::LR;

double resint::alpqr;
double resint::alpqf;
double resint::aass;

double resint::alpqfac;
double resint::alpqren;
double resint::alpqres;

double resint::bc;

//using namespace resint;
//#pragma omp threadprivate(_qt,_m,_y,_costh,_mode,muren,mufac,mures,muren2,mufac2,mures2,a,loga,rloga,logmuf2q2,logq2muf2,logq2mur2,rlogq2mur2,alpqr,alpqf,aass,alpqfac,alpqren,alpqres)

extern "C"
{
  double besselint_besche_(double &x)
  {
    return real(besselint::bint(x));
  }

  double bint_mp_up_real_(double &x)
  {
    complex <double> jac(cos(M_PI/opts.phibr),sin(M_PI/opts.phibr));
    complex <double> b = resint::bc + jac*x;
    fcomplex arg = fcx(resint::_qt*b);
    return real(besselint::bint(b)*jac*cx(h1_(arg)));
  }

  double bint_mp_up_imag_(double &x)
  {
    complex <double> jac(cos(M_PI/opts.phibr),sin(M_PI/opts.phibr));
    complex <double> b = resint::bc + jac*x;
    fcomplex arg = fcx(resint::_qt*b);
    return imag(besselint::bint(b)*jac*cx(h1_(arg)));
  }

  double bint_mp_dn_real_(double &x)
  {
    complex <double> jac(cos(M_PI/opts.phibr),-sin(M_PI/opts.phibr));
    complex <double> b = resint::bc + jac*x;
    fcomplex arg = fcx(resint::_qt*b);
    return real(besselint::bint(b)*jac*cx(h2_(arg)));
  }

  double bint_mp_dn_imag_(double &x)
  {
    complex <double> jac(cos(M_PI/opts.phibr),-sin(M_PI/opts.phibr));
    complex <double> b = resint::bc + jac*x;
    fcomplex arg = fcx(resint::_qt*b);
    return imag(besselint::bint(b)*jac*cx(h2_(arg)));
  }

//  double besselint_mp_real_(double &x)
//  {
//    return besselint_mp_real_dequad(x);
//  }
}

int besselint_cubature(unsigned ndim, const double x[], void *data, unsigned ncomp, double f[])
{
  f[0] = real(besselint::bint(x[0]));
  return 0;
}

double besselint_dequad(double x)
{
  return real(besselint::bint(x));
}

double besselint_hankel(double x)
{
  double qtb = x*resint::_qt;
  return real(besselint::bint(x))/x/fort_besj0_(qtb);
}

double besselint_588(double &x)
{
  double qtb = x*resint::_qt;
  return real(besselint::bint(x))/fort_besj0_(qtb);
}


double besselint_mp_real_dequad(double x)
{
  /*
  complex <double> jacu(cos(M_PI/opts.phibr),sin(M_PI/opts.phibr));
  complex <double> bu = resint::bc + jacu*x;
  fcomplex argu = fcx(resint::_qt*bu);
  complex <double> jacd(cos(M_PI/opts.phibr),-sin(M_PI/opts.phibr));
  complex <double> bd = resint::bc + jacd*x;
  fcomplex argd = fcx(resint::_qt*bd);

  complex <double> h1, h2;

//  //call modified Hankel functions
//  h1 = cx(h1_(argu));
//  h2 = cx(h2_(argd));
//  //cout << bu << "  " << h1 << "  " << bd << "  " << h2 << endl;

  //call original Hankel functions
  int n = 1;
  int nm;    
  fcomplex chf1[2], chd1[2], chf2[2], chd2[2];
  ch12n_ (n, argu, nm, chf1, chd1, chf2, chd2);
  h1 = cx(chf1[0]);
  ch12n_ (n, argd, nm, chf1, chd1, chf2, chd2);
  h2 = cx(chf2[0]);
  //cout << endl;
  //cout << " b " << bu << "  " << bd << " bint " << besselint::bint(bu) << "  " << besselint::bint(bd) << " h " << h1 << "  " << h2 << " jac " << jacu << "  " << jacd << endl;
  
  //return (real(besselint::bint(bu)*jacu*h1) + real(besselint::bint(bd)*jacd*h2))/2.;
  */
  
  //Real axis integration (phib = 0)
  double bb = resint::bc + x;
  double qtb = bb*resint::_qt;
  complex <double> bint = besselint::bint(bb);
  double JN, YN;
  if (resint::_mode == 3 || resint::_mode == 4)
    {
      JN = fort_besj1_(qtb);
      YN = fort_besy1_(qtb);
    }
  else
    {
      JN = fort_besj0_(qtb);
      YN = fort_besy0_(qtb);
    }
  //double res = real(bint)*JN-imag(bint)*YN;
  double res = real(bint)*JN;
  //cout << " b " << bb << "  "  << " bint " << bint << " J0 " << J0 << "  " << Y0 << endl;

  //cout << " b " << bb << " h " << (real(besselint::bint(bu)*jacu*h1) + real(besselint::bint(bd)*jacd*h2))/2. << " j " << res << endl;
  
  return res;

}

double besselint_mp_complex_dequad(double x)
{
  /*
  complex <double> jacu(cos(M_PI/opts.phibr),sin(M_PI/opts.phibr));
  complex <double> bu = resint::bc + jacu*x;
  fcomplex qtbu = fcx(resint::_qt*bu);
  complex <double> jacd(cos(M_PI/opts.phibr),-sin(M_PI/opts.phibr));
  complex <double> bd = resint::bc + jacd*x;
  fcomplex qtbd = fcx(resint::_qt*bd);
  */


  /*********** Butterfly contour **************/
  //Straight line from 0 to bc
  if (x < resint::bc)
    {
      double bb = x;
      double qtb = bb*resint::_qt;
      complex <double> bint = besselint::bint(bb);
      double JN, YN;
      if (resint::_mode == 3 || resint::_mode == 4)
  	JN = fort_besj1_(qtb);
      else
  	JN = fort_besj0_(qtb);
      
      return real(bint)*JN;
    }
  
  
  //Move away from real axis from bc to blim::sudakov
  complex <double> jacu, jacd;
  complex <double> bu, bd;
  jacu = complex <double> (1.,tan(M_PI/opts.phibr));
  jacd = complex <double> (1.,-tan(M_PI/opts.phibr));
  bu = resint::bc + jacu*(x-resint::bc);
  bd = resint::bc + jacd*(x-resint::bc);
  /**********************************************/

  ///*********** Parallelogram contour **************/
  // //--> make this an option?
  // //Continue parallel to the real axis after the Landau pole at bmax
  // double bmax;
  // if (opts.modlog)
  //   bmax = resconst::b0/scales::res * sqrt(exp(1./(resint::aass*resconst::beta0))-1.);
  // else
  //   bmax = resconst::b0/scales::res * exp(1./(2.*resint::aass*resconst::beta0));
  // 
  // double pole = resint::bc + tan(M_PI/opts.phibr)*(bmax-resint::bc);
  // //double pole = 0.;//resint::bc + tan(M_PI/opts.phibr)*(blim::sudakov-resint::bc);
  // if (x > pole)
  //   {
  //     jacu = 1.;
  //     jacd = 1.;
  //     bu = bmax + (x-pole) + complex <double> (0., (bmax-resint::bc)*tan(M_PI/opts.phibr));
  //     bd = bmax + (x-pole) + complex <double> (0.,-(bmax-resint::bc)*tan(M_PI/opts.phibr));
  //   }
  ///**********************************************/

  /*********** Triangular contour **************/
  // //Go back to the real axis after the Landau pole at bmax
  // double bmax;
  // if (opts.modlog)
  //   bmax = resconst::b0/scales::res * sqrt(exp(1./(resint::aass*resconst::beta0))-1.);
  // else
  //   bmax = resconst::b0/scales::res * exp(1./(2.*resint::aass*resconst::beta0));
  // 
  // if (x > bmax)
  //   {
  //     jacu = complex <double> (1.,-tan(M_PI/opts.phibr));
  //     jacd = complex <double> (1.,+tan(M_PI/opts.phibr));
  //     bu = x + complex <double> (0., (bmax-resint::bc)*tan(M_PI/opts.phibr)) + jacu*(x-bmax);
  //     bd = x + complex <double> (0.,-(bmax-resint::bc)*tan(M_PI/opts.phibr)) + jacd*(x-bmax);
  //   }
  // 
  // if (x > bmax+(bmax-resint::bc))
  //   {
  //     double bb = x;
  //     double qtb = bb*resint::_qt;
  //     complex <double> bint = besselint::bint(bb);
  //     double JN, YN;
  //     if (resint::_mode == 3 || resint::_mode == 4)
  // 	JN = fort_besj1_(qtb);
  //     else
  // 	JN = fort_besj0_(qtb);
  //     
  //     return real(bint)*JN;
  //   }
  // /**********************************************/

  /*********** Rectangular contour **************/
  // //Straight line from 0 to bc
  // if (x < resint::bc)
  //   {
  //     double bb = x;
  //     double qtb = bb*resint::_qt;
  //     complex <double> bint = besselint::bint(bb);
  //     double JN, YN;
  //     if (resint::_mode == 3 || resint::_mode == 4)
  // 	JN = fort_besj1_(qtb);
  //     else
  // 	JN = fort_besj0_(qtb);
  //     
  //     return real(bint)*JN;
  //   }
  // 
  // 
  // //Move away from real axis at b = bc
  // double bmax;
  // if (opts.modlog)
  //   bmax = resconst::b0/scales::res * sqrt(exp(1./(resint::aass*resconst::beta0))-1.);
  // else
  //   bmax = resconst::b0/scales::res * exp(1./(2.*resint::aass*resconst::beta0));
  // double bim = bmax/10.;
  // complex <double> jacu, jacd;
  // complex <double> bu, bd;
  // //Straight line from (bc,0) to (bc,bim)
  // if (x < resint::bc+bim)
  //   {
  //     jacu = complex <double> (0.,1.);
  //     jacd = complex <double> (0.,-1.);
  //     bu = resint::bc + jacu*x;
  //     bd = resint::bc + jacd*x;
  //   }
  // if (x >= resint::bc+bim)
  //   {
  //     jacu = 1.;
  //     jacd = 1.;
  //     bu = complex <double> (x-bim,bim);
  //     bd = complex <double> (x-bim,-bim);
  //   }
  // /**********************************************/

  // double bmax;
  // if (opts.modlog)
  //   bmax = resconst::b0/scales::res * sqrt(exp(1./(resint::aass*resconst::beta0))-1.);
  // else
  //   bmax = resconst::b0/scales::res * exp(1./(2.*resint::aass*resconst::beta0));
  // double bim = bmax/2.;
  // complex <double> jacu, jacd;
  // complex <double> bu, bd;
  // jacu = 1.;
  // jacd = 1.;
  // bu = complex <double> (x,bim);
  // bd = complex <double> (x,-bim);

  
  fcomplex qtbu = fcx(resint::_qt*bu);
  fcomplex qtbd = fcx(resint::_qt*bd);
  
  complex <double> h1, h2;

//  //call modified Hankel functions
//  h1 = cx(h1_(qtbu));
//  h2 = cx(h2_(qtbd));
//  //cout << bu << "  " << h1 << "  " << bd << "  " << h2 << endl;

  //call original Hankel functions
  // --> Here should use Hankel with J1 and Y1 for qt-integrated mode = 3!!!
  int n = 1;
  int nn;
  if (resint::_mode == 3 || resint::_mode == 4)
    nn = 1;
  else
    nn = 0;
  int nm;    
  //fcomplex chf1[n+1], chd1[n+1], chf2[n+1], chd2[n+1];
  fcomplex chf1[2], chd1[2], chf2[2], chd2[2];
  ch12n_ (n, qtbu, nm, chf1, chd1, chf2, chd2);
  h1 = cx(chf1[nn]);
  //ch12n_ (n, qtbd, nm, chf1, chd1, chf2, chd2);
  //h2 = cx(chf2[nn]);
  //cout << endl;
  //cout << " b " << bu << "  " << bd << " bint " << besselint::bint(bu) << "  " << besselint::bint(bd) << " h " << h1 << "  " << h2 << " jac " << jacu << "  " << jacd << endl;

  //Compute Hankel functions from fortran bessel functions (does not work, because they accept only real values...)
  //h1 = complex <double> (fort_besj0_(qtbu), fort_besy0_(qtbu));
  //h2 = complex <double> (fort_besj0_(qtbd), -fort_besy0_(qtbd));
  //h2 = conj(h1);
  
  complex <double> bint_up = besselint::bint(bu);
  //complex <double> bint_dn = besselint::bint(bd);
  //complex <double> bint_dn = conj(besselint::bint(bu));

  //cout << bu << " up-conj(dn) " << bint_up-conj(bint_dn) << endl;
  
  //return (real(bint_up*jacu*h1) + real(bint_dn*jacd*h2))/2.;
  return real(bint_up*jacu*h1); //<--

  //test pure bessel J
  double JN;
  double qtb = real(resint::_qt*bu);
  if (resint::_mode == 3 || resint::_mode == 4)
    JN = fort_besj1_(qtb);
  else
    JN = fort_besj0_(qtb);
  
  return real(bint_up)*JN;
  
  
}

double frealu_dequad(double x)
{
  complex <double> jacu(cos(M_PI/opts.phibr),sin(M_PI/opts.phibr));
  complex <double> bu = resint::bc + jacu*x;
  fcomplex argu = fcx(resint::_qt*bu);

  complex <double> h1;

  //call modified Hankel functions
  /*
  h1 = cx(h1_(argu));
  //cout << bu << "  " << h1 << endl;
  */

  //call original Hankel functions
  int n = 1;
  int nm;    
  fcomplex chf1[2], chd1[2], chf2[2], chd2[2];
  ch12n_ (n, argu, nm, chf1, chd1, chf2, chd2);
  h1 = cx(chf1[0]);
  
  return real(besselint::bint(bu)*jacu*h1);
}

double freald_dequad(double x)
{
  complex <double> jacd(cos(M_PI/opts.phibr),-sin(M_PI/opts.phibr));
  complex <double> bd = resint::bc + jacd*x;
  fcomplex argd = fcx(resint::_qt*bd);

  complex <double> h2;

  //call modified Hankel functions
  /*
  h2 = cx(h2_(argd));
  //cout << bd << "  " << h2 << endl;
  */

  //call original Hankel functions
  int n = 1;
  int nm;    
  fcomplex chf1[2], chd1[2], chf2[2], chd2[2];
  ch12n_ (n, argd, nm, chf1, chd1, chf2, chd2);
  h2 = cx(chf2[0]);
  
  return real(besselint::bint(bd)*jacd*h2);
}

void resint::init()
{
  //Set values in the common block for the sudakov
  flag1_.flag1_ = opts.order;              //order for the sudakov
  iorder_.iord_ = opts.order - 1;          //order for LL/NLL running of alphas
  //modified_.imod_ = 1;                     // normal (imod=0) or modified (imod=1) sudakov
  modified_.imod_ = opts.modlog;           // canonical (imod=0) or modified (imod=1) logarithms in the Sudakov

  // choose real axis (complex plane) integration of bstar (b) (always 0)
  if (opts.bprescription == 2 || opts.bprescription == 3)
    flagrealcomplex_.flagrealcomplex_ = 1;
  else
    flagrealcomplex_.flagrealcomplex_ = 0;

  //v parameter of the modified Hankel functions
  v_.v_=1.;

  //g-parameter of the non perturbative form factor
  np_.g_ = opts.g1;

  //a-parameter of the resummation scale, set it for the dynamic case
  if (opts.fmures > 0)
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
  //intdeoini(lenaw, tiny, 1.0e-6, awinf);
  //intdeiini(lenaw, tiny, 1.0e-6, awdei);
  //intdeini(lenaw, tiny, 1.0e-2, awfin);
}

void resint::rint(double costh, double m, double qt, double y, int mode, double f[2])
{
  f[0] = 0.;
  f[1] = 0.;

  /*
  //limit the integration to qt = m
  double jac = 1;
  double qtp = qt;
  if (qtp >= m*0.999)
    return 0.;
  else
    qt = qtp/sqrt(1-pow(qtp/m,2));
  jac = pow(m/sqrt(m*m-qtp*qtp),3);
  */
  
  //point in phase space (should use the phasespace namespace instead of storing them in resint)
  _qt = qt;
  _y = y;
  _m = m;
  _costh = costh;

  //compute x1 and x2 for the NP ff
  double exppy = exp(y);
  double expmy = 1./exppy;
  tau = m/opts.sroot;
      
  x1 = tau*exppy;
  x2 = tau*expmy;

  //mode 0: differential
  //mode 1: integrated in costh
  //mode 2: integrated in costh and y
  //mode 3: integrated in costh, y, and pt
  //mode 4: integrated in costh and pt
  _mode = mode;

  // Kinematical limit in rapidity
  double q2 = m*m;
  double x = q2/pow(opts.sroot,2);
  double ax = log(x);
  double ylim=-0.5*ax;
  if (abs(y) > ylim)
    return;
  
  //Set QCD scales
  scales::set(m);
  muren = scales::ren;
  mufac = scales::fac;
  mures = scales::res;
  //cout << "m " << m << " scales " << scales::ren << "  " << scales::fac << "  " << scales::res << endl;
  
  //for fixed resummation scale need to recompute a_param --> always recompute
  //if (opts.fmures == 0)
  //{
      a = m/mures;
      loga = log(a);
      rloga = log(a);
      a_param_.a_param_ = a;
      a_param_.b0p_ = resconst::b0*a;
      rlogs_.rloga_ = rloga;
      //}

  muren2=pow(muren,2);
  mufac2=pow(mufac,2);
  mures2=pow(mures,2);

  LQ = log(q2/mures2);
  LF = log(q2/mufac2);
  LR = log(q2/muren2);

  if (opts.melup == 2)
    {
      // -> put all this into mellinint::allocate()
      mellinint::allocate();
      ccoeff::allocate();
      //pegasus::allocate();
      pmom::allocate();
      mellinint::update();
    }
  //  else
  //    {
  pdfevol::allocate();
  pdfevol::update();
      //update PDFs in Mellin space at the starting scale, if the factorisation scale is proportional to mll
      //  if (opts.fmufac > 0)
      //    {
      //Truncate moments from xmin to 1, with xmin = m/sqrt(s) e^-y0 (currently works only at LL)
      //if (opts.mellin1d)
      //	{
      //  double xmin = phasespace::m/opts.sroot*exp(phasespace::ymin);
      //  delete[] mellinpdf::t;
      //  delete[] mellinpdf::fac;
      //  delete[] mellinpdf::kern;
      //  mellinpdf::init(xmin);
      //  //ccoeff::truncate();
      //}
      
      //pdfevol::allocate();
      //pdfevol::update();
      //    }
      //    }
  
  if (!opts.mellin1d && (mode == 2 || mode == 3))
    {
      double ylim = 0.5*log(pow(opts.sroot,2)/phasespace::m2);
      double ymn = min(max(-ylim, phasespace::ymin),ylim);
      double ymx = max(min(ylim, phasespace::ymax),-ylim);
      rapint::allocate();
      if (opts.makecuts)
	{
	  //there is a potential issue here, when lepton cuts are applied
	  //the rapidity dependent exponential are cached assuming integration between ymin and ymax
	  //for consistency, has to keep the integration between ymin and ymax
	  //rapint::integrate(phasespace::ymin,phasespace::ymax,phasespace::m);
	  rapint::numint(ymn,ymx,phasespace::m);
	}
      else
	rapint::integrate(ymn,ymx,phasespace::m);
    }
  
  //alphas at various scales (alphas convention is alphas/4/pi)
  //alpqren=LHAPDF::alphasPDF(muren)/4./M_PI;
  //alpqfac=LHAPDF::alphasPDF(mufac)/4./M_PI;
  //alpqres=LHAPDF::alphasPDF(mures)/4./M_PI;

  if (opts.alphaslha)
    {
      alpqren=pdf::alphas(muren)/4./M_PI;
      alpqfac=pdf::alphas(mufac)/4./M_PI;
      alpqres=pdf::alphas(mures)/4./M_PI;
    }
  else
    {
      alpqren=pdf::rgktalphas(muren)/4./M_PI;
      alpqfac=pdf::rgktalphas(mufac)/4./M_PI;
      alpqres=pdf::rgktalphas(mures)/4./M_PI;
    }

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
  blim::set();

//
//  double blim;
//  //dynamic blim
//  if (opts.blim < 0)
//    {
//      //if (opts.order == 2 && q2 > 2.2*muren2) //avoid large values of b in f2(y) when q2>mur2
//      if (q2 > 2.2*muren2) //avoid large values of b in f2(y) when q2>mur2 
//	blim = a_param_.b0p_*(1./q)*exp(1./(2.*aass*resconst::beta0))*sqrt(muren2/q2);
//      else
//	blim = a_param_.b0p_*(1./q)*exp(1./(2.*aass*resconst::beta0)); // avoid Landau pole
//      
//      /*
//      double lambdaqcd         = muren/(exp(1./(2.*aass*resconst::beta0))); //--> corresponds to a divergence in alphas
//      double lambdasudakov     = mures/(exp(1./(2.*aass*resconst::beta0))); //--> correspond to a divergence in the Sudakov
//
//      //should blim depend or not on a_param? In principle yes because PDFs are evolved to b0/bstar*a
//      blim = min(resconst::b0/lambdaqcd,resconst::b0/lambdasudakov);
//      //cout << _m << "  " << resconst::b0/lambdaqcd << "  " << resconst::b0/lambdasudakov << "  " << a_param_.b0p_*(1./q)*exp(1./(2.*aass*resconst::beta0)) << endl;
//      */
//
//      blim = blim/(opts.blim*(-1));
//    }
//
//  //fixed blim
//  if (opts.blim > 0)
//    blim = opts.blim;
//
//  //  blim = 100;
//  //  blim = 4.0;
//  //  blim = 1.5;
//
//  blimit_.rblim_ = blim;
//  blimit_.cblim_ = fcx(blim);
//  //  cout << q << "  " << mur << "  " << blimit_.rblim_ << "  " << a_param_.b0p_/lambdaqcd << "  " << lambdaqcd <<"  " <<  a_param_.b0p_/blim << endl;


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

  if (!(opts.order == 0 && opts.xspace))
    {
      mesq::allocate();
      mesq::setmesq_expy(mode, m, costh, y);
    }
  else //No Mellin transform case
    {
      mesq::setpropagators(m);
      
      if (mode == 0) //costh differential
	mesq::setmesq(1., costh, pow(costh,2));
      else //costh integrated
	{
	  double cthmom0, cthmom1, cthmom2;
	  omegaintegr::cthmoments(cthmom0, cthmom1, cthmom2);
	  mesq::setmesq(cthmom0, cthmom1, cthmom2);
	}
    }
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
  
  expc::allocate();
  if (opts.numsud || opts.numexpc || opts.order >= 4)
    gint::allocate();

  hcoeff::allocate();
  hcoeff::reset();
  hcoeff::calc();

  muf::allocate();
  muf::reset();

  /*
  if (opts.mellin1d)
    {

      //hcoeff_check::allocate();
      //hcoeff_check::reset();
      //hcoeff_check::calc();
    }
  else
    {
      //hcoefficients::allocate();
      //hcoefficients::reset();
      //hcoefficients::calc(aass,logmuf2q2,logq2muf2,logq2mur2,loga);
    }
  */
  
  pdfevol::allocate_fx();
  //*****************************************

  double res = calc(mode);
//        double res;
//        if (mode == 3)
//          {
//            /*
//      	Since qt appears only in the Bessel function J_0(bqt),
//      	the integration in qt could be factorised by using
//      	the following result: int_a^b J0(x) dx = Lambda0(b) - Lambda0(a)
//      	where Lambda0(x) = xJ0(x) + pi x/2 [J1(x)H0(x) - J0(x)H1(x)], with Hn Struve functions.
//      	!!! Things are actually much simpler because there is a qt factor (res *= qt/2./pow(opts.sroot,2);) !!!
//      	int_a^b x J0(x) dx = b J1(b) - a J1(a)
//      	See http://fisica.ciens.ucv.ve/~svincenz/TISPISGIMR.pdf and http://web.eah-jena.de/~rsh/Forschung/Stoer/besint.pdf (1.1.1)
//      	!!! The problem with this is the switching function, which depends on qt !!!
//      	--> The approach would be valid only for qt < m*k
//      	Now the integration call will automatically switch to mode = 2 for bins where qtmax > mmin*k
//            */
//      
//            //double qtmn = max(opts.qtcutoff,phasespace::qtmin);
//            double qtmn = phasespace::qtmin;
//      
//            //double kinqtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+phasespace::m2,2)/(4*pow(opts.sroot,2))-phasespace::m2)); //introduced max to avoid neqative argument of sqrt
//            double kinqtlim = 1e10;
//            double switchqtlim = opts.damp?phasespace::m*opts.dampk:1e10; //here use the value of qt where the switching start
//            double qtlim = min(kinqtlim, switchqtlim);
//            double qtmx = min(qtlim, phasespace::qtmax);
//      
//            if (qtmn >= qtmx || qtmx < opts.qtcutoff)
//      	res = 0;
//            else
//      	{
//      //	  double res2, err2;
//      //	  _qt = qtmx;
//      //	  intdeo(besselint_dequad, 0.0, qtmx, awinf, &res2, &err2);
//      //
//      //	  double res1, err1;
//      //	  _qt = qtmn;
//      //	  if (qtmn < opts.qtcutoff)
//      //	    res1 =  0.;
//      //	  else
//      //	    intdeo(besselint_dequad, 0.0, qtmn, awinf, &res1, &err1);
//      
//      	  double res1, res2;
//      	  _qt = qtmx;
//      	  res2 = bintegral(qtmx);
//      	  _qt = qtmn;
//      	  if (qtmn < opts.qtcutoff)
//      	    res1 =  0.;
//      	  else
//      	    res1 = bintegral(qtmn);
//      
//      	  res = qtmx*res2 - qtmn*res1;
//      
//      	  // Normalization
//      	  res *= 1./2./pow(opts.sroot,2);
//      	  //cout << " res " << res << endl;
//      	}
//          }
//        else
//          {
//      
//      //      //if (opts.fast && opts.xspace) {
//      //      //Fast function
//      //      res = exp(-qt);
//      //
//      //      // Normalization
//      //      double shad = pow(opts.sroot,2);
//      //      res *= qt/shad;
//      //      //}
//      //      //else
//      //      //{
//            
//            //*****************************************
//            //dependence on qt, m, also y, and costh unless integrated
//            //perform b integration (int_0^inf db)
//            res = bintegral(qt);
//            //  cout << "phase space point in resumm" << "  " << _qt << "  " <<  _y << "  " <<  _m << "  " << _costh << "  " << res << endl;
//            //*****************************************
//      
//            // Normalization
//            res *= qt/2./pow(opts.sroot,2);
//          }
//      
//      
//        //Do not need the Mellin transform for the LL case, because the HN coefficient is 1 --> Use PDFs in x space
//        if (opts.order == 0 && opts.xspace && (opts.evolmode == 0 || opts.evolmode == 1)) //if (opts.order == 0 && opts.xspace && (opts.evolmode == 1))
//          if (mode < 2) //rapidity differential
//            res *= mesq::loxs(x1, x2, mufac);
//          else           //rapidity integrated
//            res *= mesq::loxs(tau, mufac);      
//      
//        /*
//        if (opts.order == 0 && opts.xspace && (opts.evolmode == 0 || opts.evolmode == 1))
//          {
//            double fun = 0.;
//            if (mode < 2) //rapidity differential
//      	{
//      	  //PDFs
//      	  double fx1[2*MAXNF+1],fx2[2*MAXNF+1];
//      	  fdist_(opts.ih1,x1,mufac,fx1);
//      	  fdist_(opts.ih2,x2,mufac,fx2);
//      	  
//      	  for (int sp = 0; sp < mesq::totpch; sp++)
//      	    fun += fx1[mesq::pid1[sp]]*fx2[mesq::pid2[sp]]*real(mesq::mesqij[sp]);
//      	}
//            else //rapidity integration
//      	{
//      	  for (int i=0; i < opts.yintervals; i++)
//      	    {
//      	      double ya = phasespace::ymin+(phasespace::ymax-phasespace::ymin)*i/opts.yintervals;
//      	      double yb = phasespace::ymin+(phasespace::ymax-phasespace::ymin)*(i+1)/opts.yintervals;
//      	      double xc = 0.5*(ya+yb);
//      	      double xm = 0.5*(yb-ya);
//      	      for (int j=0; j < opts.yrule; j++)
//      		{
//      		  double y = xc+xm*gr::xxx[opts.yrule-1][j];
//      		  double exppy = exp(y);
//      		  double expmy = 1./exppy;
//      		  x1 = tau*exppy;
//      		  x2 = tau*expmy;
//      		      
//      		  //PDFs
//      		  double fx1[2*MAXNF+1],fx2[2*MAXNF+1];
//      		  fdist_(opts.ih1,x1,mufac,fx1);
//      		  fdist_(opts.ih2,x2,mufac,fx2);
//      		      
//      		  for (int sp = 0; sp < mesq::totpch; sp++)
//      		    fun += fx1[mesq::pid1[sp]]*fx2[mesq::pid2[sp]]*real(mesq::mesqij[sp])*gr::www[opts.yrule-1][j]*xm;;
//      		}
//      	    }
//      	}
//            //cout << "resint " << fun << endl;
//            res *= fun;
//          }
//        */

  double res_m = 0.;
  if (opts.helicity >= 0)
    if (!(opts.order == 0 && opts.xspace))
      {
	if (!opts.mellin1d && (mode == 2 || mode == 3))
	  {
	    double ylim = 0.5*log(pow(opts.sroot,2)/phasespace::m2);
	    double ymn = min(max(-ylim, phasespace::ymin),ylim);
	    double ymx = max(min(ylim, phasespace::ymax),-ylim);
	    rapint::integrate(ymn,ymx,phasespace::m, opts.helicity);
	    mesq::setmesq_expy(mode, m, costh, y);
	    res_m = calc(mode);
	  }
	else
	  {
	    //hcoeff::reset();
	    //hcoeff::calc(aass,logmuf2q2,logq2muf2,logq2mur2,loga);
	    
	    mesq::setmesq_expy(mode, m, costh, y, opts.helicity);
	    res_m = calc(mode);
	  }
      }

  //free allocated Local Thread Storage (LTS) memory
  expc::free();
  if (opts.numsud || opts.numexpc || opts.order >= 4)
    gint::free();
  if (opts.mellin1d)
    {
      hcoeff::free();
      muf::free();
      //hcoeff_check::free();
    }
  else
    {
      //hcoefficients::free();
      hcoeff::free();
    }

  if (!(opts.order == 0 && opts.xspace))
    mesq::free();

  if (opts.melup == 2)
    {
      mellinint::free();
      ccoeff::free();
      //pegasus::free();
      pmom::free();
    }

  
  //if (opts.melup != 2)
  //  if (opts.fmufac > 0)
  pdfevol::free();
  
  pdfevol::free_fx();
  //res *= jac;//jacobian for the change of variable qt=qtp/sqrt(1-qtp^2/m^2)

  f[0] = res;
  f[1] = res_m;

  //f[0] = 1;
  //f[1] = res;
  
  //return res;
}

double resint::calc(int mode)
{
  double res;
  if (mode == 3 || mode == 4)
    {
      /*
	Since qt appears only in the Bessel function J_0(bqt),
	the integration in qt could be factorised by using
	the following result: int_a^b J0(x) dx = Lambda0(b) - Lambda0(a)
	where Lambda0(x) = xJ0(x) + pi x/2 [J1(x)H0(x) - J0(x)H1(x)], with Hn Struve functions.
	!!! Things are actually much simpler because there is a qt factor (res *= qt/2./pow(opts.sroot,2);) !!!
	int_a^b x J0(x) dx = b J1(b) - a J1(a)
	See http://fisica.ciens.ucv.ve/~svincenz/TISPISGIMR.pdf and http://web.eah-jena.de/~rsh/Forschung/Stoer/besint.pdf (1.1.1)
	!!! The problem with this is the switching function, which depends on qt !!!
	--> The approach would be valid only for qt < m*k
	Now the integration call will automatically switch to mode = 2 for bins where qtmax > mmin*k
      */

      //double qtmn = max(opts.qtcutoff,phasespace::qtmin);
      double qtmn = phasespace::qtmin;

      //double kinqtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+phasespace::m2,2)/(4*pow(opts.sroot,2))-phasespace::m2)); //introduced max to avoid neqative argument of sqrt
      double kinqtlim = 1e10;
      double switchqtlim = opts.damp?phasespace::m*opts.dampk:1e10; //here use the value of qt where the switching start
      double qtlim = min(kinqtlim, switchqtlim);
      double qtmx = min(qtlim, phasespace::qtmax);

      double ylim = 0.5*log(pow(opts.sroot,2)/phasespace::m2);
      double ymn = min(max(-ylim, phasespace::ymin),ylim);
      double ymx = max(min(ylim, phasespace::ymax),-ylim);
      
      if (qtmn >= qtmx || qtmx < opts.qtcutoff)
	res = 0;
      else
	{
//	  double res2, err2;
//	  _qt = qtmx;
//	  intdeo(besselint_dequad, 0.0, qtmx, awinf, &res2, &err2);
//
//	  double res1, err1;
//	  _qt = qtmn;
//	  if (qtmn < opts.qtcutoff)
//	    res1 =  0.;
//	  else
//	    intdeo(besselint_dequad, 0.0, qtmn, awinf, &res1, &err1);

	  double res1, res2;
	  _qt = qtmx;

	  //if (opts.makecuts)
	  //  {
	  //    phasespace::set_qt(qtmx);
	  //    rapint::numint(ymn,ymx,phasespace::m);
	  //    mesq::setmesq_expy(mode, phasespace::m, _costh, _y);
	  //  }
	  res2 = bintegral(qtmx);
	  _qt = qtmn;
	  if (qtmn < opts.qtcutoff)
	    res1 =  0.;
	  else
	    {
	//if (opts.makecuts)
	//	{
	//	  phasespace::set_qt(qtmn);
	//	  rapint::numint(ymn,ymx,phasespace::m);
	//	  mesq::setmesq_expy(mode, phasespace::m, _costh, _y);
	//	}
	    res1 = bintegral(qtmn);
	    }
	  res = qtmx*res2 - qtmn*res1;

	  // Normalization
	  res *= 1./2./pow(opts.sroot,2);
	  //cout << " res " << res << endl;
	}
    }
  else
    {

//      //if (opts.fast && opts.xspace) {
//      //Fast function
//      res = exp(-_qt);
//
//      // Normalization
//      double shad = pow(opts.sroot,2);
//      res *= _qt/shad;
//      //}
//      //else
//      //{
      
      //*****************************************
      //dependence on qt, m, also y, and costh unless integrated
      //perform b integration (int_0^inf db)
      res = bintegral(_qt);
      //  cout << "phase space point in resumm" << "  " << _qt << "  " <<  _y << "  " <<  _m << "  " << _costh << "  " << res << endl;
      //*****************************************

      // Normalization
      res *= _qt/2./pow(opts.sroot,2);
    }


  //Do not need the Mellin transform for the LL case, because the HN coefficient is 1 --> Use PDFs in x space
  if (opts.order == 0 && opts.xspace && (opts.evolmode == 0 || opts.evolmode == 1)) //if (opts.order == 0 && opts.xspace && (opts.evolmode == 1))
    if (mode < 2 || mode == 4) //rapidity differential
      res *= mesq::loxs(x1, x2, mufac);
    else           //rapidity integrated
      res *= mesq::loxs(tau, mufac);      

  /*
  if (opts.order == 0 && opts.xspace && (opts.evolmode == 0 || opts.evolmode == 1))
    {
      double fun = 0.;
      if (mode < 2) //rapidity differential
	{
	  //PDFs
	  double fx1[2*MAXNF+1],fx2[2*MAXNF+1];
	  fdist_(opts.ih1,x1,mufac,fx1);
	  fdist_(opts.ih2,x2,mufac,fx2);
	  
	  for (int sp = 0; sp < mesq::totpch; sp++)
	    fun += fx1[mesq::pid1[sp]]*fx2[mesq::pid2[sp]]*real(mesq::mesqij[sp]);
	}
      else //rapidity integration
	{
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
		  x1 = tau*exppy;
		  x2 = tau*expmy;
		      
		  //PDFs
		  double fx1[2*MAXNF+1],fx2[2*MAXNF+1];
		  fdist_(opts.ih1,x1,mufac,fx1);
		  fdist_(opts.ih2,x2,mufac,fx2);
		      
		  for (int sp = 0; sp < mesq::totpch; sp++)
		    fun += fx1[mesq::pid1[sp]]*fx2[mesq::pid2[sp]]*real(mesq::mesqij[sp])*gr::www[opts.yrule-1][j]*xm;;
		}
	    }
	}
      //cout << "resint " << fun << endl;
      res *= fun;
    }
  */
  return res;
}

double resint::bintegral(double qt)
{
  double res, err;

  //Warning!!! when evolmode = 4, the integrand of bessel transform is discontinous at mb and mc
  //--> Need to split the integral in this case
      
  // Compute integral using double exponential quadratures for oscillatory fuctions (intde2.c)

  //bstar prescription
  if (opts.bprescription == 0 || opts.bprescription == 4)
    {
      //intdeo(besselint::bint, 0.0, qt, awinf, &res, &err);
      intdeo(besselint_dequad, 0.0, qt, awinf, &res, &err);
      if (err < 0)
      	cout << "warning: dequad abnormal termination, pt=" << _qt << " m=" << _m << " y=" << _y << " bint: " << res << " err " << err << endl;

      /*
      //test alternative integrations
      //      cout << endl;
      //      cout << "dequad result of inverse bessel transform, pt=" << _qt << " m=" << _m << " y=" << _y << " : " << setprecision(16) << res << " +- " << err << endl;

      int bqrule = 20;
      complex <double> h1 = 0;
      complex <double> h2 = 0;
      complex <double> bb;
      for (int i = 0; i < bqrule; i++)
	{
	  bb = complex <double> {0,bq::xxx[bqrule-1][i]/qt};
	  h1 += 2./M_PI/qt * bq::www[bqrule-1][i] * besselint::bint(bb);
	  //	  cout << "bb " << bb << " h1  " << h1 << endl;
	  bb = complex <double> {0,-bq::xxx[bqrule-1][i]/qt};
	  h2 += 2./M_PI/qt * bq::www[bqrule-1][i] * besselint::bint(bb);
	  //	  cout << "bb " << bb << " h2  " << h2 << endl;
	}
      res = real(1./2.*(h1+h2));

      //      cout << "bequad result " << res << endl;
      //      cout.precision(6); cout.unsetf(ios_base::floatfield);
      /*
     hankel::init(0,0,0.1);
     hankel::transform(besselint_hankel, qt, res, err);
     hankel::free();
     cout << "hankel result " << res << " +- " << err << endl;
     cout.precision(6); cout.unsetf(ios_base::floatfield);

      //input
      int nb = 1; //lagged convolutions
      int nrel = 1; //related convolutions
      int ntol = 1; 
      int nord[1] = {0}; //order of the transform
      int ijrel[2*1]; //exponents of the related convolutions (not used)
      double dwork [1*801]; //work array
      //output
      double dans[1*1]; //result
      double arg[1]; //arguments of the lagged convolutions
      int nofun1; //number of function evalutions
      int ierr; //error code
      dhankl_(qt, nb, nrel, opts.bintaccuracy, ntol, nord, besselint_588, ijrel, dwork, dans, arg, nofun1, ierr);
      res = dans[0];
      err = 0.;
      cout << "588 result " << res << " +- " << err << "  " << nofun1 << endl;
      cout.precision(6); cout.unsetf(ios_base::floatfield);
      */
	   
    }
  //Integrate up to blim with besche
  else if (opts.bprescription == 1)
    {
      double bmax;
      if (opts.modlog)
	bmax = resconst::b0/scales::res * sqrt(exp(1./(resint::aass*resconst::beta0))-1.);
      else
	bmax = resconst::b0/scales::res * exp(1./(2.*resint::aass*resconst::beta0));

      //cout << endl;
      //cout << "pt " << _qt << " m " << _m << endl;
      double C = bmax; //blim::sudakov;  //blimit_.rblim_; //upper limit of integration
      double ALFA[1] = {_qt};    //oscillation frequency
      int NUM = 1;               //number of integrals to be computed
      int NU = (_mode != 3 && _mode != 4)?0:1; //order of the bessel function
      int N = 30;                //degree of the chebyshev approximation
      double RESULT[1];          //output result
      int INFO[1];               //output status code
      besche_(besselint_besche_,C,ALFA,NUM,NU,N,RESULT,INFO);
      res = RESULT[0];
      err = 0;
      //cout << "besche result " << res << " status " << INFO[0] << endl;
    }
      
  //Minimal prescription to avoid Landau singularity (see hep-ph/9604351 and hep-ph/0002078)
  else if (opts.bprescription == 2)
    {
      res = 0;

      //Noo that is not the pole!!!
      //bc =  opts.bcf*blim::sudakov; //blimit_.rblim_;

      double bmax;
      if (opts.modlog)
	//bmax = resconst::b0/scales::res * (exp(1./(2.*resint::aass*resconst::beta0))-1.);
      //I think the correct formula is:
	bmax = resconst::b0/scales::res * sqrt(exp(1./(resint::aass*resconst::beta0))-1.);
      else
	bmax = resconst::b0/scales::res * exp(1./(2.*resint::aass*resconst::beta0));

      //Empirically seems that this adjustment is needed:
      //bmax *= scales::res/scales::ren;

      bc =  opts.bcf*bmax;

      /*
      //set bc so that the branching point is at lambda = opts.bcf
      if (opts.modlog)
	bc = resconst::b0/scales::res * sqrt(exp(opts.bcf/(resint::aass*resconst::beta0))-1.);
      else
	bc = resconst::b0/scales::res * exp(opts.bcf/(2.*resint::aass*resconst::beta0));
      */
      

      //if (opts.bcf > 0)
      if (opts.bcf > 1)
	{
	  //cout << endl;
	  //cout << "pt " << _qt << " m " << _m << endl;
	  
	  //besche from 0 to bc
	  double C = bc;             //upper limit of integration
	  double ALFA[1] = {_qt};    //oscillation frequency
	  int NUM = 1;               //number of integrals to be computed
	  int NU = (_mode != 3 && _mode != 4)?0:1; //order of the bessel function
	  int N = 30;                //degree of the chebyshev approximation
	  double RESULT[1];          //output result
	  int INFO[1];               //output status code
	  besche_(besselint_besche_,C,ALFA,NUM,NU,N,RESULT,INFO);
	  res = RESULT[0];
	  err = 0;
	  //cout << "besche result up to bc " << res << " status " << INFO[0] << endl;
	}
      
      /*
	double min = 0.;
	double max = bc+blimit_.rblim_/_qt;
	//adzint settings
	double abserr = 1e-10;
	double relerr = 1e-6;
	double errest;
	double ier;
	double iacta = 0;
	double iactb = 0;

	//res = 0.;
	complex <double> r = 0.;
	complex <double> i(0.,1.);
	  
	//upper real integral
	//r += adzint_(fureal_,min,max,abserr,relerr,errest,ier,iacta,iactb);
  
	//lower real integral
	//r += adzint_(fdreal_,min,max,abserr,relerr,errest,ier,iacta,iactb);

	r += adzint_(freal_,min,max,abserr,relerr,errest,ier,iacta,iactb);

	//res += real(r)/2.;
	cout << "adzint result " << real(r)/2. << " error " << errest << endl;
      */
	  
      double resdequad;
      double errdequad;

      /*
	complex <double> r = 0.;
	intdeo(frealu_dequad, 0, qt*cos(M_PI/opts.phibr), awinf, &resdequad, &errdequad);
	r += resdequad;	  
	cout << "dequad up result " << resdequad << " error " << errdequad << endl;
	
	intdeo(freald_dequad, 0, qt*cos(M_PI/opts.phibr), awinf, &resdequad, &errdequad);
	r += resdequad;	  
	cout << "dequad dn result " << resdequad << " error " << errdequad << endl;
	
	res += real(r)/2.;
	cout << "Split dequad result " << real(r)/2. << " error " << errdequad << endl;
      */
      
      //intdeo(besselint_mp_complex_dequad, 0, qt*cos(M_PI/opts.phibr), awinf, &resdequad, &errdequad);
      intdeo(besselint_mp_complex_dequad, 0, qt, awinf, &resdequad, &errdequad);
      //cout << "Complex dequad result " << resdequad << " error " << errdequad << endl;
      if (errdequad < 0)
	cout << "warning: dequad abnormal termination, pt=" << _qt << " m=" << _m << " y=" << _y << " bint: " << resdequad << " err " << errdequad << endl;
      res += resdequad;
      
      //cout << "Total bp = 2 " << res << endl;

      /*
      //test alternative integrations
      //      cout << endl;
      //      cout << "dequad result of inverse bessel transform, pt=" << _qt << " m=" << _m << " y=" << _y << " : " << setprecision(16) << resdequad << " +- " << errdequad << endl;

      int bqrule = 40;
      complex <double> h1 = 0;
      complex <double> h2 = 0;
      complex <double> bb;
      for (int i = 0; i < bqrule; i++)
       {
	 bb = complex <double> {0,bq::xxx[bqrule-1][i]/qt};
	 h1 += 2./M_PI/qt * bq::www[bqrule-1][i] * besselint::bint(bb);
	 //	 cout << "bb " << bb << " h1  " << h1 << endl;
	 bb = complex <double> {0,-bq::xxx[bqrule-1][i]/qt};
	 h2 += 2./M_PI/qt * bq::www[bqrule-1][i] * besselint::bint(bb);
	 //	 cout << "bb " << bb << " h2  " << h2 << endl;
       }
      res = real(1./2.*(h1+h2));

      //      cout << "bequad result " << res << endl;
      //      cout.precision(6); cout.unsetf(ios_base::floatfield);
      */
      
    }
  //Minimal prescription along the real axis, crossing the Landau singularity
  else if (opts.bprescription == 3)
    {
      //minimal prescription with real b integration (not always more efficient, but save the calculation of PDFs at complex scales)
      res = 0;
      bc = 0.;
      /*
      if (opts.bcf > 0)
	{
	  bc =  opts.bcf*blimit_.rblim_;
	  
	  //besche from 0 to bc
	  double C = bc;             //upper limit of integration
	  double ALFA[1] = {_qt};    //oscillation frequency
	  int NUM = 1;               //number of integrals to be computed
	  int NU = (_mode != 3 && _mode != 4)?0:1; //order of the bessel function
	  int N = 30;                //degree of the chebyshev approximation
	  double RESULT[1];          //output result
	  int INFO[1];               //output status code
	  besche_(besselint_besche_,C,ALFA,NUM,NU,N,RESULT,INFO);
	  res = RESULT[0];
	  err = 0;
	}
      */

      double resdequad;
      double errdequad;
      intdeo(besselint_mp_real_dequad, 0, qt, awinf, &resdequad, &errdequad);
      //cout << "Real dequad result " << resdequad << " error " << errdequad << endl;
      res += resdequad;
      
      //cout << "Total bp = 3 " << res << endl;
    }
      
  //  //vfn integration --> Not correct!!!
  //  //split the integral above and below the b and c masses
  //
  //  //scale b0/b without a_param
  //  double bstar_mb = resconst::b0/(LHAPDF::getThreshold(5)*opts.kmub);
  //  //convert bstar to b
  //  double b_mb = bstar_mb / sqrt(1-(bstar_mb*bstar_mb)/(blimit_.rblim_*blimit_.rblim_));
  //
  //  //same for charm mass
  //  double bstar_mc = resconst::b0/(LHAPDF::getThreshold(4)*opts.kmuc);
  //  double b_mc = bstar_mc / sqrt(1-(bstar_mc*bstar_mc)/(blimit_.rblim_*blimit_.rblim_));
  //  
  //  //Integrate from pt=0 to pt=mc (from b_mc to infinity)
  //  //Set NF = 3
  //  sudakov::setnf(3);
  //  double res1, err1;
  //  //intdeo(besselint::bint, b_mc, qt, awinf, &res1, &err1);
  //  intdeo(besselint_dequad, b_mc, qt, awinf, &res1, &err1);
  //
  //
  //  //Integrate from pt=mb to pt=infinity (from 0 to b_mb)
  //  //Set NF = 5
  //  sudakov::setnf(5);
  //  double res2, err2;
  //  //dequad seems not to be appropriate for this integral 
  //  //  intde(besselint::bint, 0.0, b_mb, awfin, &res2, &err2);
  //
  //  //try simple pcubature
  //  const int ndim = 1;     //dimensions of the integral
  //  const int ncomp = 1;  //components of the integrand
  //  void *userdata = NULL;
  //  double integral[1];
  //  double error[1];
  //  const int eval = 0;
  //  const double epsrel = 1e-2;//opts.bintaccuracy;
  //  const double epsabs = 0.;
  //  double xmin[1] = {0.000001};
  //  double xmax[1] = {b_mb};
  //  /*
  //    pcubature(ncomp, besselint_cubature, userdata, 
  //    ndim, xmin, xmax, 
  //    eval, epsabs, epsrel, ERROR_INDIVIDUAL, integral, error);
  //    res2 = integral[0];
  //    err2 = error[0];
  //    cout << endl;
  //    cout << "pcubature result " << res2 << " error " << err2 << endl;
  //  */
  //
  //  double C = b_mb;           //upper limit of integration
  //  double ALFA[1] = {_qt};    //oscillation frequency
  //  int NUM = 1;               //number of integrals to be computed
  //  int NU = 0;                //order of the bessel function
  //  int N = 30;                //degree of the chebyshev approximation
  //  double RESULT[1];          //output result
  //  int INFO[1];               //output status code
  //  besche_(besselint_besche_,C,ALFA,NUM,NU,N,RESULT,INFO);
  //  res2 = RESULT[0];
  //  err2 = 0;
  //
  //  //cout << "besche result " << res2 << " status " << INFO[0] << endl;
  //      
  //  //Integrate from pt=mc to pt=mb (from b_mb to b_mc)
  //  //Set NF = 4
  //  sudakov::setnf(4);
  //  double res3, err3;
  //  xmin[0] = {b_mb};
  //  xmax[0] = {b_mc};
  //  pcubature(ncomp, besselint_cubature, userdata, 
  //	    ndim, xmin, xmax, 
  //	    eval, epsabs, epsrel, ERROR_INDIVIDUAL, integral, error);
  //  res3 = integral[0];
  //  err3 = error[0];
  //      
  //  res = res1+res2+res3;
  //
  //  /*
  //    cout << "split: dequad result of inverse bessel transform, " << endl;
  //    cout << "total: "
  //    //<< setprecision(16)
  //    << res << "  " << (err1+err2+err3)/res*100 << endl;
  //    cout << " b(mc)-inf: " << res1 << "  " << err1/res1*100 << "%"
  //    << " 0-b(mb): "  << res2 << "  " << err2/res2*100 << "%"
  //    << " b(mb)-b(mc): "  << res3 << "  " << err3/res3*100 << "%"
  //    << endl;
  //    //cout.precision(6); cout.unsetf(ios_base::floatfield);
  //    //cout << "b_mb " << b_mb << " b_mc " << b_mc << " blim " << blimit_.rblim_ << endl;
  //    cout << endl;
  //  */


  /*
  //plot the b-integrand
  cout << "{" << endl;
  cout << "TGraph *g = new TGraph();" << endl;
  for (int i = 0; i < 1000; i++)
    {
      //      double b = b_mb*0.95+i*(b_mb*1.05 -b_mb*0.95)/100.;
      double b = 0.+i*(b_mb*5 -0)/1000.;
      //      cout << b << "  " << b_mb << "  " << besselint::bint(b) << endl;;
      cout << "g->SetPoint(g->GetN(), " << b << ", " << besselint::bint(b) << ");" << endl;
    }
  cout << "g->Draw();" << endl;
  cout << "}" << endl;
  */

  /*
  cout << "{" << endl;
  cout << "TGraph *g = new TGraph();" << endl;
  double bmax = bc/opts.bcf;
  cout << "//bmax = " << bmax << endl;
  for (int i = 0; i < 1000; i++)
    {
      double x = i/1000. * bmax * 2;
      complex <double> b;
      //      double b = b_mb*0.95+i*(b_mb*1.05 -b_mb*0.95)/100.;
      //double b = 0.+i*(b_mb*5 -0)/1000.;
      if (x < resint::bc)
	b = x;
      else
	{
	  complex <double> jacu = complex <double> (1.,tan(M_PI/opts.phibr));
	  b = resint::bc + jacu*(x-resint::bc);
	}
      //      cout << b << "  " << b_mb << "  " << besselint::bint(b) << endl;;
      //cout << "g->SetPoint(g->GetN(), " << b/bmax << ", " << besselint::bint(b) << ");" << endl;
      cout << "g->SetPoint(g->GetN(), " << real(b)/bmax << ", " << real(besselint::bint(b)) << ");" << endl;
    }
  cout << "g->Draw();" << endl;
  cout << "}" << endl;
  exit(0);  
  */  
  return res;
}
