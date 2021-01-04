#include "sudakovff.h"
#include "resconst.h"
#include "resint.h"
#include "blim.h"
#include "scales.h"
#include "besselint.h"
#include "pegasus.h"
#include "npff.h"
#include "isnan.h"
#include "gint.h"
#include "alphas.h"
#include <iostream>
#include <iomanip>

using namespace std;
using namespace resconst;
using namespace resint;

const double pi2 = M_PI*M_PI;
const double pi3 = pi2*M_PI;
const double pi4 = pi2*pi2;

//int sudakov::nf;
//double sudakov::beta0;
//double sudakov::beta1;
//double sudakov::beta2;
complex <double> sudakov::log1y;
//double sudakov::ry;
complex <double> sudakov::logS;
complex <double> sudakov::S;

//void sudakov::setnf(int nff)
//{
//  nf = nff;
//  beta0=(33.-2.*nf)/12.;
//  beta1=(153.-19.*nf)/24.;
//  beta2=2857./128.-5033.*nf/1152.+325.*nf*nf/3456.;
//}



complex <double> sudakov::f0(complex <double> y)
{
  //FFN beta
  //return (resconst::A1q/beta0(fabs(b)))*(y+log1y)/(y);
  //VFN beta
  //cout << "f0" << "  " << A1q << "  " << beta0 << "  " << y << " " << log1y << "  " << (A1q/beta0)*(y+log1y)/(y) << endl;

  return (A1q/beta0)*(y+log1y)/(y);
}

complex <double> sudakov::f1(complex <double> y) //Eq.23 of arXiv:hep-ph/0508068
{
  complex <double> F1 = 
    ((A1q*beta1)/(pow(beta0,3)))*((1./2.)*log1y*log1y +
					    (y)/(1.-y)+log1y/(1.-y)) -
    (A2q/(pow(beta0,2)))*(log1y+(y)/(1.-y)) + 
    (B1q/beta0)*log1y +
    (A1q/beta0)*(y/(1.-y)+log1y)*resint::rlogq2mur2; // !!! should be rlogq2mur2 -> log(mures^2/mur^2) = (logmur2q2-2.*loga) (is this a bug in DYRES?)
    //    a dependence      
  F1=F1-2.*resint::rloga*A1q/beta0*y/(1.-y);
  return F1;
}
complex <double> sudakov::f2(complex <double> y)
{
  complex <double> F2 = 
    ((A2q*beta1)/(pow(beta0,3)))*((y/2.)*((3.*y-2.)/pow(1.-y,2))-
				  ((1.-2.*y)*log1y/(1.-y)/(1.-y)))
    -(B2q/beta0)*((y)/(1.-y))
    +(B1q*beta1/pow(beta0,2))*((y)/(1.-y)+log1y/(1.-y))
    -(A3q/2./pow(beta0,2))*(y)*(y)/(1.-y)/(1.-y)
    +A1q*((pow(beta1,2)/2./pow(beta0,4))*(1.-2.*y)/(1.-y)
	  /(1.-y)*log1y*log1y +
	  log1y*((beta0*beta2-pow(beta1,2))/(pow(beta0,4))+
		 pow(beta1,2)/pow(beta0,4)/(1.-y)) +
	  // pow(beta1,4)/pow(beta0,4)/(1.-y)) +
	  (y)/(2.*pow(beta0,4)*(1.-y)*(1.-y))*
	  (beta0*beta2*(2.-3.*y)+pow(beta1,2)*y)) -
    (A1q/2.)*(y)*(y)/(1.-y)/(1.-y)*resint::rlogq2mur2*resint::rlogq2mur2 + //should be rlogq2mur2 -> log(mures^2/mur^2) = (logmur2q2-2.*loga) (is this a bug in DYRES?)
    resint::rlogq2mur2*(B1q*y/(1.-y)+A2q/beta0*y*y/(1.-y)/(1.-y)+
			A1q*beta1/pow(beta0,2)*(y/(1.-y)+(1.-2.*y)/(1.-y)/(1.-y)*log1y))
    +2.*C1qqn*((y)/(1.-y));
    
  //    a dependence (now without constant term)
  F2 += 2*A1q*y*(y-2.)/pow(1.-y,2)*pow(resint::rloga,2)-resint::rloga
    *(2.*B1q*y/(1.-y)+2.*y/beta0*A2q/pow(1.-y,2)
      -2.*A1q*beta1/pow(beta0,2)*y*log1y/pow(1.-y,2))
    + A1q*resint::rloga*resint::rlogq2mur2*y*2./pow(1.-y,2);
    
    return F2;
}

complex <double> sudakov::g1(complex <double> y)
{
  if (abs(y) < 1e-16)
    return A1q/beta0;
  
  return (A1q*(y + log1y))/(beta0*y);
  //return (A1q*(1. + log1y/y))/(beta0);
  
  //return (A1q*(ry + log1y))/(beta0*ry);
}
complex <double> sudakov::g2(complex <double> y)
{
  double B1qbar = B1q + A1q*LQ;
  //double LQR = 0.;
  double LQR = LR-LQ;
  return (A2q*(-(y/(1. - y)) - log1y))/pow(beta0,2) + (B1qbar*log1y)/beta0 + (A1q*beta1*(y/(1. - y) + log1y/(1. - y) + pow(log1y,2)/2.))/pow(beta0,3)
    + (A1q*(-y + (-1. + y)*log1y)*(LQR))/(beta0*(-1. + y));

  //return
  //  + (A2q*(-(ry/(1. - y)) - log1y))/pow(beta0,2)
  //  + (B1qbar*log1y)/beta0
  //  + (A1q*beta1*(ry/(1. - y) + log1y/(1. - y) + pow(log1y,2)/2.))/pow(beta0,3)
  //  + (A1q*(-ry + (-1. + ry)*log1y)*(LQR))/(beta0*(-1. + y));

}
complex <double> sudakov::g3(complex <double> y)
{
  double B1qbar = B1q + A1q*LQ;
  double B2qbar = B2q + A2q*LQ;
  //double LQR = 0.;
  double LQR = LR-LQ;
  return (A1q*beta2*y)/(pow(beta0,3)*pow(1. - y,2)) - (3.*A1q*beta2*pow(y,2))/(2.*pow(beta0,3)*pow(1. - y,2)) - (A3q*pow(y,2))/(2.*pow(beta0,2)*pow(-1. + y,2)) + (A2q*beta1*y*(-2. + 3.*y))/(2.*pow(beta0,3)*pow(1. - y,2)) - (B2qbar*y)/(beta0 - beta0*y) + (A1q*beta2*log1y)/(pow(beta0,3)*pow(1. - y,2)) - (2.*A1q*beta2*y*log1y)/(pow(beta0,3)*pow(1. - y,2)) + (A1q*beta2*pow(y,2)*log1y)/(pow(beta0,3)*pow(1. - y,2)) + (A2q*beta1*(-1. + 2.*y)*log1y)/(pow(beta0,3)*pow(1. - y,2)) + (B1qbar*beta1*(y + log1y))/(pow(beta0,2)*(1. - y)) - (A1q*pow(beta1,2)*(y + log1y)*(-y + (-1. + 2.*y)*log1y))/(2.*pow(beta0,4)*pow(-1. + y,2))
    - ((LQR)*(2.*A1q*beta1*(-1. + 2.*y)*log1y + y*(2.*B1qbar*pow(beta0,2)*(-1. + y) + 2.*A1q*beta1*(-1. + y) - 2.*A2q*beta0*y + A1q*pow(beta0,2)*y*(LQR))))/(2.*pow(beta0,2)*pow(-1. + y,2));

  //return
  //  +(A1q*beta2*ry)/(pow(beta0,3)*pow(1. - y,2))
  //  -(3.*A1q*beta2*pow(ry,2))/(2.*pow(beta0,3)*pow(1. - y,2))
  //  -(A3q*pow(ry,2))/(2.*pow(beta0,2)*pow(-1. + y,2))
  //  +(A2q*beta1*ry*(-2. + 3.*ry))/(2.*pow(beta0,3)*pow(1. - y,2))
  //  -(B2qbar*ry)/(beta0 - beta0*y)
  //  +(A1q*beta2*log1y)/(pow(beta0,3)*pow(1. - y,2))
  //  -(2.*A1q*beta2*ry*log1y)/(pow(beta0,3)*pow(1. - y,2))
  //  +(A1q*beta2*pow(ry,2)*log1y)/(pow(beta0,3)*pow(1. - y,2))
  //  +(A2q*beta1*(-1. + 2.*ry)*log1y)/(pow(beta0,3)*pow(1. - y,2))
  //  +(B1qbar*beta1*(ry + log1y))/(pow(beta0,2)*(1. - y))
  //  -(A1q*pow(beta1,2)*(ry + log1y)*(-ry + (-1. + 2.*ry)*log1y))/(2.*pow(beta0,4)*pow(-1. + y,2))
  //  
  //  - (
  //     (LQR)*(
  //	      + 2.*A1q*beta1*(-1. + 2.*ry)*log1y
  //	      + ry*(
  //		   + 2.*B1qbar*pow(beta0,2)*(-1. + y)
  //		   + 2.*A1q*beta1*(-1. + y)
  //		   - 2.*A2q*beta0*ry
  //		   + A1q*pow(beta0,2)*ry*(LQR)
  //		   )
  //	      )
  //     )
  //  /(2.*pow(beta0,2)*pow(-1. + y,2));
  
}
complex <double> sudakov::g4(complex <double> y)
{
  double B1qbar = B1q + A1q*LQ;
  double B2qbar = B2q + A2q*LQ;
  double B3qbar = B3q + A3q*LQ;

  //double LQR = 0.;
  double LQR = LR-LQ;
  
  bool term = true;  //switch off one term which gives divergent Sudakov
  return (B3qbar*(-2. + y)*y)/(2.*beta0*pow(1. - y,2)) + (A4q*(-3. + y)*pow(y,2))/(6.*pow(beta0,2)*pow(1. - y,3)) - (A3q*beta1*y*(6. - 15.*y + 5.*pow(y,2)))/(12.*pow(beta0,3)*pow(1. - y,3)) + (A3q*beta1*(-1. + 3.*y)*log1y)/(2.*pow(beta0,3)*pow(1. - y,3)) + (B2qbar*beta1*((2. - y)*y + 2.*log1y))/(2.*pow(beta0,2)*pow(1. - y,2)) + (B1qbar*((pow(beta1,2) - beta0*beta2)*pow(y,2) - pow(beta1,2)*pow(log1y,2)))/(2.*pow(beta0,3)*pow(1. - y,2)) - (A2q*(y*(8.*beta0*beta2*pow(y,2) + pow(beta1,2)*(-6. + (9. - 11.*y)*y)) + 6.*pow(beta1,2)*log1y*(-1. + y + (-1. + 3.*y)*log1y)))/(12.*pow(beta0,4)*pow(1. - y,3)) - (A1q*(y*(2.*pow(beta1,3)*pow(y,2) + pow(beta0,2)*beta3*(-6. + (15. - 7.*y)*y) + beta0*beta1*beta2*(6. + 5.*(-3. + y)*y)) + 2.*log1y*(3.*(-(pow(beta0,2)*beta3*pow(1. - y,3)) + pow(beta1,3)*pow(y,2)*(1. + y) + beta0*beta1*beta2*(1. + y*(-3. + 2.*(1. - y)*y))) + pow(beta1,3)*log1y*(3.*y + (1. - 3.*y)*log1y))))/(12.*pow(beta0,5)*pow(1. - y,3))
    + ((LQR)*(3.*y*(-2.*B2qbar*pow(beta0,3)*(-2. + y)*(-1. + y) - A2q*beta0*beta1*(-2. + y)*(-1. + y) + y*(A3q*pow(beta0,2)*(-3. + y) + A1q*(pow(beta1,2) - beta0*beta2)*(1. + y))) + 3.*beta1*log1y*(2.*(-(B1qbar*pow(beta0,2)*(-1. + y)) + A1q*beta1*y + A2q*beta0*(-1. + 3.*y)) + term*A1q*beta1*(1. - 3.*y)*log1y) + 3.*pow(beta0,2)*(y*(A1q*beta1 + B1qbar*pow(beta0,2)*(-2. + y)*(-1. + y) - A2q*beta0*(-3. + y)*y) + A1q*beta1*(1. - 3.*y)*log1y)*(LQR) + A1q*pow(beta0,4)*(-3. + y)*pow(y,2)*pow(LQR,2)))/(6.*pow(beta0,3)*pow(-1. + y,3));

  //return 
  //  +(B3qbar*(-2. + ry)*ry)/(2.*beta0*pow(1. - y,2))
  //  +(A4q*(-3. + ry)*pow(ry,2))/(6.*pow(beta0,2)*pow(1. - y,3))
  //  -(A3q*beta1*ry*(6. - 15.*ry + 5.*pow(ry,2)))/(12.*pow(beta0,3)*pow(1. - y,3))
  //  +(A3q*beta1*(-1. + 3.*ry)*log1y)/(2.*pow(beta0,3)*pow(1. - y,3))    
  //  +(B2qbar*beta1*((2. - ry)*ry + 2.*log1y))/(2.*pow(beta0,2)*pow(1. - y,2))
  //  +(B1qbar*((pow(beta1,2) - beta0*beta2)*pow(ry,2) - pow(beta1,2)*pow(log1y,2)))/(2.*pow(beta0,3)*pow(1. - y,2))
  //  -(A2q*(ry*(8.*beta0*beta2*pow(ry,2) + pow(beta1,2)*(-6. + (9. - 11.*ry)*ry)) + 6.*pow(beta1,2)*log1y*(-1. + y + (-1. + 3.*ry)*log1y)))/(12.*pow(beta0,4)*pow(1. - y,3))
  //  -(A1q*(ry*(2.*pow(beta1,3)*pow(ry,2) + pow(beta0,2)*beta3*(-6. + (15. - 7.*ry)*ry) + beta0*beta1*beta2*(6. + 5.*(-3. + ry)*ry)) + 2.*log1y*(3.*(-(pow(beta0,2)*beta3*pow(1. - y,3)) + pow(beta1,3)*pow(ry,2)*(1. + ry) + beta0*beta1*beta2*(1. + ry*(-3. + 2.*(1. - y)*ry))) + pow(beta1,3)*log1y*(3.*ry + (1. - 3.*ry)*log1y))))/(12.*pow(beta0,5)*pow(1. - y,3))
  //
  //  + (
  //     (LQR)*(
  //	      +3.*ry*(-2.*B2qbar*pow(beta0,3)*(-2. + ry)*(-1. + y)
  //		    - A2q*beta0*beta1*(-2. + ry)*(-1. + y)
  //		    + ry*(A3q*pow(beta0,2)*(-3. + ry)
  //			 + A1q*(pow(beta1,2) - beta0*beta2)*(1. + y)))
  //	      + 3.*beta1*log1y*(
  //				+2.*(0.
  //				     -(B1qbar*pow(beta0,2)*(-1. + y))
  //				     + A1q*beta1*ry
  //				     + A2q*beta0*(-1. + 3.*ry)
  //				    )
  //				+ term*A1q*beta1*(1. - 3.*ry)*log1y
  //				)
  //	      + 3.*pow(beta0,2)*(ry*(A1q*beta1 + B1qbar*pow(beta0,2)*(-2. + ry)*(-1. + y)
  //				    - A2q*beta0*(-3. + ry)*ry)
  //				 + A1q*beta1*(1. - 3.*ry)*log1y)*(LQR)
  //	      + A1q*pow(beta0,4)*(-3. + ry)*pow(ry,2)*pow(LQR,2)
  //	      )
  //     )/(6.*pow(beta0,3)*pow(-1. + y,3))
  //  ;
}


//Sudakov form factor
complex <double> sudakov::sff(complex <double> b)
{
  //Return 0 if b is close to the Landau pole
  double b_L;// = a_param_.b0p_*(1./resint::_m)*exp(1./(2.*resint::aass*resconst::beta0));
  if (opts.modlog)
    //b_L = b0/scales::res * (exp(1./(2.*aass*beta0))-1.);
  //I think the correct formula is:
    b_L = b0/scales::res * sqrt(exp(1./(aass*beta0))-1.);
  else
    b_L = b0/scales::res * exp(1./(2.*aass*beta0));

  if (fabs(b-b_L) < 1e-8)
    return 0;

  //In reading these formulas, notice that blog is L = log(q*bstar/b0p) = log[(q/a_param)*bstar/b0] = log[Q * bstar/b0], according to Eq. (13) and (17) of hep-ph/0508068.
      
  //choose bstar (b) for real axis (complex plane) integration
  //mass dependence in blim --> mass dependence in bstar
  //complex <double> blim = cx(blimit_.cblim_); //blim=(1/q)*exp(1/(2*aass*beta0))
  double blim = blim::sudakov;
  complex <double> bstar;
  if (opts.bprescription == 0 || opts.bprescription == 4 || opts.bstar_sudakov)
    bstar = real(b)/sqrt(1.+pow(real(b)/blim,2));
  else //if (opts.bprescription == 1 || opts.bprescription == 2 || opts.bprescription == 3)
    bstar = b;

  //cout << "bstar " << bstar << endl;
  
  //mass dependence (q) in blog
  //complex <double> blog;
  //if (opts.modlog)
  //  blog = log(pow(scaleh_.q_*bstar/a_param_.b0p_,2)+1.); //modified sudakov
  //else
  //  blog = log(pow(scaleh_.q_*bstar/a_param_.b0p_,2));    //normal sudakov

  double Q = scales::res;
  complex <double> blog;
  if (opts.modlog)
    blog = log(pow(Q*bstar/b0,2) + 1.); //modified sudakov
  else
    blog = log(pow(Q*bstar/b0,2));   //normal sudakov

  //double rblog;
  //if (opts.modlog)
  //  rblog = log(pow(Q*real(bstar)/b0,2) + 1.); //modified sudakov
  //else
  //  rblog = log(pow(Q*real(bstar)/b0,2));   //normal sudakov
  
  //cout << modified_.imod_ << endl;
  //cout << "blog " << blog << endl;

  complex <double> y = beta0*resint::aass*blog;

  //Numerical integration of the Sudakov form factor
  if (opts.numsud || opts.order >= 4)
    {
      logS = gint::calc(b);
      S = exp(logS);
    }
  //Analytical solution for the Sudakov form factor
  else
    {
      //mass dependence in f0(y), f1(y), f2(y)
      //FFN beta
      //complex <double> y = resconst::beta0*aass_.aass_*blog;
      //VFN beta
      log1y = log(1.-y);
      csud_.log1y_ = fcx(log(1.-y));

      //ry = beta0*resint::aass*rblog;
      //ry = abs(y);
      //ry = beta0*resint::aass*real(blog);
      
      fcomplex fy = fcx(y);
      
      if (opts.order == 0)
	//S = exp(blog*cx(f0_(fy)));
	//S = exp(blog*f0(y));
	S = exp(blog*g1(y));
      else if (opts.order == 1)
	//S = exp(blog*cx(f0_(fy))+cx(f1_(fy)));
	//S = exp(blog*f0(y)+f1(y));
	S = exp(blog*g1(y)+g2(y));
      else if (opts.order == 2)
	//S = exp(blog*cx(f0_(fy))+cx(f1_(fy))+aass_.aass_*cx(f2_(fy)));
	//S = exp(blog*f0(y)+f1(y)+resint::aass*f2(y));
	//S = exp(blog*g1(y)+g2(y)+resint::aass*(g3(y)+2.*C1qqn*((y)/(1.-y))));
	S = exp(blog*g1(y)+g2(y)+resint::aass*g3(y));
      else if (opts.order == 3)
	//S = exp(blog*g1(y)+g2(y)+resint::aass*g3(y));
	S = exp(blog*g1(y)+g2(y)+resint::aass*g3(y)+pow(resint::aass,2)*g4(y));

      /*
	if (opts.order == 0)
	S = exp(blog*g1(y));
	else if (opts.order == 1)
	S = exp(blog*g1(y)+g2(y));
	else if (opts.order == 2)
	S = exp(blog*g1(y)+g2(y)+resint::aass*g3(y));
	else if (opts.order == 3)
	S = exp(blog*g1(y)+g2(y)+resint::aass*g3(y)+pow(resint::aass,2)*g4(y));
      */

      if (opts.order == 0)
	logS = blog*g1(y);
      else if (opts.order == 1)
	logS = blog*g1(y)+g2(y);
      else if (opts.order == 2)
	logS = blog*g1(y)+g2(y)+resint::aass*g3(y);
      else if (opts.order == 3)
	logS = blog*g1(y)+g2(y)+resint::aass*g3(y)+pow(resint::aass,2)*g4(y);

      //cout << " blog " << blog << " g1 " << g1(y) << " g2 " << g2(y) << endl;
      //  cout << " y " << y << "  " << cx(fy) << endl;

      //cout << setprecision(16);
      //cout << "f0(y) " << cx(f0_(fy)) << "  " << f0(y) << endl;
      //cout << "f1(y) " << cx(f1_(fy)) << "  " << f1(y) << endl;
      //cout << "f2(y) " << cx(f2_(fy)) << "  " << f2(y) << endl;

      //cout << setprecision(16);
      //cout << "f0(y) " << f0(y) << " g1 " << g1(y) << endl;
      //cout << "f1(y) " << f1(y) << " g2 " << g2(y) << endl;
      //cout << "f2(y) " << f2(y) << " g3 " << g3(y) << endl;
      //cout << "f3(y) " << 0     << " g4 " << g4(y) << endl;

      //compare analytical and numerical integrations
      //cout << setprecision(16);
      //cout << " b " << b << endl;
      //cout << " logS " << logS << endl;
      //cout << " gint " << gint::calc(b) << endl;
      //cout << endl;
    }

  if (opts.bprescription == 0 && fabs(S) > 100.)
    {
      fcomplex fbb = fcx(b);
      //cout << "Warning! Large Sudakov, S(b) = " << S <<"; for bstar = " << bstar << " fortran S(b) = " << cx(s_(fbb)) << endl;
      //S = 0.;
      //cout << "C++ " << b << "  " << bstar << "  " << blog << "  " << y << "  " << blog*f0(y) << "  " << f1(y) << "  " << resint::aass*f2(y) << "  " <<  S << endl;
      //cout << "C++ " << b << "  " << bstar << "  " << blog << "  " << y << "  " << blog*g1(y) << "  " << g2(y) << "  " << resint::aass*g3(y) << "  " << pow(resint::aass,2)*g4(y) << "  " << S << endl;
    }
  if (isnan_ofast(real(S)))
    {
      fcomplex fbb = fcx(b);
      cout << "Warning! Sudakov is nan, S(b) = " << S << " log(S) = " << logS << "; for b = " << b << " bstar = " << bstar << " Landau pole " << b_L << endl;// << " fortran S(b) = " << cx(s_(fbb)) << endl;
      //cout << "C++ " << b << "  " << bstar << "  " << blog << "  " << y << "  " << blog*f0(y) << "  " << f1(y) << "  " << resint::aass*f2(y) << "  " <<  S << endl;
      S = 0;
    }
  //cout << "C++ " << b << "  " << bstar << "  " << y << "  " << f0(y) << "  " << S << endl;
  
  //S = S*exp(-opts.g_param*pow(b,2));
  //S = S * npff::S(b,resint::_m,resint::x1,resint::x2);

  if (opts.sumlogs)
    return 1.;
  else
    return S;
}
