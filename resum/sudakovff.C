#include "sudakovff.h"
#include "resconst.h"
#include "resint.h"
#include "besselint.h"
#include "pegasus.h"
#include "npff.h"
#include <iostream>
#include <iomanip>

using namespace std;
using namespace resconst;

//int sudakov::nf;
//double sudakov::beta0;
//double sudakov::beta1;
//double sudakov::beta2;
complex <double> sudakov::log1y;

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

complex <double> sudakov::f1(complex <double> y)
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



//Sudakov form factor
complex <double> sudakov::sff(complex <double> b)
{
  double b_L = a_param_.b0p_*(1./resint::_m)*exp(1./(2.*resint::aass*resconst::beta0));
  if (fabs(b-b_L) < 1e-8)
    return 0;
    
  //In reading these formulas, notice that blog is L = log(q*bstar/b0p) = log[(q/a_param)*bstar/b0] = log[Q * bstar/b0], according to Eq. (13) and (17) of hep-ph/0508068.
      
  //choose bstar (b) for real axis (complex plane) integration
  complex <double> bstar;
  if (opts.bprescription == 0) 
    {
      // mass dependence in blim --> mass dependence in bstar
      complex <double> blim = cx(blimit_.cblim_); //blim=(1/q)*exp(1/(2*aass*beta0))
      bstar=b/sqrt(1.+(pow(b,2))/(pow(blim,2)));
    }
  else if (opts.bprescription == 1 || opts.bprescription == 2 || opts.bprescription == 3)
    bstar = b;

  //cout << "bstar " << bstar << endl;
  
  //mass dependence (q) in blog
  complex <double> blog;
  if (opts.modlog)
    blog = log(pow(scaleh_.q_*bstar/a_param_.b0p_,2)+1.); //modified sudakov
  else
    blog = log(pow(scaleh_.q_*bstar/a_param_.b0p_,2));    //normal sudakov

  //cout << modified_.imod_ << endl;
  //cout << "blog " << blog << endl;
  
  //mass dependence in f0(y), f1(y), f2(y)
  //FFN beta
  //complex <double> y = resconst::beta0*aass_.aass_*blog;
  //VFN beta
  complex <double> y = beta0*aass_.aass_*blog;
  csud_.log1y_ = fcx(log(1.-y));
  log1y = log(1.-y);

  complex <double> S;
  fcomplex fy = fcx(y);
  
  if (opts.order == 0)
    //S = exp(blog*cx(f0_(fy)));
    S = exp(blog*f0(y));
  else if (opts.order == 1)
    //S = exp(blog*cx(f0_(fy))+cx(f1_(fy)));
    S = exp(blog*f0(y)+f1(y));
  else if (opts.order == 2)
    //S = exp(blog*cx(f0_(fy))+cx(f1_(fy))+aass_.aass_*cx(f2_(fy)));
    S = exp(blog*f0(y)+f1(y)+aass_.aass_*f2(y));
  //  cout << " y " << y << "  " << cx(fy) << endl;

  //cout << setprecision(16);
  //cout << "f0(y) " << cx(f0_(fy)) << "  " << f0(y) << endl;
  //cout << "f1(y) " << cx(f1_(fy)) << "  " << f1(y) << endl;
  //cout << "f2(y) " << cx(f2_(fy)) << "  " << f2(y) << endl;
  
  if (opts.bprescription == 0 && fabs(S) > 100.)
    {
      fcomplex fbb = fcx(b);
      cout << "Warning! Large Sudakov, S(b) = " << S <<"; for bstar = " << bstar << " fortran S(b) = " << cx(s_(fbb)) << endl;
      S = 0.;
    }
  //cout << "C++ " << b << "  " << bstar << "  " << y << "  " << f0(y) << "  " << S << endl;
  
  //S = S*exp(-opts.g_param*pow(b,2));
  //S = S * npff::S(b,resint::_m,resint::x1,resint::x2);

  return S;
}
