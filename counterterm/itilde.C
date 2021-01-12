#include "itilde.h"
#include "resint.h"
#include "resconst.h"
#include "besselint.h"
//#include "hankel.h"
#include "numbers.h"
#include "settings.h"
#include "dyres_interface.h"
#include <iostream>

#include "intde2_c.h"
double itilde::aw[lenaw];

double itilde::iqt;
double itilde::iq2;
double itilde::ib02;
int itilde::in;

void itilde::init()
{
  intdeoini(lenaw, tiny, 1e-8, aw);
}

double itilde::besselk(double x, int n)
{
  //fortran interface
  xmio_.xmio_ = x;
  switch (n)
    {
    case 1: return itilde_(one);   break;
    case 2: return itilde_(two);   break;
    case 3: return itilde_(three); break;
    case 4: return itilde_(four);  break;
    }
  return 0;
}

//Calculate Itilde from the integral expression:
//Itilde_n(qt/q)=q^2 int_0^inf db b/2 J0(bqt) ln^n(Q^2 b^2 / b0^2 + 1)
double itilde::calc(double qt, double q, int n)
{
  iqt = qt;
  iq2 = pow(q,2);
  in = n;
  ib02 = pow(resconst::b0,2);
  
  double res, err;
  intdeo(integrand, 0.0, iqt, aw, &res, &err);
  //  cout << "itilde " << res*iq2 << " +- " << err*iq2 << endl;
  //  hankel::init(0,0,0.015);
  //  hankel::transform(logn, iqt, res, err);
  //  hankel::free();
  //  cout << "itilde " << res*iq2 << " +- " << err*iq2 << endl;
    
  return res*iq2;
}

double itilde::integrand(double b)
{
  //cout << " iqt " << iqt << " iq2 " << sqrt(iq2) << " in " << in << endl;
  double qtb = b*iqt;
  double b2 = b*b;
  double val = b/2. * fort_besj0_(qtb);

  if (opts.modlog == 0)
    val *= pow(log(iq2*b2/ib02),in);
  else if (opts.modlog == 1)
    val *= pow(log(iq2*b2/ib02+1.),in);
  else if (opts.modlog == 2)
    val *= pow(1./opts.p*log(pow(iq2*b2/ib02,opts.p)+1.),in);

  //cout << " b " << b << " val " << val;
  return val;
}

double itilde::logn(double b)
{
  double b2 = b*b;
  double val = 1./2. * pow(log(iq2*b2/ib02+1.),in);
  return val;
}

double itilde::integral(double qt, double q, int n)
{

  iqt = qt;
  iq2 = pow(q,2);
  in = n;
  ib02 = pow(resconst::b0,2);
  
  double res, err;
  intdeo(integrand_int, 0.0, iqt, aw, &res, &err);
  //  cout << endl;
  //  cout << "qt " << qt << " n " << n << endl;
  //  cout << "itilde " << res*iq2 << " +- " << err*iq2 << endl;
  //  hankel::init(1,0,0.02);
  //  hankel::transform(lognob, iqt, res, err);
  //  hankel::free();
  //  cout << "itilde " << res*iq2 << " +- " << err*iq2 << endl;
  
  return res*iq2;
}

double itilde::integrand_int(double b)
{
  //cout << " iqt " << iqt << " iq2 " << sqrt(iq2) << " in " << in << endl;
  double qtb = b*iqt;
  double b2 = b*b;
  //double val = 1./2. * fort_besj1_(qtb) * pow(log(iq2*b2/ib02+1.),in);

  double val = 1./2. * fort_besj1_(qtb);
  if (opts.modlog == 0)
    val *= pow(log(iq2*b2/ib02),in);
  else if (opts.modlog == 1)
    val *= pow(log(iq2*b2/ib02+1.),in);
  else if (opts.modlog == 2)
    val *= pow(1./opts.p*log(pow(iq2*b2/ib02,opts.p)+1.),in);

  /*
  //(Eq. 3.8 of https://arxiv.org/pdf/1805.05916.pdf)
  double jac;
  //if (opts.modlog == 0)
  //  jac = 1.;
  //else if (opts.modlog == 1)
  //  jac = sqrt(iq2/(iq2+ib02/b2));
  //else if (opts.modlog == 2)
  //  jac = sqrt(pow(iq2*b2/ib02,opts.p)/(1.+pow(iq2*b2/ib02,opts.p)));
  if (opts.modlog == 0)
    jac = 1.;
  else if (opts.modlog == 1)
    jac = sqrt(iq2)/iqt/(1.+sqrt(iq2)/iqt);
  else if (opts.modlog == 2)
    jac = pow(sqrt(iq2)/iqt,opts.p)/(1.+pow(sqrt(iq2)/iqt,opts.p));
  val *= jac;
  */
  
  //cout << " b " << b << " val " << val;
  return val;
}

double itilde::lognob(double b)
{
  double b2 = b*b;
  double val = 1./2./b * pow(log(iq2*b2/ib02+1.),in);
  return val;
}
