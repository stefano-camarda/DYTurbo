#include "npff.h"
#include "settings.h"
#include "resconst.h"
#include "pdfevol.h"
#include "resint.h"
#include <LHAPDF/LHAPDF.h>
#include <iostream>

complex <double> npff::uvff;
complex <double> npff::usff;
complex <double> npff::dvff;
complex <double> npff::dsff;
complex <double> npff::ssff;
complex <double> npff::chff;
complex <double> npff::boff;
complex <double> npff::glff;

complex <double> npff::S(complex <double> b, double m, double x1, double x2)
{
  complex <double> ff;
    
  double a;
  double g1, g2, g3;
  double Q0;
  
  //Standard dyres 1-parameter
  a = opts.g_param;

  //BLNY (blim = 0.5)
  //g1 = 0.21;
  //g2 = 0.68;
  //g3 = -0.13;
  //Q0 = 3.2;

  //Konychev Nadolsky (blim = 1.5)
  //g1 = 0.20;
  //g2 = 0.19;
  //g3 = -0.03;
  //Q0 = 3.2;  

  //M. Hirai, H. Kawamura, K. Tanaka, https://inspirehep.net/record/1229722 (minimal prescription, MSTW08)
  //g1 = 0.330;
  //g2 = 0.066;
  //g3 = 0.;
  //Q0 = 2.*1.3;
  
  //g1 = 0.50;
  //g2 = 2.0;
  //g3 = 0.;
  //Q0 = 50;

  
  //a = g1 + g2 * log(m/Q0) + g3*log(100*x1*x2);

  
  //cout << "m " << m << " g " << a << endl;
  //cout << a1 << "  " << a2 * log(m/Q0) << "  " << a3*log(100*x1*x2) << endl;

  //gaussian form
  ff = exp(-a*pow(b,2));

  //exponential form
  //ff = exp(-a*b);


  //Collins, Rogers, large bT form (https://arxiv.org/abs/1412.3820)
  //double g0 = 0.3;
  //double gk = g0 * (1. - exp(- (resconst::Cf*LHAPDF::alphasPDF(real(pdfevol::qbstar))*pow(real(b),2))/(M_PI*g0*pow(blimit_.rblim_,2))));
  //Q0 = 1.6;
  //ff = exp(-gk*log(pow(m/Q0,2)));

  //flavor dependent form factors
  if (opts.flavour_kt)
    {
      uvff = exp(-a/2*pow(b,2));
      usff = exp(-a/2*pow(b,2));
      dvff = exp(-a/2*pow(b,2));
      dsff = exp(-a/2*pow(b,2));
      ssff = exp(-a/2*pow(b,2));
      chff = exp(-a/2*pow(b,2));
      boff = exp(-a/2*pow(b,2));
      glff = exp(-a/2*pow(b,2));
    }

  //universal form factor
  return ff;
}
