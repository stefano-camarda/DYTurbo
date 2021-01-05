#include "npff.h"
#include "settings.h"
#include "resconst.h"
#include "constants.h"
#include "pdfevol.h"
#include "resint.h"
//#include <LHAPDF/LHAPDF.h>
#include "pdf.h"
#include <iostream>

using namespace constants;

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
  complex <double> ff = 1;
    
  //Gaussian form
  //BLNY (blim = 0.5)
  //g1 = 0.21;  //g2 = 0.68;  //g3 = -0.13;  //Q0 = 3.2;
  //Konychev Nadolsky (blim = 1.5)
  //g1 = 0.20;  //g2 = 0.19;  //g3 = -0.03;  //Q0 = 3.2;  
  //M. Hirai, H. Kawamura, K. Tanaka, https://inspirehep.net/record/1229722 (minimal prescription, MSTW08)
  //g1 = 0.330;  //g2 = 0.066;  //g3 = 0.;  //Q0 = 2.*1.3;
  //g1 = 0.50;   //g2 = 2.0;  //g3 = 0.;  //Q0 = 50;
  if (opts.npff == 0)
    //ff = exp(-(opts.g1 + opts.g2 * log(m/opts.Q0) + opts.g3*log(100*x1*x2)) * pow(b,2));
    //ff = exp(-(opts.g1 + opts.g2 * log(m/opts.Q0) + opts.g3*log(100*x1*x2)) * pow(real(b),2));
    ff = exp(-(opts.g1 + opts.g2 * log(m/opts.Q0) + opts.g3*log(100*x1*x2)) * pow(abs(b),2));

  //Exponential form
  if (opts.npff == 1)
    ff = exp(-opts.e*b);

  //Collins, Rogers, large bT form (https://arxiv.org/abs/1412.3820 Eq.(39), (79))
  //g0 = 0.3; //Q0 = 1.6;
  if (opts.npff == 2)
    {
      double gk = opts.g0 * (1. - exp(- (CF*pdf::alphas(real(pdfevol::qbstar))*pow(real(b),2))/(M_PI*opts.g0*pow(blimit_.rblim_,2))));
      ff = exp(-gk*log(pow(m/opts.Q0,2)));
    }

  //Dispersive approach (DMW) https://arxiv.org/abs/hep-ph/9512336  https://iopscience.iop.org/article/10.1088/1126-6708/1999/07/012/pdf
  if (opts.npff == 3)
    ff = exp(-(opts.a2 * (log(m)+1.) - opts.a2p) * pow(real(b),2)/2. * 2.); //one factor for each leg -> multiply by 2
  
  //Additional flavor dependent form factors
  if (opts.flavour_kt)
    {
      uvff = exp(-opts.g1_uv*pow(b,2)) * ff;
      usff = exp(-opts.g1_us*pow(b,2)) * ff;
      dvff = exp(-opts.g1_dv*pow(b,2)) * ff;
      dsff = exp(-opts.g1_ds*pow(b,2)) * ff;
      ssff = exp(-opts.g1_ss*pow(b,2)) * ff;
      chff = exp(-opts.g1_ch*pow(b,2)) * ff;
      boff = exp(-opts.g1_bo*pow(b,2)) * ff;
      glff = exp(-opts.g1_gl*pow(b,2)) * ff;
    }
  
  //universal form factor
  return ff;
}
