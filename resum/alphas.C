#include "alphas.h"
#include "resconst.h"
#include "resint.h"
#include "coupling.h"
#include "pdf.h"
#include "scales.h"
#include "settings.h"
#include "constants.h"
#include "phasespace.h"
#include "pegasus.h"

using namespace resconst;
using namespace constants;
using namespace resint;

complex <double> alphas::as;
double alphas::as0;
complex <double> alphas::asLO;
complex <double> alphas::asNLO;
complex <double> alphas::asNNLO;
complex <double> alphas::asNNNLO;
complex <double> alphas::asNNNNLO;
complex <double> alphas::as1_1l;
complex <double> alphas::as1_2l;
complex <double> alphas::as1_3l;
complex <double> alphas::as1_4l;
complex <double> alphas::as1_5l;
complex <double> alphas::as2_2l;
complex <double> alphas::as2_3l;
complex <double> alphas::as2_4l;
complex <double> alphas::as2_5l;
complex <double> alphas::as3_3l;
complex <double> alphas::as3_4l;
complex <double> alphas::as3_5l;
complex <double> alphas::as4_4l;
complex <double> alphas::as4_5l;
complex <double> alphas::as5_5l;

//Truncated solution for the running of alphas up to N3LO
complex <double> alphas::calc(complex <double> q, int nloop)
{
  //cout << "alphas evolution, q = " << q << endl;
  
  //double q0 = phasespace::m;
  //double q0 = scales::ren;
  //double q0 = coupling::zmass;
  //as0 = pdf::alphas(q0);
  
  double q0 = scales::res;
  as0 = pdf::alphas(scales::ren);

  /*
  //Runge Kutta evolution
  //  if (opts.asrgkt)
  //    {
  //q = q0+complex <double>(0.,q0);
  //q = complex <double>(0.1,10.);
  rgkt(q, q0, as0);
  cout << "as Runge-Kutta C++:     " << asLO << "  " << asNLO << "  " << asNNLO << "  " << asNNNLO << "  " << asNNNNLO << endl;
      
  //double R20  = pow(q0,2);
  //double R2   = pow(q,2);
  //double AS0 = as0/4./M_PI;
  //int nf = 5;
  //aspar_.naord_ = 0;    asLO = as_(R2, R20, AS0, nf)*4*M_PI;
  //aspar_.naord_ = 1;   asNLO = as_(R2, R20, AS0, nf)*4*M_PI;
  //aspar_.naord_ = 2;  asNNLO = as_(R2, R20, AS0, nf)*4*M_PI;
  //aspar_.naord_ = 3; asNNNLO = as_(R2, R20, AS0, nf)*4*M_PI;
  //aspar_.naord_  = order_.npord_;
  //cout << "as Runge-Kutta fortran: " << asLO << "  " << asNLO << "  " << asNNLO << "  " << asNNNLO << endl;
  //    }
  */
  
  //compute asN4LO numerically
  if (opts.order >= 4)
    rgkt(q, q0, as0);

  complex <double> qlog = log(pow(q/q0,2));

  complex <double> lambda = as0*beta0/M_PI*qlog; //-as0*beta0/M_PI*qlog;
  complex <double> oneplambda = 1. + lambda;
  complex <double> log1plambda = log(oneplambda);
  //complex <double> log1plambda = log1p(lambda); //check if this is faster

  /*
    if (nloop >= 0) //0-loop
    as = as0;
    if (nloop >= 1) //1-loop (LO)
    as = as0/(1.+lambda);
    if (nloop >= 2) //2-loop (NLO)
    as += -(pow(as0,2)*beta1*log1plambda)/(beta0*pow(oneplambda,2)*M_PI);
    if (nloop >= 3) //3-loop (NNLO)
    as += (pow(as0,3)*((pow(beta1,2) - beta0*beta2)*lambda + pow(beta1,2)*(-1.+log1plambda)*log1plambda))/(pow(beta0,2)*pow(oneplambda,3)*pi2);
    if (nloop >= 4) //4-loop (NNNLO)
    as += (pow(as0,4)*(-(lambda*(pow(beta1,3)*lambda - 2.*beta0*beta1*beta2*(oneplambda) + pow(beta0,2)*beta3*(2.+lambda))) + beta1*log1plambda*(-4.*pow(beta1,2)*lambda + 2.*beta0*beta2*(-1.+2.*lambda) + pow(beta1,2)*(5. - 2.*log1plambda)*log1plambda)))/(2.*pow(beta0,3)*pow(oneplambda,4)*pi2*M_PI);
  */
      
  asLO    = as0/(1.+lambda); //Taylor series: asLO = as0*(1.-lambda+pow(lambda,2)-pow(lambda,3)+pow(lambda,4)-pow(lambda,5)+pow(lambda,6));
  asNLO   = asLO   - (pow(as0,2)*beta1*log1plambda)/(beta0*pow(oneplambda,2)*M_PI);
  asNNLO  = asNLO  + (pow(as0,3)*((pow(beta1,2) - beta0*beta2)*lambda + pow(beta1,2)*(-1.+log1plambda)*log1plambda))/(pow(beta0,2)*pow(oneplambda,3)*pi2);
  asNNNLO = asNNLO + (pow(as0,4)*(-(lambda*(pow(beta1,3)*lambda - 2.*beta0*beta1*beta2*(oneplambda) + pow(beta0,2)*beta3*(2.+lambda)))
				  + beta1*log1plambda*(-4.*pow(beta1,2)*lambda + 2.*beta0*beta2*(-1.+2.*lambda)
						       + pow(beta1,2)*(5. - 2.*log1plambda)*log1plambda)))/(2.*pow(beta0,3)*pow(oneplambda,4)*pi2*M_PI);

  //cout << "as analytic             " << asLO << "  " << asNLO << "  " << asNNLO << "  " << asNNNLO << endl;

  //Approximate the exact numeric alphas with higher-order iterative analytic solutions setting beta terms to zero
  //complex <double> asiter;
  //asiter   = as0/(1.+lambda); //Taylor series: asLO = as0*(1.-lambda+pow(lambda,2)-pow(lambda,3)+pow(lambda,4)-pow(lambda,5)+pow(lambda,6));
  //asiter   += - (pow(as0,2)*beta1*log1plambda)/(beta0*pow(oneplambda,2)*M_PI);
  //asiter   += + (pow(as0,3)*((pow(beta1,2))*lambda + pow(beta1,2)*(-1.+log1plambda)*log1plambda))/(pow(beta0,2)*pow(oneplambda,3)*pi2);
  //asiter   += + (pow(as0,4)*(-(lambda*(pow(beta1,3)*lambda))
  //				  + beta1*log1plambda*(-4.*pow(beta1,2)*lambda
  //						       + pow(beta1,2)*(5. - 2.*log1plambda)*log1plambda)))/(2.*pow(beta0,3)*pow(oneplambda,4)*pi2*M_PI);
  //cout << "as iterative             " << asLO << "  " << asiter << endl;
  
  complex <double> dasNLO   = asNLO-asLO;
  complex <double> dasNNLO  = asNNLO-asNLO;
  complex <double> dasNNNLO = asNNNLO-asNNLO;

  //Loop-truncated powers before scale variations
  //To calculate loop-truncate powers write asN^kLO = asLO + dasNLO + dasNNLO + ...
  //Then expand the expression pow(asN^kLO,p), and take the terms which contribute up to k loop,
  //where asLO contributes 1 loop, dasNLO contributes 2 loops, etc...
  //alphas truncated powers of 2
  as2_2l = pow(asLO,2);
  as2_3l = pow(asLO,2)+2.*dasNLO*asLO;
  as2_4l = pow(asLO,2)+pow(dasNLO,2)+2.*asLO*dasNLO+2.*asLO*dasNNLO;
  as2_5l = pow(asLO,2)+pow(dasNLO,2)+2.*asLO*dasNLO+2.*asLO*dasNNLO+2.*dasNLO*dasNNLO+2.*asLO*dasNNNLO;
  //alphas truncated powers of 3
  as3_3l = pow(asLO,3);
  as3_4l = pow(asLO,3)+3.*pow(asLO,2)*dasNLO;
  as3_5l = pow(asLO,3)+3.*pow(asLO,2)*dasNLO+3.*pow(asLO,2)*dasNNLO+3.*asLO*pow(dasNLO,2);
  //alphas truncated powers of 4
  as4_4l = pow(asLO,4);
  as4_5l = pow(asLO,4)+4.*pow(asLO,3)*dasNLO;
  //alphas truncated powers of 5
  as5_5l = pow(asLO,5);

  //QCD coupling scale dependence

  //double LQR = log(pow(q0/scales::ren,2));
  double LQR = log(pow(scales::res/scales::ren,2));

  //if (nloop >= 2) as += -pow(as,2)*beta0/M_PI*LQR;
  //if (nloop >= 3) as += -pow(as,3)*LQR*(beta1 - pow(beta0,2)*LQR)/pi2;
  //if (nloop >= 4) as += -pow(as,4)*LQR*(beta2 - 5./2.*beta0*beta1*LQR + pow(beta0,3)*pow(LQR,2))/(pi2*M_PI);
  //if (nloop >= 5) as += -pow(as,5)*LQR*(beta3 - 3./2.*pow(beta1,2)*LQR - 3.*beta0*beta2*LQR + 13./3.*pow(beta0,2)*beta1*pow(LQR,2) - pow(beta0,4)*pow(LQR,3))/pi4;

  //Runge Kutta evolution
  if (opts.asrgkt)
    rgkt(q, q0, as0);
  
  //Truncate alphas at the exact order
  asNLO    += -as2_2l * beta0/M_PI*LQR;
  asNNLO   += -as2_3l * beta0/M_PI*LQR;
  asNNLO   += -as3_3l * LQR*(beta1 - pow(beta0,2)*LQR)/pi2;
  asNNNLO  += -as2_4l * beta0/M_PI*LQR;
  asNNNLO  += -as3_4l * LQR*(beta1 - pow(beta0,2)*LQR)/pi2;
  asNNNLO  += -as4_4l * LQR*(beta2 - 5./2.*beta0*beta1*LQR + pow(beta0,3)*pow(LQR,2))/(pi2*M_PI);
  asNNNNLO += -as2_5l * beta0/M_PI*LQR;
  asNNNNLO += -as3_5l * LQR*(beta1 - pow(beta0,2)*LQR)/pi2;
  asNNNNLO += -as4_5l * LQR*(beta2 - 5./2.*beta0*beta1*LQR + pow(beta0,3)*pow(LQR,2))/(pi2*M_PI);
  asNNNNLO += -as5_5l * LQR*(beta3 - 3./2.*pow(beta1,2)*LQR - 3.*beta0*beta2*LQR + 13./3.*pow(beta0,2)*beta1*pow(LQR,2) - pow(beta0,4)*pow(LQR,3))/pi4;
  
//  //Do not truncate alphas
//  asNLO   += -pow(asNLO,2)   * beta0/M_PI*LQR;
//  asNNLO  += -pow(asNNLO,2)  * beta0/M_PI*LQR;
//  asNNLO  += -pow(asNNLO,3)  * LQR*(beta1 - pow(beta0,2)*LQR)/pi2;
//  asNNNLO += -pow(asNNNLO,2) * beta0/M_PI*LQR;
//  asNNNLO += -pow(asNNNLO,3) * LQR*(beta1 - pow(beta0,2)*LQR)/pi2;
//  asNNNLO += -pow(asNNNLO,4) * LQR*(beta2 - 5./2.*beta0*beta1*LQR + pow(beta0,3)*pow(LQR,2))/(pi2*M_PI);

//  //Do not evolve alphas  
//  asNLO   += -pow(as0,2) * beta0/M_PI*LQR;
//  asNNLO  += -pow(as0,2) * beta0/M_PI*LQR;
//  asNNLO  += -pow(as0,3) * LQR*(beta1 - pow(beta0,2)*LQR)/pi2;
//  asNNNLO += -pow(as0,2) * beta0/M_PI*LQR;
//  asNNNLO += -pow(as0,3) * LQR*(beta1 - pow(beta0,2)*LQR)/pi2;
//  asNNNLO += -pow(as0,4) * LQR*(beta2 - 5./2.*beta0*beta1*LQR + pow(beta0,3)*pow(LQR,2))/(pi2*M_PI);
  
  
  //1 over pi normalisation
  as0      /= M_PI;
  asLO     /= M_PI;
  asNLO    /= M_PI;
  asNNLO   /= M_PI;
  asNNNLO  /= M_PI;
  asNNNNLO /= M_PI;
    
  //Recompute truncated powers including scale variations and 1 over pi normalisation
  dasNLO   = asNLO-asLO;
  dasNNLO  = asNNLO-asNLO;
  dasNNNLO = asNNNLO-asNNLO;

  as1_1l = asLO;
  as1_2l = asNLO;
  as1_3l = asNNLO;
  as1_4l = asNNNLO;
  as1_5l = asNNNNLO;
  //alphas truncated powers of 2
  as2_2l = pow(asLO,2);
  as2_3l = pow(asLO,2)+2.*dasNLO*asLO;
  as2_4l = pow(asLO,2)+pow(dasNLO,2)+2.*asLO*dasNLO+2.*asLO*dasNNLO;
  as2_5l = pow(asLO,2)+pow(dasNLO,2)+2.*asLO*dasNLO+2.*asLO*dasNNLO+2.*dasNLO*dasNNLO+2.*asLO*dasNNNLO;
  //alphas truncated powers of 3
  as3_3l = pow(asLO,3);
  as3_4l = pow(asLO,3)+3.*pow(asLO,2)*dasNLO;
  as3_5l = pow(asLO,3)+3.*pow(asLO,2)*dasNLO+3.*pow(asLO,2)*dasNNLO+3.*asLO*pow(dasNLO,2);
  //alphas truncated powers of 4
  as4_4l = pow(asLO,4);
  as4_5l = pow(asLO,4)+4.*pow(asLO,3)*dasNLO;
  //alphas truncated powers of 5
  as5_5l = pow(asLO,5);

  if      (nloop == 0) as = as0;
  else if (nloop == 1) as = asLO;
  else if (nloop == 2) as = asNLO;
  else if (nloop == 3) as = asNNLO;
  else if (nloop == 4) as = asNNNLO;
  else if (nloop == 5) as = asNNNNLO;

  //do not truncate alphas
  if (opts.asrgkt)
    {
      as1_1l = as1_2l = as1_3l = as1_4l = as1_5l = as;
      as2_2l = as2_3l = as2_4l = as2_5l = pow(as,2);
      as3_3l = as3_4l = as3_5l =  pow(as,3);
      as4_4l = as4_5l =  pow(as,4);
      as5_5l =  pow(as,5);
    }
 
  return as;
}

////interface to the VFN alphas of pegasus
//double alphas::pegasus(double q)
//{
//  double ASI;
//  int NF;  
//  double R2 = pow(q,2);
//  double ASF = 0.;
//  if (R2 > asfthr_.m2t_)
//    {
//      NF = 6;
//      double R2T = asfthr_.m2t_;
//      ASI = asfthr_.ast_;
//      ASF = as_(R2, R2T, asfthr_.ast_, NF);
//    }
//  else if (M2 > asfthr_.m2b_)
//    {
//      NF = 5;
//      double R2B = asfthr_.m2b_;
//      ASI = asfthr_.asb_;
//      ASF =  as_(R2, R2B, asfthr_.asb_, NF);
//    }
//  else if (M2 > asfthr_.m2c_)
//    {
//      NF = 4;
//      double R2C = asfthr_.m2c_;
//      ASI = asfthr_.asc_;
//      ASF =  as_(R2, R2C, asfthr_.asc_, NF);
//    }
//  else
//    {
//      NF = 3;
//      double R20 = asinp_.m20_;
//      ASI = asinp_.as0_;
//      ASF =  as_(R2, R20, asinp_.as0_, NF);
//    }
//  
//  //ASF = LHAPDF::alphasPDF(sqrt(M2)) / (4.* M_PI);
//  return ASF*4.*M_PI;
//}

// code for Runge-Kutta solutions from Pegasus asrgkt.f extended to five-loop with beta4
//
// ..File: asrgkt.f      (requires a previous call of BETAFCT)
//
//
// ..The running coupling of QCD,  
//
//         AS  =  a_s  =  alpha_s(mu_r^2)/(4 pi),
//
//    obtained by integrating the evolution equation for a fixed number
//    of massless flavours  NF.  Except at leading order (LO),  AS  is 
//    obtained using a fourth-order Runge-Kutta integration. 
//
// ..The initial and final scales  R20  and  R2,  the value  AS0  at
//    R20, and  NF  are passed as function arguments.  The coefficients 
//    of the beta function up to  a_s^5 (N^3LO)  are provided by the 
//    common-block  BETA.  The order of the expansion  NAORD  (defined 
//    as the 'n' in N^nLO) and the number of steps  NASTPS  for the 
//    integration beyond LO are given by the common-block  ASPAR.
//
// =====================================================================

void alphas::rgkt(complex <double> q, double q0, double aass0)
{
  complex <double> lrrat = log(pow(q/q0,2));
  complex <double> dlr = lrrat / double(nastps);

  asLO     = aass0/M_PI / (1.+ beta0*aass0/M_PI*lrrat);
  asNLO    = aass0/M_PI;
  asNNLO   = aass0/M_PI;
  asNNNLO  = aass0/M_PI;
  asNNNNLO = aass0/M_PI;

  //Solve the evolution equation with a fourth-order Runge-Kutta.
  complex <double> xk0, xk1, xk2, xk3;  
  for (int k1 = 1; k1 <= nastps; k1++)
    {
      //2-loop
      xk0 = dlr * fbeta1(asNLO);
      xk1 = dlr * fbeta1(asNLO+0.5*xk0);
      xk2 = dlr * fbeta1(asNLO+0.5*xk1);
      xk3 = dlr * fbeta1(asNLO+xk2);
      asNLO += (xk0+2.*xk1+2.*xk2+xk3)/6.;
  
      //3-loop
      xk0 = dlr * fbeta2(asNNLO);
      xk1 = dlr * fbeta2(asNNLO + 0.5 * xk0);
      xk2 = dlr * fbeta2(asNNLO + 0.5 * xk1);
      xk3 = dlr * fbeta2(asNNLO + xk2);
      asNNLO = asNNLO + (xk0 + 2.* xk1 + 2.* xk2 + xk3)/6.;

      //4-loop
      xk0 = dlr * fbeta3(asNNNLO);
      xk1 = dlr * fbeta3(asNNNLO + 0.5 * xk0);
      xk2 = dlr * fbeta3(asNNNLO + 0.5 * xk1);
      xk3 = dlr * fbeta3(asNNNLO + xk2);
      asNNNLO = asNNNLO + (xk0 + 2.* xk1 + 2.* xk2 + xk3)/6.;

      //5-loop
      xk0 = dlr * fbeta4(asNNNNLO);
      xk1 = dlr * fbeta4(asNNNNLO + 0.5 * xk0);
      xk2 = dlr * fbeta4(asNNNNLO + 0.5 * xk1);
      xk3 = dlr * fbeta4(asNNNNLO + xk2);
      asNNNNLO = asNNNNLO + (xk0 + 2.* xk1 + 2.* xk2 + xk3)/6.;
    }

  asLO     *= M_PI;
  asNLO    *= M_PI;
  asNNLO   *= M_PI;
  asNNNLO  *= M_PI;
  asNNNNLO *= M_PI;
}
