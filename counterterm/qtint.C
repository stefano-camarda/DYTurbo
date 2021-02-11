#include "qtint.h"

#include "scales.h"
#include "mesq.h"
#include "settings.h"
#include "phasespace.h"
#include "omegaintegr.h"
#include "gaussrules.h"
#include "switch.h"
#include "itilde.h"
#include "dyres_interface.h"
#include "numbers.h"
#include <iostream>

double qtint::LL1_mesqij[12];
double qtint::LL2_mesqij[12];
double qtint::LL3_mesqij[12];
double qtint::LL4_mesqij[12];

//perform qt integration for the counterterm
void qtint::calc(double m, double qtmin, double qtmax, int mode)
{
  //Call scales::set in order to set a_param (this is a duplicate call to scales::set)
  scales::set(m);
  scales::dyres(m);
    //a_param_.a_param_ = m/scales::res;

  for (int sp = 0; sp < mesq::totpch; sp++)
    {
      LL1_mesqij[sp] = 0.;
      LL2_mesqij[sp] = 0.;
      LL3_mesqij[sp] = 0.;
      LL4_mesqij[sp] = 0.;
    }

  double q2 = m*m;
  mesq::setpropagators(m);
  if (mode == 0 || mode == 1)
    {
      double qt = qtmin;

      /*
      //limit the integration to qt = m
      double qtp = qt;
      qt = qtp/sqrt(1-pow(qtp/m,2));
      if (qtp >= m*0.999)
        LL1=LL2=LL3=LL4=0;
      */
      
      phasespace::set_qt(qt);
      omegaintegr::genV4p();
      double cthmom0,cthmom1,cthmom2;
      omegaintegr::cthmoments(cthmom0, cthmom1, cthmom2);
      mesq::setmesq(cthmom0, cthmom1, cthmom2);
      
      double qt2 = qt*qt;
      double swtch = switching::swtch(qt, m);
      //xmio is used in besselkfast for Itilde
      xmio_.xmio_ = sqrt(qt2/(q2/pow(a_param_.a_param_,2)));

      /*
      cout << endl;
      cout << " qt " << qt  << " m " << m << endl;
      cout << " fortran itilde(1) " << itilde_(one) << endl;
      itilde::calc(qt,m/a_param_.a_param_,1);
      //cout << " C++ itilde(1) " << itilde::calc(qt,m/a_param_.a_param_,1)
      cout << " fortran itilde(2) " << itilde_(two) << endl;
      itilde::calc(qt,m/a_param_.a_param_,2);
      //cout << " C++ itilde(2) " << itilde::calc(qt,m/a_param_.a_param_,2)
      cout << " fortran itilde(3) " << itilde_(three) << endl;
      itilde::calc(qt,m/a_param_.a_param_,3);
      //<< " C++ itilde(3) " << itilde::calc(qt,m/a_param_.a_param_,3)
      cout << " fortran itilde(4) " << itilde_(four) << endl;
      itilde::calc(qt,m/a_param_.a_param_,4);
      //<< " C++ itilde(4) " << itilde::calc(qt,m/a_param_.a_param_,4)
      */

      double LL1,LL2,LL3,LL4;
      
      //LL1 = itilde_(one)  /pow(q2,2)*pow(a_param_.a_param_,2)*swtch;
      //LL2 = itilde_(two)  /pow(q2,2)*pow(a_param_.a_param_,2)*swtch;
      //LL3 = itilde_(three)/pow(q2,2)*pow(a_param_.a_param_,2)*swtch;
      //LL4 = itilde_(four) /pow(q2,2)*pow(a_param_.a_param_,2)*swtch;
      
      double Q = scales::res;
      double x = qt/Q;
      LL1 = itilde::besselk(x,1)/pow(Q,2)/q2*swtch;
      LL2 = itilde::besselk(x,2)/pow(Q,2)/q2*swtch;
      LL3 = itilde::besselk(x,3)/pow(Q,2)/q2*swtch;
      LL4 = itilde::besselk(x,4)/pow(Q,2)/q2*swtch;
      
      for (int sp = 0; sp < mesq::totpch; sp++)
	{
	  LL1_mesqij[sp] = LL1*real(mesq::mesqij[sp]);
	  LL2_mesqij[sp] = LL2*real(mesq::mesqij[sp]);
	  LL3_mesqij[sp] = LL3*real(mesq::mesqij[sp]);
	  LL4_mesqij[sp] = LL4*real(mesq::mesqij[sp]);
	}
      return;
    }

  // mode == 2: qt-integrated
  
  //For fixed order predictions the counterterm is always evaluated at qt=0
  //Also when no cuts on the leptons are required, the cthmom are independent of pt, and so are the mesqij
  //Hence calculate moments and store amplitudes once for all for qt points
  if (opts.fixedorder || !opts.makecuts)
    {
      phasespace::set_qt(0.);
      omegaintegr::genV4p();
      double cthmom0,cthmom1,cthmom2;
      omegaintegr::cthmoments(cthmom0, cthmom1, cthmom2);
      mesq::setmesq(cthmom0, cthmom1, cthmom2);
    }

  //In the fixed order case we just need to calculate the integral over qt from qtcut to infinity
  //of Eq.(112) in https://arxiv.org/pdf/hep-ph/0508068.pdf:
  //int_qtcut^inf dqt Itilde_n(qt/Q) = Q^2*qtcut/2 int_0^inf db J1(b*qtcut) ln^n(Q^2 b^2 / b0^2 + 1)
  if (opts.fixedorder)
    {
      double LL1int, LL2int, LL3int, LL4int;
      //LL1int = -2.*qtmin*itilde::integral(qtmin,m/a_param_.a_param_,1)/pow(q2,2)*pow(a_param_.a_param_,2);
      //LL2int = -2.*qtmin*itilde::integral(qtmin,m/a_param_.a_param_,2)/pow(q2,2)*pow(a_param_.a_param_,2);
      //LL3int = -2.*qtmin*itilde::integral(qtmin,m/a_param_.a_param_,3)/pow(q2,2)*pow(a_param_.a_param_,2);
      //LL4int = -2.*qtmin*itilde::integral(qtmin,m/a_param_.a_param_,4)/pow(q2,2)*pow(a_param_.a_param_,2);

      double Q = scales::res;
      LL1int = -2.*qtmin*itilde::integral(qtmin,Q,1)/pow(Q,2)/q2;
      LL2int = -2.*qtmin*itilde::integral(qtmin,Q,2)/pow(Q,2)/q2;
      LL3int = -2.*qtmin*itilde::integral(qtmin,Q,3)/pow(Q,2)/q2;
      LL4int = -2.*qtmin*itilde::integral(qtmin,Q,4)/pow(Q,2)/q2;
      
      //loop on born subprocesses, i.e. born incoming partons ij
      for (int sp = 0; sp < mesq::totpch; sp++)
	{
	  LL1_mesqij[sp] += LL1int*real(mesq::mesqij[sp]);
	  LL2_mesqij[sp] += LL2int*real(mesq::mesqij[sp]);
	  LL3_mesqij[sp] += LL3int*real(mesq::mesqij[sp]);
	  LL4_mesqij[sp] += LL4int*real(mesq::mesqij[sp]);
	}
      return;
    }
  /*  //In the resummation case the integral can be performed if qtmax < dampk*m (the result is unstable for qt >~ 500)
  else if (!opts.makecuts && qtmax < m*opts.dampk && qtmax < m*5)
    {
      double LL1int = 2.*(qtmax*itilde::integral(qtmax,m/a_param_.a_param_,1) - qtmin*itilde::integral(qtmin,m/a_param_.a_param_,1))/pow(q2,2)*pow(a_param_.a_param_,2);
      double LL2int = 2.*(qtmax*itilde::integral(qtmax,m/a_param_.a_param_,2) - qtmin*itilde::integral(qtmin,m/a_param_.a_param_,2))/pow(q2,2)*pow(a_param_.a_param_,2);
      double LL3int = 2.*(qtmax*itilde::integral(qtmax,m/a_param_.a_param_,3) - qtmin*itilde::integral(qtmin,m/a_param_.a_param_,3))/pow(q2,2)*pow(a_param_.a_param_,2);
      double LL4int = 2.*(qtmax*itilde::integral(qtmax,m/a_param_.a_param_,4) - qtmin*itilde::integral(qtmin,m/a_param_.a_param_,4))/pow(q2,2)*pow(a_param_.a_param_,2);

      //loop on born subprocesses, i.e. born incoming partons ij
      for (int sp = 0; sp < mesq::totpch; sp++)
	{
	  LL1_mesqij[sp] += LL1int*real(mesq::mesqij[sp]);
	  LL2_mesqij[sp] += LL2int*real(mesq::mesqij[sp]);
	  LL3_mesqij[sp] += LL3int*real(mesq::mesqij[sp]);
	  LL4_mesqij[sp] += LL4int*real(mesq::mesqij[sp]);
	}
      return;
    }
  */


  //modlog with p > 1
  if (!opts.makecuts && opts.p > 1)
    {
      double LL1int = 2.*(qtmax*itilde::integral(qtmax,m/a_param_.a_param_,1) - qtmin*itilde::integral(qtmin,m/a_param_.a_param_,1))/pow(q2,2)*pow(a_param_.a_param_,2);
      double LL2int = 2.*(qtmax*itilde::integral(qtmax,m/a_param_.a_param_,2) - qtmin*itilde::integral(qtmin,m/a_param_.a_param_,2))/pow(q2,2)*pow(a_param_.a_param_,2);
      double LL3int = 2.*(qtmax*itilde::integral(qtmax,m/a_param_.a_param_,3) - qtmin*itilde::integral(qtmin,m/a_param_.a_param_,3))/pow(q2,2)*pow(a_param_.a_param_,2);
      double LL4int = 2.*(qtmax*itilde::integral(qtmax,m/a_param_.a_param_,4) - qtmin*itilde::integral(qtmin,m/a_param_.a_param_,4))/pow(q2,2)*pow(a_param_.a_param_,2);

      //loop on born subprocesses, i.e. born incoming partons ij
      for (int sp = 0; sp < mesq::totpch; sp++)
	{
	  LL1_mesqij[sp] += LL1int*real(mesq::mesqij[sp]);
	  LL2_mesqij[sp] += LL2int*real(mesq::mesqij[sp]);
	  LL3_mesqij[sp] += LL3int*real(mesq::mesqij[sp]);
	  LL4_mesqij[sp] += LL4int*real(mesq::mesqij[sp]);
	}
      return;
    }
  
  //In the general case (resummation, either with cuts or at high pT) use a gaussian quadrature integration

  /*
  //limit the integration to qt = m
  double Q = m;//2.*scales::res; //m;
  if (qtmin >= Q*0.999999)
    return;
  else
    //qtmin = qtmin/sqrt(1-pow(qtmin/Q,2));
    qtmin = Q*atanh(qtmin/Q);
  if (qtmax >= Q*0.999999)
    qtmax = 1e10; //ctlimit
  else
    //qtmax = qtmax/sqrt(1-pow(qtmax/Q,2));
    qtmax = Q*atanh(qtmax/Q);
  */
  
  double qtmin2 = pow(qtmin,2);
  double qtmax2 = pow(qtmax,2);
  double tiny = 1e-5; //be carefull, if qtmin^2 < tiny the phase space generation is screwed up
  double qta = 1./(1.+log(qtmax2/tiny));
  double qtb = 1./(1.+log(qtmin2/tiny));
  double min = 0;
  double max = 1;
  for (int i = 0; i < opts.qtintervals; i++)
    {
      double xa = min+(max-min)*i/opts.qtintervals;
      double xb = min+(max-min)*(i+1)/opts.qtintervals;
      double xc = 0.5*(xa+xb);
      double xm = 0.5*(xb-xa);
      for (int j = 0; j < opts.qtrule; j++)
	{
	  double x = xc+xm*gr::xxx[opts.qtrule-1][j];
	  double qtx = qta + (qtb-qta) * x;
	  double qt2 = tiny*exp(1./qtx - 1.);
	  double jac = (qtb-qta)*qt2/pow(qtx,2);
	  double w = gr::www[opts.qtrule-1][j]*xm*jac;

	  double qt = sqrt(qt2);
	  //Here there should be a jacobian: jac=jac/qt/2.;
	  //However since the original formulation of the integral is in dqt2,
	  //the integral in dqt gets a factor 2*qt from the change of variable which cancels exactly the Jacobian

	  double swtch = switching::swtch(qt, m);
	  //if(switch.le.0.01d0) cycle ! do not apply this cut to avoid discontinuities. Instead the phase space is limited to qt and m switching limits

	  //xmio is used in besselkfast for Itilde
	  xmio_.xmio_ = sqrt(qt2/(q2/pow(a_param_.a_param_,2)));

	  double LL1 = itilde_(one)/pow(q2,2)*pow(a_param_.a_param_,2)*w*swtch;
	  double LL2 = itilde_(two)/pow(q2,2)*pow(a_param_.a_param_,2)*w*swtch;
	  double LL3 = itilde_(three)/pow(q2,2)*pow(a_param_.a_param_,2)*w*swtch;
	  double LL4 = itilde_(four)/pow(q2,2)*pow(a_param_.a_param_,2)*w*swtch;

	  if (opts.makecuts)
	    {
	      phasespace::set_qt(qt);
	      omegaintegr::genV4p();
	      double cthmom0,cthmom1,cthmom2;
	      omegaintegr::cthmoments(cthmom0, cthmom1, cthmom2);
	      mesq::setmesq(cthmom0, cthmom1, cthmom2);
	    }

	  //loop on born subprocesses, i.e. born incoming partons ij
	  for (int sp = 0; sp < mesq::totpch; sp++)
	    {
	      LL1_mesqij[sp] += LL1*real(mesq::mesqij[sp]);
	      LL2_mesqij[sp] += LL2*real(mesq::mesqij[sp]);
	      LL3_mesqij[sp] += LL3*real(mesq::mesqij[sp]);
	      LL4_mesqij[sp] += LL4*real(mesq::mesqij[sp]);
	    }
	}
    }

}
