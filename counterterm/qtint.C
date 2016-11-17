#include "qtint.h"

#include "mesq.h"
#include "settings.h"
#include "phasespace.h"
#include "omegaintegr.h"
#include "gaussrules.h"
#include "switch.h"
#include "dyres_interface.h"

int one = 1;
int two = 2;
int three = 3;
int four = 4;

int qtintervals = 1; //--> make a setting
int qtrule = 64;     //--> make a setting

double qtint::LL1_mesqij[12];
double qtint::LL2_mesqij[12];
double qtint::LL3_mesqij[12];
double qtint::LL4_mesqij[12];

//perform qt integration for the counterterm
void qtint::calc(double m, double qtmin, double qtmax, int mode)
{
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

      phasespace::set_qt(qt);
      omegaintegr::genV4p();
      double cthmom0,cthmom1,cthmom2;
      omegaintegr::cthmoments(cthmom0, cthmom1, cthmom2);
      mesq::setmesq(cthmom0, cthmom1, cthmom2);
      
      double qt2 = qt*qt;
      double swtch = switching::swtch(qt, m);
      //xmio is used in besselkfast for Itilde
      xmio_.xmio_ = sqrt(qt2/(q2/pow(a_param_.a_param_,2)));
      
      double LL1 = itilde_(one)/pow(q2,2)*pow(a_param_.a_param_,2)*swtch;
      double LL2 = itilde_(two)/pow(q2,2)*pow(a_param_.a_param_,2)*swtch;
      double LL3 = itilde_(three)/pow(q2,2)*pow(a_param_.a_param_,2)*swtch;
      double LL4 = itilde_(four)/pow(q2,2)*pow(a_param_.a_param_,2)*swtch;
      for (int sp = 0; sp < mesq::totpch; sp++)
	{
	  LL1_mesqij[sp] = LL1*real(mesq::mesqij[sp]);
	  LL2_mesqij[sp] = LL2*real(mesq::mesqij[sp]);
	  LL3_mesqij[sp] = LL3*real(mesq::mesqij[sp]);
	  LL4_mesqij[sp] = LL4*real(mesq::mesqij[sp]);
	}
      return;
    }
  
  
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
      
  double qtmin2 = pow(qtmin,2);
  double qtmax2 = pow(qtmax,2);
  double tiny = 1e-5;
  double qta = 1./(1.+log(qtmax2/tiny));
  double qtb = 1./(1.+log(qtmin2/tiny));
  double min = 0;
  double max = 1;
  for (int i = 0; i < qtintervals; i++)
    {
      double xa = min+(max-min)*i/qtintervals;
      double xb = min+(max-min)*(i+1)/qtintervals;
      double xc = 0.5*(xa+xb);
      double xm = 0.5*(xb-xa);
      for (int j = 0; j < qtrule; j++)
	{
	  double x = xc+xm*gr::xxx[qtrule-1][j];
	  double qtx = qta + (qtb-qta) * x;
	  double qt2 = tiny*exp(1./qtx - 1.);
	  double jac = (qtb-qta)*qt2/pow(qtx,2);
	  double w = gr::www[qtrule-1][j]*xm*jac;

	  double qt = sqrt(qt2);

	  double swtch = switching::swtch(qt, m);
	  //if(switch.le.0.01d0) cycle ! do not apply this cut to avoid discontinuities. Instead the phase space is limited to qt and m switching limits

	  //xmio is used in besselkfast for Itilde
	  xmio_.xmio_ = sqrt(qt2/(q2/pow(a_param_.a_param_,2)));

	  double LL1 = itilde_(one)/pow(q2,2)*pow(a_param_.a_param_,2)*w*swtch;
	  double LL2 = itilde_(two)/pow(q2,2)*pow(a_param_.a_param_,2)*w*swtch;
	  double LL3 = itilde_(three)/pow(q2,2)*pow(a_param_.a_param_,2)*w*swtch;
	  double LL4 = itilde_(four)/pow(q2,2)*pow(a_param_.a_param_,2)*w*swtch;

	  if (!opts.fixedorder && opts.makecuts)
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
