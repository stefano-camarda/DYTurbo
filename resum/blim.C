#include "blim.h"
#include "resint.h"
#include "resconst.h"
#include "scales.h"
#include "settings.h"

double blim::pdf;
double blim::sudakov;
double blim::expc;

void blim::set()
{
  double Q = scales::res;
  //if (opts.mufevol)
  //Q = scales::fac;
  
  pdf     = calc(opts.blim_pdf,Q);
  sudakov = calc(opts.blim_sudakov,scales::res);
  expc    = calc(opts.blim_expc,scales::res);
}

double blim::calc(double bl, double Q)
{
  double bmax;
  //dynamic blim: set bmax to b0/lambda_QCD, where lambda_QCD = Q/(exp(1./(2.*as*beta0))) is the Landau pole
  //more precisely, the pole is at lambda = 1 <==> beta0*as*blog = 1
  if (bl < 0)
    {
      if (!opts.modlog)
	bmax = resconst::b0/Q * exp(1./(2.*resint::aass*resconst::beta0));
      else if (opts.p == 1)
	//bmax = resconst::b0/Q * (exp(1./(2.*resint::aass*resconst::beta0))-1.);
	//I think the correct formula is:
	bmax = resconst::b0/Q * sqrt(exp(1./(resint::aass*resconst::beta0))-1.);
      else
	bmax = resconst::b0/Q * pow(sqrt(exp(opts.p/(resint::aass*resconst::beta0))-1.),1./opts.p);

      //scale by factor (-bl)
      bmax = bmax/(bl*(-1));
    }

  //fixed blim
  if (bl > 0)
    bmax = bl;

  return bmax;

  /*
  //In DYRes blim was reduced when muren/q > sqrt(2.2)
  //if (opts.order == 2 && q2 > 2.2*muren2) //avoid large values of b in f2(y) when q2>mur2
  if (q2 > 2.2*muren2) //avoid large values of b in f2(y) when q2>mur2 
    blim = a_param_.b0p_*(1./q)*exp(1./(2.*aass*resconst::beta0))*sqrt(muren2/q2);
  else
    blim = a_param_.b0p_*(1./q)*exp(1./(2.*aass*resconst::beta0)); // avoid Landau pole
  */

  /*
  //Old code (I do not think there is any other divergence than that corresponding to beta0*alphas*L = 1)
  double lambdaqcd         = muren/(exp(1./(2.*aass*resconst::beta0))); //--> corresponds to a divergence in alphas
  double lambdasudakov     = mures/(exp(1./(2.*aass*resconst::beta0))); //--> correspond to a divergence in the Sudakov
    
  //should blim depend or not on a_param? In principle yes because PDFs are evolved to b0/bstar*a
  blim = min(resconst::b0/lambdaqcd,resconst::b0/lambdasudakov);
  //cout << _m << "  " << resconst::b0/lambdaqcd << "  " << resconst::b0/lambdasudakov << "  " << a_param_.b0p_*(1./q)*exp(1./(2.*aass*resconst::beta0)) << endl;
  */
}
