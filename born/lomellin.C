#include "lomellin.h"
#include "mesq.h"
#include "settings.h"
#include "muf.h"
#include "scales.h"
#include "mellinint.h"
#include "ccoeff.h"
#include "pmom.h"
#include "pegasus.h"
#include "evolnative.h"
#include "pdfevol.h"
#include "rapint.h"
#include "resint.h"
#include "pdf.h"
#include "expc.h"
#include "hcoeff.h"

#include "phasespace.h"
// 
// #include <iostream>
// #include <string.h>

using namespace std;

void lomellin::lint(double costh, double m, double y, int mode, double f[2])
{

  f[0] = 0.;
  f[1] = 0.;

  //double m2 = m*m;
  //double exppy = exp(y);
  //double expmy = 1./exppy;
  //double tau = sqrt(m2/pow(opts.sroot,2));
  
  //mesq::setpropagators(m);

  // Kinematical limit in rapidity
  double q2 = m*m;
  double x = q2/pow(opts.sroot,2);
  double ax = log(x);
  double ylim=-0.5*ax;
  if (abs(y) > ylim)
    return;
  

  /*
  double cthmom0, cthmom1, cthmom2;
  if (mode == 0)
    {
      //apply here lepton cuts
      cthmom0 = 1.;
      cthmom1 = costh;
      cthmom2 = pow(costh,2);
    }
  else if (mode == 1)
    {
      phasespace::set_mqtyphi(m, 0., y);
      omegaintegr::genV4p();
      omegaintegr::cthmoments(cthmom0,cthmom1,cthmom2);
    }

  mesq::setmesq(cthmom0, cthmom1, cthmom2);
  */

  //cout << cthmom0 << "  " << cthmom1 << "  " << cthmom2 << endl;

  //Set scales
  scales::set(m);
  //scales::mcfm();
  //double muf = scales::fac;
  //double mur = scales::ren;
  resint::LQ = 2.*log(m/scales::res);
  resint::LF = 2.*log(m/scales::fac);
  resint::LR = 2.*log(m/scales::ren);

  if (opts.melup == 2)
    {
      //anomalous::calc(); //--> should switch to ccoeff for all C coefficients!!!

      mellinint::allocate();
      ccoeff::allocate();
      //pegasus::allocate();
      pmom::allocate();
      mellinint::update();

      /*
      mellinint::updategauss();
      if (opts.mellin1d)
	ccoeff::calc1d();
      else
	ccoeff::calc2d();

      pegasus::calc_mellin();
      pmom::calc();
      */
    }

  evolnative::allocate();
  evolnative::update();
  pdfevol::allocate_fx();

  //retrieve PDFs in Mellin space at the factorisation scale
  for (int n = 0; n < mellinint::mdim; n++)
    if (opts.mellin1d)
      pdfevol::retrievemuf_1d(n);
    else
      pdfevol::retrievemuf_2d(n);
  //evolnative::evolve();
 
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
	{
	  //rapint::cache(ymn, ymx);
	  rapint::integrate(ymn,ymx,phasespace::m);
	}
    }

  mesq::allocate();
  mesq::setmesq_expy(mode, m, costh, y);
  
  
  if (opts.alphaslha)
    resint::aass = pdf::alphas(scales::ren)/M_PI;
  else
    resint::aass = pdf::rgktalphas(scales::ren)/M_PI;

  hcoeff::allocate();
  hcoeff::reset();
  hcoeff::calc();

  if (opts.mellin1d)
    {
      muf::allocate();
      muf::reset();
    }  
  mellinint::reset();
  
  expc::allocate();
  expc::reset();

  complex <double> xmsq = 0.;
  if (opts.mellin1d)
    xmsq = mellinint::calc1d();
  else
    xmsq = mellinint::calc2d();

  xmsq = xmsq * 3./8. /2./M_PI;

  //m^2/s factor
  //xmsq = xmsq * phasespace::m2/pow(opts.sroot,2);
  xmsq = xmsq /pow(opts.sroot,2);

  //phiV integration  
  xmsq = xmsq * 2.*M_PI;
  
  expc::free();
  hcoeff::free();
  if (opts.mellin1d)
    muf::free();
  mesq::free();

  if (opts.melup == 2)
    {
      mellinint::free();
      ccoeff::free();
      //pegasus::free();
      pmom::free();
    }
  //if (opts.fmufac > 0)
    evolnative::free();

    //  pdfevol::free();
    
  pdfevol::free_fx();
  if (!opts.mellin1d && (mode == 2 || mode == 3))
    rapint::free;

  f[0] = real(xmsq);
}
