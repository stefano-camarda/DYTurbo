#include "propagator.h"

#include "coupling.h"
#include "settings.h"

double prop::m;
double prop::w;
double prop::m2;
double prop::w2;

double prop::W;
double prop::Z;
double prop::G;
double prop::ZG;

//fortran interface
fcomplex bosprop_(double &q2)
{
  return fcx(q2*prop::bos(q2));
}


void prop::init()
{
  if (opts.nproc == 3)
    {
      m = coupling::zmass;
      w = coupling::zwidth;
      m2 = pow(coupling::zmass,2);
      w2 = pow(coupling::zwidth,2);
    }
  else
    {
      m = coupling::wmass;
      w = coupling::wwidth;
      m2 = pow(coupling::wmass,2);
      w2 = pow(coupling::wwidth,2);
    }
}

//Complex vector-boson propagator
complex <double> prop::bos(double q2)
{
  //running width correction
  double runw = 1.;
  if (opts.runningwidth)
    runw = q2/m2;

  return 1./complex <double> (q2-m2, m*w*runw);
}

//photon propagator
double prop::gamma(double q2)
{
  return 1./q2;
}

//Propagator's contribution to squared amplitudes

//Complex values algebra
void prop::set(double q2)
{
  Z = 0;
  W = 0;
  G = 0;
  ZG = 0;

  if (opts.nproc == 3)
    {
      Z = q2*pow(abs(bos(q2)),2);
      if (opts.useGamma)
	{
	  G = 1./q2;           // --> G = q2*pow(gamma(q2),2);
	  ZG = bos(q2).real(); //no factor of 2 --> ZG = 2*q2*gamma(q2)*bos(q2).real();
	}
    }
  else
    W = q2*pow(abs(bos(q2)),2);

  //
}


//Real values algebra (fully equivalent to prop::set(q2))
void prop::set_real(double q2)
{
  Z = 0;
  W = 0;
  G = 0;
  ZG = 0;

  double mZ2 = pow(coupling::zmass, 2);
  double wZ2 = pow(coupling::zwidth,2);
  double mW2 = pow(coupling::wmass, 2);
  double wW2 = pow(coupling::wwidth,2);

  //running width correction
  double runw = 1.;
  if (opts.runningwidth)
    {
      if (opts.nproc == 3)
	runw = pow(q2/mZ2,2);
      else
	runw = pow(q2/mW2,2);
    }
  
  if (opts.nproc == 3)
    {
      Z = q2/(pow(q2-mZ2,2)+mZ2*wZ2*runw);
      if (opts.useGamma)
	{
	  G = 1./q2;
	  ZG = (q2-mZ2)/(pow(q2-mZ2,2)+mZ2*wZ2*runw);
	}
    }
  else
    W = q2/(pow(q2-mW2,2)+mW2*wW2*runw);
}

//running width propagators using offsets as in Eqs (2.3) (2.4) (2.5) of http://inspirehep.net/record/263569 (fully equivalent to prop::set(q2))
void prop::set_offset(double q2)
{
  Z = 0;
  W = 0;
  G = 0;
  ZG = 0;

  double gamma;

  //fixed width
  if (!opts.runningwidth)
    gamma = 0.;

  //running width
  else
    if (opts.nproc == 1 || opts.nproc == 2)
      gamma = coupling::wwidth/coupling::wmass;
    else if (opts.nproc == 3)
      gamma = coupling::zwidth/coupling::zmass;

  //Eqs (2.3) (2.4) (2.5) of http://inspirehep.net/record/263569
  double zmassp = coupling::zmass / sqrt(1+pow(gamma,2));
  double zwidthp = coupling::zwidth / sqrt(1+pow(gamma,2));
  double wmassp = coupling::wmass / sqrt(1+pow(gamma,2));
  double wwidthp = coupling::wwidth / sqrt(1+pow(gamma,2));
  complex <double> Gfp = coupling::Gf / (1.+complex<double> (0,gamma));

  double mZp2 = pow(zmassp, 2);
  double wZp2 = pow(zwidthp,2);
  double mWp2 = pow(wmassp, 2);
  double wWp2 = pow(wwidthp,2);

  if (opts.nproc == 3)
    {
      Z = q2/(pow(q2-mZp2,2)+mZp2*wZp2) * pow(abs(Gfp)/coupling::Gf,2);
      if (opts.useGamma)
	{
	  G = 1./q2;
	  ZG = (q2-mZp2-gamma*sqrt(mZp2*wZp2))/(pow(q2-mZp2,2)+mZp2*wZp2) * pow(abs(Gfp)/coupling::Gf,2);
	  //propZG = (q2-mZp2                      )/(pow(q2-mZp2,2)+mZp2*wZp2) * pow(abs(Gfp)/coupling::Gf,2); //without gamma term
	}
    }
  else
    W = q2/(pow(q2-mWp2,2)+mWp2*wWp2) * pow(abs(Gfp)/coupling::Gf,2);
}
