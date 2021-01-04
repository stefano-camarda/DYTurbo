#include "scales.h"
#include "coupling.h"
#include "settings.h"
#include "interface.h"
#include "dyres_interface.h"
#include "vjint.h"
#include "pdf.h"

double scales::ren;
double scales::fac;
double scales::res;
double scales::ren2;
double scales::fac2;
double scales::res2;

string scales::func(int ff)
{
  string f;
  switch (ff)
  {
  case 0:
    f = to_string(opts.rmass);
    break;
  case 1:
    f = "m_ll";
    break;
  case 2:
    f = "sqrt(m_ll^2+p_T^2)";
    break;
  case 3:
    f = "sqrt(m_ll^2+p_T^2+m_jj^2)";
    break;
  case 4: //see slide 27 of https://gsalam.web.cern.ch/gsalam/talks/repo/2016-03-SB+SLAC-Munich-precision.pdf
    //f = "1/2(p_T + sqrt(m_ll^2+p_T^2))";
    f = "p_T + sqrt(m_ll^2+p_T^2)";
    break;
  }
  return f;
}

void scales::form(double &scale2, int ff, double m, double pt, double mjj)
{
  switch (ff)
  {
  case 0:
    scale2 = pow(opts.rmass,2);
    break;
  case 1:
    scale2 = pow(m,2);
    break;
  case 2:
    scale2 =  pow(m,2) + pow(pt,2);
    break;
  case 3:
    scale2 =  pow(m,2) + pow(pt,2) + pow(mjj,2);
    break;
  case 4: //see slide 27 of https://gsalam.web.cern.ch/gsalam/talks/repo/2016-03-SB+SLAC-Munich-precision.pdf
    scale2 = pow(pt + sqrt(m*m+pt*pt),2);
    break;
  }
}

void scales::set(double m, double pt, double mjj)
{
  //set scales according to the chosen functional form
  form(ren2, opts.fmuren, m, pt, mjj);
  form(fac2, opts.fmufac, m, pt, mjj);
  form(res2, opts.fmures, m);

  ren = sqrt(ren2);
  fac = sqrt(fac2);
  res = sqrt(res2);
  
  //Rescale by factors
  ren *= opts.kmuren;
  fac *= opts.kmufac;
  res *= opts.kmures;

  //catch absurdly large scales --> no good reason to do it
  //double scalemax = 10000.;
  //
  //if (ren > scalemax) {ren = scalemax; ren2 = pow(ren,2);}
  //if (fac > scalemax) {fac = scalemax; fac2 = pow(fac,2);}
  //if (res > scalemax) {res = scalemax; res2 = pow(res,2);}
}

//set mcfm scales and couplings
void scales::mcfm()
{
  //Set up QCD scales in the fortran common blocks
  scale_.scale_ = ren;
  scale_.musq_ = pow(ren,2);
  facscale_.facscale_ = fac;

  /*
  //Set all factorization scales of the dipole contributions to facscale
  //to avoid problems when dynamicscale=.false.
  for (int nd =0; nd <= 40; nd++)
    dipolescale_.dipscale_[nd] = fac;
  */
}

void scales::dyres(double m)
{
  a_param_.a_param_ = m/res;
}

//set scales and alpha strong in vjet analytical calculation
void scales::vjet()
{
  scales2_.xmur_ = ren;
  scales2_.xmuf_ = fac;
  scales2_.xmur2_ = pow(ren,2);
  scales2_.xmuf2_ = pow(fac,2);

  //Set alphas in the fortran common blocks
  pdf::setalphas();
}

//Fortran interface
void scaleset_mcfm_(double &m, double &pt, double &mjj)
{
  scales::set(m, pt, mjj);
  scales::mcfm();
  scales::dyres(m);

  //Set alphas and strong coupling in the fortran common blocks
  pdf::setalphas();
}

int dynamic_fac_()
{
  return (opts.fmufac > 0);
}
