#include "coupling.h"
#include "interface.h"
#include "settings.h"
#include "pdf.h"

#include <iostream>
#include <math.h>

using namespace coupling;

double coupling::aemmz;
double coupling::Gf;
double coupling::xw;
double coupling::zmass;
double coupling::wmass;

//Number of colours in QCD
const double coupling::NC = 3.;

void coupling::init()
{
  //Set up all the couplings used in MCFM

  ewcharge_.Q_[MAXNF-5] = +0.333333333333333 * opts.Zbb;
  ewcharge_.Q_[MAXNF-4] = -0.666666666666667 * opts.Zcc;
  ewcharge_.Q_[MAXNF-3] = +0.333333333333333 * opts.Zss;
  ewcharge_.Q_[MAXNF-2] = -0.666666666666667 * opts.Zuu;
  ewcharge_.Q_[MAXNF-1] = +0.333333333333333 * opts.Zdd;
  ewcharge_.Q_[MAXNF]   =  0.	       ;
  ewcharge_.Q_[MAXNF+1] = -0.333333333333333 * opts.Zdd;
  ewcharge_.Q_[MAXNF+2] = +0.666666666666667 * opts.Zuu;
  ewcharge_.Q_[MAXNF+3] = -0.333333333333333 * opts.Zss;
  ewcharge_.Q_[MAXNF+4] = +0.666666666666667 * opts.Zcc;
  ewcharge_.Q_[MAXNF+5] = -0.333333333333333 * opts.Zbb;

  ewcharge_.tau_[MAXNF-5] =  1. * opts.Zbb;
  ewcharge_.tau_[MAXNF-4] = -1. * opts.Zcc;
  ewcharge_.tau_[MAXNF-3] =  1. * opts.Zss;
  ewcharge_.tau_[MAXNF-2] = -1. * opts.Zuu;
  ewcharge_.tau_[MAXNF-1] =  1. * opts.Zdd;
  ewcharge_.tau_[MAXNF]   =  0.           ;
  ewcharge_.tau_[MAXNF+1] = -1. * opts.Zdd;
  ewcharge_.tau_[MAXNF+2] =  1. * opts.Zuu;
  ewcharge_.tau_[MAXNF+3] = -1. * opts.Zss;
  ewcharge_.tau_[MAXNF+4] =  1. * opts.Zcc;
  ewcharge_.tau_[MAXNF+5] = -1. * opts.Zbb;


  //Gmu scheme, inputs: Gf, MZ, MW
  Gf = opts.Gf;
  zmass = opts.zmass;
  wmass = opts.wmass;

  if (opts.aemmz == 0. && opts.xw == 0.)
    {
      //Gmu scheme: tree level derived sin2thetaW and alphaEM(MZ)
      xw  = 1.-pow(wmass/zmass,2);
      aemmz = sqrt(2)*Gf*pow(wmass,2)*xw/M_PI;
    }
  else if (opts.aemmz != 0. && opts.xw != 0.) 
    {
      //Scheme completely determined by the user
      aemmz = opts.aemmz;
      xw = opts.xw;
    }
  else if (opts.aemmz != 0.)
  //If the input aemmz is non zero the scheme
  //used is the old MCFM default, corresponding to an effective
  //field theory approach valid for scales below the top-mass
  //(see Georgi, Nucl. Phys. B 363 (1991) 301).
    {
      aemmz = opts.aemmz;
      xw  = 4*M_PI*aemmz/(8.*pow(wmass,2)*Gf/sqrt(2));
      dymasses_.mt_  = sqrt(16.*pow(M_PI,2)/3./sqrt(2)/Gf*(pow(wmass,2)/(zmass,2)/(1.-xw)-1.)); //is the top mass actually used at all?
    }
  else if (opts.xw != 0.)
    {
      //Effective value of the weak mixing angle
      xw = opts.xw;
      //!!!!! Check this value, should use  1.-pow(wmass/zmass,2) instead of xw?
      //      aemmz = sqrt(2)*Gf*pow(wmass,2)*xw/M_PI;
      aemmz = sqrt(2)*Gf*pow(wmass,2)*(1.-pow(wmass/zmass,2))/M_PI;
    }


  /*  if (ewscheme_.ewscheme_ == -1)
    {
      //This is the old MCFM default, corresponding to an effective
      //field theory approach valid for scales below the top-mass
      //(see Georgi, Nucl. Phys. B 363 (1991) 301).

      //inputs: Gf, alphaEM(MZ), MW, MZ
      Gf = ewinput_.Gf_inp_;
      aemmz = ewinput_.aemmz_inp_;
      wmass = ewinput_.wmass_inp_;
      zmass = ewinput_.zmass_inp_;
        
      //Derived: sin2thetaW
      xw  = 4*M_PI*aemmz/(8.*pow(wmass,2)*Gf/sqrt(2));

      //is the top mass actually used at all?
      dymasses_.mt_  = sqrt(16.*pow(M_PI,2)/3./sqrt(2)/Gf*(pow(wmass,2)/(zmass,2)/(1.-xw)-1.));
    }
  else if (ewscheme_.ewscheme_ = 1)
    {
      //Gmu scheme, inputs: MZ, MW, Gf
      zmass  = ewinput_.zmass_inp_;
      wmass  = ewinput_.wmass_inp_;
      Gf = ewinput_.Gf_inp_;
	   
      //Derived: sin2thetaW, alphaEM(MZ)
      xw  = 1.-pow(wmass/zmass,2);
      aemmz = sqrt(2)*Gf*pow(wmass,2)*xw/M_PI;
    }
  else
    { 
      cout << "ewscheme = " << ewscheme_.ewscheme_ << " is not implemented." << endl;
      exit(0);
    }
  */

  ewcouple_.Gf_ = Gf;
  ewcouple_.xw_ = xw;
  dymasses_.wmass_ = wmass;
  dymasses_.zmass_ = zmass;

  //Now set up the other derived parameters

  //W coupling
  ewcouple_.gwsq_= 4 * M_PI * aemmz/xw; //= 4*sqrt(2)*Gf*pow(wmass,2);
  ewcouple_.gw_=sqrt(ewcouple_.gwsq_);

  //photon coupling (used also for Z)
  ewcouple_.esq_= ewcouple_.gwsq_* xw; //= 4 * M_PI * aemmz;

  //calculate the couplings as given in Kunszt and Gunion
  //Modified to notation of DKS (ie divided by 2*sw*cw)
  //xw=sin^2 theta_w
  zcouple_.sin2w_=2.*sqrt(xw*(1.-xw));
  for (int j=0; j < MAXNF; j++)
    {
      zcouple_.l_[j]=(ewcharge_.tau_[j+MAXNF+1]-2.*ewcharge_.Q_[j+MAXNF+1]*xw)/zcouple_.sin2w_;
      zcouple_.r_[j]=(-2*ewcharge_.Q_[j+MAXNF+1]*xw)/zcouple_.sin2w_;
    }

  zcouple_.le_=(-1.-2.*(-1.)*xw)/zcouple_.sin2w_;
  zcouple_.re_=(-2.*(-1.)*xw)/zcouple_.sin2w_;

  //are ln and rn ever used?
  zcouple_.ln_=(+1.-2.*(+0.)*xw)/zcouple_.sin2w_;
  zcouple_.rn_=0.;

  //switch off the gamma* contribution if required
  zcouple_.q1_ = (opts.useGamma ? -1 :  0 );

  //******************* this coupling is not used ****************
  //Calculate the appropriate Higgs vacuum expectation value.
  //This vevsq is defined so that gwsq/(4*wmass**2)=Gf*rt2=1/vevsq
  //(ie differs from definition in ESW)
  ewcouple_.vevsq_=1./sqrt(2.)/ewcouple_.Gf_;
  //****************************************************************

  //set up the beta-function
  //b0 is defined in mcfm/b0.f and is used in virtint.f and mcfm/dipoles.f
  b0_.b0_=(NC*11.-2.*MAXNF)/6.;

  //initialize the pdf set
  pdfini_();

  //take the cmass and b mass from the PDF
  //      cmass=dsqrt(mcsq)
  //      bmass=dsqrt(mbsq)
  scale_.musq_ = pow(scale_.scale_,2);
 
  //initialize alpha_s
  setalphas();
}      
