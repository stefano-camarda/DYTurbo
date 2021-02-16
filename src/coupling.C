#include "coupling.h"
#include "dyres_interface.h"
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
double coupling::zwidth;
double coupling::wwidth;

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

  //electroweak schemes as in MCFM 6.8
  if (opts.ewscheme == -1)
    {
      // This is the MCFM default, corresponding to an effective
      // field theory approach valid for scales below the top-mass
      // (see Georgi, Nucl. Phys. B 363 (1991) 301).
      // There are 4 inputs here instead of the usual 3 ... --> I guess alpha(mz) is set to an effective value
      Gf = opts.Gf;
      wmass = opts.wmass;
      zmass = opts.zmass;
      aemmz  = opts.aemmz;

      // xw is derived using wmass=zmass*dsqrt(rho)*cos(theta_w)
      xw  = 4*M_PI*aemmz/(8.*pow(wmass,2)*Gf/sqrt(2));

      //is the top mass actually used at all?
      //dymasses_.mt_  = sqrt(16.*pow(M_PI,2)/3./sqrt(2)/Gf*(pow(wmass,2)/(zmass,2)/(1.-xw)-1.));
    }
  else if (opts.ewscheme == 0)
    {
      //------------------------------------------------------------
      //     option=0 : MadEvent default (= AlpGen with iewopt=2)
      //------------------------------------------------------------

      //input values: sin2thetaW, alphaEM(mZ), mZ
      xw  = opts.xw;
      aemmz  = opts.aemmz;
      zmass  = opts.zmass;
      
      //derived: mW, Gf
      wmass  = zmass * sqrt(1. - xw);
      Gf = aemmz * 4*M_PI/xw/(8.*pow(wmass,2)/sqrt(2));
    }
  else if (opts.ewscheme == 1)
    {
      //-----------------------------------------------------
      //     option=1 : LUSIFER and AlpGen (iewopt=3) default
      //-----------------------------------------------------
	
      //Gmu scheme, input values: mZ, mW, Gf
      Gf = opts.Gf;
      wmass  = opts.wmass;
      zmass  = opts.zmass;
      
      //derived: sin2thetaW, alphaEM(MZ)
      xw  = 1.-pow(wmass/zmass,2);
      aemmz  = sqrt(2)*Gf*pow(wmass,2)*xw/M_PI;
    }
  else if (opts.ewscheme == 2)
    {
      //-------------------------------------------------------------------
      //     option=2 : W and Z mass are derived from couplings
      //-------------------------------------------------------------------

      //input values: Gf, alphaEM(mZ), sin2thetaW
      Gf = opts.Gf;
      aemmz  = opts.aemmz;
      xw  = opts.xw;
      
      //derived: mW, mZ
      wmass  = sqrt(aemmz*M_PI/xw/Gf/sqrt(2));
      zmass  = wmass/sqrt(1.-xw);
    }
  else if (opts.ewscheme == 3)
    {
      //-----------------------------------------------------------------
      //     option=3 : USER choice : you should know what you're doing!!
      //-----------------------------------------------------------------
      Gf = opts.Gf;
      aemmz  = opts.aemmz;
      xw  = opts.xw;
      wmass  = opts.wmass;
      zmass  = opts.zmass;
    }
  else
    {
      cout << "opts.ewscheme = " << opts.ewscheme << " is not a valid input." << endl;
      exit (-1);
    }

  ewcouple_.Gf_ = Gf;
  ewcouple_.xw_ = xw;
  dymasses_.wmass_ = wmass;
  dymasses_.zmass_ = zmass;

  //width values
  zwidth = opts.zwidth;
  wwidth = opts.wwidth;

  //Now set up the other derived parameters

  //W coupling
  //ewcouple_.gwsq_ = 4 * M_PI * aemmz/xw;     //--> Original MCFM
  ewcouple_.gwsq_ = 4*sqrt(2)*Gf*pow(wmass,2); //--> W coupling depending only on Gf and wmass (equal to the above only in gauge invariant schemes)
  ewcouple_.gw_=sqrt(ewcouple_.gwsq_);

  //photon coupling (used also for Z)
  //ewcouple_.esq_ = ewcouple_.gwsq_* xw; //--> Original MCFM
  ewcouple_.esq_ = 4 * M_PI * aemmz;      //--> Photon coupling depending only on alpha EM (equal to the above only in gauge invariant schemes)

  //The Z coupling are evaluated as esq/(sin2w)^2 => Gf * mw^2 / cw^2 = Gf * mz^2

  //calculate the couplings as given in Kunszt and Gunion
  //Modified to notation of DKS (ie divided by 2*sw*cw)
  //xw=sin^2 theta_w
  
  //zcouple_.sin2w_=2.*sqrt(xw*(1.-xw)); //--> Original MCFM (this xw must be the on-shell xw = 1 - mW^2/mZ^2, works only for gauge invariant schemes)
  zcouple_.sin2w_=sqrt(2*sqrt(2)*M_PI*aemmz/(Gf*pow(zmass,2))); //This expression should work in all cases

  for (int j=0; j < MAXNF; j++)
    {
      zcouple_.l_[j]=(ewcharge_.tau_[j+MAXNF+1]-2.*ewcharge_.Q_[j+MAXNF+1]*xw)/zcouple_.sin2w_;
      zcouple_.r_[j]=(-2*ewcharge_.Q_[j+MAXNF+1]*xw)/zcouple_.sin2w_;
    }

  zcouple_.le_=(-1.-2.*(-1.)*xw)/zcouple_.sin2w_;
  zcouple_.re_=(-2.*(-1.)*xw)/zcouple_.sin2w_;

  //are ln and rn ever used? --> No, they are the Z -> nunu couplings
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
  //pdfini_();
  //read g from the PDF
  //setg();
  //take the cmass and b mass from the PDF
  //      cmass=dsqrt(mcsq)
  //      bmass=dsqrt(mbsq)
  //  scale_.musq_ = pow(scale_.scale_,2);
}      

void coupling::SMparameters()
{
  //None of these parameters is actually needed, check and clean up and remove this function
  //(The squared HF masses are actually used in the MCFM function for the running of alphas)
  
  //Calculational scheme for EW couplings
  //     ewscheme=-1  : MCFM default 
  //                    input values = Gf,alpha(m_Z),m_W,m_Z
  //                    output values = sin^2(theta_W),mtop
  //
  //     ewscheme=1   : New Madevent default, "G_mu scheme"
  //                    = LUSIFER and AlpGen (iewopt=3) defaults
  //                    input values = G_F,m_Z,m_W
  //                    output values = sin^2(theta_W),alpha(m_Z).

  //  ewscheme_.ewscheme_=1;

  //  ewinput_.Gf_inp_= 1.1663787e-5;
  //  ewinput_.wmass_inp_= 80.385;
  //  ewinput_.zmass_inp_= 91.1876;
  //  ewinput_.aemmz_inp_= 7.7585538055706e-03;
  //  ewinput_.xw_inp_= 0.2312;

  //  dymasses_.wwidth_ = 2.091;
  //  dymasses_.zwidth_ = 2.4950;


  // ******************************* The following parameters are not used ***************************
  //Masses, widths and initial-state flavour information
  // Masses: note that "mtausq", "mcsq" and "mbsq" are typically used
  // throughout the program to calculate couplings that depend on the
  // mass, while "mtau","mc" and "mb" are the masses that appear in
  // the rest of the matrix elements and phase space (and may be set
  // to zero in the program, depending on the process number) 
  dymasses_.mtausq_ = 3.157729;
  dymasses_.mcsq_ = 2.25;
  dymasses_.mbsq_ = 21.3444;
  dymasses_.mtau_ = 1.777;
  dymasses_.mc_ = 1.5;  //-> read HF masses from the PDF
  dymasses_.mb_ = 4.62; //-> read HF masses from the PDF
  dymasses_.mt_ = 178;

  // Widths: note that the top width is calculated in the program
  dymasses_.tauwidth_ = 2.269e-12;

  // Masses below here are currently unused      
  dymasses_.md_ = 5e-3;
  dymasses_.mu_ = 5e-3;
  dymasses_.ms_ = 1e-1;
  dymasses_.mel_ = 0.510997e-3;
  dymasses_.mmu_ = 0.105658389;

  //Dim. Reg. parameter epsilon (not used)
  epinv_.epinv_ = 1000;
  epinv2_.epinv2_= 1000;
  // ************************************************************************************************
}

void coupling::initscales()
{
  //Set up QCD scales in the fortran common blocks
  scale_.scale_               = opts.kmuren*opts.rmass;
  facscale_.facscale_         = opts.kmufac*opts.rmass;
  a_param_.a_param_           = 1./opts.kmures;

  scale_.musq_=pow(scale_.scale_,2);

  //Set all factorization scales of the dipole contributions to facscale
  //to avoid problems when dynamicscale=.false.
  for (int nd =0; nd <= 40; nd++)
    dipolescale_.dipscale_[nd]=facscale_.facscale_;

  //Set the scale pre-factors in the fortran common block
  scaleopts_.kmuren_ = opts.kmuren;
  scaleopts_.kmufac_ = opts.kmufac;
  scaleopts_.kmures_ = opts.kmures;
  
  //initialize alpha_s
  pdf::setalphas();
}
