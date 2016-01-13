#include "interface.h"
#include "settings.h"
#include "init.h"
#include "pdf.h"
#include "coupling.h"

#include <math.h>
#include <iostream>

map <int,string> plabel;

void SMparameters()
{
  //Calculational scheme for EW couplings
  //     ewscheme=-1  : MCFM default 
  //                    input values = Gf,alpha(m_Z),m_W,m_Z
  //                    output values = sin^2(theta_W),mtop
  //
  //     ewscheme=1   : New Madevent default, "G_mu scheme"
  //                    = LUSIFER and AlpGen (iewopt=3) defaults
  //                    input values = G_F,m_Z,m_W
  //                    output values = sin^2(theta_W),alpha(m_Z).

  ewscheme_.ewscheme_=1;

  ewinput_.Gf_inp_= 1.1663787e-5;
  ewinput_.aemmz_inp_= 7.7585538055706e-03;
  ewinput_.xw_inp_= 0.2312;
  ewinput_.wmass_inp_= 80.385;
  ewinput_.zmass_inp_= 91.1876;


  dymasses_.wwidth_ = 2.091;
  dymasses_.zwidth_ = 2.4950;

  //     CKM matrix entries
  cabib_.Vud_ = 0.97427;
  cabib_.Vus_ = 0.2253;
  cabib_.Vub_ = 0.00351;
  cabib_.Vcd_ = 0.2252;
  cabib_.Vcs_ = 0.97344;
  cabib_.Vcb_ = 0.0412;


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

//rewritten initialisation functions
void dyturboinit()
{
  //move here the flaq.eq.0 initialisation part of resumm() in main2 instead of using this initialisation flag
  flag_.flag_ = false;

  noglue_.noglue_=false;
  noglue_.ggonly_=false;
  noglue_.gqonly_=false;

  colc_.colourchoice_=0;
  double rtsmin = 40;
  cutoff_.cutoff_=0.001; //add to settings (what is this cutoff?)

  flags_.qflag_=true;
  flags_.gflag_=true;

  //Catani-Seymour subtraction cut-offs for initial-initial, initial-final, final-initial, and final-final dipoles
  //(add to settings)
  alfacut_.aii_=1.;//0.01;
  alfacut_.aif_=1.;//0.01;
  alfacut_.afi_=1.;//0.01;
  alfacut_.aff_=1.;//0.01;

  // Dynamic scale (if true muf=mur=q)
  dynamicscale_.dynamicscale_=false;

  bool removebr=false;

  //Set all factorization scales to facscale
  //to avoid problems when dynamicscale=.false.
  for (int nd =0; nd <= 40; nd++)
    dipolescale_.dipscale_[nd]=facscale_.facscale_;

  //Cut on qt/Q
  //(add to settings)
  qtcut_.xqtcut_=0.008;

  //Limits on invariant mass of vector boson decay products
  //(irrelevant if zerowidth=true)
  limits_.wsqmin_=pow(opts.mlow,2);
  limits_.wsqmax_=pow(opts.mhigh,2);

  //Check if the limits are compatible with sroot
  if (limits_.wsqmax_> pow(energy_.sroot_,2))
    limits_.wsqmax_=pow(energy_.sroot_,2);

  plabel[1]="pp";
  plabel[2]="pp";

  //dycoupling_();
  coupling::init();

  dymasses_.mb_=0; //probably irrelevant to set the bottom mass to 0, since it is not used anywhere
     
  if (!((density_.ih1_==1) && (density_.ih2_ == -1))
      && !((density_.ih1_==1) && (density_.ih2_ == 1)))
    {
      cout << "The selected initial state is not valid," << endl;
      cout << "please select proton-proton (1 1) or proton-antiproton (1 -1)" << endl;
      exit(0);
    }

  //the default behaviour is to remove no branching ratio
  brnrat_.brnrat_=1.;

  //branching ratios
  double brwen,brzee,brtau,brtop;
  branch_(brwen,brzee,brtau,brtop);

  //W+->l+nubar
  if (nproc_.nproc_ == 1)
    {
      plabel[3]="nl";
      plabel[4]="ea";
      plabel[5]="pp";
      plabel[6]="pp";
      nwz_.nwz_=1;
      breit_.mass3_=dymasses_.wmass_;
      breit_.width3_=dymasses_.wwidth_;
      if (removebr)
        brnrat_.brnrat_=brwen;
      
    }
  //W-=>l-nu
  else if(nproc_.nproc_ == 2)
    {
      plabel[3]="el";
      plabel[4]="na";
      plabel[5]="pp";
      plabel[6]="pp";
      nwz_.nwz_=-1;
      breit_.mass3_=dymasses_.wmass_;
      breit_.width3_=dymasses_.wwidth_;
      if (removebr)
        brnrat_.brnrat_=brwen;
    }
  else if(nproc_.nproc_==3)
    {
      plabel[3]="el";
      plabel[4]="ea";
      plabel[5]="pp";
      plabel[6]="pp";
      nwz_.nwz_=0;
      breit_.mass3_=dymasses_.zmass_;
      breit_.width3_=dymasses_.zwidth_;

      zcouple_.l1_=zcouple_.le_;
      zcouple_.r1_=zcouple_.re_;

      //with q1=0 the photon contribution is switch off
      //with q1=-1 the photon contribution is switch on

      dyphoton_.phot_=zcouple_.q1_;

      if (removebr)
	brnrat_.brnrat_=brzee;
    }
  else
    {
      cout << "Wrong process" << endl;
      exit(0);
    }

  //Breit Wigner unweighting settings
  //the jets, particles 5 and 6, are not resonant
  breit_.n2_=0; 
  breit_.mass2_=0;
  breit_.width2_=0;

  //the leptons, particles 3 and 4, are resonant
  breit_.n3_=1;

  //Decide if worth using Breit-Wigner unweighting or not for the dilepon system
  //n3=0
  //vmass=wmass
  //if(nproc.eq.3) vmass=zmass
  //if(mwmin.lt.vmass.and.mwmax.gt.vmass) n3=1

  ckmfill_(nwz_.nwz_);

  setalphas(); //probably not needed, since already set in dycoupling

  scale_.musq_=pow(scale_.scale_,2);

  // Initialize efficiency variables (get rid of these)
  efficiency_.njetzero_ = 0;
  efficiency_.ncutzero_ = 0;
  efficiency_.ntotzero_ = 0;
  efficiency_.ntotshot_ = 0;
  
  // Set-up incoming beams and PS integration cut-offs
  rtsmin = min (rtsmin, sqrt(limits_.wsqmin_ + cutoff_.cutoff_));

  if (zerowidth_.zerowidth_)
    {
      rtsmin = dymasses_.wmass_;
      if (nproc_.nproc_ == 3) rtsmin = dymasses_.zmass_;
    }

  taumin_.taumin_=pow((rtsmin/energy_.sroot_),2);
  taumin_.logtaumin_=log(taumin_.taumin_);
  xmin_.xmin_=taumin_.taumin_;

  pext_.p1ext_[3]=-0.5*energy_.sroot_;
  pext_.p1ext_[0]=0.;
  pext_.p1ext_[1]=0.;
  pext_.p1ext_[2]=-0.5*energy_.sroot_;

  pext_.p2ext_[3]=-0.5*energy_.sroot_;
  pext_.p2ext_[0]=0.;
  pext_.p2ext_[1]=0.;
  pext_.p2ext_[2]=+0.5*energy_.sroot_;
}
