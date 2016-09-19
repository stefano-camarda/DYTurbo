#include "mcfm_interface.h"
#include "settings.h"
#include "coupling.h"
#include "resconst.h"
#include "phasespace.h"

#include <cstring>
#include <iostream>
#include <math.h>

map <int,string> plabel;

void mcfm::init()
{
  energy_.sroot_              = opts.sroot;
  density_.ih1_               = opts.ih1;
  density_.ih2_               = opts.ih2;
  nproc_.nproc_               = opts.nproc;
  zerowidth_.zerowidth_       = opts.zerowidth;
  dynamicscale_.dynamicscale_ = opts.dynamicscale;   // Dynamic scale (if true mufac=muren=q)

  //    CKM matrix entries
  cabib_.Vud_ = opts.Vud;
  cabib_.Vus_ = opts.Vus;
  cabib_.Vub_ = opts.Vub;
  cabib_.Vcd_ = opts.Vcd;
  cabib_.Vcs_ = opts.Vcs;
  cabib_.Vcb_ = opts.Vcb;

  dymasses_.wwidth_ = opts.wwidth;
  dymasses_.zwidth_ = opts.zwidth;
  
  //initialise MCFM settings
  noglue_.noglue_=false;
  noglue_.ggonly_=false;
  noglue_.gqonly_=false;

  string part = "tota";
  strncpy(part_.part_ , part.c_str(), part.size());

  colc_.colourchoice_=0;
  double rtsmin = 40;
  cutoff_.cutoff_=0.001; //add to settings (this is a cut off in mcfm on invariant mass pairs between emitted and radiator)

  flags_.qflag_=true;
  flags_.gflag_=true;

  //Catani-Seymour subtraction cut-offs for initial-initial, initial-final, final-initial, and final-final dipoles
  //(add to settings)
  alfacut_.aii_=1.;//0.01;
  alfacut_.aif_=1.;//0.01;
  alfacut_.afi_=1.;//0.01;
  alfacut_.aff_=1.;//0.01;

  //Limits on invariant mass of vector boson decay products
  //(irrelevant if zerowidth=true)
  limits_.wsqmin_=pow(phasespace::mmin,2);
  limits_.wsqmax_=pow(phasespace::mmax,2);

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
  bool removebr=false;
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
      breit_.mass3_= dymasses_.wmass_;
      breit_.width3_= dymasses_.wwidth_;
      opts.rmass = dymasses_.wmass_;
      opts.rwidth = dymasses_.wwidth_;
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
      opts.rmass = dymasses_.wmass_;
      opts.rwidth = dymasses_.wwidth_;
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
      opts.rmass = dymasses_.zmass_;
      opts.rwidth = dymasses_.zwidth_;

      zcouple_.l1_=zcouple_.le_;
      zcouple_.r1_=zcouple_.re_;

      //with q1=0 the photon contribution is switched off
      //with q1=-1 the photon contribution is switched on

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

  // Set-up PS integration cut-offs
  rtsmin = min (rtsmin, sqrt(limits_.wsqmin_ + cutoff_.cutoff_));

  if (zerowidth_.zerowidth_)
    {
      rtsmin = dymasses_.wmass_;
      if (nproc_.nproc_ == 3) rtsmin = dymasses_.zmass_;
    }

  taumin_.taumin_=pow((rtsmin/energy_.sroot_),2);
  taumin_.logtaumin_=log(taumin_.taumin_);
  xmin_.xmin_=taumin_.taumin_;

  //Set-up incoming beams (used only in realint.f)
  pext_.p1ext_[3]=-0.5*energy_.sroot_;
  pext_.p1ext_[0]=0.;
  pext_.p1ext_[1]=0.;
  pext_.p1ext_[2]=-0.5*energy_.sroot_;

  pext_.p2ext_[3]=-0.5*energy_.sroot_;
  pext_.p2ext_[0]=0.;
  pext_.p2ext_[1]=0.;
  pext_.p2ext_[2]=+0.5*energy_.sroot_;

  //set NF to 5 (it is used in H2calc)
  nf_.nf_ = resconst::NF;
}
