#include "interface.h"
#include "settings.h"
#include "init.h"
#include "pdf.h"

#include <math.h>
#include <iostream>

map <int,string> plabel;

//rewrite initialisation functions
void dyturboinit()
{
  //move here the flaq.eq.0 initialisation part of resumm() in main2 instead of using this initialisation flag
  flag_.flag_ = false;

  noglue_.noglue_=false;
  noglue_.ggonly_=false;
  noglue_.gqonly_=false;

  colc_.colourchoice_=0;
  double rtsmin = 40;
  cutoff_.cutoff_=0.001;

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
  limits_.wsqmin_=pow(opts.M_min,2);
  limits_.wsqmax_=pow(opts.M_max,2);

  //Check if the limits are compatible with sroot
  if (limits_.wsqmax_> pow(energy_.sroot_,2))
    limits_.wsqmax_=pow(energy_.sroot_,2);

  plabel[1]="pp";
  plabel[2]="pp";

  dycoupling_();

  dymasses_.mb_=0; //probably irrelevant to set the bottom mass to 0, since it is not used anywhere
     
  if ((density_.ih1_==1) && (density_.ih2_ == -1))
    cout << "Colliding proton-antiproton" << endl;
  else if((density_.ih1_==1) && (density_.ih2_ == 1))
    cout << "Colliding proton-proton" << endl;
  else
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
      cout << "W+ ->  nu(p3)+e+(p4) production " << endl;
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
      cout << "W- ->  e^-(p3)+nu~(p4) production " << endl;
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
      cout << "Z ->  e-(p3)+e+(p4) production " << endl;
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
  breit_.n2_=0; //the jets, particles 5 and 6, are not resonant
  breit_.n3_=1; //the leptons, particles 3 and 4, are resonant

  //Decide if worth using Breit-Wigner unweighting or not for the dilepon system
  //n3=0
  //vmass=wmass
  //if(nproc.eq.3) vmass=zmass
  //if(mwmin.lt.vmass.and.mwmax.gt.vmass) n3=1

  //Print out some information on settings
  if (nnlo_.order_ == 1)
    cout << " Computing NLO(+NLL) cross section" << endl;
  else if (nnlo_.order_ == 2)
    cout << " Computing NNLO(+NNLL) cross section" << endl;
  else
    {
      cout << " Invalid order, please select 1 (NLO) or 2 (NNLO)" << endl;
      exit(0);
    }
  cout << "Center-of-mass energy " << energy_.sroot_ << endl;
    
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
