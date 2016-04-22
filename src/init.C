#include "banner.h"
#include "interface.h"
#include "settings.h"
#include "init.h"
#include "pdf.h"
#include "pdfevol.h"
#include "pegasus.h"
#include "coupling.h"
#include "gaussrules.h"
#include "mellinint.h"
#include "rapint.h"
#include "mesq.h"
#include "resconst.h"
#include "hcoefficients.h"
#include "anomalous.h"
#include "resint.h"
#include "switch.h"
#include "plotter.h"
#include "printsettings.h"
#include "cubacall.h"

#include <cuba.h>
#include <math.h>
#include <iostream>
#include <cstring>

map <int,string> plabel;

//rewritten initialisation functions
void dyturboinit(string conf_file)
{
  banner();
  coupling::SMparameters();
  opts.readfromfile(conf_file.c_str());
  //  opts.initDyresSettings();

  //Initialise some DYRES and MCFM settings
  energy_.sroot_              = opts.sroot;
  density_.ih1_               = opts.ih1;
  density_.ih2_               = opts.ih2;
  nproc_.nproc_               = opts.nproc;
  g_param_.g_param_           = opts.g_param;
  nnlo_.order_                = opts.order;
  zerowidth_.zerowidth_       = opts.zerowidth;
  dynamicscale_.dynamicscale_ = opts.dynamicscale;

  dofill_.doFill_ = 0;
  
  gaussinit_();
  iniflavreduce_();
  
  //move here the flaq.eq.0 initialisation part of resumm() in main2 instead of using this initialisation flag

  //Cut on qt/Q (add to settings)
  qtcut_.xqtcut_= opts.xqtcut; //0.008;

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
  flag_.flag_ = false;

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

  // Dynamic scale (if true mufac=muren=q)
  dynamicscale_.dynamicscale_=false;

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

  //set NF to 5 (it is used in H2calc)
  nf_.nf_ = resconst::NF;

  coupling::initscales();
  //Set all factorization scales to facscale
  //to avoid problems when dynamicscale=.false.
  for (int nd =0; nd <= 40; nd++)
    dipolescale_.dipscale_[nd]=facscale_.facscale_;

  //C++ resum
  //initialise all the C modules
  gr::init(); //nodes and weights of gaussian quadrature rules
  mellinint::initgauss(); //gaussian quadrature for mellin inversion
  mesq::init(); //EW couplings for born amplitudes
  rapint::init(); //allocate memory for the rapidity quadrature
  resconst::init(); //calculate beta, A and B coefficients
  anomalous::init(); //calculate anomalous dimensions, C1, C2 and gamma coefficients
  hcoefficients::init(); //allocate memory for the H coefficients
  pdfevol::init(); //allocate memory for the pdf in N-space
  pegasus::init(); //initialise Pegasus QCD
  resint::init(); //initialise dequad integration
  //end C++ resum

  switching::init();
  rescinit_();
  //bins.init();
  bins.readfromfile(conf_file.c_str());
  cubacores(opts.cubacores,1000000);   //< set number of cores (move this to cubainit)
  cubaexit((void (*)()) exitfun,NULL); //< merge at the end of the run
  // histogram output
  hists.Init();

  /***********************************/
  //print out EW and QCD parameters and other settings
  if (opts.verbose)
    opts.dumpAll();
  printsettings();
  /***********************************/
}
