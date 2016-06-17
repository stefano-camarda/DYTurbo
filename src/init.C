#include "banner.h"
#include "interface.h"
#include "mcfm_interface.h"
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
#include "anomalous.h"
#include "resint.h"
#include "vjint.h"
#include "switch.h"
#include "plotter.h"
#include "printsettings.h"
#include "cubacall.h"

#include <cuba.h>
#include <math.h>
#include <iostream>
#include <cstring>

//rewritten initialisation functions
void dyturboinit(int argc, char * argv[])
{
  banner();
  gaussinit_();             //initialisation of fortran gaussian quadrature nodes and weights
  coupling::SMparameters(); //initialisation of unused MCFM parameters
  
  opts.parse_options(argc,argv);

  dofill_.doFill_ = 0;
  
  //Initialise some DYRES settings
  g_param_.g_param_ = opts.g_param;
  nnlo_.order_ = opts.order;            //order (1=NLO, 2=NNLO)
  opts_.fixedorder_  = opts.fixedorder; //fixed order/resummation switch
  qtsub_.xqtcut_= opts.xqtcut;          //Cut on qt/Q
  qtsub_.qtcut_= opts.qtcut;            //Cut on qt
  //move here the flaq.eq.0 initialisation part of resumm() in main2 instead of using this initialisation flag
  
  mcfm::init();
  
  iniflavreduce_(); //need to call this after nproc_.nproc_ is set

  coupling::initscales();

  //C++ resum
  //initialise all the C modules
  gr::init(); //nodes and weights of gaussian quadrature rules
  mellinint::initgauss(); //gaussian quadrature for mellin inversion
  mesq::init(); //EW couplings for born amplitudes
  rapint::init(); //allocate memory for the rapidity quadrature
  resconst::init(); //calculate beta, A and B coefficients
  anomalous::init(); //calculate anomalous dimensions, C1, C2 and gamma coefficients
  pdfevol::init(); //transform the PDF from x- to N-space at the factorisation scale
  pegasus::init(); //initialise Pegasus QCD and transform the PDF from x- to N-space at the starting scale
  resint::init(); //initialise dequad integration for the bessel integral
  //end C++ resum

  //V+j fixed order initialisation
  vjint::init();
  
  switching::init(); //switching function initialisation
  rescinit_();
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
