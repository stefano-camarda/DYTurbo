#include "dyturbo.h"
#include "banner.h"
#include "settings.h"
#include "cubacall.h"
#include "interface.h"
#include "coupling.h"
#include "propagator.h"
#include "pdf.h"
#include "switch.h"
#include "itilde.h"

#include "HistoHandler.h"
#include "Kinematics.h"
#include "dyres_interface.h"
#include "mcfm_interface.h"
#include "gaussrules.h"
#include "clenshawcurtisrules.h"
#include "bequad.h"
//#include "chebyshev.h"
#include "pdfevol.h"
#include "mellinint.h"
#include "mesq.h"
#include "rapint.h"
#include "resint.h"
#include "pegasus.h"
#include "alphas.h"
#include "anomalous.h"
#include "pmom.h"
#include "ccoeff.h"
//#include "ifunc.h"
#include "resconst.h"
#include "alphas.h"
#include "vjint.h"
#include "vjloint.h"
#include "abint.h"

#include <cuba.h>
#include <math.h>

void DYTurbo::Init( int argc, char * argv[])
{
  //Should improve the flow:
  // 1) init_const
  // 2) parse options and input file to know the value of the silent option. Exit with output if options are wrong, print banner if --help is requested
  // 3) print banner
  // 4) any other printout from option parsing
  // ...
  
  //Print DYTurbo banner
  //if (!opts.silent) banner();
  banner();  

  //Initialisation
  init_const();                                       //Initialisation which does not depend on input file
  opts.parse_options(argc,argv);                      //parse options from command line (and parse input file)
  init_params();                                      //Initialisation which depends on input file settings
  if (opts.makehistos) HistoHandler::Init();          //Book hisograms and set root output file name

  //Print out parameters and other settings
  if (opts.verbose) opts.dumpAll();
  if (!opts.silent) PrintTable::Settings();
}

//Initialize constants
void DYTurbo::init_const()
{
  clock_t begin_time, end_time;
  begin_time = clock();
  
  //Integration initialisation
  gaussinit_();             //initialisation of fortran gaussian quadrature nodes and weights
  cc::init();               //nodes and weights of Clenshaw-Curtis quadrature rules
  gr::init();               //nodes and weights of gaussian quadrature rules
  bq::init();               //nodes and weights of Bessel quadrature rules
  //cheb::init();           //Chebyshev nodes for Lagrange interpolation
  itilde::init();           //itilde dequad initialisation
  if (!opts.ctcpp)
    ctquadinit_();          //alfa beta integration initialisation (fortran CT)
  //  cuba::init();

  //Physics constants initialisation
  coupling::SMparameters(); //initialisation of unused MCFM parameters (HF masses are used in the MCFM alphas function)
  resconst::init();         //calculate beta, A and B coefficients
  rescinit_();              //Resummation coefficients (fortran common blocks, still used in the Sudakov)
  pegasus::init_const();    //initialisation of constants in Pegasus QCD
  alphas::init();           //beta function
  
  //Output initialisation (depends on unicode setting)
  //PrintTable::Init();

  end_time = clock();

  //if (opts.timeprofile)
  //cout << "Constant initialisation time: " << float( end_time - begin_time ) /  CLOCKS_PER_SEC << "s" << endl;
}

/**
 * @brief  Initialize parameters of submodules
 *
 * rewritten initialisation functions
 *
 * @param _inArg Description of param
 * @return Description of return.
 */
void DYTurbo::init_params()
{
  // init filling
  dofill_.doFill_ = 0;

  //Output initialisation (depends on unicode setting)
  PrintTable::Init();

  //Set Cuba cores
  cuba::init();
  
  //Initialise parameters
  dyres::init();
  mcfm::init();                          //This functions calls coupling::init()
  pdf::init();                           //Set up PDFs from LHAPDF, set alphas, and read g from the PDF
  iniflavreduce_();                      //need to call this after nproc_.nproc_ is set
  coupling::initscales();                //Set up QCD scales in the fortran common blocks and recompute alphas(mu)
  prop::init();                          //Initialise mass and width used in the propagator

  //C++ resum: initialise all the C modules
  mellinint::initgauss();                //gaussian quadrature for mellin inversion
  mesq::init();                          //EW couplings for born amplitudes
  rapint::init();                        //allocate memory for the rapidity quadrature
  anomalous::init();                     //calculate anomalous dimensions, C1, C2 and gamma coefficients (need to be called after mellinint::initgauss())
  ccoeff::init();                        //calculate, C1, C2, and C3 coefficients
  //ifunc::test();
  pdfevol::init();                       //transform the PDF from x- to N-space at the factorisation scale
  resint::init();                        //initialise dequad integration for the bessel integral, and resummation parameters
  pmom::init();                          //calculate gamma from Pegasus (need to be called after mellinint::initgauss())

  //V+j fixed order initialisation
  vjint::init();
  vjloint::init();

  abint::init();                         //alfa beta integration initialisation
  switching::init();                     //switching function initialisation
  
  Kinematics::init();                    //precompute kinematic cuts transformations
}

//Free memory allocated by init_params
void DYTurbo::release()
{
  mellinint::release();
  rapint::release();
  anomalous::release();
  ccoeff::release();
  pdfevol::release();
  pmom::release();
  vjint::release();
  abint::release();
}
