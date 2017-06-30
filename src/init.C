#include "dyturbo.h"
#include "banner.h"
#include "settings.h"
#include "cubacall.h"
#include "interface.h"
#include "coupling.h"
#include "switch.h"

#include "HistoHandler.h"
#include "dyres_interface.h"
#include "mcfm_interface.h"
#include "gaussrules.h"
#include "clenshawcurtisrules.h"
//#include "chebyshev.h"
#include "pdfevol.h"
#include "mellinint.h"
#include "mesq.h"
#include "rapint.h"
#include "resint.h"
#include "pegasus.h"
#include "anomalous.h"
#include "resconst.h"
#include "vjint.h"
#include "vjloint.h"
#include "abint.h"

#include <cuba.h>
#include <math.h>

void DYTurbo::Init( int argc, char * argv[]){
    banner();
    gaussinit_();             //initialisation of fortran gaussian quadrature nodes and weights
    coupling::SMparameters(); //initialisation of unused MCFM parameters
    // parsing options from input file
    opts.parse_options(argc,argv);
    PrintTable::Init();
    init_params();
    HistoHandler::Init();
    /***********************************/
    //print out EW and QCD parameters and other settings
    if (opts.verbose) opts.dumpAll();
    PrintTable::Settings();
    /***********************************/
}

/**
 * @brief  Initialize parameters of submodules
 *
 * rewritten initialisation functions
 *
 * @param _inArg Description of param
 * @return Description of return.
 */
void DYTurbo::init_params(){
    // init filling
    dofill_.doFill_ = 0;
    dyres::init();
    mcfm::init();
    iniflavreduce_(); //need to call this after nproc_.nproc_ is set
    coupling::initscales();
    cc::init(); //nodes and weights of Clenshaw-Curtis quadrature rules
    //cheb::init(); //Chebyshev nodes for Lagrange interpolation
    //C++ resum
    //initialise all the C modules
    gr::init(); //nodes and weights of gaussian quadrature rules
    mellinint::initgauss(); //gaussian quadrature for mellin inversion
    mesq::init(); //EW couplings for born amplitudes
    rapint::init(); //allocate memory for the rapidity quadrature
    resconst::init(); //calculate beta, A and B coefficients
    if (!opts.fixedorder)
      {
	anomalous::init(); //calculate anomalous dimensions, C1, C2 and gamma coefficients
	pdfevol::init(); //transform the PDF from x- to N-space at the factorisation scale
	pegasus::init(); //initialise Pegasus QCD and transform the PDF from x- to N-space at the starting scale
	resint::init(); //initialise dequad integration for the bessel integral
      }
    //end C++ resum
    //V+j fixed order initialisation
    vjint::init();
    vjloint::init();
    //
    abint::init(); //alfa beta integration initialisation
    switching::init(); //switching function initialisation
    rescinit_();
    // cuba init
    cubacores(opts.cubacores,1000000);   //< set number of cores (move this to cubainit)
    cubainit((void (*)()) initfun,NULL); //< merge at the end of the run
    cubaexit((void (*)()) exitfun,NULL); //< merge at the end of the run
    // histogram output
}

