#include "mcfm/mcfm_interface.h"
#include "resum/gaussrules.h"
#include "resum/pdfevol.h"
#include "resum/mellinint.h"
#include "resum/mesq.h"
#include "resum/rapint.h"
#include "resum/resint.h"
#include "resum/pegasus.h"
#include "resum/anomalous.h"
#include "resum/resconst.h"
#include "vjfo/vjint.h"
#include "vjfo/vjloint.h"
#include "interface.h"
#include "coupling.h"
#include "switch.h"
#include "cubacall.h"
#include "settings.h"
#include "dyturbo.h"

#include <cuba.h>
#include <math.h>

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
    //Initialise some DYRES settings
    g_param_.g_param_ = opts.g_param;
    nnlo_.order_ = opts.order;            //order (0=LO, 1=NLO, 2=NNLO)
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
    vjloint::init();
    //
    switching::init(); //switching function initialisation
    rescinit_();
    // cuba init
    cubacores(opts.cubacores,1000000);   //< set number of cores (move this to cubainit)
    cubainit((void (*)()) initfun,NULL); //< merge at the end of the run
    cubaexit((void (*)()) exitfun,NULL); //< merge at the end of the run
    // histogram output
}

