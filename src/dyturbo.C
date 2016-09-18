#ifndef dyturbo_C
#define dyturbo_C
/**
 * @file dyturbo.C
 * Steering class for Drell-Yan calculation.
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-09-17
 */

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
#include "vjloint.h"
#include "switch.h"
#include "plotter.h"
#include "printsettings.h"
#include "cubacall.h"
#include "cuba.h"
#include "phasespace.h"
#include "ctintegr.h"
#include "HistoHandler.h"

#include "dyturbo.h"

bool DYTurbo::HasOnlyVegas = false;

namespace DYTurbo {

    struct Boundaries{
        std::string name;
        VecDbl data;

        inline BoundariesListItr begin() { return data.begin();};
        inline BoundariesListItr end()   { return data.end(); };
    };

    BoundariesList ActiveBoundaries;

    enum BoundaryIndex {
        b_M=0,
        b_QT,
        b_Y,
        //b_CosTh,
        //b_Phi,
        N_boundaries
    };


    BoundIterator::BoundIterator(){
        // set first boundary
        current.clear();
        for (size_t i = 0; i < N_boundaries; ++i) {
            current.push_back(ActiveBoundaries[i].begin());
        }
        isFirst=true;
    }

    bool BoundIterator::IsEnd(){
        if (isFirst) {
            isFirst=false;
            return false;
        }
        for (size_t i = 0; i < N_boundaries; ++i)
            if (current[i] != ActiveBoundaries[i].begin() )
                return false;
        return true;
    }

    BoundIterator & BoundIterator::operator++(){
        for (int i = N_boundaries-1; i >= 0; --i) { // start from end
            current[i]++; // increase
            // check if it's equal to previous to last (need two numbers as boundaries)
            if(current[i]!=ActiveBoundaries[i].end()-1) break; // go to return
            else current[i]=ActiveBoundaries[i].begin(); // is previous to last 
        }
        return (*this);
    }

    void BoundIterator::Print(){
        printf ("Boundaries ");
        printf ( "| M : %f -- %f " , *current [ b_M  ]  , * ( current [ b_M  ] +1 )  ) ;
        printf ( "| Qt: %f -- %f " , *current [ b_QT ]  , * ( current [ b_QT ] +1 )  ) ;
        printf ( "| Y : %f -- %f " , *current [ b_Y  ]  , * ( current [ b_Y  ] +1 )  ) ;
        printf ( "\n");
    }


    void Term::RunIntegration(){
        VecDbl val;
        double err=0;
        // specialized (need to recorporate)
        if (opts.doBORN && !opts.fixedorder && opts.resint2d ){
            if (opts.resumcpp) rapint::cache(phasespace::ymin, phasespace::ymax);
            else cacheyrapint_(phasespace::ymin, phasespace::ymax);
        }
        integrate(val,err);
        total_int+=val[0];
        total_err2+=err*err;
    }
    Term subtotal;
    TermList ActiveTerms;

    TermIterator::TermIterator() : icurrent(0) { }

    bool TermIterator::IsEnd(){ 
        return icurrent==ActiveTerms.size(); 
    }

    TermIterator & TermIterator::operator++(){
        icurrent++;
        return (*this);
    }

    Term & TermIterator::operator*(){
        return ActiveTerms[icurrent];
    }

    void Init( int argc, char * argv[]){
        banner();
        gaussinit_();             //initialisation of fortran gaussian quadrature nodes and weights
        coupling::SMparameters(); //initialisation of unused MCFM parameters
        // parsing options from input file
        opts.parse_options(argc,argv);
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
        hists.Init();
        /***********************************/
        //print out EW and QCD parameters and other settings
        if (opts.verbose) opts.dumpAll();
        printsettings();
        /***********************************/
    }

    void AddTermIfActive(bool isActive, string name, void (* fun)(VecDbl &val,double &err), bool is_vegas ){
        if (!isActive) return;
        ActiveTerms.push_back({name, fun, 0., 0.});
        HasOnlyVegas&=is_vegas;
    }


    void AddTerms(){
        ActiveTerms.clear();
        HasOnlyVegas=true;
        // born
        bool fixed_born = opts.doBORN && opts.fixedorder;
        AddTermIfActive ( fixed_born && opts.bornint2d      , "fixed born" , bornintegr2d   , false ) ;
        AddTermIfActive ( fixed_born && opts.bornintvegas4d , "fixed born" , bornintegrMC4d , true  ) ;
        AddTermIfActive ( fixed_born && opts.bornintvegas6d , "fixed born" , bornintegrMC6d , true  ) ;
        // resummation
        bool resum_born = opts.doBORN && !opts.fixedorder;
        AddTermIfActive ( resum_born && opts.resint2d     , "resummation" , resintegr2d , false ) ;
        AddTermIfActive ( resum_born && opts.resint3d     , "resummation" , resintegr3d , false ) ;
        AddTermIfActive ( resum_born && opts.resintvegas  , "resummation" , resintegrMC , true ) ;
        // CT
        AddTermIfActive ( opts.doCT && opts.ctint2d      , "counter term" , ctintegr2d , false ) ;
        AddTermIfActive ( opts.doCT && opts.ctint3d      , "counter term" , ctintegr2d , false ) ;
        AddTermIfActive ( opts.doCT && opts.ctintvegas6d , "counter term" , ctintegrMC , true ) ;
        AddTermIfActive ( opts.doCT && opts.ctintvegas8d , "counter term" , ctintegr   , true ) ;
        // VJ finite
        bool vj_finite = opts.doVJ && !opts.doVJREAL && !opts.doVJVIRT;
        AddTermIfActive   ( vj_finite && opts.vjint3d                         , "V+J"    , vjintegr3d   , false ) ;
        AddTermIfActive   ( vj_finite && opts.vjint5d && opts.order == 1      , "V+J LO" , vjlointegr5d , false ) ;
        AddTermIfActive   ( vj_finite && opts.vjintvegas7d && opts.order == 1 , "V+J LO" , vjlointegr7d , true ) ;
        //AddTermIfActive ( vj_finite && opts.vjintvegas7d && opts.order == 1 , "V+J LO" , vjlointegr   , true ) ;
        // VJ NLO
        AddTermIfActive ( opts.doVJREAL , "V+J real"    , vjrealintegr , true ) ;
        AddTermIfActive ( opts.doVJVIRT , "V+J virtual" , vjrealintegr , true ) ;
    }

    void AddBoundary(size_t ib, string name, VecDbl &bins ){
        ActiveBoundaries[ib].name = name;
        ActiveBoundaries[ib].data.clear();
        if (HasOnlyVegas){
            ActiveBoundaries[ib].data.push_back(bins.front());
            ActiveBoundaries[ib].data.push_back(bins.back());
        } else {
            ActiveBoundaries[ib].data = bins;
        }
    }

    void AddBoundaries(){
        ActiveBoundaries.resize(N_boundaries);
        AddBoundary ( b_M  , "q2" , bins.mbins  ) ;
        AddBoundary ( b_Y  , "y"  , bins.ybins  ) ;
        AddBoundary ( b_QT , "qT" , bins.qtbins ) ;
    }

    void WarmUpResummation(){
        /*****************************************/
        //If using the DYRES approximation for PDFs, make sure that the PDF fit is initialised in the same way
        //Need to throw a random point according to a breit wigner, which is used to determine xtauf in the PDF fit
        double costh, m, qt, y;
        int mode = 0;
        if (opts.doBORN || opts.doCT){
            if (opts_.approxpdf_ == 1) {
                srand(opts.rseed);
                double wsqmin = pow(phasespace::mmin ,2);
                double wsqmax = pow(phasespace::mmax ,2);
                double x1=((double)rand()/(double)RAND_MAX);
                double m2,wt;
                breitw_(x1,wsqmin,wsqmax,opts.rmass,opts.rwidth,m2,wt);
                printf( "Initialise PDF fit with mass = %f xtauf = %f \n", sqrt(m2), m2/opts.sroot);
                costh = 0.; m = sqrt(m2); qt = 1; y = 0;
                for (int ipdf=0; ipdf<opts.totpdf; ipdf++){
                    //setpdf_(&ipdf);
                    //setmellinpdf_(&ipdf);
                    resumm_(costh,m,qt,y,mode);
                }
            }
            else {
                costh = 0.; m = opts.rmass; qt = 1; y = 0;
                for (int ipdf=0; ipdf<opts.totpdf; ipdf++){
                    //setpdf_(&ipdf);
                    //setmellinpdf_(&ipdf);
                    resumm_(costh,m,qt,y,mode);
                }
            }
        }
        double f[opts.totpdf];
        ctint_(costh,m,qt,y,mode,f);
        /****************************************/
    }


    void WarmUp(){
        // create term list for calculation
        // check if bounds should be simple (MC) or per each bin (cubature mode)
        AddTerms();
        AddBoundaries();
        // check if all histograms are integrable if necessary (make warning)
        HistoHandler::Book();
        // check if we need to warm up CT integration or Resummation
        WarmUpResummation();
        // clear subtotal
        subtotal = {"TOATL", 0, 0., 0.};
    }

    void SetBounds(BoundIterator bounds){
    }

    void SetTotalBounds(){
    }

    void Terminate(){
    }

    namespace PrintTable {

        void BeginOfRow(){};
        void EndOfRow(){};
        void Hline(){};

        void BoundsAllLooping(bool printNames=false){
            // loop over bounds
        }

        void Header() {
            //BoundsNonLooping();
            // count line width
            // bounds non looping message
            // bounds non looping message
            Hline();
            BeginOfRow();
            BoundsAllLooping(true); // print only names
            for( TermIterator term; !term.IsEnd(); ++term){
                // term->name
            }
            EndOfRow();
            Hline();
        };

        void Footer() {};


        void Bounds() {
            BeginOfRow();
            BoundsAllLooping();
        };

        void Result(const Term &term, bool printGrandTotal) {
        };

        void ResultSubTotal() {
            Result(subtotal);
            EndOfRow();
        };


        void ResultGrandTotal() {
            Hline();
            BeginOfRow();
            SetTotalBounds();
            for( TermIterator term; !term.IsEnd(); ++term){
                Result((*term),true);
            }
            EndOfRow();
            Hline();
        };

    }

};

#endif /* ifndef dyturbo_C */
