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
#include "clock_real.h"

#include "dyturbo.h"

#include<iostream>
#include<iomanip>
using std::setw;
using std::cout;
using std::endl;

bool DYTurbo::HasOnlyVegas = false;
bool DYTurbo::isDryRun = false;

namespace DYTurbo {

    Term subtotal;

    // Boundaries
    struct Boundaries{
        std::string name ="";
        VecDbl data;
        inline BoundariesListItr begin() { return data.begin();};
        inline BoundariesListItr end()   { return data.end(); };
        inline double front() { return data.front();};
        inline double back()  { return data.back(); };
        inline size_t size()  { return data.size(); };
    };

    BoundariesList ActiveBoundaries;

    enum BoundaryIndex {
        b_M=0,
        b_Y,
        b_QT,
        b_CsTh,
        //b_Phi,
        N_boundaries
    };


    BoundIterator::BoundIterator(){
        // set first boundary
        isFirst=true;
        current.clear();
        if (!ActiveBoundaries.empty())
            for (size_t i = 0; i < N_boundaries; ++i) {
                current.push_back(ActiveBoundaries[i].begin());
            }
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
        printf ( "| M    : %f -- %f " , *current [ b_M    ] , * ( current [ b_M    ] +1 ) ) ;
        printf ( "| Qt   : %f -- %f " , *current [ b_QT   ] , * ( current [ b_QT   ] +1 ) ) ;
        printf ( "| Y    : %f -- %f " , *current [ b_Y    ] , * ( current [ b_Y    ] +1 ) ) ;
        printf ( "| CsTh : %f -- %f " , *current [ b_CsTh ] , * ( current [ b_CsTh ] +1 ) ) ;
        printf ( "\n");
    }

    void Term::last_reset() {
        last_int.assign(opts.totpdf,0.);
        last_err2=0;
        last_time = clock_real();
    }


    void Term::RunIntegration(){
        double err;
        last_reset();
        // TODO: specialized (need to reincorporate)
        if ( integrate == resintegr2d ){
            if (opts.resumcpp) rapint::cache(phasespace::ymin, phasespace::ymax);
            else cacheyrapint_(phasespace::ymin, phasespace::ymax);
        }
        // run
        if (isDryRun){
            // this is for testing interface
            last_int[0] = 1; err = 1;
        } else {
            integrate(last_int,err);
        }
        //
        last_time = clock_real()-last_time;
        last_err2 += err*err;
        // cumulate
        total_time+=last_time;
        total_int+=last_int[0];
        total_err2+=last_err2;
        // cumulate time to subtotal
        subtotal.last_time   += last_time;
        subtotal.total_time  += last_time;
        // cumulate integral to subtotal
        subtotal.last_int[0] += last_int[0];
        subtotal.total_int   += last_int[0];
        // cumulate error to subtotal
        subtotal.last_err2   += last_err2;
        subtotal.total_err2  += last_err2;
    }
    void Term::Print(){
        cout << setw(24)  << name.c_str() << ":";
        cout << description.c_str();
    }
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


    namespace PrintTable{
        struct Col4 {
            String data;
            template<class S1, class S2, class S3, class S4 > Col4( S1 col1, S2 col2, S3 col3, S4 col4){
                SStream tmp;
                tmp << setw(25) << col1;
                tmp << setw(20) << col2;
                tmp << setw(12) << col3;
                tmp << setw(12) << col4;
                tmp << '\n';
                data = tmp.str();
            }
        };
        inline std::ostream & operator<< (std::ostream & strm, const Col4 &col){ strm << col.data; return strm; }

        struct Col3{
            String data;
            template<class S2, class S3, class S4 > Col3(S2 col2, S3 col3, S4 col4){
                SStream tmp;
                tmp << setw(20) << col2;
                tmp << setw(12) << col3;
                tmp << setw(12) << col4;
                tmp << '\n';
                data = tmp.str();
            }
        };
        inline std::ostream & operator<< (std::ostream & strm, const Col3 &col){ strm << col.data; return strm; }
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

    // this is hidden only for test purposes
    bool TestAllTerms=false;

    Term & AddTermIfActive(const bool &isActive, void (* fun)(VecDbl &val,double &err), const String &name, const bool &is_vegas ){
        subtotal.description = "";
        if (!isActive && !TestAllTerms) return subtotal;
        ActiveTerms.push_back(Term());
        ActiveTerms.back().name = name;
        ActiveTerms.back().description = "";
        ActiveTerms.back().integrate = fun;
        HasOnlyVegas&=is_vegas;
        return ActiveTerms.back();
    }

    template<class Streamable> Term & Term::operator<<(const Streamable &data){
        SStream strm;
        strm << data;
        description += strm.str();
        return (*this);
    }


    void AddTerms(){
        ActiveTerms.clear();
        HasOnlyVegas=true;
        const bool isVegas=true;
        const bool isNotVegas=false;
        int w0=25;
        int w1=30;
        int w2=12;
        int w3=12;
        string name;
        // born
        bool fixed_born = opts.doBORN && opts.fixedorder;
        if (fixed_born || TestAllTerms) {
            name="Fixed born";
            AddTermIfActive ( opts.bornint2d      , bornintegr2d   , name, isNotVegas ) << PrintTable::Col3( "cuhre (dm, dpt)" , "iter ="   , opts.niterBORN       );
            AddTermIfActive ( opts.bornintvegas4d , bornintegrMC4d , name, isVegas    ) << PrintTable::Col3( "vegas 4D"        , "ncalls =" , opts.vegasncallsBORN );
            AddTermIfActive ( opts.bornintvegas6d , bornintegrMC6d , name, isVegas    ) << PrintTable::Col3( "vegas 6D"        , "ncalls =" , opts.vegasncallsBORN );
        }
        // resummation
        bool resum_born = opts.doBORN && !opts.fixedorder;
        if (resum_born || TestAllTerms) {
            name="Resummation";
            AddTermIfActive ( opts.resint2d    , resintegr2d  , name , isNotVegas ) << PrintTable::Col3 ( "cuhre (dm, dpt)"     , "iter ="      , opts.niterBORN       )
                                                                                    << PrintTable::Col4 ( "","gauss (dy)"       , "nodes ="     , opts.yrule           )
                                                                                    << PrintTable::Col4 ( "",""                 , "intervals =" , opts.yintervals      );
            AddTermIfActive ( opts.resint3d    , resintegr3d  , name , isNotVegas ) << PrintTable::Col3 ( "cuhre (dm, dpt, dy)" , "iter ="      , opts.niterBORN       );
            AddTermIfActive ( opts.resintvegas , resintegrMC  , name , isVegas    ) << PrintTable::Col3 ( "vegas"               , "ncalls ="    , opts.vegasncallsBORN );
        }
        // CT
        if (opts.doCT || TestAllTerms) {
            name="Counter term";
            AddTermIfActive ( opts.ctint2d      , ctintegr2d , name , isNotVegas  ) << PrintTable::Col3 ( "cuhre (dm, dy)"  , "iter ="      , opts.niterCT )
                                                                                    << PrintTable::Col4 ( "","gauss (dpt)"  , "nodes ="     , 20           )
                                                                                    << PrintTable::Col4 ( "",""             , "intervals =" , 5            );
            AddTermIfActive ( opts.ctint3d      , ctintegr3d , name , isNotVegas  ) << PrintTable::Col3 ( "cuhre (dm, dpt, dy)" , "iter ="   , opts.niterCT       );
            AddTermIfActive ( opts.ctintvegas6d , ctintegrMC , name , isVegas     ) << PrintTable::Col3 ( "vegas 6D"            , "ncalls =" , opts.vegasncallsCT );
            AddTermIfActive ( opts.ctintvegas8d , ctintegr   , name , isVegas     ) << PrintTable::Col3 ( "vegas 8D"            , "ncalls =" , opts.vegasncallsCT );
        }
        // VJ finite
        bool vj_finite = opts.doVJ && !opts.doVJREAL && !opts.doVJVIRT;
        name="V+J LO";
        if (vj_finite || TestAllTerms){
            AddTermIfActive   ( opts.vjint3d                         , vjintegr3d   , name , isNotVegas )  << PrintTable::Col3 ( "cuhre (dm, dpt, dy)" , "iter ="   , opts.niterVJ         );
            AddTermIfActive   ( opts.vjint5d && opts.order == 1      , vjlointegr5d , name , isVegas    )  << PrintTable::Col3 ( "vegas 5D"            , "ncalls =" , opts.vegasncallsVJLO );
            AddTermIfActive   ( opts.vjintvegas7d && opts.order == 1 , vjlointegr7d , name , isVegas    )  << PrintTable::Col3 ( "vegas 7D"            , "ncalls =" , opts.vegasncallsVJLO );
            //AddTermIfActive ( opts.vjintvegas7d && opts.order == 1 , vjlointegr   , name , isVegas    )  << PrintTable::Col3 ( "vegas 7D"            , "ncalls =" , opts.vegasncallsVJLO );
        }
        // VJ NLO
        AddTermIfActive ( opts.doVJREAL  , vjrealintegr , "V+J Real"    , isVegas) << PrintTable::Col3 ( "vegas" , "ncalls =" , opts.vegasncallsVJREAL );
        AddTermIfActive ( opts.doVJVIRT  , vjvirtintegr , "V+J Virtual" , isVegas) << PrintTable::Col3 ( "vegas" , "ncalls =" , opts.vegasncallsVJVIRT );
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
        VecDbl costh{opts.costhmin , opts.costhmax};
        AddBoundary ( b_M    , "q2"    , bins.mbins     ) ;
        AddBoundary ( b_Y    , "y"     , bins.ybins     ) ;
        AddBoundary ( b_QT   , "qT"    , bins.qtbins    ) ;
        AddBoundary ( b_CsTh , "costh" , costh );
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
        subtotal = Term();
        subtotal.name = "TOTAL";
        subtotal.last_int.assign(opts.totpdf,0);
        PrintTable::IntegrationSettings();
    }

    BoundIterator last_bounds;

    void SetBounds(BoundIterator bounds){
        phasespace::setcthbounds(
                bounds.loBound ( b_CsTh  )  , bounds.hiBound ( b_CsTh  ) 
                );
        phasespace::setbounds(
                bounds.loBound ( b_M  )  , bounds.hiBound ( b_M  )  ,
                bounds.loBound ( b_QT )  , bounds.hiBound ( b_QT )  ,
                bounds.loBound ( b_Y  )  , bounds.hiBound ( b_Y  )  );
        last_bounds = bounds;
    }

    void Terminate(){
        ActiveTerms.clear();
    }

    namespace PrintTable {

        size_t eoc=2; // bound length
        size_t wb=5; // bound length
        size_t wr=5; // result length
        size_t bound_width = 2*wb+3;
        size_t time_width = 0;
        size_t term_width = 2*wr+4 + time_width;

        inline void EndOfCell()  { cout << " |" << flush; }
        inline void BeginOfRow() { cout << "|"  << flush; }
        inline void EndOfRow()   { cout << endl;          }
        void Hline(){
            size_t N_loopingBounds =0;
            for (auto &a: ActiveBoundaries ) if (a.size()>2) N_loopingBounds++;
            size_t wtotal = 1; // begin of line
            wtotal += (bound_width +eoc)*N_loopingBounds;
            wtotal += (term_width  +eoc)*(ActiveTerms.size() + 1);
            string hline ( wtotal, '-');
            cout << hline << endl;
        };

        void BoundsAllLooping(bool printNames=false, bool useFullBound=false){
            // loop over bounds
            for (size_t ibound=0; ibound < N_boundaries; ++ibound){
                if (ActiveBoundaries[ibound].size()<3) continue;
                if (printNames){
                    cout << setw(bound_width) << ActiveBoundaries[ibound].name ;
                } else {
                    double lo = useFullBound ? ActiveBoundaries[ibound].front() : last_bounds.loBound(ibound);
                    double hi = useFullBound ? ActiveBoundaries[ibound].back()  : last_bounds.hiBound(ibound);
                    cout << setw(wb) << lo ;
                    cout << " - ";
                    cout << setw(wb) << hi ;
                }
                EndOfCell();
            }
        }

        void IntegrationSettings(){
            bool fixed_born = opts.doBORN && opts.fixedorder;
            bool fixed_born_cubature = fixed_born && opts.bornint2d ;
            bool resum_born = opts.doBORN && !opts.fixedorder;
            bool resum_born_cubature = resum_born && (opts.resint3d || opts.resint2d);
            bool ct_cubature = opts.doCT && (opts.ctint2d || opts.ctint3d);
            cout << endl;
            cout << "========================  Integration settings =======================" << endl;
            cout << endl;
            cout << Col4 ( "Random seeds"   , opts.rseed     , "" , "" ) ;
            cout << Col4 ( "Parallel cores" , opts.cubacores , "" , "" ) ;
            cout << endl;
            for (TermIterator iterm;!iterm.IsEnd();++iterm){
                (*iterm).Print();
            }
            cout << endl;
            if (resum_born && opts_.approxpdf_==0){
                cout << Col4 ( "Mellin inverse transform:" , "gaussian "      , "nodes ="     , opts.mellinrule      ) ;
                cout << Col4 ( ""                          , ""               , "intervals =" , opts.mellinintervals ) ;
                cout << Col4 ( ""                          , "imaginary axis" , "zmax ="      , opts.zmax            ) ;
            }
            if ( resum_born_cubature || ct_cubature || fixed_born_cubature ){
                if (opts.cubaint)
                    cout << Col4( "Angular variables:" , "Suave in dcosth,dphi" , "ncalls =" , opts.suavepoints );
                else if (opts.quadint) {
                    cout << Col4( "Angular variable costh:" , "semi-analytical", "ncstart ="   , opts.ncstart      );
                    cout << Col4( "Angular variables phi:"  , "gaussian"       , "intervals =" , opts.phiintervals );
                }
                else if (opts.trapezint) {
                    cout << Col4( "Angular variables costh:" , "semi-analytical" , "ncstart =" , opts.ncstart   );
                    cout << Col4( "Angular variables phi:"   , "trapezoidal"     , "points ="  , opts.nphitrape );
                }
            }
            cout << endl;
            string col1 = "Constant boundaries";
            for (size_t ibound = 0; ibound < N_boundaries; ++ibound){
                if (ActiveBoundaries[ibound].size() > 2) continue;
                string col2 = ActiveBoundaries[ibound].name;
                cout << Col4 ( col1 , col2 , "low ="  , ActiveBoundaries [ ibound ] .front (  )  ) ;
                cout << Col4 ( ""   , ""   , "high =" , ActiveBoundaries [ ibound ] .back  (  )  ) ;
                col1 = "";
            }
            cout << endl;
        }

        inline void TermName(Term &term){
            cout << setw(term_width) << term.name;
            EndOfCell();
        }

        void Header() {
            Hline();
            BeginOfRow();
            BoundsAllLooping(true); // print only names
            for( TermIterator term; !term.IsEnd(); ++term){
                TermName(*term);
            }
            TermName(subtotal);
            EndOfRow();
            Hline();
        };

        void Footer() { Hline(); };


        void Bounds(bool use_full_bound) {
            BeginOfRow();
            BoundsAllLooping(false,use_full_bound);
        };

        void Result(const Term &term, bool printGrandTotal) {
            double val  = ( printGrandTotal ) ? term.total_int  : term.last_int[0] ;
            double err2 = ( printGrandTotal ) ? term.total_err2 : term.last_err2   ;
            cout << setw(wr) <<  val;
            cout << " +- ";
            cout << setw(wr) <<  sqrt(err2);
            EndOfCell();
        };

        void ResultSubTotal(bool is_grandtotal) {
            Result(subtotal,is_grandtotal);
            EndOfRow();
            subtotal.last_reset();
        };


        void ResultGrandTotal() {
            Hline();
            Bounds(true);
            for( TermIterator term; !term.IsEnd(); ++term){
                Result((*term),true);
            }
            ResultSubTotal(true);
            Hline();
        };

    }

};

#endif /* ifndef dyturbo_C */
