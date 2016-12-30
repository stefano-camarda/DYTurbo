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

#include "dyturbo.h"

#include "cubacall.h"
#include "clock_real.h"
#include "settings.h"
#include "interface.h"
#include "ctintegr.h"

#include "phasespace.h"
#include "rapint.h"
#include "HistoHandler.h"

using DYTurbo::PrintTable::Col3;
using DYTurbo::PrintTable::Col4;

bool DYTurbo::HasOnlyVegas = false;
bool DYTurbo::isDryRun = false;

namespace DYTurbo {

    // definition of data member
    Term subtotal;
    TermList ActiveTerms;
    BoundariesList ActiveBoundaries;
    BoundIterator last_bounds;

    // Term

    void Term::last_reset() {
        last_int.assign(opts.totpdf,0.);
        last_err2=0;
        last_time = clock_real();
    }

    void Term::RunIntegration(){
        double err;
        last_reset();
        /// @todo specialized (need to reincorporate)
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
        // save results to histograms
        if (!isVegas){
            for (size_t ivar = 0; ivar < last_int.size(); ++ivar) {
                HistoHandler::SetVariation(ivar);
                HistoHandler::FillResult(last_int[ivar],sqrt(last_err2));
            }
        }
        // cumulate
        total_time+=last_time;
        total_int+=last_int[0];
        total_err2+=last_err2;
        // cumulate integral to subtotal
        subtotal.last_int[0] += last_int[0];
        subtotal.total_int   += last_int[0];
        // cumulate error to subtotal
        subtotal.last_err2   += last_err2;
        subtotal.total_err2  += last_err2;
    }

    void Term::Print(){
        printf("%24s:%s",name.c_str(), description.c_str());
    }

    // Term iterator

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

    // Boundary iterator

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
        if (isFirst) return false;
        for (size_t i = 0; i < N_boundaries; ++i)
            if (current[i] != ActiveBoundaries[i].begin() )
                return false;
        return true;
    }

    BoundIterator & BoundIterator::operator++(){
        for (int i = N_boundaries-1; i >= 0; --i) { // start from last 
            current[i]++; // increase
            // check if it's equal to previous to last (need two numbers as boundaries)
            if(current[i]!=ActiveBoundaries[i].end()-1) break; // go to return
            else current[i]=ActiveBoundaries[i].begin(); // is previous to last 
        }
        /**
         * @attention After looping through all of boundary bins it points back
         * to beginning. Therefore it is necessary to set `isFirst=false` 
         */
        if (isFirst) isFirst=false;
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



    //! Internal flag only for test purposes
    bool TestAllTerms=false;

    /**
     * @brief Add term to active terms if is requested.
     *
     * Also it is checked if we are only using Vegas for integration. TYhi
     *
     * @param isActive If this false then skip all stuff.
     * @param fun The pointer to integration function, which should be called for calculation.
     * @param name Set the name of term
     * @param is_vegas For checking if we need Integration mode.
     *
     * @return It returns newly added term so we can define description by using stream operator.
     */
    Term & AddTermIfActive(const bool &isActive, void (* fun)(VecDbl &val,double &err), const String &name, const bool &is_vegas ){
        subtotal.description = "";
        if (!isActive && !TestAllTerms) return subtotal;
        ActiveTerms.push_back(Term());
        ActiveTerms.back().name = name;
        ActiveTerms.back().description = "";
        ActiveTerms.back().integrate = fun;
        ActiveTerms.back().isVegas = is_vegas;
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
        // finite born
        bool fixed_born = opts.doBORN && opts.fixedorder;
        if (fixed_born || TestAllTerms) {
            name="Fixed born";
            AddTermIfActive ( opts.bornint2d      , bornintegr2d   , name, isNotVegas ) << Col3( "cuhre (dm, dpt)" , "iter ="   , opts.niterBORN       );
            AddTermIfActive ( opts.bornintvegas4d , bornintegrMC4d , name, isVegas    ) << Col3( "vegas 4D"        , "ncalls =" , opts.vegasncallsBORN );
            AddTermIfActive ( opts.bornintvegas6d , bornintegrMC6d , name, isVegas    ) << Col3( "vegas 6D"        , "ncalls =" , opts.vegasncallsBORN );
        }
        // resummation
        bool resum_born = opts.doBORN && !opts.fixedorder;
        if (resum_born || TestAllTerms) {
            name="Resummation";
            AddTermIfActive ( opts.resint2d    , resintegr2d  , name , isNotVegas ) << Col3 ( "cuhre (dm, dpt)"     , "iter ="      , opts.niterBORN       )
                                                                                    << Col4 ( "","gauss (dy)"       , "nodes ="     , opts.yrule           )
                                                                                    << Col4 ( "",""                 , "intervals =" , opts.yintervals      );
            AddTermIfActive ( opts.resint3d    , resintegr3d  , name , isNotVegas ) << Col3 ( "cuhre (dm, dpt, dy)" , "iter ="      , opts.niterBORN       );
            AddTermIfActive ( opts.resintvegas , resintegrMC  , name , isVegas    ) << Col3 ( "vegas"               , "ncalls ="    , opts.vegasncallsBORN );
        }
        // CT
        if (opts.doCT || TestAllTerms) {
            name="Counter term";
            AddTermIfActive ( opts.ctint2d      , ctintegr2d , name , isNotVegas  ) << Col3 ( "cuhre (dm, dy)"  , "iter ="      , opts.niterCT )
                                                                                    << Col4 ( "","gauss (dpt)"  , "nodes ="     , 20           )
                                                                                    << Col4 ( "",""             , "intervals =" , 5            );
            AddTermIfActive ( opts.ctint3d      , ctintegr3d , name , isNotVegas  ) << Col3 ( "cuhre (dm, dpt, dy)" , "iter ="   , opts.niterCT       );
            AddTermIfActive ( opts.ctintvegas6d , ctintegrMC , name , isVegas     ) << Col3 ( "vegas 6D"            , "ncalls =" , opts.vegasncallsCT );
            AddTermIfActive ( opts.ctintvegas8d , ctintegr   , name , isVegas     ) << Col3 ( "vegas 8D"            , "ncalls =" , opts.vegasncallsCT );
        }
        // VJ finite
        bool vj_finite = opts.doVJ && !opts.doVJREAL && !opts.doVJVIRT;
        name="V+J LO";
        if (vj_finite || TestAllTerms){
            AddTermIfActive   ( opts.vjint3d                         , vjintegr3d   , name , isNotVegas )  << Col3 ( "cuhre (dm, dpt, dy)" , "iter ="   , opts.niterVJ );
            AddTermIfActive   ( opts.vjint5d && opts.order == 1      , vjlointegr5d , name , isNotVegas )  << Col3 ( "cuhre 5D???"         , "iter ="   , opts.niterVJ );
            AddTermIfActive   ( opts.vjintvegas7d && opts.order == 1 , vjlointegr7d , name , isVegas    )  << Col3 ( "vegas 7D"            , "ncalls =" , opts.vegasncallsVJLO ); //phase space improved MCFM integration
            //AddTermIfActive ( opts.vjintvegas7d && opts.order == 1 , vjlointegr   , name , isVegas    )  << Col3 ( "vegas 7D"            , "ncalls =" , opts.vegasncallsVJLO ); //original MCFM integration
        }
        // VJ NLO
        AddTermIfActive ( opts.doVJREAL  , vjrealintegr , "V+J Real"    , isVegas) << Col3 ( "vegas" , "ncalls =" , opts.vegasncallsVJREAL );
        AddTermIfActive ( opts.doVJVIRT  , vjvirtintegr , "V+J Virtual" , isVegas) << Col3 ( "vegas" , "ncalls =" , opts.vegasncallsVJVIRT );
    }

    /**
     * @brief Check for Integration mode and add boundary.
     *
     * If we have only Vegas calculations we can run with simple boundaries.
     * This means we will only pick up one bin with lowest and hightest value
     * as boundaries.
     *
     * @param ib Index of variable (\ref BoundaryIndex can be used as input)
     * @param name Name of boundary. Will be used for printing.
     * @param bins All boundaries requested from user.
     * @return Description of return.
     */
    void AddBoundary(size_t ib, string name, VecDbl &bins ){
        ActiveBoundaries[ib].name = name;
        ActiveBoundaries[ib].data.clear();
        if (HasOnlyVegas && !opts.force_binsampling){ // simple boundary mode
            ActiveBoundaries[ib].data.push_back(bins.front());
            ActiveBoundaries[ib].data.push_back(bins.back());
        } else { // integration mode
            ActiveBoundaries[ib].data = bins;
        }
    }

    /**
     * @brief Adding all boundaries
     *
     *  This have to be done after adding all terms (becaus of check for Vegas
     *  only).
     *
     *  If adding new boundary definition please also add new enum item into \ref BoundaryIndex.
     *
     * @param _inArg Description of param
     * @return Description of return.
     */
    void AddBoundaries(){
        ActiveBoundaries.resize(N_boundaries);
        VecDbl costh{opts.costhmin , opts.costhmax};
        AddBoundary ( b_M    , "q2"    , bins.mbins     ) ;
        AddBoundary ( b_Y    , "y"     , bins.ybins     ) ;
        AddBoundary ( b_QT   , "qT"    , bins.qtbins    ) ;
        AddBoundary ( b_CsTh , "costh" , costh );
    }

    /**
     * @brief First run of resummation and counter term.
     *
     * If using the DYRES approximation for PDFs, make sure that the PDF fit
     * is initialised in the same way Need to throw a random point according
     * to a breit wigner, which is used to determine xtauf in the PDF fit
     *
     */
    void WarmUpResummation(){
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
    }

    // forward declaration;
    namespace PrintTable{
        void IntegrationSettings();
    }

    void WarmUp(){
        /// - Create term list for calculation
        /// - Check if bounds should be simple (MC) or per each bin (cubature mode)
        AddTerms();
        AddBoundaries();
        /// - Check if all histograms are integrable if necessary (make warning)
        //HistoHandler::Book();
        /// - Check if we need to warm up CT integration or Resummation
        WarmUpResummation();
        /// - Clear subtotal
        subtotal = Term();
        subtotal.name = "TOTAL";
        subtotal.last_int.assign(opts.totpdf,0);
        if (!HasOnlyVegas) HistoHandler::DeleteNonIntegrableHists();
        PrintTable::IntegrationSettings();
    }


    //! Set current boundaries to phasespace.
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

    //! Close, clear, delete and say bye bye..
    void Terminate(){
        PrintTable::Footer();
        ActiveTerms.clear();
        ActiveBoundaries.clear();
        subtotal.last_reset();
        HistoHandler::Terminate();
    }

};

#endif /* ifndef dyturbo_C */
