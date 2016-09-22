#ifndef DYTurbo_unittest_CXX
#define DYTurbo_unittest_CXX
/**
 * @file DYTurbo_unittest.cxx
 * Description of cpp file
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-09-16
 */


#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "src/dyturbo.h"
#include "src/settings.h"
#include "src/cubacall.h"

typedef std::vector<string> VecStr;
typedef std::vector<double> VecDbl;
typedef std::stringstream SStream;



TEST(DYTurbo,Initialization){
    int argc = 2;
    char *argv[] =  {
        ((char *) "DYTurbo_unittest"),
        ((char *) "../input/test.in"),
    };
    DYTurbo::Init(argc,argv);
}

TEST(DYTurbo,TermIteration){
    DYTurbo::WarmUp();
    std::vector <void (*)(std::vector<double>&, double&)> term_funs;
    term_funs.push_back(resintegr2d);
    term_funs.push_back(ctintegr2d);
    term_funs.push_back(vjrealintegr);
    term_funs.push_back(vjvirtintegr);
    ASSERT_EQ(term_funs.size(), DYTurbo::ActiveTerms.size() )
        << "Incorrect number of terms. Check your input file." ;
    for ( DYTurbo::TermIterator it_term; !it_term.IsEnd(); ++it_term) {
        ASSERT_EQ( term_funs[it_term.icurrent], (*it_term).integrate )
             << "Wrong term in active list index: " << it_term.icurrent << endl
             << "name: " << (*it_term).name << endl
             << "decr: " << (*it_term).description << endl
             << "Check your inputfile." << endl;
    }
}

TEST(DYTurbo,BoundIteration){
    size_t expected;
    size_t counter;
    // simple binner
    bins.qtbins = {0.,30.};
    bins.ybins  = {-3.,3.};
    bins.mbins  = {50.,500. };
    expected = (bins.qtbins.size() - 1 ) * (bins.ybins.size() - 1 ) * (bins.mbins.size() - 1 ) ;
    counter = 0 ;
    DYTurbo::WarmUp();
    for ( DYTurbo::BoundIterator it_bound; !it_bound.IsEnd(); ++it_bound) {
        counter ++;
    }
    ASSERT_EQ ( expected , counter  ) << " Incorrect number of bins from boundaries.";
    // every step
    bins.qtbins = {0., 10., 20., 30.};
    bins.ybins  = {-3., 0. ,3.};
    expected = (bins.qtbins.size() - 1 ) * (bins.ybins.size() - 1 ) * (bins.mbins.size() - 1 ) ;
    counter = 0 ;
    DYTurbo::WarmUp();
    for ( DYTurbo::BoundIterator it_bound; !it_bound.IsEnd(); ++it_bound) {
        counter ++;
    }
    ASSERT_EQ ( expected , counter  ) << " Incorrect number of bins from boundaries.";
}

TEST(DYTurbo,DryLooping){
    DYTurbo::isDryRun = true;
    DYTurbo::WarmUp();
    DYTurbo::PrintTable::Header();
    for ( DYTurbo::BoundIterator bounds; !bounds.IsEnd(); ++bounds) {
        DYTurbo::SetBounds(bounds);
        DYTurbo::PrintTable::Bounds();
        for (DYTurbo::TermIterator term; !term.IsEnd(); ++term){
            (*term).RunIntegration();
            DYTurbo::PrintTable::Result((*term));
        }
        DYTurbo::PrintTable::ResultSubTotal();
    }
    DYTurbo::PrintTable::ResultGrandTotal();
}

namespace DYTurbo{ extern bool TestAllTerms; }
struct ResTab {
    int ord;
    void (* fun)(VecDbl &val,double &err);
    double integ;
    double unc;
};
typedef vector<ResTab> VecResTab;

VecResTab results_list = { 
    { 0 , bornintegr2d   , 1.000000 , 1.000000} ,
    { 0 , bornintegrMC4d , 1.000000 , 1.000000} ,
    { 0 , bornintegrMC6d , 1.000000 , 1.000000} ,
    { 0 , resintegr2d    , 1.000000 , 1.000000} ,
    { 0 , resintegr3d    , 1.000000 , 1.000000} ,
    { 0 , resintegrMC    , 1.000000 , 1.000000} ,
    { 0 , ctintegr2d     , 1.000000 , 1.000000} ,
    { 0 , ctintegr3d     , 1.000000 , 1.000000} ,
    { 0 , ctintegrMC     , 1.000000 , 1.000000} ,
    { 0 , ctintegr       , 1.000000 , 1.000000} ,
    { 0 , vjintegr3d     , 1.000000 , 1.000000} ,
    { 0 , vjlointegr5d   , 1.000000 , 1.000000} ,
    { 0 , vjlointegr7d   , 1.000000 , 1.000000} ,
    { 0 , vjrealintegr   , 1.000000 , 1.000000} ,
    { 0 , vjvirtintegr   , 1.000000 , 1.000000} ,
    { 1 , bornintegr2d   , 1.000000 , 1.000000} ,
    { 1 , bornintegrMC4d , 1.000000 , 1.000000} ,
    { 1 , bornintegrMC6d , 1.000000 , 1.000000} ,
    { 1 , resintegr2d    , 1.000000 , 1.000000} ,
    { 1 , resintegr3d    , 1.000000 , 1.000000} ,
    { 1 , resintegrMC    , 1.000000 , 1.000000} ,
    { 1 , ctintegr2d     , 1.000000 , 1.000000} ,
    { 1 , ctintegr3d     , 1.000000 , 1.000000} ,
    { 1 , ctintegrMC     , 1.000000 , 1.000000} ,
    { 1 , ctintegr       , 1.000000 , 1.000000} ,
    { 1 , vjintegr3d     , 1.000000 , 1.000000} ,
    { 1 , vjlointegr5d   , 1.000000 , 1.000000} ,
    { 1 , vjlointegr7d   , 1.000000 , 1.000000} ,
    { 1 , vjrealintegr   , 1.000000 , 1.000000} ,
    { 1 , vjvirtintegr   , 1.000000 , 1.000000} ,
    { 2 , bornintegr2d   , 1.000000 , 1.000000} ,
    { 2 , bornintegrMC4d , 1.000000 , 1.000000} ,
    { 2 , bornintegrMC6d , 1.000000 , 1.000000} ,
    { 2 , resintegr2d    , 1.000000 , 1.000000} ,
    { 2 , resintegr3d    , 1.000000 , 1.000000} ,
    { 2 , resintegrMC    , 1.000000 , 1.000000} ,
    { 2 , ctintegr2d     , 1.000000 , 1.000000} ,
    { 2 , ctintegr3d     , 1.000000 , 1.000000} ,
    { 2 , ctintegrMC     , 1.000000 , 1.000000} ,
    { 2 , ctintegr       , 1.000000 , 1.000000} ,
    { 2 , vjintegr3d     , 1.000000 , 1.000000} ,
    { 2 , vjlointegr5d   , 1.000000 , 1.000000} ,
    { 2 , vjlointegr7d   , 1.000000 , 1.000000} ,
    { 2 , vjrealintegr   , 1.000000 , 1.000000} ,
    { 2 , vjvirtintegr   , 1.000000 , 1.000000} ,
};

const char * print_term(int &ord, DYTurbo::TermIterator &iterm){
    SStream strm;
    strm << endl;
    strm << "name: " << (*iterm).name.c_str() << endl;
    strm << "ord: " << ord << endl;
    strm << "result: " << (*iterm).last_int[0] << " +- " << sqrt((*iterm).last_err2) << endl;
    strm << "Description: " << (*iterm).description.c_str() << endl;
    return strm.str().c_str();
    
}

void CheckResult(int &ord, DYTurbo::TermIterator &iterm){
    size_t ires= ord*DYTurbo::ActiveTerms.size();
    ires+= iterm.icurrent;
    ASSERT_EQ(ord, results_list[ires].ord) 
        << "Wrong order! Expected " << ord << "But have this:" << print_term(ord,iterm);
    ASSERT_EQ((*iterm).integrate, results_list[ires].fun) 
        << "Wrong function! But have this:" << print_term(ord,iterm);
}



TEST(DYTurbo,ResultEveryTerm){
    // turn off dryrun
    DYTurbo::isDryRun = true;
    // turn on all terms 
    DYTurbo::TestAllTerms = true;
    // set boundaries
    bins.qtbins = {10. , 30.  };
    bins.ybins  = {0.  , 1.   };
    bins.mbins  = {50. , 100. };
    opts.costhmin=-1;
    opts.costhmax=+1;
    // warm up
    DYTurbo::WarmUp();
    DYTurbo::BoundIterator bound;
    DYTurbo::SetBounds(bound);
    // for all orders and for all terms
    for (int ord = 0; ord < 3; ++ord) {
        for (DYTurbo::TermIterator iterm; !iterm.IsEnd(); ++iterm ){
            // run
            (*iterm).RunIntegration();
            // check according to function and order
            CheckResult(ord,iterm);
        }
    }
}

TEST(DYTurbo,Termination){
    DYTurbo::PrintTable::Footer();
    DYTurbo::Terminate();
}

#endif /* ifndef DYTurbo_unittest_CXX */
