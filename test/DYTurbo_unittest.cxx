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

TEST(DYTurbo,Integration){
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

TEST(DYTurbo,Termination){
    DYTurbo::PrintTable::Footer();
    DYTurbo::Terminate();
}

#endif /* ifndef DYTurbo_unittest_CXX */
