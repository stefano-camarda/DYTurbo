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


TEST(DYTurbo,Initialization){
    int argc = 1;
    char *argv[] =  {
        ((char *) "DYTurbo_unittest"),
        ((char *) "input/test.in"),
    };
    DYTurbo::Init(argc,argv);
}

TEST(DYTurbo,BoundIteration){
    DYTurbo::WarmUp();
    for ( DYTurbo::BoundIterator it_bound; !it_bound.IsEnd(); ++it_bound) {
        it_bound.Print();
    }
}

TEST(DYTurbo,DISABLED_Integration){
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
