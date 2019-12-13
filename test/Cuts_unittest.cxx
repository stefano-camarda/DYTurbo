#ifndef Cuts_unittest_CXX
#define Cuts_unittest_CXX
/**
 * @file Cuts_unittest.cxx
 * Description of cpp file
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-09-07
 */

#include "gtest/gtest.h"
#include "src/settings.h"
#include "histo/Kinematics.h"
#include "histo/KinematicCuts.h"

#include "old_cuts.h"
#include "old_cuts.C"

TEST(CutInterface,StandardCuts){
    double l1[] = {1., 1., 1., sqrt(3)};
    opts.makecuts = true;
    opts.lptcut = 20;
    opts.lycut = 1000;
    opts.lepptcut = 10;
    opts.alpptcut = 10;
    opts.lepycut = 1000;
    opts.alpycut = 1000;
    ASSERT_EQ(Kinematics::Cuts::SkipEvent, Kinematics::Cuts::KeepThisEvent(l1,l1));
    double l2[] = {100., 100., 100., 100.*sqrt(3)};
    ASSERT_EQ(Kinematics::Cuts::KeepEvent, Kinematics::Cuts::KeepThisEvent(l2,l2));
}


TEST(CutInterface,CrossCheckWithPrevious){
    opts.nproc=1;
    opts.lptcut = 20;
    opts.lycut = 2.5;
    double l1[] = { -38.674565, 2.399558, 378.431728, 380.410374};
    double l2[] = { 38.674565, -2.399558, 229.881660, 233.124554};
    ASSERT_EQ(cuts::lep(l1,l2), Kinematics::Cuts::KeepThisEvent(l1,l2));
}



#endif /* ifndef Cuts_unittest_CXX */
