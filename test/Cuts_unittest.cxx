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

#include "src/old_cuts.h"
//#include "src/old_cuts.C"

TEST(CutInterface,StandardCuts){
    double l1[] = {1., 1., 1., sqrt(3)};
    opts.makecuts = true;
    opts.lptcut = 0;
    opts.lycut = 1000;
    opts.l1ptcut = 10;
    opts.l2ptcut = 10;
    opts.l1ycut = 1000;
    opts.l2ycut = 1000;
    ASSERT_EQ(Kinematics::Cuts::SkipEvent, Kinematics::Cuts::KeepThisEvent(l1,l1));
    double l2[] = {100., 100., 100., 100.*sqrt(3)};
    ASSERT_EQ(Kinematics::Cuts::KeepEvent, Kinematics::Cuts::KeepThisEvent(l2,l2));
}


TEST(CutInterface,CrossCheckWithPrevious){
    double l1[] = {1., 1., 1., sqrt(3)};
    ASSERT_EQ(cuts::lep(l1,l1), Kinematics::Cuts::KeepThisEvent(l1,l1));
    double l2[] = {100., 100., 100., 100.*sqrt(3)};
    ASSERT_EQ(cuts::lep(l2,l2), Kinematics::Cuts::KeepThisEvent(l2,l2));
}



#endif /* ifndef Cuts_unittest_CXX */
