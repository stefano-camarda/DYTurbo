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

TEST(CutInterface,StandartCuts){
    double l1[] = {1., 1., 1., sqrt(3)};
    double l2[] = {1., 1., 1., sqrt(3)};
    opts.makelepcuts = true;
    opts.l1ptcut = 10;
    opts.l2ptcut = 10;
    ASSERT_EQ(Kinematics::Cuts::SkipEvent, Kinematics::Cuts::KeepThisEvent(l1,l2));
}

#endif /* ifndef Cuts_unittest_CXX */
