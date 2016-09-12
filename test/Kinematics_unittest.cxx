#ifndef Kinematics_unittest_CXX
#define Kinematics_unittest_CXX
/**
 * @file Kinematics_unittest.cxx
 * Description of cpp file
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-09-07
 */


#include "gtest/gtest.h"

#include "Kinematics.h"
#include "KinematicDefinitions.h"

TEST(Kinemantic, CalculateOnlyOnce){
    double s22 = 2.*sqrt(2.);
    double s24 = 4.*sqrt(2.);
    Kinematics::BosPT pt1;
    Kinematics::BosPT pt2;
    ASSERT_FALSE( pt1.IsCalculated() );
    ASSERT_FALSE( pt2.IsCalculated() );

    // check calculated only once
    double l1[] = {0.,2.,0.,2.};
    double l2[] = {2.,0.,0.,2.};
    Kinematics::SetKinematics( l1,l2, 1.0);
    ASSERT_FALSE( pt1.IsCalculated() );
    ASSERT_FALSE( pt2.IsCalculated() );
    ASSERT_DOUBLE_EQ(s22, pt1());
    ASSERT_TRUE( pt1.IsCalculated() );
    ASSERT_TRUE( pt2.IsCalculated() );
    ASSERT_DOUBLE_EQ(pt2(), pt1());

    // check reset on new kinematics
    l1[0]=2.0;
    l2[1]=2.0;
    l1[3]=s22;
    l2[3]=s22;
    Kinematics::SetKinematics(l1, l2, 1.0 );
    ASSERT_FALSE( pt1.IsCalculated() );
    ASSERT_FALSE( pt2.IsCalculated() );
    ASSERT_DOUBLE_EQ(s24, pt1());
    ASSERT_TRUE( pt1.IsCalculated() );
    ASSERT_TRUE( pt2.IsCalculated() );
    ASSERT_DOUBLE_EQ(pt2(), pt1());
}

TEST(Kinemantic, SetMiddlePoint){
}

TEST(Kinemantic, CheckWithOldKinematics){
}

#endif /* ifndef Kinematics_unittest_CXX */
