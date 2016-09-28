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
#include "gmock/gmock.h"

#include "histo/Kinematics.h"
#include "histo/KinematicDefinitions.h"
#include "phasespace/phasespace.h"

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


TEST(KinemanticInterface, IntegratorVariableIsCalculated){
    double l1[] = {0.,2.,0.,2.};
    Kinematics::SetKinematics(l1, l1, 1.0 );
    phasespace::setbounds(0.,2., 2.,4., 4.,6.);
    phasespace::setcthbounds(-1.,0.);
    Kinematics::BosPT pt;
    Kinematics::BosY y;
    Kinematics::BosM m;
    Kinematics::CosThCS costh;
    Kinematics::SetMiddlePoint();
    // Check that integrable variable is calculated
    ASSERT_DOUBLE_EQ(1., m());
    ASSERT_DOUBLE_EQ(3., pt());
    ASSERT_DOUBLE_EQ(5., y());
    ASSERT_DOUBLE_EQ(-.5, costh());
}

TEST(KinemanticInterface, NonIntegratorVariableIsNotCalculated){
    double l1[] = {0.,2.,0.,2.};
    Kinematics::SetKinematics(l1, l1, 1.0 );
    phasespace::setbounds(0.,2., 2.,4., 4.,6.);
    Kinematics::SetMiddlePoint();
    // Check that non-integrable variable is not calculated
    Kinematics::LepPX lepPX;
    ASSERT_DOUBLE_EQ(666., lepPX());
}


#endif /* ifndef Kinematics_unittest_CXX */
