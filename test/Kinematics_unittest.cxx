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

#include "old_cuts.C"
#include "old_kinem.C"

double s22 = 2.*sqrt(2.);
double s24 = 4.*sqrt(2.);

TEST(Kinemantic, CalculateOnlyOnce){
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


TEST(KinemanticInterface, CrossCheckWithPrevious){
    double l1[] = { -38.674565, 2.399558, 378.431728, 380.410374};
    double l2[] = { 38.674565, -2.399558, 229.881660, 233.124554};

    Kinematics::SetKinematics(l1,l2,1.0);

    kinematic::set(l1,l2);
    kinematic::calc_vb();
    kinematic::calc_angles();


   ASSERT_DOUBLE_EQ(kinematic::lp  [0] , Kinematics::ALpPX()() );
   ASSERT_DOUBLE_EQ(kinematic::lm  [0] , Kinematics::LepPX()() );
   ASSERT_DOUBLE_EQ(kinematic::v   [0] , Kinematics::BosPX()() );

   ASSERT_DOUBLE_EQ(kinematic::lp  [1] , Kinematics::ALpPY()() );
   ASSERT_DOUBLE_EQ(kinematic::lm  [1] , Kinematics::LepPY()() );
   ASSERT_DOUBLE_EQ(kinematic::v   [1] , Kinematics::BosPY()() );

   ASSERT_DOUBLE_EQ(kinematic::lp  [2] , Kinematics::ALpPZ()() );
   ASSERT_DOUBLE_EQ(kinematic::lm  [2] , Kinematics::LepPZ()() );
   ASSERT_DOUBLE_EQ(kinematic::v   [2] , Kinematics::BosPZ()() );

   ASSERT_DOUBLE_EQ(kinematic::lp  [3] , Kinematics::ALpE()() );
   ASSERT_DOUBLE_EQ(kinematic::lm  [3] , Kinematics::LepE()() );
   ASSERT_DOUBLE_EQ(kinematic::v   [3] , Kinematics::BosE()() );

   ASSERT_DOUBLE_EQ(cuts::getY(kinematic::v) , Kinematics::BosY()() );
   ASSERT_DOUBLE_EQ(cuts::getY(l1)           , Kinematics::LepRap()() );
   ASSERT_DOUBLE_EQ(cuts::getY(l2)           , Kinematics::ALpRap()() );

   ASSERT_DOUBLE_EQ(kinematic::m2      , Kinematics::BosM2()() );
   ASSERT_DOUBLE_EQ(kinematic::m       , Kinematics::BosM()() );
   ASSERT_DOUBLE_EQ(kinematic::qt2     , Kinematics::BosPT2()() );
   ASSERT_DOUBLE_EQ(kinematic::qt      , Kinematics::BosPT()() );
   ASSERT_DOUBLE_EQ(kinematic::y       , Kinematics::BosY()() );
   ASSERT_DOUBLE_EQ(kinematic::phiZ    , Kinematics::BosPhi()() );
   ASSERT_DOUBLE_EQ(kinematic::costh   , Kinematics::CosThCS()() );
   ASSERT_DOUBLE_EQ(kinematic::phi_lep , Kinematics::PhiCS()() );
   ASSERT_DOUBLE_EQ(kinematic::mtrans  , Kinematics::BosMT()() );
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
