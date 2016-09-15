#ifndef Kinematics_CXX
#define Kinematics_CXX
/**
 * @file Kinematics.cxx
 * Description of this macro
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-08-26
 */

#include "Kinematics.h"

namespace Kinematics{

    // input values definitions
    double p3[4];
    double p4[4];
    double event_weight;

    // all flags are stored here
    VariableFlags flags;

    void SetKinematics(double _p3[4], double _p4[4], double wgt){
        for (auto i=0; i<4;i++){
            p3[i]=_p3[i];
            p4[i]=_p4[i];
        }
        event_weight=wgt;
        // Set all flags to uncalculated
        for (auto ptr_isCalculated : flags) *ptr_isCalculated=false;
    }

    void SetMiddlePoint(){
        // NOTE: There are only few observables which can be evaluated in
        // integrator mode. We need to receive middle point of integration and recalculate
        isIntegratorMode=true;
        for (auto ptr_isCalculated : flags) *ptr_isCalculated=false;
    }

    bool isIntegratorMode=false;


}


#endif /* ifndef Kinematics_CXX */
