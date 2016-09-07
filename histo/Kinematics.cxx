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
        // TODO: Implentation to DYTURBO: read current bin boundaries and set variables.
        //double qt_val = ( phasespace::qtmax + phasespace::qtmin )/2.;
        //double y_val  = ( phasespace::ymax  + phasespace::ymin  )/2.;
        //double Q_val  = ( phasespace::mmax  + phasespace::mmin  )/2.;
        // TODO: This needs to go cxx
        // BosPT qt;
        // qt.SetValue(qt_val); // set value and mark as "isCalculated"
    }


}


#endif /* ifndef Kinematics_CXX */
