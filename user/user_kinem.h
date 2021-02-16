#ifndef user_kinem_H
#define user_kinem_H
/**
 * @file user_kinem.h
 * User definition of kinematic observables.
 *
 * For example definition check \ref CosTh in `histo/KinematicDefinitions.h`.
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-09-07
 */

#include "KinematicDefinitions.h"

namespace Kinematics {
    // You can define custom observables and use them afterwards for custom
    // histograms or cuts.
    //
     NEWKIN ( BigAnswer ) { double calc(){ return 42;}  };
    //
    // Put your Observable definition here:
}

#endif /* ifndef user_kinem_H */
