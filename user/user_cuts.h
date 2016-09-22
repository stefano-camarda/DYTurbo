#ifndef user_cuts_H
#define user_cuts_H
/**
 * @file user_cuts.h
 * User definitions of fiducial cuts.
 *
 * @brief User definitions will be included in `histo/KinematicCuts.h`.
 * For example of definition check `histo/KinematicCuts.h` file.
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-09-07
 */

#include "histo/Kinematics.h"
#include "histo/KinematicCuts.h"

namespace Kinematics {
    namespace Cuts {
        struct UserCuts : public CutBase {

            // define variables here (e.g. Observables, cut values)
            // bool IDontLikeEvent = false;

            bool operator()(){
                // put your cut here:
                // if (IDontLikeEvent) return SkipEvent;
                return KeepEvent;
            }
        } user_cut;
    }
}

#endif /* ifndef user_cuts_H */
