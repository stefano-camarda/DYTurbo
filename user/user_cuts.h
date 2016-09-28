#ifndef user_cuts_H
#define user_cuts_H
/**
 * @file user_cuts.h
 *
 * @brief User definitions of cuts.
 *
 * For example of definition see class \ref StandartCuts in `histo/Kinematics.cxx`
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-09-07
 */

#include "histo/Kinematics.h"
#include "histo/KinematicCuts.h"
#include "user_kinem.h"

namespace Kinematics {
    namespace Cuts {
        struct UserCuts : public CutBase {

            // You can define variables here (e.g. Observables, cut values, etc.)
            //
             bool IDontLikeEvent = false;
             BigAnswer answer;

            bool operator()(){
                // Fastest way to cut is to return SkipEvent after each
                // condition:
                //
                 if (IDontLikeEvent) return SkipEvent;
                 if (answer()!=42) return SkipEvent;
                //
                // Put your cut here:
                return KeepEvent;
            }
        };
    }
}

#endif /* ifndef user_cuts_H */
