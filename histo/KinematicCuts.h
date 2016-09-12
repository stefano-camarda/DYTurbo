#ifndef KinematicCuts_H
#define KinematicCuts_H
/**
 * @file KinematicCuts.h
 * Description of cpp file
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-09-07
 */


namespace Kinematics {
    template<> struct Cut<"Lep:_pT>15_|eta|<2.5"> {
        BosCh ch;
        LepPt pt1;
        ALpPt pt2;
        LepAEta aeta1;
        ALpAEta aeta2;
        bool keepEvent() {
            if (pt1   () < 15. ) return false;
            if (pt2   () < 15. ) return false;
            if (aeta1 () > 2.5 ) return false;
            if (aeta2 () > 2.5 ) return false;
            return true;
        }
    }
}


#include "user_cuts.h"



#endif /* ifndef KinematicCuts_H */
