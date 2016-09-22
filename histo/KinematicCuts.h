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
    namespace Cuts {

        extern const bool KeepEvent;
        extern const bool SkipEvent;

        // Cut base class
        struct CutBase{
            virtual bool operator()()=0;
        };

        bool KeepThisEvent(double p3[4], double p4[4]);
    }
}





#endif /* ifndef KinematicCuts_H */
