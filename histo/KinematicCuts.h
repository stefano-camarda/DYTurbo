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
    /**
     * @brief Functionality around rejecting events
     */
    namespace Cuts {

        /// Constant bool to keep event.
        extern const bool KeepEvent;
        /// Constant bool to skip event.
        extern const bool SkipEvent;

        /**
         * @brief Base class for cuts.
         *
         * Contains only call operator, which makes cut decision. It is
         * implemented by standard cuts or user cuts.
         *
         * Class was chosen instead of simple function for better control of
         * variables used in cut decision (i.e. user and standard can have
         * same or different names to variables and also add functions).
         *
         * @todo Test user defined cut with user defined observable.
         */
        struct CutBase{
            virtual bool operator()()=0;
        };

        ///  function
        bool KeepThisEvent(double p3[4], double p4[4]);
    }
}





#endif /* ifndef KinematicCuts_H */
