#ifndef Kinematics_H
#define Kinematics_H

/**
 * @file Kinematics.h
 * 
 *
 * @brief Declarations for Kinematics namespace.
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-08-26
 */

#include <set>
#include <cmath>



/** @brief Kinematic calculations and cuts.
 *
 *  This namespace is used to calculate kinematic variables. And provide values
 *  for histograms and cuts.
 *
 *  For input to integration is used @ref phasespace.
 */
namespace Kinematics {

    /** 
     * @brief Call this to set new kinematic values.
     *
     * It will update lepton four-momenta and send signal to all variables
     * to recalculate their values.
     */
    void SetKinematics(double l1[4],double l2[4],double wgt);

    /**
     * @brief Update Observables in Integraion mode.
     *
     * Only selected Observable does have value in Integration mode. For those
     * it is necessary to update value.  Kinematic point will be taken from
     * \ref phasespace as middle of bin. This is implemented per each Observable in function \ref middlePoint .
     */
    void SetMiddlePoint();

    //! Input four-momenta of first lepton.
    extern double p3[4];
    //! Input four-momenta of second lepton.
    extern double p4[4];
    //! The weight of event. This will be used to fill histograms.
    extern double event_weight;

    //! Flag whether using Integration mode.
    extern bool isIntegratorMode;

    /** 
     * @brief Holder of flags (isCalculated)
     *
     * Is implemented as `std::set` of `bool` pointers. Every \ref Observable
     * adds its pointer in constructor. The `std::set` assures that every
     * pointer is there only once.
     *
     * @note This is assuming that pointer to static class member is constant.
     */
    typedef std::set<bool *> ObservableFlags;
    typedef ObservableFlags::iterator ObservableFlagsItr;
    extern ObservableFlags flags; ///< Storing pointers to all Observable flags.

    
    /** @brief Main class of Kinematics. All observables are inherited from this.
     *
     * ### Main idea 
     *
     * > Calculate only when kinematics are updated.
     *
     * To provide this functionality there is boolean flag, which stores
     * whether observable was already calculated. When kinematics are updated
     * all flags are set to false.
     *
     * ### Structure of class: Curiously recurring template (CRTP)
     *
     * This will allows us to define static variable per each class.  Make sure
     * you not recalculate variable only once.  To read more on
     * [wiki](https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern).
     *
     *
     * Well commented example of implementation can be found in \ref CosTh .
     */
    template <class T> class Observable {
        public :
            /// @brief Default constructor.
            Observable() {
                /// Adds pointer to flag into flag container. This is one flag per class.
                flags.insert(&isCalculated);
            };

            //! @brief Call operator returns actual value of observable.
            double operator()(){
                /// If available in Integration mode return middle point.
                if(isIntegratorMode) return IsIntegrableObservable() ? middlePoint() : 666.;
                /// If not calculated yet it runs \ref calc.
                if (!isCalculated){
                    value=calc();
                    isCalculated=true;
                }
                return value;
            }

            //! Wrapper around flag isCalculated.
            inline bool IsCalculated() const {return isCalculated;}

        public :
            /** @brief Is observable available in integration mode.
             *
             * Default return value is false. If observable is available in
             * integration mode you can reimplement this function and
             * reimplement calculation in \ref middlePoint .
             */
            virtual inline bool IsIntegrableObservable() const {return false;}

        protected :
            /// Body of calculation. Must be implemented per each child.
            virtual double calc()=0;

            /** @brief Return middle point value.
             *
             * If observable is available in integration mode obtain value of
             * bin from \ref phasespace and return middle point.
             *
             * @attention Dont forget to reimplement IsIntegrableObservable to return true.
             */
            virtual double middlePoint(){return 0.;};

            //! Value of observable. Is static => every instance will be updated with first call.
            static double value;

            //! Flag saying if its needed to recalculate. Is static => every instance will be updated with first call.
            static bool isCalculated;
    };

    template<class T> double Observable<T>::value = -999.;
    template<class T> bool Observable<T>::isCalculated = false;
}


#endif /* ifndef Kinematics_H */
