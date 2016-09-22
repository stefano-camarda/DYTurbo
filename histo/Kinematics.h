#ifndef Kinematics_H
#define Kinematics_H

/**
 * @file Kinematics.h
 * Description of this macro
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-08-26
 */

#include <set>
#include <cmath>


// TODO: Implementation to DYTURBO. Use Kinematics for cuts.

namespace Kinematics {
    void SetKinematics(double l1[4],double l2[4],double wgt);
    void SetMiddlePoint();

    // input values declaration
    extern double p3[4];
    extern double p4[4];
    extern double event_weight;
    extern bool isIntegratorMode;

    // Struct for flag (isCalculated)
    typedef std::set<bool *> VariableFlags;
    typedef VariableFlags::iterator VariableFlagsItr;
    extern VariableFlags flags;

    // Curiously reccurring template (CRTP)
    // This will allows us to define static variable per each class.
    // Make sure you not recalculate variable only once.
    // read more: https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
    template <class T> class Variable {
        public :
            virtual double calc()=0;
            virtual double middlePoint(){return 0.;};

            Variable() {
                // add pointer to flag, this is common flag per class
                flags.insert(&isCalculated);
            };

            double operator()(){
                // Integration mode
                if(isIntegratorMode) return IsIntegratorVariable() ? middlePoint() : 666.;
                // Filler Mode
                // calculate if needed
                if (!isCalculated){
                    value=calc();
                    isCalculated=true;
                }
                return value;
            }

            inline bool IsCalculated() const {return isCalculated;}
            virtual inline bool IsIntegratorVariable() const {return false;}

        protected :
            // Static will be per derived class
            static double value;
            static bool isCalculated;
    };

    template<class T> double Variable<T>::value = -999.;
    template<class T> bool Variable<T>::isCalculated = false;
}


#endif /* ifndef Kinematics_H */
