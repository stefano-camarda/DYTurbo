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
#include <typeinfo>
#include <iostream>

using std::set;
using std::type_info;


// TODO: Implementation to DYTURBO. Use Kinematics for cuts.

namespace Kinematics {
    void SetKinematics(double l1[4],double l2[4],double wgt);
    void SetMiddlePoint();

    // input values declaration
    extern double p3[4];
    extern double p4[4];
    extern double event_weight;

    // Struct for flag (isCalculated)
    typedef set<bool *> VariableFlags;
    typedef VariableFlags::iterator VariableFlagsItr;
    extern VariableFlags flags;

    // Curiously reccurring template (CRTP)
    // This will allows us to define static variable per each class.
    // Make sure you not recalculate variable only once.
    // read more: https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
    template <class T> class Variable {
        public :
            virtual double calc()=0;
            Variable() {
                flags.insert(&isCalculated);
                //printf ("%s %p : Creating variable with pointer %p\n", typeid(*this).name(), this, &(this->isCalculated));
            };

            double operator()(){
                if (!isCalculated){
                    value=calc();
                    isCalculated=true;
                    //printf("%s %p : Calc was eval to %f\n",typeid(*this).name(),this, this->value);
                } else {
                    //printf("%s %p : Calc is already calculated\n", typeid(*this).name(), this);
                }
                return value;
            }

            void SetValue(double val){ // For Integrator mode: set middle of the bin.
                value = val;
                isCalculated = true;
            }

            inline bool IsCalculated() const {return isCalculated;}

        protected :
            // Static will be per derived class
            static double value;
            static bool isCalculated;
    };

    template<class T> double Variable<T>::value = -999.;
    template<class T> bool Variable<T>::isCalculated = false;
}


//#include "KinematicDefinitions.h"

#endif /* ifndef Kinematics_H */
