#ifndef HistoObjects_H
#define HistoObjects_H

/**
 * @file HistoObjects.h
 * Description of this macro
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-08-26
 */

#include <string>
using std::string;

namespace HistoHandler {
    // Base Histogram class: make possible to store in one container
    class HistoBase{
        public :
            virtual void FillEvent(){};
            virtual void FillDipole(){};
            virtual void FillRealEvent(){};
            virtual void SetVariation(const KeySuffix){};
            virtual void Save(){};
            virtual void AddToBin(double int_val,double int_err){};
        protected :
            string name;
            string title;
    };
}







#endif /* ifndef HistoObjects_H */
