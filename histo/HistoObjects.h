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

#include "HistoHandler.h"
#include "config.h"

using std::string;

namespace HistoHandler {
    // Base Histogram class: make possible to store in one container
    class HistoBase{
        public :
            virtual void FillEvent()=0;
            virtual void FillDipole()=0;
            virtual void FillRealEvent()=0;
            virtual void SetVariation(const KeySuffix)=0;
            virtual void Save()=0;
            virtual void Delete()=0;
            virtual void Clear()=0;
            virtual void AddToBin(double int_val,double int_err)=0;
            virtual double GetEntries() const =0; //{ return 666.;};
            virtual const char* GetName() const =0; // {return "Belzeebos";};
        protected :
            string name;
            string title;
    };
}

#ifdef USEROOT
#include "HistoObjectsROOT.h"
#else
#include "HistoObjectsSTL.h"
#endif

#include "HistoBook.h"





#endif /* ifndef HistoObjects_H */
