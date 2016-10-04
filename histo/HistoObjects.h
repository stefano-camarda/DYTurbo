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

#include "HistoHandler.h"


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
            virtual void Reset()=0;
            virtual void AddToBin(double int_val,double int_err)=0;
            virtual double GetEntries() const =0; //{ return 666.;};
            virtual const char* GetName() const =0; // {return "Belzeebos";};
            virtual bool IsIntegrationSafe() const =0; //{ return 666.;};
        protected :
            String name;
            String title;
    };


    // Object class: covering common functionality
    template <class TH>
    class HistoObject : public HistoBase {
    };
}

// need to be include config to decide whether to use ROOT
#include "config.h"

#ifdef USEROOT
#include "HistoObjectsROOT.h"
#else
#include "HistoObjectsSTL.h"
#endif

#include "HistoBook.h"





#endif /* ifndef HistoObjects_H */
