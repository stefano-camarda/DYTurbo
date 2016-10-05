#ifndef TurboHist_P1_H
#define TurboHist_P1_H
/**
 * @file TurboHist_P1.h
 * @brief A brief description
 *
 * Detailed description of file
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-10-05
 */

#include "TurboHist_HBase.h"
#include "TurboHist_Counter.h"

namespace TurboHist {
    struct P1 : public HBase<P1,Averager> {
        P1() { SetBins(1,0.,1.);};

        P1(const String &_name, const String &_title,
                const size_t &NX, const OBS &minX, const OBS &maxX
                ){
            name = _name;
            SetBins( NX,minX,maxX, _title);
        };

        inline size_t FindBin(const OBS &valX) const {
            return binsX.FindBin(valX);
        };

        void Fill (const OBS &valX, const OBS &valY, const PRE &weight=1.0){
            data[FindBin(valX)].increment(valY,weight);
            entries++;
        };


        void SetBins(size_t NX, OBS minX, OBS maxX,
                const String &tit="X;Y"){
            VecObs newbins=Binning::GetEquidistantVector(NX,minX,maxX);
            SetBinsAxis(X,newbins,tit);
        };
    };
}

#endif /* ifndef TurboHist_P1_H */
