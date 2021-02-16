#ifndef TurboHist_P2_H
#define TurboHist_P2_H
/**
 * @file TurboHist_P2.h
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
    struct P2 : public HBase<P2,Averager> {
        P2() { SetBins(1,0.,1., 1,0.,1.);};

        P2(const String &_name, const String &_title,
                const size_t &NX, const OBS &minX, const OBS &maxX,
                const size_t &NY, const OBS &minY, const OBS &maxY
                ){
            name = _name;
            SetBins( NX,minX,maxX, NY,minY,maxY, _title);
        };


        inline size_t FindBin(const OBS &valX, const OBS &valY) const {
            return GetBin(binsX.FindBin(valX), binsY.FindBin(valY));
        };

        void Fill (const OBS &valX, const OBS &valY, const OBS &valZ, const PRE &weight=1.0){
            data[FindBin(valX,valY)].increment(valZ,weight);
            entries++;
        };


        void SetBins(size_t NX, OBS minX, OBS maxX,
                size_t NY, OBS minY, OBS maxY,
                const String &tit="X;Y"){
            VecObs newbins=Binning::GetEquidistantVector(NX,minX,maxX);
            SetBinsAxis(X,newbins,tit);
            newbins=Binning::GetEquidistantVector(NY,minY,maxY);
            SetBinsAxis(Y,newbins,tit);
        };
    };
}

#endif /* ifndef TurboHist_P2_H */
