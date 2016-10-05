#ifndef TurboHist_H2_H
#define TurboHist_H2_H
/**
 * @file TurboHist_H2.h
 * Description of cpp file
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-09-07
 */

#include "TurboHist_HBase.h"
#include "TurboHist_Counter.h"

namespace TurboHist {
    struct H2 : public HBase<H2,Counter> {
        H2() { SetBins(1,0.,1., 1,0.,1.);};

        H2(const String &_name, const String &_title,
                const size_t &NX, const OBS &minX, const OBS &maxX,
                const size_t &NY, const OBS &minY, const OBS &maxY
                ){
            name = _name;
            SetBins( NX,minX,maxX, NY,minY,maxY, _title);
        };


        inline size_t FindBin(const OBS &valX, const OBS &valY) const {
            return GetBin(binsX.FindBin(valX), binsY.FindBin(valY));
        };

        void Fill (const OBS &valX, const OBS &valY, const PRE &weight=1.0){
            data[FindBin(valX,valY)]+=weight;
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

#endif /* ifndef TurboHist_H2_H */
