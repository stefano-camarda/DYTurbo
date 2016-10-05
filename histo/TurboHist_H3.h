#ifndef TurboHist_H3_H
#define TurboHist_H3_H
/**
 * @file TurboHist_H3.h
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
    struct H3 : public HBase<H3,Counter> {
        H3() { SetBins(1,0.,1., 1,0.,1., 1,0.,1.);};

        H3(const String &_name, const String &_title,
                const size_t &NX, const OBS &minX, const OBS &maxX,
                const size_t &NY, const OBS &minY, const OBS &maxY,
                const size_t &NZ, const OBS &minZ, const OBS &maxZ
                ){
            name = _name;
            SetBins( NX,minX,maxX, NY,minY,maxY, NZ,minY,maxY, _title);
        };

        inline size_t FindBin(const OBS &valX, const OBS &valY, const OBS &valZ) const {
            return GetBin(binsX.FindBin(valX), binsY.FindBin(valY), binsZ.FindBin(valZ));
        };

        void Fill (const OBS &valX, const OBS &valY, const OBS &valZ, const PRE &weight=1.0){
            data[FindBin(valX,valY,valZ)]+=weight;
            entries++;
        };


        void SetBins(size_t NX, OBS minX, OBS maxX,
                size_t NY, OBS minY, OBS maxY,
                size_t NZ, OBS minZ, OBS maxZ,
                const String &tit="X;Y;Z"){
            VecObs newbins=Binning::GetEquidistantVector(NX,minX,maxX);
            SetBinsAxis(X,newbins,tit);
            newbins=Binning::GetEquidistantVector(NY,minY,maxY);
            SetBinsAxis(Y,newbins,tit);
            newbins=Binning::GetEquidistantVector(NZ,minZ,maxZ);
            SetBinsAxis(Z,newbins,tit);
        };
    };
}
#endif /* ifndef TurboHist_H3_H */
