#ifndef TurboHist_H1_H
#define TurboHist_H1_H
/**
 * @file TurboHist_H1.h
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
    struct H1 : public  HBase<H1,Counter> {
        H1() {SetBins(1,0.,1.);};

        H1(const String &_name, const String &_title,const size_t &N, const OBS &min, const OBS &max){
            name = _name;
            SetBins(N,min,max,_title);
        };

        inline size_t FindBin(const OBS &val) const {
            return binsX.FindBin(val);
        };

        void Fill (const OBS &val,const PRE &weight=1.0){
            data[FindBin(val)]+=weight;
            entries++;
        };

        void SetBins(size_t N, OBS min, OBS max, const String &tit="X"){
            VecObs newbins=Binning::GetEquidistantVector(N,min,max);
            SetBinsAxis(X,newbins,tit);
        };

        void SetBins(const VecObs &newbins, const String &tit="X"){
            SetBinsAxis(X,newbins,tit);
        };

    };
}


#endif /* ifndef TurboHist_H1_H */
