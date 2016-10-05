#ifndef TurboHist_H
#define TurboHist_H
/**
 * @file TurboHist.h
 * Simple histograming with STL
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-09-01
 */

#include "src/handy_typdefs.h"

namespace TurboHist {
    // OBS: observable ususaly double but I think float is fine
    typedef double OBS;
    typedef std::vector<OBS> VecObs;
    typedef VecObs::iterator VecObsItr;
    typedef VecObs::const_iterator VecObsCItr;

    // PRE: precission of data -- here we need doubles
    typedef double PRE;
    typedef std::vector<PRE> VecPre;

    enum AxisName {
        X=0, Y=1, Z=2
    };

    struct OBase {
        virtual OBase * Clone(String name) const = 0;
    };

    class File;
    class Counter;
    class Averager;
    class Binning;
    template<class HistoType, class CountType> class HBase;

};

#endif /* ifndef TurboHist_H */
