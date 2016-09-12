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


#include <vector>
using std::vector;

#include <string>
using std::string;

namespace TurboHist {
    // OBS: observable ususaly double but I think float is fine
    typedef double OBS;
    // PRE: precission of data -- here we need doubles
    typedef double PRE;

    typedef vector<string> VecStr;
    typedef vector<OBS> VecObs;
    typedef VecObs::iterator VecObsItr;
    typedef VecObs::const_iterator VecObsCItr;
    typedef vector<PRE> VecPre;

    enum AxisName {
        X=0, Y=1, Z=2
    };

    struct OBase {
        // Empty, but mother class for storing in containers
    };
    class File;
    class Counter;
    class Binning;
    template<class HistoType, class CountType> class HBase;

};

#endif /* ifndef TurboHist_H */
