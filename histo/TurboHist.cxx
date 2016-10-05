#ifndef TurboHist_CXX
#define TurboHist_CXX
/**
 * @file TurboHist.cxx
 * Description of cpp file
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-09-05
 */

#include "TurboHist.h"
#include "TurboHist_File.h"

#include "TurboHist_HBase.h"
template<class H, class C> int  TurboHist::HBase<H,C>::dim = -1;
template<class H, class C> char TurboHist::HBase<H,C>::type = 'o';

#include "TurboHist_H1.h"
template<> int  TurboHist::HBase<TurboHist::H1,TurboHist::Counter>::dim = 1;
template<> char TurboHist::HBase<TurboHist::H1,TurboHist::Counter>::type = 'h';

#include "TurboHist_H2.h"
template<> int  TurboHist::HBase<TurboHist::H2,TurboHist::Counter>::dim = 2;
template<> char TurboHist::HBase<TurboHist::H2,TurboHist::Counter>::type = 'h';

// #include "TurboHist_H3.h"
// template<> int  TurboHist::HBase<TurboHist::H3,TurboHist::Counter>::dim = 3;
// template<> char TurboHist::HBase<TurboHist::H3,TurboHist::Counter>::type = 'h';
// 
// #include "TurboHist_P1.h"
// template<> int  TurboHist::HBase<TurboHist::P1,TurboHist::Counter>::dim = 1;
// template<> char TurboHist::HBase<TurboHist::P1,TurboHist::Counter>::type = 'p';
// 
// #include "TurboHist_P2.h"
// template<> int  TurboHist::HBase<TurboHist::P2,TurboHist::Counter>::dim = 2;
// template<> char TurboHist::HBase<TurboHist::P2,TurboHist::Counter>::type = 'p';

#include <algorithm>
using std::transform;
using std::ios;

namespace TurboHist {


    void File::Open(string name, string method){
        fname=name;
        transform(method.begin(), method.end(), method.begin(), ::tolower);
        if( method.compare("recreate")==0 ){
            //! @todo ios::binary
            fst.open(name.c_str(), ios::out|ios::trunc);
        } else { // READ
            fst.open(name.c_str(), ios::in);
        }
    };


    Counter operator* ( const Counter &a , const double &b ) {
        Counter out (a);
        out*=b;
        return out;
    }

    Counter operator* ( const double &b, const Counter &a  ) {
        return a*b;
    }

    void MergeFiles(string outname,VecStr inNames){
        H1 hout;
        File f;
        bool isFirst = true;
        for (auto inname: inNames){
            H1 htmp;
            f.Open(inname);
            if (isFirst){
                f >> hout;
                isFirst=false;
            } else {
                f >> htmp;
                hout.Add(htmp);
            }
            f.Close();
        }
        f.Open(outname,"RECREATE");
        f << hout;
        f.Close();
    };
}




#endif /* ifndef TurboHist_CXX */
