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

using std::transform;
using std::ios;

namespace TurboHist {
    template<class T> int HBase<T>::dim = -1;
    template<class T> char HBase<T>::type = 'o';

    template<> int HBase<H1>::dim = 1;
    template<> char HBase<H1>::type = 'h';


    void File::Open(string name, string method){
        fname=name;
        transform(method.begin(), method.end(), method.begin(), ::tolower);
        if( method.compare("recreate")==0 ){
            // todo ios::binary
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
