#ifndef TurboHist_Counter_H
#define TurboHist_Counter_H
/**
 * @file TurboHist_Counter.h
 * Description of cpp file
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-09-07
 */

#include "TurboHist.h"

namespace TurboHist {
    struct Counter{
        PRE sum_w;
        PRE sum_w2;
        Counter() : sum_w(0.), sum_w2(0.) { };

        Counter(const Counter& rhs) : sum_w(rhs.sum_w), sum_w2(rhs.sum_w2) { };

        inline Counter operator = (const Counter &rhs){
            sum_w=rhs.sum_w;
            sum_w2=rhs.sum_w2;
            return *this;
        };

        inline Counter operator += (const PRE &weight){
            sum_w+=weight;
            sum_w2+=weight*weight;
            return *this;
        };

        inline Counter & operator += (const Counter &rhs){
            sum_w+=rhs.sum_w;
            sum_w2+=rhs.sum_w2;
            return *this;
        };

        inline Counter & operator *= (const double &weight){
            sum_w*=weight;
            sum_w2*=weight*weight;
            return *this;
        };

        inline void SaveToFile(File &file) const {
            file << sum_w;
            file << sum_w2;
        }

        inline void LoadFromFile(File &file){
            file >> sum_w;
            file >> sum_w2;
        }

    };
    Counter operator* ( const Counter &a , const double &b );
    Counter operator* ( const double &b, const Counter &a  );
}

#endif /* ifndef TurboHist_Counter_H */
