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
#include <cmath>

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

        inline void value(PRE val){ sum_w=val; };
        inline void error(PRE err){ sum_w2=err*err; };

        inline PRE value() const { return sum_w;}
        inline PRE error() const { return sqrt(sum_w2);}

        inline void reset(){
            sum_w=0.;
            sum_w2=0.;
        }

    };
    Counter operator* ( const Counter &a , const double &b );
    Counter operator* ( const double &b, const Counter &a  );

    struct Averager{
        PRE sum_w;
        PRE sum_w2;
        PRE sum_aw;
        PRE sum_a2w;
        Averager() : sum_w(0.), sum_w2(0.), sum_aw(0.), sum_a2w(0.) { };

        Averager(const Averager& rhs) :
            sum_w(rhs.sum_w),
            sum_w2(rhs.sum_w2),
            sum_aw(rhs.sum_aw),
            sum_a2w(rhs.sum_a2w)
        { };

        inline Averager operator = (const Averager &rhs){
            sum_w=rhs.sum_w;
            sum_w2=rhs.sum_w2;
            sum_aw=rhs.sum_aw;
            sum_a2w=rhs.sum_a2w;
            return *this;
        };

        inline Averager & increment (OBS a, PRE w){
            sum_w+=w;
            sum_w2+=w*w;
            sum_aw+=a*w;
            sum_a2w+=a*a*w;
            return *this;
        };

        inline Averager & operator += (const Averager &rhs){
            sum_w+=rhs.sum_w;
            sum_w2+=rhs.sum_w2;
            sum_aw+=rhs.sum_aw;
            sum_a2w+=rhs.sum_a2w;
            return *this;
        };

        inline Averager & operator *= (const double &weight){
            sum_w*=weight;
            sum_w2*=weight;
            sum_aw*=weight;
            sum_a2w*=weight;
            return *this;
        };

        inline void SaveToFile(File &file) const {
            file << sum_w;
            file << sum_w2;
            file << sum_aw;
            file << sum_a2w;
        }

        inline void LoadFromFile(File &file){
            file >> sum_w;
            file >> sum_w2;
            file >> sum_aw;
            file >> sum_a2w;
        }

        void set_suma2w_keep_sumw_sumw2 (PRE val,PRE err){
            // change sum_a2w while sumw and sumw2 keep same and new val and err
            //err*err = fabs(oNeff*( sum_a2w / sum_w - val*val) )
            //err*err*Neff=fabs(sum_a2w / sumw - val*val)
            //err*err*Neff*fabs(sum_w) = fabs(sum_a2w - val*val*sum_w)
            double C = err*err*fabs(sum_w)*sum_w*sum_w/sum_w2;
            double B = val*val*sum_w;
            sum_a2w = C<2*B ?  B-C : C-B;
        };

        inline void value(PRE val){
            double err = error();
            //if(err==0) err=sqrt(val); // else keep same error
            if(sum_w==0) sum_w=1; // else keep same sum_w
            if(sum_w2==0) sum_w2=1; // else keep same sum_w2
            sum_aw = val*sum_w;
            set_suma2w_keep_sumw_sumw2(val,err);
        };

        inline void error(PRE err){
            err = fabs(err);
            double val = value(); // keep same value
            if(sum_w==0) sum_w=1; // else keep same sum_w
            if(sum_w2==0) sum_w2=1; // else keep same sum_w2
            set_suma2w_keep_sumw_sumw2(val,err);
        };

        inline PRE value() const { return sum_w!=0 ? sum_aw/sum_w : 0;}
        inline PRE error() const {
            if (sum_w==0) return 0;
            double mu = value();
            double oNeff = sum_w2 / (sum_w*sum_w);
            double Var = sum_a2w/sum_w - mu*mu;
            return sqrt(fabs(oNeff*Var));
        }

        inline void reset(){
            sum_w=0.;
            sum_w2=0.;
            sum_aw=0.;
            sum_a2w=0.;
        }

    };
    Averager operator* ( const Averager &a , const double &b );
    Averager operator* ( const double &b, const Averager &a  );
}

#endif /* ifndef TurboHist_Counter_H */
