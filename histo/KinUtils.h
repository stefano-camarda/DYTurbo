#ifndef KinUtils_H
#define KinUtils_H
/**
 * @file KinUtils.h
 * Description of cpp file
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-09-22
 */

#include <cmath>

namespace Kinematics {
    namespace Util {

        inline double pT(double x,double y){
            return sqrt(x*x + y*y);
        }

        inline double pT2(double x,double y){
            return (x*x + y*y);
        }

        inline double mom( double x,  double y,  double z){
            return sqrt(x*x + y*y + z*z);
        }

        inline double rap( double e, double z){
            return 0.5 *log((e+z) / (e-z));
        }

        inline double eta( double x,  double y,  double z){
            double p = mom(x,y,z);
            return rap(p,z);
        }

        inline double mass2( double x,  double y,  double z, double e){
            return (e*e - (x*x + y*y + z*z));
        }
    }
}


#endif /* ifndef KinUtils_H */
