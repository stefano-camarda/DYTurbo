#include "constants.h"
#include <math.h>

//QCD Constants
const double constants::NC = 3.;    //Number of colours
const double constants::TF = 1./2.; //Trace
const double constants::CF = 4./3.; //Casimir
const double constants::CA = 3.;    //Casimir

//Mathematical constants

//Gamma_Euler
const double constants::euler = 0.57721566490153286061; //https://oeis.org/A001620
                                
//Values of the Riemann zeta function
const double constants::zeta2 = 1.64493406684822643647; //M_PI*M_PI/6.;
const double constants::zeta3 = 1.20205690315959428540; //https://oeis.org/A002117
const double constants::zeta4 = 1.08232323371113819152; //https://oeis.org/A013662
const double constants::zeta5 = 1.03692775514336992633; //https://oeis.org/A013663
const double constants::zeta6 = 1.01734306198444913971; //https://oeis.org/A013664

//Powers of Pi
const double constants::pi  = M_PI;
const double constants::pi2 = pow(M_PI,2);
const double constants::pi3 = pow(M_PI,3);
const double constants::pi4 = pow(M_PI,4);
const double constants::pi5 = pow(M_PI,5);
const double constants::pi6 = pow(M_PI,6);
