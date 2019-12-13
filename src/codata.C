#include <math.h>
#include "codata.h"
//Provides conversion from GeV^-2 to fb by evaluating (h/2pi * c / GeV)^2
//1 GeV = e[j] * 1V[j/C] *10^9 = 1.6021766208e-19 * 10^9 j
//1 fb = 10^-15 * 10^-28 m^2 = 10^-43 m^2
//evaluate (h/2pi)^2
//values from CODATA 2014 (arXiv:1507.07956)
const double cc = 299792458; // speed of light [m*s^-1]
const double e = 1.6021766208e-19; //elementary charge [C]
const double h = 6.626070040e-34; //planck constant [j*s^-1]
const double hred = h/(2.*M_PI);  //reduced planck constant [j*s^-1]
const double gevfb = pow(cc*hred/e,2) * 1e-18 * 1e43; //[GeV^-2 * fb]

//const double gevfb = 3.893793656596451e+11;
//const double gevfb = 3.8937966e11;
