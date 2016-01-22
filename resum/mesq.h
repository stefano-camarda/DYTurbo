#ifndef mesq_h
#define mesq_h

#include "mellinint.h"

extern const double eequ;
extern const double eeqd;
extern const double gevfb;

//Born level amplitudes for the resummed integrand multiplied by the x1^-z1 * x2^-z2 piece of the Mellin inverse transform
//(1/(2pi*i))^2 * x1^-z1 * x2^-z2 = (CCp/pi)^2 * (M/s)^(-(z1+z2)) * exp(-y(z1-z2))
namespace mesq
{
  //constants
  extern double fac;
  extern double mZ2;   //Z couplings
  extern double wZ2;
  extern double gLZu;
  extern double gLZd;
  extern double gRZu;
  extern double gRZd;
  extern double fLZ;
  extern double fRZ;
  extern double fLpfR;
  extern double fLmfR;
  extern double ugLpgR;
  extern double ugLmgR;
  extern double dgLpgR;
  extern double dgLmgR;
  extern double mW2;   //W coupling
  extern double wW2;
  extern double gLWfLW;
  extern double aem2pi;   //gamma* coupling
  extern double aem2pi2;

  //Initialisation of constants
  extern void init();

  //mass dependent variables
  extern double q2;
  extern double propZ;
  extern double propW;
  extern double propG;
  extern double propZG;
  
  //Amplitudes times mellin inverse transform piece x1^-z1 * x2^-z2
  //Depending on the mode 0, 1, 2 they can be differential, cos theta integrated, or cos theta and y integrated
  extern complex <double> *mesqij_expy;
  extern void setmesq_expy(int mode, double m, double costh, double y);

  //amplitudes
  template <class T>
  extern void setmesq(T one, T costh1, T costh2);
  extern complex <double> mesqij[12];

  enum sign {positive=0, negative=1};

  //Return index in the array based on partonic channel, i1 and i2 indexes in the z Mellin space, and sign of the branch
  inline int index(int pch, int i1, int i2, bool sign)
  {return i2 + mellinint::mdim*(i1 + mellinint::mdim*(sign + 2*pch));}

  //Number of partonic channels
  //extern int totpch;
  
  //function that returns the partonic channel index given f1 and f2 (not implemented yet)
  inline int pchindx(int f1, int f2)
  {
    //PDG number scheme:
    //d u s c b t
    //1 2 3 4 5 6
    //Z
    //   0  ->    2 , -2
    //   1  ->    -2, 2 
    //   2  ->    1 , -1
    //   3  ->    -1, 1 
    //   4  ->    3 , -3
    //   5  ->    -3, 3 
    //   6  ->    4 , -4
    //   7  ->    -4, 4 
    //   8  ->    5 , -5
    //   9  ->    -5, 5
    //W+
    //   0  ->     2	-1 
    //   1  ->     -1	2  
    //   2  ->     2	-3 
    //   3  ->     -3	2  
    //   4  ->     2	-5 
    //   5  ->     -5	2  
    //   6  ->     4	-3 
    //   7  ->     -3	4  
    //   8  ->     4	-1 
    //   9  ->     -1	4 
    //   10 ->     4	-5
    //   11 ->     -5	4 
    //W-
    //   0  ->    1	-2	
    //   1  ->    -2	1	
    //   2  ->    3	-2	
    //   3  ->    -2	3	
    //   4  ->    5	-2	
    //   5  ->    -2	5	
    //   6  ->    3	-4	
    //   7  ->    -4	3	
    //   8  ->    1	-4	
    //   9  ->    -4	1	
    //   10 ->    5	-4	
    //   11 ->    -4	5	
    return 0;
  }

  inline int index(int f1, int f2, int i1, int i2, bool sign)
  {return i2 + mellinint::mdim*(i1 + mellinint::mdim*(sign + 2*(pchindx(f1,f2))));}
}

#endif
