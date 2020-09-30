#ifndef hs_h
#define hs_h

#include "fcomplex.h"
#include <complex>

using namespace std;

extern "C"
{
  extern fcomplex hs_(int &k1,int &k2,int &k3,int &k4,int &k5,fcomplex &n,double &eta);
}

extern complex <double> HS(int k1,int k2,int k3,int k4,int k5,complex<double> N,int eta);

#endif
