#include "hs.h"

complex <double> HS(int k1,int k2,int k3,int k4,int k5,complex<double> N,int eta)
{
  fcomplex n = fcx(N);
  double h = double(eta);
  return cx(hs_(k1,k2,k3,k4,k5,n,h));
}

