#include "psi.h"
#include "fcomplex.h"

extern "C"
{
  //psi function and derivatives from Pegasus
  fcomplex psi_ (fcomplex &z);
  fcomplex dpsi_ (fcomplex &z, int &m);
}

complex <double> cpsi0(complex <double> z)
{
  fcomplex zz = fcx(z);
  return cx(psi_(zz));
}

complex <double> cpsi(int n, complex <double> z)
{
  fcomplex zz = fcx(z);
  if (n == 0)
    return cx(psi_(zz));
  else
    return cx(dpsi_(zz, n));
}
