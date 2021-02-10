#ifndef psi_h
#define psi_h

#include <complex>
using namespace std;

//Interface to psi function and derivatives from Pegasus
extern complex <double> cpsi0(complex <double> z);
extern complex <double> cpsi(int n, complex <double> z);

//psi functions are also available in melfun (they should be faster)

#endif
