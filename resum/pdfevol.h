#ifndef pdfevol_h
#define pdfevol_h

#include <complex>
using namespace std;

namespace pdfevol
{
  extern complex <double> *fn1;
  extern complex <double> *fn2;
  
  extern void init();
  extern void evolve(int i1, int i2, int sign);
}

#endif
