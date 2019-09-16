#ifndef npff_h
#define npff_h
#include <complex>

using namespace std;

namespace npff
{
  complex <double> S(complex <double> b, double m, double x1, double x2);

  //Flavour dependent form factors
  extern complex <double> uvff;
  extern complex <double> usff;
  extern complex <double> dvff;
  extern complex <double> dsff;
  extern complex <double> ssff;
  extern complex <double> chff;
  extern complex <double> boff;
  extern complex <double> glff;
#pragma omp threadprivate(uvff,usff,dvff,dsff,ssff,chff,boff,glff)

}
#endif
