#ifndef fcomplex_h
#define fcomplex_h

#include <complex>

using namespace std;

extern "C" {

  struct fcomplex {double real; double imag;};
  //function to convert from fortran to C complex type
  inline complex <double> cx(fcomplex fcx) {return complex<double>(fcx.real, fcx.imag);}
  //function to convert from C to fortran complex type
  //  inline fcomplex fcx(complex <double> cx) {fcomplex f; f.real = real(cx); f.imag = imag(cx); return f;};
  inline fcomplex fcx(complex <double> cx) {return (fcomplex) {real(cx), imag(cx)};};

}
#endif
