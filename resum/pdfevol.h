#ifndef pdfevol_h
#define pdfevol_h

#include "interface.h"

#include <complex>
using namespace std;

extern "C"
{
  void pdfmoments_(int &beam, double &scale, fcomplex &N, fcomplex &UV, fcomplex &DV, fcomplex &US, fcomplex &DS, fcomplex &SP, fcomplex &SM, fcomplex &GL, fcomplex &CH, fcomplex &BO);
}

namespace pdfevol
{
  //PDFs mellin moments at the factorisation scale
  extern complex <double> *UVP;
  extern complex <double> *DVP;
  extern complex <double> *USP;
  extern complex <double> *DSP;
  extern complex <double> *SSP;
  extern complex <double> *GLP;
  extern complex <double> *CHP;
  extern complex <double> *BOP;

  //moments for the PDFs convolution
  extern complex <double> *fn1;
  extern complex <double> *fn2;

  //scales
  extern complex <double> bscale;
  extern complex <double> bstarscale;
  extern complex <double> XL;
  extern complex <double> XL1;
  extern complex <double> SALP;
  extern complex <double> alpr;
  
  //initialise and compute Mellin moments of PDFs at the starting scale (factorisation scale)
  extern void init();
  //evolve Mellin moments of PDFs from the factorisation scale to the scale ~1/b, set in bscale
  extern void evolution(int i);
  //calculate Mellin moments of PDFs at all scales by direct Mellin transform
  extern void calculate(int i);
  //store the moments in the Fortran common block
  extern void storemoments(int i, complex <double> fx[11]);
  //retrieve moments corresponding to i1, i2, and sign from the Fortran common block and store them in fn1 and fn2
  extern void retrieve(int i1, int i2, int sign);
}

#endif
