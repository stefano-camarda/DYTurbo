#ifndef pdfevol_h
#define pdfevol_h

#include "interface.h"

//The pdfevol::evolve(i1, i2, sign) function
//take the PDF moments corresponding to i1, i2, and sign
//from the fortran common block, and put it in fn1 fn2.

#include <complex>
using namespace std;

extern "C"
{
  void pdfmoments_(int &beam, fcomplex &N, fcomplex &UV, fcomplex &DV, fcomplex &US, fcomplex &DS, fcomplex &SS, fcomplex &GL, fcomplex &CH, fcomplex &BO);
}

namespace pdfevol
{
  //PDFs mellin moments at the factorisation scale
  //Beam 1
  extern complex <double> *UVP1;
  extern complex <double> *DVP1;
  extern complex <double> *USP1;
  extern complex <double> *DSP1;
  extern complex <double> *SSP1;
  extern complex <double> *GLP1;
  extern complex <double> *CHP1;
  extern complex <double> *BOP1;
  extern complex <double> *UVM1;
  extern complex <double> *DVM1;
  extern complex <double> *USM1;
  extern complex <double> *DSM1;
  extern complex <double> *SSM1;
  extern complex <double> *GLM1;
  extern complex <double> *CHM1;
  extern complex <double> *BOM1;
  //Beam 2
  extern complex <double> *UVP2;
  extern complex <double> *DVP2;
  extern complex <double> *USP2;
  extern complex <double> *DSP2;
  extern complex <double> *SSP2;
  extern complex <double> *GLP2;
  extern complex <double> *CHP2;
  extern complex <double> *BOP2;
  extern complex <double> *UVM2;
  extern complex <double> *DVM2;
  extern complex <double> *USM2;
  extern complex <double> *DSM2;
  extern complex <double> *SSM2;
  extern complex <double> *GLM2;
  extern complex <double> *CHM2;
  extern complex <double> *BOM2;

  extern complex <double> *fn1;
  extern complex <double> *fn2;

  extern complex <double> XL;
  extern complex <double> XL1;
  extern complex <double> SALP;

  extern complex <double> alpr;
  
  extern void init();
  extern void evolve(int i1, int i2, int sign);
  extern void evolution(int i, int sign, int beam);
}

#endif
