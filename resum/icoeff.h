#ifndef icoeff_h
#define icoeff_h

#include <complex>
#include "fcomplex.h"

using namespace std;

extern "C"
{
  void psi0_(fcomplex &zz, fcomplex &res);
}

namespace icoeff
{
  const int rule = 500;

  extern double tlog[rule],tlin[rule];                //gauss nodes
  extern double faclog[rule],faclin[rule];            //overall factor

  //kernel of the Mellin transform
  extern complex <double> *kernlog,*kernlin;
  extern complex <double> *kernlog_1,*kernlin_1;
  extern complex <double> *kernlog_2,*kernlin_2;
  //kernlog and kernlin has to be threadprivate!!!
  
  //C coefficients in z-space (stored on disk)
  extern double c3qqpzlog[rule];
  extern double c3qqbzlog[rule];
  extern double c3qqbpzlog[rule];
  extern double c3qgzlin[rule];
  extern double c3qqzlin[rule];
  
  //C coefficients in Mellin space
  extern double C1qq_delta;
  extern double C2qq_delta;
  extern double C3qq_delta;

  extern double C2qq_plus;
  extern double C3qq_plus;
  
  extern complex <double> *c1qg;
  extern complex <double> *c1qq;
  extern complex <double> *c1qqb;
  extern complex <double> *c1qqp;
  extern complex <double> *c1qqbp;

  extern complex <double> *c2qg;
  extern complex <double> *c2qq;
  extern complex <double> *c2qqb;
  extern complex <double> *c2qqp;
  extern complex <double> *c2qqbp;

  extern complex <double> *c3qg;
  extern complex <double> *c3qq;
  extern complex <double> *c3qqb;
  extern complex <double> *c3qqp;
  extern complex <double> *c3qqbp;

  extern complex <double> *c3qg_1;
  extern complex <double> *c3qq_1;
  extern complex <double> *c3qqb_1;
  extern complex <double> *c3qqp_1;
  extern complex <double> *c3qqbp_1;
  extern complex <double> *c3qg_2;
  extern complex <double> *c3qq_2;
  extern complex <double> *c3qqb_2;
  extern complex <double> *c3qqp_2;
  extern complex <double> *c3qqbp_2;
  
  extern void allocate();
  extern void init();
  extern void delta();
  extern void calc1d();
  extern void calc2d();
  extern void release();
}

#endif
