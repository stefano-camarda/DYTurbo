#ifndef ccoeff_h
#define ccoeff_h

#include <complex>
#include "fcomplex.h"

using namespace std;

extern "C"
{
  void psi0_(fcomplex &zz, fcomplex &res);
}

namespace ccoeff
{
  extern double C1qq_delta;
  extern double C2qq_delta;
  extern double C3qq_delta;

  extern complex <double> *C1qg;
  extern complex <double> *C1qq;
  extern complex <double> *C1qqb;
  extern complex <double> *C1qqp;
  extern complex <double> *C1qqbp;

  extern complex <double> *C2qg;
  extern complex <double> *C2qq;
  extern complex <double> *C2qqb;
  extern complex <double> *C2qqp;
  extern complex <double> *C2qqbp;

  extern complex <double> *C3qg;
  extern complex <double> *C3qq;
  extern complex <double> *C3qqb;
  extern complex <double> *C3qqp;
  extern complex <double> *C3qqbp;
#pragma omp threadprivate(C1qg,C1qq,C1qqb,C1qqp,C1qqbp,C2qg,C2qq,C2qqb,C2qqp,C2qqbp,C3qg,C3qq,C3qqb,C3qqp,C3qqbp)

  extern complex <double> *C1qg_1;
  extern complex <double> *C1qq_1;
  extern complex <double> *C1qqb_1;
  extern complex <double> *C1qqp_1;
  extern complex <double> *C1qqbp_1;
  extern complex <double> *C2qg_1;
  extern complex <double> *C2qq_1;
  extern complex <double> *C2qqb_1;
  extern complex <double> *C2qqp_1;
  extern complex <double> *C2qqbp_1;
  extern complex <double> *C3qg_1;
  extern complex <double> *C3qq_1;
  extern complex <double> *C3qqb_1;
  extern complex <double> *C3qqp_1;
  extern complex <double> *C3qqbp_1;
#pragma omp threadprivate(C1qg_1,C1qq_1,C1qqb_1,C1qqp_1,C1qqbp_1,C2qg_1,C2qq_1,C2qqb_1,C2qqp_1,C2qqbp_1,C3qg_1,C3qq_1,C3qqb_1,C3qqp_1,C3qqbp_1)

  extern complex <double> *C1qg_2;
  extern complex <double> *C1qq_2;
  extern complex <double> *C1qqb_2;
  extern complex <double> *C1qqp_2;
  extern complex <double> *C1qqbp_2;
  extern complex <double> *C2qg_2;
  extern complex <double> *C2qq_2;
  extern complex <double> *C2qqb_2;
  extern complex <double> *C2qqp_2;
  extern complex <double> *C2qqbp_2;
  extern complex <double> *C3qg_2;
  extern complex <double> *C3qq_2;
  extern complex <double> *C3qqb_2;
  extern complex <double> *C3qqp_2;
  extern complex <double> *C3qqbp_2;
#pragma omp threadprivate(C1qg_2,C1qq_2,C1qqb_2,C1qqp_2,C1qqbp_2,C2qg_2,C2qq_2,C2qqb_2,C2qqp_2,C2qqbp_2,C3qg_2,C3qq_2,C3qqb_2,C3qqp_2,C3qqbp_2)
  
  extern void init();
  extern void allocate();
  extern void delta();
  extern void calc1d();
  extern void calc2d();
  extern void num_calc();
  extern void truncate();
  extern void free();
  extern void release();
}

#endif
