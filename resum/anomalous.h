#ifndef anomalous_h
#define anomalous_h
#include "mellinint.h"

namespace anomalous
{

  extern void init();
  extern void allocate();
  extern void calc();
  extern void free();
  extern void release();

  //anomalous dimensions
  extern complex <double> *ans,*am,*ap,*al,*be,*ab;
  extern complex <double> *rmin,*rplus,*rqq,*rqg,*rgq,*rgg;
  extern complex <double> *RMMQQ,*RMMQG,*RMMGQ,*RMMGG,*RMPQQ,*RMPQG,*RMPGQ,*RMPGG,*RPMQQ,*RPMQG,*RPMGQ,*RPMGG,*RPPQQ,*RPPQG,*RPPGQ,*RPPGG;
  //#pragma omp threadprivate()
  
  //gamma1 and gamma2
  extern complex <double> *gamma1qq,*gamma1qg,*gamma1gq,*gamma1gg;
  extern complex <double> *gamma2qq,*gamma2qqV,*gamma2qqbV,*gamma2qqS,*gamma2qqbS,*gamma2qg,*gamma2gq,*gamma2gg;
  //with Mellin update all this should become threadprivate!!! --> Also need to allocate and free dynamically at each iteration
  //#pragma omp threadprivate()

  //C1 and C2 coefficients
  extern complex <double> *C1QQ,*C1QG,*C1GQ;
  extern complex <double> C1GG;
  extern complex <double> *C2qgM,*C2NSqqM,*C2SqqbM,*C2NSqqbM;
  //#pragma omp threadprivate()
  
  //mellin2d
  extern complex <double> *C1QQ_1,*C1QG_1,*C1GQ_1;
  extern complex <double> *C1QQ_2,*C1QG_2,*C1GQ_2;
  extern complex <double> *C2qgM_1,*C2NSqqM_1,*C2SqqbM_1,*C2NSqqbM_1;
  extern complex <double> *C2qgM_2,*C2NSqqM_2,*C2SqqbM_2,*C2NSqqbM_2;
  //#pragma omp threadprivate()
  
  inline int index(int i, int sign)
  {return i + sign*mellinint::mdim;}
}

#endif
