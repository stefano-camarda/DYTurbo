#ifndef pmom_h
#define pmom_h
#include "mellinint.h"

namespace pmom
{

  extern void init();
  extern void allocate();
  extern void calc();
  extern void release();
  extern void free();
  extern complex <double> *gamma1qq,*gamma1qqb,*gamma1qqp,*gamma1qqbp,*gamma1qg,*gamma1gq,*gamma1gg;
  extern complex <double> *gamma2qq,*gamma2qqb,*gamma2qqp,*gamma2qqbp,*gamma2qg,*gamma2gq,*gamma2gg;
  extern complex <double> *gamma3qq,*gamma3qqb,*gamma3qqp,*gamma3qqbp,*gamma3qg,*gamma3gq,*gamma3gg;
#pragma omp threadprivate(gamma1qq,gamma1qqb,gamma1qqp,gamma1qqbp,gamma1qg,gamma1gq,gamma1gg)
#pragma omp threadprivate(gamma2qq,gamma2qqb,gamma2qqp,gamma2qqbp,gamma2qg,gamma2gq,gamma2gg)
#pragma omp threadprivate(gamma3qq,gamma3qqb,gamma3qqp,gamma3qqbp,gamma3qg,gamma3gq,gamma3gg)
  
  extern complex <double> *gamma1qq_1,*gamma1qqb_1,*gamma1qqp_1,*gamma1qqbp_1,*gamma1qg_1,*gamma1gq_1,*gamma1gg_1;
  extern complex <double> *gamma2qq_1,*gamma2qqb_1,*gamma2qqp_1,*gamma2qqbp_1,*gamma2qg_1,*gamma2gq_1,*gamma2gg_1;
  extern complex <double> *gamma3qq_1,*gamma3qqb_1,*gamma3qqp_1,*gamma3qqbp_1,*gamma3qg_1,*gamma3gq_1,*gamma3gg_1;
  extern complex <double> *gamma1qq_2,*gamma1qqb_2,*gamma1qqp_2,*gamma1qqbp_2,*gamma1qg_2,*gamma1gq_2,*gamma1gg_2;
  extern complex <double> *gamma2qq_2,*gamma2qqb_2,*gamma2qqp_2,*gamma2qqbp_2,*gamma2qg_2,*gamma2gq_2,*gamma2gg_2;
  extern complex <double> *gamma3qq_2,*gamma3qqb_2,*gamma3qqp_2,*gamma3qqbp_2,*gamma3qg_2,*gamma3gq_2,*gamma3gg_2;
#pragma omp threadprivate(gamma1qq_1,gamma1qqb_1,gamma1qqp_1,gamma1qqbp_1,gamma1qg_1,gamma1gq_1,gamma1gg_1)
#pragma omp threadprivate(gamma2qq_1,gamma2qqb_1,gamma2qqp_1,gamma2qqbp_1,gamma2qg_1,gamma2gq_1,gamma2gg_1)
#pragma omp threadprivate(gamma3qq_1,gamma3qqb_1,gamma3qqp_1,gamma3qqbp_1,gamma3qg_1,gamma3gq_1,gamma3gg_1)
#pragma omp threadprivate(gamma1qq_2,gamma1qqb_2,gamma1qqp_2,gamma1qqbp_2,gamma1qg_2,gamma1gq_2,gamma1gg_2)
#pragma omp threadprivate(gamma2qq_2,gamma2qqb_2,gamma2qqp_2,gamma2qqbp_2,gamma2qg_2,gamma2gq_2,gamma2gg_2)
#pragma omp threadprivate(gamma3qq_2,gamma3qqb_2,gamma3qqp_2,gamma3qqbp_2,gamma3qg_2,gamma3gq_2,gamma3gg_2)
  
  inline int index(int i, int sign)
  {return i + sign*mellinint::mdim;}
}

#endif
