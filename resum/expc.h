#ifndef expc_h
#define expc_h

//#include "mellinint.h"

#include <complex>

using namespace std;

namespace expc
{
  extern complex <double> *aexpqq;
  extern complex <double> *aexpqg;
  extern complex <double> *aexpqqb;
  extern complex <double> *aexpqqp;
  extern complex <double> *aexpqqbp;
  extern complex <double> *aexpqq_1;
  extern complex <double> *aexpqg_1;
  extern complex <double> *aexpqqb_1;
  extern complex <double> *aexpqqp_1;
  extern complex <double> *aexpqqbp_1;
  extern complex <double> *aexpqq_2;
  extern complex <double> *aexpqg_2;
  extern complex <double> *aexpqqb_2;
  extern complex <double> *aexpqqp_2;
  extern complex <double> *aexpqqbp_2;
#pragma omp threadprivate(aexpqq,aexpqg,aexpqqb,aexpqqp,aexpqqbp,aexpqq_1,aexpqg_1,aexpqqb_1,aexpqqp_1,aexpqqbp_1,aexpqq_2,aexpqg_2,aexpqqb_2,aexpqqp_2,aexpqqbp_2)
  
  extern complex <double> *qqb;
  extern complex <double> *qg;
  extern complex <double> *qg_1;
  extern complex <double> *qg_2;
  extern complex <double> *qq;
  extern complex <double> *qq_1;
  extern complex <double> *qq_2;
  extern complex <double> *qqp;
  extern complex <double> *qqp_1;
  extern complex <double> *qqp_2;
  extern complex <double> *qqbp;
  extern complex <double> *qqbp_1;
  extern complex <double> *qqbp_2;
  extern complex <double> *gg;
  extern complex <double> *qbg;
  extern complex <double> *qbg_1;
  extern complex <double> *qbg_2;
  extern complex <double> *qpg;
  extern complex <double> *qpg_1;
  extern complex <double> *qpg_2;
  extern complex <double> *qbpg;
  extern complex <double> *qbpg_1;
  extern complex <double> *qbpg_2;
#pragma omp threadprivate(qqb,qg,qg_1,qg_2,qq,qq_1,qq_2,qqp,qqp_1,qqp_2,qqbp,qqbp_1,qqbp_2,gg,qbg,qbg_1,qbg_2,qpg,qpg_1,qpg_2,qbpg,qbpg_1,qbpg_2)

  //  extern complex <double> aexp;
  //  extern complex <double> aexpB;
  //#pragma omp threadprivate(aexp)
  
  extern void allocate();
  extern void reset();
  extern void calc(complex <double> b);
  extern void free();
}

#endif
