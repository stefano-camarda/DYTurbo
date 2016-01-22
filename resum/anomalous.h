#ifndef anomalous_h
#define anomalous_h
#include "mellinint.h"

namespace anomalous
{

  extern void init();
  extern complex <double> *ans,*am,*ap,*al,*be,*ab;
  extern complex <double> *rmin,*rplus,*rqq,*rqg,*rgq,*rgg;
  extern complex <double> *RMMQQ,*RMMQG,*RMMGQ,*RMMGG,*RMPQQ,*RMPQG,*RMPGQ,*RMPGG,*RPMQQ,*RPMQG,*RPMGQ,*RPMGG,*RPPQQ,*RPPQG,*RPPGQ,*RPPGG;
  extern complex <double> *C1QQ,*C1QG,*C1GQ;
  extern complex <double> C1GG;
  extern complex <double> *gamma1qq,*gamma1qg,*gamma1gq,*gamma1gg;
  extern complex <double> *gamma2qq,*gamma2qqV,*gamma2qqbV,*gamma2qqS,*gamma2qqbS,*gamma2qg,*gamma2gq,*gamma2gg;
  extern complex <double> *C2qgM,*C2NSqqM,*C2SqqbM,*C2NSqqbM;

  inline int index(int i, int sign)
  {return i + sign*mellinint::mdim;}
}

#endif
