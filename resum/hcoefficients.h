#ifndef hcoefficients_h
#define hcoefficients_h

#include "mellinint.h"

#include <complex>

using namespace std;

namespace hcoefficients
{
  extern complex <double> *Hqqb;
  extern complex <double> *Hqg_1;
  extern complex <double> *Hqg_2;
  extern complex <double> *Hqq_1;
  extern complex <double> *Hqq_2;
  extern complex <double> *Hqqp_1;
  extern complex <double> *Hqqp_2;
  extern complex <double> *Hgg;
#pragma omp threadprivate(Hqqb,Hqg_1,Hqg_2,Hqq_1,Hqq_2,Hqqp_1,Hqqp_2,Hgg)
  
//  extern complex <double> *H1stqqb;
//  extern complex <double> *H1stqg;
//#pragma omp threadprivate(H1stqqb,H1stqg)
//
//  extern complex <double> *H2stqqb;
//  extern complex <double> *H2stqg;
//  extern complex <double> *H2stqq;
//  extern complex <double> *H2stqqp;
//#pragma omp threadprivate(H2stqqb,H2stqg,H2stqq,H2stqqp)
//
//  extern complex <double> *aexpqq;
//  extern complex <double> *aexpqg;
//#pragma omp threadprivate(aexpqq,aexpqg)

  extern void allocate();
  extern void reset();
  // q2-dependent quantities
  extern void calc(double aass, complex <double> logmuf2q2, complex <double> logq2muf2, complex <double> logq2mur2, complex <double> loga);

  // b-dependent quantities
  extern void calcb(double aass, complex <double> logmuf2q2, complex <double> loga, complex <double> alpq, complex <double> aexp, complex <double> aexpb);

  extern void free();

  inline int index(int i, int sign)
  {return i + mellinint::mdim*sign;}

  inline int index(int i1, int i2, int sign)
  {return i2 + mellinint::mdim*(i1 + mellinint::mdim*sign);}

}

#endif
