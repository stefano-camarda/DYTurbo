#ifndef hcoefficients_h
#define hcoefficients_h

#include "mellinint.h"

#include <complex>

using namespace std;

namespace hcoefficients
{
  extern complex <double> *Hqqb;
  extern complex <double> *Hqg;
  extern complex <double> *Hqg_1;
  extern complex <double> *Hqg_2;
  extern complex <double> *Hqq_nnll;
  extern complex <double> *Hqq;
  extern complex <double> *Hqq_1;
  extern complex <double> *Hqq_2;
  extern complex <double> *Hqqp;
  extern complex <double> *Hqqp_1;
  extern complex <double> *Hqqp_2;
  extern complex <double> *Hgg;

  extern complex <double> *H1stqqb;
  extern complex <double> *H1stqg;
  extern complex <double> H1stgg;

  extern complex <double> *H2stqq;
  extern complex <double> *H2stqqp;
  extern complex <double> *H2stqqb;
  extern complex <double> *H2stqg_1;
  extern complex <double> *H2stqg_2;
  extern complex <double> *H2stgg;

  extern complex <double> *aexpqq;
  extern complex <double> *aexpqg;
      
  extern void init();
  // q2-dependent quantities
  extern void calc(double aass, complex <double> logmuf2q2, complex <double> logq2muf2, complex <double> logq2mur2, complex <double> loga);

  // b-dependent quantities
  extern void calcb(double aass, complex <double> logmuf2q2, complex <double> loga, complex <double> alpq, complex <double> aexp, complex <double> aexpb);

  inline int index(int i, int sign)
  {return i + mellinint::mdim*sign;}

  inline int index(int i1, int i2, int sign)
  {return i1 + mellinint::mdim*(i2 + mellinint::mdim*sign);}

}

#endif
