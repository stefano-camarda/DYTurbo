#ifndef hcoeff_h
#define hcoeff_h

#include "mellinint.h"

#include <complex>

using namespace std;

namespace hcoeff
{
  extern complex <double> *Hqqb;
  extern complex <double> *Hqg;
  extern complex <double> *Hqq;
  extern complex <double> *Hqqp;
  extern complex <double> *Hgg;

  extern complex <double> *H1stqqb;
  extern complex <double> *H1stqg;
  extern complex <double> H1stgg;

  extern complex <double> *H2stqq;
  extern complex <double> *H2stqqp;
  extern complex <double> *H2stqqb;
  extern complex <double> *H2stqg;
  extern complex <double> *H2stgg;

  extern complex <double> *aexpqq;
  extern complex <double> *aexpqg;
#pragma omp threadprivate(Hqqb,Hqg,Hqq,Hqqp,Hgg,H1stgg,H1stqg,H1stqqb,H2stqq,H2stqqp,H2stqqb,H2stqg,H2stgg,aexpqq,aexpqg)


  extern complex <double> Hqqbz;
  extern complex <double> Hqgz;
  extern complex <double> Hqqz;
  extern complex <double> Hqqpz;
  extern complex <double> Hggz;
#pragma omp threadprivate(Hqqbz,Hqgz,Hqqz,Hqqpz,Hggz)

  extern void allocate();
  // q2-dependent quantities
  extern void calc(double aass, complex <double> logmuf2q2, complex <double> logq2muf2, complex <double> logq2mur2, complex <double> loga);

  // b-dependent quantities
  extern void calcb(double aass, complex <double> logmuf2q2, complex <double> loga, complex <double> alpq, complex <double> aexp, complex <double> aexpb);

  extern void truncate();
    
  extern void invert(double q2);
  
  extern void free();
}

#endif
