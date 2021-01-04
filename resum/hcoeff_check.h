#ifndef hcoeff_check_h
#define hcoeff_check_h

#include <complex>

using namespace std;

namespace hcoeff_check
{
  extern complex <double> *Hqqb;
  extern complex <double> *Hqg;
  extern complex <double> *Hqq;
  extern complex <double> *Hqqp;
  extern complex <double> *Hqqbp;
  extern complex <double> *Hgg;
  extern complex <double> *Hqbg;
  extern complex <double> *Hqpg;
  extern complex <double> *Hqbpg;
#pragma omp threadprivate(Hqqb,Hqg,Hqq,Hqqp,Hqqbp,Hgg,Hqbg,Hqpg,Hqbpg)
  
  extern void allocate();
  extern void reset();
  extern void calc();
  extern void free();
}

#endif
