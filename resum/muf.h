#ifndef muf_h
#define muf_h

#include <complex>

using namespace std;

namespace muf
{
//  extern complex <double> *qqb;
//  extern complex <double> *qg;
//  extern complex <double> *qq;
//  extern complex <double> *qqp;
//  extern complex <double> *gg;
//  extern complex <double> *qqbp;
//  extern complex <double> *qbg;
//  extern complex <double> *qpg;
//  extern complex <double> *qbpg;
//#pragma omp threadprivate(qqb,qg,qq,qqp,qqbp,gg,qbg,qpg,qbpg)

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
  
  void allocate();
  void reset();
  void calc(complex <double> b = 0.);
  void free();

}

#endif
