#ifndef luminosity_h
#define luminosity_h

#include "parton.h"

extern "C"
{
  void luminosity_calc_();
  void luminosity_pdf1_(double & x);
  void luminosity_pdf2_(double & x);
}

namespace luminosity
{
  extern double sw,cw;

  extern double fh1[2*MAXNF+1];
  extern double fh2[2*MAXNF+1];
  extern double sumckm[MAXNF];
  
  extern double flqg[MAXNF];
  extern double flgq[MAXNF];

  extern double mesq[MAXNF][MAXNF];
  extern double mesq_nointerf[MAXNF][MAXNF];

  extern double sum_mesq;
  extern double sum_mesq_nointerf;
  
#pragma omp threadprivate(fh1,fh2,mesq,mesq_nointerf,sum_mesq,sum_mesq_nointerf,flqg,flgq)
  
  extern void init();
  extern void cache();
  extern void pdf1(double x);
  extern void pdf2(double x);
  extern void calc();
}

#endif
