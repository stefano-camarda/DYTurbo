#ifndef blim_h
#define blim_h

namespace blim
{
  extern double pdf, sudakov, expc;
#pragma omp threadprivate (pdf,sudakov,expc)

  void set();
  double calc(double blm, double Q);
}

#endif
