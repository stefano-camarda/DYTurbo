#ifndef qtint_h
#define qtint_h


namespace qtint
{
  extern double LL1_mesqij[12];
  extern double LL2_mesqij[12];
  extern double LL3_mesqij[12];
  extern double LL4_mesqij[12];
#pragma omp threadprivate(LL1_mesqij,LL2_mesqij,LL3_mesqij,LL4_mesqij)

  extern void calc(double m, double qtmin, double qtmax, int mode);
}

#endif
