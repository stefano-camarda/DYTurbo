#ifndef hankel_h
#define hankel_h

namespace hankel
{

  extern int nu;
  extern int n;
  extern double h;

  extern double *zeros;
  extern double *w;
  extern double *x;
  extern double *j;
  extern double *dp;

  extern double *series;

  extern int x_power;
  extern int k_power;

  extern void init(int Nu, int N, double H=0.1);
  extern void roots();
  extern double psi(double t);
  extern double dpsi(double t);
  extern void fillx();
  extern double jn(double x);
  extern double fillj();
  extern double filldpsi();
  extern double jn1(double x);
  extern double yn(double x);
  extern double weight();
  extern double get_series(double (*f)(double), double k);
  extern void transform(double (*f)(double), double k, double &res, double &err);
  extern void free();
}

#endif
