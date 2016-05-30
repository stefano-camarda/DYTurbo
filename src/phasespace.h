#ifndef phasespace_h
#define phasespace_h

//fortran interface
extern "C" {
  void sety_(double &yy);
  void setqt_(double &qtt);
  int binner_(double p3[4], double p4[4]);
}

namespace phasespace
{
  //boson integration boundaries
  extern double mmin;
  extern double mmax;
  extern double qtmin;
  extern double qtmax;
  extern double ymin;
  extern double ymax;
  extern void setbounds(double m1, double m2, double qt1, double qt2, double y1, double y2);

  //point in phase space
  extern double costh, m, qt, y;
#pragma omp threadprivate(costh, m, qt, y)

  extern void set_mqtycth(double M, double Qt, double Y, double Costh);
  extern void set_mqty(double M, double Qt, double Y);
  extern void set_m(double M);
  extern void set_qt(double Qt);
  extern void set_y(double Y);
  extern void set_cth(double Costh);
}

#endif
