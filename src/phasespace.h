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

  extern double cthmin;
  extern double cthmax;
  extern void setcthbounds(double cth1, double cth2);
  
  //point in phase space
  extern double m, qt, y, phiV;
#pragma omp threadprivate(m, qt, y, phiV)
  extern double costh;
#pragma omp threadprivate(costh)

  extern void set_mqtyphi(double M, double Qt, double Y, double PhiV = 0.);
  extern void set_m(double M);
  extern void set_qt(double Qt);
  extern void set_y(double Y);
  extern void set_phiV(double PhiV);
  
  extern void set_cth(double Costh);
}

#endif
