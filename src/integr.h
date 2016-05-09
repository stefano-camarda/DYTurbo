#ifndef integr_h
#define integr_h

extern "C" {
  void cthmoments_(double &cthmom0, double &cthmom1, double &cthmom2);
  void sety_(double &yy);
  void setqt_(double &qtt);
  void genv4p_();
  int binner_(double p3[4], double p4[4]);
}

const int last_iter=4;
static int ICALL=0; // number of integrand calls
void tell_to_grid_we_are_alive();
extern double _costh, _m, _qt, _y;
#pragma omp threadprivate(_costh, _m, _qt, _y)

extern void genV4p(double m, double qt, double y, double phi);

extern void setcosth(double costh);
extern void setm(double mm);
extern void setqt(double qtt);
extern void sety(double yy);
extern void setmqty(double mm, double qtt, double yy);
extern void setcthmqty(double costh, double mm, double qtt, double yy);

//integration boundaries
extern void setbounds(double m1, double m2, double qt1, double qt2, double y1, double y2);
extern double mmin, mmax;
extern double qtmin, qtmax;
extern double ymin, ymax;

#endif
