#ifndef integr_h
#define integr_h
#include <iostream>
#include "cuba.h"


extern "C" {
  void cthmoments_(double &cthmom0, double &cthmom1, double &cthmom2);
  void sety_(double &yy);
  void genv4p_();
}


extern double _costh, _m, _qt, _y;
extern double resy(double y, void* param = NULL);
extern double resth(double costh, void* param = NULL);
extern double resm(double mm, void* param = NULL);
extern void set(double costh, double mm, double qtt, double yy);
extern void set(double mm, double qtt, double yy);
extern void genV4p(double m, double qt, double y, double phi);

extern void setcosth(double costh);
extern void setm(double mm);
extern void sety(double yy);
extern void setqt(double qtt);
/*
inline void setcosth(double costh) {_costh = costh;}
inline void setm(double mm) {_m = mm;}
inline void sety(double yy) {_y = yy;}
inline void setqt(double qtt) {_qt = qtt;}
*/

//integration boundaries
extern void setbounds(double m1, double m2, double qt1, double qt2, double y1, double y2);
extern double mmin, mmax;
extern double qtmin, qtmax;
extern double ymin, ymax;

integrand_t resintegrand2d(const int &ndim, const double x[], const int &ncomp, double f[]);
integrand_t resintegrand3d(const int &ndim, const double x[], const int &ncomp, double f[]);
integrand_t resintegrand4d(const int &ndim, const double x[], const int &ncomp, double f[],
			   void* userdata, const int &nvec, const int &core,
			   double &weight, const int &iter);

using namespace std;

/*
class FunctResm
{
public:
	double operator()(double x) const
	{
	  return resm(x);
	}
};
*/

#endif
