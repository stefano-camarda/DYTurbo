#ifndef omegaintegr_h
#define omegaintegr_h

#include "cuba.h"
#include "phasespace.h"

#include <vector>

using namespace std;

//fortran interface
extern "C" {
  void cthmoments_(double &cthmom0, double &cthmom1, double &cthmom2);
  void genv4p_();
  void genrfaxes_();
}

namespace omegaintegr
{
  extern void genV4p();
  extern void genRFaxes();
  extern void genzaxisCS();

  //  extern void genl4p(float costh, float phi_lep);
  extern void genl4p(double costh, double phi_lep);
  
  //CS framework
  extern double kap1[4];
  extern double xax[3];
  extern double yax[3];
  extern double zax[3];
#pragma omp threadprivate(kap1,xax,yax,zax)

  double costhCS();

  double costh_qtrec();

  //integration functions
  void costhbound(double phi_lep, vector<double> &min, vector<double> &max);
  void cthmoments(double &cthmom0, double &cthmom1, double &cthmom2);
  integrand_t thphiintegrand(const int &ndim, const double x[], const int &ncomp, double f[]);
}

#endif
