#ifndef omegaintegr_h
#define omegaintegr_h

#include "cuba.h"

#include <vector>

using namespace std;

//fortran interface
extern "C" {
  void cthmoments_(double &cthmom0, double &cthmom1, double &cthmom2);
  void genv4p_();
}

namespace omegaintegr
{
  //Vector boson 4-momentum and boost
  extern double pV[4];
  extern double gam;
  extern double beta[3];
#pragma omp threadprivate(pV,gam,beta)

  extern void genV4p();
  
  //leptons 4-momenta
  extern double p4[4];
  extern double p3[4];
#pragma omp threadprivate(p4,p3)

  //  extern void genl4p(float costh, float phi_lep);
  extern void genl4p(double costh, double phi_lep);
  extern void getp3(double p[4]);
  extern void getp4(double p[4]);
  
  //CS framework
  extern double kap1[4];
  extern double xax[3];
  extern double yax[3];
  extern double zax[3];
#pragma omp threadprivate(kap1,xax,yax,zax)

  double costhCS();
  
  //integration functions
  void costhbound(double phi_lep, vector<double> &min, vector<double> &max);
  void cthmoments(double &cthmom0, double &cthmom1, double &cthmom2);
  integrand_t thphiintegrand(const int &ndim, const double x[], const int &ncomp, double f[]);
}

#endif
