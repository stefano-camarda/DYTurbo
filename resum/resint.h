#ifndef resint_h
#define resint_h

#include "fcomplex.h"

#include <complex>

using namespace std;

//workspace for intdeo DEQUAD integration
extern double tiny;
const int lenaw = 8000;
extern double aw[lenaw];

extern "C"
{

  //real or complex b integration
  extern struct {
    int flagrealcomplex_;
  } flagrealcomplex_;

  //QCD order in the sudakov
  extern struct {
    int flag1_;
  } flag1_;

  //QCD order in the alphasl evolution
  extern struct {
    int iord_;
  } iorder_;
  
  //normal or modified sudakov
  extern struct {
    int imod_;
  } modified_;

  /*
  extern struct {
    double a_param_;
    double b0p_;
  } a_param_;
  */

  extern struct {
    double aass_;
  } aass_;

  extern struct {
      double rblim_;
      fcomplex cblim_;
  } blimit_;

  extern struct {
    double qt_;
    double qt2_;
    double xtau_;
    double q_;
    double q2_;
    double shad_;
    double sroot_;
    double mur_;
    double mur2_;
    double muf_;
    double muf2_;
  } scaleh_;

    extern struct {
      double rloga_;
      double rlogq2mur2_;
  } rlogs_;
}
#pragma omp threadprivate(blimit_,rlogs_,aass_,scaleh_)

namespace resint
{
  //point in phase space
  extern double _qt, _y, _m, _costh;
  extern double tau, x1, x2;
  extern int _mode;
#pragma omp threadprivate(_qt,_m,_y,_costh,_mode,tau,x1,x2)

  //scales
  extern double muren, mufac, mures;
  extern double muren2, mufac2, mures2;
  extern double a;
#pragma omp threadprivate(muren,mufac,mures,muren2,mufac2,mures2,a)

  //log of scales
  extern complex <double> loga;
  extern double rloga;
  extern complex <double> logmuf2q2;
  extern complex <double> logq2muf2;
  extern complex <double> logq2mur2;
  extern double rlogq2mur2;
#pragma omp threadprivate(loga,rloga,logmuf2q2,logq2muf2,logq2mur2,rlogq2mur2)
    
  //alphas
  extern double alpqr;
  extern double alpqf;
  extern double aass;

  extern double alpqfac;
  extern double alpqren;
  extern double alpqres;
#pragma omp threadprivate(alpqr,alpqf,aass,alpqfac,alpqren,alpqres)

  //integration contour
  extern double bc;
#pragma omp threadprivate(bc)
  
  extern void init();
  extern double rint(double costh, double m, double qt, double y, int mode);
  extern double bintegral(double qt);
}

#endif
