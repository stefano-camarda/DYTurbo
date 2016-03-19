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

  //non perturbative g
  extern struct {
    double g_;
  } np_;

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

namespace resint
{
  //point in phase space
  extern double _qt, _y, _m, _costh;
  extern int _mode;

  //scales
  extern double mur, muf;
  extern double mur2, muf2;

  //log of scales
  extern complex <double> loga;
  extern double rloga;
  extern complex <double> logmuf2q2;
  extern complex <double> logq2muf2;
  extern complex <double> logq2mur2;
  extern double rlogq2mur2;
    
  //alphas
  extern double alpqr;
  extern double alpqf;
  extern double aass;
  
  extern void init();
  extern double rint(double costh, double m, double qt, double y, int mode);
  extern double bintegral(double qt);
}

#endif
