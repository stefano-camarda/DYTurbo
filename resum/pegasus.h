#ifndef pegasus_h
#define pegasus_h

#include "interface.h"

namespace pegasus
{
  //settings of the interface
  const int evolution_mode = 1; //mode of evolution: there are three available schemes for solving the evolution equations at NNLO. mode = 1 reproduces evolution in x-space
  const int alphas_steps = 20;  //number of steps for the Runge-Kutta integration of alphas beyond LO
 
  extern void init();
  extern void evolve();

  extern int nff;
  extern int ivfns;
}


//fortran interface to pegasus common blocks and functions
const int ndim = 512; //maximum dimension of the moments array

extern "C" {

  extern struct {
    int imodev_;
  } evmod_;

  extern struct {
    int npord_;
  } order_;

  extern struct {
    double logfr_;
  } frrat_;

  extern struct {
    int naord_;
    int nastps_;
  } aspar_;

  extern struct {
    int nflow_;
    int nfhigh_;
  } nfused_;

  extern struct {
    int nuord_;
  } itord_;
  
  extern struct {
    int nmax_;
  } nnused_;

  extern struct {
    fcomplex na_[ndim];
  } moms_;

  extern struct {
    fcomplex s_[6][ndim];
  } hsums_;
  
  extern struct {
    double zeta_[6];
  } rzeta_;
       
  extern struct {
    fcomplex d_[2][2];
  } kron2d_;

  extern struct {
    double cf_;
    double ca_;
    double tr_;
  } colour_;


  extern struct {
    double as0_;
    double m20_;
  } asinp_;

  extern struct {
    double asc_;
    double m2c_;
    double asb_;
    double m2b_;
    double ast_;
    double m2t_;
  } asfthr_;

  fcomplex psi_ (fcomplex &z);
  fcomplex dpsi_ (fcomplex &z, int &m);

  extern void betafct_();
  extern void pns0mom_();
  extern void psg0mom_();
  extern void lsgmom_();
  extern void pns1mom_();
  extern void psg1mom_();
  extern void usg1mom_();
  extern void usg1hmom_();
  extern void pns2mom_();
  extern void psg2mom_();
  extern void uns2mom_();
  extern void usg2mom_();
  extern void usg2hmom_();
  extern void ans2mom_();
  extern void asg2mom_();

  extern double as_ (double &R2, double &R20, double &AS0, int &NF);
  
  extern void evnfthr_(double &MC2, double &MB2, double &MT2);
  extern void evnvfn_(fcomplex PDFN[13][ndim], double &ASI, double &ASF, int &NF, int &NLOW, int &NHIGH, int &IPSTD);
  extern void evnffn_(fcomplex PDFN[13][ndim],  double &ASI,  double &ASF, int &NF, int &NLOW, int &NHIGH, int &IPSTD);
  extern void dyevnffn_(fcomplex PDFN[13][ndim],  double &ASI,  double &ASF, int &NF, int &NLOW, int &NHIGH, int &IPSTD);

  extern struct {
    fcomplex vai_[ndim];
    fcomplex m3i_[ndim];
    fcomplex m8i_[ndim];
    fcomplex sgi_[ndim];
    fcomplex p3i_[ndim];
    fcomplex p8i_[ndim];
    fcomplex gli_[ndim];
  } painp_;

  extern struct {
    fcomplex m15i_[ndim];
    fcomplex m24i_[ndim];
    fcomplex p15i_[ndim];
    fcomplex p24i_[ndim];
  } hfpainp_;

}


#endif
