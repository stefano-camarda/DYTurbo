#ifndef pegasus_h
#define pegasus_h

#include "interface.h"

namespace pegasus
{
  //settings of the interface
  const int evolution_mode = 1; //mode of evolution: there are three available schemes for solving the evolution equations at NNLO. mode = 1 reproduces evolution in x-space
  const int alphas_steps = 10;  //number of steps for the Runge-Kutta integration of alphas beyond LO

  const int NFMIN = 3;
  const int NFMAX = 6;
  
  extern void init_const();  //Initialise mathematical and QCD constants

  extern void init();        //Main init which calls allocate(), init_alphas(), calc_mellin(), and init_pdf()
  extern void allocate();
  extern void init_alphas();
  extern void calc_mellin();
  extern void init_pdf();
  extern void free();
  extern void release();
  
  extern void update();
  extern void store();
  extern void retrieve();
  extern void evolve();
  extern double alphas(double M2, double R2, double &ASI, int &NF);

  extern int nff;
  extern int ivfns;
  extern int dim;

  //Decomposed PDFs at the starting scale
  extern complex <double> *gli;
  extern complex <double> *vai;
  extern complex <double> *m3i;
  extern complex <double> *m8i;
  extern complex <double> *m15i;
  extern complex <double> *m24i;
  extern complex <double> *sgi;
  extern complex <double> *p3i;
  extern complex <double> *p8i;
  extern complex <double> *p15i;
  extern complex <double> *p24i;
#pragma omp threadprivate(gli,vai,m3i,m8i,m15i,m24i,sgi,p3i,p8i,p15i,p24i) //--> make thread safe!, but copyin for evolmode = 3
}


//fortran interface to pegasus common blocks and functions //--> move inside the namespace
const int ndim = 128; //512; //maximum dimension of the moments array
const int nfl = pegasus::NFMAX-pegasus::NFMIN+1;    //= 4: values of nf (3,4,5,6)
const int numax = 20;

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
#pragma omp threadprivate(moms_)
  
  extern struct {
    fcomplex s_[6][ndim];
  } hsums_;
#pragma omp threadprivate(hsums_)
  
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
    double pgbeta0_[nfl];
    double pgbeta1_[nfl];
    double pgbeta2_[nfl];
    double pgbeta3_[nfl];
  } pgbeta_;

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
  extern void asg2mom_form_();

  extern double as_ (double &R2, double &R20, double &AS0, int &NF);
  
  extern void evnfthr_(double &MC2, double &MB2, double &MT2);
  extern void evnasthr_(double &MC2, double &MB2, double &MT2);
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

  extern struct {
    fcomplex vac_[ndim];
    fcomplex m3c_[ndim];
    fcomplex m8c_[ndim];
    fcomplex m15c_[ndim];
    fcomplex sgc_[ndim];
    fcomplex p3c_[ndim];
    fcomplex p8c_[ndim];
    fcomplex p15c_[ndim];
    fcomplex glc_[ndim];
  } pacthr_;

  extern struct {
    fcomplex vab_[ndim];
    fcomplex m3b_[ndim];
    fcomplex m8b_[ndim];
    fcomplex m15b_[ndim];
    fcomplex sgb_[ndim];
    fcomplex p3b_[ndim];
    fcomplex p8b_[ndim];
    fcomplex p15b_[ndim];
    fcomplex m24b_[ndim];
    fcomplex p24b_[ndim];
    fcomplex glb_[ndim];
  } pabthr_;

  extern struct {
    fcomplex vat_[ndim];
    fcomplex m3t_[ndim];
    fcomplex m8t_[ndim];
    fcomplex m15t_[ndim];
    fcomplex sgt_[ndim];
    fcomplex p3t_[ndim];
    fcomplex p8t_[ndim];
    fcomplex p15t_[ndim];
    fcomplex m24t_[ndim];
    fcomplex p24t_[ndim];
    fcomplex m35t_[ndim];
    fcomplex p35t_[ndim];
    fcomplex glt_[ndim];
  } patthr_;
#pragma omp threadprivate(painp_,hfpainp_,pacthr_,pabthr_,patthr_)

  //LO QCD splitting functions
  extern struct {
    fcomplex p0ns_[nfl][ndim];
  } pns0_;

  extern struct {
    fcomplex p0sg_[2][2][nfl][ndim];
  } psg0_;

  //NLO QCD splitting functions
  extern struct {
    fcomplex p1ns_[3][nfl][ndim];
  } pns1_;

  extern struct {
    fcomplex p1sg_[2][2][nfl][ndim];
  } psg1_;

  //NNLO QCD splitting functions
  extern struct {
    fcomplex p2ns_[3][nfl][ndim];
  } pns2_;

  extern struct {
    fcomplex p2sg_[2][2][nfl][ndim];
  } psg2_;
#pragma omp threadprivate(pns0_,psg0_,pns1_,psg1_,pns2_,psg2_)

  extern struct {
    fcomplex a2sg_[2][2][ndim];
  } asg2_;
#pragma omp threadprivate(asg2_)

  extern struct {
    fcomplex sschlp_[ndim];
    fcomplex sstr2p_[ndim];
    fcomplex sstr3p_[ndim];
  } spsums_;
#pragma omp threadprivate(spsums_)
  
  extern struct {
    fcomplex r_[2][nfl][ndim];
    fcomplex e_[2][2][2][nfl][ndim];
  } lsg_;
#pragma omp threadprivate(lsg_)

  extern struct {
    fcomplex u1_[2][2][nfl][ndim];
  } u1sg_;
  extern struct {
    fcomplex r1_[2][2][nfl][ndim];
  } r1sg_;
  extern struct {
    fcomplex u1h_[2][2][nfl][ndim][numax];
  } u1hsg_;
  extern struct {
    fcomplex u2_[2][2][nfl][ndim];
  } u2sg_;
  extern struct {
    fcomplex r2_[2][2][nfl][ndim];
  } r2sg_;
  extern struct {
    fcomplex u2h_[2][2][nfl][ndim][numax];
  } u2hsg_;
  extern struct {
    fcomplex uns2_[3][nfl][ndim][numax];
  } u2ns_;
#pragma omp threadprivate(u1sg_,r1sg_,u1hsg_,u2sg_,r2sg_,u2hsg_,u2ns_)
  
}



#endif
