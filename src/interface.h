#ifndef interface_h
#define interface_h

#include "mcfm_interface.h"
#include "fcomplex.h"

#include <complex>

using namespace std;

extern "C" {
  // dyres rewritten functions
  double resumm_(double &costh, double &mm, double &qtt, double &yy, int& mode);
  //void setup_();
  void dyinit_();
  void pdfini_();
  void iniflavreduce_();
  void gaussinit_();
  double dyalphas_mcfm_(double &q, double &amz, int &nloop);
  double dyalphas_lhapdf_(double &q);
  int cuts_(double p[4][12], int &njet);
  //int cutsold_(double p[4][12], int &njet);
  void rescinit_();

  //fortran interface for C++ rewritten resummation
  void setmesq_expy_(int& mode, double& m, double& costh, double& y);
  void pdfevol_(int& i1, int& i2, int& sign);
  void mellinint_pdf_mesq_expy_(int& i1, int& i2, int& sign);
  fcomplex mellinint_integrand_(int& i1, int& i2, int& sign);
  void hcoeff_calc_(double& aass, double& logmuf2q2, double& logq2muf2, double& logq2mur2, double& loga);
  void hcoeff_calcb_(double& aass, double& logmuf2q2, double& loga, double& alpq, double &aexp, double &aexpb);
  
  void dycoupling_();

  void rapintegrals_(double &ymin,double &ymax, double& mass, int& nocuts);
  void cacheyrapint_(double &ymin,double &ymax);

  void ctqtint_(double &m, double &y, double &qtmin, double &qtmax);

  void initmoments_();
  // fortran common spaces

  //Interface for anomalous dimensions and Wilson coefficients
  void ancalc_(fcomplex &QQI, fcomplex &QGF, fcomplex &GQI, fcomplex &GGI, fcomplex &GGF, fcomplex &NS1MI, fcomplex &NS1PI, fcomplex &NS1F,
	       fcomplex &QQ1F, fcomplex &QG1F, fcomplex &GQ1I, fcomplex &GQ1F, fcomplex &GG1I, fcomplex &GG1F, fcomplex &xn);
  void anom_(fcomplex &ANS, fcomplex &AM, fcomplex &AP, fcomplex &AL, fcomplex &BE, fcomplex &AB, fcomplex &RMIN, fcomplex &RPLUS, fcomplex &RQQ, fcomplex &RQG,
	     fcomplex &RGQ, fcomplex &RGG, fcomplex &C2Q, fcomplex &C2G, fcomplex &CDYQ, fcomplex &CDYG, fcomplex &xn, int &FR,
	     fcomplex &QQI, fcomplex &QGF, fcomplex &GQI, fcomplex &GGI, fcomplex &GGF, fcomplex &NS1MI, fcomplex &NS1PI, fcomplex &NS1F,
	     fcomplex &QQ1F, fcomplex &QG1F, fcomplex &GQ1I, fcomplex &GQ1F, fcomplex &GG1I, fcomplex &GG1F, fcomplex &C2QI, fcomplex &C2GF,
	     fcomplex &CDYQI, fcomplex &CDYGI);
  void h2calc_(fcomplex &C2qg, fcomplex &C2NSqqb, fcomplex &C2NSqq, fcomplex &C2Sqqb,fcomplex &xn);

  // controls
  extern struct {
    int approxpdf_;
    int pdfintervals_;
    int fixedorder_;
  } opts_;

  extern struct {
    int pdferr_;
    int totpdf_;
  } pdferropts_;

  extern struct {
    double kmuren_;
    double kmufac_;
    double kmures_;
  } scaleopts_;

  extern struct {
    double amz_;
  } couple_;

  /*extern struct {
    int nlooprun_;
    } nlooprun_;*/

  /*  extern struct {
    int lhapdfs_;
    } lhapdfs_;*/

  //initialization flag
  extern struct {
    int flag_;
  } flag_;

  extern struct {
    double a_param_;
    double b0p_;
  } a_param_;
  
  //non perturbative g
  extern struct {
    double g_param_;
  } g_param_;

  extern struct {
    double g_;
  } np_;

  extern struct {
    int order_;
  } nnlo_;

  extern struct {
    double qtcut_;
    double xqtcut_;
  } qtsub_;

  extern struct {
    int doFill_;
  } dofill_;

  extern struct {
    double sigmaij_[11][11];
  } sigmaij_;

  void hists_setpdf_(int * npdf);
  void hists_fill_(double p3[4], double p4[4], double *weight);
  void hists_AiTest_(double pjet[4][12], double p4cm[4],double *Q,double *qt,double *y,double* pcosthCS, double* pphiCS, double *pphiVB, double *wt, double *loHst );
  void hists_real_dipole_(double p3[4], double p4[4], double *weight,int *nd);
  void hists_real_event_();
  void hists_fill_pdf_(double p3[4], double p4[4], double *weight, int *npdf);
  void hists_real_dipole_pdf_(double p3[4], double p4[4], double *weight,int *nd, int *npdf);
  void hists_real_event_pdf_(int* npdf);
}
#pragma omp threadprivate(a_param_,scale_,facscale_,qcdcouple_,sigmaij_)

#endif
