#ifndef interface_h
#define interface_h

#include "rapint.h"
#include "fcomplex.h"

#include <complex>

#define MAXNF 5

using namespace std;

extern "C" {

  // rewritten functions
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

  //C++ rewritten resummation
  inline void rapint_cache_(double& ymin, double& ymax) {rapint::cache(ymin, ymax);};
  inline void rapint_integrate_(double& ymin, double& ymax, double& m) {rapint::integrate(ymin, ymax, m);};
  void setmesq_expy_(int& mode, double& m, double& costh, double& y);
  void pdfevol_(int& i1, int& i2, int& sign);
  void mellinint_pdf_mesq_expy_(int& i1, int& i2, int& sign);
  fcomplex mellinint_integrand_(int& i1, int& i2, int& sign);
  void hcoeff_calc_(double& aass, double& logmuf2q2, double& logq2muf2, double& logq2mur2, double& loga);
  void hcoeff_calcb_(double& aass, double& logmuf2q2, double& loga, double& alpq, double &aexp, double &aexpb);
  
  void breitw_(double& x1, double& mminsq, double& mmaxsq, double& rmass, double& rwidth, double& msq, double& wt);
  void boost_(double& mass, double p1[],double p_in[], double p_out[]);
  void branch_(double& brwen, double& brzee, double& brtau, double& brtop);
  void ckmfill_(int& nwz);

  void dycoupling_();

  void rapintegrals_(double &ymin,double &ymax, double& mass, int& nocuts);
  void cacheyrapint_(double &ymin,double &ymax);

  void ctqtint_(double &m, double &y, double &qtmin, double &qtmax);

  void initmoments_();
  // fortran common spaces

  //Anomalous dimensions and Wilson coefficients
  
  void ancalc_(fcomplex &QQI, fcomplex &QGF, fcomplex &GQI, fcomplex &GGI, fcomplex &GGF, fcomplex &NS1MI, fcomplex &NS1PI, fcomplex &NS1F,
	       fcomplex &QQ1F, fcomplex &QG1F, fcomplex &GQ1I, fcomplex &GQ1F, fcomplex &GG1I, fcomplex &GG1F, fcomplex &xn);
  void anom_(fcomplex &ANS, fcomplex &AM, fcomplex &AP, fcomplex &AL, fcomplex &BE, fcomplex &AB, fcomplex &RMIN, fcomplex &RPLUS, fcomplex &RQQ, fcomplex &RQG,
	     fcomplex &RGQ, fcomplex &RGG, fcomplex &C2Q, fcomplex &C2G, fcomplex &CDYQ, fcomplex &CDYG, fcomplex &xn, int &FR,
	     fcomplex &QQI, fcomplex &QGF, fcomplex &GQI, fcomplex &GGI, fcomplex &GGF, fcomplex &NS1MI, fcomplex &NS1PI, fcomplex &NS1F,
	     fcomplex &QQ1F, fcomplex &QG1F, fcomplex &GQ1I, fcomplex &GQ1F, fcomplex &GG1I, fcomplex &GG1F, fcomplex &C2QI, fcomplex &C2GF,
	     fcomplex &CDYQI, fcomplex &CDYGI);
  void h2calc_(fcomplex &C2qg, fcomplex &C2NSqqb, fcomplex &C2NSqq, fcomplex &C2Sqqb,fcomplex &xn);

  //access dyres PDF in N-space
  extern struct {
    fcomplex cfx1_[136][11];
    fcomplex cfx2p_[136][11];
    fcomplex cfx2m_[136][11];
  } creno_;
  
  //Catani-Seymour subtraction cut-offs for initial-initial, initial-final, final-initial, and final-final dipoles
  extern struct {
    double aii_;
    double aif_;
    double afi_;
    double aff_;
  } alfacut_;

  //initialization flag
  extern struct {
    int flag_;
  } flag_;

  extern struct {
    int qflag_;
    int gflag_;
  } flags_;

  extern struct {
    int noglue_;
    int ggonly_;
    int gqonly_;
  } noglue_;

  extern struct {
    int colourchoice_;
  } colc_;


  // z coupling
  extern struct {
    double l_[MAXNF];
    double r_[MAXNF];
    double q1_;
    double l1_;
    double r1_;
    double q2_;
    double l2_;
    double r2_;
    double le_;
    double ln_;
    double re_;
    double rn_;
    double sin2w_;
  } zcouple_;

  extern struct {
    int phot_;
  } dyphoton_;

  // masses
  extern struct {
    double md_;
    double mu_;
    double ms_;
    double mc_;
    double mb_;
    double mt_;
    double mel_;
    double mmu_;
    double mtau_;
    double hmass_;
    double hwidth_;
    double wmass_;
    double wwidth_;
    double zmass_;
    double zwidth_;
    double twidth_;
    double tauwidth_;
    double mtausq_;
    double mcsq_;
    double mbsq_;
  } dymasses_;

  // ewinput
  extern struct {
    double Gf_inp_;
    double aemmz_inp_;
    double xw_inp_;
    double wmass_inp_;
    double zmass_inp_;
  } ewinput_;

  //EW scheme
  extern struct {
    int ewscheme_;
  } ewscheme_;

  // ewcouple
  extern struct {
    double Gf_;
    double gw_;
    double xw_;
    double gwsq_;
    double esq_;
    double vevsq_;
  } ewcouple_;

  //alpha EM (MZ)
  /*  extern struct {
    double aemmz_;
    } em_;*/

  extern struct {
    int nf_;
  } nf_;

  extern struct {
    double Q_[2*MAXNF+1];
    double tau_[2*MAXNF+1];
  } ewcharge_;

  //branching ratio
  extern struct {
    double brnrat_;
  } brnrat_;

  //boson charge (used for the CKM matrix)
  extern struct {
    int nwz_;
  } nwz_;

  extern struct {
    int n2_;
    int n3_;
    double mass2_;
    double width2_;
    double mass3_;
    double width3_;
  } breit_;

  // CKM
  extern struct {
    double Vud_;
    double Vus_;
    double Vub_;
    double Vcd_;
    double Vcs_;
    double Vcb_;
  } cabib_;

  //QCD coupling
  extern struct {
    double gsq_;
    double as_;
    double ason2pi_;
    double ason4pi_;
  } qcdcouple_;

  extern struct {
    double b0_;
  } b0_;

  /*  // H+b mb msbar value
  extern struct {
    double mb_msbar_;
    } mb_msbar_;*/

  // dimensional regularization parameters
  extern struct {
    double epinv_;
  } epinv_;

  extern struct {
    double epinv2_;
  } epinv2_;

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
    double amz_;
  } couple_;

  /*extern struct {
    int nlooprun_;
    } nlooprun_;*/

  /*  extern struct {
    int lhapdfs_;
    } lhapdfs_;*/


  //input file related variables
  extern struct {
    double sroot_;
  } energy_;

  extern struct {
    int ih1_;
    int ih2_;
  } density_;

  extern struct {
    int nproc_;
  }  nproc_;

  extern struct {
    int dynamicscale_;
  } dynamicscale_;

  extern struct {
    double scale_;
    double musq_;
  } scale_;

  extern struct {
    double facscale_;
  } facscale_;

  //dynamic scale for each dipole
  extern struct {
    double dipscale_[41];
  } dipolescale_;

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
    char part_[4];
  } part_;

  extern struct {
    int zerowidth_;
  } zerowidth_;

  /*  extern struct {
    double Mwmin_;
    double Mwmax_;
    } mwminmax_;*/

  /*extern struct {
    int itmx1_;
    int ncall1_;
    int itmx2_;
    int ncall2_;
    } iterat_;*/

  /*  extern struct {
    int rseed_;
    } rseed_;*/

  /*  extern struct {
    int iset_;
    } pdfiset_;*/

  /*extern struct {
    double nset_;
    char prefix_[50];
    } prefix_;*/

  /*  extern struct {
    char PDFname_[30];
  } lhapdf_char_;

  extern struct {
    int PDFmember_;
    } lhapdf_int_;*/

  /*extern struct {
    char runstring_[30];
    } runstring_;*/

  /*  extern struct {
    int pr_;
    } pr_;*/

  /*  extern struct {
    double rtsmin_;
    } rtsmin_;*/

  extern struct {
    double wsqmin_;
    double wsqmax_;
  } limits_;

  extern struct {
    double taumin_;
    double logtaumin_;
  } taumin_;

  extern struct {
    double xmin_;
  } xmin_;

  extern struct {
    double p1ext_[4];
    double p2ext_[4];
  } pext_;

  //mcfm cut off on invariant mass pairs between emitted and radiator
  extern struct {
    double cutoff_;
  } cutoff_;

  extern struct {
    double xqtcut_;
  } qtcut_;

  extern struct {
    int doFill_;
  } dofill_;

  //efficiency variables (get rid of these)
  extern struct {
    int njetzero_;
    int ncutzero_;
    int ntotzero_;
    int ntotshot_;
  } efficiency_;

  void hists_setpdf_(int * npdf);
  void hists_fill_(double p3[4], double p4[4], double *weight);
  void hists_AiTest_(double pjet[4][12], double p4cm[4],double *Q,double *qt,double *y,double* pcosthCS, double* pphiCS, double *pphiVB, double *wt, double *loHst );
  void hists_real_dipole_(double p3[4], double p4[4], double *weight,int *nd);
  void hists_real_event_();
  void hists_fill_pdf_(double p3[4], double p4[4], double *weight, int *npdf);
  void hists_real_dipole_pdf_(double p3[4], double p4[4], double *weight,int *nd, int *npdf);
  void hists_real_event_pdf_(int* npdf);
}

#endif
