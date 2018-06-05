#ifndef mcfm_interface_h
#define mcfm_interface_h

#include "parton.h"
#include "fcomplex.h"

extern "C"
{
  //interface to MCFM fortran functions and common blocks
  
  void breitw_(double& x1, double& mminsq, double& mmaxsq, double& rmass, double& rwidth, double& msq, double& wt);
  void boost_(double& mass, double p1[],double p_in[], double p_out[]);
  void branch_(double& brwen, double& brzee, double& brtau, double& brtop);
  void ckmfill_(int& nwz);
  void scaleset_(double& q2);
  void qqb_z_(double p[4][12], double msqc[11][11]);
  void qqb_w_(double p[4][12], double msqc[11][11]);
  void qqb_z_g_(double p[4][12], double msqc[11][11]);
  void qqb_w_g_(double p[4][12], double msqc[11][11]);
  void spinoru_(int &N, double p[4][12], fcomplex za[12][12], fcomplex zb[12][12]);
  
  //Catani-Seymour subtraction cut-offs for initial-initial, initial-final, final-initial, and final-final dipoles
  extern struct {
    double aii_;
    double aif_;
    double afi_;
    double aff_;
  } alfacut_;

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
  /*
  extern struct {
    double Gf_inp_;
    double aemmz_inp_;
    double xw_inp_;
    double wmass_inp_;
    double zmass_inp_;
  } ewinput_;
  */

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
  
  extern struct {
    double vsq_[11][11];
    double vsum_[11][11];
  } mcfmckm_;
  
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
}  
#pragma omp threadprivate(scale_,facscale_,qcdcouple_)

namespace mcfm
{
  void init();
  void set_mass_bounds();
}
#endif
