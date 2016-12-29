#ifndef vjint_h
#define vjint_h
#include <math.h>
extern "C"
{
  void qtdy_(double &res, double &err, double &chi2a, double &y1, double &y2, int &ord);
  double xdelta_(const double *x);
  double sing_(double z[2]);

  void flavour_(void);
  void utilities_(double & sh, double & th, double & uh, double & q2, double & ss2);
  void utilities2_(double & uh, double & q2);
  void utilities3_(double & q2);

  double cgg1_(double& sh, double& th, double& uh, double& q2);
  double cgg1x_(double& sh, double& th, double& uh, double& q2);
  double cqg3_(double& sh, double& th, double& uh, double& q2, int& flag);
  double cgq3_(double& sh, double& th, double& uh, double& q2, int& flag);
  double cqqb2_(double& sh, double& th, double& uh, double& q2, int& flag);
  double cqqb2x_(double& sh, double& th, double& uh, double& q2, int& flag);
  double daa_(double& sh, double& th, double& uh, double& q2, int& flag);
  double daax_(double& sh, double& th, double& uh, double& q2, int& flag);
  double dab_(double& sh, double& th, double& uh, double& q2);
  double dabx_(double& sh, double& th, double& uh, double& q2);
  double dbb_(double& sh, double& th, double& uh, double& q2);
  double dbbx_(double& sh, double& th, double& uh, double& q2);
  double dac_(double& sh, double& th, double& uh, double& q2);
  double dad_(double& sh, double& th, double& uh, double& q2);
  double dbc_(double& sh, double& th, double& uh, double& q2);
  double dbd_(double& sh, double& th, double& uh, double& q2);
  double dcc_(double& sh, double& th, double& uh, double& q2);
  double ddd_(double& sh, double& th, double& uh, double& q2);
  double dcdll_(double& sh, double& th, double& uh, double& q2);
  double dcdllx_(double& sh, double& th, double& uh, double& q2);
  double dcdlr_(double& sh, double& th, double& uh, double& q2);
  double dcdlrx_(double& sh, double& th, double& uh, double& q2);
  double eac_(double& sh, double& th, double& uh, double& q2);
  double ebd_(double& sh, double& th, double& uh, double& q2);
  double ead_(double& sh, double& th, double& uh, double& q2);
  double ebc_(double& sh, double& th, double& uh, double& q2);

  
  extern struct {
    double gevpb_;
  } gevpb_;

  extern struct {
    int flagch_;
  } flagch_;

  extern struct {
    int iter_;
    int locall_;
    int nlocall_;
  } vegint_;

  extern struct {
    int ih1_;
    int ih2_;
  } dypdf_;

  extern struct {
    int iord_;
  } vjorder_;

  extern struct {
    int prodflag_;
  } prodflag_;

  extern struct {
    double vud_;
    double vus_;
    double vub_;
    double vcd_;
    double vcs_;
    double vcb_;
    double vtd_;
    double vts_;
    double vtb_;
  } ckm_;

  //alpha EM (MZ)
  extern struct {
    double aemmz_;
  } em_;

  extern struct {
    double eq_[5];
    double alq_[5];
    double arq_[5];
    double ckm_[6][6];
    double delta_[5][5];
    double tau3_[5][5];
  } quarks_;

  extern struct {
    double pi_;
    double cf_;
    double ca_;
    double tr_;
    double xnc_;
  } const2_;

  extern struct {
    double xw_;
    double cw_;
    double sw_;
    double alpha0_;
  } dycouplings_;

  extern struct {
    double s_;
    double ss_;
  } dypara_;

  //variables
  /*
  extern struct {
    double amv_;
    double y1_;
    double y2_;
    double qtbis_;
    double gf_;
    double ppi_;
    double ssroot_;
    double sw2_;
    double aem_;
    double ic_;
  } cdyqt_;
  */

  extern struct {
    double yv_;
    double expyp_;
    double expym_;
  } yv_;

  extern struct {
    double tm_;
  } tm_;
  
  extern struct {
    double xmur_;
    double xmuf_;
    double xmur2_;
    double xmuf2_;
  } scales2_;

  extern struct {
    double as_;
  } asnew_;

  extern struct {
    double asp_;
  } asp_;
  
  extern struct {
    double siggamma_;
    double sigint_;
    double sigz_;
    double sigw_;
  } sigs_;

  extern struct {
    double qt_;
    double q_;
    double q2_;
  } internal_;

  extern struct {
    double x1_;
    double x2_;
  } fractions_;
  
}
#pragma omp threadprivate(yv_,tm_,scales2_,asnew_,asp_,sigs_,internal_,fractions_)

namespace vjint
{
  extern const double zmin;
  extern const double zmax;
  extern const double lz;

  extern int z1rule;
  extern int z2rule;
  extern double *t1;
  extern double *t2;

  extern double brz;
  extern double brw;
  
  extern void init();
  extern double vint(double m, double pt, double y);
  extern double calc(double m, double pt, double y);
  extern double sing();

}
#endif
