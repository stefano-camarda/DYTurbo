#ifndef vjint_h
#define vjint_h

extern "C"
{
  void qtdy_(double &res, double &err, double &chi2a, double &y1, double &y2, int &ord);
  double xdelta_(const double *x);
  double sing_(double z[2]);

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

}
#pragma omp threadprivate(yv_,tm_,scales2_,asnew_,asp_,sigs_,internal_)

namespace vjint
{
  extern double brz;
  extern double brw;
  
  extern void init();
  extern double vint(double m, double pt, double y);
  extern double calc(double m, double pt, double y);

}
#endif
