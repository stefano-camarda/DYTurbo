#ifndef dyfin_h
#define dyfin_h
extern "C" {
  //  double resumm_(double &costh, double &mm, double &qtt, double &yy);
  double realvirt2_(double r[22], double &wgt);
  double countint_(double r[22], double &wgt);
  double realint_(double r[22], double &wgt);
  double virtint_(double r[22], double &wgt);
  void setup_();
  void dyinit_();
  
  double countterm_(double &costh, double &mm, double &qtt, double &yy, double &alfa, double &beta, double &cthmom0, double &cthmom1, double &cthmom2);

  extern struct {
    double amz_;
  } couple_;
  
  extern struct {
    int nlooprun_;
  } nlooprun_;

  extern struct {
    double facscale_;
  } facscale_;

  extern struct {
    double scale_;
    double musq_;
  } scale_;

  extern struct {
    int lhapdfs_;
  } lhapdfs_;


}
#endif
