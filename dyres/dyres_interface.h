#ifndef dyres_interface_h
#define dyres_interface_h

extern "C"
{
  //interface to DYRES fortran functions and common blocks
  double  pqq_(double& z);
  double  dypqg_(double& z);
  double  cqq_(double& z);
  double  cqg_(double& z);
  double  pqqint_(double& z);
  double  d0int_(double& z);
  double  d1int_(double& z);
  double  pqqqq_(double& z);
  double  pqqqg_(double& z); 
  double  pqggq_(double& z);
  double  pqggg_(double& z);
  double  cqqpqq_(double& z);
  double  cqqpqg_(double& z);
  double  cqgpgq_(double& z); 
  double  cqgpgg_(double& z);
  double  p2qqv_(double& x);
  double  p2qqbv_(double& x);
  double  p2qg_(double& x);
  double  p2qqs_(double& x);
  double  s2_(double& x);

  double itilde_(int& m);

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

  extern struct {
    double xmio_;
  } xmio_;

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
    double sigmaij_[11][11];
  } sigmaij_;


}
#pragma omp threadprivate(a_param_,sigmaij_,xmio_)

namespace dyres
{
  void init();
}
#endif
