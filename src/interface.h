#ifndef interface_h
#define interface_h

extern "C" struct {
  double a_param_;
  double b0p_;
} a_param_;

extern "C" struct {
  double g_param_;
} g_param_;


extern "C" {
  double resumm_(double &costh, double &mm, double &qtt, double &yy, int& mode);
  void setup_();
  void dyinit_();
  bool cuts_(double p[4][12], int &njet);

  void breitw_(double& x1, double& mminsq, double& mmaxsq, double& rmass, double& rwidth, double& msq, double& wt);
  void boost_(double& mass, double p1[],double p_in[], double p_out[]);

  
  void rapintegrals_(double &ymin,double &ymax, double& mass, int& nocuts);

  extern struct {
    int approxpdf_;
  } opts_;

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
