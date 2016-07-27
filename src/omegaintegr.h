#ifndef omegaintegr_h
#define omegaintegr_h

extern "C" {
  void cthmoments_(double &cthmom0, double &cthmom1, double &cthmom2);
  void genv4p_();
}

extern void genV4p(double m, double qt, double y, double phi);
//extern void genl4p(float costh, float phi_lep);
extern void genl4p(double costh, double phi_lep);
extern void getp3(double p[4]);
extern void getp4(double p[4]);

#endif
