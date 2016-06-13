#ifndef kinematic_h
#define kinematic_h
namespace kinematic
{
  void set(double p3[4], double p4[4]);
  void calc_vb();
  void calc_angles();

  // kinematic and angular variables calculations
  extern double lp[4], lm[4];
  extern double v[4];
      
  extern double m,m2,qt,qt2,y;
  extern double phiZ;
  extern double costh,phi_lep;
  extern double mtrans;

  void calcM();
  void calcQt();
  void calcY();
  void calcPhiZ();
  void calcMt();

  void calcCosThCS();
  void calcPhiCS();

  void rotate(double vin[3], double phi, double ax[3], double vout[]);
  double Vplus(double p[4]);
  double Vminus(double p[4]);
}
#endif
