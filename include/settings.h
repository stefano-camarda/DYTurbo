#ifndef settings_h
#define settings_h

#include <vector>
using namespace std;

class settings
{
public:
  settings() {};
  void init();

  // private:

  //resonance mass and width (used for breit wigner unweighting)
  double rmass, rwidth;

  //integration boundaries
  double ylow;
  double yhigh;
  double mlow;
  double mhigh;

  //type of integration
  bool int2d, int3d, int4d;
  
  //Cuba settings
  int cubaverbosity;
  int niter;
  int vegasncalls;

  bool makelepcuts;

  //integration types and settings for costh phi_lep phase space
  bool cubaint;
  int suavepoints;

  bool trapezint;
  int nphitrape;
  int ncstart;

  bool quadint;
  int quadnphi;

  //qt-recoil prescriptions
  bool qtrec_naive, qtrec_cs, qtrec_kt0;
  
  //debug settings
  bool timeprofile;
  bool verbose;
};

class binning
{
 public:
  binning() {};
  void init();
  // private:
  vector <double> qtbins;
 
};
extern settings opts;
extern binning bins;

extern bool cuts(double p3[4], double p4[4]);

#endif
