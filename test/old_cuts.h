#ifndef cuts_h
#define cuts_h

namespace cuts {
  extern double getPt(double p[4]);
  extern double getY(double p[4]);
  extern double getEta(double p[4]);
  extern double getMt(double p3[4], double p4[4]);
  extern double getM(double p3[4], double p4[4]);
  extern double getEtMiss(double p3[4], double p4[4]);
  extern double getLPt(double p3[4], double p4[4]);
  extern double getLY(double p3[4], double p4[4]);

  extern bool fiducial_D0(double p3[4], double p4[4]);
  extern bool fiducial_CDF(double p3[4], double p4[4]);
  extern bool fiducial_ATLAS(double p3[4], double p4[4]);
  extern bool fiducial_CMS7(double p3[4], double p4[4]);
  extern bool fiducial_CMS8(double p3[4], double p4[4]);
  extern bool user_cuts(double p3[4], double p4[4]);
  extern bool decide_fiducial(double p3[4], double p4[4]);

  extern bool lep(double p3[4], double p4[4]);
}

#endif
