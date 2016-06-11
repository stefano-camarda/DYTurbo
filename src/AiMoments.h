#ifndef AIMOMENTS_H
#define AIMOMENTS_H
#include "config.h"
#ifdef USEROOT

#include <vector>
#include <sstream>
#include <iostream>

#include "TH2D.h"
#include "TH3D.h"
#include "TH1D.h"
#include "TFile.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TLorentzVector.h"

#include "TLVUtils.h"

using namespace std;

class AiMoments {

 public:

  AiMoments(string AlgoName="AiMoments");
  ~AiMoments();

  bool Initialize(); 
  bool Finalize();
  bool Execute(double pt1, double eta1, double phi1, double m1, int type1,
	       double pt2, double eta2, double phi2, double m2, int type2,
	       double ebeam, double weight=1);
  inline const char * GetFileName() {return (m_algoname+".root").c_str(); }

  double pt_low;
  double pt_hgh;
  
 private:

  vector<double> m_Corr;
  
  vector<TProfile*>   aimomTruthCS_vs_pt, aimomTruthCS_vs_y, aimomTruthZD_vs_pt, aimomTruthZD_vs_y;  
  vector<TProfile2D*> aimomTruthCS_vs_pty, aimomTruthZD_vs_pty;  
  TH2D *m_yvspt, *m_yvsptS, *m_cosCSvsy, *m_cosHvsy, *m_phiHvsy, *m_phiCSvsy, *m_phiHvspt, *m_phiCSvspt, *m_cosHvspt, *m_cosCSvspt, *m_phiHvscosH, *m_phiCSvscosCS;
  TH3D *m_cosCSvsyvspt, *m_phiCSvsyvspt;
  TH1D *m_pt;

  TH1D *m_ptfid_inc, *m_ptfid_yb1, *m_ptfid_yb2, *m_ptfid_yb3;

  TH3D *m_yvsptvsm;

  TFile* fOut;
  string m_algoname;
  float A0shift;

};

#endif // USEROOT
#endif
