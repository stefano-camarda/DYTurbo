#ifndef settings_h
#define settings_h

#include "cuts.h"

#include <vector>
using namespace std;
using namespace cuts;

#include <algorithm>
#include <map>
#include <string>
#include <fstream>

class InputParser {
    public:
        // constructor
        InputParser( string _filename = "", string _charset="#=[ ]", string _white=" \t");
        ~InputParser();
        // getters
        double GetNumber(string name);
        string GetString(string name);
        bool GetBool(string name);
        void GetVectorDouble(string name, vector<double> &vec);
    private :
        // functions
        void parse_file(const string fname);
        void trim(string & str);
        void has_key(const string key);
        // data members
        map<string,string> data; ///< Data table.
        string filename; ///< Input file.
        // parser characters
        char Ccommnt; 
        char Cassign;
        char CopenAr;
        char CclosAr;
        char CdeliAr;
        string Swhite;
};

class settings
{
public:
  settings() {};
  void readfromfile(const string fname);

  // private:
  void dumpAll();
  void dumpI(string var, int    val );
  void dumpD(string var, double val );
  void dumpS(string var, string val );
  void dumpB(string var, bool   val );

  //process settings
  double sroot;
  int    ih1;
  int    ih2;
  int    nproc;
  double g_param;
  int    order;

  //PDF settings
  string LHAPDFset    ;
  int    LHAPDFmember ;

  //fixed or dynamic QCD scales
  bool dynamicscale, dynamicresscale;
  
  //scale factors for the QCD scales
  double kmures;
  double kmuren;
  double kmufac;

  //  double a_param;
  
  //IR cut-off
  double blim;
  
  //EW parameters
  double Gf, zmass, wmass;
  double xw, aemmz;
  double zwidth, wwidth;
  bool zerowidth;

  //CKM matrix
  double Vud, Vus, Vub;
  double Vcd, Vcs, Vcb;

  //Z/gamma* coupling
  double Zuu, Zdd, Zss, Zcc, Zbb;

  //resonance mass and width (used for breit wigner unweighting)
  double rmass, rwidth;

  // photon switch
  bool useGamma;

  //integration boundaries
  double ylow;
  double yhigh;
  double mlow;
  double mhigh;

  //Resummation damping
  double dampk, dampdelta;
  int dampmode;

  //Resummation cutoff
  double qtcutoff;

  //qtcut
  double xqtcut, qtcut;

  //integration settings
  int rseed;

  //dimension of integration
  int intDimRes;
  bool resint2d, resint3d, resintvegas;

  //type of integration for the counterterm
  int intDimCT;
  bool ctint2d, ctint3d, ctintvegas6d, ctintvegas8d;

  //  //type of integration for the V+j at LO
  //  int intDimVJ;
  //  bool vjint3d, vjintvegas7d;
  
  //term switches
  bool doRES, doVV, doCT, doREAL, doVIRT, doLO, doVJ;

  //resummation or fixed order switch
  bool fixedorder;

  //Cuba settings
  int cubaverbosity;
  int cubacores;
  int niterRES;
  int niterCT;
  int niterVJ;
  int vegasncallsRES  ;
  int vegasncallsVV   ;
  int vegasncallsCT   ;
  int vegasncallsLO   ;
  int vegasncallsREAL ;
  int vegasncallsVIRT ;

  //lepton fiducial cuts
  bool makelepcuts;
  double lptcut, lycut;
  double l1ptcut, l1ycut, l2ptcut, l2ycut;

  //integration types and settings for costh phi_lep phase space
  bool cubaint;
  int suavepoints;

  bool trapezint;
  int nphitrape;
  int ncstart;

  bool quadint;
  int quadnphi;

  //settings for Mellin integration
  int mellinintervals;
  int mellinrule;
  double zmax;
  
  //settings for rapidity integration
  int yintervals;
  int yrule;
  
  //qt-recoil prescriptions
  bool qtrec_naive, qtrec_cs, qtrec_kt0;

  // PDF errors
  bool PDFerrors;
  int totpdf;
  
  //debug settings
  bool timeprofile;
  bool verbose;

  bool resumcpp;

  int evolmode;

  // fiducial switches
  cuts::DetFiducial fiducial;
};

class binning
{
 public:
  binning() {};
  void readfromfile(const string fname);
  // private:
  vector <double> qtbins;
  vector <double> ybins;
  vector <double> hist_qt_bins;
  vector <double> hist_y_bins;
};


extern settings opts;
extern binning bins;

#endif
