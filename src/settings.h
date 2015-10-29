#ifndef settings_h
#define settings_h

#include <vector>
using namespace std;

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
  void init();
  void readfromfile(const string fname);
  void initDyresSettings();

  // private:
  void dumpAll();
  void dumpI(string var, int    val );
  void dumpD(string var, double val );
  void dumpS(string var, string val );
  void dumpB(string var, bool   val );


  //original DYRES settings
  double sroot        ;
  int    ih1          ;
  int    ih2          ;
  int    nproc        ;
  double mur          ;
  double muf          ;
  double a_param      ;
  double g_param      ;
  int    order        ;
  string part         ;
  bool   zerowidth    ;
  double M_min        ;
  double M_max        ;
  int    itmx1        ;
  int    ncall1       ;
  int    itmx2        ;
  int    ncall2       ;
  int    rseed        ;
  int    PDFset       ;
  int    PDFmember    ;
  string LHAPDFset    ;
  int    LHAPDFmember ;
  string outputfile   ;
  int    itmxToFile   ;

  //resonance mass and width (used for breit wigner unweighting)
  double rmass, rwidth;

  // photon switch
  bool useGamma;

  //integration boundaries
  double ylow;
  double yhigh;
  double mlow;
  double mhigh;

  //dimension of integration
  int intDimRes;
  bool int2d, int3d, int4d;

  //type of integration for the counterterm
  int intDimCT;
  bool ctint2d, ctint3d, ctintvegas;

  //term switch
  bool doRES  , doVV , doCT   , doREAL , doVIRT , doLO   ;

  //resummation or fixed order switch
  bool fixedorder;

  //Cuba settings
  int cubaverbosity;
  int cubacores;
  int niterRES;
  int niterCT;
  int vegasncallsRES  ;
  int vegasncallsVV   ;
  int vegasncallsCT   ;
  int vegasncallsLO   ;
  int vegasncallsREAL ;
  int vegasncallsVIRT ;

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
  bool HackBinnerToFiller;

  // fiducial switches
  enum DetFiducial { GENEXP=0, D0=1, CDF=2, ATLAS=3, CMS7=4, CMS8=5};
  DetFiducial fiducial;
};

class binning
{
 public:
  binning() {};
  void init();
  void readfromfile(const string fname);
  // private:
  vector <double> qtbins;
  vector <double> ybins;
  vector <double> hist_qt_bins;
  vector <double> hist_y_bins;
};



extern settings opts;
extern binning bins;

extern bool cuts(double p3[4], double p4[4]);
bool decide_fiducial(double p3[4], double p4[4]);

#endif
