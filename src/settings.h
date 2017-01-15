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
        void parse_file(const string fname);
    private :
        // functions
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

// Simple program quiting exception
struct QuitProgram : public std::runtime_error { QuitProgram(string msg) : std::runtime_error (msg) {}; };

// Forward declaration
namespace cxxopts{
    class Options;
}
namespace po=cxxopts; // inspired by po = boost::program_options

class settings
{
public:
  settings() {};
  void parse_options(int argc, char * argv[]);
  void readfromfile(const string fname);
  void check_consitency();
  void parse_binning(string name, vector<double> &vec, po::Options &args);

  // private:
  void dumpAll();
  void dumpI(string var, int    val );
  void dumpD(string var, double val );
  void dumpS(string var, string val );
  void dumpB(string var, bool   val );

  // string helpers
  void ToLower(string &val){std::transform(val.begin(), val.end(), val.begin(), ::tolower);}
  void ToUpper(string &val){std::transform(val.begin(), val.end(), val.begin(), ::toupper);}
  vector<string> Tokenize(string val, char Delim=',');
  // is number: http://stackoverflow.com/a/16575025
  bool IsNumber(const string &s);

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

  //scale factors for the matching scales
  double kmuc;
  double kmub;
  double kmut;
  
  //Additional resummation scales
  double C1,C3;
  
  //  double a_param;
  
  //IR cut-off
  double blim;
  
  //EW parameters
  int ewscheme;
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
  //double ylow;
  //double yhigh;
  //double mlow;
  //double mhigh;

  //Resummation damping
  double dampk, dampdelta;
  int dampmode;

  //Resummation cutoff
  double qtcutoff;

  //qtcut
  double xqtcut, qtcut;

  //integration settings
  int rseed;

  //dimension of integration for the resummed part
  int intDimRes;
  bool resint2d, resint3d, resintvegas;

  //dimension of integration for the born configuration
  int intDimBorn;
  bool bornint2d, bornintvegas4d, bornintvegas6d;
  
  //type of integration for the counterterm
  int intDimCT;
  bool ctint2d, ctint3d, ctintvegas6d, ctintvegas8d;

  //  //type of integration for the V+j at LO
  int intDimVJ;
  bool vjint3d, vjint5d, vjintvegas7d;
  
  //term switches
  bool doBORN;
  bool doCT;
  bool doVJ;
  
  bool doVJREAL, doVJVIRT;

  //resummation or fixed order switch
  bool fixedorder;

  //Cuba settings
  int cubaverbosity;
  int cubacores;
  int cubanbatch;
  int niterBORN;
  int niterCT;
  int niterVJ;
  int vegasncallsBORN  ;
  int vegasncallsCT   ;
  int vegasncallsVJLO   ;
  int vegasncallsVJREAL ;
  int vegasncallsVJVIRT ;

  //cubature settings
  bool pcubature;
  double pcubaccuracy;

  //costh boundaries
  double costhmin, costhmax;
  
  //lepton fiducial cuts
  bool makecuts;
  double lptcut, lycut;
  double mtcut, etmisscut;
  double l1ptcut, l1ycut, l2ptcut, l2ycut;

  //integration types and settings for costh phi_lep phase space
  bool cubaint;
  int suavepoints;

  bool trapezint;
  int nphitrape;
  int ncstart;

  //quadrature rule in phi_lep
  bool quadint;
  int phiintervals;
  int phirule;

  //settings for Bessel integration
  double bintaccuracy;

  //settings for Mellin integration
  int mellinintervals;
  int mellinrule;
  double zmax;
  double cpoint;
  int mellincores;
  bool mellin1d;
  
  //settings for rapidity integration in 2D resummed piece
  int yintervals;
  int yrule;

  //settings for qt integration in 2D counter term
  int qtintervals;
  int qtrule;

  //settings for alfa beta scaled-PDF integration in counter term and born fixed order
  int abintervals;
  int abrule;

  //settings for the phi integration in the V+J 5d LO term when makecuts is false
  int vjphirule;

  //settings for the z1, z2 integration in the V+J 3d NLO singular term
  int zrule;

  //settings for the x integration in the V+J 3d delta term
  int xrule;
  
  //qt-recoil prescriptions
  bool qtrec_naive, qtrec_cs, qtrec_kt0;

  // PDF errors
  bool PDFerrors;
  int totpdf;
  
  //debug settings
  bool timeprofile;
  bool verbose;
  bool texttable;

  //resummed code in C++
  bool resumcpp;

  //counter term code in C++
  bool ctcpp;
  
  //dyres or pegasus PDF evolution
  int evolmode;

  //bin width normalisation
  bool ptbinwidth, ybinwidth;

  // Force to loop over all bins even you have all Vegas integrands
  bool force_binner_mode = false;

  // Calculate helicity cross sections
  int helicity;
};

class binning
{
 public:
  binning() {};
  void readfromfile(const string fname);
  void GetBins(string name,vector<double> &bins);
  // private:
  vector <double> qtbins;
  vector <double> ybins;
  vector <double> mbins;
  vector <double> hist_qt_bins;
  vector <double> hist_y_bins;
  vector <double> hist_m_bins;
 private:
  InputParser in;
};


extern settings opts;
extern binning bins;

#endif
