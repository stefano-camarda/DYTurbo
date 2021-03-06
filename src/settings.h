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
  void check_consistency();
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

  //resummation or fixed order switch
  bool fixedorder;

  //order
  int    order;       //Main order
  int    order_sudak; //Order of the Sudakov
  int    order_hcoef; //Order of the H coefficients
  int    order_evol;  //Order of the PDF evolution
  int    order_expc;   //Order of the C exponentiation
  //  int    order_ford;   //Order of the finite order part
  //  int    order_alphas; //Order of the alphas running

  bool qbox;
  
  //Non-perturbative form factor
  int npff;
  double g1,g2,g3; //Gaussian
  double e;        //Exponential
  double g0;       //Collins-Rogers
  double Q0;       //reference mass
  double a2,a2p;   //Dokshitzer, Marchesini, Webber
  
  //Flavour dependent g1
  bool flavour_kt;
  double g1_uv = 0.5;
  double g1_us = 0.5;
  double g1_dv = 0.5;
  double g1_ds = 0.5;
  double g1_ss = 0.5;
  double g1_ch = 0.5;
  double g1_bo = 0.5;
  double g1_gl = 0.5;
  
  //PDF settings
  string LHAPDFset    ;
  int    LHAPDFmember ;
  bool   externalpdf;

  //alphas running
  bool alphaslha;
  
  //functional forms of the QCD scales
  int fmures;
  int fmuren;
  int fmufac;
  
  //scale factors for the QCD scales
  double kmures;
  double kmuren;
  double kmufac;

  //scale factors for the matching scales
  double kmuc;
  double kmub;
  double kmut;
  
  //  double a_param;
  
  //EW parameters
  int ewscheme;
  double Gf, zmass, wmass;
  double xw, aemmz;
  double zwidth, wwidth;

  //Fixed width to running width translation
  bool runningwidth;
  
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
  bool damp;
  double dampk, dampdelta;
  int dampmode;

  //Resummation cutoff
  double qtcutoff;

  //Modified logarithms
  bool modlog;
  int p;
  
  //qtcut
  double xqtcut, qtcut;

  //integration settings
  int rseed;

  //dimension of integration for the resummed part
  int intDimRes;
  bool resint1d, resint2d, resint3d, resintvegas;

  //dimension of integration for the born configuration
  int intDimBorn;
  bool bornint1d, bornint2d, bornintvegas4d, bornintvegas6d;
  
  //type of integration for the counterterm
  int intDimCT;
  bool ctint1d, ctint2d, ctint3d, ctintvegas6d, ctintvegas8d;

  //type of integration for the V+j at LO
  int intDimVJ;
  bool vjint3d, vjint5d, vjintvegas7d;

  //type of integration, automatic selector
  bool BORNquad, CTquad, VJquad;
  
  //term switches
  bool doBORN;
  bool doCT;
  bool doVJ;
  
  bool doVJREAL, doVJVIRT;

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
  long long int vegasncallsVJREAL ;
  int vegasncallsVJVIRT ;
  bool vegascollect;
  
  //cubature settings
  bool pcubature;
  double relaccuracy;
  double absaccuracy;

  //costh boundaries
  double costhmin, costhmax;
  
  //lepton fiducial cuts
  bool makecuts;
  double lptcut, lycut; //charged leptons
  double mtcut, etmisscut;
  double lepptcut, lepycut, alpptcut, alpycut; //lepton and antilepton
  double lcptcut, lcymin, lcymax, lfptcut, lfymin, lfymax; //lc and lf are absolute-rapidity-ordered leptons
  double cthCSmin, cthCSmax;

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

  //numerical integration of the Sudakov
  bool numsud;
  //use a numerical solution for the running of alphas
  bool asrgkt;

  //numerical integration of the C exponentiation
  bool numexpc;
  
  //b-space prescription
  int bprescription;
  
  //blim parameter of the bstar prescription (acts as an IR cut-off)
  double blim;
  double blim_pdf, blim_sudakov, blim_expc;

  //force bstar prescription
  bool bstar_pdf, bstar_sudakov, bstar_expc;
  
  //arg(z) in the complex plane for the minimal prescription
  double phibr;
  //select the point bc, where the integration contour is bended in the complex plane, as a fraction of b_L = ... (Landau singularity)
  double bcf;
  
  //settings for the Mellin transform
  int mellintr;

  //settings for the Mellin inversion
  int mellininv;
  int mellinintervals;
  int mellinrule;
  double zmax;
  int ncycle;
  double cpoint;
  double cshift;
  double phi;
  int mellincores;
  bool mellin1d;
  int melup;
  bool xspace;

  //settings for x-to-N Mellin transform
  int pdfrule;
  
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
  bool gridverbose;

  //output settings
  bool texttable;
  bool redirect;
  bool unicode;
  bool silent;
  bool makehistos;
  string output_filename; // Output Filenames
  
  //resummed code in C++
  bool resumcpp = true; //use C++ code for resummation     

  //counter term code in C++
  bool ctcpp = true;   //use C++ code for the counter term
  
  //dyres, pegasus,lhapdf PDF evolution
  int evolmode;

  //Evolve PDFs from mufac instead of mures
  //bool mufevol;

  //Do not factorize muf variations with expc
  bool mufvar;
  
  //switches for C exponentiation
  int expc, ntaylor;

  //sum all logs before exponentiation
  bool sumlogs;
  
  //bin width normalisation
  bool ptbinwidth, ybinwidth, mbinwidth;

  // Force to loop over all bins even you have all Vegas integrands
  bool force_binsampling = false;

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
