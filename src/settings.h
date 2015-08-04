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
        InputParser( string _filename = "config.ini", string _charset="#=[ ]", string _white=" \t");
        ~InputParser();
        // getters
        double GetNumber(string name);
        string GetString(string name);
        bool GetBool(string name);
        void GetVectorDouble(string name, vector<double> &vec);
    private :
        // functions
        void parse_file();
        void trim(string & str);
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
  void readfromfile();

  // private:

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
  void readfromfile();
  // private:
  vector <double> qtbins;
};



extern settings opts;
extern binning bins;

extern bool cuts(double p3[4], double p4[4]);

#endif
