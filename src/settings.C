#include <cmath>
#include <cstring>

#include "settings.h"
#include "interface.h"

settings opts;
binning bins;

void settings::init()
{
  
  //Z
  //resonance mass and width (used for breit wigner unweighting)
  rmass = 91.1876;
  rwidth = 2.495;
  //Set boundaries of integration
  ylow = 2;
  yhigh = 2.4;
  mlow = 66.;
  mhigh = 116.;

  /*
  //W
  //resonance mass and width (used for breit wigner unweighting)
  rmass = 80.385;
  rwidth = 2.091;
  //Set boundaries of integration
  ylow = -5.;
  yhigh = 5.;
  mlow = 10.;
  mhigh = 1000.;
  */
  
  //type of integration
  int2d = false;
  int3d = true;
  int4d = false;
  
  //Cuba settings
  cubaverbosity = 0;   //Cuba info messsages, from 0 to 3
  niter = 0;           //only for 2d and 3d cuhre integration
  vegasncalls = 10000; //only for 4d vegas integration
  
  //total or with lepton cuts
  makelepcuts = true;

  //integration types and settings for costh phi_lep phase space
  cubaint = true;    //integration with Cuba Suave
  trapezint = false; //trapezoidal rule for the phi_lep integration and semi-analytical for costh
  quadint  = false;  //quadrature rule for the phi_lep integration and semi-analytical for costh
  
  suavepoints = 1000000;  //number of points for suave integration, newpoints is set to suavepoints/10;
  nphitrape = 1000;       //number of steps for trapezoidal rule of phi_lep integration
  quadnphi = 20;          //number of segments for quadrature rule of phi_lep integration
  ncstart = 1000;         //starting sampling for the costh semi-analytical integration (common settings for the trapezoidal and quadrature rules)

  //qt-recoil prescriptions
  qtrec_naive = false;
  qtrec_cs = false;
  qtrec_kt0 = true;
  
  //debug settings
  timeprofile = false; //debug and time profile resummation integration
  verbose = false; //debug and time profile costh phi_lep integration

  //Set to 1 to use the dyres approximation of PDFs and integration contour in the complex plane for the Mellin inversion
  //Set to 0 to use exact PDFs and straight line contour in the complex plane for the Mellin inversion
  opts_.approxpdf_ = 0;
}

void settings::readfromfile(const char * fname){
    //read input settings from file
    InputParser in(fname);
    sroot          = in.GetNumber ( "sroot"          ); //7e3
    ih1            = in.GetNumber ( "ih1"            ); //1
    ih2            = in.GetNumber ( "ih2"            ); //1              # ih1,        ih2
    nproc          = in.GetNumber ( "nproc"          ); //3              # nproc
    mur            = in.GetNumber ( "mur"            ); //91.1876e0
    muf            = in.GetNumber ( "muf"            ); //91.1876e0      # mur,        muf
    a_param        = in.GetNumber ( "a_param"        ); //2.0e0          # a_param
    g_param        = in.GetNumber ( "g_param"        ); //1.0e0          # g_param
    order          = in.GetNumber ( "order"          ); //1              # order
    part           = in.GetString ( "part"           ); //virt           # part
    zerowidth      = in.GetBool   ( "zerowidth"      ); //false          # zerowidth
    M_min          = in.GetNumber ( "M_min"          ); //66d0
    M_max          = in.GetNumber ( "M_max"          ); //116d0          # M_min,      M_max
    itmx1          = in.GetNumber ( "itmx1"          ); //1
    ncall1         = in.GetNumber ( "ncall1"         ); //100000         # itmx1,      ncall1
    itmx2          = in.GetNumber ( "itmx2"          ); //1
    ncall2         = in.GetNumber ( "ncall2"         ); //100            # itmx2,      ncall2
    rseed          = in.GetNumber ( "rseed"          ); //123456         # rseed
    PDFset         = in.GetNumber ( "PDFset"         ); //92
    PDFmember      = in.GetNumber ( "PDFmember"      ); //0              # set,member  (native  PDFs)
    LHAPDFset      = in.GetString ( "LHAPDFset"      ); //CT10nlo.LHgrid
    LHAPDFmember   = in.GetNumber ( "LHAPDFmember"   ); //0              # set,        member   (LHAPDFs)
    outputfile     = in.GetString ( "outputfile"     ); //'LHC7-Z-nnlo'  # outputfile
    itmxToFile     = in.GetNumber ( "itmxToFile"     ); //0              # number      of       last       itmx    to          write           on          file

    rmass          = in.GetNumber ( "rmass"          ); //91.1876
    rwidth         = in.GetNumber ( "rwidth"         ); //2.495
    ylow           = in.GetNumber ( "ylow"           ); //2
    yhigh          = in.GetNumber ( "yhigh"          ); //2.4
    mlow           = in.GetNumber ( "mlow"           ); //66.
    mhigh          = in.GetNumber ( "mhigh"          ); //116.
    int2d          = in.GetBool   ( "int2d"          ); //false
    int3d          = in.GetBool   ( "int3d"          ); //true
    int4d          = in.GetBool   ( "int4d"          ); //false
    cubaverbosity  = in.GetNumber ( "cubaverbosity"  ); //0              # Cuba        info     messsages, from    0           to              3
    niter          = in.GetNumber ( "niter"          ); //0              # only        for      2d         and     3d          cuhre           integration
    vegasncalls    = in.GetNumber ( "vegasncalls"    ); //10000          # only        for      4d         vegas   integration
    makelepcuts    = in.GetBool   ( "makelepcuts"    ); //true
    cubaint        = in.GetBool   ( "cubaint"        ); //true           # integration with     Cuba       Suave
    trapezint      = in.GetBool   ( "trapezint"      ); //false          # trapezoidal rule     for        the     phi_lep     integration     and         semi-analytical for         costh
    quadint        = in.GetBool   ( "quadint"        ); //false          # quadrature  rule     for        the     phi_lep     integration     and         semi-analytical for         costh
    suavepoints    = in.GetNumber ( "suavepoints"    ); //1000000        # number      of       points     for     suave       integration,    newpoints   is              set         to    suavepoints/10;
    nphitrape      = in.GetNumber ( "nphitrape"      ); //1000           # number      of       steps      for     trapezoidal rule            of          phi_lep         integration
    quadnphi       = in.GetNumber ( "quadnphi"       ); //20             # number      of       segments   for     quadrature  rule            of          phi_lep         integration
    ncstart        = in.GetNumber ( "ncstart"        ); //1000           # starting    sampling for        the     costh       semi-analytical integration (common         settings    for   the             trapezoidal and quadrature rules)
    qtrec_naive    = in.GetBool   ( "qtrec_naive"    ); //false
    qtrec_cs       = in.GetBool   ( "qtrec_cs"       ); //false
    qtrec_kt0      = in.GetBool   ( "qtrec_kt0"      ); //true
    timeprofile    = in.GetBool   ( "timeprofile"    ); //false          # debug       and      time       profile resummation integration
    verbose        = in.GetBool   ( "verbose"        ); //false          # debug       and      time       profile costh       phi_lep         integration
    opts_.approxpdf_ = in.GetNumber ( "opts_approxpdf" ); //0


    return ;
}


void settings::initDyresSettings(){
    energy_      . sroot_     = sroot        ;         //7e3
    density_     . ih1_       = ih1          ;         //1
    density_     . ih2_       = ih2          ;         //1              # ih1,       ih2
    nproc_       . nproc_     = nproc        ;         //3              # nproc
    scale_       . scale_     = mur          ;         //91.1876e0
    facscale_    . facscale_  = muf          ;         //91.1876e0      # mur,       muf
    a_param_     . a_param_   = a_param      ;         //2.0e0          # a_param
    g_param_     . g_param_   = g_param      ;         //1.0e0          # g_param
    nnlo_       . order_     = order        ;         //1              # order
    zerowidth_   . zerowidth_ = zerowidth    ;         //false          # zerowidth
    mwminmax_    . Mwmin_     = M_min        ;         //66d0
    mwminmax_    . Mwmax_     = M_max        ;         //116d0          # M_min,     M_max
    iterat_      . itmx1_     = itmx1        ;         //1
    iterat_      . ncall1_    = ncall1       ;         //100000         # itmx1,     ncall1
    iterat_      . itmx2_     = itmx2        ;         //1
    iterat_      . ncall2_    = ncall2       ;         //100            # itmx2,     ncall2
    rseed_       . rseed_     = rseed        ;         //123456         # rseed
    pdfiset_     . iset_      = PDFset       ;         //92
    prefix_      . nset_      = PDFmember    ;         //0              # set,member (native PDFs)
    lhapdf_int_  . PDFmember_ = LHAPDFmember ;         //0              # set,       member  (LHAPDFs)
    pr_          . pr_        = itmxToFile   ;         //0              # number     of      last      itmx to write on file

    strncpy( part_        . part_      , part         .c_str(), 4); //virt           # part
    strncpy( lhapdf_char_ . PDFname_   , LHAPDFset    .c_str(), 30); //CT10nlo.LHgrid
    strncpy( runstring_   . runstring_ , outputfile   .c_str(), 30); //'LHC7-Z-nnlo'  # outputfile
}

bool cuts(double p3[4], double p4[4])
{
  if (!opts.makelepcuts)
    return true;
  double pt3 = sqrt((float)pow(p3[0],2)+pow(p3[1],2));
  if (pt3 < 20)
    return false;
  double pt4 = sqrt((float)pow(p4[0],2)+pow(p4[1],2));
  if (pt4 < 20)
    return false;
  double y3 = 0.5 *log((p3[3] + p3[2]) / (p3[3] - p3[2]));
  if (fabs(y3) > 2.4)
    return false;
  double y4 = 0.5 *log((p4[3] + p4[2]) / (p4[3] - p4[2]));
  if (fabs(y4) > 2.4)
    return false;

  //  if (y3-y4 < 0.2)
  //    return false;

  //  if (66 < m < 116
  //      && pt3 > 20 && pt4 > 20
  //      && fabs (y3) < 2.4 && fabs(y4) < 2.4)
  //    cut = true;
  return true;
}

void binning::readfromfile(const char * fname){
    InputParser in(fname);
    qtbins.clear();
    in.GetVectorDouble("qtbins",qtbins);
    return;
}

void binning::init()
{
  qtbins.push_back(0  ); 
  qtbins.push_back(2  ); 
  qtbins.push_back(4  ); 
  qtbins.push_back(6  ); 
  qtbins.push_back(8  ); 
  qtbins.push_back(10 ); 
  qtbins.push_back(12 ); 
  qtbins.push_back(14 ); 
  qtbins.push_back(16 ); 
  qtbins.push_back(18 ); 
  qtbins.push_back(22 ); 
  qtbins.push_back(26 ); 
  qtbins.push_back(30 ); 
  qtbins.push_back(34 ); 
  qtbins.push_back(38 ); 
  qtbins.push_back(42 ); 
  qtbins.push_back(46 ); 
  qtbins.push_back(50 ); 
  qtbins.push_back(54 ); 
  qtbins.push_back(60 ); 
  qtbins.push_back(70 ); 
  qtbins.push_back(80 ); 
  qtbins.push_back(100); 
  qtbins.push_back(150); 
  qtbins.push_back(200); 
  qtbins.push_back(300);
  qtbins.push_back(800);
  /*
  qtbins.push_back(0  );       qtbins.push_back(0.5  );  
  qtbins.push_back(1  );       qtbins.push_back(1.5  );  
  qtbins.push_back(2  );       qtbins.push_back(2.5  );  
  qtbins.push_back(3  );       qtbins.push_back(3.5  );  
  qtbins.push_back(4  );       qtbins.push_back(4.5  );  
  qtbins.push_back(5  );       qtbins.push_back(5.5  );  
  qtbins.push_back(6  );       qtbins.push_back(6.5  );  
  qtbins.push_back(7  );       qtbins.push_back(7.5  );  
  qtbins.push_back(8  );       qtbins.push_back(8.5  );  
  qtbins.push_back(9  );       qtbins.push_back(9.5  );  
  qtbins.push_back(10 );       qtbins.push_back(10.5 ); 
  qtbins.push_back(11 );       qtbins.push_back(11.5 ); 
  qtbins.push_back(12 );       qtbins.push_back(12.5 ); 
  qtbins.push_back(13 );       qtbins.push_back(13.5 ); 
  qtbins.push_back(14 );       qtbins.push_back(14.5 ); 
  qtbins.push_back(15 );       qtbins.push_back(15.5 ); 
  qtbins.push_back(16 );       qtbins.push_back(16.5 ); 
  qtbins.push_back(17 );       qtbins.push_back(17.5 ); 
  qtbins.push_back(18 );       qtbins.push_back(18.5 ); 
  qtbins.push_back(19 );       qtbins.push_back(19.5 ); 
  qtbins.push_back(20 );       qtbins.push_back(20.5 ); 
  qtbins.push_back(21 );       qtbins.push_back(21.5 ); 
  qtbins.push_back(22 );       qtbins.push_back(22.5 ); 
  qtbins.push_back(23 );       qtbins.push_back(23.5 ); 
  qtbins.push_back(24 );       qtbins.push_back(24.5 ); 
  qtbins.push_back(25 );       qtbins.push_back(25.5 ); 
  qtbins.push_back(26 );       qtbins.push_back(26.5 ); 
  qtbins.push_back(27 );       qtbins.push_back(27.5 ); 
  qtbins.push_back(28 );       qtbins.push_back(28.5 ); 
  qtbins.push_back(29 );       qtbins.push_back(29.5 ); 
  qtbins.push_back(30 );       qtbins.push_back(30.5 ); 
  qtbins.push_back(31 );       qtbins.push_back(31.5 ); 
  qtbins.push_back(32 );       qtbins.push_back(32.5 ); 
  qtbins.push_back(33 );       qtbins.push_back(33.5 ); 
  qtbins.push_back(34 );       qtbins.push_back(34.5 ); 
  qtbins.push_back(35 );       qtbins.push_back(35.5 ); 
  qtbins.push_back(36 );       qtbins.push_back(36.5 ); 
  qtbins.push_back(37 );       qtbins.push_back(37.5 ); 
  qtbins.push_back(38 );       qtbins.push_back(38.5 ); 
  qtbins.push_back(39 );       qtbins.push_back(39.5 ); 
  qtbins.push_back(40 );       qtbins.push_back(40.5 ); 
  qtbins.push_back(41 );       qtbins.push_back(41.5 ); 
  qtbins.push_back(42 );       qtbins.push_back(42.5 ); 
  qtbins.push_back(43 );       qtbins.push_back(43.5 ); 
  qtbins.push_back(44 );       qtbins.push_back(44.5 ); 
  qtbins.push_back(45 );       qtbins.push_back(45.5 ); 
  qtbins.push_back(46 );       qtbins.push_back(46.5 ); 
  qtbins.push_back(47 );       qtbins.push_back(47.5 ); 
  qtbins.push_back(48 );       qtbins.push_back(48.5 ); 
  qtbins.push_back(49 );       qtbins.push_back(49.5 ); 
  qtbins.push_back(50 );       qtbins.push_back(50.5 ); 
  */
}


// InputParser definitions
//

#include <stdexcept>

InputParser::InputParser( string _filename, string _charset, string _white):
    filename ( _filename   ),
    Ccommnt  ( _charset[0] ),
    Cassign  ( _charset[1] ),
    CopenAr  ( _charset[2] ),
    CdeliAr  ( _charset[3] ),
    CclosAr  ( _charset[4] ),
    Swhite   ( _white      )
{
    parse_file();
}


InputParser::~InputParser(){
}

double InputParser::GetNumber(string name){
    has_key(name);
    string val = data[name];
    return stod(val);
}

string InputParser::GetString(string name){
    has_key(name);
    string val = data[name];
    trim(val);
    return val;
}


bool InputParser::GetBool(string name){
    has_key(name);
    string val = data[name];
    // lower case
    std::transform(val.begin(), val.end(), val.begin(), ::tolower);
    if (val.compare(0,4,"true") == 0 ) return true;
    return false;
}


void InputParser::GetVectorDouble(string name, vector<double> &vec){
    has_key(name);
    string val = data[name];
    // has open/close array
    size_t strBegin = val.find(CopenAr,0);
    size_t strEnd   = val.find(CclosAr,0);
    // we need both open and close otherwise is something wrong
    if (strBegin==string::npos || strEnd==string::npos) throw invalid_argument("Missing open/close character.");
    // parse what is between them
    val = val.substr(strBegin+1,strEnd-strBegin-1);
    while (!val.empty()){
        size_t pos = 0;
        // find delim
        pos = val.find(CdeliAr);
        // get substring
        string tmp = val.substr(0,pos);
        trim(tmp);
        // expecting number otherwise exception
        if (!tmp.empty()) vec.push_back(stod(tmp));
        val = val.substr(pos+1);
    }
    return;
}


void InputParser::parse_file(){
    ifstream fstrm(filename.c_str());
    string line;
    // line by line
    while(getline(fstrm,line))   {
        size_t pos = 0;
        // check for comments
        pos = line.find(Ccommnt,pos);
        if (pos != string::npos) line = line.substr(0,pos);
        // split string by delim
        pos = line.find(Cassign,0);
        if (pos == string::npos || line.empty()) continue;
        string key  = line.substr(0,pos);
        string val = line.substr(pos+1);
        // remove leading / trailing whitespace
        trim(key);
        trim(val);
        /// @todo: Multiline setting. If has CopenAr load until CclosAr
        // save to map
        data[key]=val;
    }
    fstrm.close();
    return;
}


void InputParser::trim(string & str){
    size_t strBegin = str.find_first_not_of(Swhite);
    if (strBegin == std::string::npos){
        str="";
        return ; // no content
    }
    size_t strEnd = str.find_last_not_of(Swhite);
    size_t strRange = strEnd - strBegin + 1;
    str = str.substr(strBegin, strRange);
    return;
}


void InputParser::has_key(const string key){
    if ( data.count(key) == 0 ){
        string msg = "No setting with name '";
        msg += key;
        msg += "'";
        throw invalid_argument(msg.c_str());
    }
    return;
}




