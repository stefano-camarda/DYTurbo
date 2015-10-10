#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>

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

  //resummation or fixed order switch
  fixedorder = false;
  
  //term switch
  doRES  = false;
  doVV   = false;
  doCT   = false;
  doREAL = false;
  doVIRT = false;
  doLO   = false;
  
  //Cuba settings
  cubaverbosity = 0;   //Cuba info messsages, from 0 to 3
  cubacores = 0;   //parallelization (0 = turn off)
  niterRES = 0;           //only for 2d and 3d cuhre integration
  vegasncallsRES  = 10000; // only for res 4d vegas integration
  vegasncallsVV   = 10000; // only for double virtual 6d vegas integration
  vegasncallsCT   = 10000; // only for lo 8d vegas integration
  vegasncallsLO   = 10000; // only for lo 7d vegas integration
  vegasncallsREAL = 10000; // only for real 10d vegas integration
  vegasncallsVIRT = 10000; // only for virt 8d vegas integration
  
  //total or with lepton cuts
  makelepcuts = true;
  fiducial = static_cast<settings::DetFiducial>(0);

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
  HackBinnerToFiller = false; //debug distributions

  //Set to 1 to use the dyres approximation of PDFs and integration contour in the complex plane for the Mellin inversion
  //Set to 0 to use exact PDFs and straight line contour in the complex plane for the Mellin inversion
  opts_.approxpdf_ = 0;
}

void settings::readfromfile(const string fname){
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

    rmass              = in.GetNumber ( "rmass"           ); //91.1876
    rwidth             = in.GetNumber ( "rwidth"          ); //2.495
    ylow               = in.GetNumber ( "ylow"            ); //2
    yhigh              = in.GetNumber ( "yhigh"           ); //2.4
    mlow               = in.GetNumber ( "mlow"            ); //66.
    mhigh              = in.GetNumber ( "mhigh"           ); //116.
    intDimRes          = in.GetNumber ( "intDimRes"       ); //3
    //int2d              = in.GetBool   ( "int2d"           ); //false
    //int3d              = in.GetBool   ( "int3d"           ); //true
    //int4d              = in.GetBool   ( "int4d"           ); //false
    intDimCT           = in.GetNumber ( "intDimCT"       ); //3
    //ctint2d            = in.GetBool   ( "ctint2d"         ); //true
    //ctint3d            = in.GetBool   ( "ctint3d"         ); //false
    //ctintvegas         = in.GetBool   ( "ctintvegas"      ); //false
    fixedorder         = in.GetBool   ( "fixedorder"      ); //false
    doRES              = in.GetBool   ( "doRES"           ); //false
    doVV               = in.GetBool   ( "doVV"            ); //false
    doCT               = in.GetBool   ( "doCT"            ); //false
    doREAL             = in.GetBool   ( "doREAL"          ); //false
    doVIRT             = in.GetBool   ( "doVIRT"          ); //false
    doLO               = in.GetBool   ( "doLO"            ); //false
    cubaverbosity      = in.GetNumber ( "cubaverbosity"   ); //0       # Cuba        info     messsages, from    0           to              3
    cubacores          = in.GetNumber ( "cubacores"       ); //0
    niterRES           = in.GetNumber ( "niterRES"        ); //0       # only        for      2d         and     3d          cuhre           integration
    niterCT            = in.GetNumber ( "niterCT"         ); //0       # only        for      3d          cuhre           integration
    vegasncallsRES     = in.GetNumber ( "vegasncallsRES"  ); //10000
    vegasncallsVV      = in.GetNumber ( "vegasncallsVV"   ); //10000
    vegasncallsCT      = in.GetNumber ( "vegasncallsCT"   ); //10000
    vegasncallsLO      = in.GetNumber ( "vegasncallsLO"   ); //10000
    vegasncallsREAL    = in.GetNumber ( "vegasncallsREAL" ); //10000
    vegasncallsVIRT    = in.GetNumber ( "vegasncallsVIRT" ); //10000
    makelepcuts        = in.GetBool   ( "makelepcuts"     ); //true
    cubaint            = in.GetBool   ( "cubaint"         ); //true    # integration with     Cuba       Suave
    trapezint          = in.GetBool   ( "trapezint"       ); //false   # trapezoidal rule     for        the     phi_lep     integration     and         semi-analytical for         costh
    quadint            = in.GetBool   ( "quadint"         ); //false   # quadrature  rule     for        the     phi_lep     integration     and         semi-analytical for         costh
    suavepoints        = in.GetNumber ( "suavepoints"     ); //1000000 # number      of       points     for     suave       integration,    newpoints   is              set         to    suavepoints/10;
    nphitrape          = in.GetNumber ( "nphitrape"       ); //1000    # number      of       steps      for     trapezoidal rule            of          phi_lep         integration
    quadnphi           = in.GetNumber ( "quadnphi"        ); //20      # number      of       segments   for     quadrature  rule            of          phi_lep         integration
    ncstart            = in.GetNumber ( "ncstart"         ); //1000    # starting    sampling for        the     costh       semi-analytical integration (common         settings    for   the             trapezoidal and quadrature rules)
    qtrec_naive        = in.GetBool   ( "qtrec_naive"     ); //false
    qtrec_cs           = in.GetBool   ( "qtrec_cs"        ); //false
    qtrec_kt0          = in.GetBool   ( "qtrec_kt0"       ); //true
    timeprofile        = in.GetBool   ( "timeprofile"     ); //false   # debug       and      time       profile resummation integration
    verbose            = in.GetBool   ( "verbose"         ); //false   # debug       and      time       profile costh       phi_lep         integration
    HackBinnerToFiller = in.GetBool   ( "HackBinnerToFiller"         ); //false   # debug       and      time       profile costh       phi_lep         integration
    useGamma           = in.GetBool ( "useGamma" );//
    fiducial           = static_cast<settings::DetFiducial>( in.GetNumber( "fiducial" )); //0
    opts_.approxpdf_    = in.GetNumber ( "opts_approxpdf" ); //0
    opts_.pdfintervals_ = in.GetNumber ( "opts_pdfintervals" ); //100
    opts_.fixedorder_ = fixedorder;

    // additional conditions
    // finite order (NLO vs NNLO)
    if (opts.doLO     && opts.order != 1) throw invalid_argument( "You are trying to run LO term calculation, but order is not 1. Check your input file.") ;
    if (opts.doREAL   && opts.order != 2) throw invalid_argument( "You are trying to run REAL term calculation, but order is not 2. Check your input file.") ;
    if (opts.doVIRT   && opts.order != 2) throw invalid_argument( "You are trying to run VIRT term calculation, but order is not 2. Check your input file.") ;
    // resummation term integration dimension
    if (intDimRes<5 && intDimRes>1){
        int2d = (intDimRes == 2);
        int3d = (intDimRes == 3);
        int4d = (intDimRes == 4);
    } else {
        int2d = false;
        int3d = false;
        int4d = true;
    }

    // counter term integration dimension
    if (intDimCT<4 && intDimCT>1){
        ctint2d = (intDimCT == 2);
        ctint3d = (intDimCT == 3);
        ctintvegas = false;
    } else {
        ctint2d = false;
        ctint3d = false;
        ctintvegas = true;
    }

    if (opts_.approxpdf_ == 1)
      {
	cout << "DYRES-style approximate PDF requested, enforce vegas integration for the resummed cross section" << endl;
	int2d = false;
	int3d = false;
	int4d = true;
      }

    if (fixedorder == true && doRES == true)
      {
	cout << "Asked for fixed order predictions, switching off resummation term" << endl;
	doRES = false;
      }
    else if (doVV == true)
      {
	cout << "Asked for resummed predictions, switching off double virtual term" << endl;
	doVV = false;
      }
    if (fixedorder == true)
      {
	cout << "Asked for fixed order predictions, enforce a_param = 1.0" << endl;
	a_param = 1.0;
      }
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
    nnlo_        . order_     = order        ;         //1              # order
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

    strncpy( part_        . part_      , part         .c_str(), part       .size() ); //virt           # part
    strncpy( lhapdf_char_ . PDFname_   , LHAPDFset    .c_str(), LHAPDFset  .size() ); //CT10nlo.LHgrid
    strncpy( runstring_   . runstring_ , outputfile   .c_str(), outputfile .size() ); //'LHC7-Z-nnlo'  # outputfile

    zcouple_ . q1_ = (useGamma ? -1 :  0 );

    dofill_.doFill_ = 0;
}

void settings::dumpAll(){
    printf("==Listing settings==\n");
    bool print_inputs = true;
    bool print_inputsDYRES = true;
    bool print_masses = true;

    if (print_inputsDYRES) {
        printf("Input DYRES settings:\n");
        dumpD ( "sroot       ",  energy_      . sroot_      ) ;
        dumpI ( "ih1         ",  density_     . ih1_        ) ;
        dumpI ( "ih2         ",  density_     . ih2_        ) ;
        dumpI ( "nproc       ",  nproc_       . nproc_      ) ;
        dumpD ( "mur         ",  scale_       . scale_      ) ;
        dumpD ( "muf         ",  facscale_    . facscale_   ) ;
        dumpD ( "a_param     ",  a_param_     . a_param_    ) ;
        dumpD ( "g_param     ",  g_param_     . g_param_    ) ;
        dumpI ( "order       ",  nnlo_        . order_      ) ;
        dumpS ( "part        ",  part_        . part_       ) ;
        dumpB ( "zerowidth   ",  zerowidth_   . zerowidth_  ) ;
        dumpD ( "M_min       ",  mwminmax_    . Mwmin_      ) ;
        dumpD ( "M_max       ",  mwminmax_    . Mwmax_      ) ;
        dumpI ( "itmx1       ",  iterat_      . itmx1_      ) ;
        dumpI ( "ncall1      ",  iterat_      . ncall1_     ) ;
        dumpI ( "itmx2       ",  iterat_      . itmx2_      ) ;
        dumpI ( "ncall2      ",  iterat_      . ncall2_     ) ;
        dumpI ( "rseed       ",  rseed_       . rseed_      ) ;
        dumpI ( "PDFset      ",  pdfiset_     . iset_       ) ;
        dumpI ( "PDFmember   ",  prefix_      . nset_       ) ;
        dumpS ( "LHAPDFset   ",  lhapdf_char_ . PDFname_    ) ;
        dumpI ( "LHAPDFmember",  lhapdf_int_  . PDFmember_  ) ;
        dumpS ( "outputfile  ",  runstring_   . runstring_  ) ;
        dumpI ( "itmxToFile  ",  pr_          . pr_         ) ;
    }

    if (print_inputs) {
        printf("Input settings:\n");
        dumpD("rmass              ", rmass               );
        dumpD("rwidth             ", rwidth              );
        dumpD("ylow               ", ylow                );
        dumpD("yhigh              ", yhigh               );
        dumpD("mlow               ", mlow                );
        dumpD("mhigh              ", mhigh               );
        dumpI("intDimRes          ", intDimRes           );
        dumpB("int2d              ", int2d               );
        dumpB("int3d              ", int3d               );
        dumpB("int4d              ", int4d               );
        dumpI("intDimCT           ", intDimCT            );
        dumpB("ctint2d            ", ctint2d             );
        dumpB("ctint3d            ", ctint3d             );
        dumpB("ctintvegas         ", ctintvegas          );
        dumpB("fixedorder         ", fixedorder          );
	dumpB("doRES              ", doRES               );
	dumpB("doVV               ", doVV                );
        dumpB("doCT               ", doCT                );
        dumpB("doREAL             ", doREAL              );
        dumpB("doVIRT             ", doVIRT              );
        dumpB("doLO               ", doLO                );
        dumpI("cubaverbosity      ", cubaverbosity       );
        dumpI("cubacores          ", cubacores           );
        dumpI("niterRES           ", niterRES            );
        dumpI("niterCT            ", niterCT             );
        dumpD("vegasncallsRES     ", vegasncallsRES      );
	dumpD("vegasncallsVV      ", vegasncallsVV       );
        dumpD("vegasncallsCT      ", vegasncallsCT       );
        dumpD("vegasncallsLO      ", vegasncallsLO       );
        dumpD("vegasncallsREAL    ", vegasncallsREAL     );
        dumpD("vegasncallsVIRT    ", vegasncallsVIRT     );
        dumpB("makelepcuts        ", makelepcuts         );
        dumpI("fiducial           ", fiducial            );
        dumpB("cubaint            ", cubaint             );
        dumpB("trapezint          ", trapezint           );
        dumpB("quadint            ", quadint             );
        dumpI("suavepoints        ", suavepoints         );
        dumpI("nphitrape          ", nphitrape           );
        dumpI("quadnphi           ", quadnphi            );
        dumpI("ncstart            ", ncstart             );
        dumpB("qtrec_naive        ", qtrec_naive         );
        dumpB("qtrec_cs           ", qtrec_cs            );
        dumpB("qtrec_kt0          ", qtrec_kt0           );
        dumpB("timeprofile        ", timeprofile         );
        dumpB("verbose            ", verbose             );
        dumpB("HackBinnerToFiller ", HackBinnerToFiller  );
        dumpI("approxpdf          ", opts_.approxpdf_    );
        dumpI("pdfintervals       ", opts_.pdfintervals_ );
    }

    if (print_masses){
        printf("Masses and EW constants:\n");
        dumpD ( "md        " , dymasses_   . md_        );
        dumpD ( "mu        " , dymasses_   . mu_        );
        dumpD ( "ms        " , dymasses_   . ms_        );
        dumpD ( "mc        " , dymasses_   . mc_        );
        dumpD ( "mb        " , dymasses_   . mb_        );
        dumpD ( "mt        " , dymasses_   . mt_        );
        dumpD ( "mel       " , dymasses_   . mel_       );
        dumpD ( "mmu       " , dymasses_   . mmu_       );
        dumpD ( "mtau      " , dymasses_   . mtau_      );
        dumpD ( "hmass     " , dymasses_   . hmass_     );
        dumpD ( "hwidth    " , dymasses_   . hwidth_    );
        dumpD ( "wmass     " , dymasses_   . wmass_     );
        dumpD ( "wwidth    " , dymasses_   . wwidth_    );
        dumpD ( "zmass     " , dymasses_   . zmass_     );
        dumpD ( "zwidth    " , dymasses_   . zwidth_    );
        dumpD ( "twidth    " , dymasses_   . twidth_    );
        dumpD ( "mtausq    " , dymasses_   . mtausq_    );
        dumpD ( "mcsq      " , dymasses_   . mcsq_      );
        dumpD ( "mbsq      " , dymasses_   . mbsq_      );
        dumpD ( "Gf_inp    " , ewinput_  . Gf_inp_    );
        dumpD ( "aemmz_inp " , ewinput_  . aemmz_inp_ );
        dumpD ( "xw_inp    " , ewinput_  . xw_inp_    );
        dumpD ( "wmass_inp " , ewinput_  . wmass_inp_ );
        dumpD ( "zmass_inp " , ewinput_  . zmass_inp_ );
        dumpI ( "ewscheme  " , ewscheme_ . ewscheme_  );
        dumpD ( "Vud       " , cabib_    . Vud_       );
        dumpD ( "Vus       " , cabib_    . Vus_       );
        dumpD ( "Vub       " , cabib_    . Vub_       );
        dumpD ( "Vcd       " , cabib_    . Vcd_       );
        dumpD ( "Vcs       " , cabib_    . Vcs_       );
        dumpD ( "Vcb       " , cabib_    . Vcb_       );
        dumpD ( "mb_msbar  " , mb_msbar_ . mb_msbar_  );
        dumpD ( "epinv     " , epinv_    . epinv_     );
        dumpD ( "epinv2    " , epinv2_   . epinv2_    );
    }

    printf("end of setting dump.\n\n");
}

void settings::dumpI( string var,int    val){
    printf( " %s = %d\n", var.c_str(), val);
}
void settings::dumpD( string var,double val){
    printf( " %s = %f\n", var.c_str(), val);
}
void settings::dumpS( string var,string val){
    printf( " %s = %s\n", var.c_str(), val.c_str());
}
void settings::dumpB( string var,bool val){
    printf( " %s = %s\n", var.c_str(), val ? "true" : "false" );
}

int cuts_(double p[4][12], int &njet){
    double p3[4];
    double p4[4];
    for (int i=0; i<4; i++){
        p3[i] = p[i][2];
        p4[i] = p[i][3];
    }
    // because dyres expects oposite logic false=accept event
    return !cuts(p3,p4);
}


/// fiducial cuts
double getPt(double p[4]){
    return sqrt((float)pow(p[0],2)+pow(p[1],2));
}
double getY(double p[4]){
    return 0.5 *log((p[3] + p[2]) / (p[3] - p[2]));
}
double getEta(double p[4]){
    return atanh(p[2]/sqrt((float) p[0]*p[0] + p[1]*p[1] + p[2]*p[2]));
}
double getMt(double p3[4], double p4[4]){
    double p[4];
    p[2]=p3[2]+p4[2];
    p[3]=p3[3]+p4[3];
    return sqrt((float)pow(p[3],2)-pow(p[2],2));
}

bool fiducial_D0(double p3[4], double p4[4]){
    double pt3 = getPt(p3);
    if (pt3<20) return false;
    double aeta3 = fabs(getEta(p3));
    if (aeta3>2.5 || (aeta3>1.1 && aeta3<1.5) ) return false;
    double pt4 = getPt(p4);
    if (opts.nproc==3){ // z
        if (pt4<20) return false;
        double aeta4 = fabs(getEta(p4));
        if (aeta4>2.5 || (aeta4>1.1 && aeta4<1.5) ) return false;
    } else { // W
        if (pt4<25) return false;
    }
    return true; // keep
}

bool fiducial_CDF(double p3[4], double p4[4]){
    double pt3 = getPt(p3);
    if (pt3<20) return false;
    double aeta3 = fabs(getEta(p3));
    if (aeta3>2.0 || (aeta3>1.0 && aeta3<1.2) ) return false;
    if (aeta3<1.0 && pt3<25) return false;
    double pt4 = getPt(p4);
    if (opts.nproc==3){ // z
        if (pt4<20) return false;
        double aeta4 = fabs(getEta(p4));
        if (aeta4>2.0 || (aeta4>1.0 && aeta3<1.2) ) return false;
        if (aeta4<1.0 && pt4<25) return false;
    } else { // W
        if (pt4<25) return false;
    }
    return true; // keep
}

bool fiducial_ATLAS(double p3[4], double p4[4]){
    double pt3 = getPt(p3);
    if (pt3<20) return false;
    double aeta3 = fabs(getEta(p3));
    if (aeta3>2.5 ) return false;
    double pt4 = getPt(p4);
    if (opts.nproc==3){ // z
        if (pt4<20) return false;
        double aeta4 = fabs(getEta(p4));
        if (aeta4>2.5 ) return false;
    } else { // W
        if (pt4<25) return false;
        double mt = getMt(p3,p4);
        if (mt<40) return false;
    }
    return true; // keep
}

bool fiducial_CMS(double p3[4], double p4[4]){
    return fiducial_CDF(p3,p4);
}

bool cuts(double p3[4], double p4[4])
{
  if (!opts.makelepcuts)
    return true;
  switch (opts.fiducial){
      case settings::D0    : fiducial_D0    (p3,p4); break;
      case settings::CDF   : fiducial_CDF   (p3,p4); break;
      case settings::ATLAS : fiducial_ATLAS (p3,p4); break;
      case settings::CMS   : fiducial_CMS   (p3,p4); break;
      case settings::GENEXP :
      default:
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
  }



  //printf("c++ cuts %f %f %f %f  \n", pt3, pt4, y3, y4);

  //if (pt3 < 20)
    //return false;
  //if (pt4 < 20)
    //return false;
  //if (fabs(y3) > 2.4)
    //return false;
  //if (fabs(y4) > 2.4)
    //return false;

  //  if (y3-y4 < 0.2)
  //    return false;

  //  if (66 < m < 116
  //      && pt3 > 20 && pt4 > 20
  //      && fabs (y3) < 2.4 && fabs(y4) < 2.4)
  //    cut = true;
  return true;
}

void binning::readfromfile(const string fname){
    InputParser in(fname);
    qtbins       .clear(); in.GetVectorDouble( "qtbins"      , qtbins       );
    ybins        .clear(); in.GetVectorDouble( "ybins"       , ybins        );
    hist_qt_bins .clear(); in.GetVectorDouble( "plot_qtbins" , hist_qt_bins );
    hist_y_bins  .clear(); in.GetVectorDouble( "plot_ybins"  , hist_y_bins  );
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


InputParser::InputParser( string _filename, string _charset, string _white):
    filename ( _filename   ),
    Ccommnt  ( _charset[0] ),
    Cassign  ( _charset[1] ),
    CopenAr  ( _charset[2] ),
    CdeliAr  ( _charset[3] ),
    CclosAr  ( _charset[4] ),
    Swhite   ( _white      )
{
    // parse default
    try {
      parse_file((string)SHAREDIR+"/default.in");
    } catch (invalid_argument &e1){
        try {
            parse_file("default.in");
        } catch (invalid_argument &e2){
            throw invalid_argument( string( "Can not load default settings:\n")
                    + string(e1.what()) + string("\n")
                    + string(e2.what()) + string("\n")
                    );
        }
    }
    // parse custom settings
    if (filename.compare("")!=0) parse_file(filename);
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
    bool result = (val.compare(0,4,"true") == 0);
    //printf(" InputParser::GetBool processing key `%s` with value `%s` as result `%d`\n", name.c_str(), val.c_str(), result);
    return result;
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


void InputParser::parse_file(const string fname){
    ifstream fstrm(fname.c_str());
    if ( ! fstrm.good() ) throw invalid_argument( string("Uknown file name ") + fname );
    string line;
    // line by line
    while(getline(fstrm,line)) {
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




