#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>

#include "settings.h"
#include "interface.h"

// CXX option parser: https://raw.githubusercontent.com/jarro2783/cxxopts/master/src/cxxopts.hpp
#include "cxxopts.hpp"


settings opts;
binning bins;

void settings::parse_options(int argc, char* argv[]){
    // Declare the supported options.
    po::Options args(argv[0], " [config.in] \n\n"
            " Fast Drell-Yan Monte Carlo and quadrature integrator. \n\n"
            " NOTE: Command line options are overiding the default and config file settings."
            );
    // Hidden arguments
    args.add_options("Hidden")
        ("conf_file"     , "Name of output file. ", po::value<string>()->default_value("")  )
    ;
    // Program options
    args.add_options("")
        ("h,help"            , "Print this help and die."       )
        ("v,verbose"         , "Be verbose"                     )
        ("q,small-stat"      , "Set quick run with small stat." )
        ("p,proc"            , "Set process [z0/wp/wm]"                              , po::value<string>() )
        ("c,collider"        , "Set beam conditions [tev2/lhc7/lhc8/lhc13/lhc14]"    , po::value<string>() )
        ("t,term"            , "Set term [REAL,VIRT,CT,..]"                          , po::value<string>() )
        ("r,seed"            , "Set random seed [integer]"                           , po::value<int>()    )
        ("s,pdfset"          , "Set PDF set [LHAPDF name]"                           , po::value<string>() )
        ("m,pdfvar"          , "Set PDF member [integer/all]"                        , po::value<string>() )
        ("qtbins"            , "Set equdistan binning for mass [N,lo,hi]"            , po::value<string>() )
        ("ybins"             , "Set equdistan binning for qt [N,lo,hi]"              , po::value<string>() )
        ("mbins"             , "Set equdistan binning for y [N,lo,hi]"               , po::value<string>() )

    ;
    // Parse
    try {
        args.parse_positional( std::vector<string>({"conf_file"}) );
        args.parse(argc,argv);
    }
    catch (po::OptionException &e){
        printf("%s\n", args.help().c_str());
        printf("Bad arguments: %s \n",e.what());
        throw e;
    }
    // Print help and die
    if (args.count("help")) {
        throw QuitProgram(args.help().c_str());
    }

    // load config file (or default settings)
    readfromfile      ( args["conf_file"].as<string>() );
    bins.readfromfile ( args["conf_file"].as<string>() );

    // NOTE: Command line options are overiding the default and config file settings.
    // verbose
    if (args.count("verbose")) verbose=true;
    // rseed
    if (args.count("seed")) rseed=args["seed"].as<int>();
    // small stat
    if (args.count("small-stat")){
        niterRES           = 1;
        niterCT            = 1;
        niterVJ            = 1;
        vegasncallsRES     = 1e3;
        vegasncallsVV      = 1e5;
        vegasncallsCT      = 1e4;
        vegasncallsLO      = 1e5;
        vegasncallsREAL    = 1e5;
        vegasncallsVIRT    = 1e5;
        pcubaccuracy       = 0.5;
    } 

    // proc
    if (args.count("proc")) {
        // FIXME: is it only needed to set proc or also mass and scale values?
        string val=args["proc"].as<string>(); 
        ToLower(val);
        if        (val == "z0"){ nproc=3;
        } else if (val == "wp"){ nproc=1;
        } else if (val == "wm"){ nproc=2;
        } else {
            throw QuitProgram("Unsupported value of proc: "+val);
        }
    }
    
    // PDF
    if (args.count("pdfset")) LHAPDFset=args["pdfset"].as<string>();
    if (args.count("pdfvar")) {
        // set default values
        PDFerrors=false;
        LHAPDFmember=0;
        string val=args["pdfvar"].as<string>();
        ToLower(val);
        if   (val=="all") { PDFerrors=true;
        } else if (IsNumber(val)) { LHAPDFmember=stod(val);
        } else {
            throw QuitProgram("Unsupported value of pdfvar: "+val);
        }
    }

    // Collider
    if (args.count("collider")) {
        string val=args["collider"].as<string>();
        ToLower(val);
        if        (val == "tev1"  ){ sroot=1.80e3; ih1=1; ih2=-1;
        } else if (val == "tev2"  ){ sroot=1.96e3; ih1=1; ih2=-1;
        } else if (val == "lhc7"  ){ sroot=7.00e3; ih1=1; ih2=1;
        } else if (val == "lhc8"  ){ sroot=8.00e3; ih1=1; ih2=1;
        } else if (val == "lhc13" ){ sroot=13.0e3; ih1=1; ih2=1;
        } else if (val == "lhc14" ){ sroot=14.0e3; ih1=1; ih2=1;
        } else {
            throw QuitProgram("Unsupported value of collider: "+val);
        }
    }

    // term
    if (args.count("term")) {
        // first turn off all terms
        doRES = doVV = doCT = doREAL = doVIRT = doLO = doVJ = false ;
        string val=args["term"].as<string>();
        ToUpper(val);
        for (auto piece : Tokenize(val)) {
            if        ( piece == "REAL"     ) { doREAL=true;
            } else if ( piece ==  "VIRT"    ) { doVIRT=true;
            } else if ( piece ==  "VV"      ) { doVV=true;   fixedorder=true;
            } else if ( piece ==  "LO"      ) { doLO=true;
            } else if ( piece ==  "VJ"      ) { doVJ=true;
            } else if ( piece ==  "RES"     ) { doRES=true;  intDimRes=8;
            } else if ( piece ==  "RES3D"   ) { doRES=true;  intDimRes=3;
            } else if ( piece ==  "RES2D"   ) { doRES=true;  intDimRes=2;
            } else if ( piece ==  "CT"      ) { doCT=true;   intDimCT=6; //< FIXME: ia 6D save ?
            } else if ( piece ==  "CT3D"    ) { doCT=true;   intDimCT=3;
            } else if ( piece ==  "CT2D"    ) { doCT=true;   intDimCT=2;
            } else if ( piece ==  "FIXCT"   ) { doCT=true;   intDimCT=6; fixedorder=true; 
            } else if ( piece ==  "FIXCT3D" ) { doCT=true;   intDimCT=3; fixedorder=true; 
            } else if ( piece ==  "FIXCT2D" ) { doCT=true;   intDimCT=2; fixedorder=true; 
            } else {
                throw QuitProgram("Unsupported value of term : "+piece);
            }
        }
    }

    //binning
    parse_binning("qtbins" , bins.qtbins ,args);
    parse_binning("ybins"  , bins.ybins  ,args);
    parse_binning("mbins"  , bins.mbins  ,args);


    // check consistency of settings
    check_consitency();
}

void settings::readfromfile(const string fname){
    //read input settings from file
    InputParser in(fname);
    sroot          = in.GetNumber ( "sroot"          ); //7e3
    ih1            = in.GetNumber ( "ih1"            ); //1
    ih2            = in.GetNumber ( "ih2"            ); //1              # ih1,        ih2
    nproc          = in.GetNumber ( "nproc"          ); //3              # nproc
    kmuren         = in.GetNumber ( "kmuren"            ); //91.1876e0
    kmufac         = in.GetNumber ( "kmufac"            ); //91.1876e0      # mur,        muf
    kmures         = in.GetNumber ( "kmures"            ); //91.1876e0      # mur,        muf
    C1             = in.GetNumber ( "C1"            ); //91.1876e0      # mur,        muf
    C3             = in.GetNumber ( "C3"            ); //91.1876e0      # mur,        muf
    //    a_param        = in.GetNumber ( "a_param"        ); //2.0e0          # a_param
    g_param        = in.GetNumber ( "g_param"        ); //1.0e0          # g_param
    order          = in.GetNumber ( "order"          ); //1              # order
    zerowidth      = in.GetBool   ( "zerowidth"      ); //false          # zerowidth
    dynamicscale   = in.GetBool   ( "dynamicscale"   );
    dynamicresscale= in.GetBool   ( "dynamicresscale"   );
    rseed          = in.GetNumber ( "rseed"          ); //123456         # rseed
    blim           = in.GetNumber ( "blim"           );
    LHAPDFset      = in.GetString ( "LHAPDFset"      ); //CT10nlo.LHgrid
    LHAPDFmember   = in.GetNumber ( "LHAPDFmember"   ); //0              # set,        member   (LHAPDFs)
    Gf             = in.GetNumber ( "Gf"        );
    zmass          = in.GetNumber ( "zmass"        );
    wmass          = in.GetNumber ( "wmass"        );
    xw             = in.GetNumber ( "xw"        );
    aemmz          = in.GetNumber ( "aemmz"        );
    zwidth         = in.GetNumber ( "zwidth"        );
    wwidth         = in.GetNumber ( "wwidth"        );
    Vud            = in.GetNumber ( "Vud"        );
    Vus            = in.GetNumber ( "Vus"        );
    Vub            = in.GetNumber ( "Vub"        );
    Vcd            = in.GetNumber ( "Vcd"        );
    Vcs            = in.GetNumber ( "Vcs"        );
    Vcb            = in.GetNumber ( "Vcb"        );
    Zuu            = in.GetNumber ( "Zuu"        );
    Zdd            = in.GetNumber ( "Zdd"        );
    Zcc            = in.GetNumber ( "Zcc"        );
    Zss            = in.GetNumber ( "Zss"        );
    Zbb            = in.GetNumber ( "Zbb"        );
    //ylow           = in.GetNumber ( "ylow"            ); //2
    //yhigh          = in.GetNumber ( "yhigh"           ); //2.4
    //mlow           = in.GetNumber ( "mlow"            ); //66.
    //mhigh          = in.GetNumber ( "mhigh"           ); //116.
    dampk          = in.GetNumber ( "dampk"           );
    dampdelta      = in.GetNumber ( "dampdelta"           );
    dampmode       = in.GetNumber ( "dampmode"           );
    qtcutoff       = in.GetNumber ( "qtcutoff"           );
    xqtcut         = in.GetNumber ( "xqtcut"           );
    qtcut          = in.GetNumber ( "qtcut"           );
    intDimRes      = in.GetNumber ( "intDimRes"       ); //3
    intDimCT           = in.GetNumber ( "intDimCT"       ); //3
    //    intDimVJ           = in.GetNumber ( "intDimVJ"        );
    fixedorder         = in.GetBool   ( "fixedorder"      ); //false
    doRES              = in.GetBool   ( "doRES"           ); //false
    doVV               = in.GetBool   ( "doVV"            ); //false
    doCT               = in.GetBool   ( "doCT"            ); //false
    doREAL             = in.GetBool   ( "doREAL"          ); //false
    doVIRT             = in.GetBool   ( "doVIRT"          ); //false
    doLO               = in.GetBool   ( "doLO"            ); //false
    doVJ               = in.GetBool   ( "doVJ"            ); //false
    cubaverbosity      = in.GetNumber ( "cubaverbosity"   ); //0       # Cuba        info     messsages, from    0           to              3
    cubacores          = in.GetNumber ( "cubacores"       ); //0
    niterRES           = in.GetNumber ( "niterRES"        ); //0       # only        for      2d         and     3d          cuhre           integration
    niterCT            = in.GetNumber ( "niterCT"         ); //0       # only        for      3d          cuhre           integration
    niterVJ            = in.GetNumber ( "niterVJ"         ); //0       # only        for      3d          cuhre           integration
    vegasncallsRES     = in.GetNumber ( "vegasncallsRES"  ); //10000
    vegasncallsVV      = in.GetNumber ( "vegasncallsVV"   ); //10000
    vegasncallsCT      = in.GetNumber ( "vegasncallsCT"   ); //10000
    vegasncallsLO      = in.GetNumber ( "vegasncallsLO"   ); //10000
    vegasncallsREAL    = in.GetNumber ( "vegasncallsREAL" ); //10000
    vegasncallsVIRT    = in.GetNumber ( "vegasncallsVIRT" ); //10000
    pcubature          = in.GetBool   ( "pcubature" );
    pcubaccuracy       = in.GetNumber ( "pcubaccuracy" );
    makelepcuts        = in.GetBool   ( "makelepcuts"     ); //true
    lptcut             = in.GetNumber ( "lptcut"          );
    lycut              = in.GetNumber ( "lycut"          );
    mtcut              = in.GetNumber ( "mtcut"          );
    etmisscut          = in.GetNumber ( "etmisscut"          );
    l1ptcut            = in.GetNumber ( "l1ptcut"          );
    l1ycut             = in.GetNumber ( "l1ycut"          );
    l2ptcut            = in.GetNumber ( "l2ptcut"          );
    l2ycut             = in.GetNumber ( "l2ycut"          );
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
    resumcpp           = in.GetBool   ( "resumcpp"        );
    useGamma           = in.GetBool ( "useGamma" );//
    fiducial           = static_cast<cuts::DetFiducial>( in.GetNumber( "fiducial" )); //0
    PDFerrors           = in.GetBool ( "PDFerrors" );//
    opts_.approxpdf_    = in.GetNumber ( "opts_approxpdf" ); //0
    opts_.pdfintervals_ = in.GetNumber ( "opts_pdfintervals" ); //100
    evolmode           = in.GetNumber  ("evolmode");
    opts_.fixedorder_  = fixedorder;
    bintaccuracy       = in.GetNumber ( "bintaccuracy" );
    mellinintervals    = in.GetNumber ( "mellinintervals" );
    mellinrule         = in.GetNumber ( "mellinrule" );
    zmax               = in.GetNumber ( "zmax" );
    mellincores        = in.GetNumber ( "mellincores" );
    yintervals         = in.GetNumber ( "yintervals" );
    yrule              = in.GetNumber ( "yrule" );
    ptbinwidth         = in.GetBool ( "ptbinwidth" );
    ybinwidth          = in.GetBool ( "ybinwidth" );

    return ;
}

void settings::check_consitency(){

    // additional conditions
    if (order != 1 && order != 2)
      throw invalid_argument("Invalid order, please select 1 (NLO) or 2 (NNLO)");
    if (nproc != 1 && nproc != 2 && nproc != 3)
      throw invalid_argument("Wrong process, please select nproc = 1 (W+), 2 (W-), or 3(Z)");

    // finite order (NLO vs NNLO)
    /*
    if (opts.doLO     && opts.order != 1) throw invalid_argument( "You are trying to run LO term calculation, but order is not 1. Check your input file.") ;
    if (opts.doREAL   && opts.order != 2) throw invalid_argument( "You are trying to run REAL term calculation, but order is not 2. Check your input file.") ;
    if (opts.doVIRT   && opts.order != 2) throw invalid_argument( "You are trying to run VIRT term calculation, but order is not 2. Check your input file.") ;
    */

    if (order != 1)
      opts.doLO = false;
    if (order != 2)
      {
	opts.doREAL = false;
	opts.doVIRT = false;
      }
    //Asked for fixed order predictions, switching off resummation term
    if (fixedorder == true)
      doRES = false;
    //Asked for resummed predictions, switching off double virtual term
    if (fixedorder == false)
      doVV = false;
    //In fixed order mode, a_param must be one
    if (fixedorder == true)
      {
	cout << "Asked for fixed order predictions, enforce kmures = 1.0" << endl;
	kmures = 1.0;
      }

    if (evolmode > 4 || evolmode < 0)
      {
	cout << "wrong value for evolmode: available evolmodes: 0,1,2,3,4" << endl;
	exit (-1);
      }

    if (dynamicscale == true && evolmode < 3)
      {
	//cannot use a dynamic muren, mufac, when the PDFs are converted from x-space to N-space at the factorisation scale
	cout << "dynamicscale possible only with evolmode = 3 or 4" << endl;
	exit (-1);
      }

    if (PDFerrors == true && LHAPDFmember != 0)
      {
	cout << "Asked for PDFerrors, enforce LHAPDFmember  = 0" << endl;
	LHAPDFmember = 0;
      }

    if (qtcut <= 0 && xqtcut <= 0)
      {
	cout << "At least one between qtcut and xqtcut must be > 0" << endl;
	exit (-1);
      }

    // resummation term integration dimension
    if (intDimRes<4 && intDimRes>1){
        resint2d = (intDimRes == 2);
        resint3d = (intDimRes == 3);
        resintvegas = false;
    } else {
        resint2d = false;
        resint3d = false;
        resintvegas = true;
    }

    // counter term integration dimension
    if (intDimCT<4 && intDimCT>1){
        ctint2d = (intDimCT == 2);
        ctint3d = (intDimCT == 3);
        ctintvegas6d = false;
	ctintvegas8d = false;
    } else {
        ctint2d = false;
        ctint3d = false;
        ctintvegas6d = (intDimCT <= 6);
	ctintvegas8d = (intDimCT > 6);
    }

    /*
    // V+J integration dimension
    if (intDimVJ < 4 && intDimVJ > 2){
        vjint3d = (intDimVJ == 3);
        vjintvegas7d = false;
    } else {
        vjint3d = false;
        vjintvegas7d = true;
      }
   if (makelepcuts)
      {
	cout << "Required cuts on the final state leptons, enforce vegas integration for V+J fixed order cross section" << endl;
	vjint3d = false;
	vjintvegas7d = true;
      }
    */


    if (opts_.approxpdf_ == 1)
      {
	cout << "DYRES-style approximate PDF requested, enforce vegas integration for the resummed cross section" << endl;
	resint2d = false;
	resint3d = false;
	resintvegas = true;
      }

    // -- binning
    // check bins size
    if ( bins.qtbins .size() < 2) throw QuitProgram("Option `qtbins` needs at least 2 items ");
    if ( bins.ybins  .size() < 2) throw QuitProgram("Option `ybins`  needs at least 2 items ");
    if ( bins.mbins  .size() < 2) throw QuitProgram("Option `mbins`  needs at least 2 items ");
    // check sorting
    sort( bins.qtbins .begin () , bins.qtbins .end () ) ;
    sort( bins.ybins  .begin () , bins.ybins  .end () ) ;
    sort( bins.mbins  .begin () , bins.mbins  .end () ) ;
    // plotmode
    ToLower(bins.plotmode);
    // set histogram bins
    bins.hist_qt_bins = bins.qtbins ;
    bins.hist_y_bins  = bins.ybins ;
    bins.hist_m_bins  = bins.mbins ;
    // integration boundaries
    // FIXME: where is this used ?
    ylow  = bins.ybins.front();
    yhigh = bins.ybins.back();
    mlow  = bins.mbins.front();
    mhigh = bins.mbins.back();
    // plot mode consitency with integration
    if ( bins.plotmode=="fill" && 
         ( resint2d || resint3d || ctint2d || ctint3d )
       ) {
        printf ("Warning: plotmode: Filling not work for quadrature integration. I am switching to integrate.\n");
        bins.plotmode="integrate";
    }
    if (bins.plotmode=="fill"){ 
        // integration boundaries are max and minimum
        bins.qtbins = { bins.qtbins .front(), bins.qtbins .back() };
        bins.ybins  = { bins.ybins  .front(), bins.ybins  .back() };
        bins.mbins  = { bins.mbins  .front(), bins.mbins  .back() };
    } else if (bins.plotmode=="integrate"){
        // keep integration bins same as plotting (quadrature)
    } else {
        throw QuitProgram("Unsupported value of plotmode : "+bins.plotmode);
    }

    return ;
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
	//        dumpD ( "a_param     ",  a_param_     . a_param_    ) ;
        dumpD ( "g_param     ",  g_param_     . g_param_    ) ;
        dumpI ( "order       ",  nnlo_        . order_      ) ;
        dumpB ( "zerowidth   ",  zerowidth_   . zerowidth_  ) ;
    }

    if (print_inputs) {
        printf("Input settings:\n");
	dumpB ( "dynamicscale"      , dynamicscale        );
	dumpB ( "dynamicresscale"      , dynamicresscale        );
        dumpD ( "kmuren      ",  kmuren                     ) ;
        dumpD ( "kmufac      ",  kmufac                     ) ;
	dumpD ( "kmures       ",  kmures   ) ;
	dumpD ( "C1       ",  C1   ) ;
	dumpD ( "C3       ",  C3   ) ;
        dumpD( "blim              ",  blim    ) ;
        dumpS("LHAPDFset          ", LHAPDFset           );
        dumpI("LHAPDFmember       ", LHAPDFmember        );
        dumpI("rseed              ", rseed               );
	dumpD("Gf"                 , Gf);
	dumpD("zmass"              , zmass);
	dumpD("wmass"              , wmass   );
	dumpD("xw"                 , xw);
	dumpD("aemmz"              , aemmz);
	dumpD("zwidth"             , zwidth);
	dumpD("wwidth"             , wwidth);
	dumpD( "Vud",        Vud);
	dumpD( "Vus",        Vus);
	dumpD( "Vub",        Vub);
	dumpD( "Vcd",        Vcd);
	dumpD( "Vcs",        Vcs);
	dumpD( "Vcb",        Vcb);
	dumpD( "Zuu",        Zuu);
	dumpD( "Zdd",        Zdd);
	dumpD( "Zss",        Zss);
	dumpD( "Zcc",        Zcc);
	dumpD( "Zbb",        Zbb);
        dumpD("ylow               ", ylow                );
        dumpD("yhigh              ", yhigh               );
        dumpD("mlow               ", mlow                );
        dumpD("mhigh              ", mhigh               );
        dumpD("dampk",             dampk       );
        dumpD("dampdelta",     dampdelta      );
        dumpI("dampmode",       dampmode     );
	dumpD("qtcutoff",       qtcutoff     );
	dumpD("xqtcut",       xqtcut     );
	dumpD("qtcut",       qtcut     );
        dumpB("useGamma           ", useGamma            );
        dumpI("intDimRes          ", intDimRes           );
        dumpB("resint2d           ", resint2d            );
        dumpB("resint3d           ", resint3d            );
        dumpB("resintvegas        ", resintvegas         );
        dumpI("intDimCT           ", intDimCT            );
        dumpB("ctint2d            ", ctint2d             );
        dumpB("ctint3d            ", ctint3d             );
        dumpB("ctintvegas6d         ", ctintvegas6d          );
	dumpB("ctintvegas8d         ", ctintvegas8d          );
        dumpB("fixedorder         ", fixedorder          );
	dumpB("doRES              ", doRES               );
	dumpB("doVV               ", doVV                );
        dumpB("doCT               ", doCT                );
        dumpB("doREAL             ", doREAL              );
        dumpB("doVIRT             ", doVIRT              );
        dumpB("doLO               ", doLO                );
        dumpB("doVJ               ", doVJ                );
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
	dumpB("pcubature          ", pcubature     );
	dumpD("pcubaccuracy       ", pcubaccuracy     );
        dumpB("makelepcuts        ", makelepcuts         );
	dumpD("lptcut        ", lptcut         );
	dumpD("lycut        ", lycut         );
	dumpD("mtcut        ", mtcut         );
	dumpD("etmisscut        ", etmisscut         );
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
	dumpB("resumcpp           ", resumcpp            );
        dumpI("approxpdf          ", opts_.approxpdf_    );
        dumpI("pdfintervals       ", opts_.pdfintervals_ );
	dumpI("evolmode           ", evolmode            );
        dumpB("PDFerrors          ", PDFerrors           );
        dumpI("mellinintervals    ", mellinintervals     );
        dumpI("mellinrule         ", mellinrule          );
	dumpD("zmax               ", zmax          );
	dumpD("mellincores        ", mellincores          );
        dumpI("yintervals         ", yintervals     );
        dumpI("yrule              ", yrule          );
        dumpB("ptbinwidth         ", ptbinwidth     );
        dumpB("ybinwidth          ", ybinwidth         );
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

vector<string> settings::Tokenize(string val,char Delim){
    vector<string> vec;
    size_t pos = 0;
    string tmp;
    while (!val.empty()){
        // find delim
        pos = val.find_first_of(Delim);
        if (pos==string::npos){
            // last item
            tmp = val;
            val.clear();
        } else {
            // get substring
            tmp = val.substr(0,pos);
            val = val.substr(pos+1);
        }
        // add to vector
        if (!tmp.empty()) vec.push_back(tmp);
    }
    return vec;
}

bool settings::IsNumber(const string &s) {
    double dummy;
    try {
        dummy = stod(s);
        return true;
    } catch (const std::exception &e){
        return false;
    }
}

void settings::parse_binning(string name, vector<double> &bins, po::Options &args){
    if (args.count(name)) {
        string e("Unsupported value of "+name+" : need 3 numbers seperated by comma: 'N,lo,hi' ");
        string val=args[name.c_str()].as<string>();
        ToLower(val);
        vector<string> vec = Tokenize(val);
        // check value
        if (vec.size()!=3) throw QuitProgram(e+" not 3 numbers.");
        for (auto s : vec) if (!IsNumber(s))  throw QuitProgram(e+" not numbers.");
        // retrieve N,lo,hi
        int N = stod(vec[0]);
        double lo = stod(vec[1]);
        double hi = stod(vec[2]);
        if (lo>hi)  throw QuitProgram(e+" lo is more then hi.");
        if (N<1)  throw QuitProgram(e+" N is not at least 1.");
        // make binning
        bins.clear();
        for (double loedge=lo; loedge<=hi; loedge+=(hi-lo)/double(N)) bins.push_back(loedge);
    }
}


void binning::readfromfile(const string fname){
    InputParser in(fname);
    plotmode =  in.GetString( "plotmode");
    qtbins       .clear(); in.GetVectorDouble( "qt_bins"      , qtbins       );
    ybins        .clear(); in.GetVectorDouble( "y_bins"       , ybins        );
    mbins        .clear(); in.GetVectorDouble( "m_bins"       , mbins        );
    hist_qt_bins .clear();
    hist_y_bins  .clear();
    hist_m_bins  .clear();
    return;
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
    try { return stod(val);
    } catch (const std::exception &e){
        printf("Cannot read option '%s' with value '%s' as number.\n",name.c_str(), val.c_str() );
        throw e;
    }
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
    size_t pos = 0;
    string tmp;
    while (!val.empty()){
        // find delim
        pos = val.find(CdeliAr);
        if (pos==string::npos){
            // last item
            tmp = val;
            val.clear();
        } else {
            // get substring
            tmp = val.substr(0,pos);
            val = val.substr(pos+1);
        }
        // get substring
        trim(tmp);
        // expecting number otherwise exception
        if (!tmp.empty()) vec.push_back(stod(tmp));
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
