#ifndef print_table_C
#define print_table_C
/**
 * @file print_table.C
 * @brief A brief description
 *
 * Detailed description of file
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>, Stefano Camarda <Stefano.Camarda@cern.ch>
 * @date 2016-09-29
 */

#include "dyturbo.h"
#include "settings.h"
#include "dyres_interface.h"
#include "interface.h"
#include "coupling.h"
#include "scales.h"
#include "HistoHandler.h"
#include "clock_real.h"

#include <LHAPDF/LHAPDF.h>

#include<iostream>
#include<limits>
#include<fstream>
#include<iomanip>

using std::setw;
using std::flush;
using std::cout;
using std::endl;

namespace DYTurbo {
    namespace PrintTable {

      SStream textfile("empty", std::ios_base::ate);

      bool useUnicode;
      String c_hor;
      String c_ver;
      String c_ul;
      String c_uc;
      String c_ur;
      String c_ml;
      String c_mc;
      String c_mr;
      String c_bl;
      String c_bc;
      String c_br;

      size_t eoc;
      size_t wb;
      size_t wr;
      size_t we;
      size_t wx;
      size_t wt;
      size_t bound_width;
      size_t time_width;
      size_t term_width;

      void Init()
      {
	useUnicode=opts.unicode;
	c_hor = useUnicode ? "\u2501" : "-";
	c_ver = useUnicode ? "\u2503" : "|";
	c_ul = useUnicode  ? "\u250f" : "/";
	c_uc = useUnicode  ? "\u2530" : "-";
	c_ur = useUnicode  ? "\u2513" : "\\";
	c_ml = useUnicode  ? "\u2523" : "|";
	c_mc = useUnicode  ? "\u254b" : "+";
	c_mr = useUnicode  ? "\u252b" : "|";
	c_bl = useUnicode  ? "\u2517" : "\\";
	c_bc = useUnicode  ? "\u253b" : "-";
	c_br = useUnicode  ? "\u251b" : "/";

	eoc=2; //!< end of cell length
	wb=5; //!< bound length
	wr=8; //!< result length
	we=6; //!< error length
	wx=(useUnicode ? 7 : 5 ); //!< exponent length
	wt=6; //!< time length
	bound_width = 2*wb+3;
	time_width = 2+ // parenthesis
	  wt+
	  1; // parenthesis
	term_width = wr + // result
	  (useUnicode ? 3:4) + // plus minus
	  we + // error
	  wx + // exponenent
	  time_width;
      }
      
        inline void EndOfCell()  { std::cout << " " << c_ver << flush; }
        inline void BeginOfRow() { std::cout << c_ver << flush; }
        inline void EndOfRow()   { cout << endl;          }

        void coutN(size_t N, String s) { for (size_t i = 0; i < N; ++i) { cout << s; }}

        void Hline(char p='m'){
            String left   = c_ml;
            String center = c_mc;
            String right  = c_mr;
            if (p=='u') {
                left   = c_ul;
                center = c_uc;
                right  = c_ur;
            }
            if (p=='b') {
                left   = c_bl;
                center = c_bc;
                right  = c_br;
            }
            size_t N_loopingBounds =0;
            for (auto &a: ActiveBoundaries ) if (a.size()>2) N_loopingBounds++;
            cout << left;
            size_t wtotal = 1; // begin of line
            wtotal += (bound_width +eoc)*N_loopingBounds;
            for (size_t i = 0; i < N_loopingBounds; ++i) {
                coutN(bound_width+1,c_hor);
                cout << center;
            }
            for (size_t i = 0; i < ActiveTerms.size(); ++i) {
                coutN(term_width+1,c_hor);
                cout << center;
            }
            // total column
            coutN(term_width+1,c_hor);
            cout << right << endl;
        };

        void BoundsAllLooping(bool printNames=false, bool useFullBound=false){
            // loop over bounds
            for (size_t ibound=0; ibound < N_boundaries; ++ibound){
                if (ActiveBoundaries[ibound].size()<3) continue;
                if (printNames){
                    cout << setw(bound_width) << ActiveBoundaries[ibound].name ;
                    if (opts.texttable) {
                        textfile << ActiveBoundaries[ibound].name << "lo ";
                        textfile << ActiveBoundaries[ibound].name << "hi ";
                    } 
                } else {
                    double lo = useFullBound ? ActiveBoundaries[ibound].front() : last_bounds.loBound(ibound);
                    double hi = useFullBound ? ActiveBoundaries[ibound].back()  : last_bounds.hiBound(ibound);
                    cout << setw(wb) << lo ;
                    cout << " - ";
                    cout << setw(wb) << hi ;
                    if (opts.texttable) textfile << lo << " " << hi << " ";
                }
                EndOfCell();
            }
        }

        void SettingsHeader(String title="");

        void IntegrationSettings(){
            SettingsHeader("Integration settings");
            bool fixed_born = opts.doBORN && opts.fixedorder;
            bool fixed_born_cubature = fixed_born && opts.bornint2d ;
            bool resum_born = opts.doBORN && !opts.fixedorder;
            bool resum_born_cubature = resum_born && (opts.resint3d || opts.resint2d);
            bool ct_cubature = opts.doCT && (opts.ctint2d || opts.ctint3d);
            cout << Col4 ( "Random seeds"   , opts.rseed     , "" , "" ) ;
            cout << Col4 ( "Parallel cores" , opts.cubacores , "" , "" ) ;
            cout << endl;
            for (TermIterator iterm;!iterm.IsEnd();++iterm){
                (*iterm).Print();
            }
            cout << endl;
            if (resum_born && opts_.approxpdf_==0){
                cout << Col4 ( "Mellin inverse transform:" , "gaussian "      , "nodes ="     , opts.mellinrule      ) ;
                cout << Col4 ( ""                          , ""               , "intervals =" , opts.mellinintervals ) ;
                cout << Col4 ( ""                          , "imaginary axis" , "zmax ="      , opts.zmax            ) ;
            }
            if ( resum_born_cubature || ct_cubature || fixed_born_cubature ){
                if (opts.cubaint)
                    cout << Col4( "Angular variables:" , "Suave in dcosth,dphi" , "ncalls =" , opts.suavepoints );
                else if (opts.quadint) {
                    cout << Col4( "Angular variable costh:" , "semi-analytical", "ncstart ="   , opts.ncstart      );
                    cout << Col4( "Angular variables phi:"  , "gaussian"       , "intervals =" , opts.phiintervals );
                }
                else if (opts.trapezint) {
                    cout << Col4( "Angular variables costh:" , "semi-analytical" , "ncstart =" , opts.ncstart   );
                    cout << Col4( "Angular variables phi:"   , "trapezoidal"     , "points ="  , opts.nphitrape );
                }
            }
            cout << endl;
            string col1 = "Constant boundaries";
            for (size_t ibound = 0; ibound < N_boundaries; ++ibound){
                if (ActiveBoundaries[ibound].size() > 2) continue;
                string col2 = ActiveBoundaries[ibound].name;
                cout << Col4 ( col1 , col2 , "low ="  , ActiveBoundaries [ ibound ] .front (  )  ) ;
                cout << Col4 ( ""   , ""   , "high =" , ActiveBoundaries [ ibound ] .back  (  )  ) ;
                col1 = "";
            }
            cout << endl;
        }

        inline void TermName(Term &term){
            cout << setw(term_width) << term.name;
            EndOfCell();
        }

        void Header() {
            Hline('u');
            BeginOfRow();
            if (opts.texttable){
                textfile.str("#"); // start of text file
                textfile.precision(std::numeric_limits<double>::digits10+1);
            }
            BoundsAllLooping(true); // print only names
            for( TermIterator term; !term.IsEnd(); ++term){
                TermName(*term);
            }
            TermName(subtotal);
            if (opts.texttable){
                for (int ipdf = 0; ipdf < opts.totpdf; ++ipdf) textfile << "PDF" << ipdf << " ";
                textfile << " uncertainty\n";
            }
            EndOfRow();
            Hline();
            subtotal.last_reset();
        };



        void Bounds(bool use_full_bound) {
            BeginOfRow();
            BoundsAllLooping(false,use_full_bound);
        };


        String engineering_base(double value, int dec_max, bool sign){
            // Decimal order where we print (==exponent)
            int decimal_print = (dec_max/3)*3;
            // Prepare number for printing
            double fac = pow (10, decimal_print);
            value/=fac;
            // Prepare number for printing
            SStream tmp;
            // leading sign
            if (sign && value >=0) tmp << "+";
            tmp << value;
            // if (decimal_print!=0 && value!=0){
            //     tmp << "e";
            //     if (decimal_print>0) tmp << "+"; // always print sign on exponent
            //     tmp <<setw(2)<<setfill('0')<<decimal_print; // two digits leading with zero
            // }
            return tmp.str();
        }

        String unicode_exponent(int exp){
            char ch_exp[10];
            sprintf(ch_exp,"%d",exp+100);
            String tmp;
            for ( size_t ich = 1; ich<3; ++ich){
                if      (ch_exp[ich] == '0' ) tmp += "\u2070";
                else if (ch_exp[ich] == '1' ) tmp += "\u00B9";
                else if (ch_exp[ich] == '2' ) tmp += "\u00B2";
                else if (ch_exp[ich] == '3' ) tmp += "\u00B3";
                else if (ch_exp[ich] == '4' ) tmp += "\u2074";
                else if (ch_exp[ich] == '5' ) tmp += "\u2075";
                else if (ch_exp[ich] == '6' ) tmp += "\u2076";
                else if (ch_exp[ich] == '7' ) tmp += "\u2077";
                else if (ch_exp[ich] == '8' ) tmp += "\u2078";
                else if (ch_exp[ich] == '9' ) tmp += "\u2079";
                else tmp += " ";
            }
            return tmp;
        }

        String engineering_exponent(int dec_max,bool unit=false){
            // Decimal order where we print (==exponent)
            int decimal_print = (dec_max/3)*3;
            SStream tmp;
            if (unit && decimal_print > 0 && decimal_print < 12 ){
                if (decimal_print < 3)  tmp << "fb";
                else if (decimal_print < 6)  tmp << "pb";
                else if (decimal_print < 9)  tmp << "nb";
                else if (decimal_print < 12) tmp << "ub";
            } else{
                if (decimal_print!=0){
                    if (useUnicode){
                        //tmp << " \u00d7"; // mult
                        tmp << " \u00b7"; // dot
                        tmp << "10";
                        if (decimal_print>0) tmp << "\u207A"; // always print sign on exponent
                        if (decimal_print<0) tmp << "\u207B"; // always print sign on exponent
                        tmp << unicode_exponent(abs(decimal_print));
                    } else {
                        tmp << " e";
                        if (decimal_print>0) tmp << "+"; // always print sign on exponent
                        if (decimal_print<0) tmp << "-"; // always print sign on exponent
                        tmp <<setw(2)<<setfill('0')<<abs(decimal_print); // two digits leading with zero
                    }
                } else tmp << " ";
            }
            return tmp.str();
        }

        // Rounding to significant digits
        void round_by_error(double &value, double &error, int nsd = 2 /*number of significant digits*/){
	  if (value == 0 && error == 0)
	    return;
            // Decimal order where we round
            //   If no error, take significant digits of value
            int decimal_round = ceil(log10(fabs(  (error==0) ? value : error ))) - nsd;
            double fac = pow (10, decimal_round);
            // round numbers
            value = roundf(value/fac)*fac;
            error = roundf(error/fac)*fac;
            // Put correct decimal base
        }


        VecStr result_format(double value, double error=0, bool unit=false) {
            // make sure error always positive
            error=fabs(error);
            // round by significant digits in error
            round_by_error(value,error);
            // Highest decimal part from value
            int decimal_max;
	    if (value == 0 && error == 0)
	      decimal_max = 0;
	    else
	      decimal_max = floor(log10(fabs(  (value==0 ? error : value )  )));
            //
            VecStr result;
            result.push_back(engineering_base(value,decimal_max,false)); // value
            result.push_back(engineering_base(error,decimal_max,false)); // error
            result.push_back(engineering_exponent(decimal_max,unit)); // exponent
            return result;
        }

        String  pm () {
            SStream out;
            if (useUnicode){
                out << ' ' << "\u00b1" << ' ' ;
            } else {
                out << " +- ";
            }
            return out.str();
        }

        VecStr time_format (double time_in_sec, bool useLong=false) {
            double time = time_in_sec;
            string time_unit = "s";
            string time_unit_long = " seconds";
            if      ( time_in_sec < 1e-3  ) { time *= 1e6   ; time_unit="us"; time_unit_long=" microseconds" ; }
            else if ( time_in_sec < 1.    ) { time *= 1000. ; time_unit="ms"; time_unit_long=" miliseconds"  ; }
            else if ( time_in_sec < 60.   ) { time /= 1.    ; time_unit="s"; time_unit_long=" seconds"      ; }
            else if ( time_in_sec < 3600. ) { time /= 60.   ; time_unit="mi"; time_unit_long=" minutes"      ; }
            else                            { time /= 3600. ; time_unit="hr"; time_unit_long=" hours"        ; }
            SStream out;
            VecStr result;
            out << setprecision(3) << time;
            result.push_back(out.str());
            if (useLong) result.push_back(time_unit_long);
            else result.push_back(time_unit);
            return result;
        }

        void Result(const Term &term, bool printGrandTotal) {
            double time = ( printGrandTotal ) ? term.total_time  : term.last_time ;
            double val  = ( printGrandTotal ) ? term.total_int   : term.last_int[0] ;
            double err  = ( printGrandTotal ) ? term.total_err2  : term.last_err2   ;
            err = sqrt(err);
            VecStr result = result_format(val,err);
            // integration
            cout << setw(wr) <<  result[0];
            cout << pm();
            cout << setw(we) <<  result[1];
            cout << setw(wx) <<  result[2];
            // time
            result = time_format(time);
            cout << " (";
            cout << setw(wt) << result[0] + result[1];
            cout << ")";
            EndOfCell();
        };

        void ResultSubTotal(bool is_grandtotal) {
            // cumulate time to subtotal
            subtotal.last_time   = clock_real()-subtotal.last_time;
            subtotal.total_time  += subtotal.last_time;
            Result(subtotal,is_grandtotal);
            if (opts.texttable) {
                if (is_grandtotal){
                    textfile  << subtotal.total_int << " ";
                    textfile  << sqrt(subtotal.total_err2);
                } else{
                    for (double val : subtotal.last_int) textfile << val << " ";
                    textfile  << sqrt(subtotal.last_err2);
                }
                textfile  << "\n";
            }
            EndOfRow();
            subtotal.last_reset();
        };


        void ResultGrandTotal() {
            Hline();
            Bounds(true);
            for( TermIterator term; !term.IsEnd(); ++term){
                Result((*term),true);
            }
            ResultSubTotal(true);
            Hline('b');
        };

        void SettingsHeader(String title){
            size_t t_length =title.size();
            size_t head_width = 72;
            String symb = useUnicode ? "\u254d" : "=";
            cout << endl;
            if (t_length == 0 ) {
                coutN(head_width, symb);
                cout << endl;
            }
            else {
                t_length+=2;
                int length = (head_width-t_length)/2;
                int add = (head_width-t_length)%2;
                coutN(length+add,symb);
                cout << " " << title << " ";
                coutN(length,symb);
            }
            cout << endl;
        }

        void ProcessSettings(){
            SettingsHeader("Process settings");
            // Hadron collision
            int den =  density_.ih1_ * density_.ih2_  ;
            if (den<0) cout << "   proton-antiproton collisions";
            if (den>0) cout << "   proton-proton collisions";
            cout << " at sqrt(s) = " << energy_.sroot_ << " GeV " << endl;
            // Order
            if (opts.fixedorder) {
                if (nnlo_.order_ == 0)
                    cout << "   Computing LO fixed order cross section" << endl;
                else if (nnlo_.order_ == 1)
                    cout << "   Computing NLO fixed order cross section" << endl;
                else if (nnlo_.order_ == 2)
                    cout << "   Computing NNLO fixed order cross section" << endl;
                else if (nnlo_.order_ == 2)
                    cout << "   Computing NNNLO fixed order cross section" << endl;
            } else {
                if (nnlo_.order_ == 0)
                    cout << "   Computing LO+LL resummed cross section" << endl;
                else if (nnlo_.order_ == 1)
                    cout << "   Computing NLO+NLL resummed cross section" << endl;
                else if (nnlo_.order_ == 2)
                    cout << "   Computing NNLO+NNLL resummed cross section" << endl;
                else if (nnlo_.order_ == 3)
                    cout << "   Computing NNNLO+NNNLL resummed cross section" << endl;
            }
            // Process
            if (nproc_.nproc_ == 1)
                cout << "   W+ ->  nu(p3)+e+(p4) production " << endl;
            else if(nproc_.nproc_ == 2)
                cout << "   W- ->  e^-(p3)+nu~(p4) production " << endl;
            else if(nproc_.nproc_== 3 && !opts.useGamma)
                cout << "   Z ->  e-(p3)+e+(p4) production " << endl;
            else if(nproc_.nproc_== 3 && opts.useGamma)
                cout << "   Z/gamma* ->  e-(p3)+e+(p4) production " << endl;
        }

        void CKMmatrix(){
            if (nproc_.nproc_ == 1 || nproc_.nproc_ == 2)
            {
                SettingsHeader("CKM matrix");
                cout << "     | Vud Vus Vub |   | " << setw(9) << cabib_.Vud_ << setw(9) <<  cabib_.Vus_ << setw(9)  << cabib_.Vub_ << " | " << endl;
                cout << "     | Vcd Vcs Vcb | = | " << setw(9) << cabib_.Vcd_ << setw(9) <<  cabib_.Vcs_ << setw(9)  << cabib_.Vcb_ << " | " << endl;
                cout << "     | Vtd Vts Vtb |   | " << setw(9) <<          0  << setw(9) <<            0 << setw(9)  <<           0 << " | " << endl;
            }
        }

        void EWparameters(){
            SettingsHeader("EW parameters");
            cout << Col4( "" , "Gf ="          , ewcouple_.Gf_     , "GeV^-2" , 12);
            cout << Col4( "" , "1/alpha ="     , 1/coupling::aemmz , ""       , 12);
            cout << Col4( "" , "sin(thW)^2 ="  , ewcouple_.xw_     , ""       , 12);
            cout << Col4( "" , "W mass ="      , dymasses_.wmass_  , "GeV"    , 12);
            cout << Col4( "" , "Z mass ="      , dymasses_.zmass_  , "GeV"    , 12);
            cout << Col4( "" , "W width ="     , dymasses_.wwidth_ , "GeV"    , 12);
            cout << Col4( "" , "Z width ="     , dymasses_.zwidth_ , "GeV"    , 12);
        }

        void QCDsettings(){
            SettingsHeader("QCD settings");
            cout << Col4(  "Renormalization scale:" , "mur ="    , scale_.scale_       , "GeV" );
            cout << Col4(  "Factorization scale:"   , "muf ="    , facscale_.facscale_ , "GeV" );
            cout << Col4(  "Resummation scale:"     , "m_ll/Q =" , a_param_.a_param_   ,  ""   );

	    //            cout << Col4(  "Renormalization scale:" , String(opts.fmuren    ? "dynamic" : "fixed") + " muren =" , opts.kmuren , String("* ") + (opts.dynamicscale ? "m_ll" : to_string(opts.rmass)) );
	    //            cout << Col4(  "Factorization scale:"   , String(opts.fmufac    ? "dynamic" : "fixed") + " mufac =" , opts.kmufac , String("* ") + (opts.dynamicscale ? "m_ll" : to_string(opts.rmass)) );
	    //            cout << Col4(  "Resummation scale:"     , String(opts.fmures    ? "dynamic" : "fixed") + " mures =" , opts.kmures , String("* ") + (opts.dynamicscale ? "m_ll" : to_string(opts.rmass)) );

            cout << Col4(  "Renormalization scale:" , String(opts.fmuren    ? "dynamic" : "fixed") + " muren =" , opts.kmuren , String("* ") + scales::func(opts.fmuren) );
            cout << Col4(  "Factorization scale:"   , String(opts.fmufac    ? "dynamic" : "fixed") + " mufac =" , opts.kmufac , String("* ") + scales::func(opts.fmufac) );
            cout << Col4(  "Resummation scale:"     , String(opts.fmures    ? "dynamic" : "fixed") + " mures =" , opts.kmures , String("* ") + scales::func(opts.fmures) );
	    
            cout << Col4( "Strong coupling"  , "alpha_s(MZ) ="   , couple_.amz_                 , ""      );
            cout << Col4( ""                 , "alpha_s(mur) ="  , qcdcouple_.as_               , ""      );
            cout << Col4( ""                 , "running order =" , (LHAPDF::getOrderAlphaS()+1) , "-loop" );
            cout << Col4( "Non-perturbative" , "form factor g =" , g_param_.g_param_            , "GeV^2" );
        }

        void PDFsettings(){
            SettingsHeader("PDF settings");
            cout << Col4( "" , "PDF set:"    , opts.LHAPDFset                      , "" );
            cout << Col4( "" , "PDF member:" , opts.LHAPDFmember                   , "" );
            cout << Col4( "" , "PDF errors:" , (opts.PDFerrors ? "true" : "false") , "" );
            if (opts.doBORN && !opts.fixedorder)
            {
                if (opts_.approxpdf_ == 1) cout << Col4( "Mellin moments of PDFs:" , "" , "Polynomial interpolation" , "" );
                else
                {
                    cout << Col4( "Mellin moments of PDFs:" , ""                              , "Numerical integration" , "" );
                    cout << Col4( ""                        , "nodes for the quadrature:"     , 64                      , "" );
                    cout << Col4( ""                        , "Intervals for the quadrature:" , opts_.pdfintervals_     , "" );
                }
            }
        }

        void QTrecoil(){
            SettingsHeader("qT recoil");
            if (opts.qtrec_cs)         cout << Col4(  "qt-recoil prescription:" , "Collins-Soper" , " kt1 = p_V(x)/2; kt2 = p_V(y)/2"   , "" );
            else if (opts.qtrec_kt0)   cout << Col4(  "qt-recoil prescription:" , "kt0"           , " kt1 = kt2 = 0"                    , "" );
            else if (opts.qtrec_naive) cout << Col4(  "qt-recoil prescription:" , "Naive"         , " kt = kt_CS * (1+p_V(z)/(m+E_V))"  , "" );
        }

        void AppliedCuts(){
            SettingsHeader("Cuts");
            cout << Col4( "", "apply cuts:" ,   (opts.makecuts ? "true" : "false") , "" );
            if (opts.makecuts){
                cout << Col4("" , "pt_l > "         , opts.lptcut  , "" );
                cout << Col4("" , "|eta_l| < "      , opts.lycut   , "" );
                cout << Col4("" , "pt_lep > "       , opts.lepptcut , "" );
                cout << Col4("" , "|eta_lep| < "    , opts.lepycut  , "" );
                cout << Col4("" , "pt_alep > "      , opts.alpptcut , "" );
                cout << Col4("" , "|eta_alep| < "   , opts.alpycut  , "" );
                cout << Col4("" , "user cuts: "     , "true"       , "" );
            }
        }

        void ResummationDamp(){
            if (!opts.fixedorder && ((opts.doBORN && !opts.fixedorder) || opts.doCT)) {
                SettingsHeader("ResummationDamp");
		cout << Col4( "" , "damping :"           , opts.damp                 , ""    );
                cout << Col4( "" , "damping mode:"       , opts.dampmode             , ""    );
                cout << Col4( "" , "damp above:"         , opts.dampk*opts.rmass     , "GeV" );
                cout << Col4( "" , "damping width:"      , opts.dampdelta*opts.rmass , "GeV" );
                cout << Col4( "" , "resummation cutoff:" , opts.qtcutoff*1000        , "MeV" );
            }
        }

        void DebugSettings(){
            SettingsHeader("Debug settings");
            cout << Col4( "" , "timeprofile:"   , (opts.timeprofile ? "true" : "false") , "" );
            cout << Col4( "" , "verbose:"       , (opts.verbose     ? "true" : "false") , "" );
            cout << Col4( "" , "cubaverbosity:" , opts.cubaverbosity                    , "" );
        }


        void Settings(){
            ProcessSettings();
            CKMmatrix();
            EWparameters();
            QCDsettings();
            PDFsettings();
            QTrecoil();
            AppliedCuts();
            ResummationDamp();
            DebugSettings();
            SettingsHeader(); // empty header will create line
        }

        void Footer() {
            SettingsHeader("Summary");
            VecStr s_res = result_format(subtotal.total_int, sqrt(subtotal.total_err2), true);
            cout << Col4 ("Total cross section" , "" , s_res[0]+pm()+s_res[1] , s_res[2]  );
            s_res = time_format(subtotal.total_time, true);
            cout << Col4 ("Calculation time"    , "" , s_res[0] , s_res[1] );
            cout << Col4 ("Output file"         , "" , HistoHandler::result_filename + HistoHandler::file_suffix , "");
            if (opts.texttable){
                String name = HistoHandler::result_filename + ".txt";
                cout << Col4 ("Result in text file" , "" , name , "");
                std::ofstream f (name.c_str(), std::ofstream::trunc);
                f << textfile.str();
                f.close();
            }
            SettingsHeader();
        }
    }
}
#endif /* ifndef print_table_C */
