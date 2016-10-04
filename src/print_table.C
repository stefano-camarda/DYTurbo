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
#include "interface.h"
#include "coupling.h"
#include "histo/HistoHandler.h"

#include <LHAPDF/LHAPDF.h>

#include<iostream>
#include<iomanip>

using std::setw;
using std::flush;
using std::cout;
using std::endl;

namespace DYTurbo {
    namespace PrintTable {

        size_t eoc=2; //!< end of cell length
        size_t wb=5; //!< bound length
        size_t wr=8; //!< result length
        size_t wt=5; //!< time length
        size_t bound_width = 2*wb+3;
        size_t time_width = 2+wt+2;
        size_t term_width = 2*wr+4 + time_width;

        inline void EndOfCell()  { cout << " |" << flush; }
        inline void BeginOfRow() { cout << "|"  << flush; }
        inline void EndOfRow()   { cout << endl;          }

        void Hline(){
            size_t N_loopingBounds =0;
            for (auto &a: ActiveBoundaries ) if (a.size()>2) N_loopingBounds++;
            size_t wtotal = 1; // begin of line
            wtotal += (bound_width +eoc)*N_loopingBounds;
            wtotal += (term_width  +eoc)*(ActiveTerms.size() + 1/*total column*/ );
            String hline ( wtotal, '-');
            cout << hline << endl;
        };

        void BoundsAllLooping(bool printNames=false, bool useFullBound=false){
            // loop over bounds
            for (size_t ibound=0; ibound < N_boundaries; ++ibound){
                if (ActiveBoundaries[ibound].size()<3) continue;
                if (printNames){
                    cout << setw(bound_width) << ActiveBoundaries[ibound].name ;
                } else {
                    double lo = useFullBound ? ActiveBoundaries[ibound].front() : last_bounds.loBound(ibound);
                    double hi = useFullBound ? ActiveBoundaries[ibound].back()  : last_bounds.hiBound(ibound);
                    cout << setw(wb) << lo ;
                    cout << " - ";
                    cout << setw(wb) << hi ;
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
            Hline();
            BeginOfRow();
            BoundsAllLooping(true); // print only names
            for( TermIterator term; !term.IsEnd(); ++term){
                TermName(*term);
            }
            TermName(subtotal);
            EndOfRow();
            Hline();
        };

        void Footer() {
            cout << "Total cross-section:        "; Result(subtotal,true); cout << endl;
            cout << "Result was written in file: " << HistoHandler::result_filename << HistoHandler::file_suffix <<endl;
        };


        void Bounds(bool use_full_bound) {
            BeginOfRow();
            BoundsAllLooping(false,use_full_bound);
        };

        vector<string> Round(double value, double error=0, bool sign=false) {
            vector <string> result;
            int decimal = 0;
            //If no error, value is rounded to two significant digits
            if (error == 0) error = value;
            if (error != 0) decimal = -log10(fabs(error)) + 2; //-> this 2 sets the number of digits
            decimal = max(0, decimal);

            char Dec[10];
            sprintf (Dec, "%d", decimal);
            string D = Dec;

            char Numb[50];
            // write number
            if(sign) sprintf (Numb, ((string)"%+." + D + "f").c_str(), value);
            else     sprintf (Numb, ((string)"%."  + D + "f").c_str(), value);
            result.push_back(Numb);
            // write error
            sprintf (Numb, ((string)"%." + D + "f").c_str(), error);
            result.push_back(Numb);

            return result;
        }

        void Result(const Term &term, bool printGrandTotal) {
            /// @todo set correct width of table cell
            double time = ( printGrandTotal ) ? term.total_time  : term.last_time ;
            double val  = ( printGrandTotal ) ? term.total_int  : term.last_int[0] ;
            double err  = ( printGrandTotal ) ? term.total_err2 : term.last_err2   ;
            err = sqrt(err);
            VecStr result = Round(val,err);
            cout << setw(wr) <<  result[0];
            cout << " +- ";
            cout << setw(wr) <<  result[1];
            //
            cout << " (" << setprecision(3) << setw(wt) << time << "s)";
            cout << setprecision(6);
            EndOfCell();
        };

        void ResultSubTotal(bool is_grandtotal) {
            Result(subtotal,is_grandtotal);
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
            Hline();
        };

        void SettingsHeader(String title){
            size_t t_length =title.size();
            size_t head_width = 70;
            char symb = '=';
            cout << endl;
            if (t_length == 0 ) cout << String(head_width, symb) <<endl;
            else {
                String line ( int((head_width-t_length)/2), symb); 
                cout << line << ((head_width-t_length)%2 ? "=" : "") ;
                cout << " " << title << " ";
                cout << line << endl;
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
            } else {
                if (nnlo_.order_ == 1)
                    cout << "   Computing NLO+NLL resummed cross section" << endl;
                else if (nnlo_.order_ == 2)
                    cout << "   Computing NNLO+NNLL resummed cross section" << endl;
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
            cout << Col4( "" , "Narrow width:" , (opts.zerowidth ? "true" : "false")    , ""       , 16);
        }

        void QCDsettings(){
            SettingsHeader("QCD settings");
            cout << Col4(  "Renormalization scale:" , "mur ="    , scale_.scale_       , "GeV" );
            cout << Col4(  "Factorization scale:"   , "muf ="    , facscale_.facscale_ , "GeV" );
            cout << Col4(  "Resummation scale:"     , "m_ll/Q =" , a_param_.a_param_   ,  ""   );

            cout << Col4(  "Renormalization scale:" , String(opts.dynamicscale    ? "dynamic" : "fixed") + " muren =" , opts.kmuren , String("* ") + (opts.dynamicscale ? "m_ll" : to_string(opts.rmass)) );
            cout << Col4(  "Factorization scale:"   , String(opts.dynamicscale    ? "dynamic" : "fixed") + " mufac =" , opts.kmufac , String("* ") + (opts.dynamicscale ? "m_ll" : to_string(opts.rmass)) );
            cout << Col4(  "Resummation scale:"     , String(opts.dynamicresscale ? "dynamic" : "fixed") + " mures =" , opts.kmures , String("* ") + (opts.dynamicscale ? "m_ll" : to_string(opts.rmass)) );

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
            SettingsHeader("QCD settings");
            cout << Col4( "", "apply cuts:" ,   (opts.makecuts ? "true" : "false") , "" );
            if (opts.makecuts){
                cout << Col4("" , "pt_l > "         , opts.lptcut  , "" );
                cout << Col4("" , "|eta_l| < "      , opts.lycut   , "" );
                cout << Col4("" , "pt_l(1st) > "    , opts.l1ptcut , "" );
                cout << Col4("" , "|eta_l|(1st) < " , opts.l1ycut  , "" );
                cout << Col4("" , "pt_l(2nd) > "    , opts.l1ptcut , "" );
                cout << Col4("" , "|eta_l|(2nd) < " , opts.l1ycut  , "" );
                cout << Col4("" , "user cuts: "     , "true"       , "" );
            }
        }

        void ResummationDamp(){
            if (!opts.fixedorder && ((opts.doBORN && !opts.fixedorder) || opts.doCT)) {
                SettingsHeader("ResummationDamp");
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
    }
}
#endif /* ifndef print_table_C */
