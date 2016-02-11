#ifndef plotter_h
#define plotter_h

#include "config.h"
//#include "AiMoments.h"

#ifdef USEROOT
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TProfile2D.h>
#endif // USEROOT

#define NMOM 9

class plotter {
    public :
        plotter();
        ~plotter();

        enum TermType { Resum, CT, LO, Real, Virt, Total, None=999 };

        void Init();
        bool IsInitialized();
        void FillQuadrature(double int_val, double int_error); ///< Adding integration to proper bin
        void FillEvent(double p3[4], double p4[4], double wgt); ///<Normal filling of histograms.
        void FillRealDipole(double p3[4], double p4[4], double wgt,int nd); ///<Collect dipole kinematics and weights. Fill ai profiles.
        void FillRealEvent(TermType term = None); ///< for real filling without correlations. Need to FillRealDipole before.
        void FillResult(TermType term, double int_val, double int_error, double time); ///< Fill the result of integration.
        void SetPDF(int npdf); ///< Set histograms for pdf memeber npdf.
        void CumulateResult(TermType term, double wgt); ///< Fill the result of integration.
        void Merge();
        void Dump();
        void Finalise(double xsection=0);

        double Q2,qt,y,a[NMOM],c[NMOM];
        double costh,phi,phi_lep;

    protected :
#ifdef USEROOT
        void CalculateKinematics(double p3[4], double p4[4]);
        void addToBin(TH1* h, double int_val, double int_err);

        std::vector<double> v_wgt;

        /// @todo: use one object instead
        double N;
        // kinematic histograms
        TH1D * h_qt;
        TH1D * h_y ;
        TH1D * h_costh;
        TH1D * h_phi;
        TH1D * h_phi_lep;
        TH2D * h_qtVy;
        // profiles
        struct AiProf2D { TProfile2D* A[NMOM]; };
        struct AiProf   { TProfile*   A[NMOM]; };
        bool doAiMoments;
        //TProfile2D * p_qtVy_A[NMOM];
        //TProfile * p_qt_A[NMOM];
        //TProfile * p_y_A[NMOM];
        AiProf2D pa_qtVy;
        AiProf pa_qt;
        AiProf pa_y;
        // final results
        TH2D* qt_y_resum ;
        TH2D* qt_y_ct    ;
        TH2D* qt_y_lo    ;
        TH2D* qt_y_real  ;
        TH2D* qt_y_virt  ;
        TH2D* qt_y_total ;
        // kinematics
        //double Q2,qt,y,a[NMOM],c[NMOM];
        //double costh,phi;
        // dipole variables
        struct XsecPoint {
            int ibin;
            double qt,y,costh,phi,wgt;
            bool fid;
        } point;
        std::vector<XsecPoint> dipole_points;
        void print_dipole(XsecPoint pt);
        void print_dipoleVec(std::vector<XsecPoint> vec );

        // PDF hists
        int last_npdf;
        std::vector<TH1D     *> h_qt_PDF      ;
        std::vector<TH1D     *> h_y_PDF       ;
        std::vector<TH2D     *> h_qtVy_PDF    ;
        std::vector<TH1D     *> h_costh_PDF   ;
        std::vector<TH1D     *> h_phi_PDF     ;
        std::vector<TH1D     *> h_phi_lep_PDF ;
        std::vector<AiProf2D >  pa_qtVy_PDF   ;
        std::vector<AiProf   >  pa_qt_PDF     ;
        std::vector<AiProf   >  pa_y_PDF      ;
        //
        TH1 * clone_PDF( TH1 *h, int npdf);
        template<typename T>
            void clone_Array_PDF( std::vector<T> &v_ha, int npdf);

        //AiMoments ai_maarten;


        double verbose;
#endif // USEROOT

};

#ifdef USEROOT
template<typename T>
void plotter::clone_Array_PDF( std::vector<T> &v_ha, int npdf){
    v_ha.resize(npdf);
    for(int i=0;i<NMOM;i++) v_ha.at(npdf).A[i] = clone_PDF(v_ha.at(0).A[i],npdf);
}
#endif // USEROOT


extern plotter hists;


#endif //plotter_h

