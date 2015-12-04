#ifndef plotter_h
#define plotter_h

#include "config.h"

#ifdef USEROOT
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile2D.h>
#endif // USEROOT

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
        TH2D * h_qtVy;
        // profiles
        bool doAiMoments;
        TProfile2D * p_qtVy_A[8];
        // final results
        TH2D* qt_y_resum ;
        TH2D* qt_y_ct    ;
        TH2D* qt_y_lo    ;
        TH2D* qt_y_real  ;
        TH2D* qt_y_virt  ;
        TH2D* qt_y_total ;
        // kinematics
        double Q2,qt,y,a[8];
        double costh,phi;
        // dipole variables
        struct XsecPoint {
            int ibin;
            double qt,y,wgt;
            bool fid;
        } point;
        std::vector<XsecPoint> dipole_points;
        void print_dipole(XsecPoint pt);
        void print_dipoleVec(std::vector<XsecPoint> vec );

        // PDF hists
        int last_npdf;
        std::vector<TH1D *> h_qt_PDF   ;
        std::vector<TH1D *> h_y_PDF    ;
        std::vector<TH2D *> h_qtVy_PDF ;
        TH1 * clone_PDF( TH1 *h, int npdf);

        double verbose;
#endif // USEROOT

};

extern plotter hists;


#endif //plotter_h

