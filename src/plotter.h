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

        // term definition
        enum TermType { Resum, CT, LO, Real, Virt, VV, VJ, Total, None=999 };
        // kinematic and angular variables
        double Q2,Q,qt,y,a[NMOM],c[NMOM];
        double costh,phi,phi_lep;

        // Interface
        void Init();
        bool IsInitialized();
        void FillEvent(double p3[4], double p4[4], double wgt); ///<Normal filling of histograms.
        void FillRealDipole(double p3[4], double p4[4], double wgt,int nd); ///<Collect dipole kinematics and weights. Fill ai profiles.
        void FillRealEvent(TermType term = None); ///< for real filling without correlations. Need to FillRealDipole before.
        void FillQuadrature(double int_val, double int_error); ///< Adding integration to proper bin
        void FillResult(TermType term, double int_val, double int_error, double time); ///< Fill the result of integration.
        void CumulateResult(TermType term, double wgt); ///< Fill the result of integration.
        void SetPDF(int npdf); ///< Set histograms for pdf memeber npdf.
        void Merge();
        void Dump();
        void Finalise(double xsection=0);


    protected :
#ifdef USEROOT
        // Data members
        //
        // Histogramming
        /// @todo: use one object instead, be more variable
        // boson kinematic histograms
        TH1D * h_qt;
        TH1D * h_y ;
        TH2D * h_qtVy;
        TH1D * h_Q;
        // angular histograms
        bool doAiMoments;
        struct AiProf2D { TProfile2D* A[NMOM]; };
        struct AiProf   { TProfile*   A[NMOM]; };
        AiProf2D pa_qtVy;
        AiProf pa_qt;
        AiProf pa_y;
        TH1D * h_costh;
        TH1D * h_phi;
        TH1D * h_phi_lep;
        // final results
        TH2D* qt_y_resum ;
        TH2D* qt_y_ct    ;
        TH2D* qt_y_lo    ;
        TH2D* qt_y_real  ;
        TH2D* qt_y_virt  ;
        TH2D* qt_y_vv    ;
        TH2D* qt_y_vj    ;
        TH2D* qt_y_total ;
        // dipole helper variables
        struct XsecPoint {
            int ibin;
            double qt,y,costh,phi,phi_lep,wgt;
            bool fid;
            double Q,A[NMOM];
        } point;
        std::vector<XsecPoint> dipole_points;

        // PDF hists
        int last_npdf;
        std::vector<TH1D     *> h_qt_PDF      ;
        std::vector<TH1D     *> h_y_PDF       ;
        std::vector<TH2D     *> h_qtVy_PDF    ;
        std::vector<TH1D     *> h_Q_PDF       ;
        std::vector<TH1D     *> h_costh_PDF   ;
        std::vector<TH1D     *> h_phi_PDF     ;
        std::vector<TH1D     *> h_phi_lep_PDF ;
        std::vector<AiProf2D >  pa_qtVy_PDF   ;
        std::vector<AiProf   >  pa_qt_PDF     ;
        std::vector<AiProf   >  pa_y_PDF      ;

        // helpers
        double N;
        std::vector<double> v_wgt;
        double verbose;
        //AiMoments ai_maarten;

        // Function members
        //
        // helper function
        TH1 * clone_PDF( TH1 *h, int npdf);
        template<typename T,typename P>
            void clone_Array_PDF( std::vector<T> &v_ha, int npdf);
        void addToBin(TH1* h, double int_val, double int_err);
        // calculation of kinematics and angular variables
        void calculate_kinematics(double p3[4], double p4[4]);
        double calcQt(double p[4]);
        double calcY(double p[4]);
        double calcQ2(double p[4]);
        double Vplus(double p[4]);
        double Vminus(double p[4]);
        double calcCosThCS(double Q2,double qt,double p3[4],double p4[4]);
        double calcPhiCS(double p3[4],double p4[4],double &phi_lep);
        // Debug printing
        void print_dipole(XsecPoint pt);
        void print_dipoleVec(std::vector<XsecPoint> vec );


#endif // USEROOT

};

#ifdef USEROOT
template<typename T,typename P>
void plotter::clone_Array_PDF( std::vector<T> &v_ha, int npdf){
    v_ha.push_back(T());
    for(int i=0;i<NMOM;i++) {
            v_ha.at(npdf).A[i] = (P *) clone_PDF(v_ha.at(0).A[i],npdf);
    } 
}
#endif // USEROOT


extern plotter hists;


#endif //plotter_h

