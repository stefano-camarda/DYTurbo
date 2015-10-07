#ifndef plotter_h
#define plotter_h

#include "config.h"

#ifdef USEROOT
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile2D.h>
#endif // USEROOT

#include<memory>
#include<mutex>

class plotter {
    public :
        plotter();
        ~plotter();

        enum TermType { Resum, CT, LO, Real, Virt, Total };

        void Init();
        void FillEvent(double p3[4], double p4[4], double wgt); ///<Normal filling of histograms.
        void FillRealDipole(double p3[4], double p4[4], double wgt,int nd); ///<Collect dipole kinematics and weights. Fill ai profiles.
        void FillRealEvent(); ///< for real filling without correlations. Need to FillRealDipole before.
        void FillResult(TermType term, double int_val, double int_error, double time); ///< Fill the result of integration.
        void Merge();
        void Dump();
        void Finalise(double xsection=0);


    private :
#ifdef USEROOT
        void CalculateKinematics(double p3[4], double p4[4]);

        /// shared space
        std::shared_ptr<int> sh_N;
        std::mutex m;

        std::vector<double> v_wgt;

        /// @todo: use one object instead
        double N;
        // kinematic histograms
        //TH1D * h_l1_pt;
        TH1D * h_qt;
        TH1D * h_y ;
        TH2D * h_qtVy;
        // profiles
        TProfile2D * p_qtVy_A4;
        // final results
        TH2D* qt_y_resum ;
        TH2D* qt_y_ct    ;
        TH2D* qt_y_lo    ;
        TH2D* qt_y_real  ;
        TH2D* qt_y_virt  ;
        TH2D* qt_y_total ;
        // kinematics
        double qt,y,a4;
        // dipole variables
        struct XsecPoint {
            int ibin;
            double qt,y,wgt;
        } point;
        std::vector<XsecPoint> dipole_points;
        void print_dipole(XsecPoint pt);
        void print_dipoleVec(std::vector<XsecPoint> vec );

#endif // USEROOT

};

extern plotter hists;


#endif //plotter_h

