#ifndef KinematicDefinitions_H
#define  KinematicDefinitions_H
/**
 * @file KinematicDefinitions.h
 * @brief Definition of all Observables inside Kinematics namespace.
 * Only developers should edit this file. User-defined observables should be in
 * `user/user_kinem.h`
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>, Stefano Camarda <Stefano.Camarda@cern.ch>
 * @date 2016-08-29
 */

#include "Kinematics.h"
#include "KinUtils.h"

#define NEWKIN(CLS) class CLS : public Observable< CLS >

#include <cmath>
#include <algorithm>
using std::max;

// forward phasespace
#include "phasespace/phasespace.h"
#include "src/settings.h"

namespace Kinematics{



    /// @todo Decide which lepton is charged. Should be done. Please write test.
    /// @defgroup LepVar Lepton related observables.
    /// @{
    NEWKIN( LepPX  ) { double calc(){ return p3[0]; } };
    NEWKIN( LepPY  ) { double calc(){ return p3[1]; } };
    NEWKIN( LepPZ  ) { double calc(){ return p3[2]; } };
    NEWKIN( LepE   ) { double calc(){ return p3[3]; } };
    NEWKIN( LepPT  ) { double calc(){ return Util::pT(p3[0],p3[1]); } };
    NEWKIN( LepCh  ) { double calc(){ return opts.nproc==1 ? 0 : -1 ; } };
    NEWKIN( LepEta ) { double calc(){ return Util::eta(p3[0],p3[1],p3[2]); } };
    NEWKIN( LepAbsEta ) { LepEta eta; double calc(){ return fabs(eta()); } };
    NEWKIN( LepRap ) { double calc(){ return Util::rap(p3[3],p3[2]); } };
    NEWKIN( LepAbsRap ) { LepRap rap; double calc(){ return fabs(rap()); } };
    /// @}

    /// @defgroup AlpVar Anti-lepton related observables.
    /// @{
    NEWKIN( ALpPX  ) { double calc(){ return p4[0]; } };
    NEWKIN( ALpPY  ) { double calc(){ return p4[1]; } };
    NEWKIN( ALpPZ  ) { double calc(){ return p4[2]; } };
    NEWKIN( ALpE   ) { double calc(){ return p4[3]; } };
    NEWKIN( ALpPT  ) { double calc(){ return Util::pT(p4[0],p4[1]); } };
    NEWKIN( ALpCh  ) { double calc(){ return opts.nproc==2 ? 0 : +1 ; } };
    NEWKIN( ALpEta ) { double calc(){ return Util::eta(p4[0],p4[1],p4[2]); } };
    NEWKIN( ALpAbsEta ) { ALpEta eta; double calc(){ return fabs(eta()); } };
    NEWKIN( ALpRap ) { double calc(){ return Util::rap(p4[3],p4[2]); } };
    NEWKIN( ALpAbsRap ) { ALpRap rap; double calc(){ return fabs(rap()); } };
    /// @}


    /// @defgroup BosVar Vector boson related observables.
    /// @{
    NEWKIN( BosPX ) { double calc(){ return p3[0]+p4[0]; } };
    NEWKIN( BosPY ) { double calc(){ return p3[1]+p4[1]; } };
    NEWKIN( BosPZ ) { double calc(){ return p3[2]+p4[2]; } };
    NEWKIN( BosE  ) { double calc(){ return p3[3]+p4[3]; } };

    class BosM2 : public Observable< BosM2 > {
        BosPX px;
        BosPY py;
        BosPZ pz;
        BosE  e;
        double calc(){
            return Util::mass2(px(),py(),pz(),e());
        }
        double middlePoint(){
            double m = ( phasespace::mmax + phasespace::mmin )/2. ;
            return m*m;
        }
        public :
        inline bool IsIntegrableObservable() const {return true;}
    };

    class BosM : public Observable< BosM > {
        BosM2 m2;
        double calc(){
            return sqrt(m2());
        }
        double middlePoint(){
            return ( phasespace::mmax + phasespace::mmin )/2. ;
        }
        public :
        inline bool IsIntegrableObservable() const {return true;}
    };

    /**
     * @brief Boson transverse mass (experimentalist definition)
     *
     * Implemetation follows definition:
     *
     * \f$m_T \equiv \sqrt{ 2 p_{T,\ell} p_{T,\nu} (1-\cos(\phi) }\f$
     *
     * For saving calculation time cos is calculated as dot product.
     */
    class BosMT : public Observable< BosMT > {
        LepPX lmPX;
        LepPY lmPY;
        ALpPX lpPX;
        ALpPY lpPY;

        double calc(){
            double mtrans = 1;
            mtrans*=sqrt(lmPX()*lmPX() + lmPY()*lmPY());
            mtrans*=sqrt(lpPX()*lpPX() + lpPY()*lpPY());
            mtrans-=    (lmPX()*lpPX() + lmPY()*lpPY());
            mtrans*=2;
            return sqrt(max(mtrans,0.));
        }
    };

    class BosPT2 : public Observable< BosPT2 > {
        BosPX px;
        BosPY py;
        double calc(){
            return Util::pT2(px(),py())  ;
        }
        double middlePoint(){
            double qt = ( phasespace::qtmax + phasespace::qtmin )/2. ; 
            return  qt*qt;
        }
        public :
        inline bool IsIntegrableObservable() const {return true;}
    };

    class BosPT : public Observable< BosPT > {
        BosPT2 qt2;
        double calc(){
            return sqrt(qt2());
        }
        double middlePoint(){
            return ( phasespace::qtmax + phasespace::qtmin )/2. ; 
        }
        public :
        inline bool IsIntegrableObservable() const {return true;}
    };

    NEWKIN( BosPhi ) {
        BosPX px;
        BosPY py;
        double calc(){
            return atan2(py(),px())  ;
        }
    };

    class BosY : public Observable<BosY> {
        BosPZ pz;
        BosE  e;
        double calc(){
            return 0.5*log((e()+pz())/(e()-pz()));
        }
        double middlePoint(){
            return ( phasespace::ymax + phasespace::ymin )/2. ;
        }
        public :
        inline bool IsIntegrableObservable() const {return true;}
    };
    /// @}

    NEWKIN ( MET ) {
        LepCh ch1;
        LepPT pt1;
        ALpCh ch2;
        ALpPT pt2;
        double calc(){
            if (ch1()) return pt1();
            if (ch2()) return pt2();
            return 0;
        }
    };

    /// @defgroup AngulVar Final state kinematic observables.
    /// @{

    /**
     * @brief Cosine of longitudinal angle \f$\theta\f$ in Collin-Soper frame.
     *
     * This class was choosen as example implementation of Observable.
     *
     * @note Writing twice the class name is intentional. Check \ref
     * Observable description for more info. Or use preprocessor macro NEWKIN.
     */
    class CosThCS : public Observable< CosThCS > {
        //! You can define additional functions within this class.
        double Vplus  (double p[4]) { return (p[3]+p[2]); };
        double Vminus (double p[4]) { return (p[3]-p[2]); };
        //! Pointers to global arrays -- stays updated.
        double * lm = p3;
        double * lp = p4;
        /**
         * @brief Observable objects are updated automatically.
         *
         * Keeping them as class memebers is less time consuming.
         */
        BosPZ pz;
        BosM2 m2;
        BosPT2 pt2;
        /**
         * Mandatory function `calc`: is called when kinematics is updated.
         *
         * This where actual calculation of observable goes.
         *
         * This definition is consistent with Collins-Sopper-Sterman definition 
         */
        double calc(){
            double costh=0;
            costh = (Vplus(lm)*Vminus(lp) - Vplus(lp)*Vminus(lm));
            costh /= sqrt(m2() * (m2() + pt2()));
            costh *=  pz() < 0. ? -1 : 1; //sign flip according to boson rapidity
            return costh;
        }
        /**
         * Optional reimplementation of `middlePoint` for Integrable Observables. 
         * 
         * @note Dont forget also to set `IsIntegrableObservable` to true !
         */
        double middlePoint(){
            return ( phasespace::cthmax + phasespace::cthmin )/2. ;
        }
        //! Optional reimplementation of `IsIntegrableObservable` for Integrable Observables. Otherwise false.
        public :
        inline bool IsIntegrableObservable() const {return true;}
    };

    //! Azimuthal angle phi in Collin-Soper frame.
    NEWKIN(PhiCS){
        double * lm = p3;
        double * lp = p4;
        BosPX px;
        BosPY py;
        BosPZ pz;
        BosM2 m2;
        BosM m;
        BosPT pt;
        BosPT2 pt2;
        double calc(){
            double plxCS, plyCS;
            if (pt() == 0) //if qt = 0, use the original x and y axis to determine phiCS in the boson rest frame
            {
                plxCS = lm[0];
                plyCS = lm[1];
            }
            else // if qt > 0, use the boson direction in the transverse plane as x axis in the boson rest frame
            {
                /*******************************************************************/
                //Mirkes definition as in Nucl.Phys. B387 (1992) 385 Eq.(22)
                //first needs to rotate the lab frame so that the x axis lies in the event plane defined by the boson and proton directions
                double c = px()/pt();//cos(-phiZ);
                double s = -py()/pt();//sin(-phiZ);
                double plx = c*lm[0] - s*lm[1];
                double ply = s*lm[0] + c*lm[1];

                //Now apply formulas (22) of Nucl.Phys. B387 (1992) 385
                plxCS = 0.5 * m() / sqrt(m2()+pt2()) * (2.*plx - pt());
                plyCS = ply;

                /*******************************************************************/
                //Original formula, as in Eq.(2.1) of Phys.Rev. D16 (1977) 2219 
                double delta[2];
                delta[0] = lm[0]-lp[0];
                delta[1] = lm[1]-lp[1];

                //unit vector in the transverse plane from the projection of the boson momentum
                double qht[2];
                qht[0] = px()/pt();
                qht[1] = py()/pt();

                //unit vector in the transverse plane perpendicular to the boson momentum and beam 1
                double rht[2];
                rht[0] = -py()/pt();
                rht[1] = px()/pt();

                plxCS = 0.5 * m()/sqrt(m2()+pt2()) * (delta[0]*qht[0]+delta[1]*qht[1]);
                plyCS = 0.5 * (delta[0]*rht[0]+delta[1]*rht[1]);
                /*******************************************************************/
            }
            double sign = pz() > 0. ? 1 : -1; //sign flip according to boson rapidity
            return atan2(-sign*plyCS,-plxCS); //rotate by M_PI (just a convention, it means that the x axis is opposite to the boson direction in the transverse plane)
        }
    };

    //Theta goniometrics
    NEWKIN( SinThCS      ){ CosThCS  costh; /*       */     double calc(){ return sqrt(max(0.,1.-costh()*costh()));                } };
    NEWKIN( Sin2ThCS     ){ CosThCS  costh; SinThCS  sinth; double calc(){ return 2*costh()*sinth();                               } };
    //Phi   goniometrics
    NEWKIN( CosPhiCS     ){ PhiCS    phi;   /*       */     double calc(){ return cos(phi());                                      } };
    NEWKIN( Cos2PhiCS    ){ CosPhiCS cosph; /*       */     double calc(){ return 2*cosph()*cosph()-1.;                            } };
    NEWKIN( SinPhiCS     ){ CosPhiCS cosph; PhiCS    phi;   double calc(){ return sqrt(max(0.,1.-cosph()*cosph()))*(phi()>0?1:-1); } };
    NEWKIN( Sin2PhiCS    ){ CosPhiCS cosph; SinPhiCS sinph; double calc(){ return 2*cosph()*sinph();                               } };
    //Angular polynomials
    NEWKIN( A0 ){ CosThCS  costh;  /*        */      double calc(){ return 20./3. * (0.5-1.5*costh()*costh()  ) +2./3.; } };
    NEWKIN( A1 ){ Sin2ThCS sin2th; CosPhiCS  cosph;  double calc(){ return 5.     * (sin2th()*cosph()         ) ;       } };
    NEWKIN( A2 ){ SinThCS  sinth;  Cos2PhiCS cos2ph; double calc(){ return 10.    * (sinth()*sinth()*cos2ph() ) ;       } };
    NEWKIN( A3 ){ SinThCS  sinth;  CosPhiCS  cosph;  double calc(){ return 4.     * (sinth()*cosph()          ) ;       } };
    NEWKIN( A4 ){ CosThCS  costh;  /*        */      double calc(){ return 4.     * (costh()                  ) ;       } };
    NEWKIN( A5 ){ SinThCS  sinth;  Sin2PhiCS sin2ph; double calc(){ return 5.     * (sinth()*sinth()*sin2ph() ) ;       } };
    NEWKIN( A6 ){ Sin2ThCS sin2th; SinPhiCS  sinph;  double calc(){ return 4.     * (sin2th()*sinph()         ) ;       } };
    NEWKIN( A7 ){ SinThCS  sinth;  SinPhiCS  sinph;  double calc(){ return 4.     * (sinth()*sinph()          ) ;       } };
    /// @}
}

#include "user/user_kinem.h"

#endif /* ifndef KinematicDefinitions_H */
