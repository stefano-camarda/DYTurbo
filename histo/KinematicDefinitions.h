#ifndef KinematicDefinitions_H
#define  KinematicDefinitions_H
/**
 * @file KinematicDefinitions.h
 * Description of cpp file
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-08-29
 */

#include "Kinematics.h"

#define NEWKIN(CLS) class CLS : public Variable< CLS >

#include <cmath>
#include <algorithm>
using std::max;

// forward phasespace
namespace phasespace {
    extern double mmin, mmax, qtmin, qtmax, ymin, ymax;
}


namespace Kinematics{
    // NOTE: Writing twice the class name is intentional (CTRP)! Check
    // `Kinematic.h` for more info. Or use preprocessor macro NEWKIN.
    // Description of posibilities can be found in `CosThCS`

    // TODO: Decide which lepton
    // lepton
    NEWKIN( LepPX ) { double calc(){ return p3[0]; } };
    NEWKIN( LepPY ) { double calc(){ return p3[1]; } };
    NEWKIN( LepPZ ) { double calc(){ return p3[2]; } };
    NEWKIN( LepE  ) { double calc(){ return p3[3]; } };

    // antilepton
    NEWKIN( ALpPX ) { double calc(){ return p4[0]; } };
    NEWKIN( ALpPY ) { double calc(){ return p4[1]; } };
    NEWKIN( ALpPZ ) { double calc(){ return p4[2]; } };
    NEWKIN( ALpE  ) { double calc(){ return p4[3]; } };

    // boson
    class BosPX : public Variable< BosPX > { double calc(){ return p3[0]+p4[0]; } };
    class BosPY : public Variable< BosPY > { double calc(){ return p3[1]+p4[1]; } };
    class BosPZ : public Variable< BosPZ > { double calc(){ return p3[2]+p4[2]; } };
    class BosE  : public Variable< BosE  > { double calc(){ return p3[3]+p4[3]; } };


    class BosM : public Variable< BosM > {
        BosPX px;
        BosPY py;
        BosPZ pz;
        BosE  e;
        // Calc is 
        double calc(){
            return sqrt( e()*e() - (px()*px() + py()*py() + pz()*pz()) )  ;
        }
        double middlePoint(){
            return ( phasespace::mmax + phasespace::mmin )/2. ;
        }
    };
    template<> const bool Variable<BosM>::isIntegratorVariable=true;

    class BosMT : public Variable< BosMT > {
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

    class BosPT : public Variable< BosPT > {
        BosPX px;
        BosPY py;
        double calc(){
            return sqrt(px()*px()+py()*py())  ;
        }
        double middlePoint(){
            return ( phasespace::qtmax + phasespace::qtmin )/2. ; 
        }
    };
    template<> const bool Variable<BosPT>::isIntegratorVariable=true;

    NEWKIN( BosPhi ) {
        BosPX px;
        BosPY py;
        double calc(){
            return atan2(py(),px())  ;
        }
    };

    class BosY : public Variable<BosY> {
        BosPZ pz;
        BosE  e;
        double calc(){
            return 0.5*log((e()+pz())/(e()-pz()));
        }
        double middlePoint(){
            return ( phasespace::ymax + phasespace::ymin )/2. ;
        }
    };
    template<> const bool Variable<BosY>::isIntegratorVariable=true;

    // Longitudinal angle theta in Collin-Soper frame
    NEWKIN(CosThCS){
        // You can define additional functions for this class
        double Vplus  (double p[4]) { return (p[3]+p[2]); };
        double Vminus (double p[4]) { return (p[3]-p[2]); };
        // Pointers to global arrays stays updated
        double * lm = p3;
        double * lp = p4;
        // Variable objects updated automatically. Keeping as class memebers is
        // less time consuming.
        BosPZ pz;
        BosM m;
        BosPT pt;
        // Mandatory function `calc`: is called when kinematics is updated
        double calc(){
            double costh=0;
            costh = (Vplus(lm)*Vminus(lp) - Vplus(lp)*Vminus(lm));
            costh /= sqrt(m()*m() * (m()*m() + pt()*pt()));
            costh *=  pz() < 0. ? -1 : 1; //sign flip according to boson rapidity
            return costh;
        }
    };

    // Azimuthal angle phi in Collin-Soper frame
    NEWKIN(PhiCS){
        double * lm = p3;
        double * lp = p4;
        BosPX px;
        BosPY py;
        BosPZ pz;
        BosM m;
        BosPT pt;
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
                plxCS = 0.5 * m() / sqrt(m()*m()+pt()*pt()) * (2.*plx - pt());
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

                plxCS = 0.5 * m()/sqrt(m()*m()+pt()*pt()) * (delta[0]*qht[0]+delta[1]*qht[1]);
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
}

#include "user_kinem.h"

#endif /* ifndef KinematicDefinitions_H */