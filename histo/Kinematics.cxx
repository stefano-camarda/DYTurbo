#ifndef Kinematics_CXX
#define Kinematics_CXX
/**
 * @file Kinematics.cxx
 * Description of this macro
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-08-26
 */

#include "Kinematics.h"
#include "KinematicDefinitions.h"
#include "user/user_cuts.h"

namespace Kinematics{

    // input values definitions
    double p3[4];
    double p4[4];
    double event_weight;

    // All flags are stored here
    ObservableFlags flags;

    void SetKinematics(double _p3[4], double _p4[4], double wgt=1.){
        for (auto i=0; i<4;i++){
            p3[i]=_p3[i];
            p4[i]=_p4[i];
        }
        event_weight=wgt;
        // Set all variables to re-calculate
        for (auto ptr_isCalculated : flags) *ptr_isCalculated=false;
    }

    void SetMiddlePoint(){
        isIntegratorMode=true;
        for (auto ptr_isCalculated : flags) *ptr_isCalculated=false;
    }

    bool isIntegratorMode=false;

    namespace Cuts {

        const bool KeepEvent=true;
        const bool SkipEvent=false;

        struct StandardCuts : public CutBase {
            // define here all observables and variables
            LepCh      ch1;
            LepPT      pt1;
            LepAbsEta  aeta1;
            AlpCh      ch2;
            AlpPT      pt2;
            AlpAbsEta  aeta2;
            BosMT      mt;
            MET        etmiss;
            // define here cut decision
            bool operator()(){
                if (opts.lptcut > 0){
                    if (ch1()!=0 && pt1() < opts.lptcut) return SkipEvent;
                    if (ch2()!=0 && pt2() < opts.lptcut) return SkipEvent;
                }
                if (opts.lycut < 100){
                    if (ch1()!=0 && aeta1() > opts.lycut) return SkipEvent;
                    if (ch2()!=0 && aeta2() > opts.lycut) return SkipEvent;
                }
                if (opts.mtcut     >0 && mt()     < opts.mtcut     ) return SkipEvent;
                if (opts.etmisscut >0 && etmiss() < opts.etmisscut ) return SkipEvent;
                if (opts.l1ptcut   >0 && pt1()    < opts.l1ptcut   ) return SkipEvent;
                if (opts.l2ptcut   >0 && pt2()    < opts.l2ptcut   ) return SkipEvent;
                if (opts.l1ycut    >0 && aeta1()  > opts.l1ycut    ) return SkipEvent;
                if (opts.l2ycut    >0 && aeta2()  > opts.l2ycut    ) return SkipEvent;
                return KeepEvent;
            }
        } standard_cuts;

        bool KeepThisEvent(double p3[4], double p4[4]){
            SetKinematics(p3,p4);
            if (opts.makecuts){
                if ( standard_cuts() == SkipEvent ) return SkipEvent ;
                if ( user_cut() == SkipEvent ) return SkipEvent ;
            }
            return KeepEvent;
        }
    }
}


#endif /* ifndef Kinematics_CXX */
