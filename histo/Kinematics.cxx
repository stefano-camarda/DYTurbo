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
#include "KinematicCuts.h"

namespace Kinematics{

    // input values definitions
    double p3[4];
    double p4[4];
    double event_weight;

    // All flags are stored here
    ObservableFlags flags;

    void SetKinematics(double _p3[4], double _p4[4], double wgt=1.){
        isIntegratorMode=false;
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
            LepAbsRap  absy1;
            ALpCh      ch2;
            ALpPT      pt2;
            ALpAbsRap  absy2;
            BosMT      mt;
            MET        etmiss;
	    //absolute-rapidity-ordered leptons
  	    LepCPT lcpt;
	    LepFPT lfpt;
	    LepCAbsEta lcy;
	    LepFAbsEta lfy;
            // define here cut decision
            bool operator()(){
	      
                if (opts.lptcut > 0){
                    if (ch1()!=0 && pt1() < opts.lptcut) return SkipEvent;
                    if (ch2()!=0 && pt2() < opts.lptcut) return SkipEvent;
                }
                if (opts.lycut < 100){
                    if (ch1()!=0 && absy1() > opts.lycut) return SkipEvent;
                    if (ch2()!=0 && absy2() > opts.lycut) return SkipEvent;
                }
                if (opts.mtcut     >0   && mt()     < opts.mtcut     ) return SkipEvent;
                if (opts.etmisscut >0   && etmiss() < opts.etmisscut ) return SkipEvent;
                if (opts.lepptcut   >0   && pt1()    < opts.lepptcut   ) return SkipEvent;
                if (opts.alpptcut   >0   && pt2()    < opts.alpptcut   ) return SkipEvent;
                if (opts.lepycut    <100 && absy1()  > opts.lepycut    ) return SkipEvent;
                if (opts.alpycut    <100 && absy2()  > opts.alpycut    ) return SkipEvent;
                if ((opts.lcymax    <100 || opts.lcymin    > 0) && (lcy()  > opts.lcymax || lcy()  < opts.lcymin)    ) return SkipEvent;
                if ((opts.lfymax    <100 || opts.lfymin    > 0) && (lfy()  > opts.lfymax || lfy()  < opts.lfymin)    ) return SkipEvent;
                if (opts.lcptcut   >0   && lcpt()    < opts.lcptcut   ) return SkipEvent;
                if (opts.lfptcut   >0   && lfpt()    < opts.lfptcut   ) return SkipEvent;

		return KeepEvent;
            }
        } standard_cuts;

        UserCuts user_cuts;

        bool KeepThisEvent(double p3[4], double p4[4]){
            SetKinematics(p3,p4);
            if (opts.makecuts){
                if ( standard_cuts() == SkipEvent ) return SkipEvent ;
                if ( user_cuts() == SkipEvent ) return SkipEvent ;
            }
            return KeepEvent;
        }
    }
}


#endif /* ifndef Kinematics_CXX */
