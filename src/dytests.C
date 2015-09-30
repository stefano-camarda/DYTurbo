#ifndef DYTESTS
#define DYTESTS

/**
 * @file dytests.C
 * Description of this file
 *
 * @brief A brief description
 *
 * @author cuto <Jakub.Cuth@cern.ch>
 * @date 2015-09-30
 */

#include "interface.h"
#include "settings.h"
#include <cmath>
#include <cstdio>


void test_CentralLeptonCut(){
    // force cuts on
    opts.makelepcuts=true;
    nproc_.nproc_=3;
    // particle kinematics
    double p[4][12];
    int njets = 0;
    // testing point
    for ( double pt3 : { 10., 30. }) for ( double pt4 : { 10., 30. }) for ( double y3 : { 1., 5. }) for ( double y4 : { 1., 5. }){
        // set lep kinem
        double e2y3= exp(-2*y3);
        double e2y4= exp(-2*y4);
        double sqY3 = pow ( (1-e2y3)/(1+e2y3), 2);
        double sqY4 = pow ( (1-e2y4)/(1+e2y4), 2);
        double pz3sq = pt3*pt3 * sqY3/(1-sqY3);
        double pz4sq = pt4*pt4 * sqY4/(1-sqY4);
        //
        p[0][2]= pt3;
        p[1][2]= 0;
        p[2][2]= sqrt(pz3sq);
        p[3][2]= sqrt(pt3*pt3 + pz3sq);
        //
        p[0][3]= pt4;
        p[1][3]= 0;
        p[2][3]= sqrt(pz4sq);
        p[3][3]= sqrt(pt4*pt4 + pz4sq);
        // use cuts
        bool cutC = cuts_    (p,njets);
        bool cutF = cutsold_ (p,njets);
        if (cutC!=cutF) printf (" NOT SAME decission cutC %d cutF %d  pt3 %g pt4 %g y3 %g y4 %g \n", cutC, cutF, pt3, pt4, y3, y4);
        else printf (" same decission cutC %d cutF %d  pt3 %g pt4 %g y3 %g y4 %g \n", cutC, cutF, pt3, pt4, y3, y4);
    }
}

/**
 * Description of main program
 *
 */
int main(int argc, const char * argv[]){

    test_CentralLeptonCut();

    return 0;
}


#endif // DYTEST
