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

// nasty hack -- to get acces to privite stuff
#define private public
#define protected public

#include "interface.h"
#include "settings.h"
#include <cmath>
#include <cstdio>

#include <TLorentzVector.h>
#include <TRandom.h>
#include "plotter.h"


void test_CentralLeptonCut(){
    // force cuts on
    opts.makelepcuts=true;
    nproc_.nproc_=3;
    opts.nproc=nproc_.nproc_;
    // particle kinematics
    double p[4][12];
    int njets = 0;
    // testing point
    for (int i : {1,2,3,4,5} ){
        opts.fiducial = static_cast<settings::DetFiducial> (i);
        printf ("fiducial %d\n",i);
        for ( double pt3 : { 10., 22., 30. }) for ( double pt4 : { 10., 22.,  30. }) for ( double y3 : { 0.5, 1.1, 1.3, 2.2 }) for ( double y4 : { 0.5, 1.1, 1.3, 2.2 }){
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
            double pvb[4];
            for ( int j : {0,1,2,3}) pvb[j] = p[j][2]+p[j][3];
            double mt = sqrt( pvb[3]*pvb[3] - pvb[2]*pvb[2]);
            const char * throw_away = (cuts_(p,njets) ? "THROW " : "KEEP  ");
            printf ("%s :   pt3 %g pt4 %g y3 %g y4 %g mt %g \n", throw_away, pt3, pt4, y3, y4, mt);
            // use cuts
            //bool cutC = cuts_    (p,njets);
            //bool cutF = cutsold_ (p,njets);
            //if (cutC!=cutF) printf (" NOT SAME decission cutC %d cutF %d  pt3 %g pt4 %g y3 %g y4 %g \n", cutC, cutF, pt3, pt4, y3, y4);
            //else printf (" same decission cutC %d cutF %d  pt3 %g pt4 %g y3 %g y4 %g \n", cutC, cutF, pt3, pt4, y3, y4);
        }
    }
}


  
void getCSFAngles(const TLorentzVector & lep1, const int &charge1, const TLorentzVector & lep2, double ebeam, double &costh, double &phi, std::vector<double> *aimom=0)
{
    TLorentzVector boson = lep1+lep2;
    double Lplus  = (charge1 < 0) ? lep1.E()+lep1.Pz() : lep2.E()+lep2.Pz();
    double Lminus = (charge1 < 0) ? lep1.E()-lep1.Pz() : lep2.E()-lep2.Pz();
    double Pplus  = (charge1 < 0) ? lep2.E()+lep2.Pz() : lep1.E()+lep1.Pz();
    double Pminus = (charge1 < 0) ? lep2.E()-lep2.Pz() : lep1.E()-lep1.Pz();

    costh  = (Lplus*Pminus - Lminus*Pplus);
    costh *= TMath::Abs(boson.Pz());
    costh /= (boson.Mag()*boson.Pz());
    costh /= TMath::Sqrt(boson.Mag2() + boson.Pt()*boson.Pt());

    TVector3 boostV = -boson.BoostVector();
    TLorentzVector lep1_boosted = (charge1 < 0 ) ? lep1 : lep2;
    lep1_boosted.Boost(boostV);

    TVector3 CSAxis, xAxis, yAxis;
    TLorentzVector p1, p2;
    double sign = +1.;
    if (boson.Z() < 0)
        sign = -1.;
    p1.SetXYZM(0., 0., sign*ebeam, 0.938); // quark (?)
    p2.SetXYZM(0., 0., -sign*ebeam, 0.938); // antiquark (?)

    p1.Boost(boostV);
    p2.Boost(boostV);
    CSAxis = (p1.Vect().Unit()-p2.Vect().Unit()).Unit();
    yAxis = (p1.Vect().Unit()).Cross(p2.Vect().Unit());
    yAxis = yAxis.Unit();
    xAxis = yAxis.Cross(CSAxis);
    xAxis = xAxis.Unit();

    phi = TMath::ATan2(lep1_boosted.Vect()*yAxis, lep1_boosted.Vect()*xAxis);
}
void getCSFAngles(double p3[4], double p4[4], double &costh, double &phi){
    TLorentzVector lep1(p3);
    TLorentzVector lep2(p4);
    const int ch1 = -1;
    const int ch2 = 1;
    double ebeam = opts.sroot/2.;
    getCSFAngles(lep1,ch1,lep2, ebeam, costh,phi,0 );
}

void test_CalculationCollinsSopper(){
    opts.sroot=7e3;
    double costh_maarten, phi_maarten;
    double p3[4],p4[4];
    plotter hists;
    TRandom rand;

    for (int atry=0; atry<1000; atry++){
        double p3s = 0;
        double p4s = 0;
        for (int i=0; i<3; i++){
            //p3[i] = rand.Uniform(-100, 100);
            //p4[i] = rand.Uniform(-100, 100);
            p3[i] = rand.Gaus(0, 10);
            p4[i] = rand.Gaus(0, 10);
            p3s+=p3[i]*p3[i];
            p4s+=p4[i]*p4[i];
        }
        p3[3] = TMath::Sqrt(p3s);
        p4[3] = TMath::Sqrt(p4s);
        hists.CalculateKinematics(p3,p4);
        //
        TLorentzVector lep1(p3);
        TLorentzVector lep2(p4);
        TLorentzVector VB(lep1+lep2);
        const int ch1 = -1;
        const int ch2 = 1;
        double ebeam = opts.sroot/2.;
        getCSFAngles(lep1,ch1,lep2, ebeam, costh_maarten,phi_maarten,0 );
        if (TMath::Abs(costh_maarten) > 1.0) continue;
        double Q2 = VB.M2();
        double qt = VB.Pt();
        double y = VB.Rapidity();
        //getCSFAngles(p3,p4,costh_maarten,phi_maarten);

        //printf ( " p3 is %g %g %g %g\n",p3[0],p3[1],p3[2],p3[3]);
        //printf ( " p4 is %g %g %g %g\n",p4[0],p4[1],p4[2],p4[3]);
        //printf( " Q2 %g qt %g y %g z %g \n" , hists.Q2, hists.qt, hists.y, VB.Z() );
        printf("\n\n\nEvent\n");
        lep1.Print();
        lep2.Print();
        VB.Print();
        if (hists.phi   != phi_maarten   ) printf ("  PHI   is not same: me %g maarten %g\n",   hists.phi   , phi_maarten   );
        if (hists.costh != costh_maarten ) printf ("  COSTH is not same: me %g maarten %g\n",   hists.costh , costh_maarten );
        if (hists.Q2    != Q2            ) printf ("  Q2    is not same: me %g maarten %g\n", hists.Q2    , Q2            );
        if (hists.qt    != qt            ) printf ("  QT    is not same: me %g maarten %g\n", hists.qt    , qt            );
        if (hists.y     != y             ) printf ("  Y     is not same: me %g maarten %g\n", hists.y     , y             );
    }
}

/**
 * Description of main program
 *
 */
int main(int argc, const char * argv[]){

    //test_CentralLeptonCut();
    test_CalculationCollinsSopper();

    return 0;
}


#endif // DYTEST
