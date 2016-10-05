#ifndef interface_C
#define interface_C
/**
 * @file interface.C
 * @brief A brief description
 *
 * Detailed description of file
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-10-05
 */

#include "interface.h"
#include "histo/KinematicCuts.h"
#include "histo/HistoHandler.h"


#include <cassert>


int cuts_(double p[4][12], int &njet){
    double p3[4];
    double p4[4];
    for (int i=0; i<4; i++){
        p3[i] = p[i][2];
        p4[i] = p[i][3];
    }
    // MCFM expects opposite logic false=accept event
    return !Kinematics::Cuts::KeepThisEvent(p3,p4);
}


void hists_setpdf_(int * npdf){
    HistoHandler::SetVariation(*npdf);
}

void hists_fill_(double p3[4], double p4[4], double *weight){
    HistoHandler::FillEvent(p3,p4,*weight);
}

void hists_fill_pdf_(double p3[4], double p4[4], double *weight, int *npdf){
    HistoHandler::SetVariation(*npdf);
    HistoHandler::FillEvent(p3,p4,*weight);
}

void hists_real_dipole_(double p3[4], double p4[4], double *weight, int * nd){
    HistoHandler::FillDipole(p3,p4,*weight);
}

void hists_real_dipole_pdf_(double p3[4], double p4[4], double *weight, int * nd, int* npdf){
    HistoHandler::SetVariation(*npdf);
    HistoHandler::FillDipole(p3,p4,*weight);
}

void hists_real_event_(){
    HistoHandler::FillRealEvent();
}

void disabled_hists_real_event_pdf_(int *npdf){
    /// @attention This is very dangerous !!!
    HistoHandler::SetVariation(*npdf);
    HistoHandler::FillRealEvent();
}

void disabled_hists_finalize_(){
    HistoHandler::Terminate();
}

#endif /* ifndef interface_C */
