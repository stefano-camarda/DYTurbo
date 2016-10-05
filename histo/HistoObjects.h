#ifndef HistoObjects_H
#define HistoObjects_H

/**
 * @file HistoObjects.h
 * Description of this macro
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-08-26
 */

#include "HistoHandler.h"
#include "HistoBase.h"

#include "config.h" // need to be include config to decide whether to use ROOT

#ifdef USEROOT
#include "HistoObjectsROOT.h"
#define HistoImp HistoROOT
namespace HistoHandler{
    // c++11: template<typename T> using HistoImp = HistoROOT<T>;
    typedef TH1D        H1;
    typedef TH2D        H2;
    typedef TH3D        H3;
    typedef TProfile    P1;
    typedef TProfile2D  P2;
}

#else // STL
#include "HistoObjectsSTL.h"
#define HistoImp HistoSTL
namespace HistoHandler{
    typedef TurboHist::H1 H1;
    typedef TurboHist::H2 H2;
    typedef TurboHist::H3 H3;
    typedef TurboHist::P1 P1;
    typedef TurboHist::P2 P2;
}
#endif

#include "HistoSpecialization.h"
#include "HistoBook.h" // move this to HistoBook cxx and compile



#endif /* ifndef HistoObjects_H */
