#ifndef TurboHist_H2_H
#define TurboHist_H2_H
/**
 * @file TurboHist_H2.h
 * Description of cpp file
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-09-07
 */

#include "TurboHist_HBase.h"
#include "TurboHist_Counter.h"

namespace TurboHist {
    struct H2 : public HBase<H2,Counter> { };
}

#endif /* ifndef TurboHist_H2_H */
