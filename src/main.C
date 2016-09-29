#ifndef main_C
#define main_C
/**
 * @file main.C
 * @brief Main program of DYTURBO
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-09-29
 */

#include "dyturbo.h"

 /**
  * @brief Main program.
  *
  * This is using the DYTurbo namespace to run your calculation.
  */
int main(int argc, char *argv[])
{
    DYTurbo::Init(argc,argv);
    DYTurbo::WarmUp();
    DYTurbo::PrintTable::Header();
    for ( DYTurbo::BoundIterator bounds; !bounds.IsEnd(); ++bounds) {
        DYTurbo::SetBounds(bounds);
        DYTurbo::PrintTable::Bounds();
        for (DYTurbo::TermIterator term; !term.IsEnd(); ++term){
            (*term).RunIntegration();
            DYTurbo::PrintTable::Result((*term));
        }
        DYTurbo::PrintTable::ResultSubTotal();
    }
    DYTurbo::PrintTable::ResultGrandTotal();
    DYTurbo::Terminate();
    return 0;
}

#endif /* ifndef main_C */
