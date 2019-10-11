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
#include "settings.h"

 /**
  * @brief Main program.
  *
  * This is using the DYTurbo namespace to run your calculation.
  */
int main(int argc, char *argv[])
{
  DYTurbo::Init(argc,argv);              //Init, read config file
  DYTurbo::WarmUp();                     //Set up integration terms and bins boundaries
  if (!opts.silent) DYTurbo::PrintTable::Header();

    //Loop on phase space bins
    for ( DYTurbo::BoundIterator bounds; !bounds.IsEnd(); ++bounds)
      {
        DYTurbo::SetBounds(bounds);
	if (!opts.silent) DYTurbo::PrintTable::Bounds();

	//Loop on active terms
        for (DYTurbo::TermIterator term; !term.IsEnd(); ++term)
	  {
            (*term).RunIntegration();
	    if (!opts.silent) DYTurbo::PrintTable::Result((*term));
	  }
	
	if (!opts.silent) DYTurbo::PrintTable::ResultSubTotal();
      }
    
    if (!opts.silent) DYTurbo::PrintTable::ResultGrandTotal();
    DYTurbo::Terminate();
    return 0;
}

#endif /* ifndef main_C */
