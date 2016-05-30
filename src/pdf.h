#ifndef pdf_h
#define pdf_h

#include "mcfm_interface.h"

extern void setalphas();
extern void setg();

extern "C" {
  void fdist_(int& ih, double& x, double& xmu, double fx[2*MAXNF+1]);
  void dysetpdf_(int& member);
  void setmellinpdf_(int& member);
}



#endif
