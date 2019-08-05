#ifndef pdf_h
#define pdf_h

#include <LHAPDF/LHAPDF.h>

#include "mcfm_interface.h"

namespace pdf
{
  extern LHAPDF::PDF* lhapdf;

  extern void init();
  extern void setalphas();
  extern void setg();
}

extern "C" {
  void fdist_(int& ih, double& x, double& xmu, double fx[2*MAXNF+1]);
  void dysetpdf_(int& member);
  void setmellinpdf_(int& member);
}



#endif
