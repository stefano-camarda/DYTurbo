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

  extern void (*xfxq)(const double &x, const double &Q, double *fPDF);
  extern void lhaxfxq(const double &x, const double &Q, double *fPDF);

  extern double (*alphas)(const double &Q);
  extern double lhaalphas(const double &Q);
}

extern "C" {
  void fdist_(int& ih, double& x, double& xmu, double fx[2*MAXNF+1]);
  void dysetpdf_(int& member);
  void setmellinpdf_(int& member);
}



#endif
