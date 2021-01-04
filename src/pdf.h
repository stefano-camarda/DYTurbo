#ifndef pdf_h
#define pdf_h

#include "parton.h"

#include <LHAPDF/LHAPDF.h>

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
  extern double rgktalphas(double q);

  extern int order;
  extern double xmin;
  extern double qmin;
  extern double mc;
  extern double mb;
  extern double mt;
  extern double g;
}

extern "C" {
  void fdist_(int& ih, double& x, double& xmu, double fx[2*MAXNF+1]);
  void dysetpdf_(int& member);
  void setmellinpdf_(int& member);
  double dyalphas_mcfm_(double &q, double &amz, int &nloop);
  //  double dyalphas_lhapdf_(double &q); //no need for this, as there is a C++ LHAPDF function
}



#endif
