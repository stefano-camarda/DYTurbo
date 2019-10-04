#ifndef mellinpdf_h
#define mellinpdf_h

#include <complex>
using namespace std;

namespace mellinpdf
{
  void allocate();
  void free();

  void init(double xmin);
  void evalpdfs(double scale);
  void gauss_quad();
  void laguerre_ipol();

  extern complex <double> *UV;
  extern complex <double> *DV;
  extern complex <double> *US;
  extern complex <double> *DS;
  extern complex <double> *SP;
  extern complex <double> *SM;
  extern complex <double> *GL;
  extern complex <double> *CP;
  extern complex <double> *CM;
  extern complex <double> *BP;
  extern complex <double> *BM;
#pragma omp threadprivate(UV,DV,US,DS,SP,SM,GL,CP,CM,BP,BM)
  
  extern double *fuv;
  extern double *fdv;
  extern double *fus;
  extern double *fds;
  extern double *fsp;
  extern double *fsm;
  extern double *fgl;
  extern double *fcp;
  extern double *fcm;
  extern double *fbp;
  extern double *fbm;
#pragma omp threadprivate(fuv,fdv,fus,fds,fsp,fsm,fgl,fcp,fcm,fbp,fbm)
}
#endif
