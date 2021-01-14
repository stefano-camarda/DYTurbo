#ifndef gint_h
#define gint_h

#include <complex.h>

using namespace std;

namespace gint
{
  void init();
  void intg(complex <double> q, complex<double> jac, double blim);
  void intbeta(complex <double> q, complex <double> jac, double blim);
  void intexpc(complex <double> q, complex <double> jac, double blim);
  void calc(complex <double> b);
  void allocate();
  void reset();
  void free();

  extern complex <double> logS;
  
  extern complex <double> logasl;
  extern complex <double> logasl_pdf;
  extern complex <double> logasl_expc;
  extern complex <double> *alogqq;
  extern complex <double> *alogqg;
  extern complex <double> *alogqqb;
  extern complex <double> *alogqqp;
  extern complex <double> *alogqqbp;
#pragma omp threadprivate(logS,logasl,logasl_pdf,logasl_expc,alogqq,alogqg,alogqqb,alogqqp,alogqqbp)
}

#endif
