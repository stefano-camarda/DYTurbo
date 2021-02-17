#ifndef evolnative_h
#define evolnative_h

#include "mellinint.h"

#include <complex>
using namespace std;

namespace evolnative
{
  extern int ndim;
  
  //PDFs mellin moments at the factorisation scale (only positive branch is needed as PDFs at the starting scale are real valued)
  extern complex <double> *UVP;
  extern complex <double> *DVP;
  extern complex <double> *SVP;
  extern complex <double> *CVP;
  extern complex <double> *BVP;
  extern complex <double> *USP;
  extern complex <double> *DSP;
  extern complex <double> *SSP;
  extern complex <double> *CHP;
  extern complex <double> *BOP;
  extern complex <double> *GLP;
#pragma omp threadprivate(UVP,DVP,SVP,CVP,BVP,USP,DSP,SSP,CHP,BOP,GLP)
  
  //Singlet/non-singlet decomposition
  extern complex <double> *SIP;
  extern complex <double> *NS3P;
  extern complex <double> *NS8P;
  extern complex <double> *NS15P;
  extern complex <double> *NS24P;
  extern complex <double> *NS35P;
#pragma omp threadprivate(SIP,NS3P,NS8P,NS15P,NS24P,NS35P)

  //alphas for the evolution
  extern complex <double> alpqf;
  extern complex <double> alpq;
  extern complex <double> alpr;
#pragma omp threadprivate(alpqf,alpq,alpr)

  //scales for the PDF evolution used in evolnative::evolve()
  extern complex <double> XL;
  extern complex <double> XL1;
  extern complex <double> SALP;
#pragma omp threadprivate(XL,XL1,SALP)

  //anomalous dimensions
  extern complex <double> *ans,*am,*ap,*al,*be,*ab;
  extern complex <double> *rmin,*rplus,*rqq,*rqg,*rgq,*rgg;
  extern complex <double> *rmmqq,*rmmqg,*rmmgq,*rmmgg,*rmpqq,*rmpqg,*rmpgq,*rmpgg,*rpmqq,*rpmqg,*rpmgq,*rpmgg,*rppqq,*rppqg,*rppgq,*rppgg;
#pragma omp threadprivate(ans,am,ap,al,be,ab,rmin,rplus,rqq,rqg,rgq,rgg,rmmqq,rmmqg,rmmgq,rmmgg,rmpqq,rmpqg,rmpgq,rmpgg,rpmqq,rpmqg,rpmgq,rpmgg,rppqq,rppqg,rppgq,rppgg)
  
  extern void init();
  extern void allocate(); //allocate dynamic memory
  extern void allocate_pdfs();
  extern void allocate_engine();
  extern void update();
  extern void update_pdfs();
  extern void update_engine();
  extern void scales();
  extern void evolve();
  extern void free();     //free dynamic memory
  extern void free_pdfs();
  extern void free_engine();
  extern void release();

  inline int index(int i, int beam, int sign)  {return i + mellinint::mdim*((beam-1) + 2*sign);};
  inline int index(int i, int sign)            {return i + ndim*sign;};
  
  //obsolete
  extern void init_fortran();
}

#endif
