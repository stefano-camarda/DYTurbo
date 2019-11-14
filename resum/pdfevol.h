#ifndef pdfevol_h
#define pdfevol_h

#include "interface.h"

#include <complex>
using namespace std;

extern "C"
{
  void pdfmoments_(int &beam, double &scale, fcomplex &N, fcomplex &UV, fcomplex &DV, fcomplex &US, fcomplex &DS, fcomplex &SP, fcomplex &SM, fcomplex &GL, fcomplex &CH, fcomplex &BO, double &xmin);

  //access dyres PDF in N-space
  extern struct {
    //    fcomplex cfx1_[136][11];
    //    fcomplex cfx2p_[136][11];
    //    fcomplex cfx2m_[136][11];
    fcomplex cfx1_[512][11];
    fcomplex cfx2p_[512][11];
    fcomplex cfx2m_[512][11];
  } creno_;
}
#pragma omp threadprivate(creno_)


namespace pdfevol
{
  //PDFs mellin moments at the factorisation scale
  extern complex <double> *UVP;
  extern complex <double> *DVP;
  extern complex <double> *USP;
  extern complex <double> *DSP;
  extern complex <double> *SSP;
  extern complex <double> *GLP;
  extern complex <double> *CHP;
  extern complex <double> *BOP;
#pragma omp threadprivate(UVP,DVP,USP,DSP,SSP,GLP,CHP,BOP)

  //Singlet/non-singlet decomposition
  extern complex <double> *SIP;
  extern complex <double> *NS3P;
  extern complex <double> *NS8P;
  extern complex <double> *NS15P;
  extern complex <double> *NS24P;
  extern complex <double> *NS35P;
#pragma omp threadprivate(SIP,NS3P,NS8P,NS15P,NS24P,NS35P)
  
  //evolved PDFs mellin moments
  extern complex <double> *fx1;
  extern complex <double> *fx2;
#pragma omp threadprivate(fx1,fx2)
    
  //moments for the PDFs convolution
  extern complex <double> fn1[2*MAXNF+1];
  extern complex <double> fn2[2*MAXNF+1];
#pragma omp threadprivate(fn1,fn2)

  //scales
  extern complex <double> bscale;
  extern complex <double> bstarscale;
  extern complex <double> bstartilde;
  extern complex <double> qbstar;
  extern complex <double> bcomplex;
  extern complex <double> XL;
  extern complex <double> XL1;
  extern complex <double> SALP;
  extern complex <double> alpr;
#pragma omp threadprivate(bscale,bstarscale,bstartilde,qbstar,bcomplex,XL,XL1,SALP,alpr)

  extern void allocate(); //allocate dynamic memory
  extern void free();     //free dynamic memory

  extern void allocate_fx(); //allocate dynamic memory
  extern void free_fx();     //free dynamic memory
  
  //initialise and compute Mellin moments of PDFs at the starting scale (factorisation scale)
  extern void init();
  //release memory if fmufac = 0
  extern void release();
  //update Mellin moments of PDFs at the starting scale with a dynamic factorisation scale
  extern void update();
  //evolve Mellin moments of PDFs from the factorisation scale to the scale ~1/b, set in bscale
  extern void evolution();
  //calculate Mellin moments of PDFs at all scales by direct Mellin transform
  extern void calculate(int i);
  extern void calculate();
  //store the moments in the fx1 and fx2 arrays
  extern void storemoments(int i, complex <double> fx[11]);
  //store the moments in the Fortran common block
  extern void storemoments_fortran(int i, complex <double> fx[11]);
  //retrieve moments corresponding to i1, i2, and sign from the fx1 and fx2 arrays and store them in fn1 and fn2
  extern void retrieve(int i1, int i2, int sign);
  //retrieve moments corresponding to i1, i2, and sign from the Fortran common block and store them in fn1 and fn2
  extern void retrieve_fortran(int i1, int i2, int sign);
  //retrieve only fn1
  extern void retrieve_beam1(int i1);
  //retrieve only fn2 positive
  extern void retrieve_beam2_pos(int i2);
  //retrieve only fn2 negative
  extern void retrieve_beam2_neg();

  //flavour dependent form factors
  extern void flavour_kt();
  
  //retrieve moments corresponding to i, and positive branch from fx1 and fx2 and store them in fn1 and fn2, to be used with the mellin1d option
  extern void retrieve1d_pos(int i);
  //retrieve moments corresponding to i, and negative branch by complex conjugation of the positive branch moments into fn1 fn2, to be used with the mellin1d option
  extern void retrieve1d_neg();
  //retrieve moments corresponding to i, and sign from the Fortran common block and store them in fn1 and fn2, to be used with the mellin1d option
  extern void retrieve1d_fortran(int i, int sign);

  extern void retrievemuf(int i);

  extern void truncate();
  extern void uppertruncate();
  
}

#endif
