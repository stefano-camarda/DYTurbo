#ifndef pdfevol_h
#define pdfevol_h

#include "interface.h"

#include <complex>
using namespace std;

extern "C"
{
  void pdfmoments_(int &beam, double &scale, fcomplex &N, fcomplex &UV, fcomplex &DV, fcomplex &US, fcomplex &DS, fcomplex &SP, fcomplex &SM, fcomplex &GL, fcomplex &CH, fcomplex &BO, double &xmin);
  fcomplex alphasl_(fcomplex &q2);
  
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
  //evolved PDFs mellin moments
  extern complex <double> *fx1;
  extern complex <double> *fx2;
#pragma omp threadprivate(fx1,fx2)
    
  //moments for the PDFs convolution
  extern complex <double> fn1[2*MAXNF+1];
  extern complex <double> fn2[2*MAXNF+1];
#pragma omp threadprivate(fn1,fn2)

  //scales
  //extern complex <double> bb;
  extern complex <double> bstar;
  extern complex <double> mub_a;
  extern complex <double> mub;
  extern complex <double> mubstar_a;
  extern complex <double> mubstar;
  extern complex <double> mubstartilde;
#pragma omp threadprivate(bstar,mub_a,mub,mubstar_a,mubstar,mubstartilde)

  //alphas for the evolution
  extern complex <double> logasl, asl;
#pragma omp threadprivate(logasl, asl)
  
  //obsolete
  extern complex <double> bscale;
  extern complex <double> bstarscale;
  extern complex <double> bstartilde;
  extern complex <double> qbstar;
  extern complex <double> bcomplex;
#pragma omp threadprivate(bscale,bstarscale,bstartilde,qbstar,bcomplex)

  extern void allocate_fx(); //allocate dynamic memory
  extern void free_fx();     //free dynamic memory

  extern void allocate(); //allocate dynamic memory
  extern void free();     //free dynamic memory
  
  //initialise and compute Mellin moments of PDFs at the starting scale (factorisation scale)
  extern void init();
  extern void init_fortran();
  //release memory if fmufac = 0
  extern void release();
  //update Mellin moments of PDFs at the starting scale with a dynamic factorisation scale
  extern void update();
  //set various scales of order b0/b
  extern void scales(complex <double> b);
  //compute the LL,NLL,NNLL,NNNLL evolution of alphas
  extern void alphasl(complex <double> b);
  //evolution main selector
  extern void evolution();
  //evolve Mellin moments of PDFs from the factorisation scale to the scale ~b0/b, set in bscale
  extern void evolve();
  //calculate Mellin moments of PDFs at all scales by direct Mellin transform
  extern void calculate(int i);
  extern void calculate();
  //store the moments in the fx1 and fx2 arrays
  extern void storemoments(int i, complex <double> fx[11]);
  extern void storemoments_1(int i, complex <double> fx[11]);
  extern void storemoments_2(int i, complex <double> fx[11]);
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

  extern void retrievemuf_1d(int i);
  extern void retrievemuf_2d(int i);

  extern void truncate();
  extern void uppertruncate();
  
}

#endif
