#ifndef mellinpdf_h
#define mellinpdf_h

#include <complex>
using namespace std;

//Transform PDF from x to N space
namespace mellinpdf
{
  void allocate();
  void free();

  extern double *t;              //Gauss nodes
  extern double *fac;            //overall factor (jacobian times Gauss weights)
  extern complex <double> *kern; //kernel of the Mellin transform
#pragma omp threadprivate(t,fac,kern)

  extern double *t_1;
  extern double *fac_1;
  extern complex <double> *kern_1;
#pragma omp threadprivate(t_1,fac_1,kern_1)
  extern double *t_2;
  extern double *fac_2;
  extern complex <double> *kern_2;
#pragma omp threadprivate(t_2,fac_2,kern_2)

  void init();
  void evalpdfs(double scale = -1, double m = -1, double y = -1);
  void update_mellin();
  void transform();
  void release();

  void gauss_quad();
  void laguerre_ipol();
  void ceres_pdf();

  extern complex <double> *UP;
  extern complex <double> *DO;
  extern complex <double> *ST;
  extern complex <double> *CH;
  extern complex <double> *BO;
  extern complex <double> *GL;
  extern complex <double> *UB;
  extern complex <double> *DB;
  extern complex <double> *SB;
  extern complex <double> *CB;
  extern complex <double> *BB;
#pragma omp threadprivate(UP,DO,ST,CH,BO,GL,UB,DB,SB,CB,BB)

  extern complex <double> *UP_1;
  extern complex <double> *DO_1;
  extern complex <double> *ST_1;
  extern complex <double> *CH_1;
  extern complex <double> *BO_1;
  extern complex <double> *GL_1;
  extern complex <double> *UB_1;
  extern complex <double> *DB_1;
  extern complex <double> *SB_1;
  extern complex <double> *CB_1;
  extern complex <double> *BB_1;
#pragma omp threadprivate(UP_1,DO_1,ST_1,CH_1,BO_1,GL_1,UB_1,DB_1,SB_1,CB_1,BB_1)

  extern complex <double> *UP_2;
  extern complex <double> *DO_2;
  extern complex <double> *ST_2;
  extern complex <double> *CH_2;
  extern complex <double> *BO_2;
  extern complex <double> *GL_2;
  extern complex <double> *UB_2;
  extern complex <double> *DB_2;
  extern complex <double> *SB_2;
  extern complex <double> *CB_2;
  extern complex <double> *BB_2;
#pragma omp threadprivate(UP_2,DO_2,ST_2,CH_2,BO_2,GL_2,UB_2,DB_2,SB_2,CB_2,BB_2)
  
  extern double *fup;
  extern double *fdo;
  extern double *fst;
  extern double *fch;
  extern double *fbo;
  extern double *fgl;
  extern double *fub;
  extern double *fdb;
  extern double *fsb;
  extern double *fcb;
  extern double *fbb;
#pragma omp threadprivate(fup,fdo,fst,fch,fbo,fgl,fub,fdb,fsb,fcb,fbb)

  extern double *fup_1;
  extern double *fdo_1;
  extern double *fst_1;
  extern double *fch_1;
  extern double *fbo_1;
  extern double *fgl_1;
  extern double *fub_1;
  extern double *fdb_1;
  extern double *fsb_1;
  extern double *fcb_1;
  extern double *fbb_1;
#pragma omp threadprivate(fup_1,fdo_1,fst_1,fch_1,fbo_1,fgl_1,fub_1,fdb_1,fsb_1,fcb_1,fbb_1)

  extern double *fup_2;
  extern double *fdo_2;
  extern double *fst_2;
  extern double *fch_2;
  extern double *fbo_2;
  extern double *fgl_2;
  extern double *fub_2;
  extern double *fdb_2;
  extern double *fsb_2;
  extern double *fcb_2;
  extern double *fbb_2;
#pragma omp threadprivate(fup_2,fdo_2,fst_2,fch_2,fbo_2,fgl_2,fub_2,fdb_2,fsb_2,fcb_2,fbb_2)
}
#endif
