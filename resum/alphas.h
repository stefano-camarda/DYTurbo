#ifndef alphas_h
#define alphas_h

//#include "resconst.h"

#include <complex>
using namespace std;
//using namespace resconst;

namespace alphas
{
  const int nastps = 10;

  const int NFMIN = 3;
  const int NFMAX = 6;

  extern double bet0[NFMAX-NFMIN+1];
  extern double bet1[NFMAX-NFMIN+1];
  extern double bet2[NFMAX-NFMIN+1];
  extern double bet3[NFMAX-NFMIN+1];
  extern double bet4[NFMAX-NFMIN+1];

  //NF = 5 beta function
  //extern double beta0;
  //extern double beta1;
  //extern double beta2;
  //extern double beta3;
  //extern double beta4;
  
  inline complex <double> fbeta1(complex <double> a, int nf) {return -pow(a,2)*(bet0[nf-NFMIN]+a*bet1[nf-NFMIN]);};
  inline complex <double> fbeta2(complex <double> a, int nf) {return -pow(a,2)*(bet0[nf-NFMIN]+a*bet1[nf-NFMIN]+pow(a,2)*bet2[nf-NFMIN]);};
  inline complex <double> fbeta3(complex <double> a, int nf) {return -pow(a,2)*(bet0[nf-NFMIN]+a*bet1[nf-NFMIN]+pow(a,2)*bet2[nf-NFMIN]+pow(a,3)*bet3[nf-NFMIN]);};
  inline complex <double> fbeta4(complex <double> a, int nf) {return -pow(a,2)*(bet0[nf-NFMIN]+a*bet1[nf-NFMIN]+pow(a,2)*bet2[nf-NFMIN]+pow(a,3)*bet3[nf-NFMIN]+pow(a,4)*bet4[nf-NFMIN]);};
  extern void rgkt(complex <double> q, double q0, double aass0, int nf = 5);
  extern void iter(complex <double> q, double q0, double aass0, int nf = 5);

  //alphas at order nloop
  extern complex <double> as;
#pragma omp threadprivate(as)
  
  //alphas at all computed orders (up to 4-loop)
  extern double as0;      //0-loop
  extern complex <double> asLO;     //1-loop
  extern complex <double> asNLO;    //2-loop
  extern complex <double> asNNLO;   //3-loop
  extern complex <double> asNNNLO;  //4-loop
  extern complex <double> asNNNNLO; //5-loop
#pragma omp threadprivate(as0,asLO,asNLO,asNNLO,asNNNLO,asNNNNLO)

  //powers of alphas at truncated orders
  extern complex <double> as1_1l; //alphas truncated at 1-loop
  extern complex <double> as1_2l; //alphas truncated at 2-loop
  extern complex <double> as1_3l; //alphas truncated at 3-loop
  extern complex <double> as1_4l; //alphas truncated at 4-loop
  extern complex <double> as1_5l; //alphas truncated at 5-loop
  extern complex <double> as2_2l; //alphas truncated powers of 2 at 2-loop
  extern complex <double> as2_3l; //alphas truncated powers of 2 at 3-loop
  extern complex <double> as2_4l; //alphas truncated powers of 2 at 4-loop
  extern complex <double> as2_5l; //alphas truncated powers of 2 at 5-loop
  extern complex <double> as3_3l; //alphas truncated powers of 3 at 3-loop
  extern complex <double> as3_4l; //alphas truncated powers of 3 at 4-loop
  extern complex <double> as3_5l; //alphas truncated powers of 3 at 5-loop
  extern complex <double> as4_4l; //alphas truncated powers of 4 at 4-loop
  extern complex <double> as4_5l; //alphas truncated powers of 4 at 5-loop
  extern complex <double> as5_5l; //alphas truncated powers of 5 at 5-loop
#pragma omp threadprivate(as1_1l,as1_2l,as1_3l,as1_4l,as1_5l,as2_2l,as2_3l,as2_4l,as2_5l,as3_3l,as3_4l,as3_5l,as4_4l,as4_5l,as5_5l)

  extern void init();
  extern complex <double> calc(complex <double> q, int nloop);
}

#endif
