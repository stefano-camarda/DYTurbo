#ifndef scales_h
#define scales_h

#include <string>

using namespace std;

extern "C" {
  //  double dyalphas_mcfm_(double &q, double &amz, int &nloop);
  //  double dyalphas_lhapdf_(double &q); //no need for this, as there is a C++ LHAPDF function
  void scaleset_mcfm_(double &m, double &pt, double &mjj);
  int dynamic_fac_();
}

namespace scales
{
  extern double ren, ren2;
  extern double fac, fac2;
  extern double res, res2;
  //  extern double alphasmz;
  //  extern double alphas;
#pragma omp threadprivate(ren,fac,res,ren2,fac2,res2)
  
  //  void init();
  string func(int ff);
  void set(double m, double pt = 0, double mjj = 0);
  void form(double &scale2, int ff, double m, double pt = 0, double mjj = 0);
  void mcfm();
  void dyres(double m);
  void vjet();
}

#endif
