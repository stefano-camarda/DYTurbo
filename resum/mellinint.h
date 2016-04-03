#ifndef mellinint_h
#define mellinint_h
#include "fcomplex.h"

#include <map>
#include <complex>

//fortran interfaces
extern "C"
{
  void mellinint_pdf_mesq_expy_(int& i1, int& i2, int& sign);
  fcomplex mellinint_integrand_(int& i1, int& i2, int& sign);
}

using namespace std;
namespace mellinint
{
  //DYRES convention
  // bb cb sb db ub g u d s c b
  // -5 -4 -3 -2 -1 0 1 2 3 4 5
  enum quark {bb=0, cb=1, sb=2, db=3, ub=4, g=5, u=6, d=7, s=8, c=9, b=10};
    
  extern int mdim;
  extern double *wn; //weights of the gaussian quadrature rule
  extern complex <double> *Np; //nodes on the positive branch of the contour
  extern complex <double> *Nm; //nodes on the negative branch of the contour
  extern complex <double> CCp,CCm; //angles of the positive and negative branches in the unitary complex circle
  extern void initgauss();

  extern void pdf_mesq_expy(int i1, int i2, int sign);
  extern complex <double> integrand(int i1, int i2, int sign);

  inline int index(int i1, int i2)
  {return i2 + mellinint::mdim*i1;}
  
  //luminosity times mesqij
  extern complex<double> GGN;
  extern complex<double> QGN_1;
  extern complex<double> QGN_2;
  extern complex<double> QQBN;
  extern complex<double> QQN;
  extern complex <double> QQN_1;
  extern complex <double> QQN_2;
  extern complex<double> QQPN_1;
  extern complex<double> QQPN_2;
}

#endif
