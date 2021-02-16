#define CHEBNMAX 500

#include <complex>

using namespace std;
namespace cheb
{
  extern double xxx[CHEBNMAX][CHEBNMAX];
  
  //legendre polynomial j for an interpolation rule of order n;
  extern complex <double> Lnj(complex <double> a, complex <double> b, int n, int j, complex <double> z);
  
  //polynomial interpolation order n
  extern complex <double> ipol(complex <double> a, complex <double> b, int n, complex <double> f[], complex <double> z);

  extern void init();
}
