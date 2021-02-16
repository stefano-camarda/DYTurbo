#include "chebyshev.h"
#include "gaussrules.h"

double cheb::xxx[CHEBNMAX][CHEBNMAX];

//legendre polynomial j for an interpolation of order n;
complex <double> cheb::Lnj(complex <double> a, complex <double> b, int n, int j, complex <double> z)
{
  complex <double> c = 0.5*(a+b);
  complex <double> m = 0.5*(b-a);

  complex <double> res = 1;
  for (int k = 1; k <= n; k++)
    if (k != j)
      //res *= (z-(c+m*gr::xxx[n-1][k-1]))/((c+m*gr::xxx[n-1][j-1])-(c+m*gr::xxx[n-1][k-1]));
      res *= (z-(c+m*xxx[n-1][k-1]))/((c+m*xxx[n-1][j-1])-(c+m*xxx[n-1][k-1]));

  return res;
}

//polynomial interpolation of order n
complex <double> cheb::ipol(complex <double> a, complex <double> b, int n, complex <double> f[], complex <double> z)
{
  complex <double> res = 0.;
  for (int j = 1; j <= n; j++)
    res += f[j-1]*Lnj(a,b,n,j,z);

  return res;
}

void cheb::init()
{
  for (int n = 1; n <= CHEBNMAX; n++)
    for ( int i = 1; i <= n; i++ )
      xxx[n-1][i-1] = cos((2.*i-1.)/(2.*n) * M_PI);
}
