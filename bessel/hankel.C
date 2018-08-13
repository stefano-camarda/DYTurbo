#include "hankel.h"
#include "specialfunctions_interface.h"
#include "besselint.h"
#include <iostream>

//rewrite steven-murray python code from https://github.com/steven-murray/hankel/blob/master/hankel/hankel.py

//General quadrature method for Hankel transformations.
//
//Based on the algorithm provided in
//H. Ogata, A Numerical Integration Formula Based on the Bessel Functions,
//Publications of the Research Institute for Mathematical Sciences,
//vol. 41, no. 4, pp. 949-970, 2005.


//    The basis of the Hankel Transformation algorithm by Ogata 2005.
//
//    This algorithm is used to solve the equation :math:`\int_0^\infty f(x) J_\nu(x) dx`
//    where :math:`J_\nu(x)` is a Bessel function of the first kind of order
//    :math:`nu`, and :math:`f(x)` is an arbitrary (slowly-decaying) function.
//
//    The algorithm is presented in
//    H. Ogata, A Numerical Integration Formula Based on the Bessel Functions,
//    Publications of the Research Institute for Mathematical Sciences, vol. 41, no. 4, pp. 949-970, 2005.
//
//    This class provides a method for directly performing this integration, and also
//    for doing a Hankel Transform.
//
//    Parameters
//    ----------
//
//    nu : int or 0.5, optional, default = 0
//        The order of the bessel function (of the first kind) J_nu(x)
//
//    N : int, optional, default = 3.2/`h`
//        The number of nodes in the calculation. Generally this must increase
//        for a smaller value of the step-size h. Default value is based on where the series will truncate according
//        to the double-exponential convergence to the roots of the Bessel function.
//
//    h : float, optional, default = 0.1
//        The step-size of the integration.

int hankel::nu;
int hankel::n;
double hankel::h;

double *hankel::zeros;
double *hankel::w; //weights
double *hankel::x; //nodes
double *hankel::j; //Jn evaluated at the nodes
double *hankel::dp;

double *hankel::series;

//Some quantities only useful in the FourierTransform
int hankel::x_power = 1;
int hankel::k_power = 2;

void hankel::init(int Nu, int N, double H)
{
  nu = Nu;
  h = H;
  n = N;
  
  if (n <= 0)
    n = int(3.2/h);

  //fill the bessel zeros
  roots();

  //fill the nodes
  fillx();

  //Commpute the bessel functions at the nodes
  fillj();
  
  //compute the weights
  weight();
  
  filldpsi();
}

void hankel::roots()
{
  zeros = new double[n];

  double rj0[n];
  double rj1[n];
  double ry0[n];
  double ry1[n];
  jyzo_ (nu, n, rj0, rj1, ry0, ry1);

  for (int i = 0; i < n; i++)
    zeros[i] = rj0[i]/M_PI;
}

double hankel::psi(double t)
{
  return t*tanh(M_PI*sinh(t)/2.);
  double K = 6.;
  //return t/(1.-exp(-K*sinh(t)));
}

double hankel::dpsi(double t)
{
  return (M_PI*t*cosh(t) + sinh(M_PI*sinh(t)))/(1.0 + cosh(M_PI*sinh(t)));
  double K = 6.;
  //return (1.-(1.+K*t*cosh(t))*exp(-K*sinh(t)))/pow(1.-exp(-K*sinh(t)),2);
}

void hankel::fillx()
{
  x = new double[n];

  for (int i = 0; i < n; i++)
    x[i] = M_PI*psi(h*zeros[i])/h;
}

double hankel::jn(double x)
{
  if (nu == 0)
    return fort_besj0_(x);
  else if (nu == 1)
    return fort_besj1_(x);
  else
    return fort_besjn_(nu, x);
}

double hankel::fillj()
{
  j = new double[n];

  for (int i = 0; i < n; i++)
    j[i] = jn(x[i]);
}      


double hankel::filldpsi()
{
  dp = new double[n];
  
  for (int i = 0; i < n; i++)
    {
      //double a = (M_PI*h*zeros[i]*cosh(h*zeros[i]) + sinh(M_PI*sinh(h*zeros[i])))/(1.0 + cosh(M_PI*sinh(h*zeros[i])));
      double a = dpsi(h*zeros[i]);
      if (isnan(a))
	a = 1.0;
      dp[i] = a;
    }
}

//this is J_n+1
double hankel::jn1(double x)
{
  int nup1 = nu+1;
  if (nup1 == 0)
    return fort_besj0_(x);
  if (nup1 == 1)
    return fort_besj1_(x);
  else
    return fort_besjn_(nup1, x);
}

//this is Y_n
double hankel::yn(double x)
{
  if (nu == 0)
    return fort_besy0_(x);
  else if (nu == 1)
    return fort_besy1_(x);
  else
    return fort_besyn_(nu, x);
}

double hankel::weight()
{
  w = new double[n];

  for (int i = 0; i < n; i++)
    w[i] = yn(M_PI*zeros[i])/jn1(M_PI*zeros[i]);
}


double hankel::get_series(double (*f)(double), double k)
{
  series = new double[n];
  //formula (5.2) of http://www.kurims.kyoto-u.ac.jp/~okamoto/paper/Publ_RIMS_DE/41-4-40.pdf
  for (int i = 0; i < n; i++)
    {
      //evaluate f(x)*x at the nodes
      double fres = f(x[i]/k)*pow(x[i],x_power); //x_power is 1
    
      series[i] = M_PI*w[i]*fres*j[i]*dp[i];
    }
}

void hankel::transform(double (*f)(double), double k, double &res, double &err)
//        Do the Hankel-transform of the function f.
//
//        Parameters
//        ----------
//        f : callable
//            A function of one variable, representing :math:`f(x)`
//
//        ret_err : boolean, optional, default = True
//            Whether to return the estimated error
//
//        ret_cumsum : boolean, optional, default = False
//            Whether to return the cumulative sum
//
//        Returns
//        -------
//        ret : array-like
//            The Hankel-transform of f(x) at the provided k. If
//            `k` is scalar, then this will be scalar.
//
//        err : array-like
//            The estimated error of the approximate integral, at every `k`.
//            It is merely the last term in the sum. Only returned if `ret_err=True`.
//
//        cumsum : array-like
//            The total cumulative sum, for which the last term is itself the transform.
//            One can use this to check whether the integral is converging.
//            Only returned if `ret_cumsum=True`
//
//
//        Notes
//        -----
//        The Hankel transform is defined as
//
//        .. math:: F(k) = \int_0^\infty r f(r) J_\nu(kr) dr.
//
//        The inverse transform is identical (swapping *k* and *r* of course).
{

  // The following is the scalar normalisation of the transform
  // The basic transform has a norm of 1
  double norm = 1.;

  // The following renormalises by the fourier dual to some power
  double knorm = pow(k,k_power); //k_power is 2

  get_series(f,k);

  res = 0.;
  for (int i = 0; i < n; i++)
    res += series[i];

  res *= norm/knorm;

  //use the last term of series to estimate the error
  err = norm/knorm*series[n-1];

  //  for (int i = 0; i < n; i++)
  //    cout << norm/knorm*series[i] << endl;
}

void hankel::free()
{
  delete[] zeros;
  delete[] x;
  delete[] j;
  delete[] dp;
  delete[] w;
  delete[] series;
}
