#include "mellinpdf.h"
#include "gaussrules.h"
#include "settings.h"
#include "mellinint.h"
#include "pdf.h"
#include "parton.h"
#include <complex>
#include <cmath>
#include <iostream>
#include <string.h>

using namespace parton;

//Mellin moments of PDFs
complex <double> *mellinpdf::UV;
complex <double> *mellinpdf::DV;
complex <double> *mellinpdf::US;
complex <double> *mellinpdf::DS;
complex <double> *mellinpdf::SP;
complex <double> *mellinpdf::SM;
complex <double> *mellinpdf::GL;
complex <double> *mellinpdf::CP;
complex <double> *mellinpdf::CM;
complex <double> *mellinpdf::BP;
complex <double> *mellinpdf::BM;

//Values of the PDFs at the nodes of the Gauss rule
double *mellinpdf::fuv;
double *mellinpdf::fdv;
double *mellinpdf::fus;
double *mellinpdf::fds;
double *mellinpdf::fsp;
double *mellinpdf::fsm;
double *mellinpdf::fgl;
double *mellinpdf::fcp;
double *mellinpdf::fcm;
double *mellinpdf::fbp;
double *mellinpdf::fbm;

void mellinpdf::allocate()
{
  UV = new complex <double> [mellinint::mdim];
  DV = new complex <double> [mellinint::mdim];
  US = new complex <double> [mellinint::mdim];
  DS = new complex <double> [mellinint::mdim];
  SP = new complex <double> [mellinint::mdim];
  SM = new complex <double> [mellinint::mdim];
  GL = new complex <double> [mellinint::mdim];
  CP = new complex <double> [mellinint::mdim];
  CM = new complex <double> [mellinint::mdim];
  BP = new complex <double> [mellinint::mdim];
  BM = new complex <double> [mellinint::mdim];

  fuv = new double [opts.pdfrule];
  fdv = new double [opts.pdfrule];
  fus = new double [opts.pdfrule];
  fds = new double [opts.pdfrule];
  fsp = new double [opts.pdfrule];
  fsm = new double [opts.pdfrule];
  fgl = new double [opts.pdfrule];
  fcp = new double [opts.pdfrule];
  fcm = new double [opts.pdfrule];
  fbp = new double [opts.pdfrule];
  fbm = new double [opts.pdfrule];
}

void mellinpdf::free()
{
  delete [] UV;
  delete [] DV;
  delete [] US;
  delete [] DS;
  delete [] SP;
  delete [] SM;
  delete [] GL;
  delete [] CP;
  delete [] CM;
  delete [] BP;
  delete [] BM;

  delete [] fuv;
  delete [] fdv;
  delete [] fus;
  delete [] fds;
  delete [] fsp;
  delete [] fsm;
  delete [] fgl;
  delete [] fcp;
  delete [] fcm;
  delete [] fbp;
  delete [] fbm;
}

//Cache PDFs, t values, jacobian, kernels
double *t;              //gauss nodes
double *fac;            //overall factor
complex <double> *kern; //kernel of the Mellin transform
void mellinpdf::init(double xmin)
{
  t = new double [opts.pdfrule];
  fac = new double [opts.pdfrule];
  kern = new complex <double> [opts.pdfrule*mellinint::mdim];

  // boundaries of integration
  double xmax = 1;
  double ll = log(xmax/xmin);
  double cx = 0.5;
  double mx = 0.5;
  for (int i = 0; i < opts.pdfrule; i++)
    {
      double x = cx+mx*gr::xxx[opts.pdfrule-1][i];
      t[i] = xmin*pow(xmax/xmin,x);
      //t[i] = xmin*exp(ll*x);
      double jac = mx * t[i] * ll;
      fac[i] = jac * gr::www[opts.pdfrule-1][i];
      for (int m = 0; m < mellinint::mdim; m++)
	kern[i*mellinint::mdim+m] = pow(t[i], mellinint::Np[m]-1.);
    }
}

//Evaluate PDFs in x-space at the required factorisation scale
void mellinpdf::evalpdfs(double scale)
{
  //factorization scale
  double muf = scale;
  //muf=opts.facscale
  //muf=2D0

  //Call the PDFs and store values at the knots
  double fx[11];
  for (int i = 0; i < opts.pdfrule; i++)
    {
      fdist_(opts.ih1,t[i],muf,fx);
      fuv[i] = fac[i] * (fx[U]-fx[Ub]);
      fdv[i] = fac[i] * (fx[D]-fx[Db]);
      fus[i] = fac[i] * fx[Ub];
      fds[i] = fac[i] * fx[Db];
      fsp[i] = fac[i] * fx[S];
      fsm[i] = fac[i] * fx[Sb];
      fgl[i] = fac[i] * fx[G];
      fcp[i] = fac[i] * fx[C];
      fcm[i] = fac[i] * fx[Cb];
      fbp[i] = fac[i] * fx[B];
      fbm[i] = fac[i] * fx[Bb];
      //cout << i << "  " << t[i] << "  " << fac[i] << "  " << fx[G] << endl;
    }
}

//Computes the complex N Mellin moments of pdfs by means of Gauss-Legendre quadrature
void mellinpdf::gauss_quad()
{
  //initialise the moments
  for (int m = 0; m < mellinint::mdim; m++)
    {
      UV[m] = 0.;
      DV[m] = 0.;
      US[m] = 0.;
      DS[m] = 0.;
      SP[m] = 0.;
      SM[m] = 0.;
      GL[m] = 0.;
      CP[m] = 0.;
      CM[m] = 0.;
      BP[m] = 0.;
      BM[m] = 0.;
    }
  
  //Calculate Mellin moments as:
  //integral_0^1{ x^(N-1) fx dx}
  for (int i = 0; i < opts.pdfrule; i++)
    for (int m = 0; m < mellinint::mdim; m++)
      {
	UV[m] += fuv[i] * kern[i*mellinint::mdim+m];
	DV[m] += fdv[i] * kern[i*mellinint::mdim+m];
	US[m] += fus[i] * kern[i*mellinint::mdim+m];
	DS[m] += fds[i] * kern[i*mellinint::mdim+m];
	SP[m] += fsp[i] * kern[i*mellinint::mdim+m];
	SM[m] += fsm[i] * kern[i*mellinint::mdim+m];
	GL[m] += fgl[i] * kern[i*mellinint::mdim+m];
	CP[m] += fcp[i] * kern[i*mellinint::mdim+m];
	CM[m] += fcm[i] * kern[i*mellinint::mdim+m];
	BP[m] += fbp[i] * kern[i*mellinint::mdim+m];
	BM[m] += fbm[i] * kern[i*mellinint::mdim+m];
	//cout << i << "  " << kern[i*mellinint::mdim+m] << endl;
	//if (m == 0) cout << m << "  " << i << "  " << GL[m] << endl;
      }
}

void laguerre(int n, double x, double *L)
{
  L[0] = 1.;
  double lkm2 = 0.; //--> L[-1]
  double lkm1 = L[0];// = 1.;
  for (int k = 1; k <= n; k++)
    {
      L[k] = ((2*k-1.-x)*lkm1-(k-1.)*lkm2)/k;
      lkm2=lkm1;
      lkm1=L[k];
    }
  return;
}

//Computes the complex N Mellin moments of pdfs by means of Laguerre interpolation and Laplace transform
void mellinpdf::laguerre_ipol()
{
  //order of the Laguerre approximation --> make this an option
  int N = 150;

  //scaling x -> x^p
  double p = 2.;

  //Calculate Laguerre projection
  double luv[N+1] = {0.};
  double ldv[N+1] = {0.};
  double lus[N+1] = {0.};
  double lds[N+1] = {0.};
  double lsp[N+1] = {0.};
  double lsm[N+1] = {0.};
  double lgl[N+1] = {0.};
  double lcp[N+1] = {0.};
  double lcm[N+1] = {0.};
  double lbp[N+1] = {0.};
  double lbm[N+1] = {0.};
  
  double cx = 0.5;
  double mx = 0.5;

  //calculate the x^alpha*(1-x)^beta approximation
  //calculate the 0, 1, 2 Mellin moments
  double fnuv[3] = {0.};
  double fndv[3] = {0.};
  double fnus[3] = {0.};
  double fnds[3] = {0.};
  double fnsp[3] = {0.};
  double fnsm[3] = {0.};
  double fngl[3] = {0.};
  double fncp[3] = {0.};
  double fncm[3] = {0.};
  double fnbp[3] = {0.};
  double fnbm[3] = {0.};


  for (int i = 0; i < opts.pdfrule; i++)
    {
      //double x = cx+mx*xxx[rule-1][i];
      //double t = a*pow(1./a,x);
      //double jac = 0.5* t * log(1./a);
      //f0 += fknots[i]      * fac[i];
      //f1 += fknots[i] *t   * fac[i];
      //f2 += fknots[i] *t*t * fac[i];
      for (int n = 0; n < 3; n++)
	{
	  fnuv[n] += fuv[i] * pow(t[i],n) *fac[i];
	  fndv[n] += fdv[i] * pow(t[i],n) *fac[i];
	  fnus[n] += fus[i] * pow(t[i],n) *fac[i];
	  fnds[n] += fds[i] * pow(t[i],n) *fac[i];
	  fnsp[n] += fsp[i] * pow(t[i],n) *fac[i];
	  fnsm[n] += fsm[i] * pow(t[i],n) *fac[i];
	  fngl[n] += fgl[i] * pow(t[i],n) *fac[i];
	  fncp[n] += fcp[i] * pow(t[i],n) *fac[i];
	  fncm[n] += fcm[i] * pow(t[i],n) *fac[i];
	  fnbp[n] += fbp[i] * pow(t[i],n) *fac[i];
	  fnbm[n] += fbm[i] * pow(t[i],n) *fac[i];
	}
    }

  //This is used to regularise the PDFs with the same x^alpha behaviour at low x
  double alpha = 2;
  double powuv = alpha - fnuv[1]*(fnuv[2]-fnuv[1])/(fnuv[1]*fnuv[1]-fnuv[0]*fnuv[2]);
  double powdv = alpha - fndv[1]*(fndv[2]-fndv[1])/(fndv[1]*fndv[1]-fndv[0]*fndv[2]);
  double powus = alpha - fnus[1]*(fnus[2]-fnus[1])/(fnus[1]*fnus[1]-fnus[0]*fnus[2]);
  double powds = alpha - fnds[1]*(fnds[2]-fnds[1])/(fnds[1]*fnds[1]-fnds[0]*fnds[2]);
  double powsp = alpha - fnsp[1]*(fnsp[2]-fnsp[1])/(fnsp[1]*fnsp[1]-fnsp[0]*fnsp[2]);
  double powsm = alpha - fnsm[1]*(fnsm[2]-fnsm[1])/(fnsm[1]*fnsm[1]-fnsm[0]*fnsm[2]);
  double powgl = alpha - fngl[1]*(fngl[2]-fngl[1])/(fngl[1]*fngl[1]-fngl[0]*fngl[2]);
  double powcp = alpha - fncp[1]*(fncp[2]-fncp[1])/(fncp[1]*fncp[1]-fncp[0]*fncp[2]);
  double powcm = alpha - fncm[1]*(fncm[2]-fncm[1])/(fncm[1]*fncm[1]-fncm[0]*fncm[2]);
  double powbp = alpha - fnbp[1]*(fnbp[2]-fnbp[1])/(fnbp[1]*fnbp[1]-fnbp[0]*fnbp[2]);
  double powbm = alpha - fnbm[1]*(fnbm[2]-fnbm[1])/(fnbm[1]*fnbm[1]-fnbm[0]*fnbm[2]);

  //The largest poles of the Mellin transform are at these values of z:
  //cout << powuv << "  "
  //     << powdv << "  "
  //     << powus << "  "
  //     << powds << "  "
  //     << powsp << "  "
  //     << powsm << "  "
  //     << powgl << "  "
  //     << powch << "  "
  //     << powbo << endl;

  
  //double alpha = 2;
  for (int i = 0; i < opts.pdfrule; i++)
    {
      double L[N+1];
      double u = -log(t[i])/p;
      laguerre(N, u, L);// --> cache it as L[n][i]
      //      cout << endl << i << "  " << u << endl;
      //double po = pow(t,alpha);
      for (int n = 0; n <= N; n++)
	{
	  luv[n] += fuv[i] * L[n] * pow(t[i],1./p-1.+powuv);
	  ldv[n] += fdv[i] * L[n] * pow(t[i],1./p-1.+powdv);
	  lus[n] += fus[i] * L[n] * pow(t[i],1./p-1.+powus);
	  lds[n] += fds[i] * L[n] * pow(t[i],1./p-1.+powds);
	  lsp[n] += fsp[i] * L[n] * pow(t[i],1./p-1.+powsp);
	  lsm[n] += fsm[i] * L[n] * pow(t[i],1./p-1.+powsm);
	  lgl[n] += fgl[i] * L[n] * pow(t[i],1./p-1.+powgl);
	  lcp[n] += fcp[i] * L[n] * pow(t[i],1./p-1.+powcp);
	  lcm[n] += fcm[i] * L[n] * pow(t[i],1./p-1.+powcm);
	  lbp[n] += fbp[i] * L[n] * pow(t[i],1./p-1.+powbp);
	  lbm[n] += fbm[i] * L[n] * pow(t[i],1./p-1.+powbm);

//	  cout << i << "  " << n << "  " << L[n] << endl;
	}
    }
  //for (int n = 0; n <= N; n++)
  //printf("Laguerre projection %d \t (%e) \n", n, lgl[n]);
  
  //initialise the moments
  for (int m = 0; m < mellinint::mdim; m++)
    {
      UV[m] = 0.;
      DV[m] = 0.;
      US[m] = 0.;
      DS[m] = 0.;
      SP[m] = 0.;
      SM[m] = 0.;
      GL[m] = 0.;
      CP[m] = 0.;
      CM[m] = 0.;
      BP[m] = 0.;
      BM[m] = 0.;
    }
  //Calculate Mellin moments from the
  //Mellin transform of Laguerre polynomials
  for (int m = 0; m < mellinint::mdim; m++)
    {
      //complex <double> z = (mellinint::Np[m] - alpha)*p;
      complex <double> zuv = (mellinint::Np[m] - powuv)*p;
      complex <double> zdv = (mellinint::Np[m] - powdv)*p;
      complex <double> zus = (mellinint::Np[m] - powus)*p;
      complex <double> zds = (mellinint::Np[m] - powds)*p;
      complex <double> zsp = (mellinint::Np[m] - powsp)*p;
      complex <double> zsm = (mellinint::Np[m] - powsm)*p;
      complex <double> zgl = (mellinint::Np[m] - powgl)*p;
      complex <double> zcp = (mellinint::Np[m] - powcp)*p;
      complex <double> zcm = (mellinint::Np[m] - powcm)*p;
      complex <double> zbp = (mellinint::Np[m] - powbp)*p;
      complex <double> zbm = (mellinint::Np[m] - powbm)*p;
      for (int n = 0; n <= N; n++)
	{
	  //complex <double> lapl = 1./z * pow(((z-1.)/z),n);
	  //complex <double> lapl = 1./z * pow(1.-1./z,n); //cache this as lapl[m][n]

	  complex <double> lapluv = 1./zuv * pow(1.-1./zuv,n); //cache this as lapluv[m][n]
	  complex <double> lapldv = 1./zdv * pow(1.-1./zdv,n); //cache this as lapldv[m][n]
	  complex <double> laplus = 1./zus * pow(1.-1./zus,n); //cache this as laplus[m][n]
	  complex <double> laplds = 1./zds * pow(1.-1./zds,n); //cache this as laplds[m][n]
	  complex <double> laplsp = 1./zsp * pow(1.-1./zsp,n); //cache this as laplsp[m][n]
	  complex <double> laplsm = 1./zsm * pow(1.-1./zsm,n); //cache this as laplsm[m][n]
	  complex <double> laplgl = 1./zgl * pow(1.-1./zgl,n); //cache this as laplgl[m][n]
	  complex <double> laplcp = 1./zcp * pow(1.-1./zcp,n); //cache this as laplcp[m][n]
	  complex <double> laplcm = 1./zcm * pow(1.-1./zcm,n); //cache this as laplcm[m][n]
	  complex <double> laplbp = 1./zbp * pow(1.-1./zbp,n); //cache this as laplbp[m][n]
	  complex <double> laplbm = 1./zbm * pow(1.-1./zbm,n); //cache this as laplbm[m][n]

	  UV[m] += luv[n] * lapluv;
	  DV[m] += ldv[n] * lapldv;
	  US[m] += lus[n] * laplus;
	  DS[m] += lds[n] * laplds;
	  SP[m] += lsp[n] * laplsp;
	  SM[m] += lsm[n] * laplsm;
	  GL[m] += lgl[n] * laplgl;
	  CP[m] += lcp[n] * laplcp;
	  CM[m] += lcm[n] * laplcm;
	  BP[m] += lbp[n] * laplbp;
	  BM[m] += lbm[n] * laplbm;
	  //cout << m << "  " << n << "  " << real(fli[n] * lapl) << "  " << z << "  " << lapl << endl;
	}
    }
  
}
