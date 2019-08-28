#include "mellinint.h"
#include "mesq.h"
#include "settings.h"
#include "gaussrules.h"
#include "clenshawcurtisrules.h"
#include "hcoefficients.h"
#include "hcoeff.h"
#include "pdfevol.h"
#include "interface.h"
#include "anomalous.h"
#include "parton.h"
#include "string.h"

#include <complex>
#include <iostream>
#include <iomanip>

using namespace parton;

//double *mellinint::wn;
complex <double> *mellinint::wn;
complex <double> *mellinint::Np;
complex <double> *mellinint::Nm;
complex <double> mellinint::CCp;
complex <double> mellinint::CCm;

complex <double> mellinint::GGN;
complex <double> mellinint::QGN_1;
complex <double> mellinint::QGN_2;
complex <double> mellinint::QQBN;
complex <double> mellinint::QQN;
complex <double> mellinint::QQN_1;
complex <double> mellinint::QQN_2;
complex <double> mellinint::QQPN_1;
complex <double> mellinint::QQPN_2;

int mellinint::mdim;

//fortran interfaces
void mellinint_pdf_mesq_expy_(int& i1, int& i2, int& sign)
{
  mellinint::pdf_mesq_expy(i1-1, i2-1, sign-1);
};
fcomplex mellinint_integrand_(int& i1, int& i2, int& sign)
{
  return fcx(mellinint::integrand2d(i1-1, i2-1, sign-1));
};

void mellinint::initgauss()
{
  //set up weights and nodes for gaussian quadrature
  //cpoint is the starting point on the real axis for the positive and negative part of the integration path
  //phi is the angle in the complex plane of the positive part of the integration path
       
  double cpoint = opts.cpoint;
  //double phi = M_PI * 1./2.;
  //double phi = M_PI * 3./4.;
  double phi = M_PI * opts.phi;

  double zmin = 0;
  //upper limit for the mellin integration in the complex plane (z)
  //Above 50 the integral becomes unstable, especially when m_ll << mur,
  //due to large logs in the Sudakov, which can be reduced with smaller blim.
  //The issue can be avoided by using dynamicscale
  //Also, larger values of zmax requires more support points for the
  //Mellin inverse transform, i.e. higher mellinintervals or mellinrule.
  double zmax = opts.zmax;

  mdim = opts.mellinintervals*opts.mellinrule;
  
  //allocate memory
  //wn = new double [mdim];
  wn = new complex <double> [mdim];
  Np = new complex <double> [mdim];
  Nm = new complex <double> [mdim];
  
  //initialise to 0
  for (int i=0; i < mdim; i++)
    {
      Np[i]=cpoint+1.;
      Nm[i]=cpoint+1.;
      wn[i]=0;
    }

  //imaginary unit
  complex <double> ii = complex <double>(0.,1.);
  
  //Gauss-Legendre quadrature along a linear contour
  if (opts.mellininv == 0)
    {
      // positive branch      
      CCp = complex <double> (cos(phi), sin(phi));
      for (int i=0; i < opts.mellinintervals; i++)
	{
	  double a = 0.+(1.-0.)*i/opts.mellinintervals;
	  double b = 0.+(1.-0.)*(i+1)/opts.mellinintervals;
	  double c = 0.5*(a+b);
	  double m = 0.5*(b-a);
	  for (int j=0; j < opts.mellinrule; j++)
	    {
	      double x = c+m*gr::xxx[opts.mellinrule-1][j];
	      //double x = c+m*cc::xxx[opts.mellinrule-1][j];
	      double t = zmin+(zmax-zmin)*x;
	      double jac = zmax-zmin;
	      Np[j+i*opts.mellinrule]=complex <double> (cpoint+cos(phi)*t+1.,sin(phi)*t);
	      wn[j+i*opts.mellinrule]=gr::www[opts.mellinrule-1][j]*m*jac;
	      //wn[j+i*opts.mellinrule]=cc::www[opts.mellinrule-1][j]*m*jac;
	      //cout << setprecision(16) <<  t << " " << Np[j+i*opts.mellinrule] << "  " << wn[j+i*opts.mellinrule] << endl;
	    }
	}      
      // negative branch
      CCm = complex <double> (cos(phi), -sin(phi));
      for (int i=0; i < opts.mellinintervals; i++)
	{
	  double a = 0.+(1.-0.)*i/opts.mellinintervals;
	  double b = 0.+(1.-0.)*(i+1)/opts.mellinintervals;
	  double c = 0.5*(a+b);
	  double m = 0.5*(b-a);
	  for (int j=0; j < opts.mellinrule; j ++)
	    {
	      double x = c+m*gr::xxx[opts.mellinrule-1][j];
	      //double x = c+m*cc::xxx[opts.mellinrule-1][j];
	      double t = zmin+(zmax-zmin)*x;
	      Nm[j+i*opts.mellinrule]=complex <double> (cpoint+cos(phi)*t+1.,-sin(phi)*t);
	      //cout << setprecision(16) <<  t << " " << Nm[j+i*opts.mellinrule] << endl;
	    }
	}
    }

  //Talbot contour
  else if (opts.mellininv == 1)
    {
      //calculate the range of the corresponding Laplace transform
      double xmin, xmax;
      double t0, t1;

      //could restrict xmin and xmax further for limited y range
      xmin = pow(bins.mbins.front()/opts.sroot,2);
      xmax = 1;

      t0 = -log(xmax);
      t1 = -log(xmin);

      double tmid = (t0+t1)/2.;

      
      double Not, sigma, lambda, nu, alpha, mu;

      //fixed talbot of Peter Valko'
      //lambda = 2*opts.mellinrule/(5*tmid);
      //nu = 1.;

      //Modified Talbot of Rizzardi
      //lambda = 4.8/tmid;
      //nu = 1.;

      //Empirically good for Talbot
      sigma = 0.6;
      lambda = 0.6;
      nu = 2.5;

      //Weideman
      //Not    = opts.mellinrule/tmid;
      //sigma  = -0.6122;
      //mu     = 0.5017;
      //alpha  = 0.6407;
      //nu     = 0.2645;

      //Empirically good for Weideman
      //Not = 1.;
      //sigma = 0.6;
      //mu = 0.8;
      //alpha = 1.;
      //nu = 2.;

      bool midpoint = false;
      bool weideman = false;
      
      CCp = 1;
      for (int j=0; j < opts.mellinrule; j++)
	{
	  double x;
	  if (midpoint)
	    x = (double(j)+0.5)/opts.mellinrule;   //midpoint
	  else
	    x = double(j)/opts.mellinrule;         //trapezoidal

	  double thmax = M_PI;
	  double theta = thmax*x;

	  complex <double> s, jac;
	  if (weideman)
	    {
	      //Weideman contour
	      s = Not*(sigma+(theta==0?mu/alpha:mu*theta/tan(alpha*theta)+theta*nu*ii));
	      jac = thmax * Not * (ii * nu + (theta==0?0 : mu * (1./tan(alpha*theta) - alpha*theta*(1. + 1./tan(alpha*theta)/tan(alpha*theta)))));
	    }
	  else
	    {
	      //Talbot contour
	      s = sigma+(theta==0?lambda:lambda*theta*(1./tan(theta)+nu*ii));
	      jac = thmax * (theta==0? ii*lambda*nu : ii * lambda *( nu + ii * (theta + (theta /tan(theta) - 1.) / tan(theta))));
	    }


	  
	  Np[j] = s+1.;
	  if (midpoint)
	    wn[j] = jac/double(opts.mellinrule);                                //midpoint
	  else
	    wn[j] = (j==0 ? jac/2.:jac)/double(opts.mellinrule);                //trapezoidal
	  
	  //cout << setprecision(16) <<  theta << " " << Np[j] << "  " << wn[j] << endl;

	  /*
	  double x = 0.5+0.5*gr::xxx[opts.mellinrule-1][j];
	  double t = zmin+(zmax-zmin)*x;
	  complex <double> jac = (zmax-zmin) * (cos(phi)+ii*sin(phi));
	  complex <double> s = cpoint + t * (cos(phi)+ii*sin(phi));
	  Np[j] = s + 1.;
	  wn[j] = gr::www[opts.mellinrule-1][j]*0.5*jac;
	  cout << setprecision(16) <<  t << " " << Np[j] << "  " << wn[j] << endl;
	  */
	}
      // negative branch
      CCm = 1;
      for (int j=0; j < opts.mellinrule; j ++)
	{
	  double x;
	  if (midpoint)
	    x = (double(j)+0.5)/opts.mellinrule;     //midpoint
	  else
	    x = double(j)/opts.mellinrule;           //trapezoidal
	  
	  double thmax = -M_PI;
	  double theta = thmax*x;
	  
	  complex <double> s;
	  if (weideman) //Weideman contour
	    s = Not*(sigma+(theta==0?mu/alpha:mu*theta/tan(alpha*theta)+theta*nu*ii));
	  else//Talbot contour
	    s = sigma+(theta==0?lambda:lambda*theta*(1./tan(theta)+nu*ii));

	  Nm[j] = s+1.;
	  //cout << setprecision(16) <<  theta << " " << Nm[j] << "  " << wn[j] << endl;

	  /*
	  double x = 0.5+0.5*gr::xxx[opts.mellinrule-1][j];
	  double t = zmin+(zmax-zmin)*x;
	  complex <double> s = cpoint + t * (cos(phi)-ii*sin(phi));
	  Nm[j] = s + 1.;
	  cout << setprecision(16) <<  t << " " << Nm[j]  << endl;
	  */
	}
    }
  else
    {
      //Implement here other mellin inversions
      cout << "Not valid option for mellininv (should be 0 or 1) " << endl;
      exit (-1);
    }
}


//This function performs the product of PDF, born level amplitudes and expy piece, and sums over partonic channels (i,j)
//It is a function of z1 and z2 in Mellin space

// The integrand of the inverse Mellin transform is composed of
// A) parton luminosities in Mellin space: fn1 * fn2 -> GGN, QGN, QQBN, QQN, QQPN which are functions of the quark flavours (i,j) and Mellin indices (i1, i2)
// B) the EW born level squared amplitudes mesq::mesqij, which are functions of m, costh and quark flavours (i,j) (if integrated in costh, and in the presence of cuts on the leptons, the costh moments are functions of y, pt, m)
// C) the x1^-z1 * x2^-z2 piece of the Mellin inverse transform, which depends explicitly on rapidity and on the Mellin indices (i1, i2)
// D) the Wilson coefficients Hgg, Hqqb, etc.. which are functions of the Mellin indices (i1, i2)

//In mode = 1, the integration in costh and phi enters only in B
//In mode = 2, the integration in costh, phi and rapidity enters in B and C, hence B and C are calculated together in mesq::mesqij_expy
//The product with the Wilson coefficients (D) can be performed at the very end, in mellinint::integrand

//input : pdfevol::fn1, pdfevol::fn2, mesq::mesqij_expy
//output: qqbn, qgn, ggn

//#define fn1(x) pdfevol::fx1[i1*11+x]
//#define fn2(x) pdfevol::fx2[i2*11+x]

void mellinint::pdf_mesq_expy(int i1, int i2, int sign)
{
  //QQBN=0;
  if (opts.order < 1)
    {
      QGN_1=0;
      QGN_2=0;
    }
      
  if (opts.order < 2)
    {
      GGN=0;
      //QQN=0;
      QQN_1=0;
      QQN_2=0;
      QQPN_1=0;
      QQPN_2=0;
    }

  complex<double>* fn1 = pdfevol::fn1;
  complex<double>* fn2 = pdfevol::fn2;

  //complex<double> fn1[11];
  //complex<double> fn2[11];
  //memcpy(fn1, &(pdfevol::fx1[i1*11]), 11*sizeof(complex<double>));
  //if (sign == mesq::positive)
  //memcpy(fn2, &(pdfevol::fx2[i2*11]), 11*sizeof(complex<double>));
  //else
  //memcpy(fn2, &(pdfevol::fx2[mellinint::mdim*11+i2*11]), 11*sizeof(complex<double>));


  //set b to 0
  //  fn2[bb] = 0;  fn1[bb] = 0;
  //  fn2[b ] = 0;  fn1[b ]  = 0;
  
  //set s and c to 0
  //  fn2[cb] = 0;  fn1[cb] = 0;
  //  fn2[sb] = 0;  fn1[sb] = 0;
  //  fn2[s ] = 0;  fn1[s ]  = 0;
  //  fn2[c ] = 0;  fn1[c ]  = 0;

  //set u and d to 0
  //  fn2[db] = 0; fn1[db] = 0;
  //  fn2[ub] = 0; fn1[ub] = 0;
  //  fn2[u ] = 0; fn1[u ]  = 0;
  //  fn2[d ] = 0; fn1[d ]  = 0;

  //set gluon to 0
  //  fn2[g] = 0;  fn1[g] = 0;
  
  //DYRES convention
  // bb cb sb db ub g u d s c b
  // -5 -4 -3 -2 -1 0 1 2 3 4 5
    
  //mesqij[n] -> mesqij_expy[mesq::index(n,i1,i2,sign)]
  if (opts.nproc == 3)
    {
      complex <double> mesq_uub = mesq::mesqij_expy[mesq::index(0,i1,i2,sign)];
      complex <double> mesq_ubu = mesq::mesqij_expy[mesq::index(1,i1,i2,sign)];
      complex <double> mesq_ddb = mesq::mesqij_expy[mesq::index(2,i1,i2,sign)];
      complex <double> mesq_dbd = mesq::mesqij_expy[mesq::index(3,i1,i2,sign)];
      complex <double> mesq_ssb = mesq::mesqij_expy[mesq::index(4,i1,i2,sign)];
      complex <double> mesq_sbs = mesq::mesqij_expy[mesq::index(5,i1,i2,sign)];
      complex <double> mesq_ccb = mesq::mesqij_expy[mesq::index(6,i1,i2,sign)];
      complex <double> mesq_cbc = mesq::mesqij_expy[mesq::index(7,i1,i2,sign)];
      complex <double> mesq_bbr = mesq::mesqij_expy[mesq::index(8,i1,i2,sign)];
      complex <double> mesq_brb = mesq::mesqij_expy[mesq::index(9,i1,i2,sign)];
      
      //      cout << "pdf_mesq_expy " << i1 << "  " << i2 << "  " << mesq_uub << "  " << fn1[u] << "  " << fn2[ub]  << endl;
      //LL part
      QQBN = fn1[u ]*fn2[ub]*mesq_uub
	    +fn1[c ]*fn2[cb]*mesq_ccb
	    +fn1[d ]*fn2[db]*mesq_ddb
	    +fn1[s ]*fn2[sb]*mesq_ssb
	    +fn1[b ]*fn2[bb]*mesq_bbr
	    +fn1[ub]*fn2[u ]*mesq_ubu
	    +fn1[cb]*fn2[c ]*mesq_cbc
	    +fn1[db]*fn2[d ]*mesq_dbd
	    +fn1[sb]*fn2[s ]*mesq_sbs
	    +fn1[bb]*fn2[b ]*mesq_brb;

      //NLL part
      if(opts.order >= 1)
	{
	  QGN_1 = fn1[g]*(fn2[ub]*mesq_uub
			  +fn2[cb]*mesq_ccb
			  +fn2[db]*mesq_ddb
			  +fn2[sb]*mesq_ssb
			  +fn2[bb]*mesq_bbr
			  +fn2[u]*mesq_ubu
			  +fn2[c]*mesq_cbc
			  +fn2[d]*mesq_dbd
			  +fn2[s]*mesq_sbs
			  +fn2[b]*mesq_brb);

	  QGN_2 = fn2[g]*(fn1[u]*mesq_uub
			  +fn1[c]*mesq_ccb
			  +fn1[d]*mesq_ddb
			  +fn1[s]*mesq_ssb
			  +fn1[b]*mesq_bbr
			  +fn1[ub]*mesq_ubu
			  +fn1[cb]*mesq_cbc
			  +fn1[db]*mesq_dbd
			  +fn1[sb]*mesq_sbs
			  +fn1[bb]*mesq_brb);
	}

      //NNLL
      if(opts.order >= 2)
	{
	  GGN = fn1[g]*fn2[g]*
	    (mesq_uub + mesq_ccb
	     + mesq_ddb + mesq_ssb + mesq_bbr
	     + mesq_ubu + mesq_cbc
	     + mesq_dbd + mesq_sbs + mesq_brb);
            
	  //I suspect a mistake here: there are different costh couplings for u-ubar and ubar-u,
	  //the contribution u u -> ub u and u u -> u ub should be separated (need to split hqq_1 and hqq_2)

	  /*
	  QQN = fn1[u]*fn2[u]*mesq_uub
	    +fn1[c]*fn2[c]*mesq_ccb
	    +fn1[d]*fn2[d]*mesq_ddb
	    +fn1[s]*fn2[s]*mesq_ssb
	    +fn1[b]*fn2[b]*mesq_bbr
	    +fn1[ub]*fn2[ub]*mesq_ubu
	    +fn1[cb]*fn2[cb]*mesq_cbc
	    +fn1[db]*fn2[db]*mesq_dbd
	    +fn1[sb]*fn2[sb]*mesq_sbs
	    +fn1[bb]*fn2[bb]*mesq_brb;
	  */

	  QQN_1 = fn1[u]*fn2[u]*mesq_ubu
	    +fn1[c]*fn2[c]*mesq_cbc
	    +fn1[d]*fn2[d]*mesq_dbd
	    +fn1[s]*fn2[s]*mesq_sbs
	    +fn1[b]*fn2[b]*mesq_brb
	    +fn1[ub]*fn2[ub]*mesq_uub
	    +fn1[cb]*fn2[cb]*mesq_ccb
	    +fn1[db]*fn2[db]*mesq_ddb
	    +fn1[sb]*fn2[sb]*mesq_ssb
	    +fn1[bb]*fn2[bb]*mesq_bbr;
	  
	  QQN_2 = fn1[u]*fn2[u]*mesq_uub
	    +fn1[c]*fn2[c]*mesq_ccb
	    +fn1[d]*fn2[d]*mesq_ddb
	    +fn1[s]*fn2[s]*mesq_ssb
	    +fn1[b]*fn2[b]*mesq_bbr
	    +fn1[ub]*fn2[ub]*mesq_ubu
	    +fn1[cb]*fn2[cb]*mesq_cbc
	    +fn1[db]*fn2[db]*mesq_dbd
	    +fn1[sb]*fn2[sb]*mesq_sbs
	    +fn1[bb]*fn2[bb]*mesq_brb;

	  QQPN_1 = fn1[u]*(fn2[db]*mesq_ddb
			   +fn2[sb]*mesq_ssb
			   +fn2[cb]*mesq_ccb
			   +fn2[bb]*mesq_bbr
			   +fn2[d]*mesq_dbd
			   +fn2[s]*mesq_sbs
			   +fn2[c]*mesq_cbc
			   +fn2[b]*mesq_brb)
	    + fn1[d]*(fn2[ub]*mesq_uub
		      +fn2[sb]*mesq_ssb
		      +fn2[cb]*mesq_ccb
		      +fn2[bb]*mesq_bbr
		      +fn2[u]*mesq_ubu
		      +fn2[s]*mesq_sbs
		      +fn2[c]*mesq_cbc
		      +fn2[b]*mesq_brb)
	    + fn1[s]*(fn2[ub]*mesq_uub
		      +fn2[db]*mesq_ddb
		      +fn2[cb]*mesq_ccb
		      +fn2[bb]*mesq_bbr
		      +fn2[u]*mesq_ubu
		      +fn2[d]*mesq_dbd
		      +fn2[c]*mesq_cbc
		      +fn2[b]*mesq_brb)
	    + fn1[c]*(fn2[ub]*mesq_uub
		      +fn2[db]*mesq_ddb
		      +fn2[sb]*mesq_ssb
		      +fn2[bb]*mesq_bbr
		      +fn2[u]*mesq_ubu
		      +fn2[d]*mesq_dbd
		      +fn2[s]*mesq_sbs
		      +fn2[b]*mesq_brb)
	    + fn1[b]*(fn2[ub]*mesq_uub
		      +fn2[db]*mesq_ddb
		      +fn2[sb]*mesq_ssb
		      +fn2[cb]*mesq_ccb
		      +fn2[u]*mesq_ubu
		      +fn2[d]*mesq_dbd
		      +fn2[s]*mesq_sbs
		      +fn2[c]*mesq_cbc)
	    + fn1[ub]*(fn2[db]*mesq_ddb
		       +fn2[sb]*mesq_ssb
		       +fn2[cb]*mesq_ccb
		       +fn2[bb]*mesq_bbr
		       +fn2[d]*mesq_dbd
		       +fn2[s]*mesq_sbs
		       +fn2[c]*mesq_cbc
		       +fn2[b]*mesq_brb)
	    + fn1[db]*(fn2[ub]*mesq_uub
		       +fn2[sb]*mesq_ssb
		       +fn2[cb]*mesq_ccb
		       +fn2[bb]*mesq_bbr
		       +fn2[u]*mesq_ubu
		       +fn2[s]*mesq_sbs
		       +fn2[c]*mesq_cbc
		       +fn2[b]*mesq_brb)
	    + fn1[sb]*(fn2[ub]*mesq_uub
		       +fn2[db]*mesq_ddb
		       +fn2[cb]*mesq_ccb
		       +fn2[bb]*mesq_bbr
		       +fn2[u]*mesq_ubu
		       +fn2[d]*mesq_dbd
		       +fn2[c]*mesq_cbc
		       +fn2[b]*mesq_brb)
	    + fn1[cb]*(fn2[ub]*mesq_uub
		       +fn2[db]*mesq_ddb
		       +fn2[sb]*mesq_ssb
		       +fn2[bb]*mesq_bbr
		       +fn2[u]*mesq_ubu
		       +fn2[d]*mesq_dbd
		       +fn2[s]*mesq_sbs
		       +fn2[b]*mesq_brb)
	    + fn1[bb]*(fn2[ub]*mesq_uub
		       +fn2[db]*mesq_ddb
		       +fn2[sb]*mesq_ssb
		       +fn2[cb]*mesq_ccb
		       +fn2[u]*mesq_ubu
		       +fn2[d]*mesq_dbd
		       +fn2[s]*mesq_sbs
		       +fn2[c]*mesq_cbc);

	  QQPN_2 = fn2[u]*(fn1[d]*mesq_ddb
			   +fn1[s]*mesq_ssb
			   +fn1[c]*mesq_ccb
			   +fn1[b]*mesq_bbr
			   +fn1[db]*mesq_dbd
			   +fn1[sb]*mesq_sbs
			   +fn1[cb]*mesq_cbc
			   +fn1[bb]*mesq_brb)
	    + fn2[d]* (fn1[u]*mesq_uub
		       +fn1[s]*mesq_ssb
		       +fn1[c]*mesq_ccb
		       +fn1[b]*mesq_bbr
		       +fn1[ub]*mesq_ubu
		       +fn1[sb]*mesq_sbs
		       +fn1[cb]*mesq_cbc
		       +fn1[bb]*mesq_brb)
	    + fn2[s]* (fn1[u]*mesq_uub
		       +fn1[d]*mesq_ddb
		       +fn1[c]*mesq_ccb
		       +fn1[b]*mesq_bbr
		       +fn1[ub]*mesq_ubu
		       +fn1[db]*mesq_dbd
		       +fn1[cb]*mesq_cbc
		       +fn1[bb]*mesq_brb)
	    + fn2[c]* (fn1[u]*mesq_uub
		       +fn1[d]*mesq_ddb
		       +fn1[s]*mesq_ssb
		       +fn1[b]*mesq_bbr
		       +fn1[ub]*mesq_ubu
		       +fn1[db]*mesq_dbd
		       +fn1[sb]*mesq_sbs
		       +fn1[bb]*mesq_brb)
	    + fn2[b]* (fn1[u]*mesq_uub
		       +fn1[d]*mesq_ddb
		       +fn1[s]*mesq_ssb
		       +fn1[c]*mesq_ccb
		       +fn1[ub]*mesq_ubu
		       +fn1[db]*mesq_dbd
		       +fn1[sb]*mesq_sbs
		       +fn1[cb]*mesq_cbc)
	    + fn2[ub]*(fn1[d]*mesq_ddb
		       +fn1[s]*mesq_ssb
		       +fn1[c]*mesq_ccb
		       +fn1[b]*mesq_bbr
		       +fn1[db]*mesq_dbd
		       +fn1[sb]*mesq_sbs
		       +fn1[cb]*mesq_cbc
		       +fn1[bb]*mesq_brb)
	    + fn2[db]*(fn1[u]*mesq_uub
		       +fn1[s]*mesq_ssb
		       +fn1[c]*mesq_ccb
		       +fn1[b]*mesq_bbr
		       +fn1[ub]*mesq_ubu
		       +fn1[sb]*mesq_sbs
		       +fn1[cb]*mesq_cbc
		       +fn1[bb]*mesq_brb)
	    + fn2[sb]*(fn1[u]*mesq_uub
		       +fn1[d]*mesq_ddb
		       +fn1[c]*mesq_ccb
		       +fn1[b]*mesq_bbr
		       +fn1[ub]*mesq_ubu
		       +fn1[db]*mesq_dbd
		       +fn1[cb]*mesq_cbc
		       +fn1[bb]*mesq_brb)
	    + fn2[cb]*(fn1[u]*mesq_uub
		       +fn1[d]*mesq_ddb
		       +fn1[s]*mesq_ssb
		       +fn1[b]*mesq_bbr
		       +fn1[ub]*mesq_ubu
		       +fn1[db]*mesq_dbd
		       +fn1[sb]*mesq_sbs
		       +fn1[bb]*mesq_brb)
	    + fn2[bb]*(fn1[u]*mesq_uub
		       +fn1[d]*mesq_ddb
		       +fn1[s]*mesq_ssb
		       +fn1[c]*mesq_ccb
		       +fn1[ub]*mesq_ubu
		       +fn1[db]*mesq_dbd
		       +fn1[sb]*mesq_sbs
		       +fn1[cb]*mesq_cbc);
	}
    }
  if (opts.nproc == 1)
    {
      //DYRES convention
      // bb cb sb db ub g u d s c b
      // -5 -4 -3 -2 -1 0 1 2 3 4 5
      complex <double> mesq_udb = mesq::mesqij_expy[mesq::index(0,i1,i2,sign)];
      complex <double> mesq_dbu = mesq::mesqij_expy[mesq::index(1,i1,i2,sign)];
      complex <double> mesq_usb = mesq::mesqij_expy[mesq::index(2,i1,i2,sign)];
      complex <double> mesq_sbu = mesq::mesqij_expy[mesq::index(3,i1,i2,sign)];
      complex <double> mesq_ubb = mesq::mesqij_expy[mesq::index(4,i1,i2,sign)];
      complex <double> mesq_bbu = mesq::mesqij_expy[mesq::index(5,i1,i2,sign)];
      complex <double> mesq_csb = mesq::mesqij_expy[mesq::index(6,i1,i2,sign)];
      complex <double> mesq_sbc = mesq::mesqij_expy[mesq::index(7,i1,i2,sign)];
      complex <double> mesq_cdb = mesq::mesqij_expy[mesq::index(8,i1,i2,sign)];
      complex <double> mesq_dbc = mesq::mesqij_expy[mesq::index(9,i1,i2,sign)];
      complex <double> mesq_cbb = mesq::mesqij_expy[mesq::index(10,i1,i2,sign)];
      complex <double> mesq_bbc = mesq::mesqij_expy[mesq::index(11,i1,i2,sign)];

      //LL part
      QQBN = fn1[u ]*fn2[db]*mesq_udb
	   + fn1[u ]*fn2[sb]*mesq_usb
	   + fn1[u ]*fn2[bb]*mesq_ubb
	   + fn1[c ]*fn2[db]*mesq_cdb
	   + fn1[c ]*fn2[sb]*mesq_csb
	   + fn1[c ]*fn2[bb]*mesq_cbb
	   + fn1[db]*fn2[u ]*mesq_dbu
	   + fn1[db]*fn2[c ]*mesq_dbc
	   + fn1[sb]*fn2[u ]*mesq_sbu
	   + fn1[sb]*fn2[c ]*mesq_sbc
	   + fn1[bb]*fn2[u ]*mesq_bbu
	   + fn1[bb]*fn2[c ]*mesq_bbc;

      // NLL part
      if(opts.order >= 1)
	{
	  QGN_1 = fn1[g]*(fn2[db]*mesq_udb
			  + fn2[sb]*mesq_usb
			  + fn2[bb]*mesq_ubb
			  + fn2[db]*mesq_cdb
			  + fn2[sb]*mesq_csb
			  + fn2[bb]*mesq_cbb
			  + fn2[u]*mesq_dbu
			  + fn2[c]*mesq_dbc
			  + fn2[u]*mesq_sbu
			  + fn2[c]*mesq_sbc
			  + fn2[u]*mesq_bbu
			  + fn2[c]*mesq_bbc);
      
	  QGN_2 = fn2[g]*(fn1[u]*mesq_udb
			  + fn1[u]*mesq_usb
			  + fn1[u]*mesq_ubb
			  + fn1[c]*mesq_cdb
			  + fn1[c]*mesq_csb
			  + fn1[c]*mesq_cbb
			  + fn1[db]*mesq_dbu
			  + fn1[db]*mesq_dbc
			  + fn1[sb]*mesq_sbu
			  + fn1[sb]*mesq_sbc
			  + fn1[bb]*mesq_bbu
			  + fn1[bb]*mesq_bbc);
	}
      
      //NNLL
      if(opts.order >= 2)
	{
	  GGN = fn1[g]*fn2[g]*
	    (mesq_udb
             +   mesq_usb
             +   mesq_ubb
             +   mesq_cdb
             +   mesq_csb
             +   mesq_cbb
             +   mesq_dbu
             +   mesq_dbc
             +   mesq_sbu
             +   mesq_sbc
             +   mesq_bbu
	     +   mesq_bbc);

	  /*
	  QQN = fn1[u]*fn2[d]*mesq_udb
	    + fn1[u]*fn2[s]*mesq_usb
	    + fn1[u]*fn2[b]*mesq_ubb
	    + fn1[c]*fn2[d]*mesq_cdb
	    + fn1[c]*fn2[s]*mesq_csb
	    + fn1[c]*fn2[b]*mesq_cbb
	    + fn1[db]*fn2[ub]*mesq_dbu
	    + fn1[db]*fn2[cb]*mesq_dbc
	    + fn1[sb]*fn2[ub]*mesq_sbu
	    + fn1[sb]*fn2[cb]*mesq_sbc
	    + fn1[bb]*fn2[ub]*mesq_bbu
	    + fn1[bb]*fn2[cb]*mesq_bbc;
	  */

	  QQN_1 = fn1[ub]*fn2[db]*mesq_udb
	    + fn1[ub]*fn2[sb]*mesq_usb
	    + fn1[ub]*fn2[bb]*mesq_ubb
	    + fn1[cb]*fn2[db]*mesq_cdb
	    + fn1[cb]*fn2[sb]*mesq_csb
	    + fn1[cb]*fn2[bb]*mesq_cbb
	    + fn1[d]*fn2[u]*mesq_dbu
	    + fn1[d]*fn2[c]*mesq_dbc
	    + fn1[s]*fn2[u]*mesq_sbu
	    + fn1[s]*fn2[c]*mesq_sbc
	    + fn1[b]*fn2[u]*mesq_bbu
	    + fn1[b]*fn2[c]*mesq_bbc;
	  
	  QQN_2 = fn1[u]*fn2[d]*mesq_udb
	    + fn1[u]*fn2[s]*mesq_usb
	    + fn1[u]*fn2[b]*mesq_ubb
	    + fn1[c]*fn2[d]*mesq_cdb
	    + fn1[c]*fn2[s]*mesq_csb
	    + fn1[c]*fn2[b]*mesq_cbb
	    + fn1[db]*fn2[ub]*mesq_dbu
	    + fn1[db]*fn2[cb]*mesq_dbc
	    + fn1[sb]*fn2[ub]*mesq_sbu
	    + fn1[sb]*fn2[cb]*mesq_sbc
	    + fn1[bb]*fn2[ub]*mesq_bbu
	    + fn1[bb]*fn2[cb]*mesq_bbc;

	  QQPN_1 =
	    (fn1[u]+fn1[ub])*
	    ( fn2[db]*mesq_cdb
	      + fn2[sb]*mesq_csb
	      + fn2[bb]*mesq_cbb
	      + fn2[u]*mesq_dbu
	      + fn2[c]*mesq_dbc
	      + fn2[u]*mesq_sbu
	      + fn2[c]*mesq_sbc
	      + fn2[u]*mesq_bbu
	      + fn2[c]*mesq_bbc)
	    + (fn1[d]+fn1[db])*
	    ( fn2[db]*mesq_udb
	      + fn2[sb]*mesq_usb
	      + fn2[bb]*mesq_ubb
	      + fn2[db]*mesq_cdb
	      + fn2[sb]*mesq_csb
	      + fn2[bb]*mesq_cbb
	      + fn2[u]*mesq_sbu
	      + fn2[c]*mesq_sbc
	      + fn2[u]*mesq_bbu
	      + fn2[c]*mesq_bbc)
	    + (fn1[s]+fn1[sb])*
	    ( fn2[db]*mesq_udb
	      + fn2[sb]*mesq_usb
	      + fn2[bb]*mesq_ubb
	      + fn2[db]*mesq_cdb
	      + fn2[sb]*mesq_csb
	      + fn2[bb]*mesq_cbb
	      + fn2[u]*mesq_dbu
	      + fn2[c]*mesq_dbc
	      + fn2[u]*mesq_bbu
	      + fn2[c]*mesq_bbc)
	    + (fn1[c]+fn1[cb])*
	    ( fn2[db]*mesq_udb
	      + fn2[sb]*mesq_usb
	      + fn2[bb]*mesq_ubb
	      + fn2[u]*mesq_dbu
	      + fn2[c]*mesq_dbc
	      + fn2[u]*mesq_sbu
	      + fn2[c]*mesq_sbc
	      + fn2[u]*mesq_bbu
	      + fn2[c]*mesq_bbc)
	    + (fn1[b]+fn1[bb])*
	    ( fn2[db]*mesq_udb
	      + fn2[sb]*mesq_usb
	      + fn2[bb]*mesq_ubb
	      + fn2[db]*mesq_cdb
	      + fn2[sb]*mesq_csb
	      + fn2[bb]*mesq_cbb
	      + fn2[u]*mesq_dbu
	      + fn2[c]*mesq_dbc
	      + fn2[u]*mesq_sbu
	      + fn2[c]*mesq_sbc);

	  QQPN_2 = 
	    (fn2[u]+fn2[ub])*
            (fn1[u]*mesq_udb
	     + fn1[u]*mesq_usb
	     + fn1[u]*mesq_ubb
	     + fn1[c]*mesq_cdb
	     + fn1[c]*mesq_csb
	     + fn1[c]*mesq_cbb
	     + fn1[db]*mesq_dbc
	     + fn1[sb]*mesq_sbc
	     + fn1[bb]*mesq_bbc)
	    + (fn2[d]+fn2[db])*
            (fn1[u]*mesq_usb
	     + fn1[u]*mesq_ubb
	     + fn1[c]*mesq_csb
	     + fn1[c]*mesq_cbb
	     + fn1[db]*mesq_dbu
	     + fn1[db]*mesq_dbc
	     + fn1[sb]*mesq_sbu
	     + fn1[sb]*mesq_sbc
	     + fn1[bb]*mesq_bbu
	     + fn1[bb]*mesq_bbc)
	    + (fn2[s]+fn2[sb])*
            (fn1[u]*mesq_udb
	     + fn1[u]*mesq_ubb
	     + fn1[c]*mesq_cdb
	     + fn1[c]*mesq_cbb
	     + fn1[db]*mesq_dbu
	     + fn1[db]*mesq_dbc
	     + fn1[sb]*mesq_sbu
	     + fn1[sb]*mesq_sbc
	     + fn1[bb]*mesq_bbu
	     + fn1[bb]*mesq_bbc)
	    + (fn2[c]+fn2[cb])*
            (fn1[u]*mesq_udb
	     + fn1[u]*mesq_usb
	     + fn1[u]*mesq_ubb
	     + fn1[c]*mesq_cdb
	     + fn1[c]*mesq_csb
	     + fn1[c]*mesq_cbb
	     + fn1[db]*mesq_dbu
	     + fn1[sb]*mesq_sbu
	     + fn1[bb]*mesq_bbu)
	    + (fn2[b]+fn2[bb])*
            (fn1[u]*mesq_udb
	     + fn1[u]*mesq_usb
	     + fn1[c]*mesq_cdb
	     + fn1[c]*mesq_csb
	     + fn1[db]*mesq_dbu
	     + fn1[db]*mesq_dbc
	     + fn1[sb]*mesq_sbu
	     + fn1[sb]*mesq_sbc
	     + fn1[bb]*mesq_bbu
	     + fn1[bb]*mesq_bbc);
	}
    }
    if (opts.nproc == 2)
    {
      //DYRES convention
      // bb cb sb db ub g u d s c b
      // -5 -4 -3 -2 -1 0 1 2 3 4 5
      complex <double> mesq_dub = mesq::mesqij_expy[mesq::index(0,i1,i2,sign)];
      complex <double> mesq_ubd = mesq::mesqij_expy[mesq::index(1,i1,i2,sign)];
      complex <double> mesq_sub = mesq::mesqij_expy[mesq::index(2,i1,i2,sign)];
      complex <double> mesq_ubs = mesq::mesqij_expy[mesq::index(3,i1,i2,sign)];
      complex <double> mesq_bub = mesq::mesqij_expy[mesq::index(4,i1,i2,sign)];
      complex <double> mesq_ubb = mesq::mesqij_expy[mesq::index(5,i1,i2,sign)];
      complex <double> mesq_scb = mesq::mesqij_expy[mesq::index(6,i1,i2,sign)];
      complex <double> mesq_cbs = mesq::mesqij_expy[mesq::index(7,i1,i2,sign)];
      complex <double> mesq_dcb = mesq::mesqij_expy[mesq::index(8,i1,i2,sign)];
      complex <double> mesq_cbd = mesq::mesqij_expy[mesq::index(9,i1,i2,sign)];
      complex <double> mesq_bcb = mesq::mesqij_expy[mesq::index(10,i1,i2,sign)];
      complex <double> mesq_cbb = mesq::mesqij_expy[mesq::index(11,i1,i2,sign)];

      //LL part
      QQBN = fn1[ub]*fn2[d]*mesq_ubd
	+ fn1[ub]*fn2[s]*mesq_ubs
	+ fn1[ub]*fn2[b]*mesq_ubb
	+ fn1[cb]*fn2[d]*mesq_cbd
	+ fn1[cb]*fn2[s]*mesq_cbs
	+ fn1[cb]*fn2[b]*mesq_cbb
	+ fn1[d]*fn2[ub]*mesq_dub
	+ fn1[d]*fn2[cb]*mesq_dcb
	+ fn1[s]*fn2[ub]*mesq_sub
	+ fn1[s]*fn2[cb]*mesq_scb
	+ fn1[b]*fn2[ub]*mesq_bub
	+ fn1[b]*fn2[cb]*mesq_bcb;
      
      // NLL part
      if(opts.order >= 1)
	{
	  QGN_1 = fn1[g]*(fn2[d]*mesq_ubd
			  + fn2[s]*mesq_ubs
			  + fn2[b]*mesq_ubb
			  + fn2[d]*mesq_cbd
			  + fn2[s]*mesq_cbs
			  + fn2[b]*mesq_cbb
			  + fn2[ub]*mesq_dub
			  + fn2[cb]*mesq_dcb
			  + fn2[ub]*mesq_sub
			  + fn2[cb]*mesq_scb
			  + fn2[ub]*mesq_bub
			  + fn2[cb]*mesq_bcb);
      
	  QGN_2 = fn2[g]*(fn1[ub]*mesq_ubd
			  + fn1[ub]*mesq_ubs
			  + fn1[ub]*mesq_ubb
			  + fn1[cb]*mesq_cbd
			  + fn1[cb]*mesq_cbs
			  + fn1[cb]*mesq_cbb
			  + fn1[d]*mesq_dub
			  + fn1[d]*mesq_dcb
			  + fn1[s]*mesq_sub
			  + fn1[s]*mesq_scb
			  + fn1[b]*mesq_bub
			  + fn1[b]*mesq_bcb);
	}
      
      //NNLL
      if(opts.order >= 2)
	{
	  GGN = fn1[g]*fn2[g]*
	    (mesq_ubd
             +   mesq_ubs
             +   mesq_ubb
             +   mesq_cbd
             +   mesq_cbs
             +   mesq_cbb
             +   mesq_dub
             +   mesq_dcb
             +   mesq_sub
             +   mesq_scb
             +   mesq_bub
	     +   mesq_bcb);

	  /*
	  QQN = fn1[ub]*fn2[db]*mesq_ubd
	    + fn1[ub]*fn2[sb]*mesq_ubs
	    + fn1[ub]*fn2[bb]*mesq_ubb
	    + fn1[cb]*fn2[db]*mesq_cbd
	    + fn1[cb]*fn2[sb]*mesq_cbs
	    + fn1[cb]*fn2[bb]*mesq_cbb
	    + fn1[d]*fn2[u]*mesq_dub
	    + fn1[d]*fn2[c]*mesq_dcb
	    + fn1[s]*fn2[u]*mesq_sub
	    + fn1[s]*fn2[c]*mesq_scb
	    + fn1[b]*fn2[u]*mesq_bub
	    + fn1[b]*fn2[c]*mesq_bcb;
	  */

	  QQN_1 = fn1[u]*fn2[d]*mesq_ubd
	    + fn1[u]*fn2[s]*mesq_ubs
	    + fn1[u]*fn2[b]*mesq_ubb
	    + fn1[c]*fn2[d]*mesq_cbd
	    + fn1[c]*fn2[s]*mesq_cbs
	    + fn1[c]*fn2[b]*mesq_cbb
	    + fn1[db]*fn2[ub]*mesq_dub
	    + fn1[db]*fn2[cb]*mesq_dcb
	    + fn1[sb]*fn2[ub]*mesq_sub
	    + fn1[sb]*fn2[cb]*mesq_scb
	    + fn1[bb]*fn2[ub]*mesq_bub
	    + fn1[bb]*fn2[cb]*mesq_bcb;
	  
	  QQN_2 = fn1[ub]*fn2[db]*mesq_ubd
	    + fn1[ub]*fn2[sb]*mesq_ubs
	    + fn1[ub]*fn2[bb]*mesq_ubb
	    + fn1[cb]*fn2[db]*mesq_cbd
	    + fn1[cb]*fn2[sb]*mesq_cbs
	    + fn1[cb]*fn2[bb]*mesq_cbb
	    + fn1[d]*fn2[u]*mesq_dub
	    + fn1[d]*fn2[c]*mesq_dcb
	    + fn1[s]*fn2[u]*mesq_sub
	    + fn1[s]*fn2[c]*mesq_scb
	    + fn1[b]*fn2[u]*mesq_bub
	    + fn1[b]*fn2[c]*mesq_bcb;

	  QQPN_1 =
	    (fn1[ub]+fn1[u])*
	    ( fn2[d]*mesq_cbd
	      + fn2[s]*mesq_cbs
	      + fn2[b]*mesq_cbb
	      + fn2[ub]*mesq_dub
	      + fn2[cb]*mesq_dcb
	      + fn2[ub]*mesq_sub
	      + fn2[cb]*mesq_scb
	      + fn2[ub]*mesq_bub
	      + fn2[cb]*mesq_bcb)
	    + (fn1[db]+fn1[d])*
	    ( fn2[d]*mesq_ubd
	      + fn2[s]*mesq_ubs
	      + fn2[b]*mesq_ubb
	      + fn2[d]*mesq_cbd
	      + fn2[s]*mesq_cbs
	      + fn2[b]*mesq_cbb
	      + fn2[ub]*mesq_sub
	      + fn2[cb]*mesq_scb
	      + fn2[ub]*mesq_bub
	      + fn2[cb]*mesq_bcb)
	    + (fn1[sb]+fn1[s])*
	    ( fn2[d]*mesq_ubd
	      + fn2[s]*mesq_ubs
	      + fn2[b]*mesq_ubb
	      + fn2[d]*mesq_cbd
	      + fn2[s]*mesq_cbs
	      + fn2[b]*mesq_cbb
	      + fn2[ub]*mesq_dub
	      + fn2[cb]*mesq_dcb
	      + fn2[ub]*mesq_bub
	      + fn2[cb]*mesq_bcb)
	    + (fn1[cb]+fn1[c])*
	    ( fn2[d]*mesq_ubd
	      + fn2[s]*mesq_ubs
	      + fn2[b]*mesq_ubb
	      + fn2[ub]*mesq_dub
	      + fn2[cb]*mesq_dcb
	      + fn2[ub]*mesq_sub
	      + fn2[cb]*mesq_scb
	      + fn2[ub]*mesq_bub
	      + fn2[cb]*mesq_bcb)
	    + (fn1[bb]+fn1[b])*
	    ( fn2[d]*mesq_ubd
	      + fn2[s]*mesq_ubs
	      + fn2[b]*mesq_ubb
	      + fn2[d]*mesq_cbd
	      + fn2[s]*mesq_cbs
	      + fn2[b]*mesq_cbb
	      + fn2[ub]*mesq_dub
	      + fn2[cb]*mesq_dcb
	      + fn2[ub]*mesq_sub
	      + fn2[cb]*mesq_scb);

	  QQPN_2 = 
	    (fn2[ub]+fn2[u])*
            (fn1[ub]*mesq_ubd
	     + fn1[ub]*mesq_ubs
	     + fn1[ub]*mesq_ubb
	     + fn1[cb]*mesq_cbd
	     + fn1[cb]*mesq_cbs
	     + fn1[cb]*mesq_cbb
	     + fn1[d]*mesq_dcb
	     + fn1[s]*mesq_scb
	     + fn1[b]*mesq_bcb)
	    + (fn2[db]+fn2[d])*
            (fn1[ub]*mesq_ubs
	     + fn1[ub]*mesq_ubb
	     + fn1[cb]*mesq_cbs
	     + fn1[cb]*mesq_cbb
	     + fn1[d]*mesq_dub
	     + fn1[d]*mesq_dcb
	     + fn1[s]*mesq_sub
	     + fn1[s]*mesq_scb
	     + fn1[b]*mesq_bub
	     + fn1[b]*mesq_bcb)
	    + (fn2[sb]+fn2[s])*
            (fn1[ub]*mesq_ubd
	     + fn1[ub]*mesq_ubb
	     + fn1[cb]*mesq_cbd
	     + fn1[cb]*mesq_cbb
	     + fn1[d]*mesq_dub
	     + fn1[d]*mesq_dcb
	     + fn1[s]*mesq_sub
	     + fn1[s]*mesq_scb
	     + fn1[b]*mesq_bub
	     + fn1[b]*mesq_bcb)
	    + (fn2[cb]+fn2[c])*
            (fn1[ub]*mesq_ubd
	     + fn1[ub]*mesq_ubs
	     + fn1[ub]*mesq_ubb
	     + fn1[cb]*mesq_cbd
	     + fn1[cb]*mesq_cbs
	     + fn1[cb]*mesq_cbb
	     + fn1[d]*mesq_dub
	     + fn1[s]*mesq_sub
	     + fn1[b]*mesq_bub)
	    + (fn2[bb]+fn2[b])*
            (fn1[ub]*mesq_ubd
	     + fn1[ub]*mesq_ubs
	     + fn1[cb]*mesq_cbd
	     + fn1[cb]*mesq_cbs
	     + fn1[d]*mesq_dub
	     + fn1[d]*mesq_dcb
	     + fn1[s]*mesq_sub
	     + fn1[s]*mesq_scb
	     + fn1[b]*mesq_bub
	     + fn1[b]*mesq_bcb);
	}
    }
}

double mellinint::integrand2d(int i1, int i2, int sign)
{
  //cout << "C++ " << i1 << "  " << i2 << "  " << QQBN << endl;
  if (opts.order == 0)
    return real(QQBN);

  int i = hcoefficients::index(i1,i2,sign);
  //cout << "C++ " << i1 << "  " << i2 << "  " << GGN << "  " << hcoefficients::Hgg[i] << endl;
  
  return
    //contribution starting at LL
    real(QQBN*hcoefficients::Hqqb[i])

    //contribution starting at NLL
    + real(QGN_1*hcoefficients::Hqg_1[i]) + real(QGN_2*hcoefficients::Hqg_2[i])
    
    //contributions starting at NNLL
    + real(GGN*hcoefficients::Hgg[i])
    //+ QQN*hcoefficients::Hqq[i]
    + real(QQN_1*hcoefficients::Hqq_1[i]) + real(QQN_2*hcoefficients::Hqq_2[i]) //Bug fix in DYRES (I believe this is more correct, since it accounts for which leg undergoes the q -> qb or qb -> q transformation)
    + real(QQPN_1*hcoefficients::Hqqp_1[i]) + real(QQPN_2*hcoefficients::Hqqp_2[i]);
  
}

double mellinint::integrand1d(int i)
{
  //  cout << i << "  " << QQBN << endl;
  if (opts.order == 0)
    return real(QQBN);

  return
    //contribution starting at LL
    real(QQBN*hcoeff::Hqqb[i])

    //contribution starting at NLL
    + real((QGN_1+QGN_2)*hcoeff::Hqg[i])
    
    //contributions starting at NNLL
    + real(GGN*hcoeff::Hgg[i])
    //+ QQN*hcoeff::Hqq[i]
    + real((QQN_1+QQN_2)*hcoeff::Hqq[i]) //Bug fix in DYRES (I believe this is more correct, since it accounts for which leg undergoes the q -> qb or qb -> q transformation)
    + real((QQPN_1+QQPN_2)*hcoeff::Hqqp[i]);
}

//Failed attempt to perform N to z Mellin inversion before PDF convolution
/*
complex <double> mellinint::integrand()
{
  if (opts.order == 0)
    return QQBN;

  //  cout << "mellinint:integrand()" << QQBN << "  " << hcoeff::Hqqbz << endl;
  return
    //contribution starting at LL
    QQBN*hcoeff::Hqqbz

    //contribution starting at NLL
    + (QGN_1+QGN_2)*hcoeff::Hqgz
    
    //contributions starting at NNLL
    + GGN*hcoeff::Hggz
    + QQN*hcoeff::Hqqz
    //+ (QQN_1+QQN_2)*hcoeff::Hqqz //I believe this is more correct, since it accounts for which leg undergoes the q -> qb or qb -> q transformation
    + (QQPN_1+QQPN_2)*hcoeff::Hqqpz;  
}
*/
