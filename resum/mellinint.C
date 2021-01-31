#include "mellinint.h"
#include "mesq.h"
#include "settings.h"
#include "gaussrules.h"
#include "clenshawcurtisrules.h"
#include "hcoefficients.h"
#include "hcoeff.h"
#include "expc.h"
#include "muf.h"
#include "pdfevol.h"
#include "interface.h"
#include "anomalous.h"
#include "parton.h"
#include "ccoeff.h"
#include "mellinpdf.h"
#include "pegasus.h"
#include "pmom.h"
#include "phasespace.h"

#include <complex>
#include <iostream>
#include <iomanip>
#include <string>

using namespace parton;

//hard-coded Talbot options
double mellinint::Not;
double mellinint::Not_1;
double mellinint::Not_2;
double mellinint::sigma;
double mellinint::lambda;
double mellinint::nu;
double mellinint::alpha;
double mellinint::mu;

bool mellinint::midpoint;
bool mellinint::weideman;

//mellin 1d case
complex <double> *mellinint::wn;
complex <double> *mellinint::Np;
complex <double> *mellinint::Nm;
complex <double> mellinint::CCp; //--> Can drop and put into the Jacobian (jac)
complex <double> mellinint::CCm; //--> Can drop and put into the Jacobian (jac)
double mellinint::cphi;
double mellinint::sphi;

//mellin 2d case
complex <double> *mellinint::Np_1;
complex <double> *mellinint::Nm_1;
complex <double> *mellinint::wn_1;
complex <double> *mellinint::Np_2;
complex <double> *mellinint::Nm_2;
complex <double> *mellinint::wn_2;

complex <double> mellinint::QQBN;
complex <double> mellinint::QGN_1;
complex <double> mellinint::QGN_2;
complex <double> mellinint::QQN;
complex <double> mellinint::QQN_1;
complex <double> mellinint::QQN_2;
complex <double> mellinint::QQPN_1;
complex <double> mellinint::QQPN_2;
complex <double> mellinint::QQBN_nfz;
complex <double> mellinint::QQBPN_1;
complex <double> mellinint::QQBPN_2;
complex <double> mellinint::GGN;
complex <double> mellinint::QBGN_1;
complex <double> mellinint::QBGN_2;
complex <double> mellinint::QPGN_1;
complex <double> mellinint::QPGN_2;
complex <double> mellinint::QBPGN_1;
complex <double> mellinint::QBPGN_2;

int mellinint::mdim;

//imaginary unit
const complex <double> ii = complex <double>(0.,1.);

//fortran interfaces
void mellinint_pdf_mesq_expy_(int& i1, int& i2, int& sign)
{
  mellinint::pdf_mesq_expy(i1-1, i2-1, sign-1);
};
fcomplex mellinint_integrand_(int& i1, int& i2, int& sign)
{
  return fcx(mellinint::integrand2d(i1-1, i2-1, sign-1));
};

void mellinint::allocate()
{
  //wn = new double [mdim];
  if (opts.mellin1d)
    {
      wn = new complex <double> [mdim];
      Np = new complex <double> [mdim];
      Nm = new complex <double> [mdim];
  
      //initialise to 0 --> not needed
      //fill(Np , Np+mdim, cpoint+1.);
      //fill(Nm , Nm+mdim, cpoint+1.);
      //fill(wn , wn+mdim, 0.);
    }
  else
    {
      wn_1 = new complex <double> [mdim];
      Np_1 = new complex <double> [mdim];
      Nm_1 = new complex <double> [mdim];
      wn_2 = new complex <double> [mdim];
      Np_2 = new complex <double> [mdim];
      Nm_2 = new complex <double> [mdim];
  
      //initialise to 0 --> not needed
      //fill(Np_1 , Np_1+mdim, cpoint+1.);
      //fill(Nm_1 , Nm_1+mdim, cpoint+1.);
      //fill(wn_1 , wn_1+mdim, 0.);
      //fill(Np_2 , Np_2+mdim, cpoint+1.);
      //fill(Nm_2 , Nm_2+mdim, cpoint+1.);
      //fill(wn_2 , wn_2+mdim, 0.);
    }
}

void mellinint::release()
{
  if (opts.melup <= 1)
    free();
}

void mellinint::free()
{
  if (opts.mellin1d)
    {
      delete[] wn;
      delete[] Np;
      delete[] Nm;
    }
  else
    {
      delete[] wn_1;
      delete[] Np_1;
      delete[] Nm_1;
      delete[] wn_2;
      delete[] Np_2;
      delete[] Nm_2;
    }
}

void mellinint::initgauss()
{
  //allocate memory
  mdim = opts.mellinintervals*opts.mellinrule;

  //if (!opts.silent) cout << "start initgauss " << endl;
  //fixed talbot of Peter Valko'
  //lambda = 2./5.; (*opts.mellinrule/mlogz);
  //nu = 1.;

  //Modified Talbot of Rizzardi
  //lambda = 4.8; (*1./mlogz);
  //nu = 1.;

  //Empirically good for Talbot
  //sigma = 0.6;
  //lambda = 0.6;
  //nu = 2.5;

  //Weideman (https://arxiv.org/pdf/1304.2505.pdf)
  if (opts.mellintr == 2) // --> google ceres to approximate PDFs
    {
      Not    = 1.; //opts.mellinrule/tmid;  //Global scaling of the contour
      sigma  = -0.6122;               //Offset along the real axis
      mu     = 0.5017;                //Correspond to lambda in the Talbot contour
      alpha  = 0.6407;                //Scale the extrema at minus infinity
      nu     = 0.2645;                //Scaling on the imaginary axis
    }

  //Empirically good for Weideman
  if (opts.mellintr == 1) // --> Laguerre polynomials to approximate PDFs
    {
      Not = 1.;
      sigma = 0.6;
      mu = 0.8;
      alpha = 1.;
      nu = 2.;
    }

  //Not = 1.;
  //sigma = 1.;
  //mu = 0.8;
  //alpha = 1;
  //nu = 4.;
      
  //straight line
  if (opts.mellintr == 0) // --> Gauss Legendre quadrature for the Mellin transform of PDFs
    {
      Not = 1.;
      sigma = opts.cpoint/opts.mellinrule; //opts.cpoint;
      mu = 0.;
      alpha = 1.;
      nu = opts.zmax/M_PI; //-> nu should be a multiple of pi, to cancel oscillations
      //nu = double(opts.ncycle)/double(mdim);
    }
      
  midpoint = true;
  weideman = true;

  if (opts.phi == 0.5)
    {
      cphi = 0;
      sphi = 1;
    }
  else
    {
      cphi = cos(M_PI * opts.phi);
      sphi = sin(M_PI * opts.phi);
    }

  if (opts.mellininv == 0)
    CCp = complex <double> (cphi, sphi);
  else if (opts.mellininv == 1)
    CCp = 1;

  CCm = conj (CCp);
  
  if (opts.melup <= 1)
    allocate();
  if (opts.melup == 0)
    updategauss();
}

//setting m and y to phasespace::m phasespace::y, average of bin boundaries, or average of first-to-last bin boundaries (depending on input file options)
void mellinint::updategauss()
{
  //Update quadrature nodes according to the N/t scaling (N is in the number of nodes, t = -log(x))
  //and recalculate all the Mellin transforms
  //--> should distinguish between N1 and N2 for the double Mellin inversion

  //set up weights and nodes for gaussian quadrature
  //cpoint is the starting point on the real axis for the positive and negative part of the integration path
  //opts.phi*M_PI is the angle in the complex plane of the positive part of the integration path
  //opts.phi = 1/2 -> straight line; opts.phi = 3/4 -> Pegasus default

  double zmin = 0;

  //mellin1d case
  if (opts.mellin1d)
    {
      //zmax is the upper limit for the mellin integration in the complex plane (z)
      //Above 50 the integral becomes unstable, especially when m_ll << mur,
      //due to large logs in the Sudakov, which can be reduced with smaller blim.
      //The issue can be avoided by using dynamicscale
      //Also, larger values of zmax requires more support points for the
      //Mellin inverse transform, i.e. higher mellinintervals or mellinrule.

      //zmax is set to opts.ncycle * M_PI,
      //and scaled by N/t

      //cpoint is the intersection with the real axis of the integration contour
      //cpoint is set to opts.cpoint scaled by 1/t, and shifted to the right by opts.cshift

      double m0;
      //Compute once for all
      if (opts.melup == 0)
	//m0 = (bins.mbins.front() + bins.mbins.back())/2.; //arithmetic mean
        m0 = sqrt(bins.mbins.front() * bins.mbins.back()); //geometric mean
      
      //Recompute in each bin
      if (opts.melup == 1)
	//m0 = (phasespace::mmin + phasespace::mmax)/2.; //arithmetic mean
	m0 = sqrt(phasespace::mmin * phasespace::mmax); //geometric mean
      
      //Recompute in each phasespace point
      if (opts.melup == 2)
	m0 = phasespace::m;

      //Recompute in each phasespace point
      double mlogz = -log(m0/opts.sroot);
      
      double cpoint = opts.cpoint       / mlogz;
      //double zmax   = opts.zmax   *mdim / mlogz;
      double zmax   = opts.ncycle *M_PI / mlogz;
      cpoint += opts.cshift;

      //Gauss-Legendre quadrature along a linear contour
      if (opts.mellininv == 0)
	{
	  // positive branch      
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
		  Np[j+i*opts.mellinrule]=complex <double> (cpoint+cphi*t+1.,sphi*t);
		  wn[j+i*opts.mellinrule]=gr::www[opts.mellinrule-1][j]*m*jac;
		  //wn[j+i*opts.mellinrule]=cc::www[opts.mellinrule-1][j]*m*jac;
		  //cout << setprecision(16) <<  t << " " << Np[j+i*opts.mellinrule] << "  " << wn[j+i*opts.mellinrule] << endl;
		}
	    }      
	}
      //Talbot contour
      else if (opts.mellininv == 1)
	{
	  //fixed talbot of Peter Valko'
	  //Not = mdim/mlogz;
	  
	  //Modified Talbot of Rizzardi
	  //Not = 1./mlogz;
	  
	  //Empirically good for Talbot
	  //Not = 1.;
	  
	  //Weideman (https://arxiv.org/pdf/1304.2505.pdf)
	  if (opts.mellintr == 2)
	    Not    = mdim / mlogz;          //Global scaling of the contour
	  
	  //Empirically good for Weideman
	  if (opts.mellintr == 1)
	    Not = 5. / mlogz;
	  
	  //straight line
	  if (opts.mellintr == 0)
	    {
	      Not    = mdim / mlogz;          //Global scaling of the contour
	      //Not    = 1. / mlogz;          //Global scaling of the contour
	      //Set mu to a value that curves the contour so that the last point has real part equal to c (and never becomes negative)
	      //double c = 1.;                
	      /*sigma + mu (1-2N) = c --> as in https://people.math.ethz.ch/~hiptmair/Seminars/CONVQUAD/Articles/WEI06.pdf */
	      /*sigma + mu (1-N) + 1 = c --> Empirically correct (for the trapezoidal rule) */
	      //mu = (sigma-(c-1.)/Not)/(mdim-1.);
	      //cout << Not*(sigma + mu* (1.-mdim))+1. << endl;
	    }
	  
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

	      Np[j] = s+1.+opts.cshift;
	      //cout << j << "  " << Np[j] << endl;
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
	}
    }
  else
    //mellin 2d case
    {
      //Scale zmax by N/t and cpoint by 1/t (mellin1d case)
      double mlogz_1,mlogz_2;
      
      //Recompute in each bin
      double m0,y0;

      //Compute once for all
      if (opts.melup == 0)
	{
	  //m0 = (bins.mbins.front() + bins.mbins.back())/2.; //arithmetic mean
	  //m0 = sqrt(bins.mbins.front()*bins.mbins.back()); //geometric mean
	  m0 = bins.mbins.front(); //take the minimum to avoid z < 0
	  double ylim = log(opts.sroot/bins.mbins.back());  //here take bins.mbins.back() so that this is the smallest ylim --> avoid z < 0
	  double ymn = min(max(-ylim, bins.ybins.front()),ylim);
	  double ymx = max(min(ylim, bins.ybins.back()),-ylim);
	  y0 = (ymn+ymx)/2.;
	}

      //Recompute in each bin
      if (opts.melup == 1)
	{
	  //m0 = (phasespace::mmin + phasespace::mmax)/2.;    //arithmetic mean
	  //m0 = sqrt(phasespace::mmin * phasespace::mmax); //geometric mean
	  m0 = phasespace::mmin; //take the minimum to avoid log(z) < 0 (--> x > 1)
	  double ylim = log(opts.sroot/phasespace::mmax); //here take mmax so that this is the smallest ylim --> avoid log(z) < 0 (--> x > 1)
	  double ymn = min(max(-ylim, phasespace::ymin),ylim);
	  double ymx = max(min(ylim, phasespace::ymax),-ylim);
	  y0 = (ymn+ymx)/2.;
	}
      
      //Recompute in each phasespace point
      if (opts.melup == 2)
	{
	  m0 = phasespace::m;

	  double ylim = 0.5*log(pow(opts.sroot,2)/phasespace::m2);
	  double ymn = min(max(-ylim, phasespace::ymin),ylim);
	  double ymx = max(min(ylim, phasespace::ymax),-ylim);
	  //if (y-differential mode...)
	  //y0 = phasespace::y;
	  //else //y-integrated mode
	  //{
	  y0 = (ymn+ymx)/2.;
	  if (phasespace::ymax >= ylim && phasespace::ymin > -ylim)
	    y0 = ymn;
	  if (phasespace::ymin <= -ylim && phasespace::ymax < ylim)
	    y0 = ymx;
	  //}
	  
//      if (ymx >= ylim-1e-4)
//	{
//	  mlogz_1 = -(log(phasespace::m/opts.sroot)+ymn);
//	  mlogz_2 = -(log(phasespace::m/opts.sroot)-ymn);
//	}
//      else
//	{
//	  mlogz_1 = 2./ (1./(-(log(phasespace::m/opts.sroot)+ymn)) + 1./(-(log(phasespace::m/opts.sroot)+ymx))    );
//	  mlogz_2 = 2./ (1./(-(log(phasespace::m/opts.sroot)-ymn)) + 1./(-(log(phasespace::m/opts.sroot)-ymx))    );
//	}

	}
      mlogz_1 = -(log(m0/opts.sroot)+y0);
      mlogz_2 = -(log(m0/opts.sroot)-y0);
      //cout << " z1 " << m0/opts.sroot*exp(y0) << endl;
      //cout << " z2 " << m0/opts.sroot*exp(-y0) << endl;

      //mlogz_1 = -log(m0/opts.sroot);
      //mlogz_2 = -log(m0/opts.sroot);

      
      //cout << ymx << "  " << ylim- (1e-4) << "  " << (-mlogz_1-log(phasespace::m/opts.sroot)) << endl;
      //cout << " x1 " << phasespace::m/opts.sroot*exp(phasespace::y) << " mlogz1 " << mlogz_1 << endl;
      //cout << " x2 " << phasespace::m/opts.sroot*exp(-phasespace::y) << " mlogz2 " << mlogz_2 << endl;
      
      
      double cpoint_1 = opts.cpoint       / mlogz_1;
      //double zmax_1   = opts.zmax   *mdim / mlogz_1;
      double zmax_1   = opts.ncycle *M_PI / mlogz_1;

      double cpoint_2 = opts.cpoint       / mlogz_2;
      //double zmax_2   = opts.zmax   *mdim / mlogz_2;
      double zmax_2   = opts.ncycle *M_PI / mlogz_2;

      //Patch to fix issues with NNPDF, likely due to the position of the PDF poles in N-space
      cpoint_1 += opts.cshift;
      cpoint_2 += opts.cshift;
      
      //Gauss-Legendre quadrature along a linear contour
      if (opts.mellininv == 0)
	{
	  // positive branch      
	  for (int i=0; i < opts.mellinintervals; i++)
	    {
	      double a = 0.+(1.-0.)*i/opts.mellinintervals;
	      double b = 0.+(1.-0.)*(i+1)/opts.mellinintervals;
	      double c = 0.5*(a+b);
	      double m = 0.5*(b-a);
	      for (int j=0; j < opts.mellinrule; j++)
		{
		  double x = c+m*gr::xxx[opts.mellinrule-1][j];

		  double t_1 = zmin+(zmax_1-zmin)*x;
		  complex <double> jac_1 = (zmax_1-zmin);//*CCp;
		  Np_1[j+i*opts.mellinrule]=cpoint_1+1.+CCp*t_1;//  complex <double> (cpoint_1+cphi*t_1+1.,sphi*t_1);
		  wn_1[j+i*opts.mellinrule]=gr::www[opts.mellinrule-1][j]*m*jac_1;
		  //cout << setprecision(16) <<  t_1 << " " << Np_1[j+i*opts.mellinrule] << "  " << wn_1[j+i*opts.mellinrule] << endl;

		  double t_2 = zmin+(zmax_2-zmin)*x;
		  complex <double> jac_2 = (zmax_2-zmin);//*CCp;
		  Np_2[j+i*opts.mellinrule]=cpoint_2+1.+CCp*t_2; //complex <double> (cpoint_2+cphi*t_2+1.,sphi*t_2);
		  wn_2[j+i*opts.mellinrule]=gr::www[opts.mellinrule-1][j]*m*jac_2;
		  //cout << setprecision(16) <<  t_2 << " " << Np_2[j+i*opts.mellinrule] << "  " << wn_2[j+i*opts.mellinrule] << endl;

		  //double t = zmin+(zmax_1-zmin)*x;
		  //double jac = zmax_1-zmin;
		  //Np[j+i*opts.mellinrule]=complex <double> (cpoint_1+cphi*t+1.,sphi*t);
		  //wn[j+i*opts.mellinrule]=gr::www[opts.mellinrule-1][j]*m*jac;
		}
	    }
	}

      //Talbot contour
      else if (opts.mellininv == 1)
	{
	  double mu_1 = mu;
	  double mu_2 = mu;
	  
	  //fixed talbot of Peter Valko'
	  //Not = mdim/mlogz;
	  
	  //Modified Talbot of Rizzardi
	  //Not = 1./mlogz;
	  
	  //Empirically good for Talbot
	  //Not = 1.;
	  
	  //Weideman (https://arxiv.org/pdf/1304.2505.pdf)
	  if (opts.mellintr == 2)
	    {
	      Not_1    = mdim / mlogz_1;          //Global scaling of the contour
	      Not_2    = mdim / mlogz_2;          //Global scaling of the contour
	    }
	  
	  //Empirically good for Weideman
	  if (opts.mellintr == 1)
	    {
	      Not_1 = 5. / mlogz_1;
	      Not_2 = 5. / mlogz_2;
	    }
	  
	  //straight line
	  if (opts.mellintr == 0)
	    {
	      Not_1    = mdim / mlogz_1;          //Global scaling of the contour
	      Not_2    = mdim / mlogz_2;          //Global scaling of the contour

	      //Not_1    = 1. / mlogz_1;          //Global scaling of the contour
	      //Not_2    = 1. / mlogz_2;          //Global scaling of the contour

	      //Set mu to a value that curves the contour so that the last point has real part equal to c (and never becomes negative)
	      //double c = 1.;
	      /*sigma + mu (1-2N) = c --> as in https://people.math.ethz.ch/~hiptmair/Seminars/CONVQUAD/Articles/WEI06.pdf */
	      /*sigma + mu (1-N) + 1 = c --> Empirically correct (for the trapezoidal rule) */
	      //mu_1 = (sigma-(c-1.)/Not_1)/(mdim-1.);
	      //mu_2 = (sigma-(c-1.)/Not_2)/(mdim-1.);
	      //cout << Not*(sigma + mu* (1.-mdim))+1. << endl;

	      //Not_1 = 2;
	      //Not_2 = 1;
	      //sigma = 1.;
	      //nu = 22./M_PI;

	      //Not = 1.;
	      //sigma = opts.cpoint/opts.mellinrule;
	      //mu = 0.;
	      //alpha = 1.;
	      //nu = opts.zmax/M_PI;
	      //
	      //cout << endl;
	      //cout << sigma * Not_1 << endl;
	      //cout << sigma * Not_2 << endl;
	    }
	  
	  for (int j=0; j < opts.mellinrule; j++)
	    {
	      double x;
	      if (midpoint)
		x = (double(j)+0.5)/opts.mellinrule;   //midpoint
	      else
		x = double(j)/opts.mellinrule;         //trapezoidal
	      
	      double thmax = M_PI;
	      double theta = thmax*x;
	      
	      complex <double> s_1,s_2, jac_1,jac_2;
	      if (weideman)
		{
		  //Weideman contour
		  s_1 = Not_1*(sigma+(theta==0?mu_1/alpha:mu_1*theta/tan(alpha*theta)+theta*nu*ii));
		  jac_1 = thmax * Not_1 * (ii * nu + (theta==0?0 : mu_1 * (1./tan(alpha*theta) - alpha*theta*(1. + 1./tan(alpha*theta)/tan(alpha*theta)))));
		  s_2 = Not_2*(sigma+(theta==0?mu_2/alpha:mu_2*theta/tan(alpha*theta)+theta*nu*ii));
		  jac_2 = thmax * Not_2 * (ii * nu + (theta==0?0 : mu_2 * (1./tan(alpha*theta) - alpha*theta*(1. + 1./tan(alpha*theta)/tan(alpha*theta)))));
		}
	      else
		{
		  //Talbot contour
		  s_1 = sigma+(theta==0?lambda:lambda*theta*(1./tan(theta)+nu*ii));
		  jac_1 = thmax * (theta==0? ii*lambda*nu : ii * lambda *( nu + ii * (theta + (theta /tan(theta) - 1.) / tan(theta))));
		  s_2 = sigma+(theta==0?lambda:lambda*theta*(1./tan(theta)+nu*ii));
		  jac_2 = thmax * (theta==0? ii*lambda*nu : ii * lambda *( nu + ii * (theta + (theta /tan(theta) - 1.) / tan(theta))));
		}

	      Np_1[j] = s_1+1.+opts.cshift;
	      Np_2[j] = s_2+1.+opts.cshift;
	      //cout << j << " Np_1 " << Np_1[j] << " Np_2 " << Np_2[j] << endl;
	      if (midpoint)
		{
		  wn_1[j] = jac_1/double(opts.mellinrule);                                //midpoint
		  wn_2[j] = jac_2/double(opts.mellinrule);                                //midpoint
		}
	      else
		{
		  wn_1[j] = (j==0 ? jac_1/2.:jac_1)/double(opts.mellinrule);                //trapezoidal
		  wn_2[j] = (j==0 ? jac_2/2.:jac_2)/double(opts.mellinrule);                //trapezoidal
		}
	  
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
	}
    }

  //Compute negative branch by complex conjugation of the positive branch
  if (opts.mellin1d)
    for (int i=0; i < opts.mellinintervals; i++)
      for (int j=0; j < opts.mellinrule; j ++)
	Nm[j+i*opts.mellinrule] = conj(Np[j+i*opts.mellinrule]);
  else
    for (int i=0; i < opts.mellinintervals; i++)
      for (int j=0; j < opts.mellinrule; j ++)
	{
	  Nm_1[j+i*opts.mellinrule] = conj(Np_1[j+i*opts.mellinrule]);
	  Nm_2[j+i*opts.mellinrule] = conj(Np_2[j+i*opts.mellinrule]);
	}
  
  
  //cout << endl;
  //cout << " Mellin inversion support points" << endl;
  //for (int i=0; i < opts.mellinintervals; i++)
  //  for (int j=0; j < opts.mellinrule; j ++)
  //    {
  //	if (opts.mellin1d)
  //	  cout << j << "  " << Np[j] << endl;
  //	else
  //	  {
  //	    cout << j << "  " << Np_1[j] << "  " << Np_2[j] << endl;
  //	  }
  //    }
}

void mellinint::update()
{
  //Update quadrature nodes according to the N/t scaling (N is in the number of nodes, t = -log(x))
  //cout << "Update mellin points... " << flush;

  clock_t begin_time_tot, end_time_tot;
  begin_time_tot = clock();  
  
  updategauss();

  //Recalculate all the Mellin transforms, i.e. update all other initialisations which depends on the Mellin quadrature nodes
  clock_t begin_time, end_time;

  /*
  //calculate anomalous dimensions, C1, C2 and gamma coefficients
  begin_time = clock();  
  anomalous::calc();
  end_time = clock();  
  //cout << "anomalous done in "  << float(end_time - begin_time)/CLOCKS_PER_SEC*1000. << "ms" << endl;
  */
  
  //calculate, C1, C2, and C3 coefficients
  begin_time = clock();
  if (opts.mellin1d)
    ccoeff::calc1d();
  else
    ccoeff::calc2d();
  end_time = clock();  
  //cout << "ccoeff done in "  << float(end_time - begin_time)/CLOCKS_PER_SEC*1000. << "ms" << endl;

  /*
  //transform the PDF from x- to N-space at the factorisation scale
  begin_time = clock();  
  //mellinpdf::update_mellin();
  pdfevol::allocate();
  pdfevol::update();
  pdfevol::free();
  end_time = clock();  
  //cout << "pdfevol done in "  << float(end_time - begin_time)/CLOCKS_PER_SEC*1000. << "ms" << endl;
  */
  
  begin_time = clock();  
  pegasus::calc_mellin();
  //pegasus::init_pdf(); //pegasus::update();
  end_time = clock();  
  //cout << "pegasus done in "  << float(end_time - begin_time)/CLOCKS_PER_SEC*1000. << "ms" << endl;

  begin_time = clock();  
  pmom::calc();                          //calculate gamma from Pegasus (need to be called after mellinint::initgauss())
  end_time = clock();  
  //cout << "pmom done in "  << float(end_time - begin_time)/CLOCKS_PER_SEC*1000. << "ms" << endl;

  end_time_tot = clock();  
  //cout << "Done in "  << float(end_time_tot - begin_time_tot)/CLOCKS_PER_SEC*1000. << "ms" << endl;
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

void mellinint::reset()
{
  //LL
  QQBN=0;

  //NLL
  QGN_1=0;
  QGN_2=0;

  //NNLL
  QQBN_nfz=0;
  GGN=0;
  QQN_1=0;
  QQN_2=0;
  QQPN_1=0;
  QQPN_2=0;
}

void mellinint::pdf_mesq_expy(int i1, int i2, int sign)
{
//  //QQBN=0;
//  if (opts.order < 1)
//    {
//      QGN_1=0;
//      QGN_2=0;
//    }
//      
//  if (opts.order < 2)
//    {
//      GGN=0;
//      //QQN=0;
//      QQN_1=0;
//      QQN_2=0;
//      QQPN_1=0;
//      QQPN_2=0;
//    }

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
      if (opts.order >= 1)
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
      if (opts.order >= 2)
	{
	  GGN = fn1[g]*fn2[g]*
	    (mesq_uub + mesq_ccb
	     + mesq_ddb + mesq_ssb + mesq_bbr
	     + mesq_ubu + mesq_cbc
	     + mesq_dbd + mesq_sbs + mesq_brb);
            
	  //I suspect a mistake here: there are different costh couplings for u-ubar and ubar-u,
	  //the contribution u u -> ub u and u u -> u ub should be separated (need to split hqq_1 and hqq_2)

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

	  QQPN_1 =
	    +fn1[u]*fn2[b]*mesq_brb
	    +fn1[d]*fn2[b]*mesq_brb
	    +fn1[s]*fn2[b]*mesq_brb
	    +fn1[c]*fn2[b]*mesq_brb
	    +fn1[u]*fn2[c]*mesq_cbc
	    +fn1[d]*fn2[c]*mesq_cbc
	    +fn1[s]*fn2[c]*mesq_cbc
	    +fn1[b]*fn2[c]*mesq_cbc
	    +fn1[u]*fn2[s]*mesq_sbs
	    +fn1[d]*fn2[s]*mesq_sbs
	    +fn1[c]*fn2[s]*mesq_sbs
	    +fn1[b]*fn2[s]*mesq_sbs
	    +fn1[u]*fn2[d]*mesq_dbd
	    +fn1[s]*fn2[d]*mesq_dbd
	    +fn1[c]*fn2[d]*mesq_dbd
	    +fn1[b]*fn2[d]*mesq_dbd
	    +fn1[d]*fn2[u]*mesq_ubu
	    +fn1[s]*fn2[u]*mesq_ubu
	    +fn1[c]*fn2[u]*mesq_ubu
	    +fn1[b]*fn2[u]*mesq_ubu
	    +fn1[bb]*fn2[ub]*mesq_uub
	    +fn1[cb]*fn2[ub]*mesq_uub
	    +fn1[sb]*fn2[ub]*mesq_uub
	    +fn1[db]*fn2[ub]*mesq_uub
	    +fn1[bb]*fn2[db]*mesq_ddb
	    +fn1[cb]*fn2[db]*mesq_ddb
	    +fn1[sb]*fn2[db]*mesq_ddb
	    +fn1[ub]*fn2[db]*mesq_ddb
	    +fn1[bb]*fn2[sb]*mesq_ssb
	    +fn1[cb]*fn2[sb]*mesq_ssb
	    +fn1[db]*fn2[sb]*mesq_ssb
	    +fn1[ub]*fn2[sb]*mesq_ssb
	    +fn1[bb]*fn2[cb]*mesq_ccb
	    +fn1[sb]*fn2[cb]*mesq_ccb
	    +fn1[db]*fn2[cb]*mesq_ccb
	    +fn1[ub]*fn2[cb]*mesq_ccb
	    +fn1[cb]*fn2[bb]*mesq_bbr
	    +fn1[sb]*fn2[bb]*mesq_bbr
	    +fn1[db]*fn2[bb]*mesq_bbr
	    +fn1[ub]*fn2[bb]*mesq_bbr;

	  QQPN_2 =
	    +fn1[bb]*fn2[u]*mesq_brb
	    +fn1[bb]*fn2[d]*mesq_brb
	    +fn1[bb]*fn2[s]*mesq_brb
	    +fn1[bb]*fn2[c]*mesq_brb
	    +fn1[cb]*fn2[u]*mesq_cbc
	    +fn1[cb]*fn2[d]*mesq_cbc
	    +fn1[cb]*fn2[s]*mesq_cbc
	    +fn1[cb]*fn2[b]*mesq_cbc
	    +fn1[sb]*fn2[u]*mesq_sbs
	    +fn1[sb]*fn2[d]*mesq_sbs
	    +fn1[sb]*fn2[c]*mesq_sbs
	    +fn1[sb]*fn2[b]*mesq_sbs
	    +fn1[db]*fn2[u]*mesq_dbd
	    +fn1[db]*fn2[s]*mesq_dbd
	    +fn1[db]*fn2[c]*mesq_dbd
	    +fn1[db]*fn2[b]*mesq_dbd
	    +fn1[ub]*fn2[d]*mesq_ubu
	    +fn1[ub]*fn2[s]*mesq_ubu
	    +fn1[ub]*fn2[c]*mesq_ubu
	    +fn1[ub]*fn2[b]*mesq_ubu
	    +fn1[u]*fn2[bb]*mesq_uub
	    +fn1[u]*fn2[cb]*mesq_uub
	    +fn1[u]*fn2[sb]*mesq_uub
	    +fn1[u]*fn2[db]*mesq_uub
	    +fn1[d]*fn2[bb]*mesq_ddb
	    +fn1[d]*fn2[cb]*mesq_ddb
	    +fn1[d]*fn2[sb]*mesq_ddb
	    +fn1[d]*fn2[ub]*mesq_ddb
	    +fn1[s]*fn2[bb]*mesq_ssb
	    +fn1[s]*fn2[cb]*mesq_ssb
	    +fn1[s]*fn2[db]*mesq_ssb
	    +fn1[s]*fn2[ub]*mesq_ssb
	    +fn1[c]*fn2[bb]*mesq_ccb
	    +fn1[c]*fn2[sb]*mesq_ccb
	    +fn1[c]*fn2[db]*mesq_ccb
	    +fn1[c]*fn2[ub]*mesq_ccb
	    +fn1[b]*fn2[cb]*mesq_bbr
	    +fn1[b]*fn2[sb]*mesq_bbr
	    +fn1[b]*fn2[db]*mesq_bbr
	    +fn1[b]*fn2[ub]*mesq_bbr;

	  QQBPN_1 =
	    +fn1[cb]*fn2[b]*mesq_brb
	    +fn1[sb]*fn2[b]*mesq_brb
	    +fn1[db]*fn2[b]*mesq_brb
	    +fn1[ub]*fn2[b]*mesq_brb
	    +fn1[bb]*fn2[c]*mesq_cbc
	    +fn1[sb]*fn2[c]*mesq_cbc
	    +fn1[db]*fn2[c]*mesq_cbc
	    +fn1[ub]*fn2[c]*mesq_cbc
	    +fn1[bb]*fn2[s]*mesq_sbs
	    +fn1[cb]*fn2[s]*mesq_sbs
	    +fn1[db]*fn2[s]*mesq_sbs
	    +fn1[ub]*fn2[s]*mesq_sbs
	    +fn1[bb]*fn2[d]*mesq_dbd
	    +fn1[cb]*fn2[d]*mesq_dbd
	    +fn1[sb]*fn2[d]*mesq_dbd
	    +fn1[ub]*fn2[d]*mesq_dbd
	    +fn1[bb]*fn2[u]*mesq_ubu
	    +fn1[cb]*fn2[u]*mesq_ubu
	    +fn1[sb]*fn2[u]*mesq_ubu
	    +fn1[db]*fn2[u]*mesq_ubu
	    +fn1[d]*fn2[ub]*mesq_uub
	    +fn1[s]*fn2[ub]*mesq_uub
	    +fn1[c]*fn2[ub]*mesq_uub
	    +fn1[b]*fn2[ub]*mesq_uub
	    +fn1[u]*fn2[db]*mesq_ddb
	    +fn1[s]*fn2[db]*mesq_ddb
	    +fn1[c]*fn2[db]*mesq_ddb
	    +fn1[b]*fn2[db]*mesq_ddb
	    +fn1[u]*fn2[sb]*mesq_ssb
	    +fn1[d]*fn2[sb]*mesq_ssb
	    +fn1[c]*fn2[sb]*mesq_ssb
	    +fn1[b]*fn2[sb]*mesq_ssb
	    +fn1[u]*fn2[cb]*mesq_ccb
	    +fn1[d]*fn2[cb]*mesq_ccb
	    +fn1[s]*fn2[cb]*mesq_ccb
	    +fn1[b]*fn2[cb]*mesq_ccb
	    +fn1[u]*fn2[bb]*mesq_bbr
	    +fn1[d]*fn2[bb]*mesq_bbr
	    +fn1[s]*fn2[bb]*mesq_bbr
	    +fn1[c]*fn2[bb]*mesq_bbr;

	  QQBPN_2 =
	    +fn1[bb]*fn2[cb]*mesq_brb
	    +fn1[bb]*fn2[sb]*mesq_brb
	    +fn1[bb]*fn2[db]*mesq_brb
	    +fn1[bb]*fn2[ub]*mesq_brb
	    +fn1[cb]*fn2[bb]*mesq_cbc
	    +fn1[cb]*fn2[sb]*mesq_cbc
	    +fn1[cb]*fn2[db]*mesq_cbc
	    +fn1[cb]*fn2[ub]*mesq_cbc
	    +fn1[sb]*fn2[bb]*mesq_sbs
	    +fn1[sb]*fn2[cb]*mesq_sbs
	    +fn1[sb]*fn2[db]*mesq_sbs
	    +fn1[sb]*fn2[ub]*mesq_sbs
	    +fn1[db]*fn2[bb]*mesq_dbd
	    +fn1[db]*fn2[cb]*mesq_dbd
	    +fn1[db]*fn2[sb]*mesq_dbd
	    +fn1[db]*fn2[ub]*mesq_dbd
	    +fn1[ub]*fn2[bb]*mesq_ubu
	    +fn1[ub]*fn2[cb]*mesq_ubu
	    +fn1[ub]*fn2[sb]*mesq_ubu
	    +fn1[ub]*fn2[db]*mesq_ubu
	    +fn1[u]*fn2[d]*mesq_uub
	    +fn1[u]*fn2[s]*mesq_uub
	    +fn1[u]*fn2[c]*mesq_uub
	    +fn1[u]*fn2[b]*mesq_uub
	    +fn1[d]*fn2[u]*mesq_ddb
	    +fn1[d]*fn2[s]*mesq_ddb
	    +fn1[d]*fn2[c]*mesq_ddb
	    +fn1[d]*fn2[b]*mesq_ddb
	    +fn1[s]*fn2[u]*mesq_ssb
	    +fn1[s]*fn2[d]*mesq_ssb
	    +fn1[s]*fn2[c]*mesq_ssb
	    +fn1[s]*fn2[b]*mesq_ssb
	    +fn1[c]*fn2[u]*mesq_ccb
	    +fn1[c]*fn2[d]*mesq_ccb
	    +fn1[c]*fn2[s]*mesq_ccb
	    +fn1[c]*fn2[b]*mesq_ccb
	    +fn1[b]*fn2[u]*mesq_bbr
	    +fn1[b]*fn2[d]*mesq_bbr
	    +fn1[b]*fn2[s]*mesq_bbr
	    +fn1[b]*fn2[c]*mesq_bbr;
	  
//	  QQPN_1 = fn1[u]*(fn2[db]*mesq_ddb
//			   +fn2[sb]*mesq_ssb
//			   +fn2[cb]*mesq_ccb
//			   +fn2[bb]*mesq_bbr
//			   +fn2[d]*mesq_dbd
//			   +fn2[s]*mesq_sbs
//			   +fn2[c]*mesq_cbc
//			   +fn2[b]*mesq_brb)
//	    + fn1[d]*(fn2[ub]*mesq_uub
//		      +fn2[sb]*mesq_ssb
//		      +fn2[cb]*mesq_ccb
//		      +fn2[bb]*mesq_bbr
//		      +fn2[u]*mesq_ubu
//		      +fn2[s]*mesq_sbs
//		      +fn2[c]*mesq_cbc
//		      +fn2[b]*mesq_brb)
//	    + fn1[s]*(fn2[ub]*mesq_uub
//		      +fn2[db]*mesq_ddb
//		      +fn2[cb]*mesq_ccb
//		      +fn2[bb]*mesq_bbr
//		      +fn2[u]*mesq_ubu
//		      +fn2[d]*mesq_dbd
//		      +fn2[c]*mesq_cbc
//		      +fn2[b]*mesq_brb)
//	    + fn1[c]*(fn2[ub]*mesq_uub
//		      +fn2[db]*mesq_ddb
//		      +fn2[sb]*mesq_ssb
//		      +fn2[bb]*mesq_bbr
//		      +fn2[u]*mesq_ubu
//		      +fn2[d]*mesq_dbd
//		      +fn2[s]*mesq_sbs
//		      +fn2[b]*mesq_brb)
//	    + fn1[b]*(fn2[ub]*mesq_uub
//		      +fn2[db]*mesq_ddb
//		      +fn2[sb]*mesq_ssb
//		      +fn2[cb]*mesq_ccb
//		      +fn2[u]*mesq_ubu
//		      +fn2[d]*mesq_dbd
//		      +fn2[s]*mesq_sbs
//		      +fn2[c]*mesq_cbc)
//	    + fn1[ub]*(fn2[db]*mesq_ddb
//		       +fn2[sb]*mesq_ssb
//		       +fn2[cb]*mesq_ccb
//		       +fn2[bb]*mesq_bbr
//		       +fn2[d]*mesq_dbd
//		       +fn2[s]*mesq_sbs
//		       +fn2[c]*mesq_cbc
//		       +fn2[b]*mesq_brb)
//	    + fn1[db]*(fn2[ub]*mesq_uub
//		       +fn2[sb]*mesq_ssb
//		       +fn2[cb]*mesq_ccb
//		       +fn2[bb]*mesq_bbr
//		       +fn2[u]*mesq_ubu
//		       +fn2[s]*mesq_sbs
//		       +fn2[c]*mesq_cbc
//		       +fn2[b]*mesq_brb)
//	    + fn1[sb]*(fn2[ub]*mesq_uub
//		       +fn2[db]*mesq_ddb
//		       +fn2[cb]*mesq_ccb
//		       +fn2[bb]*mesq_bbr
//		       +fn2[u]*mesq_ubu
//		       +fn2[d]*mesq_dbd
//		       +fn2[c]*mesq_cbc
//		       +fn2[b]*mesq_brb)
//	    + fn1[cb]*(fn2[ub]*mesq_uub
//		       +fn2[db]*mesq_ddb
//		       +fn2[sb]*mesq_ssb
//		       +fn2[bb]*mesq_bbr
//		       +fn2[u]*mesq_ubu
//		       +fn2[d]*mesq_dbd
//		       +fn2[s]*mesq_sbs
//		       +fn2[b]*mesq_brb)
//	    + fn1[bb]*(fn2[ub]*mesq_uub
//		       +fn2[db]*mesq_ddb
//		       +fn2[sb]*mesq_ssb
//		       +fn2[cb]*mesq_ccb
//		       +fn2[u]*mesq_ubu
//		       +fn2[d]*mesq_dbd
//		       +fn2[s]*mesq_sbs
//		       +fn2[c]*mesq_cbc);
//
//	  QQPN_2 = fn2[u]*(fn1[d]*mesq_ddb
//			   +fn1[s]*mesq_ssb
//			   +fn1[c]*mesq_ccb
//			   +fn1[b]*mesq_bbr
//			   +fn1[db]*mesq_dbd
//			   +fn1[sb]*mesq_sbs
//			   +fn1[cb]*mesq_cbc
//			   +fn1[bb]*mesq_brb)
//	    + fn2[d]* (fn1[u]*mesq_uub
//		       +fn1[s]*mesq_ssb
//		       +fn1[c]*mesq_ccb
//		       +fn1[b]*mesq_bbr
//		       +fn1[ub]*mesq_ubu
//		       +fn1[sb]*mesq_sbs
//		       +fn1[cb]*mesq_cbc
//		       +fn1[bb]*mesq_brb)
//	    + fn2[s]* (fn1[u]*mesq_uub
//		       +fn1[d]*mesq_ddb
//		       +fn1[c]*mesq_ccb
//		       +fn1[b]*mesq_bbr
//		       +fn1[ub]*mesq_ubu
//		       +fn1[db]*mesq_dbd
//		       +fn1[cb]*mesq_cbc
//		       +fn1[bb]*mesq_brb)
//	    + fn2[c]* (fn1[u]*mesq_uub
//		       +fn1[d]*mesq_ddb
//		       +fn1[s]*mesq_ssb
//		       +fn1[b]*mesq_bbr
//		       +fn1[ub]*mesq_ubu
//		       +fn1[db]*mesq_dbd
//		       +fn1[sb]*mesq_sbs
//		       +fn1[bb]*mesq_brb)
//	    + fn2[b]* (fn1[u]*mesq_uub
//		       +fn1[d]*mesq_ddb
//		       +fn1[s]*mesq_ssb
//		       +fn1[c]*mesq_ccb
//		       +fn1[ub]*mesq_ubu
//		       +fn1[db]*mesq_dbd
//		       +fn1[sb]*mesq_sbs
//		       +fn1[cb]*mesq_cbc)
//	    + fn2[ub]*(fn1[d]*mesq_ddb
//		       +fn1[s]*mesq_ssb
//		       +fn1[c]*mesq_ccb
//		       +fn1[b]*mesq_bbr
//		       +fn1[db]*mesq_dbd
//		       +fn1[sb]*mesq_sbs
//		       +fn1[cb]*mesq_cbc
//		       +fn1[bb]*mesq_brb)
//	    + fn2[db]*(fn1[u]*mesq_uub
//		       +fn1[s]*mesq_ssb
//		       +fn1[c]*mesq_ccb
//		       +fn1[b]*mesq_bbr
//		       +fn1[ub]*mesq_ubu
//		       +fn1[sb]*mesq_sbs
//		       +fn1[cb]*mesq_cbc
//		       +fn1[bb]*mesq_brb)
//	    + fn2[sb]*(fn1[u]*mesq_uub
//		       +fn1[d]*mesq_ddb
//		       +fn1[c]*mesq_ccb
//		       +fn1[b]*mesq_bbr
//		       +fn1[ub]*mesq_ubu
//		       +fn1[db]*mesq_dbd
//		       +fn1[cb]*mesq_cbc
//		       +fn1[bb]*mesq_brb)
//	    + fn2[cb]*(fn1[u]*mesq_uub
//		       +fn1[d]*mesq_ddb
//		       +fn1[s]*mesq_ssb
//		       +fn1[b]*mesq_bbr
//		       +fn1[ub]*mesq_ubu
//		       +fn1[db]*mesq_dbd
//		       +fn1[sb]*mesq_sbs
//		       +fn1[bb]*mesq_brb)
//	    + fn2[bb]*(fn1[u]*mesq_uub
//		       +fn1[d]*mesq_ddb
//		       +fn1[s]*mesq_ssb
//		       +fn1[c]*mesq_ccb
//		       +fn1[ub]*mesq_ubu
//		       +fn1[db]*mesq_dbd
//		       +fn1[sb]*mesq_sbs
//		       +fn1[cb]*mesq_cbc);
	}
      if (opts.order >= 3)
	{
	  complex <double> mesq_uub_nfz = mesq::mesqij_expy_NFZ[mesq::index(0,i1,i2,sign)];
	  complex <double> mesq_ubu_nfz = mesq::mesqij_expy_NFZ[mesq::index(1,i1,i2,sign)];
	  complex <double> mesq_ddb_nfz = mesq::mesqij_expy_NFZ[mesq::index(2,i1,i2,sign)];
	  complex <double> mesq_dbd_nfz = mesq::mesqij_expy_NFZ[mesq::index(3,i1,i2,sign)];
	  complex <double> mesq_ssb_nfz = mesq::mesqij_expy_NFZ[mesq::index(4,i1,i2,sign)];
	  complex <double> mesq_sbs_nfz = mesq::mesqij_expy_NFZ[mesq::index(5,i1,i2,sign)];
	  complex <double> mesq_ccb_nfz = mesq::mesqij_expy_NFZ[mesq::index(6,i1,i2,sign)];
	  complex <double> mesq_cbc_nfz = mesq::mesqij_expy_NFZ[mesq::index(7,i1,i2,sign)];
	  complex <double> mesq_bbr_nfz = mesq::mesqij_expy_NFZ[mesq::index(8,i1,i2,sign)];
	  complex <double> mesq_brb_nfz = mesq::mesqij_expy_NFZ[mesq::index(9,i1,i2,sign)];
      
	  //      cout << "pdf_mesq_expy " << i1 << "  " << i2 << "  " << mesq_uub << "  " << fn1[u] << "  " << fn2[ub]  << endl;
	  //LL part
	  QQBN_nfz = fn1[u ]*fn2[ub]*mesq_uub_nfz
	            +fn1[c ]*fn2[cb]*mesq_ccb_nfz
	            +fn1[d ]*fn2[db]*mesq_ddb_nfz
	            +fn1[s ]*fn2[sb]*mesq_ssb_nfz
	            +fn1[b ]*fn2[bb]*mesq_bbr_nfz
	            +fn1[ub]*fn2[u ]*mesq_ubu_nfz
	            +fn1[cb]*fn2[c ]*mesq_cbc_nfz
	            +fn1[db]*fn2[d ]*mesq_dbd_nfz
	            +fn1[sb]*fn2[s ]*mesq_sbs_nfz
	            +fn1[bb]*fn2[b ]*mesq_brb_nfz;

	  QBGN_1 = fn1[g]*(fn2[u]*mesq_uub
			  +fn2[c]*mesq_ccb
			  +fn2[d]*mesq_ddb
			  +fn2[s]*mesq_ssb
			  +fn2[b]*mesq_bbr
			  +fn2[ub]*mesq_ubu
			  +fn2[cb]*mesq_cbc
			  +fn2[db]*mesq_dbd
			  +fn2[sb]*mesq_sbs
			  +fn2[bb]*mesq_brb);

	  QBGN_2 = fn2[g]*(fn1[ub]*mesq_uub
			  +fn1[cb]*mesq_ccb
			  +fn1[db]*mesq_ddb
			  +fn1[sb]*mesq_ssb
			  +fn1[bb]*mesq_bbr
			  +fn1[u]*mesq_ubu
			  +fn1[c]*mesq_cbc
			  +fn1[d]*mesq_dbd
			  +fn1[s]*mesq_sbs
			  +fn1[b]*mesq_brb);

	  QPGN_1 =
	    +fn1[g]*fn2[u]*mesq_brb
	    +fn1[g]*fn2[d]*mesq_brb
	    +fn1[g]*fn2[s]*mesq_brb
	    +fn1[g]*fn2[c]*mesq_brb
	    +fn1[g]*fn2[u]*mesq_cbc
	    +fn1[g]*fn2[d]*mesq_cbc
	    +fn1[g]*fn2[s]*mesq_cbc
	    +fn1[g]*fn2[b]*mesq_cbc
	    +fn1[g]*fn2[u]*mesq_sbs
	    +fn1[g]*fn2[d]*mesq_sbs
	    +fn1[g]*fn2[c]*mesq_sbs
	    +fn1[g]*fn2[b]*mesq_sbs
	    +fn1[g]*fn2[u]*mesq_dbd
	    +fn1[g]*fn2[s]*mesq_dbd
	    +fn1[g]*fn2[c]*mesq_dbd
	    +fn1[g]*fn2[b]*mesq_dbd
	    +fn1[g]*fn2[d]*mesq_ubu
	    +fn1[g]*fn2[s]*mesq_ubu
	    +fn1[g]*fn2[c]*mesq_ubu
	    +fn1[g]*fn2[b]*mesq_ubu
	    +fn1[g]*fn2[bb]*mesq_uub
	    +fn1[g]*fn2[cb]*mesq_uub
	    +fn1[g]*fn2[sb]*mesq_uub
	    +fn1[g]*fn2[db]*mesq_uub
	    +fn1[g]*fn2[bb]*mesq_ddb
	    +fn1[g]*fn2[cb]*mesq_ddb
	    +fn1[g]*fn2[sb]*mesq_ddb
	    +fn1[g]*fn2[ub]*mesq_ddb
	    +fn1[g]*fn2[bb]*mesq_ssb
	    +fn1[g]*fn2[cb]*mesq_ssb
	    +fn1[g]*fn2[db]*mesq_ssb
	    +fn1[g]*fn2[ub]*mesq_ssb
	    +fn1[g]*fn2[bb]*mesq_ccb
	    +fn1[g]*fn2[sb]*mesq_ccb
	    +fn1[g]*fn2[db]*mesq_ccb
	    +fn1[g]*fn2[ub]*mesq_ccb
	    +fn1[g]*fn2[cb]*mesq_bbr
	    +fn1[g]*fn2[sb]*mesq_bbr
	    +fn1[g]*fn2[db]*mesq_bbr
	    +fn1[g]*fn2[ub]*mesq_bbr;

	  QPGN_2 =
	    +fn1[u]*fn2[g]*mesq_brb
	    +fn1[d]*fn2[g]*mesq_brb
	    +fn1[s]*fn2[g]*mesq_brb
	    +fn1[c]*fn2[g]*mesq_brb
	    +fn1[u]*fn2[g]*mesq_cbc
	    +fn1[d]*fn2[g]*mesq_cbc
	    +fn1[s]*fn2[g]*mesq_cbc
	    +fn1[b]*fn2[g]*mesq_cbc
	    +fn1[u]*fn2[g]*mesq_sbs
	    +fn1[d]*fn2[g]*mesq_sbs
	    +fn1[c]*fn2[g]*mesq_sbs
	    +fn1[b]*fn2[g]*mesq_sbs
	    +fn1[u]*fn2[g]*mesq_dbd
	    +fn1[s]*fn2[g]*mesq_dbd
	    +fn1[c]*fn2[g]*mesq_dbd
	    +fn1[b]*fn2[g]*mesq_dbd
	    +fn1[d]*fn2[g]*mesq_ubu
	    +fn1[s]*fn2[g]*mesq_ubu
	    +fn1[c]*fn2[g]*mesq_ubu
	    +fn1[b]*fn2[g]*mesq_ubu
	    +fn1[bb]*fn2[g]*mesq_uub
	    +fn1[cb]*fn2[g]*mesq_uub
	    +fn1[sb]*fn2[g]*mesq_uub
	    +fn1[db]*fn2[g]*mesq_uub
	    +fn1[bb]*fn2[g]*mesq_ddb
	    +fn1[cb]*fn2[g]*mesq_ddb
	    +fn1[sb]*fn2[g]*mesq_ddb
	    +fn1[ub]*fn2[g]*mesq_ddb
	    +fn1[bb]*fn2[g]*mesq_ssb
	    +fn1[cb]*fn2[g]*mesq_ssb
	    +fn1[db]*fn2[g]*mesq_ssb
	    +fn1[ub]*fn2[g]*mesq_ssb
	    +fn1[bb]*fn2[g]*mesq_ccb
	    +fn1[sb]*fn2[g]*mesq_ccb
	    +fn1[db]*fn2[g]*mesq_ccb
	    +fn1[ub]*fn2[g]*mesq_ccb
	    +fn1[cb]*fn2[g]*mesq_bbr
	    +fn1[sb]*fn2[g]*mesq_bbr
	    +fn1[db]*fn2[g]*mesq_bbr
	    +fn1[ub]*fn2[g]*mesq_bbr;

	  QBPGN_1=
	    +fn1[g]*fn2[cb]*mesq_brb
	    +fn1[g]*fn2[sb]*mesq_brb
	    +fn1[g]*fn2[db]*mesq_brb
	    +fn1[g]*fn2[ub]*mesq_brb
	    +fn1[g]*fn2[bb]*mesq_cbc
	    +fn1[g]*fn2[sb]*mesq_cbc
	    +fn1[g]*fn2[db]*mesq_cbc
	    +fn1[g]*fn2[ub]*mesq_cbc
	    +fn1[g]*fn2[bb]*mesq_sbs
	    +fn1[g]*fn2[cb]*mesq_sbs
	    +fn1[g]*fn2[db]*mesq_sbs
	    +fn1[g]*fn2[ub]*mesq_sbs
	    +fn1[g]*fn2[bb]*mesq_dbd
	    +fn1[g]*fn2[cb]*mesq_dbd
	    +fn1[g]*fn2[sb]*mesq_dbd
	    +fn1[g]*fn2[ub]*mesq_dbd
	    +fn1[g]*fn2[bb]*mesq_ubu
	    +fn1[g]*fn2[cb]*mesq_ubu
	    +fn1[g]*fn2[sb]*mesq_ubu
	    +fn1[g]*fn2[db]*mesq_ubu
	    +fn1[g]*fn2[d]*mesq_uub
	    +fn1[g]*fn2[s]*mesq_uub
	    +fn1[g]*fn2[c]*mesq_uub
	    +fn1[g]*fn2[b]*mesq_uub
	    +fn1[g]*fn2[u]*mesq_ddb
	    +fn1[g]*fn2[s]*mesq_ddb
	    +fn1[g]*fn2[c]*mesq_ddb
	    +fn1[g]*fn2[b]*mesq_ddb
	    +fn1[g]*fn2[u]*mesq_ssb
	    +fn1[g]*fn2[d]*mesq_ssb
	    +fn1[g]*fn2[c]*mesq_ssb
	    +fn1[g]*fn2[b]*mesq_ssb
	    +fn1[g]*fn2[u]*mesq_ccb
	    +fn1[g]*fn2[d]*mesq_ccb
	    +fn1[g]*fn2[s]*mesq_ccb
	    +fn1[g]*fn2[b]*mesq_ccb
	    +fn1[g]*fn2[u]*mesq_bbr
	    +fn1[g]*fn2[d]*mesq_bbr
	    +fn1[g]*fn2[s]*mesq_bbr
	    +fn1[g]*fn2[c]*mesq_bbr;
	  

	  QBPGN_2=
	    +fn1[cb]*fn2[g]*mesq_brb
	    +fn1[sb]*fn2[g]*mesq_brb
	    +fn1[db]*fn2[g]*mesq_brb
	    +fn1[ub]*fn2[g]*mesq_brb
	    +fn1[bb]*fn2[g]*mesq_cbc
	    +fn1[sb]*fn2[g]*mesq_cbc
	    +fn1[db]*fn2[g]*mesq_cbc
	    +fn1[ub]*fn2[g]*mesq_cbc
	    +fn1[bb]*fn2[g]*mesq_sbs
	    +fn1[cb]*fn2[g]*mesq_sbs
	    +fn1[db]*fn2[g]*mesq_sbs
	    +fn1[ub]*fn2[g]*mesq_sbs
	    +fn1[bb]*fn2[g]*mesq_dbd
	    +fn1[cb]*fn2[g]*mesq_dbd
	    +fn1[sb]*fn2[g]*mesq_dbd
	    +fn1[ub]*fn2[g]*mesq_dbd
	    +fn1[bb]*fn2[g]*mesq_ubu
	    +fn1[cb]*fn2[g]*mesq_ubu
	    +fn1[sb]*fn2[g]*mesq_ubu
	    +fn1[db]*fn2[g]*mesq_ubu
	    +fn1[d]*fn2[g]*mesq_uub
	    +fn1[s]*fn2[g]*mesq_uub
	    +fn1[c]*fn2[g]*mesq_uub
	    +fn1[b]*fn2[g]*mesq_uub
	    +fn1[u]*fn2[g]*mesq_ddb
	    +fn1[s]*fn2[g]*mesq_ddb
	    +fn1[c]*fn2[g]*mesq_ddb
	    +fn1[b]*fn2[g]*mesq_ddb
	    +fn1[u]*fn2[g]*mesq_ssb
	    +fn1[d]*fn2[g]*mesq_ssb
	    +fn1[c]*fn2[g]*mesq_ssb
	    +fn1[b]*fn2[g]*mesq_ssb
	    +fn1[u]*fn2[g]*mesq_ccb
	    +fn1[d]*fn2[g]*mesq_ccb
	    +fn1[s]*fn2[g]*mesq_ccb
	    +fn1[b]*fn2[g]*mesq_ccb
	    +fn1[u]*fn2[g]*mesq_bbr
	    +fn1[d]*fn2[g]*mesq_bbr
	    +fn1[s]*fn2[g]*mesq_bbr
	    +fn1[c]*fn2[g]*mesq_bbr;
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
      if (opts.order >= 1)
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
      if (opts.order >= 2)
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

	  /*	  
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
	  */

	  QQPN_1 =
	    +fn1[u]*fn2[u]*mesq_bbu
	    +fn1[d]*fn2[u]*mesq_bbu
	    +fn1[s]*fn2[u]*mesq_bbu
	    +fn1[c]*fn2[u]*mesq_bbu
	    +fn1[u]*fn2[c]*mesq_bbc
	    +fn1[d]*fn2[c]*mesq_bbc
	    +fn1[s]*fn2[c]*mesq_bbc
	    +fn1[c]*fn2[c]*mesq_bbc
	    +fn1[u]*fn2[u]*mesq_sbu
	    +fn1[d]*fn2[u]*mesq_sbu
	    +fn1[c]*fn2[u]*mesq_sbu
	    +fn1[b]*fn2[u]*mesq_sbu
	    +fn1[u]*fn2[c]*mesq_sbc
	    +fn1[d]*fn2[c]*mesq_sbc
	    +fn1[c]*fn2[c]*mesq_sbc
	    +fn1[b]*fn2[c]*mesq_sbc
	    +fn1[u]*fn2[u]*mesq_dbu
	    +fn1[s]*fn2[u]*mesq_dbu
	    +fn1[c]*fn2[u]*mesq_dbu
	    +fn1[b]*fn2[u]*mesq_dbu
	    +fn1[u]*fn2[c]*mesq_dbc
	    +fn1[s]*fn2[c]*mesq_dbc
	    +fn1[c]*fn2[c]*mesq_dbc
	    +fn1[b]*fn2[c]*mesq_dbc
	    +fn1[bb]*fn2[bb]*mesq_ubb
	    +fn1[cb]*fn2[bb]*mesq_ubb
	    +fn1[sb]*fn2[bb]*mesq_ubb
	    +fn1[db]*fn2[bb]*mesq_ubb
	    +fn1[bb]*fn2[sb]*mesq_usb
	    +fn1[cb]*fn2[sb]*mesq_usb
	    +fn1[sb]*fn2[sb]*mesq_usb
	    +fn1[db]*fn2[sb]*mesq_usb
	    +fn1[bb]*fn2[db]*mesq_udb
	    +fn1[cb]*fn2[db]*mesq_udb
	    +fn1[sb]*fn2[db]*mesq_udb
	    +fn1[db]*fn2[db]*mesq_udb
	    +fn1[bb]*fn2[bb]*mesq_cbb
	    +fn1[sb]*fn2[bb]*mesq_cbb
	    +fn1[db]*fn2[bb]*mesq_cbb
	    +fn1[ub]*fn2[bb]*mesq_cbb
	    +fn1[bb]*fn2[sb]*mesq_csb
	    +fn1[sb]*fn2[sb]*mesq_csb
	    +fn1[db]*fn2[sb]*mesq_csb
	    +fn1[ub]*fn2[sb]*mesq_csb
	    +fn1[bb]*fn2[db]*mesq_cdb
	    +fn1[sb]*fn2[db]*mesq_cdb
	    +fn1[db]*fn2[db]*mesq_cdb
	    +fn1[ub]*fn2[db]*mesq_cdb;
	  
	  QQPN_2 =
	    +fn1[bb]*fn2[u]*mesq_bbu
	    +fn1[bb]*fn2[d]*mesq_bbu
	    +fn1[bb]*fn2[s]*mesq_bbu
	    +fn1[bb]*fn2[c]*mesq_bbu
	    +fn1[bb]*fn2[u]*mesq_bbc
	    +fn1[bb]*fn2[d]*mesq_bbc
	    +fn1[bb]*fn2[s]*mesq_bbc
	    +fn1[bb]*fn2[c]*mesq_bbc
	    +fn1[sb]*fn2[u]*mesq_sbu
	    +fn1[sb]*fn2[d]*mesq_sbu
	    +fn1[sb]*fn2[c]*mesq_sbu
	    +fn1[sb]*fn2[b]*mesq_sbu
	    +fn1[sb]*fn2[u]*mesq_sbc
	    +fn1[sb]*fn2[d]*mesq_sbc
	    +fn1[sb]*fn2[c]*mesq_sbc
	    +fn1[sb]*fn2[b]*mesq_sbc
	    +fn1[db]*fn2[u]*mesq_dbu
	    +fn1[db]*fn2[s]*mesq_dbu
	    +fn1[db]*fn2[c]*mesq_dbu
	    +fn1[db]*fn2[b]*mesq_dbu
	    +fn1[db]*fn2[u]*mesq_dbc
	    +fn1[db]*fn2[s]*mesq_dbc
	    +fn1[db]*fn2[c]*mesq_dbc
	    +fn1[db]*fn2[b]*mesq_dbc
	    +fn1[u]*fn2[bb]*mesq_ubb
	    +fn1[u]*fn2[cb]*mesq_ubb
	    +fn1[u]*fn2[sb]*mesq_ubb
	    +fn1[u]*fn2[db]*mesq_ubb
	    +fn1[u]*fn2[bb]*mesq_usb
	    +fn1[u]*fn2[cb]*mesq_usb
	    +fn1[u]*fn2[sb]*mesq_usb
	    +fn1[u]*fn2[db]*mesq_usb
	    +fn1[u]*fn2[bb]*mesq_udb
	    +fn1[u]*fn2[cb]*mesq_udb
	    +fn1[u]*fn2[sb]*mesq_udb
	    +fn1[u]*fn2[db]*mesq_udb
	    +fn1[c]*fn2[bb]*mesq_cbb
	    +fn1[c]*fn2[sb]*mesq_cbb
	    +fn1[c]*fn2[db]*mesq_cbb
	    +fn1[c]*fn2[ub]*mesq_cbb
	    +fn1[c]*fn2[bb]*mesq_csb
	    +fn1[c]*fn2[sb]*mesq_csb
	    +fn1[c]*fn2[db]*mesq_csb
	    +fn1[c]*fn2[ub]*mesq_csb
	    +fn1[c]*fn2[bb]*mesq_cdb
	    +fn1[c]*fn2[sb]*mesq_cdb
	    +fn1[c]*fn2[db]*mesq_cdb
	    +fn1[c]*fn2[ub]*mesq_cdb;

	  QQBPN_1 =
	    +fn1[cb]*fn2[u]*mesq_bbu
	    +fn1[sb]*fn2[u]*mesq_bbu
	    +fn1[db]*fn2[u]*mesq_bbu
	    +fn1[ub]*fn2[u]*mesq_bbu
	    +fn1[cb]*fn2[c]*mesq_bbc
	    +fn1[sb]*fn2[c]*mesq_bbc
	    +fn1[db]*fn2[c]*mesq_bbc
	    +fn1[ub]*fn2[c]*mesq_bbc
	    +fn1[bb]*fn2[u]*mesq_sbu
	    +fn1[cb]*fn2[u]*mesq_sbu
	    +fn1[db]*fn2[u]*mesq_sbu
	    +fn1[ub]*fn2[u]*mesq_sbu
	    +fn1[bb]*fn2[c]*mesq_sbc
	    +fn1[cb]*fn2[c]*mesq_sbc
	    +fn1[db]*fn2[c]*mesq_sbc
	    +fn1[ub]*fn2[c]*mesq_sbc
	    +fn1[bb]*fn2[u]*mesq_dbu
	    +fn1[cb]*fn2[u]*mesq_dbu
	    +fn1[sb]*fn2[u]*mesq_dbu
	    +fn1[ub]*fn2[u]*mesq_dbu
	    +fn1[bb]*fn2[c]*mesq_dbc
	    +fn1[cb]*fn2[c]*mesq_dbc
	    +fn1[sb]*fn2[c]*mesq_dbc
	    +fn1[ub]*fn2[c]*mesq_dbc
	    +fn1[d]*fn2[bb]*mesq_ubb
	    +fn1[s]*fn2[bb]*mesq_ubb
	    +fn1[c]*fn2[bb]*mesq_ubb
	    +fn1[b]*fn2[bb]*mesq_ubb
	    +fn1[d]*fn2[sb]*mesq_usb
	    +fn1[s]*fn2[sb]*mesq_usb
	    +fn1[c]*fn2[sb]*mesq_usb
	    +fn1[b]*fn2[sb]*mesq_usb
	    +fn1[d]*fn2[db]*mesq_udb
	    +fn1[s]*fn2[db]*mesq_udb
	    +fn1[c]*fn2[db]*mesq_udb
	    +fn1[b]*fn2[db]*mesq_udb
	    +fn1[u]*fn2[bb]*mesq_cbb
	    +fn1[d]*fn2[bb]*mesq_cbb
	    +fn1[s]*fn2[bb]*mesq_cbb
	    +fn1[b]*fn2[bb]*mesq_cbb
	    +fn1[u]*fn2[sb]*mesq_csb
	    +fn1[d]*fn2[sb]*mesq_csb
	    +fn1[s]*fn2[sb]*mesq_csb
	    +fn1[b]*fn2[sb]*mesq_csb
	    +fn1[u]*fn2[db]*mesq_cdb
	    +fn1[d]*fn2[db]*mesq_cdb
	    +fn1[s]*fn2[db]*mesq_cdb
	    +fn1[b]*fn2[db]*mesq_cdb;

	  QQBPN_2 =
	    +fn1[bb]*fn2[cb]*mesq_bbu
	    +fn1[bb]*fn2[sb]*mesq_bbu
	    +fn1[bb]*fn2[db]*mesq_bbu
	    +fn1[bb]*fn2[ub]*mesq_bbu
	    +fn1[bb]*fn2[cb]*mesq_bbc
	    +fn1[bb]*fn2[sb]*mesq_bbc
	    +fn1[bb]*fn2[db]*mesq_bbc
	    +fn1[bb]*fn2[ub]*mesq_bbc
	    +fn1[sb]*fn2[bb]*mesq_sbu
	    +fn1[sb]*fn2[cb]*mesq_sbu
	    +fn1[sb]*fn2[db]*mesq_sbu
	    +fn1[sb]*fn2[ub]*mesq_sbu
	    +fn1[sb]*fn2[bb]*mesq_sbc
	    +fn1[sb]*fn2[cb]*mesq_sbc
	    +fn1[sb]*fn2[db]*mesq_sbc
	    +fn1[sb]*fn2[ub]*mesq_sbc
	    +fn1[db]*fn2[bb]*mesq_dbu
	    +fn1[db]*fn2[cb]*mesq_dbu
	    +fn1[db]*fn2[sb]*mesq_dbu
	    +fn1[db]*fn2[ub]*mesq_dbu
	    +fn1[db]*fn2[bb]*mesq_dbc
	    +fn1[db]*fn2[cb]*mesq_dbc
	    +fn1[db]*fn2[sb]*mesq_dbc
	    +fn1[db]*fn2[ub]*mesq_dbc
	    +fn1[u]*fn2[d]*mesq_ubb
	    +fn1[u]*fn2[s]*mesq_ubb
	    +fn1[u]*fn2[c]*mesq_ubb
	    +fn1[u]*fn2[b]*mesq_ubb
	    +fn1[u]*fn2[d]*mesq_usb
	    +fn1[u]*fn2[s]*mesq_usb
	    +fn1[u]*fn2[c]*mesq_usb
	    +fn1[u]*fn2[b]*mesq_usb
	    +fn1[u]*fn2[d]*mesq_udb
	    +fn1[u]*fn2[s]*mesq_udb
	    +fn1[u]*fn2[c]*mesq_udb
	    +fn1[u]*fn2[b]*mesq_udb
	    +fn1[c]*fn2[u]*mesq_cbb
	    +fn1[c]*fn2[d]*mesq_cbb
	    +fn1[c]*fn2[s]*mesq_cbb
	    +fn1[c]*fn2[b]*mesq_cbb
	    +fn1[c]*fn2[u]*mesq_csb
	    +fn1[c]*fn2[d]*mesq_csb
	    +fn1[c]*fn2[s]*mesq_csb
	    +fn1[c]*fn2[b]*mesq_csb
	    +fn1[c]*fn2[u]*mesq_cdb
	    +fn1[c]*fn2[d]*mesq_cdb
	    +fn1[c]*fn2[s]*mesq_cdb
	    +fn1[c]*fn2[b]*mesq_cdb;
	}
      if (opts.order >= 3)
	{
	  QBGN_1 =
	    +fn1[g]*fn2[ub]*mesq_bbu
	    +fn1[g]*fn2[cb]*mesq_bbc
	    +fn1[g]*fn2[ub]*mesq_sbu
	    +fn1[g]*fn2[cb]*mesq_sbc
	    +fn1[g]*fn2[ub]*mesq_dbu
	    +fn1[g]*fn2[cb]*mesq_dbc
	    +fn1[g]*fn2[b]*mesq_ubb
	    +fn1[g]*fn2[s]*mesq_usb
	    +fn1[g]*fn2[d]*mesq_udb
	    +fn1[g]*fn2[b]*mesq_cbb
	    +fn1[g]*fn2[s]*mesq_csb
	    +fn1[g]*fn2[d]*mesq_cdb;
	  
	  QBGN_2 =
	    +fn1[b]*fn2[g]*mesq_bbu
	    +fn1[b]*fn2[g]*mesq_bbc
	    +fn1[s]*fn2[g]*mesq_sbu
	    +fn1[s]*fn2[g]*mesq_sbc
	    +fn1[d]*fn2[g]*mesq_dbu
	    +fn1[d]*fn2[g]*mesq_dbc
	    +fn1[ub]*fn2[g]*mesq_ubb
	    +fn1[ub]*fn2[g]*mesq_usb
	    +fn1[ub]*fn2[g]*mesq_udb
	    +fn1[cb]*fn2[g]*mesq_cbb
	    +fn1[cb]*fn2[g]*mesq_csb
	    +fn1[cb]*fn2[g]*mesq_cdb;
	  
	  QPGN_1 =
	    +fn1[g]*fn2[u]*mesq_bbu
	    +fn1[g]*fn2[d]*mesq_bbu
	    +fn1[g]*fn2[s]*mesq_bbu
	    +fn1[g]*fn2[c]*mesq_bbu
	    +fn1[g]*fn2[u]*mesq_bbc
	    +fn1[g]*fn2[d]*mesq_bbc
	    +fn1[g]*fn2[s]*mesq_bbc
	    +fn1[g]*fn2[c]*mesq_bbc
	    +fn1[g]*fn2[u]*mesq_sbu
	    +fn1[g]*fn2[d]*mesq_sbu
	    +fn1[g]*fn2[c]*mesq_sbu
	    +fn1[g]*fn2[b]*mesq_sbu
	    +fn1[g]*fn2[u]*mesq_sbc
	    +fn1[g]*fn2[d]*mesq_sbc
	    +fn1[g]*fn2[c]*mesq_sbc
	    +fn1[g]*fn2[b]*mesq_sbc
	    +fn1[g]*fn2[u]*mesq_dbu
	    +fn1[g]*fn2[s]*mesq_dbu
	    +fn1[g]*fn2[c]*mesq_dbu
	    +fn1[g]*fn2[b]*mesq_dbu
	    +fn1[g]*fn2[u]*mesq_dbc
	    +fn1[g]*fn2[s]*mesq_dbc
	    +fn1[g]*fn2[c]*mesq_dbc
	    +fn1[g]*fn2[b]*mesq_dbc
	    +fn1[g]*fn2[bb]*mesq_ubb
	    +fn1[g]*fn2[cb]*mesq_ubb
	    +fn1[g]*fn2[sb]*mesq_ubb
	    +fn1[g]*fn2[db]*mesq_ubb
	    +fn1[g]*fn2[bb]*mesq_usb
	    +fn1[g]*fn2[cb]*mesq_usb
	    +fn1[g]*fn2[sb]*mesq_usb
	    +fn1[g]*fn2[db]*mesq_usb
	    +fn1[g]*fn2[bb]*mesq_udb
	    +fn1[g]*fn2[cb]*mesq_udb
	    +fn1[g]*fn2[sb]*mesq_udb
	    +fn1[g]*fn2[db]*mesq_udb
	    +fn1[g]*fn2[bb]*mesq_cbb
	    +fn1[g]*fn2[sb]*mesq_cbb
	    +fn1[g]*fn2[db]*mesq_cbb
	    +fn1[g]*fn2[ub]*mesq_cbb
	    +fn1[g]*fn2[bb]*mesq_csb
	    +fn1[g]*fn2[sb]*mesq_csb
	    +fn1[g]*fn2[db]*mesq_csb
	    +fn1[g]*fn2[ub]*mesq_csb
	    +fn1[g]*fn2[bb]*mesq_cdb
	    +fn1[g]*fn2[sb]*mesq_cdb
	    +fn1[g]*fn2[db]*mesq_cdb
	    +fn1[g]*fn2[ub]*mesq_cdb;

	  QPGN_2 =
	    +fn1[u]*fn2[g]*mesq_bbu
	    +fn1[d]*fn2[g]*mesq_bbu
	    +fn1[s]*fn2[g]*mesq_bbu
	    +fn1[c]*fn2[g]*mesq_bbu
	    +fn1[u]*fn2[g]*mesq_bbc
	    +fn1[d]*fn2[g]*mesq_bbc
	    +fn1[s]*fn2[g]*mesq_bbc
	    +fn1[c]*fn2[g]*mesq_bbc
	    +fn1[u]*fn2[g]*mesq_sbu
	    +fn1[d]*fn2[g]*mesq_sbu
	    +fn1[c]*fn2[g]*mesq_sbu
	    +fn1[b]*fn2[g]*mesq_sbu
	    +fn1[u]*fn2[g]*mesq_sbc
	    +fn1[d]*fn2[g]*mesq_sbc
	    +fn1[c]*fn2[g]*mesq_sbc
	    +fn1[b]*fn2[g]*mesq_sbc
	    +fn1[u]*fn2[g]*mesq_dbu
	    +fn1[s]*fn2[g]*mesq_dbu
	    +fn1[c]*fn2[g]*mesq_dbu
	    +fn1[b]*fn2[g]*mesq_dbu
	    +fn1[u]*fn2[g]*mesq_dbc
	    +fn1[s]*fn2[g]*mesq_dbc
	    +fn1[c]*fn2[g]*mesq_dbc
	    +fn1[b]*fn2[g]*mesq_dbc
	    +fn1[bb]*fn2[g]*mesq_ubb
	    +fn1[cb]*fn2[g]*mesq_ubb
	    +fn1[sb]*fn2[g]*mesq_ubb
	    +fn1[db]*fn2[g]*mesq_ubb
	    +fn1[bb]*fn2[g]*mesq_usb
	    +fn1[cb]*fn2[g]*mesq_usb
	    +fn1[sb]*fn2[g]*mesq_usb
	    +fn1[db]*fn2[g]*mesq_usb
	    +fn1[bb]*fn2[g]*mesq_udb
	    +fn1[cb]*fn2[g]*mesq_udb
	    +fn1[sb]*fn2[g]*mesq_udb
	    +fn1[db]*fn2[g]*mesq_udb
	    +fn1[bb]*fn2[g]*mesq_cbb
	    +fn1[sb]*fn2[g]*mesq_cbb
	    +fn1[db]*fn2[g]*mesq_cbb
	    +fn1[ub]*fn2[g]*mesq_cbb
	    +fn1[bb]*fn2[g]*mesq_csb
	    +fn1[sb]*fn2[g]*mesq_csb
	    +fn1[db]*fn2[g]*mesq_csb
	    +fn1[ub]*fn2[g]*mesq_csb
	    +fn1[bb]*fn2[g]*mesq_cdb
	    +fn1[sb]*fn2[g]*mesq_cdb
	    +fn1[db]*fn2[g]*mesq_cdb
	    +fn1[ub]*fn2[g]*mesq_cdb;

	  QBPGN_1 =
	    +fn1[g]*fn2[cb]*mesq_bbu
	    +fn1[g]*fn2[sb]*mesq_bbu
	    +fn1[g]*fn2[db]*mesq_bbu
	    +fn1[g]*fn2[ub]*mesq_bbu
	    +fn1[g]*fn2[cb]*mesq_bbc
	    +fn1[g]*fn2[sb]*mesq_bbc
	    +fn1[g]*fn2[db]*mesq_bbc
	    +fn1[g]*fn2[ub]*mesq_bbc
	    +fn1[g]*fn2[bb]*mesq_sbu
	    +fn1[g]*fn2[cb]*mesq_sbu
	    +fn1[g]*fn2[db]*mesq_sbu
	    +fn1[g]*fn2[ub]*mesq_sbu
	    +fn1[g]*fn2[bb]*mesq_sbc
	    +fn1[g]*fn2[cb]*mesq_sbc
	    +fn1[g]*fn2[db]*mesq_sbc
	    +fn1[g]*fn2[ub]*mesq_sbc
	    +fn1[g]*fn2[bb]*mesq_dbu
	    +fn1[g]*fn2[cb]*mesq_dbu
	    +fn1[g]*fn2[sb]*mesq_dbu
	    +fn1[g]*fn2[ub]*mesq_dbu
	    +fn1[g]*fn2[bb]*mesq_dbc
	    +fn1[g]*fn2[cb]*mesq_dbc
	    +fn1[g]*fn2[sb]*mesq_dbc
	    +fn1[g]*fn2[ub]*mesq_dbc
	    +fn1[g]*fn2[d]*mesq_ubb
	    +fn1[g]*fn2[s]*mesq_ubb
	    +fn1[g]*fn2[c]*mesq_ubb
	    +fn1[g]*fn2[b]*mesq_ubb
	    +fn1[g]*fn2[d]*mesq_usb
	    +fn1[g]*fn2[s]*mesq_usb
	    +fn1[g]*fn2[c]*mesq_usb
	    +fn1[g]*fn2[b]*mesq_usb
	    +fn1[g]*fn2[d]*mesq_udb
	    +fn1[g]*fn2[s]*mesq_udb
	    +fn1[g]*fn2[c]*mesq_udb
	    +fn1[g]*fn2[b]*mesq_udb
	    +fn1[g]*fn2[u]*mesq_cbb
	    +fn1[g]*fn2[d]*mesq_cbb
	    +fn1[g]*fn2[s]*mesq_cbb
	    +fn1[g]*fn2[b]*mesq_cbb
	    +fn1[g]*fn2[u]*mesq_csb
	    +fn1[g]*fn2[d]*mesq_csb
	    +fn1[g]*fn2[s]*mesq_csb
	    +fn1[g]*fn2[b]*mesq_csb
	    +fn1[g]*fn2[u]*mesq_cdb
	    +fn1[g]*fn2[d]*mesq_cdb
	    +fn1[g]*fn2[s]*mesq_cdb
	    +fn1[g]*fn2[b]*mesq_cdb;

	  QBPGN_2 =
	    +fn1[cb]*fn2[g]*mesq_bbu
	    +fn1[sb]*fn2[g]*mesq_bbu
	    +fn1[db]*fn2[g]*mesq_bbu
	    +fn1[ub]*fn2[g]*mesq_bbu
	    +fn1[cb]*fn2[g]*mesq_bbc
	    +fn1[sb]*fn2[g]*mesq_bbc
	    +fn1[db]*fn2[g]*mesq_bbc
	    +fn1[ub]*fn2[g]*mesq_bbc
	    +fn1[bb]*fn2[g]*mesq_sbu
	    +fn1[cb]*fn2[g]*mesq_sbu
	    +fn1[db]*fn2[g]*mesq_sbu
	    +fn1[ub]*fn2[g]*mesq_sbu
	    +fn1[bb]*fn2[g]*mesq_sbc
	    +fn1[cb]*fn2[g]*mesq_sbc
	    +fn1[db]*fn2[g]*mesq_sbc
	    +fn1[ub]*fn2[g]*mesq_sbc
	    +fn1[bb]*fn2[g]*mesq_dbu
	    +fn1[cb]*fn2[g]*mesq_dbu
	    +fn1[sb]*fn2[g]*mesq_dbu
	    +fn1[ub]*fn2[g]*mesq_dbu
	    +fn1[bb]*fn2[g]*mesq_dbc
	    +fn1[cb]*fn2[g]*mesq_dbc
	    +fn1[sb]*fn2[g]*mesq_dbc
	    +fn1[ub]*fn2[g]*mesq_dbc
	    +fn1[d]*fn2[g]*mesq_ubb
	    +fn1[s]*fn2[g]*mesq_ubb
	    +fn1[c]*fn2[g]*mesq_ubb
	    +fn1[b]*fn2[g]*mesq_ubb
	    +fn1[d]*fn2[g]*mesq_usb
	    +fn1[s]*fn2[g]*mesq_usb
	    +fn1[c]*fn2[g]*mesq_usb
	    +fn1[b]*fn2[g]*mesq_usb
	    +fn1[d]*fn2[g]*mesq_udb
	    +fn1[s]*fn2[g]*mesq_udb
	    +fn1[c]*fn2[g]*mesq_udb
	    +fn1[b]*fn2[g]*mesq_udb
	    +fn1[u]*fn2[g]*mesq_cbb
	    +fn1[d]*fn2[g]*mesq_cbb
	    +fn1[s]*fn2[g]*mesq_cbb
	    +fn1[b]*fn2[g]*mesq_cbb
	    +fn1[u]*fn2[g]*mesq_csb
	    +fn1[d]*fn2[g]*mesq_csb
	    +fn1[s]*fn2[g]*mesq_csb
	    +fn1[b]*fn2[g]*mesq_csb
	    +fn1[u]*fn2[g]*mesq_cdb
	    +fn1[d]*fn2[g]*mesq_cdb
	    +fn1[s]*fn2[g]*mesq_cdb
	    +fn1[b]*fn2[g]*mesq_cdb;
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
      if (opts.order >= 1)
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
      if (opts.order >= 2)
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
	  /*
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
	  */

	  QQPN_1 =
	    +fn1[u]*fn2[d]*mesq_cbd
	    +fn1[d]*fn2[d]*mesq_cbd
	    +fn1[s]*fn2[d]*mesq_cbd
	    +fn1[b]*fn2[d]*mesq_cbd
	    +fn1[u]*fn2[s]*mesq_cbs
	    +fn1[d]*fn2[s]*mesq_cbs
	    +fn1[s]*fn2[s]*mesq_cbs
	    +fn1[b]*fn2[s]*mesq_cbs
	    +fn1[u]*fn2[b]*mesq_cbb
	    +fn1[d]*fn2[b]*mesq_cbb
	    +fn1[s]*fn2[b]*mesq_cbb
	    +fn1[b]*fn2[b]*mesq_cbb
	    +fn1[d]*fn2[d]*mesq_ubd
	    +fn1[s]*fn2[d]*mesq_ubd
	    +fn1[c]*fn2[d]*mesq_ubd
	    +fn1[b]*fn2[d]*mesq_ubd
	    +fn1[d]*fn2[s]*mesq_ubs
	    +fn1[s]*fn2[s]*mesq_ubs
	    +fn1[c]*fn2[s]*mesq_ubs
	    +fn1[b]*fn2[s]*mesq_ubs
	    +fn1[d]*fn2[b]*mesq_ubb
	    +fn1[s]*fn2[b]*mesq_ubb
	    +fn1[c]*fn2[b]*mesq_ubb
	    +fn1[b]*fn2[b]*mesq_ubb
	    +fn1[bb]*fn2[cb]*mesq_dcb
	    +fn1[cb]*fn2[cb]*mesq_dcb
	    +fn1[sb]*fn2[cb]*mesq_dcb
	    +fn1[ub]*fn2[cb]*mesq_dcb
	    +fn1[bb]*fn2[ub]*mesq_dub
	    +fn1[cb]*fn2[ub]*mesq_dub
	    +fn1[sb]*fn2[ub]*mesq_dub
	    +fn1[ub]*fn2[ub]*mesq_dub
	    +fn1[bb]*fn2[cb]*mesq_scb
	    +fn1[cb]*fn2[cb]*mesq_scb
	    +fn1[db]*fn2[cb]*mesq_scb
	    +fn1[ub]*fn2[cb]*mesq_scb
	    +fn1[bb]*fn2[ub]*mesq_sub
	    +fn1[cb]*fn2[ub]*mesq_sub
	    +fn1[db]*fn2[ub]*mesq_sub
	    +fn1[ub]*fn2[ub]*mesq_sub
	    +fn1[cb]*fn2[cb]*mesq_bcb
	    +fn1[sb]*fn2[cb]*mesq_bcb
	    +fn1[db]*fn2[cb]*mesq_bcb
	    +fn1[ub]*fn2[cb]*mesq_bcb
	    +fn1[cb]*fn2[ub]*mesq_bub
	    +fn1[sb]*fn2[ub]*mesq_bub
	    +fn1[db]*fn2[ub]*mesq_bub
	    +fn1[ub]*fn2[ub]*mesq_bub;

	  QQPN_2 =
	    +fn1[cb]*fn2[u]*mesq_cbd
	    +fn1[cb]*fn2[d]*mesq_cbd
	    +fn1[cb]*fn2[s]*mesq_cbd
	    +fn1[cb]*fn2[b]*mesq_cbd
	    +fn1[cb]*fn2[u]*mesq_cbs
	    +fn1[cb]*fn2[d]*mesq_cbs
	    +fn1[cb]*fn2[s]*mesq_cbs
	    +fn1[cb]*fn2[b]*mesq_cbs
	    +fn1[cb]*fn2[u]*mesq_cbb
	    +fn1[cb]*fn2[d]*mesq_cbb
	    +fn1[cb]*fn2[s]*mesq_cbb
	    +fn1[cb]*fn2[b]*mesq_cbb
	    +fn1[ub]*fn2[d]*mesq_ubd
	    +fn1[ub]*fn2[s]*mesq_ubd
	    +fn1[ub]*fn2[c]*mesq_ubd
	    +fn1[ub]*fn2[b]*mesq_ubd
	    +fn1[ub]*fn2[d]*mesq_ubs
	    +fn1[ub]*fn2[s]*mesq_ubs
	    +fn1[ub]*fn2[c]*mesq_ubs
	    +fn1[ub]*fn2[b]*mesq_ubs
	    +fn1[ub]*fn2[d]*mesq_ubb
	    +fn1[ub]*fn2[s]*mesq_ubb
	    +fn1[ub]*fn2[c]*mesq_ubb
	    +fn1[ub]*fn2[b]*mesq_ubb
	    +fn1[d]*fn2[bb]*mesq_dcb
	    +fn1[d]*fn2[cb]*mesq_dcb
	    +fn1[d]*fn2[sb]*mesq_dcb
	    +fn1[d]*fn2[ub]*mesq_dcb
	    +fn1[d]*fn2[bb]*mesq_dub
	    +fn1[d]*fn2[cb]*mesq_dub
	    +fn1[d]*fn2[sb]*mesq_dub
	    +fn1[d]*fn2[ub]*mesq_dub
	    +fn1[s]*fn2[bb]*mesq_scb
	    +fn1[s]*fn2[cb]*mesq_scb
	    +fn1[s]*fn2[db]*mesq_scb
	    +fn1[s]*fn2[ub]*mesq_scb
	    +fn1[s]*fn2[bb]*mesq_sub
	    +fn1[s]*fn2[cb]*mesq_sub
	    +fn1[s]*fn2[db]*mesq_sub
	    +fn1[s]*fn2[ub]*mesq_sub
	    +fn1[b]*fn2[cb]*mesq_bcb
	    +fn1[b]*fn2[sb]*mesq_bcb
	    +fn1[b]*fn2[db]*mesq_bcb
	    +fn1[b]*fn2[ub]*mesq_bcb
	    +fn1[b]*fn2[cb]*mesq_bub
	    +fn1[b]*fn2[sb]*mesq_bub
	    +fn1[b]*fn2[db]*mesq_bub
	    +fn1[b]*fn2[ub]*mesq_bub;

	  QQBPN_1 =
	    +fn1[bb]*fn2[d]*mesq_cbd
	    +fn1[sb]*fn2[d]*mesq_cbd
	    +fn1[db]*fn2[d]*mesq_cbd
	    +fn1[ub]*fn2[d]*mesq_cbd
	    +fn1[bb]*fn2[s]*mesq_cbs
	    +fn1[sb]*fn2[s]*mesq_cbs
	    +fn1[db]*fn2[s]*mesq_cbs
	    +fn1[ub]*fn2[s]*mesq_cbs
	    +fn1[bb]*fn2[b]*mesq_cbb
	    +fn1[sb]*fn2[b]*mesq_cbb
	    +fn1[db]*fn2[b]*mesq_cbb
	    +fn1[ub]*fn2[b]*mesq_cbb
	    +fn1[bb]*fn2[d]*mesq_ubd
	    +fn1[cb]*fn2[d]*mesq_ubd
	    +fn1[sb]*fn2[d]*mesq_ubd
	    +fn1[db]*fn2[d]*mesq_ubd
	    +fn1[bb]*fn2[s]*mesq_ubs
	    +fn1[cb]*fn2[s]*mesq_ubs
	    +fn1[sb]*fn2[s]*mesq_ubs
	    +fn1[db]*fn2[s]*mesq_ubs
	    +fn1[bb]*fn2[b]*mesq_ubb
	    +fn1[cb]*fn2[b]*mesq_ubb
	    +fn1[sb]*fn2[b]*mesq_ubb
	    +fn1[db]*fn2[b]*mesq_ubb
	    +fn1[u]*fn2[cb]*mesq_dcb
	    +fn1[s]*fn2[cb]*mesq_dcb
	    +fn1[c]*fn2[cb]*mesq_dcb
	    +fn1[b]*fn2[cb]*mesq_dcb
	    +fn1[u]*fn2[ub]*mesq_dub
	    +fn1[s]*fn2[ub]*mesq_dub
	    +fn1[c]*fn2[ub]*mesq_dub
	    +fn1[b]*fn2[ub]*mesq_dub
	    +fn1[u]*fn2[cb]*mesq_scb
	    +fn1[d]*fn2[cb]*mesq_scb
	    +fn1[c]*fn2[cb]*mesq_scb
	    +fn1[b]*fn2[cb]*mesq_scb
	    +fn1[u]*fn2[ub]*mesq_sub
	    +fn1[d]*fn2[ub]*mesq_sub
	    +fn1[c]*fn2[ub]*mesq_sub
	    +fn1[b]*fn2[ub]*mesq_sub
	    +fn1[u]*fn2[cb]*mesq_bcb
	    +fn1[d]*fn2[cb]*mesq_bcb
	    +fn1[s]*fn2[cb]*mesq_bcb
	    +fn1[c]*fn2[cb]*mesq_bcb
	    +fn1[u]*fn2[ub]*mesq_bub
	    +fn1[d]*fn2[ub]*mesq_bub
	    +fn1[s]*fn2[ub]*mesq_bub
	    +fn1[c]*fn2[ub]*mesq_bub;

	  QQBPN_2 =
	    +fn1[cb]*fn2[bb]*mesq_cbd
	    +fn1[cb]*fn2[sb]*mesq_cbd
	    +fn1[cb]*fn2[db]*mesq_cbd
	    +fn1[cb]*fn2[ub]*mesq_cbd
	    +fn1[cb]*fn2[bb]*mesq_cbs
	    +fn1[cb]*fn2[sb]*mesq_cbs
	    +fn1[cb]*fn2[db]*mesq_cbs
	    +fn1[cb]*fn2[ub]*mesq_cbs
	    +fn1[cb]*fn2[bb]*mesq_cbb
	    +fn1[cb]*fn2[sb]*mesq_cbb
	    +fn1[cb]*fn2[db]*mesq_cbb
	    +fn1[cb]*fn2[ub]*mesq_cbb
	    +fn1[ub]*fn2[bb]*mesq_ubd
	    +fn1[ub]*fn2[cb]*mesq_ubd
	    +fn1[ub]*fn2[sb]*mesq_ubd
	    +fn1[ub]*fn2[db]*mesq_ubd
	    +fn1[ub]*fn2[bb]*mesq_ubs
	    +fn1[ub]*fn2[cb]*mesq_ubs
	    +fn1[ub]*fn2[sb]*mesq_ubs
	    +fn1[ub]*fn2[db]*mesq_ubs
	    +fn1[ub]*fn2[bb]*mesq_ubb
	    +fn1[ub]*fn2[cb]*mesq_ubb
	    +fn1[ub]*fn2[sb]*mesq_ubb
	    +fn1[ub]*fn2[db]*mesq_ubb
	    +fn1[d]*fn2[u]*mesq_dcb
	    +fn1[d]*fn2[s]*mesq_dcb
	    +fn1[d]*fn2[c]*mesq_dcb
	    +fn1[d]*fn2[b]*mesq_dcb
	    +fn1[d]*fn2[u]*mesq_dub
	    +fn1[d]*fn2[s]*mesq_dub
	    +fn1[d]*fn2[c]*mesq_dub
	    +fn1[d]*fn2[b]*mesq_dub
	    +fn1[s]*fn2[u]*mesq_scb
	    +fn1[s]*fn2[d]*mesq_scb
	    +fn1[s]*fn2[c]*mesq_scb
	    +fn1[s]*fn2[b]*mesq_scb
	    +fn1[s]*fn2[u]*mesq_sub
	    +fn1[s]*fn2[d]*mesq_sub
	    +fn1[s]*fn2[c]*mesq_sub
	    +fn1[s]*fn2[b]*mesq_sub
	    +fn1[b]*fn2[u]*mesq_bcb
	    +fn1[b]*fn2[d]*mesq_bcb
	    +fn1[b]*fn2[s]*mesq_bcb
	    +fn1[b]*fn2[c]*mesq_bcb
	    +fn1[b]*fn2[u]*mesq_bub
	    +fn1[b]*fn2[d]*mesq_bub
	    +fn1[b]*fn2[s]*mesq_bub
	    +fn1[b]*fn2[c]*mesq_bub;
	}
      if (opts.order >= 3)
	{
	  QBGN_1 =
	    +fn1[g]*fn2[db]*mesq_cbd
	    +fn1[g]*fn2[sb]*mesq_cbs
	    +fn1[g]*fn2[bb]*mesq_cbb
	    +fn1[g]*fn2[db]*mesq_ubd
	    +fn1[g]*fn2[sb]*mesq_ubs
	    +fn1[g]*fn2[bb]*mesq_ubb
	    +fn1[g]*fn2[c]*mesq_dcb
	    +fn1[g]*fn2[u]*mesq_dub
	    +fn1[g]*fn2[c]*mesq_scb
	    +fn1[g]*fn2[u]*mesq_sub
	    +fn1[g]*fn2[c]*mesq_bcb
	    +fn1[g]*fn2[u]*mesq_bub;

	  QBGN_2 =
	    +fn1[c]*fn2[g]*mesq_cbd
	    +fn1[c]*fn2[g]*mesq_cbs
	    +fn1[c]*fn2[g]*mesq_cbb
	    +fn1[u]*fn2[g]*mesq_ubd
	    +fn1[u]*fn2[g]*mesq_ubs
	    +fn1[u]*fn2[g]*mesq_ubb
	    +fn1[db]*fn2[g]*mesq_dcb
	    +fn1[db]*fn2[g]*mesq_dub
	    +fn1[sb]*fn2[g]*mesq_scb
	    +fn1[sb]*fn2[g]*mesq_sub
	    +fn1[bb]*fn2[g]*mesq_bcb
	    +fn1[bb]*fn2[g]*mesq_bub;

	  QPGN_1 =
	    +fn1[g]*fn2[u]*mesq_cbd
	    +fn1[g]*fn2[d]*mesq_cbd
	    +fn1[g]*fn2[s]*mesq_cbd
	    +fn1[g]*fn2[b]*mesq_cbd
	    +fn1[g]*fn2[u]*mesq_cbs
	    +fn1[g]*fn2[d]*mesq_cbs
	    +fn1[g]*fn2[s]*mesq_cbs
	    +fn1[g]*fn2[b]*mesq_cbs
	    +fn1[g]*fn2[u]*mesq_cbb
	    +fn1[g]*fn2[d]*mesq_cbb
	    +fn1[g]*fn2[s]*mesq_cbb
	    +fn1[g]*fn2[b]*mesq_cbb
	    +fn1[g]*fn2[d]*mesq_ubd
	    +fn1[g]*fn2[s]*mesq_ubd
	    +fn1[g]*fn2[c]*mesq_ubd
	    +fn1[g]*fn2[b]*mesq_ubd
	    +fn1[g]*fn2[d]*mesq_ubs
	    +fn1[g]*fn2[s]*mesq_ubs
	    +fn1[g]*fn2[c]*mesq_ubs
	    +fn1[g]*fn2[b]*mesq_ubs
	    +fn1[g]*fn2[d]*mesq_ubb
	    +fn1[g]*fn2[s]*mesq_ubb
	    +fn1[g]*fn2[c]*mesq_ubb
	    +fn1[g]*fn2[b]*mesq_ubb
	    +fn1[g]*fn2[bb]*mesq_dcb
	    +fn1[g]*fn2[cb]*mesq_dcb
	    +fn1[g]*fn2[sb]*mesq_dcb
	    +fn1[g]*fn2[ub]*mesq_dcb
	    +fn1[g]*fn2[bb]*mesq_dub
	    +fn1[g]*fn2[cb]*mesq_dub
	    +fn1[g]*fn2[sb]*mesq_dub
	    +fn1[g]*fn2[ub]*mesq_dub
	    +fn1[g]*fn2[bb]*mesq_scb
	    +fn1[g]*fn2[cb]*mesq_scb
	    +fn1[g]*fn2[db]*mesq_scb
	    +fn1[g]*fn2[ub]*mesq_scb
	    +fn1[g]*fn2[bb]*mesq_sub
	    +fn1[g]*fn2[cb]*mesq_sub
	    +fn1[g]*fn2[db]*mesq_sub
	    +fn1[g]*fn2[ub]*mesq_sub
	    +fn1[g]*fn2[cb]*mesq_bcb
	    +fn1[g]*fn2[sb]*mesq_bcb
	    +fn1[g]*fn2[db]*mesq_bcb
	    +fn1[g]*fn2[ub]*mesq_bcb
	    +fn1[g]*fn2[cb]*mesq_bub
	    +fn1[g]*fn2[sb]*mesq_bub
	    +fn1[g]*fn2[db]*mesq_bub
	    +fn1[g]*fn2[ub]*mesq_bub;

	  QPGN_2 =
	    +fn1[u]*fn2[g]*mesq_cbd
	    +fn1[d]*fn2[g]*mesq_cbd
	    +fn1[s]*fn2[g]*mesq_cbd
	    +fn1[b]*fn2[g]*mesq_cbd
	    +fn1[u]*fn2[g]*mesq_cbs
	    +fn1[d]*fn2[g]*mesq_cbs
	    +fn1[s]*fn2[g]*mesq_cbs
	    +fn1[b]*fn2[g]*mesq_cbs
	    +fn1[u]*fn2[g]*mesq_cbb
	    +fn1[d]*fn2[g]*mesq_cbb
	    +fn1[s]*fn2[g]*mesq_cbb
	    +fn1[b]*fn2[g]*mesq_cbb
	    +fn1[d]*fn2[g]*mesq_ubd
	    +fn1[s]*fn2[g]*mesq_ubd
	    +fn1[c]*fn2[g]*mesq_ubd
	    +fn1[b]*fn2[g]*mesq_ubd
	    +fn1[d]*fn2[g]*mesq_ubs
	    +fn1[s]*fn2[g]*mesq_ubs
	    +fn1[c]*fn2[g]*mesq_ubs
	    +fn1[b]*fn2[g]*mesq_ubs
	    +fn1[d]*fn2[g]*mesq_ubb
	    +fn1[s]*fn2[g]*mesq_ubb
	    +fn1[c]*fn2[g]*mesq_ubb
	    +fn1[b]*fn2[g]*mesq_ubb
	    +fn1[bb]*fn2[g]*mesq_dcb
	    +fn1[cb]*fn2[g]*mesq_dcb
	    +fn1[sb]*fn2[g]*mesq_dcb
	    +fn1[ub]*fn2[g]*mesq_dcb
	    +fn1[bb]*fn2[g]*mesq_dub
	    +fn1[cb]*fn2[g]*mesq_dub
	    +fn1[sb]*fn2[g]*mesq_dub
	    +fn1[ub]*fn2[g]*mesq_dub
	    +fn1[bb]*fn2[g]*mesq_scb
	    +fn1[cb]*fn2[g]*mesq_scb
	    +fn1[db]*fn2[g]*mesq_scb
	    +fn1[ub]*fn2[g]*mesq_scb
	    +fn1[bb]*fn2[g]*mesq_sub
	    +fn1[cb]*fn2[g]*mesq_sub
	    +fn1[db]*fn2[g]*mesq_sub
	    +fn1[ub]*fn2[g]*mesq_sub
	    +fn1[cb]*fn2[g]*mesq_bcb
	    +fn1[sb]*fn2[g]*mesq_bcb
	    +fn1[db]*fn2[g]*mesq_bcb
	    +fn1[ub]*fn2[g]*mesq_bcb
	    +fn1[cb]*fn2[g]*mesq_bub
	    +fn1[sb]*fn2[g]*mesq_bub
	    +fn1[db]*fn2[g]*mesq_bub
	    +fn1[ub]*fn2[g]*mesq_bub;

	  QBPGN_1 =
	    +fn1[g]*fn2[bb]*mesq_cbd
	    +fn1[g]*fn2[sb]*mesq_cbd
	    +fn1[g]*fn2[db]*mesq_cbd
	    +fn1[g]*fn2[ub]*mesq_cbd
	    +fn1[g]*fn2[bb]*mesq_cbs
	    +fn1[g]*fn2[sb]*mesq_cbs
	    +fn1[g]*fn2[db]*mesq_cbs
	    +fn1[g]*fn2[ub]*mesq_cbs
	    +fn1[g]*fn2[bb]*mesq_cbb
	    +fn1[g]*fn2[sb]*mesq_cbb
	    +fn1[g]*fn2[db]*mesq_cbb
	    +fn1[g]*fn2[ub]*mesq_cbb
	    +fn1[g]*fn2[bb]*mesq_ubd
	    +fn1[g]*fn2[cb]*mesq_ubd
	    +fn1[g]*fn2[sb]*mesq_ubd
	    +fn1[g]*fn2[db]*mesq_ubd
	    +fn1[g]*fn2[bb]*mesq_ubs
	    +fn1[g]*fn2[cb]*mesq_ubs
	    +fn1[g]*fn2[sb]*mesq_ubs
	    +fn1[g]*fn2[db]*mesq_ubs
	    +fn1[g]*fn2[bb]*mesq_ubb
	    +fn1[g]*fn2[cb]*mesq_ubb
	    +fn1[g]*fn2[sb]*mesq_ubb
	    +fn1[g]*fn2[db]*mesq_ubb
	    +fn1[g]*fn2[u]*mesq_dcb
	    +fn1[g]*fn2[s]*mesq_dcb
	    +fn1[g]*fn2[c]*mesq_dcb
	    +fn1[g]*fn2[b]*mesq_dcb
	    +fn1[g]*fn2[u]*mesq_dub
	    +fn1[g]*fn2[s]*mesq_dub
	    +fn1[g]*fn2[c]*mesq_dub
	    +fn1[g]*fn2[b]*mesq_dub
	    +fn1[g]*fn2[u]*mesq_scb
	    +fn1[g]*fn2[d]*mesq_scb
	    +fn1[g]*fn2[c]*mesq_scb
	    +fn1[g]*fn2[b]*mesq_scb
	    +fn1[g]*fn2[u]*mesq_sub
	    +fn1[g]*fn2[d]*mesq_sub
	    +fn1[g]*fn2[c]*mesq_sub
	    +fn1[g]*fn2[b]*mesq_sub
	    +fn1[g]*fn2[u]*mesq_bcb
	    +fn1[g]*fn2[d]*mesq_bcb
	    +fn1[g]*fn2[s]*mesq_bcb
	    +fn1[g]*fn2[c]*mesq_bcb
	    +fn1[g]*fn2[u]*mesq_bub
	    +fn1[g]*fn2[d]*mesq_bub
	    +fn1[g]*fn2[s]*mesq_bub
	    +fn1[g]*fn2[c]*mesq_bub;

	  QBPGN_2 =
	    +fn1[bb]*fn2[g]*mesq_cbd
	    +fn1[sb]*fn2[g]*mesq_cbd
	    +fn1[db]*fn2[g]*mesq_cbd
	    +fn1[ub]*fn2[g]*mesq_cbd
	    +fn1[bb]*fn2[g]*mesq_cbs
	    +fn1[sb]*fn2[g]*mesq_cbs
	    +fn1[db]*fn2[g]*mesq_cbs
	    +fn1[ub]*fn2[g]*mesq_cbs
	    +fn1[bb]*fn2[g]*mesq_cbb
	    +fn1[sb]*fn2[g]*mesq_cbb
	    +fn1[db]*fn2[g]*mesq_cbb
	    +fn1[ub]*fn2[g]*mesq_cbb
	    +fn1[bb]*fn2[g]*mesq_ubd
	    +fn1[cb]*fn2[g]*mesq_ubd
	    +fn1[sb]*fn2[g]*mesq_ubd
	    +fn1[db]*fn2[g]*mesq_ubd
	    +fn1[bb]*fn2[g]*mesq_ubs
	    +fn1[cb]*fn2[g]*mesq_ubs
	    +fn1[sb]*fn2[g]*mesq_ubs
	    +fn1[db]*fn2[g]*mesq_ubs
	    +fn1[bb]*fn2[g]*mesq_ubb
	    +fn1[cb]*fn2[g]*mesq_ubb
	    +fn1[sb]*fn2[g]*mesq_ubb
	    +fn1[db]*fn2[g]*mesq_ubb
	    +fn1[u]*fn2[g]*mesq_dcb
	    +fn1[s]*fn2[g]*mesq_dcb
	    +fn1[c]*fn2[g]*mesq_dcb
	    +fn1[b]*fn2[g]*mesq_dcb
	    +fn1[u]*fn2[g]*mesq_dub
	    +fn1[s]*fn2[g]*mesq_dub
	    +fn1[c]*fn2[g]*mesq_dub
	    +fn1[b]*fn2[g]*mesq_dub
	    +fn1[u]*fn2[g]*mesq_scb
	    +fn1[d]*fn2[g]*mesq_scb
	    +fn1[c]*fn2[g]*mesq_scb
	    +fn1[b]*fn2[g]*mesq_scb
	    +fn1[u]*fn2[g]*mesq_sub
	    +fn1[d]*fn2[g]*mesq_sub
	    +fn1[c]*fn2[g]*mesq_sub
	    +fn1[b]*fn2[g]*mesq_sub
	    +fn1[u]*fn2[g]*mesq_bcb
	    +fn1[d]*fn2[g]*mesq_bcb
	    +fn1[s]*fn2[g]*mesq_bcb
	    +fn1[c]*fn2[g]*mesq_bcb
	    +fn1[u]*fn2[g]*mesq_bub
	    +fn1[d]*fn2[g]*mesq_bub
	    +fn1[s]*fn2[g]*mesq_bub
	    +fn1[c]*fn2[g]*mesq_bub;
	}      
    }
}

double mellinint::integrand2d(int i1, int i2, int sign) //should return a complex value, for the minimal prescription! (the expc:: are not complex-conjugate between pos and neg branches when b is complex)
{
  //cout << "C++ " << i1 << "  " << i2 << "  " << QQBN << endl;
  int i = hcoefficients::index(i1,i2,sign);
  //cout << "C++ " << i1 << "  " << i2 << "  " << GGN << "  " << hcoefficients::Hgg[i] << endl;

  //cout << "C++ " << i1 << "  " << i2 << "  " << sign << "  " << i  << endl;
  //cout << expc::qqb[i] << endl;
  //cout << QGN_1+QGN_2 << " hqg 2d " << hcoefficients::Hqg_1[0] << "  " << hcoefficients::Hqg_2[0] << endl;
  
  if (opts.order == 0)
    return real(QQBN*expc::qqb[i]);

  return
    //contribution starting at LL
    real(QQBN*hcoefficients::Hqqb[i]*expc::qqb[i])

    //contribution starting at NLL
    + real(QGN_1*hcoefficients::Hqg_1[i]*expc::qg_1[i]) + real(QGN_2*hcoefficients::Hqg_2[i]*expc::qg_2[i])
    
    //contributions starting at NNLL
    + real(GGN*hcoefficients::Hgg[i]*expc::gg[i])
    //+ QQN*hcoefficients::Hqq[i]
    + real(QQN_1*hcoefficients::Hqq_1[i]*expc::qq_1[i]) + real(QQN_2*hcoefficients::Hqq_2[i]*expc::qq_2[i]) //Bug fix in DYRES (I believe this is more correct, since it accounts for which leg undergoes the q -> qb or qb -> q transformation)
    + real(QQPN_1*hcoefficients::Hqqp_1[i]*expc::qq_1[i]) + real(QQPN_2*hcoefficients::Hqqp_2[i]*expc::qq_2[i])
    + real(QQBPN_1*hcoefficients::Hqqp_1[i]*expc::qq_1[i]) + real(QQBPN_2*hcoefficients::Hqqp_2[i]*expc::qq_2[i]);
  
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

complex <double> mellinint::calc1d()
{
  complex <double> fun = 0.;
  int idx;
  for (int i = 0; i < mdim; i++)
    {
      //Positive branch
      complex <double> pos;
      idx = anomalous::index(i,mesq::positive);
      pdfevol::retrieve1d_pos(i);
      pdf_mesq_expy(i,i,mesq::positive);

      if (opts.order == 0)
	pos = QQBN*expc::qqb[idx];
      else
	{
	  //cout << i << "  " << muf::qqb[i] << endl;
	  pos =
	    //contribution starting at LL
	    QQBN             *(hcoeff::Hqqb[i]+muf::qqb[i])*expc::qqb[idx]
	    
	    //contribution starting at NLL
	    + (QGN_1+QGN_2)  *(hcoeff::Hqg[i]+muf::qg[i]) *expc::qg[idx]
    
	    //contributions starting at NNLL
	    + GGN              *(hcoeff::Hgg[i]  +muf::gg[i]  )*expc::gg[idx]
	    + (QQN_1+QQN_2)    *(hcoeff::Hqq[i]  +muf::qq[i]  )*expc::qq[idx]
	    + (QQPN_1+QQPN_2)  *(hcoeff::Hqqp[i] +muf::qqp[i] )*expc::qqp[idx]
	    + (QQBPN_1+QQBPN_2)*(hcoeff::Hqqbp[i]+muf::qqbp[i])*expc::qqbp[idx]
	    
	    //contributions starting at NNNLL
	    + QQBN_nfz         *(hcoeff::Hqqb_nfz             )*expc::qqb[idx]
	    + (QBGN_1+QBGN_2)  *(hcoeff::Hqbg[i] +muf::qbg[i] )*expc::qbg[idx]
	    + (QPGN_1+QPGN_2)  *(hcoeff::Hqpg[i] +muf::qpg[i] )*expc::qpg[idx]
	    + (QBPGN_1+QBPGN_2)*(hcoeff::Hqbpg[i]+muf::qbpg[i])*expc::qbpg[idx];


	}
      //if (opts.order == 0)
      //cout << i << " positive " << " QQBN " << QQBN
      //	   << " expc::qqb[idx] " << expc::qqb[idx]
      //	   << " hcoeff::Hqqb[i] " << hcoeff::Hqqb[i]
      //	   << endl;
      //cout << pdfevol::fn1[u ] << "  " << pdfevol::fn2[ub] << "  "
      //	   << mesq::mesqij_expy[mesq::index(0,i,i,mesq::positive)] << endl;
      //else
      //	{
      //cout << i << " positive " << setprecision(16) << " QQBN " << QQBN << " Hqqb " << hcoeff::Hqqb[i]+muf::qqb[i] << " expc::qqb[idx] " << expc::qqb[idx] << endl;
      //cout << " QGN_1 " << QGN_1 << " Hqg " << hcoeff::Hqg[i] << " expc::qg[idx] " << expc::qg[idx] << endl;
      //	}
      //complex <double> pos = 	    QQBN             *(hcoeff::Hqqb[i]+muf::qqb[i])*expc::qqb[idx]
      //+ (QGN_1+QGN_2)  *(hcoeff::Hqg[i]+muf::qg[i]) *expc::qg[idx]
      //;
      //cout << " QQBN_nfx " << QQBN_nfz
      //	   << hcoeff::Hqqb_nfz
      //	   << endl;
      //Negative branch
      complex <double> neg;
      //if (opts.mellininv == 1 && i == 0)
      //continue;
      idx = anomalous::index(i,mesq::negative);
      pdfevol::retrieve1d_neg();
      pdf_mesq_expy(i,i,mesq::negative);

      if (opts.order == 0)
	neg = QQBN*expc::qqb[idx];
      else
	{
	  neg =
	    //contribution starting at LL
	    QQBN             *conj(hcoeff::Hqqb[i]+muf::qqb[i])*expc::qqb[idx]

	    //contribution starting at NLL
	    + (QGN_1+QGN_2)  *conj(hcoeff::Hqg[i]+muf::qg[i]) *expc::qg[idx]
    
	    //contributions starting at NNLL
	    + GGN              *conj(hcoeff::Hgg[i]  +muf::gg[i]  )*expc::gg[idx]
	    + (QQN_1+QQN_2)    *conj(hcoeff::Hqq[i]  +muf::qq[i]  )*expc::qq[idx]
	    + (QQPN_1+QQPN_2)  *conj(hcoeff::Hqqp[i] +muf::qqp[i] )*expc::qqp[idx]
	    + (QQBPN_1+QQBPN_2)*conj(hcoeff::Hqqbp[i]+muf::qqbp[i])*expc::qqbp[idx]
	    
	    //contributions starting at NNNLL
	    + QQBN_nfz         *conj(hcoeff::Hqqb_nfz             )*expc::qqb[idx]
	    + (QBGN_1+QBGN_2)  *conj(hcoeff::Hqbg[i] +muf::qbg[i] )*expc::qbg[idx]
	    + (QPGN_1+QPGN_2)  *conj(hcoeff::Hqpg[i] +muf::qpg[i] )*expc::qpg[idx]
	    + (QBPGN_1+QBPGN_2)*conj(hcoeff::Hqbpg[i]+muf::qbpg[i])*expc::qbpg[idx];
	}

      //if (opts.order == 0)
      //	cout << i << " negative " << " QQBN " << -QQBN  << " expc::qqb[idx] " << expc::qqb[idx] << endl;
      //else
      //	{
      //cout << i << " negative " << " QQBN " << -QQBN << " Hqqb " << conj(hcoeff::Hqqb[i]+muf::qqb[i]) << " expc::qqb[idx] " << expc::qqb[idx] << endl;
      //cout << " QGN_1 " << -QGN_1 << " Hqg " << conj(hcoeff::Hqg[i]) << " expc::qg[idx] " << expc::qg[idx] << endl;
      //	}

      //      complex <double> neg =QQBN             *conj(hcoeff::Hqqb[i]+muf::qqb[i])*expc::qqb[idx]
      //	+ (QGN_1+QGN_2)  *conj(hcoeff::Hqg[i]+muf::qg[i]) *expc::qg[idx]
      //	;
      
      //Check carefully in the real case (bprescription = 0) that all pieces are complex conjugate between positive and negative branches!!!
      fun += pos-neg;
      //cout << i << " pos - neg " << pos-neg << " fun " << fun << endl;
    }
  //return real(fun);
  //cout << QGN_1+QGN_2 << "  " << "hqg 1d " << hcoeff::Hqg[0] << endl;
  //cout << " fun " << fun/2. << endl;

  if (opts.bprescription == 0)
    return real(fun/2.);
  else
    return fun/2.;
}

complex <double> mellinint::calc2d()
{
  complex <double> fun = 0.;

  //#pragma omp parallel for reduction(+:fun) num_threads(opts.mellincores) copyin(creno_,mesq::mesqij_expy,hcoefficients::Hqqb,hcoefficients::Hqg_1,hcoefficients::Hqg_2,hcoefficients::Hgg,hcoefficients::Hqq,hcoefficients::Hqq_1,hcoefficients::Hqq_2,hcoefficients::Hqqp_1,hcoefficients::Hqqp_2)
  for (int i1 = 0; i1 < mdim; i1++)
    {
      pdfevol::retrieve_beam1(i1);
      for (int i2 = 0; i2 < mdim; i2++)
	{
	  //Positive branch
	  pdfevol::retrieve_beam2_pos(i2);
	  pdf_mesq_expy(i1,i2,mesq::positive);
	  int ii = hcoeff::index(i1,i2,mesq::positive);

	  if (opts.order == 0)
	    fun += QQBN*expc::qqb[ii];
	  else
	    {
	      fun +=
		//contribution starting at LL
		QQBN*hcoeff::Hqqb[ii]*expc::qqb[ii]

		//contribution starting at NLL
		+ QGN_1*hcoeff::Hqg_1[ii]*expc::qg_1[ii] + QGN_2*hcoeff::Hqg_2[ii]*expc::qg_2[ii]
    
		//contributions starting at NNLL
		+ GGN*hcoeff::Hgg[ii]*expc::gg[ii]
		+ QQN_1*hcoeff::Hqq_1[ii]*expc::qq_1[ii] + QQN_2*hcoeff::Hqq_2[ii]*expc::qq_2[ii]
		+ QQPN_1*hcoeff::Hqqp_1[ii]*expc::qqp_1[ii] + QQPN_2*hcoeff::Hqqp_2[ii]*expc::qqp_2[ii]
		+ QQBPN_1*hcoeff::Hqqbp_1[ii]*expc::qqbp_1[ii] + QQBPN_2*hcoeff::Hqqbp_2[ii]*expc::qqbp_2[ii]

		//contributions starting at NNNLL
		+ QQBN_nfz*hcoeff::Hqqb_nfz*expc::qqb[ii]
		+ QBGN_1*hcoeff::Hqbg_1[ii]*expc::qbg_1[ii]+QBGN_2  *hcoeff::Hqbg_2[ii]*expc::qbg_2[ii]
		+ QPGN_1*hcoeff::Hqbg_1[ii]*expc::qbg_1[ii]+QPGN_2  *hcoeff::Hqpg_2[ii]*expc::qpg_2[ii]   //!!! bug on the first leg qbg->qpg 
		+ QBPGN_1*hcoeff::Hqbg_1[ii]*expc::qbg_1[ii]+QBPGN_2*hcoeff::Hqbpg_2[ii]*expc::qbpg_2[ii] //!!! bug on the first leg qbg->qbpg
		;
	    }	  
	      
//	  if (opts.order == 0)
//	    cout << i1 << "  " << i2 << " positive " << " QQBN " << QQBN << " expc::qqb[ii] " << expc::qqb[ii]
//		 << " mesq " << mesq::mesqij_expy[mesq::index(0,i1,i2,mesq::positive)]
//		 << " fn1 " << pdfevol::fn1[g]
//		 << " fn2 " << pdfevol::fn2[g]
//		 << endl;
//	  cout << i1 << "  " << i2 << " positive " << " QQBN " << QQBN << " expc::qqb[ii] " << expc::qqb[ii]
//	       << " mesq " << mesq::mesqij_expy[mesq::index(0,i1,i2,mesq::positive)]
//	       << " fn2 " << pdfevol::fn2[ub]
//	       << endl;
//	    cout << i1 << "  " << i2 << " positive "
//		 << " QQBN "    << hcoeff::Hqqb[ii]
//		 << " QGN_1   " << hcoeff::Hqg_1[ii]
//		 << " QQN_1   " << hcoeff::Hqq_1[ii]
//		 << " QQPN_1  " << hcoeff::Hqqp_1[ii]
//		 << " QQBPN_1 " << hcoeff::Hqqbp_1[ii]
//		 << " QBGN_1  " << hcoeff::Hqbg_1[ii]
//		 << " QPGN_1  " << hcoeff::Hqpg_1[ii]
//		 << " QBPGN_1 " << hcoeff::Hqbpg_1[ii]
//		 << " GGN     " << hcoeff::Hgg[ii]
//		 << endl;
	      
	  //Negative branch
	  pdfevol::retrieve_beam2_neg();
	  pdf_mesq_expy(i1,i2,mesq::negative);
	  ii = hcoeff::index(i1,i2,mesq::negative);

	  if (opts.order == 0)
	    fun -= QQBN*expc::qqb[ii];
	  else
	    {
	      fun -=
		//contribution starting at LL
		QQBN*hcoeff::Hqqb[ii]*expc::qqb[ii]

		//contribution starting at NLL
		+ QGN_1*hcoeff::Hqg_1[ii]*expc::qg_1[ii] + QGN_2*hcoeff::Hqg_2[ii]*expc::qg_2[ii]
    		
		//contributions starting at NNLL
		+ GGN*hcoeff::Hgg[ii]*expc::gg[ii]
		+ QQN_1*hcoeff::Hqq_1[ii]*expc::qq_1[ii] + QQN_2*hcoeff::Hqq_2[ii]*expc::qq_2[ii]
		+ QQPN_1*hcoeff::Hqqp_1[ii]*expc::qqp_1[ii] + QQPN_2*hcoeff::Hqqp_2[ii]*expc::qqp_2[ii]
		+ QQBPN_1*hcoeff::Hqqbp_1[ii]*expc::qqbp_1[ii] + QQBPN_2*hcoeff::Hqqbp_2[ii]*expc::qqbp_2[ii]
		
		//contributions starting at NNNLL
		+ QQBN_nfz*hcoeff::Hqqb_nfz*expc::qqb[ii]
		+ QBGN_1*hcoeff::Hqbg_1[ii]*expc::qbg_1[ii]+QBGN_2  *hcoeff::Hqbg_2[ii]*expc::qbg_2[ii]    
		+ QPGN_1*hcoeff::Hqbg_1[ii]*expc::qbg_1[ii]+QPGN_2  *hcoeff::Hqpg_2[ii]*expc::qpg_2[ii]	   //!!! bug on the first leg qbg->qpg
		+ QBPGN_1*hcoeff::Hqbg_1[ii]*expc::qbg_1[ii]+QBPGN_2*hcoeff::Hqbpg_2[ii]*expc::qbpg_2[ii]  //!!! bug on the first leg qbg->qbpg
		;
	      
	    }
//	  if (opts.order == 0)
//	    cout << i1 << "  " << i2 << " negative " << " QQBN " << -QQBN  << " expc::qqb[ii] " << expc::qqb[ii]
//		 << " mesq " << mesq::mesqij_expy[mesq::index(0,i1,i2,mesq::negative)]
//		 << " fn2 " << pdfevol::fn2[ub]
//		 << endl;
//	    cout << i1 << "  " << i2 << " negative " << " QQBN " << -QQBN  << " expc::qqb[ii] " << expc::qqb[ii]
//		 << " mesq " << mesq::mesqij_expy[mesq::index(0,i1,i2,mesq::negative)]
//		 << " fn2 " << pdfevol::fn2[ub]
//		 << endl;
	}
    }

  if (opts.bprescription == 0)
    return -real(fun/2.);
  else
    return -fun/2.;
}  


complex <double> mellinint::calc1d_muf()
{
  complex <double> fun = 0.;
  int idx;
  for (int i = 0; i < mdim; i++)
    {
      //Positive branch
      complex <double> pos = 0.;
      idx = anomalous::index(i,mesq::positive);
      pdfevol::retrieve1d_pos(i);
      pdf_mesq_expy(i,i,mesq::positive);

      //contribution starting at LL
      pos +=
	+ QQBN             *(hcoeff::Hqqb[i]+muf::qqb[idx])*expc::qqb[idx];
      //contribution starting at NLL
      if (opts.order_hcoef >= 1)
	pos +=
	  +(QGN_1+QGN_2)  *(hcoeff::Hqg[i]    *expc::qg[idx]   + muf::qg[idx]  *expc::qqb[idx]);
    
      //contributions starting at NNLL
      if (opts.order_hcoef >= 2)
	pos +=
	  + GGN              *(hcoeff::Hgg[i]  *expc::gg[idx]   + muf::gg[idx]  *expc::qqb[idx])
	  + (QQN_1+QQN_2)    *(hcoeff::Hqq[i]  *expc::qq[idx]   + muf::qq[idx]  *expc::qqb[idx])
	  + (QQPN_1+QQPN_2)  *(hcoeff::Hqqp[i] *expc::qqp[idx]  + muf::qqp[idx] *expc::qqb[idx])
	  + (QQBPN_1+QQBPN_2)*(hcoeff::Hqqbp[i]*expc::qqbp[idx] + muf::qqbp[idx]*expc::qqb[idx]);
	    
	    //contributions starting at NNNLL
      if (opts.order_hcoef >= 3)
	pos +=
	  + QQBN_nfz         *hcoeff::Hqqb_nfz *expc::qqb[idx]
	  + (QBGN_1+QBGN_2)  *(hcoeff::Hqbg[i] *expc::qbg[idx]  + muf::qbg[idx] *expc::qqb[idx])
	  + (QPGN_1+QPGN_2)  *(hcoeff::Hqpg[i] *expc::qpg[idx]  + muf::qpg[idx] *expc::qqb[idx])
	  + (QBPGN_1+QBPGN_2)*(hcoeff::Hqbpg[i]*expc::qbpg[idx] + muf::qbpg[idx]*expc::qqb[idx]);

      //Negative branch
      complex <double> neg = 0.;
      idx = anomalous::index(i,mesq::negative);
      pdfevol::retrieve1d_neg();
      pdf_mesq_expy(i,i,mesq::negative);

      //contribution starting at LL
      neg +=
	QQBN               *(conj(hcoeff::Hqqb[i])+muf::qqb[idx])*expc::qqb[idx];

      //contribution starting at NLL
      if (opts.order >= 1)
	neg +=
	  + (QGN_1+QGN_2)    *(conj(hcoeff::Hqg[i]  ) *expc::qg[idx]   + muf::qg[idx]*expc::qqb[idx]);

      //contributions starting at NNLL
      if (opts.order >= 2)
	neg +=
	  + GGN              *(conj(hcoeff::Hgg[i]  )*expc::gg[idx]    + muf::gg[idx]  *expc::qqb[idx])
	  + (QQN_1+QQN_2)    *(conj(hcoeff::Hqq[i]  )*expc::qq[idx]    + muf::qq[idx]  *expc::qqb[idx])
	  + (QQPN_1+QQPN_2)  *(conj(hcoeff::Hqqp[i] )*expc::qqp[idx]   + muf::qqp[idx] *expc::qqb[idx])
	  + (QQBPN_1+QQBPN_2)*(conj(hcoeff::Hqqbp[i])*expc::qqbp[idx]  + muf::qqbp[idx]*expc::qqb[idx]);
	    
      //contributions starting at NNNLL
      if (opts.order >= 3)
	neg +=
	  + QQBN_nfz         * conj(hcoeff::Hqqb_nfz)*expc::qqb[idx]
	  + (QBGN_1+QBGN_2)  *(conj(hcoeff::Hqbg[i] )*expc::qbg[idx]   + muf::qbg[idx] *expc::qqb[idx])
	  + (QPGN_1+QPGN_2)  *(conj(hcoeff::Hqpg[i] )*expc::qpg[idx]   + muf::qpg[idx] *expc::qqb[idx])
	  + (QBPGN_1+QBPGN_2)*(conj(hcoeff::Hqbpg[i])*expc::qbpg[idx]  + muf::qbpg[idx]*expc::qqb[idx]);

      fun += pos-neg;
    }

  if (opts.bprescription == 0)
    return real(fun/2.);
  else
    return fun/2.;
}

complex <double> mellinint::calc2d_muf()
{
  complex <double> fun = 0.;
  for (int i1 = 0; i1 < mdim; i1++)
    {
      pdfevol::retrieve_beam1(i1);
      for (int i2 = 0; i2 < mdim; i2++)
	{
	  //Positive branch
	  pdfevol::retrieve_beam2_pos(i2);
	  pdf_mesq_expy(i1,i2,mesq::positive);
	  int ii = hcoeff::index(i1,i2,mesq::positive);

	  //contribution starting at LL
	  fun +=
	    QQBN*(hcoeff::Hqqb[ii]+muf::qqb[ii])*expc::qqb[ii];

	  //contribution starting at NLL
	  if (opts.order_hcoef >= 1)
	    fun +=
	      + QGN_1  *(hcoeff::Hqg_1[ii]*expc::qg_1[ii]    +muf::qg_1[ii]  *expc::qqb[ii]) + QGN_2  *(hcoeff::Hqg_2[ii]  *expc::qg_2[ii]  +muf::qg_2[ii]*expc::qqb[ii]);

	  //contributions starting at NNLL
	  if (opts.order_hcoef >= 2)
	    fun +=
	      + GGN    *(hcoeff::Hgg[ii]    *expc::gg[ii]    +muf::gg[ii]    *expc::qqb[ii])
	      + QQN_1  *(hcoeff::Hqq_1[ii]  *expc::qq_1[ii]  +muf::qq_1[ii]  *expc::qqb[ii]) + QQN_2  *(hcoeff::Hqq_2[ii]  *expc::qq_2[ii]  +muf::qq_2[ii]  *expc::qqb[ii])
	      + QQPN_1 *(hcoeff::Hqqp_1[ii] *expc::qqp_1[ii] +muf::qqp_1[ii] *expc::qqb[ii]) + QQPN_2 *(hcoeff::Hqqp_2[ii] *expc::qqp_2[ii] +muf::qqp_2[ii] *expc::qqb[ii])
	      + QQBPN_1*(hcoeff::Hqqbp_1[ii]*expc::qqbp_1[ii]+muf::qqbp_1[ii]*expc::qqb[ii]) + QQBPN_2*(hcoeff::Hqqbp_2[ii]*expc::qqbp_2[ii]+muf::qqbp_2[ii]*expc::qqb[ii]);

	  //contributions starting at NNNLL
	  if (opts.order_hcoef >= 3)
	    fun +=
	      + QQBN_nfz*hcoeff::Hqqb_nfz*expc::qqb[ii]
	      + QBGN_1  *(hcoeff::Hqbg_1[ii]*expc::qbg_1[ii] +muf::qbg_1[ii]*expc::qqb[ii])+QBGN_2  *(hcoeff::Hqbg_2[ii] *expc::qbg_2[ii] +muf::qbg_2[ii] *expc::qqb[ii] )
	      + QPGN_1  *(hcoeff::Hqbg_1[ii]*expc::qpg_1[ii] +muf::qbg_1[ii]*expc::qqb[ii])+QPGN_2  *(hcoeff::Hqpg_2[ii] *expc::qpg_2[ii] +muf::qpg_2[ii] *expc::qqb[ii] )
	      + QBPGN_1 *(hcoeff::Hqbg_1[ii]*expc::qbpg_1[ii]+muf::qbg_1[ii]*expc::qqb[ii])+QBPGN_2 *(hcoeff::Hqbpg_2[ii]*expc::qbpg_2[ii]+muf::qbpg_2[ii]*expc::qqb[ii])
	      ;
	      
	  //Negative branch
	  pdfevol::retrieve_beam2_neg();
	  pdf_mesq_expy(i1,i2,mesq::negative);
	  ii = hcoeff::index(i1,i2,mesq::negative);

	  //contribution starting at LL
	  fun -=
	    QQBN*(hcoeff::Hqqb[ii]+muf::qqb[ii])*expc::qqb[ii];

	  //contribution starting at NLL
	  if (opts.order_hcoef >= 1)
	    fun -=
	      + QGN_1  *(hcoeff::Hqg_1[ii]*expc::qg_1[ii]    +muf::qg_1[ii]  *expc::qqb[ii]) + QGN_2  *(hcoeff::Hqg_2[ii]  *expc::qg_2[ii]  +muf::qg_2[ii]*expc::qqb[ii]);
    		
	  //contributions starting at NNLL
	  if (opts.order_hcoef >= 2)
	    fun -=
	      + GGN    *(hcoeff::Hgg[ii]    *expc::gg[ii]    +muf::gg[ii]    *expc::qqb[ii])
	      + QQN_1  *(hcoeff::Hqq_1[ii]  *expc::qq_1[ii]  +muf::qq_1[ii]  *expc::qqb[ii]) + QQN_2  *(hcoeff::Hqq_2[ii]  *expc::qq_2[ii]  +muf::qq_2[ii]  *expc::qqb[ii])
	      + QQPN_1 *(hcoeff::Hqqp_1[ii] *expc::qqp_1[ii] +muf::qqp_1[ii] *expc::qqb[ii]) + QQPN_2 *(hcoeff::Hqqp_2[ii] *expc::qqp_2[ii] +muf::qqp_2[ii] *expc::qqb[ii])
	      + QQBPN_1*(hcoeff::Hqqbp_1[ii]*expc::qqbp_1[ii]+muf::qqbp_1[ii]*expc::qqb[ii]) + QQBPN_2*(hcoeff::Hqqbp_2[ii]*expc::qqbp_2[ii]+muf::qqbp_2[ii]*expc::qqb[ii]);

	  //contributions starting at NNNLL
	  if (opts.order_hcoef >= 3)
	    fun -=
	      + QQBN_nfz*hcoeff::Hqqb_nfz*expc::qqb[ii]
	      + QBGN_1  *(hcoeff::Hqbg_1[ii]*expc::qbg_1[ii] +muf::qbg_1[ii]*expc::qqb[ii])+QBGN_2  *(hcoeff::Hqbg_2[ii] *expc::qbg_2[ii] +muf::qbg_2[ii] *expc::qqb[ii] )
	      + QPGN_1  *(hcoeff::Hqbg_1[ii]*expc::qpg_1[ii] +muf::qbg_1[ii]*expc::qqb[ii])+QPGN_2  *(hcoeff::Hqpg_2[ii] *expc::qpg_2[ii] +muf::qpg_2[ii] *expc::qqb[ii] )
	      + QBPGN_1 *(hcoeff::Hqbg_1[ii]*expc::qbpg_1[ii]+muf::qbg_1[ii]*expc::qqb[ii])+QBPGN_2 *(hcoeff::Hqbpg_2[ii]*expc::qbpg_2[ii]+muf::qbpg_2[ii]*expc::qqb[ii])
	      ;
	}
    }

  if (opts.bprescription == 0)
    return -real(fun/2.);
  else
    return -fun/2.;
}  
