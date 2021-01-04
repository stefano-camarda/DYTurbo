#include "rapint.h"
#include "settings.h"
#include "gaussrules.h"
#include "mellinint.h"
#include "omegaintegr.h"
#include "phasespace.h"
#include "pdfevol.h"
#include <iostream>

int rapint::ydim;

complex <double> *rapint::cfpy;
complex <double> *rapint::cfmy;

complex <double> *rapint::Ith0p;
complex <double> *rapint::Ith1p;
complex <double> *rapint::Ith2p;
complex <double> *rapint::Ith0m;
complex <double> *rapint::Ith1m;
complex <double> *rapint::Ith2m;

void rapint::init()
{
  ydim = opts.yrule*opts.yintervals;

  //variables to cache the rapidity dependent exponentials
  cfpy = new complex <double>[mellinint::mdim*mellinint::mdim*ydim];
  cfmy = new complex <double>[mellinint::mdim*mellinint::mdim*ydim];
}

void rapint::release()
{
  delete[] cfpy;
  delete[] cfmy;
}

void rapint::cache(double ymin, double ymax)
{
  mellinint::update();
  //cache the rapidity dependent exponentials
  //The caching is done for the numerical integration, but can be done also for the analytical integration (need to calculate only the ymin and ymax exp)
  for (int i=0; i < opts.yintervals; i++)
    {
      double ya = ymin+(ymax-ymin)*i/opts.yintervals;
      double yb = ymin+(ymax-ymin)*(i+1)/opts.yintervals;
      double xc=0.5*(ya+yb);
      double xm=0.5*(yb-ya);
      for (int j=0; j < opts.yrule; j++)
	{
	  double y=xc+xm*gr::xxx[opts.yrule-1][j];
	  for (int i1 = 0; i1 < mellinint::mdim; i1++)
	    for (int i2 = 0; i2 < mellinint::mdim; i2++)
	      {
		complex <double> diffnpp = mellinint::Np_1[i1]-mellinint::Np_2[i2];
		complex <double> diffnpm = mellinint::Np_1[i1]-mellinint::Nm_2[i2];
		cfpy[index(i,j,i1,i2)]=exp(-diffnpp*y)*mellinint::wn_1[i1]*     mellinint::wn_2[i2] *gr::www[opts.yrule-1][j]*xm;
		cfmy[index(i,j,i1,i2)]=exp(-diffnpm*y)*mellinint::wn_1[i1]*conj(mellinint::wn_2[i2])*gr::www[opts.yrule-1][j]*xm;
	      }
	}
    }
}

void rapint::allocate()
{
  //allocate memory
  Ith0p = new complex <double>[mellinint::mdim*mellinint::mdim];
  Ith1p = new complex <double>[mellinint::mdim*mellinint::mdim];
  Ith2p = new complex <double>[mellinint::mdim*mellinint::mdim];
  Ith0m = new complex <double>[mellinint::mdim*mellinint::mdim];
  Ith1m = new complex <double>[mellinint::mdim*mellinint::mdim];
  Ith2m = new complex <double>[mellinint::mdim*mellinint::mdim];
}
void rapint::integrate(double ymin, double ymax, double m, int helicity)
{
  //integration results are stored in Ithxx

  //Initialize integrals
  for (int i1 = 0; i1 < mellinint::mdim; i1++)
    for (int i2 = 0; i2 < mellinint::mdim; i2++)
      {
	Ith0p[mellinint::index(i1,i2)] = 0.;
	Ith1p[mellinint::index(i1,i2)] = 0.;
	Ith2p[mellinint::index(i1,i2)] = 0.;
	Ith0m[mellinint::index(i1,i2)] = 0.;
	Ith1m[mellinint::index(i1,i2)] = 0.;
	Ith2m[mellinint::index(i1,i2)] = 0.;
      }
  
  double q2 = m*m;
  double bjx= q2/pow(opts.sroot,2);
  double ax = log(bjx);
  //double ax1 = (ax+2*y)/2.; 
  //double ax2 = (ax-2*y)/2.; 

  //Notice that in all the expressions, the dependence on I1, I2 is of the type I1+I2 or I1-I2
  //Does it mean that the double loops on I1 I2 can be reduced to single loops on I1+I2 and I1-I2 to calculate all the exponentials?
  //Is it true for both analytical and numerical integration?
  //No, because Np[i1]+Np[i2] != Np[i1+i2]

  //If there are no cuts on the leptons, calculate the integrals analitically
  //Int_ymin^ymax (CCp/M_PI)^2 * exp(Np(i1)*ax1) * exp(Np(i2)*ax2)
  if (!opts.makecuts)
    //Analytical integration
    for (int i1 = 0; i1 < mellinint::mdim; i1++)
      for (int i2 = 0; i2 < mellinint::mdim; i2++)
	{
	  complex <double> yintp, yintm;

	  //split between negative and positive y, so as to allow costh boundaries, and account for y+- sign flip on costh
	  double ymn, ymx;
	  ymn = min(0., ymin);
	  ymx = min(0., ymax);
	  if (ymn < 0)
	    {
	      phasespace::set_y((ymn+ymx)/2.);
	      complex <double> sumnpp = mellinint::Np_1[i1]+mellinint::Np_2[i2];
	      complex <double> diffnpp = mellinint::Np_1[i1]-mellinint::Np_2[i2];
	      //if (i1 == i2)
	      if (mellinint::Np_1[i1] == mellinint::Np_2[i2])
		yintp=pow(mellinint::CCp/M_PI,2)*exp(-sumnpp*ax/2.)*(ymx-ymn);
	      else
		//yintp=1./(-diffnpp)*pow(mellinint::CCp/M_PI,2)*exp(-sumnpp*ax/2.)*(exp(-diffnpp*ymx)-exp(-diffnpp*ymn));
		//yintp=1./(-diffnpp)*pow(mellinint::CCp/M_PI,2)*(exp(-diffnpp*ymx-sumnpp*ax/2.)-exp(-diffnpp*ymn-sumnpp*ax/2.));
		yintp=1./(-diffnpp)*pow(mellinint::CCp/M_PI,2)
		  *(exp(-mellinint::Np_1[i1]*(ax/2.+ymx))*exp(-mellinint::Np_2[i2]*(ax/2.-ymx))
		    -exp(-mellinint::Np_1[i1]*(ax/2.+ymn))*exp(-mellinint::Np_2[i2]*(ax/2.-ymn)));

//	      double ylim = 0.5*log(pow(opts.sroot,2)/phasespace::m2);
//	      if (phasespace::ymin <= -ylim) //when ymn = -ylim, then x2_max = 1 and (ax/2.-ymn) = log(x2_max) = 0
//		yintp=1./(-diffnpp)*pow(mellinint::CCp/M_PI,2)
//		  *(exp(-mellinint::Np_1[i1]*(ax/2.+ymx))*exp(-mellinint::Np_2[i2]*(ax/2.-ymx))
//		    -exp(-mellinint::Np_1[i1]*(ax/2.+ymn)));

	      complex <double> sumnpm = mellinint::Np_1[i1]+mellinint::Nm_2[i2];
	      complex <double> diffnpm = mellinint::Np_1[i1]-mellinint::Nm_2[i2];
	      //if (i1 == 0 && i2 == 0 && opts.mellininv == 1 && !mellinint::midpoint)
	      if (mellinint::Np_1[i1] == mellinint::Nm_2[i2]) // this never happens --> It happens for Talbot!
		yintm=(mellinint::CCp/M_PI)*(mellinint::CCm/M_PI)*exp(-sumnpm*ax/2.)*(ymx-ymn);
	      else
		//yintm=1./(-diffnpm)*(mellinint::CCp/M_PI)*(mellinint::CCm/M_PI)*exp(-sumnpm*ax/2.)*(exp(-diffnpm*ymx)-exp(-diffnpm*ymn));
		//yintm=1./(-diffnpm)*(mellinint::CCp/M_PI)*(mellinint::CCm/M_PI)*(exp(-diffnpm*ymx-sumnpm*ax/2.)-exp(-diffnpm*ymn-sumnpm*ax/2.));
		yintm=1./(-diffnpm)*(mellinint::CCp/M_PI)*(mellinint::CCm/M_PI)
		  *(exp(-mellinint::Np_1[i1]*(ax/2.+ymx))*exp(-mellinint::Nm_2[i2]*(ax/2.-ymx))
		    -exp(-mellinint::Np_1[i1]*(ax/2.+ymn))*exp(-mellinint::Nm_2[i2]*(ax/2.-ymn)));

//	      if (phasespace::ymin <= -ylim)  //when ymn = -ylim, then x2_max = 1 and (ax/2.-ymn) = log(x2_max) = 0
//		yintm=1./(-diffnpm)*(mellinint::CCp/M_PI)*(mellinint::CCm/M_PI)
//		  *(exp(-mellinint::Np_1[i1]*(ax/2.+ymx))*exp(-mellinint::Nm_2[i2]*(ax/2.-ymx))
//		    -exp(-mellinint::Np_1[i1]*(ax/2.+ymn)));

	      double cthmom0, cthmom1, cthmom2;
	      omegaintegr::cthmoments(cthmom0,cthmom1,cthmom2);
	      if (helicity >= 0)
		{
		  if (helicity == 0)
		    {
		      cthmom0 = 4./3.; //20./3. * (0.5*2. - 1.5*2./3.  ) +2./3.*2.;
		      cthmom1 = 0.;
		      cthmom2 = -4./3.; //20./3. * (0.5*2./3.-1.5*2./5.  ) +2./3.*2./3.;
		    }
		  else if (helicity == 4)
		    {
		      cthmom0=0.;
		      //cthmom1=4*cthmom2*(y < 0 ? -1. : 1.);
		      cthmom1=-4*cthmom2;
		      cthmom2=0.;
		    }
		  else
		    {
		      cthmom0=0.;
		      cthmom1=0.;
		      cthmom2=0.;
		    }
		}
	      
	      Ith0p[mellinint::index(i1,i2)] = cthmom0*yintp*mellinint::wn_1[i1]*mellinint::wn_2[i2];
	      Ith1p[mellinint::index(i1,i2)] = cthmom1*yintp*mellinint::wn_1[i1]*mellinint::wn_2[i2];
	      Ith2p[mellinint::index(i1,i2)] = cthmom2*yintp*mellinint::wn_1[i1]*mellinint::wn_2[i2];
	      Ith0m[mellinint::index(i1,i2)] = cthmom0*yintm*mellinint::wn_1[i1]*conj(mellinint::wn_2[i2]);
	      Ith1m[mellinint::index(i1,i2)] = cthmom1*yintm*mellinint::wn_1[i1]*conj(mellinint::wn_2[i2]);
	      Ith2m[mellinint::index(i1,i2)] = cthmom2*yintm*mellinint::wn_1[i1]*conj(mellinint::wn_2[i2]);
	    }

	  ymn = max(0., ymin);
	  ymx = max(0., ymax);
	  if (ymx > 0.)
	    {
	      phasespace::set_y((ymn+ymx)/2.);
	      complex <double> sumnpp = mellinint::Np_1[i1]+mellinint::Np_2[i2];
	      complex <double> diffnpp = mellinint::Np_1[i1]-mellinint::Np_2[i2];
	      //if (i1 == i2)
	      if (mellinint::Np_1[i1] == mellinint::Np_2[i2])
		yintp=pow(mellinint::CCp/M_PI,2)*exp(-sumnpp*ax/2.)*(ymx-ymn);
	      else
		//yintp=1./(-diffnpp)*pow(mellinint::CCp/M_PI,2)*exp(-sumnpp*ax/2.)*(exp(-diffnpp*ymx)-exp(-diffnpp*ymn));
		//yintp=1./(-diffnpp)*pow(mellinint::CCp/M_PI,2)*(exp(-diffnpp*ymx-sumnpp*ax/2.)-exp(-diffnpp*ymn-sumnpp*ax/2.));
		yintp=1./(-diffnpp)*pow(mellinint::CCp/M_PI,2)
		  *(exp(-mellinint::Np_1[i1]*(ax/2.+ymx))*exp(-mellinint::Np_2[i2]*(ax/2.-ymx))
		    -exp(-mellinint::Np_1[i1]*(ax/2.+ymn))*exp(-mellinint::Np_2[i2]*(ax/2.-ymn)));
	      
//	      double ylim = 0.5*log(pow(opts.sroot,2)/phasespace::m2);
//	      if (phasespace::ymax >= ylim) //when ymx = ylim, then x1_max = 1 and (ax/2.+ymx) = log(x1_max) = 0
//		yintp=1./(-diffnpp)*pow(mellinint::CCp/M_PI,2)
//		  *(exp(-mellinint::Np_2[i2]*(ax/2.-ymx))
//		    -exp(-mellinint::Np_1[i1]*(ax/2.+ymn))*exp(-mellinint::Np_2[i2]*(ax/2.-ymn)));
		
	//cout << " Np_1 " << mellinint::Np_1[i1] << " Np_2 " <<  mellinint::Np_2[i2]
	//	   << " CCp " << pow(mellinint::CCp/M_PI,2) << " exp(sumnpp) " << exp(-sumnpp*ax/2.) << " deltay " << ymx-ymn
	//	   << " diffnpp " << diffnpp << " expymx " << exp(-diffnpp*ymx) << " expymn " << exp(-diffnpp*ymn) << endl;
	      
	      
	      complex <double> sumnpm = mellinint::Np_1[i1]+mellinint::Nm_2[i2];
	      complex <double> diffnpm = mellinint::Np_1[i1]-mellinint::Nm_2[i2];
	      //if (i1 == 0 && i2 == 0 && opts.mellininv == 1 && !mellinint::midpoint)
	      if (mellinint::Np_1[i1] == mellinint::Nm_2[i2]) // this never happens --> It happens for Talbot!
		yintm=(mellinint::CCp/M_PI)*(mellinint::CCm/M_PI)*exp(-sumnpm*ax/2.)*(ymx-ymn);
	      else
		//yintm=1./(-diffnpm)*(mellinint::CCp/M_PI)*(mellinint::CCm/M_PI)*exp(-sumnpm*ax/2.)*(exp(-diffnpm*ymx)-exp(-diffnpm*ymn));
		//yintm=1./(-diffnpm)*(mellinint::CCp/M_PI)*(mellinint::CCm/M_PI)*(exp(-diffnpm*ymx-sumnpm*ax/2.)-exp(-diffnpm*ymn-sumnpm*ax/2.));
		yintm=1./(-diffnpm)*(mellinint::CCp/M_PI)*(mellinint::CCm/M_PI)
		  *(exp(-mellinint::Np_1[i1]*(ax/2.+ymx))*exp(-mellinint::Nm_2[i2]*(ax/2.-ymx))
		    -exp(-mellinint::Np_1[i1]*(ax/2.+ymn))*exp(-mellinint::Nm_2[i2]*(ax/2.-ymn)));

//	      if (phasespace::ymax >= ylim)  //when ymx = ylim, then x1_max = 1 and (ax/2.+ymx) = log(x1_max) = 0
//		yintm=1./(-diffnpm)*(mellinint::CCp/M_PI)*(mellinint::CCm/M_PI)
//		  *(exp(-mellinint::Nm_2[i2]*(ax/2.-ymx))
//		    -exp(-mellinint::Np_1[i1]*(ax/2.+ymn))*exp(-mellinint::Nm_2[i2]*(ax/2.-ymn)));
	      
	      double cthmom0, cthmom1, cthmom2;
	      omegaintegr::cthmoments(cthmom0,cthmom1,cthmom2);
	      if (helicity >= 0)
		{
		  if (helicity == 0)
		    {
		      cthmom0 = 4./3.; //20./3. * (0.5*2. - 1.5*2./3.  ) +2./3.*2.;
		      cthmom1 = 0.;
		      cthmom2 = -4./3.; //20./3. * (0.5*2./3.-1.5*2./5.  ) +2./3.*2./3.;
		    }
		  else if (helicity == 4)
		    {
		      cthmom0=0.;
		      //cthmom1=4*cthmom2*(y < 0 ? -1. : 1.);
		      cthmom1=4*cthmom2;
		      cthmom2=0.;
		    }
		  else
		    {
		      cthmom0=0.;
		      cthmom1=0.;
		      cthmom2=0.;
		    }
		}
	      
	      Ith0p[mellinint::index(i1,i2)] += cthmom0*yintp*mellinint::wn_1[i1]*mellinint::wn_2[i2];
	      Ith1p[mellinint::index(i1,i2)] += cthmom1*yintp*mellinint::wn_1[i1]*mellinint::wn_2[i2];
	      Ith2p[mellinint::index(i1,i2)] += cthmom2*yintp*mellinint::wn_1[i1]*mellinint::wn_2[i2];
	      Ith0m[mellinint::index(i1,i2)] += cthmom0*yintm*mellinint::wn_1[i1]*conj(mellinint::wn_2[i2]);
	      Ith1m[mellinint::index(i1,i2)] += cthmom1*yintm*mellinint::wn_1[i1]*conj(mellinint::wn_2[i2]);
	      Ith2m[mellinint::index(i1,i2)] += cthmom2*yintm*mellinint::wn_1[i1]*conj(mellinint::wn_2[i2]);
	      //cout << "cthmom0 " << cthmom0 << " yintp " << yintp << " wn_1 " << mellinint::wn_1[i1] << " wn_2 " << mellinint::wn_2[i2] << endl;
	    }
	}
  else //Numerical integration
    {
      //cache the mass dependent part (ax) of the exponential
      //The caching is done for the numerical integration only, but can be done also for the analytical integration (there is no gain in speed, but the code would be cleaner)
      complex <double> cfpm[mellinint::mdim][mellinint::mdim];
      complex <double> cfmm[mellinint::mdim][mellinint::mdim];
      for (int i1 = 0; i1 < mellinint::mdim; i1++)
	for (int i2 = 0; i2 < mellinint::mdim; i2++)
	  {
	    //Reduce by a factor of 2 the number of exponentials (only upper diagonal triangle) the positive piece is symmetric under i1 <-> i2		  
	    //In the negative piece exchange of i1 <-> i2 needs a complex conjugate
	    if (i2 >= i1)
	      {
		complex <double> sumnpp = mellinint::Np_1[i1]+mellinint::Np_2[i2];
		complex <double> sumnpm = mellinint::Np_1[i1]+mellinint::Nm_2[i2];
		cfpm[i1][i2]=pow(mellinint::CCp/M_PI,2)*exp(-sumnpp*ax/2.);                
		cfmm[i1][i2]=mellinint::CCp*mellinint::CCm/pow(M_PI,2)*exp(-sumnpm*ax/2.);
	      }
	    else
	      {
		cfpm[i1][i2]=cfpm[i2][i1];
		cfmm[i1][i2]=conj(cfmm[i2][i1]);
	      }
	  }

      // start integration
      double xc=0.5*(ymin+ymax);
      double xm=0.5*(ymax-ymin);
      for (int i=0; i < opts.yintervals; i++)
	{
	  double ya = ymin+(ymax-ymin)*i/opts.yintervals;
	  double yb = ymin+(ymax-ymin)*(i+1)/opts.yintervals;
	  xc=0.5*(ya+yb);
	  xm=0.5*(yb-ya);
	  for (int j=0; j < opts.yrule; j++)
	    {
	      double y=xc+xm*gr::xxx[opts.yrule-1][j];

	      //calculate costheta moments as a function of y
	      phasespace::set_y(y);
	      omegaintegr::genV4p();
	      double cthmom0, cthmom1, cthmom2;
	      omegaintegr::cthmoments(cthmom0,cthmom1,cthmom2);
	      //cout << cthmom0 << "  " << cthmom1 << "  " << cthmom2 << endl;
	      
	      for (int i1 = 0; i1 < mellinint::mdim; i1++)
		for (int i2 = 0; i2 < mellinint::mdim; i2++)
		  {
		    //if (i2 >= i1) //Be careful!!! can symmetrise this only for lepton charge-symmetric cuts
		    // {
		    //functions to integrate (rapidity part is cached at initialisation)
		    complex <double> fpy=cfpm[i1][i2]*cfpy[index(i, j, i1, i2)];
		    complex <double> fmy=cfmm[i1][i2]*cfmy[index(i, j, i1, i2)];
		    //integrals
		    Ith0p[mellinint::index(i1,i2)]+=fpy*cthmom0;
		    Ith1p[mellinint::index(i1,i2)]+=fpy*cthmom1;
		    Ith2p[mellinint::index(i1,i2)]+=fpy*cthmom2;
		    Ith0m[mellinint::index(i1,i2)]+=fmy*cthmom0;
		    Ith1m[mellinint::index(i1,i2)]+=fmy*cthmom1;
		    Ith2m[mellinint::index(i1,i2)]+=fmy*cthmom2;
		    //  }
		    //else
		    //  {
		    // 	 Ith0p[mellinint::index(i1,i2)]=Ith0p[mellinint::index(i2,i1)];
		    //	 Ith1p[mellinint::index(i1,i2)]=Ith1p[mellinint::index(i2,i1)];
		    //	 Ith2p[mellinint::index(i1,i2)]=Ith2p[mellinint::index(i2,i1)];
		    //	 Ith0m[mellinint::index(i1,i2)]=conj(Ith0m[mellinint::index(i2,i1)]);
		    //	 Ith1m[mellinint::index(i1,i2)]=conj(Ith1m[mellinint::index(i2,i1)]);
		    //	 Ith2m[mellinint::index(i1,i2)]=conj(Ith2m[mellinint::index(i2,i1)]);
		    //  }
		  }
	    }
	}
      //      for (int i1 = 0; i1 < 3; i1++)
      //	for (int i2 = 0; i2 < 3; i2++)
      //	  cout << i1 << "  " << i2 << "  " << Ith1p[mellinint::index(i1,i2)] << "  " << cfpm[i1][i2] << endl;
    }
//  for (int i1 = 0; i1 < 3; i1++)
//    for (int i2 = 0; i2 < 3; i2++)
//      cout << i1 << "  " << i2 << "  " << Ith0p[mellinint::index(i1,i2)] << "  " << Ith0m[mellinint::index(i1,i2)] << endl;
}

//Numerical integration without any caching
void rapint::numint(double ymin, double ymax, double m, int helicity)
{
  //integration results are stored in Ithxx

  fill(Ith0p,Ith0p+mellinint::mdim*mellinint::mdim,0.);
  fill(Ith1p,Ith1p+mellinint::mdim*mellinint::mdim,0.);
  fill(Ith2p,Ith2p+mellinint::mdim*mellinint::mdim,0.);
  fill(Ith0m,Ith0m+mellinint::mdim*mellinint::mdim,0.);
  fill(Ith1m,Ith1m+mellinint::mdim*mellinint::mdim,0.);
  fill(Ith2m,Ith2m+mellinint::mdim*mellinint::mdim,0.);

  double q2 = m*m;
  double bjx= q2/pow(opts.sroot,2);
  double ax = log(bjx);
  //double ax1 = (ax+2*y)/2.; 
  //double ax2 = (ax-2*y)/2.; 

  //cache the mass dependent part (ax) of the exponential
  //The caching is done for the numerical integration only, but can be done also for the analytical integration (there is no gain in speed, but the code would be cleaner)
  complex <double> cfpm[mellinint::mdim][mellinint::mdim];
  complex <double> cfmm[mellinint::mdim][mellinint::mdim];
  for (int i1 = 0; i1 < mellinint::mdim; i1++)
    for (int i2 = 0; i2 < mellinint::mdim; i2++)
      {
	complex <double> sumnpp = mellinint::Np_1[i1]+mellinint::Np_2[i2];
	complex <double> sumnpm = mellinint::Np_1[i1]+mellinint::Nm_2[i2];
	cfpm[i1][i2]=pow(mellinint::CCp/M_PI,2)*exp(-sumnpp*ax/2.);                
	cfmm[i1][i2]=mellinint::CCp*mellinint::CCm/pow(M_PI,2)*exp(-sumnpm*ax/2.);
      }

  // start integration
  double xc=0.5*(ymin+ymax);
  double xm=0.5*(ymax-ymin);
  for (int i=0; i < opts.yintervals; i++)
    {
      double ya = ymin+(ymax-ymin)*i/opts.yintervals;
      double yb = ymin+(ymax-ymin)*(i+1)/opts.yintervals;
      xc=0.5*(ya+yb);
      xm=0.5*(yb-ya);
      for (int j=0; j < opts.yrule; j++)
	{
	  double y=xc+xm*gr::xxx[opts.yrule-1][j];
	  //<< endl;
	  //calculate costheta moments as a function of y
	  phasespace::set_y(y);
	  omegaintegr::genV4p();
	  double cthmom0, cthmom1, cthmom2;
	  omegaintegr::cthmoments(cthmom0,cthmom1,cthmom2);
	  //helicity option is to be supported only without cuts on the leptons
	  if (!opts.makecuts && helicity >= 0)
	    {
	      if (helicity == 0)
		{
		  cthmom0 = 4./3.; //20./3. * (0.5*2. - 1.5*2./3.  ) +2./3.*2.;
		  cthmom1 = 0.;
		  cthmom2 = -4./3.; //20./3. * (0.5*2./3.-1.5*2./5.  ) +2./3.*2./3.;
		}
	      else if (helicity == 4)
		{
		  cthmom0=0.;
		  cthmom1=4*cthmom2*(y < 0 ? -1. : 1.);
		  cthmom2=0.;
		}
	      else
		{
		  cthmom0=0.;
		  cthmom1=0.;
		  cthmom2=0.;
		}
	    }

	  //mellinint::update();
	  //pdfevol::update();

	  double z1 = ax/2.+y; 
	  double z2 = ax/2.-y; 
	  for (int i1 = 0; i1 < mellinint::mdim; i1++)
	    for (int i2 = 0; i2 < mellinint::mdim; i2++)
	      {
		//functions to integrate (rapidity part is cached at initialisation)
		//complex <double> diffnpp = mellinint::Np_1[i1]-mellinint::Np_2[i2];
		//complex <double> diffnpm = mellinint::Np_1[i1]-mellinint::Nm_2[i2];
		//complex <double> cfpy=exp(-diffnpp*y)*mellinint::wn_1[i1]*mellinint::wn_2[i2]*gr::www[opts.yrule-1][j]*xm;
		//complex <double> cfmy=exp(-diffnpm*y)*mellinint::wn_1[i1]*conj(mellinint::wn_2[i2])*gr::www[opts.yrule-1][j]*xm;
		//complex <double> fp=cfpm[i1][i2]*cfpy;
		//complex <double> fm=cfmm[i1][i2]*cfmy;

		complex <double> sumnpp = mellinint::Np_1[i1]+mellinint::Np_2[i2];
		complex <double> sumnpm = mellinint::Np_1[i1]+mellinint::Nm_2[i2];
		complex <double> diffnpp = mellinint::Np_1[i1]-mellinint::Np_2[i2];
		complex <double> diffnpm = mellinint::Np_1[i1]-mellinint::Nm_2[i2];
		complex <double> fp = pow(mellinint::CCp/M_PI,2)
		  //* exp(-diffnpp*y-sumnpp*ax/2.)*mellinint::wn_1[i1]*mellinint::wn_2[i2]*gr::www[opts.yrule-1][j]*xm;
		  * exp(-mellinint::Np_1[i1]*z1 -mellinint::Np_2[i2]*z2)*mellinint::wn_1[i1]*mellinint::wn_2[i2]*gr::www[opts.yrule-1][j]*xm;
		complex <double> fm = mellinint::CCp*mellinint::CCm/pow(M_PI,2)
		  //* exp(-diffnpm*y-sumnpm*ax/2.)*mellinint::wn_1[i1]*conj(mellinint::wn_2[i2])*gr::www[opts.yrule-1][j]*xm;
		  * exp(-mellinint::Np_1[i1]*z1 -mellinint::Nm_2[i2]*z2)*mellinint::wn_1[i1]*conj(mellinint::wn_2[i2])*gr::www[opts.yrule-1][j]*xm;
	
		//integrals
		Ith0p[mellinint::index(i1,i2)]+=fp*cthmom0;
		Ith1p[mellinint::index(i1,i2)]+=fp*cthmom1;
		Ith2p[mellinint::index(i1,i2)]+=fp*cthmom2;
		Ith0m[mellinint::index(i1,i2)]+=fm*cthmom0;
		Ith1m[mellinint::index(i1,i2)]+=fm*cthmom1;
		Ith2m[mellinint::index(i1,i2)]+=fm*cthmom2;
	      }
	}
    }
//  for (int i1 = 0; i1 < 3; i1++)
//    for (int i2 = 0; i2 < 3; i2++)
//      cout << i1 << "  " << i2 << "  " << Ith0p[mellinint::index(i1,i2)] << "  " << Ith0m[mellinint::index(i1,i2)] << endl;
//
//  cout << endl;
}

void rapint::free()
{
  //integration results
  delete[] Ith0p;
  delete[] Ith1p;
  delete[] Ith2p;
  delete[] Ith0m;
  delete[] Ith1m;
  delete[] Ith2m;
}
