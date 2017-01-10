#include "vjint.h"
#include "phasespace.h"
#include "settings.h"
#include "coupling.h"
#include "isnan.h"
#include "gaussrules.h"
#include "numbers.h"
#include "luminosity.h"
#include "mesq.h"
#include "resconst.h"

#include <iostream>
#include <iomanip>

double vjint::sing()
{
  double intsing = 0.;
  //include 'fodyqt_inc.f'
  //include 'internal_inc.f'
  //include 'luminosities_inc.f'
  //include 'scales2_inc.f'

  //....definition of quantitites for phase space integration
  //....as defined in formula B.1 of Glosser-Schmidt paper (arXiv:hep-ph/0209248)
  // zz1 <-> z1
  // zz2 <-> z2
  // to improve numerical accuracy, the change of variable z1 -> 1-z1 and z2 -> 1-z2 is applied

  //.....x1_0,x2_0 and modification by Massimiliano (dcut)
  double x10 = phasespace::mt/opts.sroot*phasespace::exppy;
  double x20 = phasespace::mt/opts.sroot*phasespace::expmy;
  double dcut = x10*pow(phasespace::qt/phasespace::mt,2)/(1.-x10*(1.-pow(phasespace::qt/phasespace::mt,2)));

  double q2 = phasespace::m2;
  
  double fac = gevfb/1000.*coupling::aemmz*asp_.asp_*resconst::Cf;
  double xjacz1 = (1.-dcut-x20);

  double sh, uh, th, s2;
  double x1,x2;
  
  //.....z1,z2,lambda_b and corresponding jacobian 
  //.....(to allow integration from 0 to 1)

  //gaussian quadrature
  
  //start zz1 integration
  double az1 = zmin;
  double bz1 = zmax;
  double cz1=0.5*(az1+bz1);
  double mz1=0.5*(bz1-az1);
  for (int j = 0; j < z1rule; j++)
    {
      double zz1 = t1[j];//zmin*pow(zmax/zmin,z1);
      double jacz1=zz1*lz;

      double z1 = dcut+(1.-dcut-x20)*zz1;
      double lb = phasespace::qt2/phasespace::mt2*(1.-z1)/z1;

      //     *********************************************************
      //     This piece depends only on z1, phase space is the same, with z2=1 (which leads to s2=0)
      //     It corresponds to fs(0) in Eq. (2.5) of [Gonsalves, Pawlowsky, Wai]
  
      //.....subtraction of z2=1 term according to plus prescription     
      x1 = x10*(1.+lb);
      x2 = x20/(1.-z1);
      luminosity::pdf1(x1);
      luminosity::pdf2(x2);
      
      double tiny = 0.;//1e-7;
      if(x1 > 1.-tiny || x2 > 1.-tiny)
	return 0.;

      //.....definition of partonic Mandelstam invariants -> bug fix in DYqT: th <-> uh
      double sh = phasespace::mt2/(1.-z1)*(1.+lb);
      double uh = -phasespace::mt2/(1.-z1)*z1*(1.+lb);
      double th = q2-phasespace::mt2*(1.+lb);
      double s2 = 0.; // set s2=0 explicitly in utilities
      //cout << sh << " " << x1*x2*s << " " << th << " " << q2-ss*x1*tm*expym << " " << uh,q2-ss*x2*tm*expyp << endl;

      //.....phase space prefactor
      double pre10 = phasespace::mt2*(1.+lb)/pow(1.-z1,2);

      //.....compute parton luminosity
      //fractions_.x1_ = x1;
      //fractions_.x2_ = x2;
      //flavour_();
      luminosity::calc();
      utils_fu_(uh,q2);
      utils_(sh,th,uh,q2,s2);
         
      //.....common factor for all the contributions
      double factor = fac/sh;

      //.....(log(s2)/s2)_A+ contributions
      double xsqg10=(factor/(pow(coupling::NC,2)-1.))*
	(cqg3_(sh,th,uh,q2,one)+cgq3_(sh,th,uh,q2,one));
      double xsqqb10=(factor/coupling::NC)*
	(cqqb2_(sh,th,uh,q2,one)+cqqb2x_(sh,th,uh,q2,one));
      double xs10=xsqg10+xsqqb10;

      //.....(1/s2)_A+ contributions
      double xsqg20=(factor/(pow(coupling::NC,2)-1.))*
	(cqg3_(sh,th,uh,q2,two)+cgq3_(sh,th,uh,q2,two));
      double xsqqb20=(factor/coupling::NC)*
	(cqqb2_(sh,th,uh,q2,two)+cqqb2x_(sh,th,uh,q2,two)+
	 daa_(sh,th,uh,q2,two)+daax_(sh,th,uh,q2,two));
      double xs20=xsqg20+xsqqb20;
      //     *********************************************************


//.....final result.................................................
//.....note: here I am passing from s2 to z2 distributions..........
//.....according to the formulas:...................................
//.....(1/s2)_A+ = (1/-t)*((1/(1-z2))_+ -1 
//.....                   + log(z2min/(1-z2min))*delta(1-z2))
//.....(log(s2/q2)/s2)_A+ = log(-t/q2)*(1/s2)_A+ - 1/t*(
//.....                   (log(1-z2)/(1-z2))_+ + log(z2/(1-z2))
//.....                   - log(z2)/(1-z2) - 1/2*delta(1-z2)*
//.....                   log^2((1-z2min)/z2min))
//.....as a consequence, there appear some "spurious" terms.........
//.....in nlo delta, 1/s2 and regular contributions.................     
//.....I do not include spurious delta terms: they cancel against...
//.....similar contributions in subroutine 'xdelta'.................
//.....provided that A=-th in Gonsalves' formulas...................
//.....I also add mismatch terms (xmism) coming from 0-z2min region.

//     bug fix in DYqT: th <-> uh
      double a = -uh;
      double logaq2 = log(a/q2);
      double z2min = x10*(1.+lb);
      double xmism = -0.5*xs10*pow(log(1.-z2min),2)-(xs10*logaq2+xs20)*log(1.-z2min);

      xmism *= pre10*xjacz1/a;
      
      //start zz2 integration
      double az2 = zmin;
      double bz2 = zmax;
      double cz2 = 0.5*(az2+bz2);
      double mz2 = 0.5*(bz2-az2);
      for (int jj = 0; jj < z2rule; jj++)
	{
	  double zz2 = t2[jj];//zmin*pow(zmax/zmin,z2);
	  double jacz2 = zz2*lz;
	      

	  double z2 = (1.-x10*(1.+lb))*zz2;
	  double xjacz2 = (1.-x10*(1.+lb));

	  double xjac = xjacz1*xjacz2;
	  
	  //.....first term of the integrand for the plus prescription
	  x1 = x10/(1.-z2)*(1.+lb);
	  luminosity::pdf1(x1);

	  //.....definition of partonic Mandelstam invariants --> bug fix in DYqT: th <-> uh
	  double sh = phasespace::mt2/(1.-z1)/(1.-z2)*(1.+lb);
	  double th = q2-phasespace::mt2/(1.-z2)*(1.+lb);
	  double s2 = phasespace::mt2/(1.-z1)/(1.-z2)*z2*z1*(1.+lb);
	  //cout << "sh " << sh << " " << x1*x2*pow(opts.sroot,2) << endl;
	  //cout << "th " << th << " " << q2-opts.sroot*x1*phasespace::mt*phasespace::expmy << endl;
	  //cout << "uh " << uh << " " << q2-opts.sroot*x2*phasespace::mt*phasespace::exppy << endl;
	  //cout << "s2 " << s2 << " " << sh+uh+th-q2 << endl;

	  //check s2 > 0
	  if (s2 == 0.) 
	    return 0.;
	  if (s2 < 0)
	    {
	      cout << " s2 < 0 ! s2 = " << s2 << endl;
	      cout << "z1 " << z1 << " z2 " << z2 << " x1 " << x1 << " x2 " << x2 << " sh " << sh << " th " << th << " uh " << uh << endl;
	      cout << " m " << phasespace::m  << " pt " << phasespace::qt << " y  " << phasespace::y << endl;
	      return 0.;
	    }
	
	  //.....phase space prefactor
	  double pre1 = phasespace::mt2*(1.+lb)/pow((1.-z1)*(1.-z2),2);

	  //check x1,x2 < 1
	  double tiny = 0.;//1e-7;
	  if(x1 > 1.-tiny || x2 > 1.-tiny)
	    {
	      cout << "x1 " << x1 << " x2 " << x2 << endl;
	      return 0.;
	    }

	  //.....compute parton luminosity
	  //fractions_.x1_ = x1;
	  //fractions_.x2_ = x2;
	  //flavour_();
	  luminosity::calc();	  
	  utils_(sh,th,uh,q2,s2);

	  //.....common factor for all the contributions
	  double factor = fac/sh;

	  //.....regular contributions (non-delta and non-singular in s2)
	  double xrgg=(factor*(coupling::NC/pow(pow(coupling::NC,2)-1.,2)))*
	    (cgg1_(sh,th,uh,q2)+cgg1x_(sh,th,uh,q2));
	  double xrqg=(factor/(pow(coupling::NC,2)-1.))*
	    (cqg3_(sh,th,uh,q2,three)+cgq3_(sh,th,uh,q2,three)+
	     cqg3_(sh,th,uh,q2,four)+cgq3_(sh,th,uh,q2,four));
	  double xrqqb=(factor/coupling::NC)*
	    (cqqb2_(sh,th,uh,q2,four)+cqqb2x_(sh,th,uh,q2,four)
	     +daa_(sh,th,uh,q2,four)+daax_(sh,th,uh,q2,four)
	     +dab_(sh,th,uh,q2)+dabx_(sh,th,uh,q2)
	     +dbb_(sh,th,uh,q2)+dbbx_(sh,th,uh,q2)
	     +dac_(sh,th,uh,q2)+dad_(sh,th,uh,q2)
	     +dbc_(sh,th,uh,q2)+dbd_(sh,th,uh,q2)
	     +dcc_(sh,th,uh,q2)+ddd_(sh,th,uh,q2)
	     +dcdll_(sh,th,uh,q2)+dcdllx_(sh,th,uh,q2)
	     +dcdlr_(sh,th,uh,q2)+dcdlrx_(sh,th,uh,q2));
	  double xrqq=(factor/coupling::NC)*(1/2.)*
	    (eac_(sh,th,uh,q2)+ebd_(sh,th,uh,q2)
	     +ead_(sh,th,uh,q2)+ebc_(sh,th,uh,q2)
	     // missing pieces in qq
	     //        + eaa_(sh,th,uh,q2)+ecc_(sh,th,uh,q2) // eaa = ecc = dcc -> eaa+ecc=2*dcc (see A28 of http://journals.aps.org/prd/pdf/10.1103/PhysRevD.40.2245)
	     //        + ebb_(sh,th,uh,q2)+edd_(sh,th,uh,q2) // ebb = edd = ddd -> ebb+edd=2*ddd (see A30 of http://journals.aps.org/prd/pdf/10.1103/PhysRevD.40.2245)
	     //        + eabll+eabllx //eabll=-dcdlr
	     //        + eablr+eablrx //
	     //        + ecdll+ecdllx //ecdll=-dcdlr
	     //        + ecdlr+ecdlrx //ecdlr=-dcdll
	     );
	  double xreg = xrgg+xrqg+xrqqb+xrqq;

	  if (isnan_ofast(xrqg))
	    {
	      cout << "nan in xrqg s2 = " << s2 << endl;
	      xrqg = 0.;
	    }

	  if (isnan_ofast(xrqqb))
	    {
	      cout << "nan in xrqqb, s2 = " << s2 << endl;
	      xrqqb = 0.;
	    }
	  //cout << xreg << "  " << factor<< " " << sh<< " " << th<< " " << uh<< " " << q2<< " " << s2 << endl;

	  //.....(log(s2)/s2)_A+ contributions
	  double xsqg1=(factor/(pow(coupling::NC,2)-1.))*
	    (cqg3_(sh,th,uh,q2,one)+cgq3_(sh,th,uh,q2,one));
	  double xsqqb1=(factor/coupling::NC)*
	    (cqqb2_(sh,th,uh,q2,one)+cqqb2x_(sh,th,uh,q2,one));
	  double xs1 = xsqg1+xsqqb1;


	  //.....(1/s2)_A+ contributions
	  double xsqg2=(factor/(pow(coupling::NC,2)-1.))*
	    (cqg3_(sh,th,uh,q2,two)+cgq3_(sh,th,uh,q2,two));
	  double xsqqb2=(factor/coupling::NC)*
	    (cqqb2_(sh,th,uh,q2,two)+cqqb2x_(sh,th,uh,q2,two)+
	     daa_(sh,th,uh,q2,two)+daax_(sh,th,uh,q2,two));
	  double xs2=xsqg2+xsqqb2;


	  double sing1 = log(z2)/(z2)*(pre1*xs1-pre10*xs10)*xjac/a;
	  double sing2 = (pre1*(xs1*logaq2+xs2)-pre10*(xs10*logaq2+xs20))/(z2)*xjac/a;
	  double reg = pre1*(xreg+(xs1*(log((1.-z2)/(z2))-(log(1.-z2))/(z2)-logaq2)-xs2)/(a))*xjac;
      
	  double sing=(sing1+sing2+reg-xmism)/pow(opts.sroot,2);
	  sing *= M_PI;

	  intsing += sing*gr::www[z1rule-1][j]*gr::www[z2rule-1][jj]
	    *jacz1*jacz2*mz1*mz2;
	  
	  if (isnan_ofast(sing))
	    {
	      cout << "nan in xsing " << sing1 << "  " << sing2 << "  " << reg << "  " << xmism << endl;
	      return 0.;
	    }
	}
    }
  return intsing;
}
