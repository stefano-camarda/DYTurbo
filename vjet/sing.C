#include "vjint.h"
#include "codata.h"
#include "phasespace.h"
#include "settings.h"
#include "coupling.h"
#include "isnan.h"
#include "gaussrules.h"
#include "numbers.h"
#include "luminosity.h"
#include "mesq.h"
#include "scales.h"
#include "resconst.h"
#include "LHAPDF/LHAPDF.h"

#include <iostream>
#include <iomanip>

double az = 0.;
double bz = 1.;
double cz = 0.5*(az+bz);
double mz = 0.5*(bz-az);

double vjint::sing()
{
  double intsing = 0.;

  //variables of phase space integration as defined in appendix B (formula B.1) of Glosser-Schmidt arXiv:hep-ph/0209248

  //x1_0,x2_0 and delta (B.2 of arXiv:hep-ph/0209248)
  double x10 = phasespace::mt/opts.sroot*phasespace::exppy;
  double x20 = phasespace::mt/opts.sroot*phasespace::expmy;
  //double delta = phasespace::qt/(phasespace::mt+phasespace::qt);                                              //delta as in arXiv:hep-ph/0209248
  //double delta = x10*phasespace::qt2/phasespace::mt2/(1.-x10*(1.-phasespace::qt2/phasespace::mt2));           //delta as in dyqt
  double delta = x10*phasespace::qt2/(x10*phasespace::qt2+(1-x10)*phasespace::mt2);                             //delta from the relation x1 < 1 <=> x10*(1.+lb) < 1

  double q2 = phasespace::m2;
  double fac = gevfb/1000.*coupling::aemmz*asp_.asp_*resconst::Cf;

  //double sh, uh, th, s2;
  //double x1,x2;

  //z1,z2,lambda_b as defined in arXiv:hep-ph/0209248
  //to improve numerical accuracy, the change of variable z1 -> 1-z1 and z2 -> 1-z2 is applied

  //start z1 integration with gaussian quadrature
  double z1min = delta;
  double z1max = (1.-x20);
  double lz1 = log(z1max/z1min);
  for (int j = 0; j < z1rule; j++)
    {
      //double z1 = delta+(1.-x20-delta)*t1[j];
      //double jacz1 = t1[j]*lz*(1.-x20-delta);

      double z1 = z1min*pow(z1max/z1min, cz+mz*gr::xxx[z1rule-1][j]);
      double jacz1 = z1*lz1;
      
      double lb = phasespace::qt2/phasespace::mt2*(1.-z1)/z1; //(B.4 of arXiv:hep-ph/0209248)

      //This piece depends only on z1, phase space is the same, with z2=1 (which leads to s2=0)
      //It corresponds to fs(0) in Eq. (2.5) of [Gonsalves, Pawlowsky, Wai]
  
      //subtraction of z2=1 term according to plus prescription     
      double x1 = x10*(1.+lb);
      double x2 = x20/(1.-z1);

      double tiny = 0.;//1e-7;
      if(x1 > 1.-tiny || x2 > 1.-tiny)
	{
	  //cout << "x1 " << x1 << " x2 " << x2 << endl;
	  continue;
	}

      //partonic Mandelstam invariants -> bug fix in DYqT: th <-> uh
      double sh = phasespace::mt2/(1.-z1)*(1.+lb);
      double uh = -phasespace::mt2/(1.-z1)*z1*(1.+lb);
      double th = q2-phasespace::mt2*(1.+lb);
      double s2 = 0.; // set s2=0 explicitly in utilities
      //cout << sh << " " << x1*x2*s << " " << th << " " << q2-ss*x1*tm*expym << " " << uh,q2-ss*x2*tm*expyp << endl;

      if (opts.fmuren > 2 || opts.fmufac > 2)
	{
	  scales::set(phasespace::m, phasespace::qt);
	  scales::vjet();
	  fac = gevfb/1000.*coupling::aemmz*asp_.asp_*resconst::Cf;
	  utils_scales_(q2);
	}

      /*      
      if (opts.kmjj_muren != 0 || opts.kmjj_mufac != 0)
	{
	  scales2_.xmur_ = sqrt(pow(opts.kmuren*opts.rmass,2) + pow(opts.kpt_muren*phasespace::qt,2) + pow(opts.kmjj_muren*sqrt(s2),2));
	  scales2_.xmuf_ = sqrt(pow(opts.kmufac*opts.rmass,2) + pow(opts.kpt_mufac*phasespace::qt,2) + pow(opts.kmjj_mufac*sqrt(s2),2));
	  scales2_.xmur2_ = pow(scales2_.xmur_,2);
	  scales2_.xmuf2_ = pow(scales2_.xmuf_,2);
	  asnew_.as_ = LHAPDF::alphasPDF(scales2_.xmur_)/M_PI;
	  asp_.asp_ = asnew_.as_*M_PI;
	  fac = gevfb/1000.*coupling::aemmz*asp_.asp_*resconst::Cf;
	  utils_scales_(q2);
	}
      */

      luminosity::pdf1(x1);
      luminosity::pdf2(x2);
      
      //phase space prefactor
      double pre10 = phasespace::mt2*(1.+lb)/pow(1.-z1,2);

      //parton luminosities
      luminosity::calc();
      utils_fu_(uh,q2);
      utils_(sh,th,uh,q2,s2);
         
      //common factor for all the contributions
      double factor = fac/sh;

      //(log(s2)/s2)_A+ contributions
      double xsqg10=(factor/(pow(coupling::NC,2)-1.))*
	(cqg3_(sh,th,uh,q2,one)+cgq3_(sh,th,uh,q2,one));
      double xsqqb10=(factor/coupling::NC)*
	(cqqb2_(sh,th,uh,q2,one)+cqqb2x_(sh,th,uh,q2,one));
      double xs10=xsqg10+xsqqb10;

      //(1/s2)_A+ contributions
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

      //bug fix in DYqT: th <-> uh
      double a = -uh;
      double logaq2 = log(a/q2);
      double z2min = x10*(1.+lb);
      double xmism = -0.5*xs10*pow(log(1.-z2min),2)-(xs10*logaq2+xs20)*log(1.-z2min);

      xmism *= pre10/a*asp_.asp_;
      
      //start z2 integration
      for (int jj = 0; jj < z2rule; jj++)
	{
	  double z2 = (1.-z2min)*t2[jj];
	  double jacz2 = (1.-z2min)*t2[jj]*lz;
	  
	  //double z2mn = 1e-13;
	  //double z2mx = 1.-z2min;
	  //double z2 = z2mn*pow(z2mx/z2mn, cz+mz*gr::xxx[z2rule-1][jj]);
	  //double jacz2 = z2*log(z2mx/z2mn);

	  //first term of the integrand for the plus prescription
	  x1 = x10/(1.-z2)*(1.+lb);

	  //check x1 > 1
	  double tiny = 0.;//1e-7;
	  if(x1 > 1.-tiny)
	    {
	      //cout << "x1 " << x1 << " x2 " << x2 << endl;
	      continue;
	    }

	  //partonic Mandelstam invariants --> bug fix in DYqT: th <-> uh
	  sh = phasespace::mt2/(1.-z1)/(1.-z2)*(1.+lb);
	  th = q2-phasespace::mt2/(1.-z2)*(1.+lb);
	  s2 = phasespace::mt2/(1.-z1)/(1.-z2)*z2*z1*(1.+lb);
	  //cout << "sh " << sh << " " << x1*x2*pow(opts.sroot,2) << endl;
	  //cout << "th " << th << " " << q2-opts.sroot*x1*phasespace::mt*phasespace::expmy << endl;
	  //cout << "uh " << uh << " " << q2-opts.sroot*x2*phasespace::mt*phasespace::exppy << endl;
	  //cout << "s2 " << s2 << " " << sh+uh+th-q2 << endl;

	  //check s2 > 0
	  if (s2 == 0.)
	    {
	      //cout << " s2 == 0 ! s2 = " << s2 << endl;
	      //continue;
	      return 0.;
	    }
	  if (s2 < 0)
	    {
	      cout << " s2 < 0 ! s2 = " << s2 << endl;
	      cout << "z1 " << z1 << " z2 " << z2 << " x1 " << x1 << " x2 " << x2 << " sh " << sh << " th " << th << " uh " << uh << endl;
	      cout << " m " << phasespace::m  << " pt " << phasespace::qt << " y  " << phasespace::y << endl;
	      return 0.;
	    }

	  if (opts.fmuren > 2 || opts.fmufac > 2)
	    {
	      scales::set(phasespace::m, phasespace::qt, sqrt(s2));
	      scales::vjet();
	      fac = gevfb/1000.*coupling::aemmz*asp_.asp_*resconst::Cf;
	      utils_scales_(q2);
	    }
	  /*
	  if (opts.kmjj_muren != 0 || opts.kmjj_mufac != 0)
	    {
	      scales2_.xmur_ = sqrt(pow(opts.kmuren*opts.rmass,2) + pow(opts.kpt_muren*phasespace::qt,2) + pow(opts.kmjj_muren*sqrt(s2),2));
	      scales2_.xmuf_ = sqrt(pow(opts.kmufac*opts.rmass,2) + pow(opts.kpt_mufac*phasespace::qt,2) + pow(opts.kmjj_mufac*sqrt(s2),2));
	      scales2_.xmur2_ = pow(scales2_.xmur_,2);
	      scales2_.xmuf2_ = pow(scales2_.xmuf_,2);
	      asnew_.as_ = LHAPDF::alphasPDF(scales2_.xmur_)/M_PI;
	      asp_.asp_ = asnew_.as_*M_PI;
	      fac = gevfb/1000.*coupling::aemmz*asp_.asp_*resconst::Cf;
	      utils_scales_(q2);
	      //printf ("s2 %f mur %f alphas %f\n", s2, scales2_.xmur_, asnew_.as_);
	    }
	  */
	  luminosity::pdf1(x1);

	  //phase space prefactor
	  double pre1 = phasespace::mt2*(1.+lb)/pow((1.-z1)*(1.-z2),2);

	  //compute parton luminosities
	  luminosity::calc();	  
	  utils_(sh,th,uh,q2,s2);

	  //common factor for all the contributions
	  double factor = fac/sh;

	  //regular contributions (non-delta and non-singular in s2)
	  double xrgg=(factor*(coupling::NC/pow(pow(coupling::NC,2)-1.,2)))*
	    (cgg1_(sh,th,uh,q2)+cgg1x_(sh,th,uh,q2));
	  double xrqg=(factor/(pow(coupling::NC,2)-1.))*
	    (cqg3_(sh,th,uh,q2,three)+cgq3_(sh,th,uh,q2,three)
	     +cqg3_(sh,th,uh,q2,four)+cgq3_(sh,th,uh,q2,four));
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

	  //(log(s2)/s2)_A+ contributions
	  double xsqg1=(factor/(pow(coupling::NC,2)-1.))*
	    (cqg3_(sh,th,uh,q2,one)+cgq3_(sh,th,uh,q2,one));
	  double xsqqb1=(factor/coupling::NC)*
	    (cqqb2_(sh,th,uh,q2,one)+cqqb2x_(sh,th,uh,q2,one));
	  double xs1 = xsqg1+xsqqb1;


	  //(1/s2)_A+ contributions
	  double xsqg2=(factor/(pow(coupling::NC,2)-1.))*
	    (cqg3_(sh,th,uh,q2,two)+cgq3_(sh,th,uh,q2,two));
	  double xsqqb2=(factor/coupling::NC)*
	    (cqqb2_(sh,th,uh,q2,two)+cqqb2x_(sh,th,uh,q2,two)+
	     daa_(sh,th,uh,q2,two)+daax_(sh,th,uh,q2,two));
	  double xs2=xsqg2+xsqqb2;

	  double sing1 = log(z2)/(z2)*(pre1*xs1-pre10*xs10)/a;
	  double sing2 = (pre1*(xs1*logaq2+xs2)-pre10*(xs10*logaq2+xs20))/(z2)/a;
	  double reg = pre1*(xreg+(xs1*(log((1.-z2)/(z2))-(log(1.-z2))/(z2)-logaq2)-xs2)/(a));
      
	  double sing = (sing1+sing2+reg)/pow(opts.sroot,2)*asp_.asp_;

	  intsing += sing*gr::www[z1rule-1][j]*gr::www[z2rule-1][jj]
	    *jacz1*jacz2*mz*mz;
	  
	  if (isnan_ofast(sing))
	    {
	      cout << "nan in xsing " << sing1 << "  " << sing2 << "  " << reg << "  " << xmism << endl;
	      return 0.;
	    }
	}
      intsing += (-xmism/pow(opts.sroot,2))*gr::www[z1rule-1][j]*jacz1*mz;
    }
  return intsing*M_PI;
}
