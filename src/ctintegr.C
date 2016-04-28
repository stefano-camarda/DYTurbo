#include "interface.h"
#include "ctintegr.h"
#include "integr.h"
#include "settings.h"
#include "switch.h"

#include <iostream>
#include <iomanip>
#include <ctime>
#include <math.h>


using namespace std;


integrand_t ctintegrand(const int &ndim, const double x[], const int &ncomp, double f[],
                        void* userdata, const int &nvec, const int &core,
                        double &weight, const int &iter)
{
  //here generate the phase space according to x[], and pass the p vector to countint_

  if (opts.PDFerrors)
    for (int i = 1; i < opts.totpdf; i++)
      f[i] = 0.;
    
  double rct[22];
  for (int i = 0; i < ndim; i++)
    rct[i]=x[i];
  rct[8] = rct[1];

  dofill_.doFill_ = int(iter==last_iter);
  f[0] = countint_(rct,weight);

  return 0;
}

integrand_t ctintegrandMC(const int &ndim, const double x[], const int &ncomp, double f[],
                        void* userdata, const int &nvec, const int &core,
                        double &weight, const int &iter)
{
  clock_t begin_time, end_time;

  begin_time = clock();

  if (opts.PDFerrors)
    for (int i = 1; i < opts.totpdf; i++)
      f[i] = 0.;
    
  //Jacobian of the change of variables from the unitary hypercube x[6] to the m, y, qt, phi, costh, philep boundaries
  double jac = 1.;

  // Generate the boson invariant mass between the integration boundaries
  double mcut = qtmax/qtcut_.xqtcut_;
  double wsqmin = pow(mmin,2);
  double wsqmax = pow(min(mmax,mcut),2);
  if (wsqmin >= wsqmax)
    {
      f[0]=0.;
      return 0;
    }
  double x1=x[0];
  double m2,wt;
  breitw_(x1,wsqmin,wsqmax,opts.rmass,opts.rwidth,m2,wt);
  double m=sqrt(m2);
  jac=jac*wt;

  //     Dynamic scale
  //      if(dynamicscale) call scaleset(m2)

  //Limit y boundaries to the kinematic limit in y
  double ylim = 0.5*log(pow(energy_.sroot_,2)/m2);
  double ymn = min(max(-ylim, ymin),ylim);
  double ymx = max(min(ylim, ymax),-ylim);
  if (ymn >= ymx)
    {
      f[0]=0.;
      return 0;
    }

  //integrate between ymin and ymax
  double y=ymn+(ymx-ymn)*x[1];
  jac=jac*(ymx-ymn);
  
  //integrate between qtmin and qtmax
  //use qt2 to get the correct jacobian!

  //In the fixed order calculation, integrate from qtcut to infinity
  double expy, expmy;
  double qtmn, qtmx;
  double qtcut = qtcut_.xqtcut_*m;
  if (opts.fixedorder)
    {
      if (qtmin > 0)
	{
	  f[0]=0.;
	  return 0;
	}
      qtmn = qtcut;
      qtmx = 1e10;
    }
  else
    {
      expy = exp(y);
      expmy = exp(-y);
      qtmn = max(qtcut, qtmin);
      double cosh2y34=pow((expy+expmy)*0.5,2);
      double kinqtlim = sqrt(pow(pow(energy_.sroot_,2)+m*m,2)/(4*pow(energy_.sroot_,2)*cosh2y34)-m*m);
      double switchqtlim = switching::qtlimit(m);
      double qtlim = min(kinqtlim, switchqtlim);
      qtmx = min(qtlim, qtmax);
      if (qtmn >= qtmx)
	{
	  f[0]=0.;
	  return 0;
	}
    }
  double qtmn2 = pow(qtmn,2);
  double qtmx2 = pow(qtmx,2);
  
  double tiny = 1E-5;
  double a = 1./(1+log(qtmx2/tiny));
  double b = 1./(1+log(qtmn2/tiny));
  double x2 = a + (b-a) * x[2];
  jac = jac * (b-a);
  double qt2=tiny*exp(1./x2 - 1);
  jac=jac*qt2/pow(x2,2);
  
  double qt=sqrt(qt2);

  // incoming quarks (are they actually used?)
  double xx1=sqrt(m2/pow(energy_.sroot_,2))*expy;
  double xx2=sqrt(m2/pow(energy_.sroot_,2))*expmy;
  double p1[4];
  p1[0]=0.;
  p1[1]=0.;
  p1[2]=xx1*0.5*energy_.sroot_;
  p1[3]=xx1*0.5*energy_.sroot_;
  double p2[4];
  p2[0]=0.;
  p2[1]=0.;
  p2[2]=-xx2*0.5*energy_.sroot_;
  p2[3]=xx2*0.5*energy_.sroot_;

  // First lepton direction: Cos of the polar angle
  double costh=-1.+2.*x[3];
  jac=jac*2.;
  
  //Generate vector boson and lepton azimuthal angles 
  double phi = 2.*M_PI*x[4];
  double phi_lep = 2.*M_PI*x[5];
      
  double p3[4];
  double p4[4];
  double costh_CS;
  if (opts.fixedorder)
    {
      double pV[4];
      pV[0]=0.;
      pV[1]=0.;
      pV[2]=0.5*m*(expy-expmy);
      pV[3]=0.5*m*(expy+expmy);

      // momentum of the first lepton 
      double p4cm[4];
      p4cm[3]=m/2.;
      p4cm[0]=p4cm[3]*sin(acos(costh))*sin(phi_lep);
      p4cm[1]=p4cm[3]*sin(acos(costh))*cos(phi_lep);
      p4cm[2]=p4cm[3]*costh;

      // Boost to go in the Lab frame
      boost_(m,pV,p4cm,p4);

      //  momentum of the second lepton
      p3[3]=pV[3]-p4[3];
      p3[0]=pV[0]-p4[0];
      p3[1]=pV[1]-p4[1];
      p3[2]=pV[2]-p4[2];

      //apply lepton cuts
      if (opts.makelepcuts)
	if (!cuts::lep(p3, p4))
	  {
	    f[0]=0.;
	    return 0;
	  }
      costh_CS = costh;
    }
  else
    {
      //  vector boson momentum: pV(4)^2-pV(1)^2-pV(2)^2-pV(3)^2=m2
      double mt2=m2+qt2;
      double pV[4];
      pV[0]=qt*cos(phi);
      pV[1]=qt*sin(phi);
      pV[2]=0.5*sqrt(mt2)*(expy-expmy);
      pV[3]=0.5*sqrt(mt2)*(expy+expmy);

      // momentum of the first lepton 
      double p4cm[4];
      p4cm[3]=m/2.;
      p4cm[0]=p4cm[3]*sin(acos(costh))*sin(phi_lep);
      p4cm[1]=p4cm[3]*sin(acos(costh))*cos(phi_lep);
      p4cm[2]=p4cm[3]*costh;

      // Boost to go in the Lab frame
      boost_(m,pV,p4cm,p4);

      //  momentum of the second lepton
      p3[3]=pV[3]-p4[3];
      p3[0]=pV[0]-p4[0];
      p3[1]=pV[1]-p4[1];
      p3[2]=pV[2]-p4[2];

      //apply lepton cuts
      if (opts.makelepcuts)
	if (!cuts::lep(p3, p4))
	  {
	    f[0]=0.;
	    return 0;
	  }
  			    
      //Calculate costh according to a qt-recoil prescription
      double kt1,kt2;
  
      //CS frame prescription
      if (opts.qtrec_cs)
	{
	  kt1 = pV[0]/2.;
	  kt2 = pV[1]/2.;
	}

      //MY (DYRES) prescription
      if (opts.qtrec_naive)
	{
	  kt1=(1.+pV[2]/(sqrt(m2)+pV[3]))*pV[0]/2;
	  kt2=(1.+pV[2]/(sqrt(m2)+pV[3]))*pV[1]/2;
	}  

      //alternative k1t = 0 prescription
      if (opts.qtrec_kt0)
	{
	  kt1 = 0;
	  kt2 = 0;
	}

      double zeta1=1./m2/2.*(m2+2.*(pV[0]*kt1+pV[1]*kt2)+sqrt(pow((m2+2.*(pV[0]*kt1+pV[1]*kt2)),2)-4.*mt2*(pow(kt1,2)+pow(kt2,2))));

      double qP1=(pV[3]-pV[2])*energy_.sroot_/2.;
      double qP2=(pV[3]+pV[2])*energy_.sroot_/2.;
  
      double kap1[4];
      kap1[3]=energy_.sroot_/2.*(zeta1*m2/2./qP1+(pow(kt1,2)+pow(kt2,2))/zeta1*qP1/m2/pow(energy_.sroot_,2)*2.);
      kap1[0]=kt1;
      kap1[1]=kt2;
      kap1[2]=energy_.sroot_/2.*(zeta1*m2/2./qP1-(pow(kt1,2)+pow(kt2,2))/zeta1*qP1/m2/pow(energy_.sroot_,2)*2.);

      costh_CS=1.-4.*(kap1[3]*p4[3]-kap1[2]*p4[2]-kap1[1]*p4[1]-kap1[0]*p4[0])/m2;
    }

  //apply resummation switching
  double swtch = switching::swtch(qt, m);
  /*
  //No need to check the switching, since the phase space is generated up to the switching qt limit
  if (swtch < switching::cutoff*switching::tolerance)
  {
    f[0]=0.;
    return 0;
  }
  */

  //Call the counterterm
  int mode = 0;
  dofill_.doFill_ = int(iter==last_iter);
  f[0] = ctint_(costh_CS,m,qt,y,mode,f);
  
  //avoid nans
  if (f[0] != f[0])
    {
      f[0]=0.;
      if (opts.PDFerrors)
	for (int i = 1; i < opts.totpdf; i++)
	  f[i] = 0.;
      return 0;
    }
	   
  //apply switching and jacobian
  f[0] = f[0]*jac*swtch;
  if (opts.PDFerrors)
    for (int i = 1; i < opts.totpdf; i++)
      f[i] = f[i]*jac*swtch;
  
  if (iter==last_iter)
    for (int i = 0; i < opts.totpdf; i++)
      {
	double wt = weight*f[i];
	hists_setpdf_(&i);
	hists_fill_(p3, p4, &wt);
      }

  end_time = clock();
  if (opts.timeprofile)
    cout << setw (3) << "m" << setw(10) << m << setw(4) << "qt" << setw(10) <<  qt
	 << setw (3) << "y" << setw(10) << y << setw(4) << "costh" << setw(10) <<  costh
	 << setw(8) << "result" << setw(10) << f[0]
	 << setw(10) << "tot time" << setw(10) << float( end_time - begin_time ) /  CLOCKS_PER_SEC
	 << endl;
  return 0;
}

integrand_t ctintegrand3d(const int &ndim, const double x[], const int &ncomp, double f[])
//Generates the phase space 4 vectors
//Calculates the ct integrand as a function of m, qt, y
//dOmega integration is factorised in the costh moments
//The integration in alpha and beta is performed inside countdy
{
  clock_t begin_time, end_time;

  begin_time = clock();

  if (opts.PDFerrors)
    for (int i = 1; i < opts.totpdf; i++)
      f[i] = 0.;

  //Jacobian of the change of variables from the unitary hypercube x[3] to the m, y, qt boundaries
  double jac = 1.;

  // Generate the boson invariant mass between the integration boundaries
  double mcut = qtmax/qtcut_.xqtcut_;
  double wsqmin = pow(mmin,2);
  double wsqmax = pow(min(mmax,mcut),2);
  if (wsqmin >= wsqmax)
    {
      f[0]=0.;
      return 0;
    }
  double x1=x[0];
  double m2,wt;
  breitw_(x1,wsqmin,wsqmax,opts.rmass,opts.rwidth,m2,wt);
  double m=sqrt(m2);
  jac=jac*wt;

  //Dynamic scale (not implemented yet)
  //if(dynamicscale) call scaleset(m2)

  //Limit y boundaries to the kinematic limit in y
  double ylim = 0.5*log(pow(energy_.sroot_,2)/m2);
  double ymn = min(max(-ylim, ymin),ylim);
  double ymx = max(min(ylim, ymax),-ylim);
  if (ymn >= ymx)
    {
      f[0]=0.;
      return 0;
    }

  //integrate between ymin and ymax
  double y=ymn+(ymx-ymn)*x[1];
  jac=jac*(ymx-ymn);
  
  //integrate between qtmin and qtmax
  //use qt2 to get the correct jacobian!
  double qtcut = qtcut_.xqtcut_*m;
  double qtmn = max(qtcut, qtmin);
  double cosh2y34=pow((exp(y)+exp(-y))*0.5,2);
  double kinqtlim = sqrt(pow(pow(energy_.sroot_,2)+m*m,2)/(4*pow(energy_.sroot_,2)*cosh2y34)-m*m);
  double switchqtlim = switching::qtlimit(m);
  double qtlim = min(kinqtlim, switchqtlim);
  double qtmx = min(qtlim, qtmax);
  if (qtmn >= qtmx)
    {
      f[0]=0.;
      return 0;
    }

  //In the fixed order calculation, integrate from qtcut to infinity
  if (opts.fixedorder)
    {
      if (qtmin > 0)
	{
	  f[0]=0.;
	  return 0;
	}
      qtmn = qtcut;
      qtmx = 1e10;
    }
  double qtmn2 = pow(qtmn,2);
  double qtmx2 = pow(qtmx,2);
  
  double tiny = 1E-5;
  double a = 1./(1+log(qtmx2/tiny));
  double b = 1./(1+log(qtmn2/tiny));
  double x2 = a + (b-a) * x[2];
  jac = jac * (b-a);
  double qt2=tiny*exp(1./x2 - 1);
  jac=jac*qt2/pow(x2,2);
  
  double qt=sqrt(qt2);
  
  //set global variables to costh, m, qt, y
  if (opts.fixedorder) //In the fixed order calculation evaluate kinematic cuts with qt=0
    setcthmqty(0, m, 0.0001, y);
  else
    setcthmqty(0, m, qt, y);

  //generate boson 4-momentum, with m, qt, y and phi=0
  if (opts.fixedorder) //In the fixed order calculation evaluate kinematic cuts with qt=0
    genV4p(m, 0.0001, y, 0.);
  else
    genV4p(m, qt, y, 0.);

  //  SWITCHING FUNCTIONS is inside countdy
  //  double swtch=1.;
  //  if (qt >= m*3/4.)  swtch=exp(-pow((m*3/4.-qt),2)/pow((m/2.),2)); // GAUSS SWITCH
  //  if (swtch <= 0.01) swtch = 0;
  double swtch = switching::swtch(qt, m);

  //In this point of phase space (m, qt, y) the costh integration is performed by 
  //calculating the 0, 1 and 2 moments of costh
  //that are the integrals int(dcosth dphi1 dphi2), int(costh dcosth dphi1 dphi2), and int(costh^2 dcosth dphi1 dphi2) convoluted with cuts
  //Then the epxressions 1, costh and costh^2 in sigmaij are substituted by these costh moments
  double costh = 0;
  int mode = 1;
  dofill_.doFill_ = 1;

  /*
  //No need to check the switching, since the phase space is generated up to the switching qt limit
  if (swtch < switching::cutoff*switching::tolerance)
  {
    f[0]=0.;
    return 0;
  }
  */

  //evaluate the fixed order expansion of the resummed cross section
  f[0]=ctint_(costh,m,qt,y,mode,f);

  //avoid nans
  if (f[0] != f[0])
    {
      f[0]=0.;
      if (opts.PDFerrors)
	for (int i = 1; i < opts.totpdf; i++)
	  f[i] = 0.;
      return 0;
    }
	   
  f[0] = f[0]*jac; //*swtch; switching function is inside ctint_
  if (opts.PDFerrors)
    for (int i = 1; i < opts.totpdf; i++)
      f[i] = f[i]*jac;

  end_time = clock();
  if (opts.timeprofile)
    cout << setw (3) << "m" << setw(10) << m << setw(4) << "qt" << setw(10) <<  qt
	 << setw(8) << "result" << setw(10) << f[0]
	 << setw(10) << "tot time" << setw(10) << float( end_time - begin_time ) /  CLOCKS_PER_SEC
	 << endl;
  return 0;
}

int ctintegrand2d_cubature_v(unsigned ndim, long unsigned npts, const double x[], void *data, unsigned ncomp, double f[])
{
  //  cout << "parallel " << npts << endl;
#pragma omp parallel for num_threads(opts.cubacores) copyin(a_param_,scale_,facscale_,qcdcouple_)
  for (unsigned i = 0; i < npts; i++)
    {
      // evaluate the integrand for npts points
      double xi[ndim];
      double fi[ncomp];
      for (unsigned j = 0; j < ndim; j++)
	xi[j] = x[i*ndim + j];

      ctintegrand2d(ndim, xi, ncomp, fi);
      
      for (unsigned k = 0; k < ncomp; ++k)
	f[i*ncomp + k] = fi[k];
    }
  return 0;
}

int ctintegrand2d_cubature(unsigned ndim, const double x[], void *data, unsigned ncomp, double f[])
{
  ctintegrand2d(ndim, x, ncomp, f);
  return 0;
}

integrand_t ctintegrand2d(const int &ndim, const double x[], const int &ncomp, double f[])
//Calculates the ct integrand as a function of m, y
//The integration in qt is factorised in LL1, LL2, LL3 and LL4 large logs
//dOmega integration is factorised in the costh moments
//The integration in alpha and beta is performed inside countdy
{
  clock_t begin_time, end_time;

  begin_time = clock();

  if (opts.PDFerrors)
    for (int i = 1; i < opts.totpdf; i++)
      f[i] = 0.;

  //Jacobian of the change of variables from the unitary hypercube x[3] to the m, y, qt boundaries
  double jac = 1.;

  // Generate the boson invariant mass between the integration boundaries
  double mcut = qtmax/qtcut_.xqtcut_;
  double wsqmin = pow(mmin,2);
  double wsqmax = pow(min(mmax,mcut),2);
  if (wsqmin >= wsqmax)
    {
      f[0]=0.;
      return 0;
    }
  double x1=x[0];
  double m2,wt;
  breitw_(x1,wsqmin,wsqmax,opts.rmass,opts.rwidth,m2,wt);
  double m=sqrt(m2);
  jac=jac*wt;

  //Dynamic scale (not implemented yet)
  //if(dynamicscale) call scaleset(m2)

  //Limit y boundaries to the kinematic limit in y
  double ylim = 0.5*log(pow(energy_.sroot_,2)/m2);
  double ymn = min(max(-ylim, ymin),ylim);
  double ymx = max(min(ylim, ymax),-ylim);
  if (ymn >= ymx)
    {
      f[0]=0.;
      return 0;
    }

  //integrate between ymin and ymax
  double y=ymn+(ymx-ymn)*x[1];
  jac=jac*(ymx-ymn);
  
  //factorised integration between qtmin and qtmax
  //use qt2 to get the correct jacobian!
  double qtcut = qtcut_.xqtcut_*m;
  double qtmn = max(qtcut, qtmin);
  double cosh2y34=pow((exp(y)+exp(-y))*0.5,2);
  double kinqtlim = sqrt(pow(pow(energy_.sroot_,2)+m*m,2)/(4*pow(energy_.sroot_,2)*cosh2y34)-m*m);
  double switchqtlim = switching::qtlimit(m);
  double qtlim = min(kinqtlim, switchqtlim);
  double qtmx = min(qtlim, qtmax);
  if (qtmn >= qtmx)
    {
      f[0]=0.;
      return 0;
    }

  //In the fixed order calculation, integrate from qtcut to infinity
  if (opts.fixedorder)
    {
      if (qtmin > 0)
	{
	  f[0]=0.;
	  return 0;
	}
      qtmn = qtcut;
      qtmx = 1e10;
    }

  //set global variables to costh, m, qt, y
  if (opts.fixedorder) //In the fixed order calculation evaluate kinematic cuts with qt=0
    setcthmqty(0, m, 0.001, y);
  else
    setcthmqty(0, m, (qtmn+qtmx)/2., y);

  clock_t qtbt, qtet;
  qtbt = clock();
  ctqtint_(m,y,qtmn,qtmx);
  qtet = clock();
  
  //generate boson 4-momentum, with m, qt, y and phi=0
  if (opts.fixedorder) //In the fixed order calculation evaluate kinematic cuts with qt=0
    genV4p(m, 0.0001, y, 0.);
  else
    genV4p(m, (qtmn+qtmx)/2., y, 0.);

  //In this point of phase space (m, qt, y) the costh integration is performed by 
  //calculating the 0, 1 and 2 moments of costh
  //that are the integrals int(dcosth dphi1 dphi2), int(costh dcosth dphi1 dphi2), and int(costh^2 dcosth dphi1 dphi2) convoluted with cuts
  //Then the epxressions 1, costh and costh^2 in sigmaij are substituted by these costh moments
  double costh = 0;
  int mode = 2;
  dofill_.doFill_ = 1;

  //evaluate the fixed order expansion of the resummed cross section
  double qt = (qtmn+qtmx)/2.;
  clock_t cbt = clock();
  f[0]=ctint_(costh,m,qt,y,mode,f);
  clock_t cet = clock();

  //avoid nans
  if (f[0] != f[0])
    {
      f[0]=0.;
      if (opts.PDFerrors)
	for (int i = 1; i < opts.totpdf; i++)
	  f[i] = 0.;
      return 0;
    }
	   
  f[0] = f[0]*jac;
  if (opts.PDFerrors)
    for (int i = 1; i < opts.totpdf; i++)
      f[i] = f[i]*jac;

  end_time = clock();
  if (opts.timeprofile)
    cout << setw (3) << "m" << setw(10) << m << setw(4) << "y" << setw(10) <<  y
	 << setw(8) << "result" << setw(10) << f[0]
	 << setw(10) << "tot time" << setw(10) << float( end_time - begin_time ) /  CLOCKS_PER_SEC
      	 << setw(10) << "qtint"  << setw(10) << float( qtet - qtbt ) /  CLOCKS_PER_SEC
	 << setw(10) << "ctint"  << setw(10) << float( cet - cbt ) /  CLOCKS_PER_SEC
	 << endl;
  return 0;
}
