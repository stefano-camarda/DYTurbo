#include "resintegr.h"
#include "integr.h"
#include "settings.h"
#include "interface.h"
#include "switch.h"
#include "resint.h"
#include "rapint.h"
//#include "plotter.h"

#include <math.h>
#include <iomanip>
#include <iostream>
#include <omp.h>

int resintegrand2d_cubature_v(unsigned ndim, long unsigned npts, const double x[], void *data, unsigned ncomp, double f[])
{
  //The current issue with openmp parallelisation, is that the resummation integrand has many
  //global variables, and the use of these variables has a race
  
  //  cout << "parallel " << npts << endl;

#pragma omp parallel for num_threads(opts.cubacores) copyin(a_param_,rlogs_,resint::a,resint::loga,resint::rloga)
  for (unsigned i = 0; i < npts; i++)
    {
      //      int nThreads = omp_get_num_threads();
      //      printf("Total number of threads: %d\n", nThreads);
      
      // evaluate the integrand for npts points
      double xi[ndim];
      double fi[ncomp];
      for (unsigned j = 0; j < ndim; j++)
	xi[j] = x[i*ndim + j];

      resintegrand2d(ndim, xi, ncomp, fi);
      
      for (unsigned k = 0; k < ncomp; ++k)
	f[i*ncomp + k] = fi[k];
    }
  return 0;
}

int resintegrand2d_cubature(unsigned ndim, const double x[], void *data, unsigned ncomp, double f[])
{
  resintegrand2d(ndim, x, ncomp, f);
  return 0;
}

integrand_t resintegrand2d(const int &ndim, const double x[], const int &ncomp, double f[])
{
  clock_t begin_time, end_time;

  begin_time = clock();

  //Jacobian of the change of variables from the unitary hypercube x[3] to the m, y, qt boundaries
  double jac = 1.;

  // Generate the boson invariant mass between the integration boundaries
  double wsqmin = pow(mmin,2);
  double wsqmax = pow(mmax,2);
  double x1=x[0];
  double m2,wt;
  breitw_(x1,wsqmin,wsqmax,opts.rmass,opts.rwidth,m2,wt);
  double m=sqrt(m2);
  jac=jac*wt;
  //Dynamic scale (not implemented yet)
  //if(dynamicscale) call scaleset(m2)

  //Limit the y boundaries to the kinematic limit in y
  double ylim = 0.5*log(pow(energy_.sroot_,2)/m2);
  double ymn = min(max(-ylim, ymin),ylim);
  double ymx = max(min(ylim, ymax),-ylim);
  if (ymn >= ymx)
    {
      f[0]=0.;
      return 0;
    }

  //integrate between qtmin and qtmax
  /*
  double qtmn = max(opts.qtcutoff, qtmin);
  double qt=qtmn+(qtmax-qtmn)*x[1];
  jac=jac*(qtmax-qtmn);
  */
    
  //integrate between qtmin and qtmax
  //limit qtmax to the qT kinematical limit, or to the switching function boundary
  double qtmn = max(opts.qtcutoff, qtmin);
  double miny;
  if (ymn * ymx <= 0)
    miny = 0;
  else miny = min(fabs(ymn),fabs(ymx));
  double cosh2y34=pow((exp(miny)+exp(-miny))*0.5,2);
  double kinqtlim = sqrt(pow(pow(energy_.sroot_,2)+m*m,2)/(4*pow(energy_.sroot_,2)*cosh2y34)-m*m);
  double switchqtlim = switching::qtlimit(m);
  double qtlim = min(kinqtlim, switchqtlim);
  double qtmx = min(qtlim, qtmax);
  if (qtmn >= qtmx)
    {
      f[0]=0.;
      return 0;
    }
  //double qt=qtmn+(qtmx-qtmn)*x[1];
  //jac=jac*(qtmx-qtmn);
  /*double qtmn2 = pow(qtmn,2);
  double qtmx2 = pow(qtmx,2);
  double tiny = 1E-5;
  double a = 1./(1+log(qtmx2/tiny));
  double b = 1./(1+log(qtmn2/tiny));
  double x2 = a + (b-a) * x[1];
  jac = jac * (b-a);
  double qt2=tiny*exp(1./x2 - 1);
  jac=jac*qt2/pow(x2,2);
  double qt=sqrt(qt2);
  jac=jac*0.5/qt;*/

  /*double tiny = 1E-5;
  double a = 1./(1+log(qtmx/tiny));
  double b = 1./(1+log(qtmn/tiny));
  double x2 = a + (b-a) * x[1];
  jac = jac * (b-a);
  double qt=tiny*exp(1./x2 - 1);
  jac=jac*qt/pow(x2,2);*/

  /*  double exp = 1./3.;
  double a = pow(qtmn,exp);
  double b = pow(qtmx,exp);
  double x2=a+(b-a)*x[1];
  jac=jac*(b-a);
  double qt = pow(x2,1./exp);
  jac=jac*pow(x2,1./exp-1)/exp;*/

  /*double a = log(qtmn);
  double b = log(qtmx);
  double x2=a+(b-a)*x[1];
  jac=jac*(b-a);
  double qt = exp(x2);
  jac=jac*qt;*/

  double esp = 1./2.;
  double tiny = 1E-3;
  double a = pow(log(qtmn/tiny),esp);
  double b = pow(log(qtmx/tiny),esp);
  double x2=a+(b-a)*x[1];
  jac=jac*(b-a);
  double qt = tiny*exp(pow(x2,1./esp));
  jac=jac*qt*pow(x2,1./esp-1)/esp;

  /*double tiny = 1E-3;
  double a = log(log(qtmn/tiny));
  double b = log(log(qtmx/tiny));
  double x2=a+(b-a)*x[1];
  jac=jac*(b-a);
  double qt = tiny*exp(exp(x2));
  jac=jac*qt*exp(x2);*/
  
  /*double base = 100000.;
  double a = log(qtmn)/log(base);
  double b = log(qtmx)/log(base);
  double x2=a+(b-a)*x[1];
  jac=jac*(b-a);
  double qt = exp(x2*log(base));
  jac=jac*qt*log(base);*/
  
  //set global variables to m, qt
  setmqty(m, qt, 0);

  //generate boson 4-momentum, with m, qt, y and phi=0
  genV4p(m, qt, 0, 0.);

  //Perform quadrature rule integration in rapidity and nested costh integration
  int nocuts = !opts.makelepcuts;

  clock_t ybt, yet;
  ybt = clock();
  //there is a potential issue here, when lepton cuts are applied
  //the rapidity dependent exponential are cached assuming integration between ymin and ymax
  //for consistency, has to keep the integration between ymin and ymax
  if (opts.makelepcuts)
    {
      if (opts.resumcpp)
	{
	  //C++ resum
	  rapint::allocate();
	  rapint::integrate(ymin,ymax,m);
	  //end C++ resum
	}
      else
	rapintegrals_(ymin,ymax,m,nocuts);
    }
  else
    {
      if (opts.resumcpp)
	{
	  //C++ rewritten resum
	  rapint::allocate();
	  rapint::integrate(ymn,ymx,m);
	  //end C++ resum
	}
      else
	rapintegrals_(ymn,ymx,m,nocuts);
    }
  yet = clock();
  
  setmqty(m, qt, 0);
  genV4p(m, qt, 0, 0.);

  //  SWITCHING FUNCTIONS
  //  double swtch=1.;
  //  if (qt >= m*3/4.)  swtch=exp(-pow((m*3/4.-qt),2)/pow((m/2.),2)); // GAUSS SWITCH
  //  if (swtch <= 0.01) swtch = 0;

  double swtch = switching::swtch(qt, m);
  
  //Call to the resummation part
  double costh = 0;
  double y = 0.5*(ymx+ymn); //needed to get the correct IFIT index of the approximate PDF
  int mode = 2;

  clock_t rbt = clock();

  /*
  //No need to check the switching, since the phase space is generated up to the switching qt limit
  if (swtch < switching::cutoff*switching::tolerance)
  {
    f[0]=0.;
    return 0;
  }
  */

  if (opts.resumcpp)
    f[0]=resint::rint(costh,m,qt,y,mode)/(8./3.);
  else
    f[0]=resumm_(costh,m,qt,y,mode)/(8./3.);
  clock_t ret = clock();

  if (f[0] != f[0])
    f[0]=0.;  //avoid nans
	   
  f[0] = f[0]*jac*swtch;
  end_time = clock();
  if (opts.timeprofile)
    cout << setw (3) << "m" << setw(10) << m << setw(4) << "qt" << setw(10) <<  qt
	 << setw(8) << "result" << setw(10) << f[0]
	 << setw(10) << "tot time" << setw(10) << float( end_time - begin_time ) /  CLOCKS_PER_SEC
	 << setw(10) << "rapint"  << setw(10) << float( yet - ybt ) /  CLOCKS_PER_SEC
	 << setw(10) << "resumm"  << setw(10) << float( ret - rbt ) /  CLOCKS_PER_SEC
	 << endl;
  if (opts.resumcpp)
    rapint::free();
  return 0;
}


integrand_t resintegrand3d(const int &ndim, const double x[], const int &ncomp, double f[])
{
  clock_t begin_time, end_time;

  begin_time = clock();

  //Jacobian of the change of variables from the unitary hypercube x[3] to the m, y, qt boundaries
  double jac = 1.;

  // Generate the boson invariant mass between the integration boundaries
  double wsqmin = pow(mmin,2);
  double wsqmax = pow(mmax,2);
  double x1=x[0];
  double m2,wt;
  breitw_(x1,wsqmin,wsqmax,opts.rmass,opts.rwidth,m2,wt);
  double m=sqrt(m2);
  jac=jac*wt;

  /*
  //split the mass in 3 pieces,
  //and unweight the breit wigner only from -5 width to +5 width
  //no sense, I have to split the mass before entering the integration
  //and set breit wigner unweighting or not
  double bwmn = opts.rmass -5*opts.rwidth;
  double bwmx = opts.rmass +5*opts.rwidth;
  double bwmnsq = pow(bwmn,2);
  double bwmxsq = pow(bwmx,2);
  double m,m2,wt;
  double xl = 0.25;
  double xu = 0.75;
  if (x[0] < xl)
    {
      double x1=x[0]/(xl-0.);
      m=mmin+(bwmn-mmin)*x[0];
      jac=jac*(bwmn-mmin)/(xl-0.);
    }
  else if (x[0] > xu)
    {
      double x1=(x[0]-xu)/(1.-xu);
      m=bwmx+(mmax-bwmx)*x1;
      jac=jac*(mmax-bwmx)/(1.-xu);
    }
  else
    {
      double xbw=(x[0]-xl)/(xu-xl);
      breitw_(xbw,bwmnsq,bwmxsq,opts.rmass,opts.rwidth,m2,wt);
      m=sqrt(m2);
      wt=wt/(xu-xl);
      jac=jac*wt;
    }
  */

  //  double m=mmin+(mmax-mmin)*x[0];
  //  jac=jac*(mmax-mmin);

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

  /*
  //Integrate the full phase space
  double y=-ylim+2.*ylim*x[1];
  jac=jac*2.*ylim;
  */


  //integrate up to the qT kinematical limit
  /*
  double dexpy=exp(y);
  double dexpmy=exp(-y);
  double cosh2y=pow(((dexpy+dexpmy)*0.5),2);
  double qtmax=sqrt(pow((pow(energy_.sroot_,2)+m2),2)/(4.*pow(energy_.sroot_,2)*cosh2y) - m2);
  double qtmin=0.1;
  double qt=qtmin+qtmax*x[2];
  jac=jac*(qtmax);
  */

  //integrate between qtmin and qtmax
  /*
  double qtmn = max(opts.qtcutoff, qtmin);
  double qt=qtmn+(qtmax-qtmn)*x[2];
  jac=jac*(qtmax-qtmn);
  */

  //integrate between qtmin and qtmax
  double qtmn = max(opts.qtcutoff, qtmin);
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
  //  double qt=qtmn+(qtmx-qtmn)*x[2];
  //  jac=jac*(qtmx-qtmn);
  double esp = 1./2.;
  double tiny = 1E-3;
  double a = pow(log(qtmn/tiny),esp);
  double b = pow(log(qtmx/tiny),esp);
  double x2=a+(b-a)*x[2];
  jac=jac*(b-a);
  double qt = tiny*exp(pow(x2,1./esp));
  jac=jac*qt*pow(x2,1./esp-1)/esp;

  
  //set global variables to m, qt, y
  setmqty(m, qt, y);

  //generate boson 4-momentum, with m, qt, y and phi=0
  genV4p(m, qt, y, 0.);

  //  SWITCHING FUNCTIONS
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

  /*
  //No need to check the switching, since the phase space is generated up to the switching qt limit
  if (swtch < switching::cutoff*switching::tolerance)
  {
    f[0]=0.;
    return 0;
  }
  */
  
  //evaluate the resummed cross section
  if (opts.resumcpp)
    f[0]=resint::rint(costh,m,qt,y,mode)/(8./3.);
  else
    f[0]=resumm_(costh,m,qt,y,mode)/(8./3.);

  if (f[0] != f[0])
    f[0]=0.;  //avoid nans

	   
  f[0] = f[0]*jac*swtch;

  end_time = clock();

  if (opts.timeprofile)
    cout << setw (3) << "m" << setw(10) << m
	 << setw(4) << "qt" << setw(10) <<  qt
      	 << setw(4) << "y" << setw(10) <<  y
	 << setw(8) << "result" << setw(10) << f[0]
	 << setw(10) << "tot time" << setw(10) << float( end_time - begin_time ) /  CLOCKS_PER_SEC
	 << endl;
  return 0;
}

//Perform the integration as in the original dyres code
integrand_t resintegrandMC(const int &ndim, const double x[], const int &ncomp, double f[],
			   void* userdata, const int &nvec, const int &core,
			   double &weight, const int &iter)
{
  clock_t begin_time, end_time;

  begin_time = clock();
  
  //Jacobian of the change of variables from the unitary hypercube x[6] to the m, y, qt, phi, costh, philep boundaries
  double jac = 1.;

  // Generate the boson invariant mass
  double wsqmin = pow(mmin,2);
  double wsqmax = pow(mmax,2);
  if (wsqmin >= wsqmax)
    {
      f[0]=0.;
      return 0;
    }
  double x1=x[0];
  double m2, wt;
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
  double expy=exp(y);
  double expmy=exp(-y);
  double cosh2y34=pow((expy+expmy)*0.5,2);
  double qtmn = max(opts.qtcutoff, qtmin);
  double kinqtlim = sqrt(pow(pow(energy_.sroot_,2)+m*m,2)/(4*pow(energy_.sroot_,2)*cosh2y34)-m*m);
  double switchqtlim = switching::qtlimit(m);
  double qtlim = min(kinqtlim, switchqtlim);
  double qtmx = min(qtlim, qtmax);
  if (qtmn >= qtmx)
    {
      f[0]=0.;
      return 0;
    }
  double qt=qtmn+(qtmx-qtmn)*x[2];
  jac=jac*(qtmx-qtmn);
  

  /********************** DYRES check ************************************/
  /*  //integrate between qtmin and qtmax
  double qtmn = 0.1;
  double qtmx = sqrt(pow(pow(energy_.sroot_,2)+m*m,2)/(4*pow(energy_.sroot_,2)*cosh2y34)-m*m);
  double qt=qtmn+qtmx*x[2];
  if (qt < qtmin || qt >= qtmax)
    {
      f[0]=0.;
      return 0;
    }
    jac=jac*qtmx;*/
  /******************************************************************/
  
  double qt2=pow(qt,2);

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
  double p4[4];
  boost_(m,pV,p4cm,p4);

  //  momentum of the second lepton
  double p3[4];
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

  double costh_CS=1.-4.*(kap1[3]*p4[3]-kap1[2]*p4[2]-kap1[1]*p4[1]-kap1[0]*p4[0])/m2;

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

  //skip PDF loop in the preconditioning phase
  int Npdf = (iter!=last_iter ? 1 : opts.totpdf);
  for(int ipdf=0; ipdf<Npdf; ipdf++){
      // Set PDF
      //setpdf_(&ipdf);
      //setmellinpdf_(&ipdf);
      //hists_setpdf_(&ipdf);
      //Call the resummation integrand
      int mode = 0;
      f[ipdf] = resumm_(costh_CS,m,qt,y,mode)/(8./3.);
      //avoid nans
      if (f[ipdf] != f[ipdf])
          f[ipdf] = 0.;
      //apply switching and jacobian
      f[ipdf] = f[ipdf]*jac*swtch;
      if (iter==last_iter){
          double wt = weight*f[ipdf];
          //hists_fill_(p4,p3,&wt);
          hists_fill_(p3,p4,&wt);
          //hists_AiTest_(pjet,p4cm,&m,&qt,&y,&costh_CS,&phi_lep,&phi,&wt,&lowintHst0);
      }
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
