#include "resintegr.h"
#include "integr.h"
#include "settings.h"
#include "interface.h"
#include "switch.h"
//#include "plotter.h"

#include <math.h>
#include <iomanip>
#include <iostream>
//#include <vector>
//#include <random>

double const qtcutoff = 0.02;

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
  double qtmn = max(qtcutoff, qtmin);
  double qt=qtmn+(qtmax-qtmn)*x[1];
  jac=jac*(qtmax-qtmn);
  */
    
  //integrate between qtmin and qtmax
  //limit qtmax to the qT kinematical limit, or to the switching function boundary
  double qtmn = max(qtcutoff, qtmin);
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
      rapintegrals_(ymin,ymax,m,nocuts);
      //C++ rewritten resum
      //rapint::integrate(ymin,ymax,m);
      //end
    }
  else
    {
      rapintegrals_(ymn,ymx,m,nocuts);
      //C++ rewritten resum
      //rapint::integrate(ymn,ymx,m);
      //end
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
  if (swtch < 0.01)
    f[0]=0.;
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
  double qtmn = max(qtcutoff, qtmin);
  double qt=qtmn+(qtmax-qtmn)*x[2];
  jac=jac*(qtmax-qtmn);
  */

  //integrate between qtmin and qtmax
  double qtmn = max(qtcutoff, qtmin);
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
  if (swtch < 0.01)
    f[0]=0.;
  else
    //evaluate the resummed cross section
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

  double r[ndim];
  for (int i = 0; i < ndim; i++)
    r[i]=x[i];

  double resumm,tempp,ran2;
  int ii,j,ij;
  double m,m2,qt2,s,wgt;
  double p[5][5],pV[4],p4cm[4],ptemp[4],kap1[5],kap1b[5];
  double pjet[4][12];
  double phi,cos_th,phi_lep;
  int n2,n3;
  double msq,wt,swtch,kk,zeta1,zeta1b,kt1,kt2,mt2,costh_CS,qP1,qP2,ritmx;

  double wsqmin = pow(mmin,2);
  double wsqmax = pow(mmax,2);
  
  double lowintHst=0;
  double lowintHst0=0;

  int jstop = 0;

  //double azloopmax = 1;
  //double azloop = 0;

  ptemp[0]=0.;
  ptemp[1]=0.;
  ptemp[2]=0.;
  ptemp[3]=0.;
  for (int i=0; i < 5; i++)
    {
      p[i][1]=0.;
      p[i][2]=0.;
      p[i][3]=0.;
      p[i][4]=0.;
    }

  // Generate the boson invariant mass

  double jac=1.;
  double x1=r[0];
  breitw_(x1,wsqmin,wsqmax,opts.rmass,opts.rwidth,msq,wt);
  m2=msq;
  m=sqrt(m2);
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
  double y=ymn+(ymx-ymn)*r[1];
  jac=jac*(ymx-ymn);

  //integrate between qtmin and qtmax
  double expy=exp(y);
  double expmy=exp(-y);
  double qtmn = max(qtcutoff, qtmin);
  double cosh2y34=pow((expy+expmy)*0.5,2);
  double kinqtlim = sqrt(pow(pow(energy_.sroot_,2)+m*m,2)/(4*pow(energy_.sroot_,2)*cosh2y34)-m*m);
  double switchqtlim = switching::qtlimit(m);
  double qtlim = min(kinqtlim, switchqtlim);
  double qtmx = min(qtlim, qtmax);
  if (qtmn >= qtmx)
    {
      f[0]=0.;
      return 0;
    }
  double qt=qtmn+(qtmx-qtmn)*r[2];
  jac=jac*(qtmx-qtmn);

  qt2=pow(qt,2);

  double xx1=sqrt(m2/pow(energy_.sroot_,2))*expy;
  double xx2=sqrt(m2/pow(energy_.sroot_,2))*expmy;


  // incoming quarks
  p[1][1]=0.;
  p[1][2]=0.;
  p[1][3]=xx1*0.5*energy_.sroot_;
  p[1][4]=xx1*0.5*energy_.sroot_;
  p[2][1]=0.;
  p[2][2]=0.;
  p[2][3]=-xx2*0.5*energy_.sroot_;
  p[2][4]=xx2*0.5*energy_.sroot_;


  // First lepton direction: Cos of the polar angle
  cos_th=-1.+2.*r[3];
  jac=jac*2.;
  
  mt2=m2+qt2;

  int mode = 0;

  // LOOP over (vector boson and lepton) azimuthal angles 
  //for (int j=0; j < azloopmax; j++)
  //    {
      // Vector boson azimuthal angle
      //      phi = 2.*M_PI*((double)rand()/(double)RAND_MAX);
      // Azimuthal angle of the first lepton in the center of mass frame 
      //phi_lep = 2.*M_PI*((double)rand()/(double)RAND_MAX);

      phi = 2.*M_PI*r[4];
      phi_lep = 2.*M_PI*r[5];
      
      //  vector boson momentum: pV(4)^2-pV(1)^2-pV(2)^2-pV(3)^2=m2
      pV[0]=qt*cos(phi);
      pV[1]=qt*sin(phi);
      pV[2]=0.5*sqrt(mt2)*(exp(y)-exp(-y));
      pV[3]=0.5*sqrt(mt2)*(exp(y)+exp(-y));

      // momentum of the first lepton 
      p4cm[3]=m/2.;
      p4cm[0]=p4cm[3]*sin(acos(cos_th))*sin(phi_lep);
      p4cm[1]=p4cm[3]*sin(acos(cos_th))*cos(phi_lep);
      p4cm[2]=p4cm[3]*cos_th;

      // Boost to go in the Lab frame
      boost_(m,pV,p4cm,ptemp);
      p[4][1]=ptemp[0];
      p[4][2]=ptemp[1];
      p[4][3]=ptemp[2];
      p[4][4]=ptemp[3];

      //  momentum of the second lepton
      p[3][4]=pV[3]-p[4][4];
      p[3][1]=pV[0]-p[4][1];
      p[3][2]=pV[1]-p[4][2];
      p[3][3]=pV[2]-p[4][3];

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

      zeta1=1./m2/2.*(m2+2.*(pV[0]*kt1+pV[1]*kt2)+sqrt(pow((m2+2.*(pV[0]*kt1+pV[1]*kt2)),2)-4.*mt2*(pow(kt1,2)+pow(kt2,2))));

      qP1=(pV[3]-pV[2])*energy_.sroot_/2.;
      qP2=(pV[3]+pV[2])*energy_.sroot_/2.;

      kap1[4]=energy_.sroot_/2.*(zeta1*m2/2./qP1+(pow(kt1,2)+pow(kt2,2))/zeta1*qP1/m2/pow(energy_.sroot_,2)*2.);
      kap1[1]=kt1;
      kap1[2]=kt2;
      kap1[3]=energy_.sroot_/2.*(zeta1*m2/2./qP1-(pow(kt1,2)+pow(kt2,2))/zeta1*qP1/m2/pow(energy_.sroot_,2)*2.);

      costh_CS=1.-4.*(kap1[4]*p[4][4]-kap1[3]*p[4][3]-kap1[2]*p[4][2]-kap1[1]*p[4][1])/m2;

       // see whether this point will pass cuts - if it will not, do not
       // bother calculating the matrix elements for it, instead bail out
      for (int i=0; i < 4; i++)
	{
	  pjet[i][0] = p[1][i+1];
	  pjet[i][1] = p[2][i+1];
	  pjet[i][2] = p[3][i+1];
	  pjet[i][3] = p[4][i+1];
	}
      int njet = 0;
       //       cout << qt << "  " <<  y << "  " <<  m << "  " << sqrt(p[3][1]*p[3][1] + p[3][2] * p[3][2]) << endl;
      if (opts.makelepcuts)
	if (cuts_(pjet,njet)) // continue;
	  {
	    f[0]=0.;
	    return 0;
	  }
			    
      //  SWITCHING FUNCTIONS
      //      swtch=1.;
      //      if (qt >= m*3/4.)  swtch=exp(-pow((m*3/4.-qt),2)/pow((m/2.),2)); // GAUSS SWITCH
      swtch = switching::swtch(qt, m);

      if (swtch <= 0.01) //break;
	{
	  f[0]=0.;
	  return 0;
	}

      if (jstop != 1)
	{
	  jstop=1;

	  //Call to the resummation part
	  if (swtch < 0.01)
	    tempp=0.;
	  else
	    tempp=resumm_(costh_CS,m,qt,y,mode)/(8./3.);  // CS PRESCRIPTION
	    //tempp=resumm_(cos_th,m,qt,y,mode)/(8./3.);
	   if (tempp != tempp)
	     tempp=0.;    //  avoid nans
	   
	   lowintHst0=jac*tempp;
	   lowintHst0=lowintHst0*swtch; // SWITCHING
	 }
      if (iter==4){
  	double wt = weight*lowintHst0;///azloopmax;
	hists_fill_(p[4]+1,p[3]+1,&wt);
        //hists_AiTest_(pjet,p4cm,&m,&qt,&y,&costh_CS,&phi_lep,&phi,&wt,&lowintHst0);
      } 
      //azloop=azloop+1;
      //}

      lowintHst=lowintHst0;//*float(azloop)/float(azloopmax);
  
  f[0] = lowintHst;

  return 0;
}
