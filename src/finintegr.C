#include <iostream>
#include <LHAPDF/LHAPDF.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <ctime>

#include "interface.h"
#include "finintegr.h"
#include "integr.h"
#include "settings.h"

using namespace std;

const int last_iter=4;

integrand_t lowintegrand(const int &ndim, const double x[], const int &ncomp, double f[],
                        void* userdata, const int &nvec, const int &core,
                        double &weight, const int &iter)
{
  
  double rlo[22];
  for (int i = 0; i < ndim; i++)
    rlo[i]=x[i];


  f[0] = lowint_(rlo,weight);
  
  return 0;
}

integrand_t realintegrand(const int &ndim, const double x[], const int &ncomp, double f[],
                        void* userdata, const int &nvec, const int &core,
                        double &weight, const int &iter)
{
  
  double rre[22];
  for (int i = 0; i < ndim; i++)
    rre[i]=x[i];

  dofill_.doFill_ = int(iter==last_iter);
  f[0] = realint_(rre,weight);
  
  return 0;
}

integrand_t virtintegrand(const int &ndim, const double x[], const int &ncomp, double f[],
                        void* userdata, const int &nvec, const int &core,
                        double &weight, const int &iter)
{
  
  double rvi[22];
  for (int i = 0; i < ndim; i++)
    rvi[i]=x[i];
  rvi[9] = rvi[7];
  
  dofill_.doFill_ = int(iter==last_iter);
  f[0]  = virtint_(rvi,weight);
  
  return 0;
}


integrand_t ctintegrand(const int &ndim, const double x[], const int &ncomp, double f[],
                        void* userdata, const int &nvec, const int &core,
                        double &weight, const int &iter)
{
    //double wgt = 1;
  //here generate the phase space according to x[], and pass the p vector to countint_
  
  double rct[22];
  for (int i = 0; i < ndim; i++)
    rct[i]=x[i];
  rct[8] = rct[1];

  dofill_.doFill_ = int(iter==last_iter);
  f[0] = countint_(rct,weight);

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

  //Jacobian of the change of variables from the unitary hypercube x[3] to the m, y, qt boundaries
  double jac = 1.;

  // Generate the boson invariant mass between the integration boundaries
  double mcut = qtmax/qtcut_.xqtcut_;
  double wsqmin = pow(mmin,2);
  double wsqmax = pow(min(mmax,mcut),2);
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

  //integrate between ymin and ymax
  double y=ymn+(ymx-ymn)*x[1];
  jac=jac*(ymx-ymn);
  
  //integrate between qtmin and qtmax
  //use qt2 to get the correct jacobian!
  double qtcut = qtcut_.xqtcut_*m;
  double qtmn = max(qtcut, qtmin);
  double cosh2y34=pow((exp(y)+exp(-y))*0.5,2);
  double qtlim = sqrt(pow(pow(energy_.sroot_,2)+m*m,2)/(4*pow(energy_.sroot_,2)*cosh2y34)-m*m);
  double qtmx = min(qtlim, qtmax);
  double qtmn2 = pow(qtmn,2);
  double qtmx2 = pow(qtmx,2);
  if (qtmn >= qtmx)
    {
      f[0]=0.;
      return 0;
    }
  
  double tiny = 1E-5;
  double a = 1./(1+log(qtmx2/tiny));
  double b = 1./(1+log(qtmn2/tiny));
  double x2 = a + (b-a) * x[2];
  jac = jac * (b-a);
  double qt2=tiny*exp(1./x2 - 1);
  jac=jac*qt2/pow(x2,2);
  
  double qt=sqrt(qt2);
  
  //set global variables to costh, m, qt, y
  set(0, m, qt, y);

  //generate boson 4-momentum, with m, qt, y and phi=0
  genV4p(m, qt, y, 0.);

  //  SWITCHING FUNCTIONS is inside countdy
  double swtch=1.;
  if (qt >= m*3/4.)  swtch=exp(-pow((m*3/4.-qt),2)/pow((m/2.),2)); // GAUSS SWITCH
  if (swtch <= 0.01) swtch = 0;

  //In this point of phase space (m, qt, y) the costh integration is performed by 
  //calculating the 0, 1 and 2 moments of costh
  //that are the integrals int(dcosth dphi1 dphi2), int(costh dcosth dphi1 dphi2), and int(costh^2 dcosth dphi1 dphi2) convoluted with cuts
  //Then the epxressions 1, costh and costh^2 in sigmaij are substituted by these costh moments
  double costh = 0;
  int mode = 1;
  if (swtch < 0.01)
    f[0]=0.;
  else
    //evaluate the fixed order expansion of the resummed cross section
    f[0]=countterm_(costh,m,qt,y,mode);

  if (f[0] != f[0])
    f[0]=0.;  //avoid nans
	   
  f[0] = f[0]*jac; //*swtch; switching function is inside countterm_

  end_time = clock();
  if (opts.timeprofile)
    cout << setw (3) << "m" << setw(10) << m << setw(4) << "qt" << setw(10) <<  qt
	 << setw(8) << "result" << setw(10) << f[0]
	 << setw(10) << "tot time" << setw(10) << float( end_time - begin_time ) /  CLOCKS_PER_SEC
	 << endl;
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

  //Jacobian of the change of variables from the unitary hypercube x[3] to the m, y, qt boundaries
  double jac = 1.;

  // Generate the boson invariant mass between the integration boundaries
  double mcut = qtmax/qtcut_.xqtcut_;
  double wsqmin = pow(mmin,2);
  double wsqmax = pow(min(mmax,mcut),2);
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

  //integrate between ymin and ymax
  double y=ymn+(ymx-ymn)*x[1];
  jac=jac*(ymx-ymn);
  
  //factorised integration between qtmin and qtmax
  //use qt2 to get the correct jacobian!
  double qtcut = qtcut_.xqtcut_*m;
  double qtmn = max(qtcut, qtmin);
  double cosh2y34=pow((exp(y)+exp(-y))*0.5,2);
  double qtlim = sqrt(pow(pow(energy_.sroot_,2)+m*m,2)/(4*pow(energy_.sroot_,2)*cosh2y34)-m*m);
  double qtmx = min(qtlim, qtmax);
  if (qtmn >= qtmx)
    {
      f[0]=0.;
      return 0;
    }
  
  //set global variables to costh, m, qt, y
  set(0, m, (qtmn+qtmx)/2., y);

  clock_t qtbt, qtet;
  qtbt = clock();
  ctqtint_(m,y,qtmn,qtmx);
  qtet = clock();
  
  //generate boson 4-momentum, with m, qt, y and phi=0
  genV4p(m, (qtmn+qtmx)/2., y, 0.);

  //In this point of phase space (m, qt, y) the costh integration is performed by 
  //calculating the 0, 1 and 2 moments of costh
  //that are the integrals int(dcosth dphi1 dphi2), int(costh dcosth dphi1 dphi2), and int(costh^2 dcosth dphi1 dphi2) convoluted with cuts
  //Then the epxressions 1, costh and costh^2 in sigmaij are substituted by these costh moments
  double costh = 0;
  int mode = 2;

  //evaluate the fixed order expansion of the resummed cross section
  double qt = (qtmn+qtmx)/2.;
  clock_t cbt = clock();
  f[0]=countterm_(costh,m,qt,y,mode);
  clock_t cet = clock();

  if (f[0] != f[0])
    f[0]=0.;  //avoid nans
	   
  f[0] = f[0]*jac;

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

int binner_(double p3[4], double p4[4])
{
    if (opts.HackBinnerToFiller) {
        return true;
    } else {

        double qt = sqrt((float)pow(p3[0]+p4[0],2)+pow(p3[1]+p4[1],2));
        if (qt < qtmin || qt > qtmax)
            return false;

        double y = 0.5 *log((p3[3] + p4[3] + p3[2] + p4[2]) / (p3[3] + p4[3] - p3[2] - p4[2]));
        if (y < ymin || y > ymax)
            return false;

        //  cout << "qt " << qt << " y " << y << endl;
    }
    return true;
}
