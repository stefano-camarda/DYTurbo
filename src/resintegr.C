#include "resintegr.h"
#include "omegaintegr.h"
#include "phasespace/phasespace.h"
#include "settings.h"
#include "interface.h"
#include "switch.h"
#include "resum/resint.h"
#include "resum/rapint.h"
#include "cubacall.h"
#include "isnan.h"

#include <math.h>
#include <iomanip>
#include <iostream>
#include <omp.h>

#include "histo/KinematicCuts.h"

int resintegrand2d_cubature_v(unsigned ndim, long unsigned npts, const double x[], void *data, unsigned ncomp, double f[])
{
  tell_to_grid_we_are_alive();
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
  tell_to_grid_we_are_alive();
  resintegrand2d(ndim, x, ncomp, f);
  return 0;
}

integrand_t resintegrand2d(const int &ndim, const double x[], const int &ncomp, double f[])
{
  tell_to_grid_we_are_alive();
  clock_t begin_time, end_time;

  begin_time = clock();

  //Jacobian of the change of variables from the unitary hypercube x[3] to the m, y, qt boundaries
  double jac = 1.;
  bool status = true;

  double r2[2] = {x[0], x[1]};
  status = phasespace::gen_mqt(r2, jac, false, true);
  if (!status)
    {
      f[0] = 0.;
      return 0;
    }
  //Limit the y boundaries to the kinematic limit in y
  double ylim = 0.5*log(pow(opts.sroot,2)/phasespace::m2);
  double ymn = min(max(-ylim, phasespace::ymin),ylim);
  double ymx = max(min(ylim, phasespace::ymax),-ylim);
  if (ymn >= ymx)
    {
      f[0]=0.;
      return 0;
    }

  phasespace::set_phiV(0);

  //Dynamic scale --> is this still needed? probably it is for the fortran version of the code
  if (opts.dynamicscale)
    scaleset_(phasespace::m2);

  //Perform quadrature rule integration in rapidity and semi-analitical costh, phi_lep integration
  int nocuts = !opts.makecuts;

  clock_t ybt, yet;
  ybt = clock();
  //there is a potential issue here, when lepton cuts are applied
  //the rapidity dependent exponential are cached assuming integration between ymin and ymax
  //for consistency, has to keep the integration between ymin and ymax
  if (opts.makecuts)
    {
      if (opts.resumcpp)
	{
	  //C++ resum
	  rapint::allocate();
	  rapint::integrate(phasespace::ymin,phasespace::ymax,phasespace::m);
	  //end C++ resum
	}
      else
	rapintegrals_(phasespace::ymin,phasespace::ymax,phasespace::m,nocuts);
    }
  else
    {
      if (opts.resumcpp)
	{
	  //C++ rewritten resum
	  rapint::allocate();
	  rapint::integrate(ymn,ymx,phasespace::m);
	  //end C++ resum
	}
      else
	rapintegrals_(ymn,ymx,phasespace::m,nocuts);
    }
  yet = clock();
  
  //Are the following three lines needed?
  //phasespace::set_m(m);
  //phasespace::set_qt(qt);
  //phasespace::set_phiV(0);

  //Are the following two lines needed?
  //phasespace::set_mqtyphi(m, qt, 0);
  //omegaintegr::genV4p(); // --> not needed, it is regenerated in rapint

  //  SWITCHING FUNCTIONS
  //  double swtch=1.;
  //  if (qt >= m*3/4.)  swtch=exp(-pow((m*3/4.-qt),2)/pow((m/2.),2)); // GAUSS SWITCH
  //  if (swtch <= 0.01) swtch = 0;

  double swtch = switching::swtch(phasespace::qt, phasespace::m);
  
  //Call the resummed cross section
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
    f[0]=resint::rint(costh,phasespace::m,phasespace::qt,y,mode)/(8./3.);
  else
    f[0]=resumm_(costh,phasespace::m,phasespace::qt,y,mode)/(8./3.);
  clock_t ret = clock();

  if (isnan_ofast(f[0]))
    f[0]=0.;  //avoid nans

  f[0] = f[0]*jac*swtch;
  end_time = clock();
  if (opts.timeprofile)
    cout << setw (3) << "m" << setw(10) << phasespace::m << setw(4) << "qt" << setw(10) <<  phasespace::qt
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
  tell_to_grid_we_are_alive();
  clock_t begin_time, end_time;

  begin_time = clock();

  //Jacobian of the change of variables from the unitary hypercube x[3] to the m, y, qt boundaries
  double jac = 1.;
  bool status = true;

  double r[3] = {x[0], x[1], x[2]};
  status = phasespace::gen_myqt(r, jac, false, true); //(qtcut = false, qtswitching = true)
  if (!status)
    {
      f[0] = 0.;
      return 0;
    }
  phasespace::set_phiV(0.);

  //generate boson 4-momentum, with m, qt, y and phi=0
  omegaintegr::genV4p();
  
  //  SWITCHING FUNCTIONS
  //  double swtch=1.;
  //  if (qt >= m*3/4.)  swtch=exp(-pow((m*3/4.-qt),2)/pow((m/2.),2)); // GAUSS SWITCH
  //  if (swtch <= 0.01) swtch = 0;
  double swtch = switching::swtch(phasespace::qt, phasespace::m);

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

  //Dynamic scale
  if (opts.dynamicscale)
    scaleset_(phasespace::m2);

  //evaluate the resummed cross section
  if (opts.resumcpp)
    f[0]=resint::rint(costh,phasespace::m,phasespace::qt,phasespace::y,mode)/(8./3.);
  else
    f[0]=resumm_(costh,phasespace::m,phasespace::qt,phasespace::y,mode)/(8./3.);

  if (isnan_ofast(f[0]))
    f[0]=0.;  //avoid nans

	   
  f[0] = f[0]*jac*swtch;

  end_time = clock();

  if (opts.timeprofile)
    cout << setw (3) << "m" << setw(10) << phasespace::m
	 << setw(4) << "qt" << setw(10) <<  phasespace::qt
      	 << setw(4) << "y" << setw(10) <<  phasespace::y
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
  tell_to_grid_we_are_alive();
  clock_t begin_time, end_time;

  begin_time = clock();

  //Jacobian of the change of variables from the unitary hypercube x[3] to the m, y, qt boundaries
  double jac = 1.;
  bool status = true;

  double r[3] = {x[0], x[1], x[2]};
  status = phasespace::gen_myqt(r, jac, false, true); //(qtcut = false, qtswitching = true)
  if (!status)
    {
      f[0] = 0.;
      return 0;
    }

  //generate boson 4-momentum, with m, qt, y, and phi
  phasespace::set_phiV(-M_PI+2.*M_PI*x[5]);

  //  phasespace::calcexpy();
  //  phasespace::calcmt();

  /////////////////////
  //This function set the RF axes according to the selected qt-recoil prescription
  omegaintegr::genV4p();

  //alternative, generate in naive RF, and calculate costh with the selected qt-recoil prescription
  //phasespace::genRFaxes(phasespace::naive);
  //phasespace::genV4p();
  /////////////////////
  
  double r2[2] = {x[3], x[4]};
  phasespace::gen_costhphi(r2, jac);
  phasespace::genl4p();

  //apply lepton cuts
  if (opts.makecuts)
    if (!Kinematics::Cuts::KeepThisEvent(phasespace::p3, phasespace::p4))
      {
	f[0]=0.;
	return 0;
      }
  
  /////////////////////
  double costh_CS = phasespace::costh;

  //alternative, generate in naive RF, and calculate costh with the selected qt-recoil prescription
  //double costh_CS = omegaintegr::costh_qtrec();
  //cout << phasespace::costh << "  " << costh_CS << endl;
  /////////////////////

  double m = phasespace::m;
  double qt = phasespace::qt;
  double y = phasespace::y;

  //Dynamic scale
  if (opts.dynamicscale)
    scaleset_(phasespace::m2);
  
  double psfac = 3./8./2./M_PI;
  
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

  int mode = 0;
  /*
  //skip PDF loop in the preconditioning phase
  int Npdf = (iter!=last_iter ? 1 : opts.totpdf);
  for(int ipdf=0; ipdf<Npdf; ipdf++){
      // Set PDF
      //setpdf_(&ipdf);
      //setmellinpdf_(&ipdf);
      //hists_setpdf_(&ipdf);
      //Call the resummation integrand
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
  */

  if (opts.resumcpp)
    f[0]=resint::rint(costh_CS,m,qt,y,mode);
  else
    f[0]=resumm_(costh_CS,m,qt,y,mode);

  if (isnan_ofast(f[0]))
    f[0]=0.;  //avoid nans
  
  //apply switching and jacobian
  f[0] = f[0]*jac*psfac*swtch;
  if (iter==last_iter)
    {
      double wt = weight*f[0];
      hists_fill_(phasespace::p3,phasespace::p4,&wt);
      //hists_AiTest_(pjet,p4cm,&m,&qt,&y,&costh_CS,&phi_lep,&phi,&wt,&lowintHst0);
    }
  
  end_time = clock();
  if (opts.timeprofile)
    cout << setw (3) << "m" << setw(10) << m << setw(4) << "qt" << setw(10) <<  qt
	 << setw (3) << "y" << setw(10) << y << setw(6) << "costh" << setw(10) <<  costh_CS
	 << setw(8) << "result" << setw(10) << f[0]
	 << setw(10) << "tot time" << setw(10) << float( end_time - begin_time ) /  CLOCKS_PER_SEC
	 << endl;
  return 0;
}
