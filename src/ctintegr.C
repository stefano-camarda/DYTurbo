#include "ctintegr.h"

#include "dyres_interface.h"
#include "interface.h"
#include "omegaintegr.h"
#include "phasespace.h"
#include "settings.h"
#include "switch.h"
#include "cubacall.h"
#include "isnan.h"
#include "qtint.h"
#include "ctint.h"
#include "ctmellin.h"
#include "evolnative.h"
#include "pegasus.h"
#include "pmom.h"
#include "ccoeff.h"

#include "KinematicCuts.h"

#include <iostream>
#include <iomanip>
#include <ctime>
#include <math.h>

//Upper limit of integration in qt (in dynnlo it is set to qtcut**2*exp(49) ~ 1e21)
const double ctlimit = 1e10;

using namespace std;


integrand_t ctintegrand(const int &ndim, const double x[], const int &ncomp, double f[],
                        void* userdata, const int &nvec, const int &core,
                        double &weight, const int &iter)
{
  tell_to_grid_we_are_alive();
  //here generate the phase space according to x[], and pass the p vector to countint_

  if (opts.PDFerrors)
    for (int i = 1; i < opts.totpdf; i++)
      f[i] = 0.;
    
  double rct[22];
  for (int i = 0; i < ndim; i++)
    rct[i]=x[i];
  rct[8] = rct[1];

  dofill_.doFill_ = int(iter==last_iter);
  f[0] = countint_(rct,weight,f);

  return 0;
}

integrand_t ctintegrandMC(const int &ndim, const double x[], const int &ncomp, double f[],
			  void* userdata, const int &nvec, const int &core,
			  double &weight, const int &iter)
{
  tell_to_grid_we_are_alive();
  clock_t begin_time, end_time;

  begin_time = clock();

  if (opts.PDFerrors)
    for (int i = 1; i < opts.totpdf; i++)
      f[i] = 0.;
    
  if (opts.fixedorder && phasespace::qtmin > 0)
    {
      f[0]=0.;
      return 0;
    }

  //Jacobian of the change of variables from the unitary hypercube x[6] to the m, y, qt, phi, costh, philep boundaries
  double jac = 1.;
  
  bool status = true;

  double r1[2] = {x[0], x[1]};
  status = phasespace::gen_my(r1, jac, !opts.fixedorder, !opts.fixedorder);  //qtcut = !opts.fixedorder, qtswitching = !opts.fixedorder
  if (!status)
    {
      f[0] = 0.;
      return 0;
    }

  double qt;
  //In the fixed order calculation, integrate from qtcut to infinity
  if (opts.fixedorder)
    {
      phasespace::gen_qt_ctfo(x[2], jac);
      qt = phasespace::qt;
      phasespace::set_qt(0.); //In the fixed order calculation evaluate kinematic cuts with qt=0
    }
  else
    {
      //Generate the boson transverse momentum between the integration boundaries
      double qtcut = max(opts.qtcut,opts.xqtcut*phasespace::m);
      double qtmn = max(qtcut, phasespace::qtmin);
      //phasespace::calcexpy();
      //double cosh2y=pow((phasespace::exppy+phasespace::expmy)*0.5,2);
      //double kinqtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+phasespace::m2,2)/(4*pow(opts.sroot,2)*cosh2y)-phasespace::m2)); //introduced max to avoid neqative argument of sqrt when y=ymax
      //double kinqtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+phasespace::m2,2)/(4*pow(opts.sroot,2))-phasespace::m2)); //introduced max to avoid negative argument of sqrt when y=ymax
      double kinqtlim = ctlimit; //There should not be any kinematic limit on qt, since the counterterm is evaluated with born level kinematic
      double switchqtlim = switching::qtlimit(phasespace::m);
      double qtlim = min(kinqtlim, switchqtlim);
      double qtmx = min(qtlim, phasespace::qtmax);
      if (qtmn >= qtmx)
	{
	  f[0] = 0.;
	  return 0;
	}
      status = phasespace::gen_qt(x[2], jac, qtlim, true);
      if (!status)
	{
	  f[0] = 0.;
	  return 0;
	}
      qt = phasespace::qt;
    }
  jac = jac *2.*qt;

  //generate boson 4-momentum, with m, qt, y, and phi
  phasespace::set_phiV(-M_PI+2.*M_PI*x[5]);

  /////////////////////
  //This function set the RF axes according to the selected qt-recoil prescription
  omegaintegr::genV4p();

  //alternative, generate in naive RF, and calculate costh with the selected qt-recoil prescription
  //phasespace::genRFaxes(phasespace::naive);
  //phasespace::genV4p();
  /////////////////////
  
  double r2[2] = {x[3], x[4]};
  phasespace::gen_costhphi(r2, jac);
  phasespace::calcphilep();
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
  double y = phasespace::y;

  jac = jac/2./M_PI;

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
  if (isnan_ofast(f[0]))
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
	hists_fill_(phasespace::p3, phasespace::p4, &wt);
      }

  end_time = clock();
  if (opts.timeprofile)
    cout << setw (3) << "m" << setw(10) << m << setw(4) << "qt" << setw(10) <<  qt
	 << setw (3) << "y" << setw(10) << y << setw(4) << "costh" << setw(10) <<  costh_CS
	 << setw(8) << "result" << setw(10) << f[0]
	 << setw(10) << "tot time" << setw(10) << float( end_time - begin_time ) /  CLOCKS_PER_SEC
	 << endl;
  return 0;
}

int ctintegrand3d_cubature_v(unsigned ndim, long unsigned npts, const double x[], void *data, unsigned ncomp, double f[])
{
  tell_to_grid_we_are_alive();
  //  cout << "parallel " << npts << endl;
#pragma omp parallel for num_threads(opts.cubacores) copyin(a_param_,scale_,facscale_,qcdcouple_)
  for (unsigned i = 0; i < npts; i++)
    {
      // evaluate the integrand for npts points
      double xi[ndim];
      double fi[ncomp];
      for (unsigned j = 0; j < ndim; j++)
	xi[j] = x[i*ndim + j];

      ctintegrand3d(ndim, xi, ncomp, fi);
      
      for (unsigned k = 0; k < ncomp; ++k)
	f[i*ncomp + k] = fi[k];
    }
  return 0;
}

int ctintegrand3d_cubature(unsigned ndim, const double x[], void *data, unsigned ncomp, double f[])
{
  tell_to_grid_we_are_alive();
  ctintegrand3d(ndim, x, ncomp, f);
  return 0;
}

integrand_t ctintegrand3d(const int &ndim, const double x[], const int &ncomp, double f[])
//Generates the phase space 4 vectors
//Calculates the ct integrand as a function of m, qt, y
//dOmega integration is factorised in the costh moments
//The integration in alpha and beta is performed inside countdy
{
  tell_to_grid_we_are_alive();
  clock_t begin_time, end_time;

  begin_time = clock();

  if (opts.fixedorder && phasespace::qtmin > 0)
    {
      f[0]=0.;
      return 0;
    }
  
  if (opts.PDFerrors)
    for (int i = 1; i < opts.totpdf; i++)
      f[i] = 0.;

  //Jacobian of the change of variables from the unitary hypercube x[3] to the m, y, qt boundaries
  double jac = 1.;
  bool status = true;

  double r2[2] = {x[0], x[1]};
  status = phasespace::gen_my(r2, jac, !opts.fixedorder, !opts.fixedorder);  //qtcut = !opts.fixedorder, qtswitching = !opts.fixedorder
  if (!status)
    {
      f[0] = 0.;
      return 0;
    }

  double qt;
  //In the fixed order calculation, integrate from qtcut to infinity
  if (opts.fixedorder)
    {
      phasespace::gen_qt_ctfo(x[2], jac);
      qt = phasespace::qt;
      phasespace::set_qt(0.); //In the fixed order calculation evaluate kinematic cuts with qt=0
    }
  else
    {
      //Generate the boson transverse momentum between the integration boundaries
      double qtcut = max(opts.qtcut,opts.xqtcut*phasespace::m);
      double qtmn = max(qtcut, phasespace::qtmin);
      //phasespace::calcexpy();
      //double cosh2y=pow((phasespace::exppy+phasespace::expmy)*0.5,2);
      //double kinqtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+phasespace::m2,2)/(4*pow(opts.sroot,2)*cosh2y)-phasespace::m2)); //introduced max to avoid neqative argument of sqrt when y=ymax
      //double kinqtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+phasespace::m2,2)/(4*pow(opts.sroot,2))-phasespace::m2)); //introduced max to avoid negative argument of sqrt when y=ymax
      double kinqtlim = ctlimit; //There should not be any kinematic limit on qt, since the counterterm is evaluated with born level kinematic
      double switchqtlim = switching::qtlimit(phasespace::m);
      double qtlim = min(kinqtlim, switchqtlim);
      double qtmx = min(qtlim, phasespace::qtmax);
      if (qtmn >= qtmx)
	{
	  f[0] = 0.;
	  return 0;
	}
      status = phasespace::gen_qt(x[2], jac, qtlim, true);
      if (!status)
	{
	  f[0] = 0.;
	  return 0;
	}
      qt = phasespace::qt;
    }
  jac = jac *2.*qt;

  phasespace::set_phiV(0.);
  
  double m = phasespace::m;
  double y = phasespace::y;

  if (opts.ctcpp)
    qtint::calc(m,qt,0,1);
  
  //generate boson 4-momentum, with m, qt, y and phi=0
  omegaintegr::genV4p();

  //switching is applied inside countdy
  //double swtch = switching::swtch(qt, m);

  //In this point of phase space (m, qt, y) the costh integration is performed by 
  //calculating the 0, 1 and 2 moments of costh
  //that are the integrals int(dcosth dphi1 dphi2), int(costh dcosth dphi1 dphi2), and int(costh^2 dcosth dphi1 dphi2) convoluted with cuts
  //Then the epxressions 1, costh and costh^2 in sigmaij are substituted by these costh moments
  double costh = 0;
  int mode = 1;
  dofill_.doFill_ = 1;

  //No need to check the switching, since the phase space is generated up to the switching qt limit
  //if (swtch < switching::cutoff*switching::tolerance)
  //{
  //  f[0]=0.;
  //  return 0;
  //}

  //evaluate the fixed order expansion of the resummed cross section
  if (opts.ctcpp)
    ctint::calc(costh,m,qt,y,mode,f);
  else
    f[0]=ctint_(costh,m,qt,y,mode,f);
  
  //avoid nans
  if (isnan_ofast(f[0]))
    {
      cout << "nan in ctintegrand3d" << endl;
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
    cout << setw (3) << "m" << setw(10) << m << setw(4) << "qt" << setw(10) << qt << setw(4) << "y" << setw(10) <<  y
	 << setw(8) << "result" << setw(10) << f[0]
	 << setw(10) << "tot time" << setw(10) << float( end_time - begin_time ) /  CLOCKS_PER_SEC
	 << endl;
  return 0;
}

int ctintegrand2d_cubature_v(unsigned ndim, long unsigned npts, const double x[], void *data, unsigned ncomp, double f[])
{
  tell_to_grid_we_are_alive();
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
  tell_to_grid_we_are_alive();
  ctintegrand2d(ndim, x, ncomp, f);
  return 0;
}

integrand_t ctintegrand2d(const int &ndim, const double x[], const int &ncomp, double f[])
//Calculates the ct integrand as a function of m, y
//The integration in qt is factorised in LL1, LL2, LL3 and LL4 large logs
//dOmega integration is factorised in the costh moments
//The integration in alpha and beta is performed inside countdy
{
  tell_to_grid_we_are_alive();
  clock_t begin_time, end_time;

  begin_time = clock();

  if (opts.fixedorder && phasespace::qtmin > 0)
    {
      f[0]=0.;
      return 0;
    }
  
  if (opts.PDFerrors)
    for (int i = 1; i < opts.totpdf; i++)
      f[i] = 0.;

  //Jacobian of the change of variables from the unitary hypercube x[3] to the m, y, qt boundaries
  double jac = 1.;
  bool status = true;

  double r2[2] = {x[0], x[1]};
  status = phasespace::gen_my(r2, jac, !opts.fixedorder, !opts.fixedorder);  //qtcut = !opts.fixedorder, qtswitching = !opts.fixedorder
  if (!status)
    {
      f[0] = 0.;
      return 0;
    }

  double qtcut = max(opts.qtcut,opts.xqtcut*phasespace::m);
  double qtmn, qtmx;
  //In the fixed order calculation, integrate from qtcut to infinity
  if (opts.fixedorder)
    {
      qtmn = qtcut;
      qtmx = ctlimit;
      phasespace::set_qt(0.); //In the fixed order calculation evaluate kinematic cuts with qt=0
    }
  else
    {
      //Generate the boson transverse momentum between the integration boundaries
      qtmn = max(qtcut, phasespace::qtmin);
      //phasespace::calcexpy();
      //double cosh2y=pow((phasespace::exppy+phasespace::expmy)*0.5,2);
      //double kinqtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+phasespace::m2,2)/(4*pow(opts.sroot,2)*cosh2y)-phasespace::m2)); //introduced max to avoid neqative argument of sqrt when y=ymax
      //double kinqtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+phasespace::m2,2)/(4*pow(opts.sroot,2))-phasespace::m2)); //introduced max to avoid negative argument of sqrt when y=ymax
      double kinqtlim = ctlimit; //There should not be any kinematic limit on qt, since the counterterm is evaluated with born level kinematic
      double switchqtlim = switching::qtlimit(phasespace::m);
      double qtlim = min(kinqtlim, switchqtlim);
      qtmx = min(qtlim, phasespace::qtmax);
      if (qtmn >= qtmx)
	{
	  f[0] = 0.;
	  return 0;
	}
      phasespace::set_qt((qtmn+qtmx)/2.); //Is this actually needed?
    }
  phasespace::set_phiV(0.);
  double m = phasespace::m;
  double y = phasespace::y;

  clock_t qtbt, qtet;
  qtbt = clock();
  if (opts.ctcpp)
    qtint::calc(m,qtmn,qtmx,2);
  else
    ctqtint_(m,y,qtmn,qtmx);
  qtet = clock();
  
  //generate boson 4-momentum, with m, qt, y and phi=0 (not actually needed)
  omegaintegr::genV4p();

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
  if (opts.ctcpp)
    ctint::calc(costh,m,qt,y,mode,f);
  else
    f[0]=ctint_(costh,m,qt,y,mode,f);
  clock_t cet = clock();

  //avoid nans
  if (isnan_ofast(f[0]))
    {
      cout << "nan in ctintegr 2d " << endl;
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

int ctintegrand1d_cubature_v(unsigned ndim, long unsigned npts, const double x[], void *data, unsigned ncomp, double f[])
{
  tell_to_grid_we_are_alive();
  //  cout << "parallel " << npts << endl;
#pragma omp parallel for num_threads(opts.cubacores) \
  copyin(a_param_,scale_,facscale_,qcdcouple_,\
	 mellinint::wn,mellinint::Np,mellinint::Nm,mellinint::wn_1,mellinint::Np_1,mellinint::Nm_1,mellinint::wn_2,mellinint::Np_2,mellinint::Nm_2,\
	 evolnative::ans,evolnative::am,evolnative::ap,evolnative::al,evolnative::be,evolnative::ab,\
	 evolnative::rmin,evolnative::rplus,evolnative::rqq,evolnative::rqg,evolnative::rgq,evolnative::rgg,\
	 evolnative::rmmqq,evolnative::rmmqg,evolnative::rmmgq,evolnative::rmmgg,evolnative::rmpqq,evolnative::rmpqg,evolnative::rmpgq,evolnative::rmpgg,\
	 evolnative::rpmqq,evolnative::rpmqg,evolnative::rpmgq,evolnative::rpmgg,evolnative::rppqq,evolnative::rppqg,evolnative::rppgq,evolnative::rppgg,\
	 pegasus::gli,pegasus::vai,pegasus::m3i,pegasus::m8i,pegasus::m15i,pegasus::m24i,pegasus::sgi,pegasus::p3i,pegasus::p8i,pegasus::p15i,pegasus::p24i,\
	 order_,painp_,hfpainp_,pacthr_,pabthr_,patthr_,moms_,hsums_,pns0_,psg0_,pns1_,psg1_,pns2_,psg2_,asg2_,lsg_,spsums_,u1sg_,r1sg_,u1hsg_,u2sg_,r2sg_,u2hsg_,u2ns_, \
	 ccoeff::C1qg,ccoeff::C1qq,ccoeff::C1qqb,ccoeff::C1qqp,ccoeff::C1qqbp,ccoeff::C2qg,ccoeff::C2qq,ccoeff::C2qqb,ccoeff::C2qqp,ccoeff::C2qqbp,ccoeff::C3qg,ccoeff::C3qq,ccoeff::C3qqb,ccoeff::C3qqp,ccoeff::C3qqbp,\
	 ccoeff::C1qg_1,ccoeff::C1qq_1,ccoeff::C1qqb_1,ccoeff::C1qqp_1,ccoeff::C1qqbp_1,ccoeff::C2qg_1,ccoeff::C2qq_1,ccoeff::C2qqb_1,ccoeff::C2qqp_1,ccoeff::C2qqbp_1,ccoeff::C3qg_1,ccoeff::C3qq_1,ccoeff::C3qqb_1,ccoeff::C3qqp_1,ccoeff::C3qqbp_1,\
	 ccoeff::C1qg_2,ccoeff::C1qq_2,ccoeff::C1qqb_2,ccoeff::C1qqp_2,ccoeff::C1qqbp_2,ccoeff::C2qg_2,ccoeff::C2qq_2,ccoeff::C2qqb_2,ccoeff::C2qqp_2,ccoeff::C2qqbp_2,ccoeff::C3qg_2,ccoeff::C3qq_2,ccoeff::C3qqb_2,ccoeff::C3qqp_2,ccoeff::C3qqbp_2,\
	 pmom::gamma1qq,pmom::gamma1qqb,pmom::gamma1qqp,pmom::gamma1qqbp,pmom::gamma1qg,pmom::gamma1gq,pmom::gamma1gg,pmom::gamma2qq,pmom::gamma2qqb,pmom::gamma2qqp,pmom::gamma2qqbp,pmom::gamma2qg,pmom::gamma2gq,pmom::gamma2gg,pmom::gamma3qq,pmom::gamma3qqb,pmom::gamma3qqp,pmom::gamma3qqbp,pmom::gamma3qg,pmom::gamma3gq,pmom::gamma3gg,\
	 pmom::gamma1qq_1,pmom::gamma1qqb_1,pmom::gamma1qqp_1,pmom::gamma1qqbp_1,pmom::gamma1qg_1,pmom::gamma1gq_1,pmom::gamma1gg_1,pmom::gamma2qq_1,pmom::gamma2qqb_1,pmom::gamma2qqp_1,pmom::gamma2qqbp_1,pmom::gamma2qg_1,pmom::gamma2gq_1,pmom::gamma2gg_1,pmom::gamma3qq_1,pmom::gamma3qqb_1,pmom::gamma3qqp_1,pmom::gamma3qqbp_1,pmom::gamma3qg_1,pmom::gamma3gq_1,pmom::gamma3gg_1,\
	 pmom::gamma1qq_2,pmom::gamma1qqb_2,pmom::gamma1qqp_2,pmom::gamma1qqbp_2,pmom::gamma1qg_2,pmom::gamma1gq_2,pmom::gamma1gg_2,pmom::gamma2qq_2,pmom::gamma2qqb_2,pmom::gamma2qqp_2,pmom::gamma2qqbp_2,pmom::gamma2qg_2,pmom::gamma2gq_2,pmom::gamma2gg_2,pmom::gamma3qq_2,pmom::gamma3qqb_2,pmom::gamma3qqp_2,pmom::gamma3qqbp_2,pmom::gamma3qg_2,pmom::gamma3gq_2,pmom::gamma3gg_2)
  //copyin(a_param_,scale_,facscale_,qcdcouple_,evolnative::UVP,evolnative::DVP,evolnative::USP,evolnative::DSP,evolnative::SSP,evolnative::GLP,evolnative::CHP,evolnative::BOP,evolnative::SVP,evolnative::CVP,evolnative::BVP)
  for (unsigned i = 0; i < npts; i++)
    {
      // evaluate the integrand for npts points
      double xi[ndim];
      double fi[ncomp];
      for (unsigned j = 0; j < ndim; j++)
	xi[j] = x[i*ndim + j];

      ctintegrand1d(ndim, xi, ncomp, fi);
      
      for (unsigned k = 0; k < ncomp; ++k)
	f[i*ncomp + k] = fi[k];
    }
  return 0;
}

int ctintegrand1d_cubature(unsigned ndim, const double x[], void *data, unsigned ncomp, double f[])
{
  tell_to_grid_we_are_alive();
  ctintegrand1d(ndim, x, ncomp, f);
  return 0;
}

integrand_t ctintegrand1d(const int &ndim, const double x[], const int &ncomp, double f[])
//Calculates the ct integrand as a function of m,
//The integration in qt is factorised in LL1, LL2, LL3 and LL4 large logs
//dOmega integration is factorised in the costh moments
//The y integration is achieved by means of a single Mellin inversion
{
  tell_to_grid_we_are_alive();
  clock_t begin_time, end_time;

  begin_time = clock();
  if (opts.fixedorder && phasespace::qtmin > 0)
    {
      f[0]=0.;
      return 0;
    }
  
  if (opts.PDFerrors)
    for (int i = 1; i < opts.totpdf; i++)
      f[i] = 0.;

  //Jacobian of the change of variables from the unitary segment x[0] to the m boundaries
  double jac = 1.;
  bool status = true;

  //This integration currently works only for the full y range, so the mass limit is sqrt(s)
  double mlim = opts.sroot;
  double r1 = {x[0]};
  status = phasespace::gen_m(r1, jac, mlim, !opts.fixedorder, !opts.fixedorder); //qtcut = !opts.fixedorder, qtswitching = !opts.fixedorder
  if (!status)
    {
      f[0] = 0.;
      return 0;
    }
  
  double qtcut = max(opts.qtcut,opts.xqtcut*phasespace::m);
  double qtmn, qtmx;
  //In the fixed order calculation, integrate from qtcut to infinity
  if (opts.fixedorder)
    {
      qtmn = qtcut;
      qtmx = ctlimit;
      phasespace::set_qt(0.); //In the fixed order calculation evaluate kinematic cuts with qt=0
    }
  else
    {
      //Generate the boson transverse momentum between the integration boundaries
      qtmn = max(qtcut, phasespace::qtmin);
      //phasespace::calcexpy();
      //double cosh2y=pow((phasespace::exppy+phasespace::expmy)*0.5,2);
      //double kinqtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+phasespace::m2,2)/(4*pow(opts.sroot,2)*cosh2y)-phasespace::m2)); //introduced max to avoid neqative argument of sqrt when y=ymax
      //double kinqtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+phasespace::m2,2)/(4*pow(opts.sroot,2))-phasespace::m2)); //introduced max to avoid negative argument of sqrt when y=ymax
      double kinqtlim = ctlimit; //There should not be any kinematic limit on qt, since the counterterm is evaluated with born level kinematic
      double switchqtlim = switching::qtlimit(phasespace::m);
      double qtlim = min(kinqtlim, switchqtlim);
      qtmx = min(qtlim, phasespace::qtmax);
      if (qtmn >= qtmx)
	{
	  f[0] = 0.;
	  return 0;
	}
      phasespace::set_qt((qtmn+qtmx)/2.); //Is this actually needed?
    }
  phasespace::set_phiV(0.);
  phasespace::set_y(0.);
  double m = phasespace::m;
  //double y = phasespace::y;

  clock_t qtbt, qtet;
  qtbt = clock();
  qtint::calc(m,qtmn,qtmx,2);
  qtet = clock();
  
  //generate boson 4-momentum, with m, qt, y and phi=0 (not actually needed)
  omegaintegr::genV4p();

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
  ctmellin::calc(m,f);
  clock_t cet = clock();

  //avoid nans
  if (isnan_ofast(f[0]))
    {
      cout << "nan in ctintegr 1d " << endl;
      f[0]=0.;
      if (opts.PDFerrors)
	for (int i = 1; i < opts.totpdf; i++)
	  f[i] = 0.;
      return 0;
    }

  //  cout << "counterterm " << f[0] << " m " << m << " jac " << jac << endl;
	   
  f[0] = f[0]*jac;
  if (opts.PDFerrors)
    for (int i = 1; i < opts.totpdf; i++)
      f[i] = f[i]*jac;

  end_time = clock();
  if (opts.timeprofile)
    cout << setw (3) << "m" << setw(10) << m 
	 << setw(8) << "result" << setw(10) << f[0]
	 << setw(10) << "tot time" << setw(10) << float( end_time - begin_time ) /  CLOCKS_PER_SEC
      	 << setw(10) << "qtint"  << setw(10) << float( qtet - qtbt ) /  CLOCKS_PER_SEC
	 << setw(10) << "ctint"  << setw(10) << float( cet - cbt ) /  CLOCKS_PER_SEC
	 << endl;
  return 0;
}
