#include "resintegr.h"
#include "omegaintegr.h"
#include "phasespace.h"
#include "settings.h"
#include "interface.h"
#include "dyres_interface.h"
#include "switch.h"
#include "resint.h"
#include "rapint.h"
#include "cubacall.h"
#include "isnan.h"
#include "scales.h"
#include "evolnative.h"
#include "pegasus.h"
#include "pmom.h"
#include "ccoeff.h"

#include <math.h>
#include <iomanip>
#include <iostream>
#include <omp.h>

#include "KinematicCuts.h"

int resintegrand1d_cubature_v(unsigned ndim, long unsigned npts, const double x[], void *data, unsigned ncomp, double f[])
{
#pragma omp parallel for num_threads(opts.cubacores) \
  copyin(a_param_,rlogs_,resint::a,resint::loga,resint::rloga,\
	 mellinint::wn,mellinint::Np,mellinint::Nm,mellinint::wn_1,mellinint::Np_1,mellinint::Nm_1,mellinint::wn_2,mellinint::Np_2,mellinint::Nm_2,\
	 evolnative::ans,evolnative::am,evolnative::ap,evolnative::al,evolnative::be,evolnative::ab,\
	 evolnative::rmin,evolnative::rplus,evolnative::rqq,evolnative::rqg,evolnative::rgq,evolnative::rgg,\
	 evolnative::rmmqq,evolnative::rmmqg,evolnative::rmmgq,evolnative::rmmgg,evolnative::rmpqq,evolnative::rmpqg,evolnative::rmpgq,evolnative::rmpgg,\
	 evolnative::rpmqq,evolnative::rpmqg,evolnative::rpmgq,evolnative::rpmgg,evolnative::rppqq,evolnative::rppqg,evolnative::rppgq,evolnative::rppgg,\
	 pegasus::gli,pegasus::vai,pegasus::m3i,pegasus::m8i,pegasus::m15i,pegasus::m24i,pegasus::sgi,pegasus::p3i,pegasus::p8i,pegasus::p15i,pegasus::p24i,\
	 painp_,hfpainp_,pacthr_,pabthr_,patthr_,moms_,hsums_,pns0_,psg0_,pns1_,psg1_,pns2_,psg2_,asg2_,lsg_,spsums_,u1sg_,r1sg_,u1hsg_,u2sg_,r2sg_,u2hsg_,u2ns_, \
	 ccoeff::C1qg,ccoeff::C1qq,ccoeff::C1qqb,ccoeff::C1qqp,ccoeff::C1qqbp,ccoeff::C2qg,ccoeff::C2qq,ccoeff::C2qqb,ccoeff::C2qqp,ccoeff::C2qqbp,ccoeff::C3qg,ccoeff::C3qq,ccoeff::C3qqb,ccoeff::C3qqp,ccoeff::C3qqbp,\
	 ccoeff::C1qg_1,ccoeff::C1qq_1,ccoeff::C1qqb_1,ccoeff::C1qqp_1,ccoeff::C1qqbp_1,ccoeff::C2qg_1,ccoeff::C2qq_1,ccoeff::C2qqb_1,ccoeff::C2qqp_1,ccoeff::C2qqbp_1,ccoeff::C3qg_1,ccoeff::C3qq_1,ccoeff::C3qqb_1,ccoeff::C3qqp_1,ccoeff::C3qqbp_1,\
	 ccoeff::C1qg_2,ccoeff::C1qq_2,ccoeff::C1qqb_2,ccoeff::C1qqp_2,ccoeff::C1qqbp_2,ccoeff::C2qg_2,ccoeff::C2qq_2,ccoeff::C2qqb_2,ccoeff::C2qqp_2,ccoeff::C2qqbp_2,ccoeff::C3qg_2,ccoeff::C3qq_2,ccoeff::C3qqb_2,ccoeff::C3qqp_2,ccoeff::C3qqbp_2,\
	 pmom::gamma1qq,pmom::gamma1qqb,pmom::gamma1qqp,pmom::gamma1qqbp,pmom::gamma1qg,pmom::gamma1gq,pmom::gamma1gg,pmom::gamma2qq,pmom::gamma2qqb,pmom::gamma2qqp,pmom::gamma2qqbp,pmom::gamma2qg,pmom::gamma2gq,pmom::gamma2gg,pmom::gamma3qq,pmom::gamma3qqb,pmom::gamma3qqp,pmom::gamma3qqbp,pmom::gamma3qg,pmom::gamma3gq,pmom::gamma3gg,\
	 pmom::gamma1qq_1,pmom::gamma1qqb_1,pmom::gamma1qqp_1,pmom::gamma1qqbp_1,pmom::gamma1qg_1,pmom::gamma1gq_1,pmom::gamma1gg_1,pmom::gamma2qq_1,pmom::gamma2qqb_1,pmom::gamma2qqp_1,pmom::gamma2qqbp_1,pmom::gamma2qg_1,pmom::gamma2gq_1,pmom::gamma2gg_1,pmom::gamma3qq_1,pmom::gamma3qqb_1,pmom::gamma3qqp_1,pmom::gamma3qqbp_1,pmom::gamma3qg_1,pmom::gamma3gq_1,pmom::gamma3gg_1,\
	 pmom::gamma1qq_2,pmom::gamma1qqb_2,pmom::gamma1qqp_2,pmom::gamma1qqbp_2,pmom::gamma1qg_2,pmom::gamma1gq_2,pmom::gamma1gg_2,pmom::gamma2qq_2,pmom::gamma2qqb_2,pmom::gamma2qqp_2,pmom::gamma2qqbp_2,pmom::gamma2qg_2,pmom::gamma2gq_2,pmom::gamma2gg_2,pmom::gamma3qq_2,pmom::gamma3qqb_2,pmom::gamma3qqp_2,pmom::gamma3qqbp_2,pmom::gamma3qg_2,pmom::gamma3gq_2,pmom::gamma3gg_2)
  //  evolnative::UVP,evolnative::DVP,evolnative::USP,evolnative::DSP,evolnative::SSP,evolnative::GLP,evolnative::CHP,evolnative::BOP,evolnative::SVP,evolnative::CVP,evolnative::BVP,evolnative::SIP,evolnative::NS3P,evolnative::NS8P,evolnative::NS15P,evolnative::NS24P,evolnative::NS35P)
  
  for (unsigned i = 0; i < npts; i++)
    {
      // evaluate the integrand for npts points
      double xi[ndim];
      double fi[ncomp];
      for (unsigned j = 0; j < ndim; j++)
	xi[j] = x[i*ndim + j];

      resintegrand1d(ndim, xi, ncomp, fi);
      
      for (unsigned k = 0; k < ncomp; ++k)
	f[i*ncomp + k] = fi[k];
    }
  return 0;
}

int resintegrand1d_cubature(unsigned ndim, const double x[], void *data, unsigned ncomp, double f[])
{
  resintegrand1d(ndim, x, ncomp, f);
  return 0;
}

//This integration works only without cuts (!opts.makecuts)
integrand_t resintegrand1d(const int &ndim, const double x[], const int &ncomp, double f[])
{
  tell_to_grid_we_are_alive();
  clock_t begin_time, end_time;

  begin_time = clock();

  //Jacobian of the change of variables from the segment [0,1] to qt boundaries
  double jac = 1.;
  bool status = true;

  //Generate the boson invariant mass between the integration boundaries
  //kinematic limit from the relation x(1,2) < 1 ==> exp(fabs(y)) < sqrt(s)/m
  //There is also a lower limit from the switching: m > qtmin/opts.dampk
  double absymn;
  if (phasespace::ymin*phasespace::ymax <= 0.)
    absymn = 0.;
  else
    absymn = min(fabs(phasespace::ymin),fabs(phasespace::ymax));
  double mlim = opts.sroot/exp(absymn);
  double r1 = {x[0]};
  status = phasespace::gen_m(r1, jac, mlim, false, true);
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
  //if (opts.dynamicscale)
  //scaleset_(phasespace::m2);
  scales::set(phasespace::m);
  scales::mcfm();

  clock_t ybt, yet;
  ybt = clock();
//  if (!opts.mellin1d)
//    {
//      rapint::allocate();
//      rapint::integrate(ymn,ymx,phasespace::m);
//    }
  yet = clock();
  
  //switching cannot be applied because we integrate analitically in pt
  //double swtch = switching::swtch(phasespace::qt, phasespace::m);
  
  //Call the resummed cross section
  double costh = 0;
  double y = 0.5*(ymx+ymn); //needed to get the correct IFIT index of the approximate PDF
  phasespace::set_y(y); //needed to get the correct mlogz factor in mellinint::update()
  double qt = 0.5*(phasespace::qtmax+phasespace::qtmin); //not needed
  phasespace::set_qt(qt); //needed to get the correct mlogz factor in mellinint::update()
  int mode = 3;

  clock_t rbt = clock();
  double ff[2];
  resint::rint(costh,phasespace::m,qt,y,mode,ff);
  clock_t ret = clock();

  f[0] = ff[0]/(8./3.);
  if (opts.helicity >= 0)
    f[1] = ff[1]/(8./3.);
  
  //  cout << "resummation " << f[0] << " m " << phasespace::m  << " jac " << jac << endl;
  
  if (isnan_ofast(f[0]))
    {
      f[0] = 0.;  //avoid nans
      if (opts.helicity >= 0)
	f[1]=0.;
    }

  f[0] = f[0]*jac;
  if (opts.helicity >= 0)
    f[1] = f[1]*jac;

  end_time = clock();
  if (opts.timeprofile)
    cout << setw (3) << "m" << setw(10) << phasespace::m << setw(4) << "qt" << setw(10) <<  qt
	 << setw(8) << "result" << setw(10) << f[0]
	 << setw(10) << "tot time" << setw(10) << float( end_time - begin_time ) /  CLOCKS_PER_SEC
	 << setw(10) << "rapint"  << setw(10) << float( yet - ybt ) /  CLOCKS_PER_SEC
	 << setw(10) << "resumm"  << setw(10) << float( ret - rbt ) /  CLOCKS_PER_SEC
	 << endl;

  if (!opts.mellin1d)
      rapint::free();

  return 0;
}

int resintegrand2d_cubature_v(unsigned ndim, long unsigned npts, const double x[], void *data, unsigned ncomp, double f[])
{
  tell_to_grid_we_are_alive();
  //The current issue with openmp parallelisation, is that the resummation integrand has many
  //global variables, and the use of these variables has a race
  
  //  cout << "parallel " << npts << endl;

#pragma omp parallel for num_threads(opts.cubacores) \
  copyin(a_param_,rlogs_,resint::a,resint::loga,resint::rloga,\
	 mellinint::wn,mellinint::Np,mellinint::Nm,mellinint::wn_1,mellinint::Np_1,mellinint::Nm_1,mellinint::wn_2,mellinint::Np_2,mellinint::Nm_2,\
	 evolnative::ans,evolnative::am,evolnative::ap,evolnative::al,evolnative::be,evolnative::ab,\
	 evolnative::rmin,evolnative::rplus,evolnative::rqq,evolnative::rqg,evolnative::rgq,evolnative::rgg,\
	 evolnative::rmmqq,evolnative::rmmqg,evolnative::rmmgq,evolnative::rmmgg,evolnative::rmpqq,evolnative::rmpqg,evolnative::rmpgq,evolnative::rmpgg,\
	 evolnative::rpmqq,evolnative::rpmqg,evolnative::rpmgq,evolnative::rpmgg,evolnative::rppqq,evolnative::rppqg,evolnative::rppgq,evolnative::rppgg,\
	 pegasus::gli,pegasus::vai,pegasus::m3i,pegasus::m8i,pegasus::m15i,pegasus::m24i,pegasus::sgi,pegasus::p3i,pegasus::p8i,pegasus::p15i,pegasus::p24i,\
	 painp_,hfpainp_,pacthr_,pabthr_,patthr_,moms_,hsums_,pns0_,psg0_,pns1_,psg1_,pns2_,psg2_,asg2_,lsg_,spsums_,u1sg_,r1sg_,u1hsg_,u2sg_,r2sg_,u2hsg_,u2ns_, \
	 ccoeff::C1qg,ccoeff::C1qq,ccoeff::C1qqb,ccoeff::C1qqp,ccoeff::C1qqbp,ccoeff::C2qg,ccoeff::C2qq,ccoeff::C2qqb,ccoeff::C2qqp,ccoeff::C2qqbp,ccoeff::C3qg,ccoeff::C3qq,ccoeff::C3qqb,ccoeff::C3qqp,ccoeff::C3qqbp,\
	 ccoeff::C1qg_1,ccoeff::C1qq_1,ccoeff::C1qqb_1,ccoeff::C1qqp_1,ccoeff::C1qqbp_1,ccoeff::C2qg_1,ccoeff::C2qq_1,ccoeff::C2qqb_1,ccoeff::C2qqp_1,ccoeff::C2qqbp_1,ccoeff::C3qg_1,ccoeff::C3qq_1,ccoeff::C3qqb_1,ccoeff::C3qqp_1,ccoeff::C3qqbp_1,\
	 ccoeff::C1qg_2,ccoeff::C1qq_2,ccoeff::C1qqb_2,ccoeff::C1qqp_2,ccoeff::C1qqbp_2,ccoeff::C2qg_2,ccoeff::C2qq_2,ccoeff::C2qqb_2,ccoeff::C2qqp_2,ccoeff::C2qqbp_2,ccoeff::C3qg_2,ccoeff::C3qq_2,ccoeff::C3qqb_2,ccoeff::C3qqp_2,ccoeff::C3qqbp_2,\
	 pmom::gamma1qq,pmom::gamma1qqb,pmom::gamma1qqp,pmom::gamma1qqbp,pmom::gamma1qg,pmom::gamma1gq,pmom::gamma1gg,pmom::gamma2qq,pmom::gamma2qqb,pmom::gamma2qqp,pmom::gamma2qqbp,pmom::gamma2qg,pmom::gamma2gq,pmom::gamma2gg,pmom::gamma3qq,pmom::gamma3qqb,pmom::gamma3qqp,pmom::gamma3qqbp,pmom::gamma3qg,pmom::gamma3gq,pmom::gamma3gg,\
	 pmom::gamma1qq_1,pmom::gamma1qqb_1,pmom::gamma1qqp_1,pmom::gamma1qqbp_1,pmom::gamma1qg_1,pmom::gamma1gq_1,pmom::gamma1gg_1,pmom::gamma2qq_1,pmom::gamma2qqb_1,pmom::gamma2qqp_1,pmom::gamma2qqbp_1,pmom::gamma2qg_1,pmom::gamma2gq_1,pmom::gamma2gg_1,pmom::gamma3qq_1,pmom::gamma3qqb_1,pmom::gamma3qqp_1,pmom::gamma3qqbp_1,pmom::gamma3qg_1,pmom::gamma3gq_1,pmom::gamma3gg_1,\
	 pmom::gamma1qq_2,pmom::gamma1qqb_2,pmom::gamma1qqp_2,pmom::gamma1qqbp_2,pmom::gamma1qg_2,pmom::gamma1gq_2,pmom::gamma1gg_2,pmom::gamma2qq_2,pmom::gamma2qqb_2,pmom::gamma2qqp_2,pmom::gamma2qqbp_2,pmom::gamma2qg_2,pmom::gamma2gq_2,pmom::gamma2gg_2,pmom::gamma3qq_2,pmom::gamma3qqb_2,pmom::gamma3qqp_2,pmom::gamma3qqbp_2,pmom::gamma3qg_2,pmom::gamma3gq_2,pmom::gamma3gg_2)
  //copyin(a_param_,rlogs_,resint::a,resint::loga,resint::rloga,painp_,hfpainp_,mellinint::wn,mellinint::Np,mellinint::Nm,mellinint::wn_1,mellinint::Np_1,mellinint::Nm_1,mellinint::wn_2,mellinint::Np_2,mellinint::Nm_2,evolnative::ans,evolnative::am,evolnative::ap,evolnative::al,evolnative::be,evolnative::ab,evolnative::rmin,evolnative::rplus,evolnative::rqq,evolnative::rqg,evolnative::rgq,evolnative::rgg,evolnative::rmmqq,evolnative::rmmqg,evolnative::rmmgq,evolnative::rmmgg,evolnative::rmpqq,evolnative::rmpqg,evolnative::rmpgq,evolnative::rmpgg,evolnative::rpmqq,evolnative::rpmqg,evolnative::rpmgq,evolnative::rpmgg,evolnative::rppqq,evolnative::rppqg,evolnative::rppgq,evolnative::rppgg)
  //#pragma omp parallel for num_threads(opts.cubacores) copyin(a_param_,rlogs_,resint::a,resint::loga,resint::rloga,evolnative::UVP,evolnative::DVP,evolnative::USP,evolnative::DSP,evolnative::SSP,evolnative::GLP,evolnative::CHP,evolnative::BOP,evolnative::SIP,evolnative::NS3P,evolnative::NS8P,evolnative::NS15P,evolnative::NS24P,evolnative::NS35P,painp_,hfpainp_)

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
  //if (opts.dynamicscale)
  //scaleset_(phasespace::m2);
  scales::set(phasespace::m);
  scales::mcfm();

  //Perform quadrature rule integration in rapidity and semi-analitical costh, phi_lep integration
  int nocuts = !opts.makecuts;

  clock_t ybt, yet;
  ybt = clock();
  //there is a potential issue here, when lepton cuts are applied
  //the rapidity dependent exponential are cached assuming integration between ymin and ymax
  //for consistency, has to keep the integration between ymin and ymax
  //if (opts.makecuts)
  //  {
  //    if (opts.resumcpp)
  //	{
  //	  //C++ resum
  //	  rapint::allocate();
  //	  rapint::integrate(phasespace::ymin,phasespace::ymax,phasespace::m);
  //	  //end C++ resum
  //	}
  //    else
  //	rapintegrals_(phasespace::ymin,phasespace::ymax,phasespace::m,nocuts);
  //  }
  //else
  //  {
  //    if (opts.resumcpp)
  //	{
  //	  //C++ rewritten resum
  //	  if (!opts.mellin1d)
  //	    {
  //	      rapint::allocate();
  //	      rapint::integrate(ymn,ymx,phasespace::m);
  //	    }
  //	  //end C++ resum
  //	}
  //    else
  //	rapintegrals_(ymn,ymx,phasespace::m,nocuts);
  //  }
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

  double ff[2];
  if (opts.resumcpp)
  //f[0]=resint::rint(costh,phasespace::m,phasespace::qt,y,mode)/(8./3.);
    resint::rint(costh,phasespace::m,phasespace::qt,y,mode,ff);
  else
    f[0]=resumm_(costh,phasespace::m,phasespace::qt,y,mode)/(8./3.);
  clock_t ret = clock();

  f[0] = ff[0]/(8./3.);
  if (opts.helicity >= 0)
    f[1] = ff[1]/(8./3.);

  if (isnan_ofast(f[0]))
    {
      f[0] = 0.;  //avoid nans
      if (opts.helicity >= 0)
	f[1]=0.;
    }

  f[0] = f[0]*jac*swtch;
  if (opts.helicity >= 0)
    f[1] = f[1]*jac*swtch;

  end_time = clock();
  if (opts.timeprofile)
    cout << setw (3) << "m" << setw(10) << phasespace::m << setw(4) << "qt" << setw(10) <<  phasespace::qt
	 << setw(8) << "result" << setw(10) << f[0]
	 << setw(10) << "tot time" << setw(10) << float( end_time - begin_time ) /  CLOCKS_PER_SEC
	 << setw(10) << "rapint"  << setw(10) << float( yet - ybt ) /  CLOCKS_PER_SEC
	 << setw(10) << "resumm"  << setw(10) << float( ret - rbt ) /  CLOCKS_PER_SEC
	 << endl;

  if (opts.resumcpp)
    if (!opts.mellin1d)
      rapint::free();

  return 0;
}
/********************************/
int resintegrand2d_my_cubature_v(unsigned ndim, long unsigned npts, const double x[], void *data, unsigned ncomp, double f[])
{
  tell_to_grid_we_are_alive();
#pragma omp parallel for num_threads(opts.cubacores) copyin(a_param_,rlogs_,resint::a,resint::loga,resint::rloga,evolnative::UVP,evolnative::DVP,evolnative::USP,evolnative::DSP,evolnative::SSP,evolnative::GLP,evolnative::CHP,evolnative::BOP,evolnative::SIP,evolnative::NS3P,evolnative::NS8P,evolnative::NS15P,evolnative::NS24P,evolnative::NS35P,painp_,hfpainp_)

  for (unsigned i = 0; i < npts; i++)
    {
      // evaluate the integrand for npts points
      double xi[ndim];
      double fi[ncomp];
      for (unsigned j = 0; j < ndim; j++)
	xi[j] = x[i*ndim + j];

      resintegrand2d_my(ndim, xi, ncomp, fi);
      
      for (unsigned k = 0; k < ncomp; ++k)
	f[i*ncomp + k] = fi[k];
    }
  return 0;
}

int resintegrand2d_my_cubature(unsigned ndim, const double x[], void *data, unsigned ncomp, double f[])
{
  tell_to_grid_we_are_alive();
  resintegrand2d_my(ndim, x, ncomp, f);
  return 0;
}

integrand_t resintegrand2d_my(const int &ndim, const double x[], const int &ncomp, double f[])
{
  tell_to_grid_we_are_alive();
  clock_t begin_time, end_time;

  begin_time = clock();

  //Jacobian of the change of variables from the unitary hypercube x[3] to the m, y, qt boundaries
  double jac = 1.;
  bool status = true;

  double r1[2] = {x[0], x[1]};
  status = phasespace::gen_my(r1, jac, false, true);
  if (!status)
    {
      f[0] = 0.;
      return 0;
    }

  phasespace::set_phiV(0);

  scales::set(phasespace::m);
  scales::mcfm();

  //Call the resummed cross section
  double costh = 0;
  double y = phasespace::y;
  int mode = 4;

  double ff[2];
  resint::rint(costh,phasespace::m,phasespace::qt,y,mode,ff);

  f[0] = ff[0]/(8./3.);
  if (opts.helicity >= 0)
    f[1] = ff[1]/(8./3.);

  if (isnan_ofast(f[0]))
    {
      f[0] = 0.;  //avoid nans
      if (opts.helicity >= 0)
	f[1]=0.;
    }

  f[0] = f[0]*jac;
  if (opts.helicity >= 0)
    f[1] = f[1]*jac;

  end_time = clock();
  if (opts.timeprofile)
    cout << setw (3) << "m" << setw(10) << phasespace::m << setw(4) << "y" << setw(10) <<  phasespace::y
	 << setw(8) << "result" << setw(10) << f[0]
	 << setw(10) << "tot time" << setw(10) << float( end_time - begin_time ) /  CLOCKS_PER_SEC
	 << endl;

  if (opts.resumcpp)
    if (!opts.mellin1d)
      rapint::free();

  return 0;
}
/*************************/

int resintegrand3d_cubature_v(unsigned ndim, long unsigned npts, const double x[], void *data, unsigned ncomp, double f[])
{
  tell_to_grid_we_are_alive();

#pragma omp parallel for num_threads(opts.cubacores) copyin(a_param_,rlogs_,resint::a,resint::loga,resint::rloga,evolnative::UVP,evolnative::DVP,evolnative::USP,evolnative::DSP,evolnative::SSP,evolnative::GLP,evolnative::CHP,evolnative::BOP,evolnative::SIP,evolnative::NS3P,evolnative::NS8P,evolnative::NS15P,evolnative::NS24P,evolnative::NS35P,painp_,hfpainp_)

  for (unsigned i = 0; i < npts; i++)
    {
      // evaluate the integrand for npts points
      double xi[ndim];
      double fi[ncomp];
      for (unsigned j = 0; j < ndim; j++)
	xi[j] = x[i*ndim + j];

      resintegrand3d(ndim, xi, ncomp, fi);
      
      for (unsigned k = 0; k < ncomp; ++k)
	f[i*ncomp + k] = fi[k];
    }
  return 0;
}

int resintegrand3d_cubature(unsigned ndim, const double x[], void *data, unsigned ncomp, double f[])
{
  tell_to_grid_we_are_alive();
  resintegrand3d(ndim, x, ncomp, f);
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
  //if (opts.dynamicscale)
  //scaleset_(phasespace::m2);
  scales::set(phasespace::m);
  scales::mcfm();


  //evaluate the resummed cross section
  double ff[2];
  if (opts.resumcpp)
    resint::rint(costh,phasespace::m,phasespace::qt,phasespace::y,mode,ff);
  else
    f[0]=resumm_(costh,phasespace::m,phasespace::qt,phasespace::y,mode)/(8./3.);

  f[0] = ff[0]/(8./3.);
  if (opts.helicity >= 0)
    f[1] = ff[1]/(8./3.);
  
  if (isnan_ofast(f[0]))
    {
      f[0] = 0.;  //avoid nans
      if (opts.helicity >= 0)
	f[1]=0.;
    }
	   
  f[0] = f[0]*jac*swtch;
  if (opts.helicity >= 0)
    f[1] = f[1]*jac*swtch;

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
  double qt = phasespace::qt;
  double y = phasespace::y;

  //Dynamic scale
  //if (opts.dynamicscale)
  //scaleset_(phasespace::m2);
  scales::set(phasespace::m);
  scales::mcfm();
  
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

  double ff[2];
  if (opts.resumcpp)
    resint::rint(costh_CS,m,qt,y,mode,ff);
  else
    f[0]=resumm_(costh_CS,m,qt,y,mode);

  f[0] = ff[0];
  if (opts.helicity >= 0)
    f[1] = ff[1];
  
  if (isnan_ofast(f[0]))
    {
      f[0] = 0.;  //avoid nans
      if (opts.helicity >= 0)
	f[1]=0.;
    }
  
  //apply switching and jacobian
  f[0] = f[0]*jac*psfac*swtch;
  if (opts.helicity >= 0)
    f[1] = f[1]*jac*psfac*swtch;
  
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
