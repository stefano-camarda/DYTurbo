#include "bornintegr.h"
#include "interface.h"
#include "phasespace.h"
#include "settings.h"
#include "loint.h"
#include "lomellin.h"
#include "cubacall.h"
#include "isnan.h"
#include "evolnative.h"
#include "pegasus.h"
#include "pmom.h"
#include "ccoeff.h"
#include "resint.h"
#include "dyres_interface.h"

#include <iostream>
#include <iomanip>
//#include <ctime>
//#include <math.h>


#include "KinematicCuts.h"

//using namespace std;

integrand_t doublevirtintegrand(const int &ndim, const double x[], const int &ncomp, double f[],
                        void* userdata, const int &nvec, const int &core,
                        double &weight, const int &iter)
{
  
  double rvv[22];
  for (int i = 0; i < ndim; i++)
    rvv[i]=x[i];

  dofill_.doFill_ = int(iter==last_iter);
  f[0] = lowinthst_dynnlo_(rvv,weight,f);

  tell_to_grid_we_are_alive();
  return 0;
}

integrand_t lointegrandMC(const int &ndim, const double x[], const int &ncomp, double f[],
			 void* userdata, const int &nvec, const int &core,
			 double &weight, const int &iter)
{
  clock_t begin_time, end_time;
  begin_time = clock();

  f[0]=0.;
  
  if (opts.fixedorder && phasespace::qtmin > 0)
    return 0;
  
  //Jacobian of the change of variables from the unitary hypercube x[3] to the m, y boundaries
  double jac = 1.;
  bool status = true;

  //generate phase space
  double r1[2] = {x[0], x[1]};
  status = phasespace::gen_my(r1, jac);
  if (!status)
    return 0;

  phasespace::set_qt(0.);
  phasespace::set_phiV(0.);
  
  double r2[2] = {x[2], x[3]};
  phasespace::gen_costhphi(r2, jac);

  //generate boson and leptons
  phasespace::calcexpy();
  phasespace::calcmt();
  phasespace::genV4p();
  phasespace::genl4p();

  //apply cuts
  if (!Kinematics::Cuts::KeepThisEvent(phasespace::p3,phasespace::p4))
    return 0;
  
  //calculate integrand
  double f2[2];
  loint::lint(phasespace::costh,phasespace::m,phasespace::y,0,f2);
  f[0] = f2[0]*jac;
  
  //fill histograms
  if (iter==last_iter)
    {
      double wt = weight*f[0];
      hists_fill_(phasespace::p3,phasespace::p4,&wt);
    }

  end_time = clock();
  if (opts.timeprofile)
    cout << setw (3) << "m" << setw(10) << phasespace::m
	 << setw(4) << "qt" << setw(10) <<  phasespace::qt
      	 << setw(4) << "y" << setw(10) <<  phasespace::y
	 << setw(8) << "result" << setw(10) << f[0]
	 << setw(10) << "tot time" << setw(10) << float( end_time - begin_time ) /  CLOCKS_PER_SEC
	 << endl;
  
  tell_to_grid_we_are_alive();
  return 0;
}

int lointegrand2d_cubature_v(unsigned ndim, long unsigned npts, const double x[], void *data, unsigned ncomp, double f[])
{
#pragma omp parallel for num_threads(opts.cubacores) copyin(scale_,facscale_,qcdcouple_)
  for (unsigned i = 0; i < npts; i++)
    {
      // evaluate the integrand for npts points
      double xi[ndim];
      double fi[ncomp];
      for (unsigned j = 0; j < ndim; j++)
	xi[j] = x[i*ndim + j];

      lointegrand2d(ndim, xi, ncomp, fi);
      
      for (unsigned k = 0; k < ncomp; ++k)
	f[i*ncomp + k] = fi[k];
    }
  return 0;
}

int lointegrand2d_cubature(unsigned ndim, const double x[], void *data, unsigned ncomp, double f[])
{
  lointegrand2d(ndim, x, ncomp, f);
  return 0;
}

double lointegrand2d_smolyak(int ndim, double x[])
{
  int ncomp = 1;
  double f[ncomp];
  lointegrand2d(ndim, x, ncomp, f);
  return f[0];
}

integrand_t lointegrand2d(const int &ndim, const double x[], const int &ncomp, double f[])
{
  clock_t begin_time, end_time;
  begin_time = clock();

  f[0]=0.;
  if (opts.helicity >= 0)
    f[1]=0.;
  
  if (opts.fixedorder && phasespace::qtmin > 0)
    return 0;

  //Jacobian of the change of variables from the unitary hypercube x[3] to the m, y boundaries
  double jac = 1.;
  bool status = true;

  double r1[2] = {x[0], x[1]};
  status = phasespace::gen_my(r1, jac);
  if (!status)
    return 0;

  phasespace::set_qt(0.);

  double m = phasespace::m;
  double y = phasespace::y;

  double f2[2];
  loint::lint(0.,phasespace::m,phasespace::y,1,f2);

  f[0] = f2[0]*jac;
  if (opts.helicity >= 0)
    f[1] = f2[1]*jac;
  
  //perform phi integration
  f[0] = f[0] * 2.*M_PI;
  if (opts.helicity >= 0)
    f[1]*=2.*M_PI;

  end_time = clock();
  if (opts.timeprofile)
    cout << setw (3) << "m" << setw(10) << phasespace::m
	 << setw(4) << "qt" << setw(10) <<  phasespace::qt
      	 << setw(4) << "y" << setw(10) <<  phasespace::y
	 << setw(8) << "result" << setw(10) << f[0]
	 << setw(10) << "tot time" << setw(10) << float( end_time - begin_time ) /  CLOCKS_PER_SEC
	 << endl;
  
  tell_to_grid_we_are_alive();
  return 0;
}


int lointegrand1d_cubature_v(unsigned ndim, long unsigned npts, const double x[], void *data, unsigned ncomp, double f[])
{
#pragma omp parallel for num_threads(opts.cubacores) \
  copyin(scale_,facscale_,qcdcouple_) \
  copyin(a_param_,rlogs_,resint::a,resint::loga,resint::rloga,\
	 mellinint::wn,mellinint::Np,mellinint::Nm,mellinint::wn_1,mellinint::Np_1,mellinint::Nm_1,mellinint::wn_2,mellinint::Np_2,mellinint::Nm_2,\
	 evolnative::ans,evolnative::am,evolnative::ap,evolnative::al,evolnative::be,evolnative::ab,\
	 evolnative::rmin,evolnative::rplus,evolnative::rqq,evolnative::rqg,evolnative::rgq,evolnative::rgg,\
	 evolnative::rmmqq,evolnative::rmmqg,evolnative::rmmgq,evolnative::rmmgg,evolnative::rmpqq,evolnative::rmpqg,evolnative::rmpgq,evolnative::rmpgg,\
	 evolnative::rpmqq,evolnative::rpmqg,evolnative::rpmgq,evolnative::rpmgg,evolnative::rppqq,evolnative::rppqg,evolnative::rppgq,evolnative::rppgg,\
	 pegasus::gli,pegasus::vai,pegasus::m3i,pegasus::m8i,pegasus::m15i,pegasus::m24i,pegasus::sgi,pegasus::p3i,pegasus::p8i,pegasus::p15i,pegasus::p24i,\
	 order_,painp_,hfpainp_,pacthr_,pabthr_,patthr_,moms_,hsums_,pns0_,psg0_,pns1_,psg1_,pns2_,psg2_,ans2_,asg2_,lsg_,spsums_,u1sg_,r1sg_,u1hsg_,u2sg_,r2sg_,u2hsg_,u2ns_, \
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

      lointegrand1d(ndim, xi, ncomp, fi);
      
      for (unsigned k = 0; k < ncomp; ++k)
	f[i*ncomp + k] = fi[k];
    }
  return 0;
}

int lointegrand1d_cubature(unsigned ndim, const double x[], void *data, unsigned ncomp, double f[])
{
  lointegrand1d(ndim, x, ncomp, f);
  return 0;
}
integrand_t lointegrand1d(const int &ndim, const double x[], const int &ncomp, double f[])
{
  tell_to_grid_we_are_alive();
  clock_t begin_time, end_time;
  begin_time = clock();

  f[0]=0.;
  if (opts.helicity >= 0)
    f[1]=0.;

  if (opts.fixedorder && phasespace::qtmin > 0)
    return 0;
  
  //Jacobian of the change of variables from the segment [0,1] to m boundaries
  double jac = 1.;
  bool status = true;

  //Generate the boson invariant mass between the integration boundaries
  //kinematic limit from the relation x(1,2) < 1 ==> exp(fabs(y)) < sqrt(s)/m
  double absymn;
  if (phasespace::ymin*phasespace::ymax <= 0.)
    absymn = 0.;
  else
    absymn = min(fabs(phasespace::ymin),fabs(phasespace::ymax));
  double mlim = opts.sroot/exp(absymn);
  double r1 = {x[0]};
  status = phasespace::gen_m(r1, jac, mlim, false, false);
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

  phasespace::set_qt(0.);
  phasespace::set_phiV(0);

  //scales::set(phasespace::m);
  //scales::mcfm();

  //Call the resummed cross section
  double costh = 0;
  double y = 0.5*(ymx+ymn);
  phasespace::set_y(y); //needed to get the correct mlogz factor in mellinint::update() ?
  int mode = 2;

  double ff[2];
  lomellin::lint(costh,phasespace::m,phasespace::y,mode,ff);

  f[0] = ff[0];///(8./3.);
  if (opts.helicity >= 0)
    f[1] = ff[1];///(8./3.);
  
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
    cout << setw (3) << "m" << setw(10) << phasespace::m << setw(4)
	 << setw(8) << "result" << setw(10) << f[0]
	 << setw(10) << "tot time" << setw(10) << float( end_time - begin_time ) /  CLOCKS_PER_SEC
	 << endl;

  return 0;
}
