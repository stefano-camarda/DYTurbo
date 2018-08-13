#include "vjloint.h"

#include "settings.h"
#include "codata.h"
#include "mcfm_interface.h"
#include "mesq.h"
#include "pdf.h"
#include "phasespace.h"
#include "omegaintegr.h"
#include "cubature.h"
#include "isnan.h"
#include "gaussrules.h"
#include "KinematicCuts.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iomanip>

const double scutoff = 1e-6;

using namespace std;

double vjloint::muf;
double vjloint::mur;

void vjloint::init()
{
  muf = opts.rmass*opts.kmufac;
  mur = opts.rmass*opts.kmures;
}

void vjloint::calc(const double x[5], double f[2])
{
  f[0] = 0.;
  f[1] = 0.;

  //generate phase space as m, qt, y, costh, phi_lep, x2
  //Jacobian of the change of variables from the unitary hypercube x[6] to the m, qt, y, costh, phi_lep, x2 boundaries
  double jac = 1.;
  bool status;

  //phase space factors (should be (1/2/pi)^3 for each free particle)
  double wtvj=1./8./M_PI; //pq -> V+j phase space
  double wtv=1./8./M_PI/2./2./M_PI; //V-> ll phase space
  double wt=wtvj*wtv/2./M_PI;
  jac = jac * wt;
  
  double r3[3] = {x[0], x[1], x[2]};
  //two alternatives for phase space generation
  /*
  status = phasespace::gen_mqty(r3, jac, true);
  phasespace::calcexpy(); //calculate exp(y) and exp(-y), since they are used in genV4p() and genx2()
  */
  status = phasespace::gen_myqt(r3, jac, true);
  phasespace::calcmt(); //calculate mt, since it is used in genV4p() and genx2()
  
  if (!status)
    return;
  
  jac = jac *2.*phasespace::qt;

  //Set factorization and renormalization scales
  if (opts.dynamicscale)
    {
      muf = phasespace::m*opts.kmufac;
      mur = phasespace::m*opts.kmuren;
      double mur2 = mur*mur;
      scaleset_(mur2);
    }

  //Generate the boson 4-momentum
  phasespace::set_phiV(0.);
  phasespace::genV4p();

  //move dx2 integration here
  //Calculate Bjorken x1 and x2 (integration is performed in dx2)
  status = phasespace::gen_x2(x[3], jac);
  if (!status)
    return;
  /*
  if (jac == 0.)
    {
      //cout << "x2 with jac = 0 " << phasespace::x2 << endl;
      return;
    }
  */

  //Generate incoming partons
  phasespace::genp12();

  //Generate jet from momentum conservation
  phasespace::genp5();

  //reject event if any s(i,j) is too small
  double s15=2.*(phasespace::p1[3]*phasespace::p5[3]-phasespace::p1[0]*phasespace::p5[0]-phasespace::p1[1]*phasespace::p5[1]-phasespace::p1[2]*phasespace::p5[2]);
  double s25=2.*(phasespace::p2[3]*phasespace::p5[3]-phasespace::p2[0]*phasespace::p5[0]-phasespace::p2[1]*phasespace::p5[1]-phasespace::p2[2]*phasespace::p5[2]);
  //if (-s15 < cutoff_.cutoff_ || -s25 < cutoff_.cutoff_)
  if (-s15 < scutoff || -s25 < scutoff)
    {
      //      cout << "failed s cut off " << s15 << "  " << s25 << endl;
      return;
    }


  /*
  //6d integration
  /////////////////////
  //perform dOmega integration (start loop on phi_lep, costh)
  double r2[2] = {x[3], x[4]};
  phasespace::gen_costhphi(r2, jac);
  phasespace::genl4p();
  //calculate V+j matrix elements
  double p[4][mxpart];
  fillp(p);
  
  double msq[11][11];
  if(opts.nproc == 3)
    qqb_z_g_(p,msq);
  else
    qqb_w_g_(p,msq);
  /////////////////////
  */

  /*
  //5d integration
  /////////////////////
  omegaintegr::genV4p(); //--> this function calls genRFaxes for qt-prescription
  phasespace::genRFaxes(phasespace::CS);  //-->overwrite RF axes
  //phasespace::genRFaxes(phasespace::naive);
  
  //generate phi phase space
  phasespace::gen_phi(x[4], jac);

  //start costh integration
  double msq_cth[11][11];
  for (int j = 0; j < 2*MAXNF+1; j++)
    for (int k = 0; k < 2*MAXNF+1; k++)
      msq_cth[j][k] = 0.;
  vector <double> cthmin;
  vector <double> cthmax;
  phasespace::setcthbounds(phasespace::getcthmin(), phasespace::getcthmax());
  omegaintegr::costhbound(phasespace::phi_lep, cthmin, cthmax); //!be carefull, this way costh is integrated according to a given boson rest frame (CS or others)
  vector<double>::iterator itmn;
  vector<double>::iterator itmx;
  itmn = cthmin.begin(); itmx = cthmax.begin();
  for (; itmn != cthmin.end(); itmn++, itmx++)
    {
      //cout << *itmn << "  " << *itmx << "  " << cthmin.size() << endl;
      double cthc=0.5*(*itmn+*itmx);
      double cthm=0.5*(*itmx-*itmn);
      
      int cthrule = 4;
      for(int i=0; i < cthrule; i++)
	{
	  double xcth = cthc+cthm*gr::xxx[cthrule-1][i];
	  
	  //Generate leptons 4-momenta: p3 is the lepton and p4 is the antilepton
	  phasespace::set_cth(xcth);
	  //phasespace::genRFaxes(CS);
	  phasespace::genl4p();

	  //calculate V+j matrix elements
	  double p[4][mxpart];
	  fillp(p);

	  double msq[11][11];
	  if(opts.nproc == 3)
	    qqb_z_g_(p,msq);
	  else
	    qqb_w_g_(p,msq);

	  double cth = xcth*(phasespace::y >= 0. ? 1 : -1);
	  double sth = sqrt(max(0.,1.-cth*cth));
	  double s2th = 2*xcth*sth;
	  double phi = phasespace::phi_lep*(phasespace::y > 0. ? 1 : -1);
	  double cphi = cos(phi);
	  double c2phi = 2*cphi*cphi-1.;
	  double sphi = sqrt(max(0.,1.-cphi*cphi))*(phi>0?1:-1);
	  double s2phi = 2*cphi*sphi;

	  //Angular polynomials
	  //NEWKIN( A5 ){ SinThCS  sinth;  Sin2PhiCS sin2ph; double calc(){ return 5.     * (sinth()*sinth()*sin2ph() ) ;       } };
	  //NEWKIN( A6 ){ Sin2ThCS sin2th; SinPhiCS  sinph;  double calc(){ return 4.     * (sin2th()*sinph()         ) ;       } };
	  //NEWKIN( A7 ){ SinThCS  sinth;  SinPhiCS  sinph;  double calc(){ return 4.     * (sinth()*sinph()          ) ;       } };

	  double pol = 1.;
	  if (opts.helicity == 0)
	    pol = 20./3.*(0.5-1.5*cth*cth) +2./3.;
	  else if (opts.helicity == 1)
	    pol = 5.*(s2th*cphi);
	  else if (opts.helicity == 2)
	    pol = 10.*(sth*sth*c2phi);
	  else if (opts.helicity == 3)
	    pol = 4.*(sth*cphi);
	  else if (opts.helicity == 4)
	    pol = 4.*cth;
	  //cout << cth << "  " << phi << endl;
	  for (int j = 0; j < 2*MAXNF+1; j++)
	    for (int k = 0; k < 2*MAXNF+1; k++)
	      if (msq[k][j] != 0.)
		msq_cth[k][j] += msq[k][j] *gr::www[cthrule-1][i]*cthm * pol;
	  
	} //end dcosth loop
    } //end loop on costh boundaries
  double msq[11][11];
  for (int j = 0; j < 2*MAXNF+1; j++)
    for (int k = 0; k < 2*MAXNF+1; k++)
      msq[j][k] = msq_cth[j][k];
  /////////////////////
  */

  //4d integration
  //dOmega integration
  omegaintegr::genV4p(); //--> this function calls genRFaxes for qt-prescription
  phasespace::genRFaxes(phasespace::CS);  //--> overwrite RF axes: RF here determines the frame for the coefficients if opts.helicity != -1, and/or the RF for the costh boundaries
  //phasespace::genRFaxes(phasespace::naive);

  double msq_omega[11][11];
  double msq_omega_m[11][11];
  for (int j = 0; j < 2*MAXNF+1; j++)
    for (int k = 0; k < 2*MAXNF+1; k++)
      {
	msq_omega[j][k] = 0.;
	msq_omega_m[j][k] = 0.;
      }
  
  //start phi integration
  int phiintervals, phirule;
  if (opts.makecuts)
    {
      phiintervals = opts.phiintervals;
      phirule = opts.phirule;
    }
  else
    {
      phiintervals = 1;
      phirule = opts.vjphirule;
    }
  double phi1 = -M_PI;
  double phi2 = M_PI;
  double hphi=(phi2-phi1)/phiintervals;
  for(int i=0; i < phiintervals; i++)
    {
      double phimin = i*hphi+phi1;
      double phimax = (i+1)*hphi+phi1;
      double phic=0.5*(phimin+phimax);
      double phim=0.5*(phimax-phimin);
      for(int iphi = 0; iphi < phirule; iphi++)
	{
	  double xphi = phic+phim*gr::xxx[phirule-1][iphi];
	  phasespace::set_philep(xphi);

	  double phi, cphi, c2phi;
	  //double sphi, s2phi;
	  if (opts.helicity >= 0)
	    {
	      phi = phasespace::phi_lep*(phasespace::y > 0. ? 1 : -1);
	      cphi = cos(phi);
	      c2phi = 2*cphi*cphi-1.;
	      //sphi = sqrt(max(0.,1.-cphi*cphi))*(phi>0?1:-1);
	      //s2phi = 2*cphi*sphi;
	    }
      
	  //start costh integration
	  vector <double> cthmin;
	  vector <double> cthmax;
	  phasespace::setcthbounds(phasespace::getcthmin(), phasespace::getcthmax()); //Allows for switching integration boundaries between positive and negative rapidity (currently not thread safe!!!)
	  omegaintegr::costhbound(phasespace::phi_lep, cthmin, cthmax); //!be carefull, this way costh is integrated according to a given boson rest frame (CS or others)
	  vector<double>::iterator itmn;
	  vector<double>::iterator itmx;
	  itmn = cthmin.begin(); itmx = cthmax.begin();
	  for (; itmn != cthmin.end(); itmn++, itmx++)
	    {
	      //cout << *itmn << "  " << *itmx << "  " << cthmin.size() << endl;
	      double cthc=0.5*(*itmn+*itmx);
	      double cthm=0.5*(*itmx-*itmn);
      
	      int cthrule = 3;
	      for(int i=0; i < cthrule; i++)
		{
		  double xcth = cthc+cthm*gr::xxx[cthrule-1][i];
	  
		  //Generate leptons 4-momenta: p3 is the lepton and p4 is the antilepton
		  phasespace::set_cth(xcth);
		  //phasespace::genRFaxes(CS);
		  phasespace::genl4p();

		  //calculate V+j matrix elements
		  double p[4][mxpart];
		  fillp(p);

		  double msq[11][11];
		  if(opts.nproc == 3)
		    //qqb_z_g_(p,msq);
		    qqb_z1jet_(p,msq);
		  else
		    qqb_w_g_(p,msq);
		  double cth, sth, s2th;
		  if (opts.helicity >= 0)
		    {
		      cth = xcth*(phasespace::y >= 0. ? 1 : -1);
		      sth = sqrt(max(0.,1.-cth*cth));
		      s2th = 2*cth*sth;
		    }

		  double poly = 1.;
		  if (opts.helicity == 0)
		    poly = 20./3.*(0.5-1.5*cth*cth) +2./3.;
		  else if (opts.helicity == 1)
		    poly = 5.*(s2th*cphi);
		  else if (opts.helicity == 2)
		    poly = 10.*(sth*sth*c2phi);
		  else if (opts.helicity == 3)
		    poly = 4.*(sth*cphi);
		  else if (opts.helicity == 4)
		    poly = 4.*cth;
		  else if (opts.helicity > 4)
		    poly = 0.;
		  //Angular polynomials
		  //NEWKIN( A5 ){ SinThCS  sinth;  Sin2PhiCS sin2ph; double calc(){ return 5.     * (sinth()*sinth()*sin2ph() ) ;       } };
		  //NEWKIN( A6 ){ Sin2ThCS sin2th; SinPhiCS  sinph;  double calc(){ return 4.     * (sin2th()*sinph()         ) ;       } };
		  //NEWKIN( A7 ){ SinThCS  sinth;  SinPhiCS  sinph;  double calc(){ return 4.     * (sinth()*sinph()          ) ;       } };

		  for (int j = 0; j < 2*MAXNF+1; j++)
		    for (int k = 0; k < 2*MAXNF+1; k++)
		      if (msq[k][j] != 0.)
			{
			  msq_omega[k][j] += msq[k][j] *gr::www[cthrule-1][i]*cthm *gr::www[phirule-1][iphi]*phim;
			  //msq_omega[k][j] += gr::www[cthrule-1][i]*cthm *gr::www[phirule-1][iphi]*phim;
			  if (opts.helicity >= 0)
			    msq_omega_m[k][j] += msq[k][j] *gr::www[cthrule-1][i]*cthm *gr::www[phirule-1][iphi]*phim * poly;
			}
	  
		} //end dcosth loop
	    } //end loop on costh boundaries
	}
    }//end philep integration
  //double msq[11][11];
  //for (int j = 0; j < 2*MAXNF+1; j++)
  //for (int k = 0; k < 2*MAXNF+1; k++)
  //msq[j][k] = msq_omega[j][k];

  /////////////////////

  // Load central PDF and QCD coupling
  //if (pdferr) then
  //int npdf = 0;
  //dysetpdf_(npdf);
  //double gsqcentral=gsq;
  //     skip PDF loop in the preconditioning phase
  //int maxpdf=0;
  //if (doFill.ne.0) maxpdf = totpdf-1
      
  //     start PDF loop
  //  for (int np=0; np <= maxpdf; np++)
  //    {
  //      dysetpdf_(np);
  //      hists_setpdf_(np);
  //     intitialise xmsq to 0
  //cout << muf << "  " << x1 << "  " << x2 << endl;
  // calculate PDFs
  double fx1[2*MAXNF+1],fx2[2*MAXNF+1];
  fdist_(opts.ih1,phasespace::x1,muf,fx1);
  fdist_(opts.ih2,phasespace::x2,muf,fx2);

  double xmsq = convolute(fx1, fx2, msq_omega);
  double xmsq_m;
  if (opts.helicity >= 0)
    xmsq_m = convolute(fx1, fx2, msq_omega_m);

  double shad = pow(opts.sroot,2);
  xmsq *= gevfb/(2.*phasespace::x1*phasespace::x2*shad);
  xmsq_m *= gevfb/(2.*phasespace::x1*phasespace::x2*shad);
  //double fbGeV2 = 0.389379e12;
  //xmsq=xmsq*fbGeV2/(2.*phasespace::x1*phasespace::x2*shad);
  
  //f[np+1]=xmsq;;

  //if (doFill.ne.0) then
  // {
  //   val=xmsq*wgt
  //   hists_fill_(p3,p4,val);
  // }

  //} end PDFs loop

  //cout << phasespace::x2 << "  " << xmsq << "  " << jac << endl;
  
  if (isnan_ofast(xmsq))
    {
      //      cout << "xmsq in vjloint is nan" << endl;
      return;
    }

  f[0] = xmsq*jac;
  f[1] = xmsq_m*jac;
  return;
}

double vjloint::calcvegas(const double x[7])
{
  clock_t begin_time, end_time;

  begin_time = clock();

  //generate phase space as m, qt, y, costh, phi_lep, x2
  //Jacobian of the change of variables from the unitary hypercube x[6] to the m, qt, y, costh, phi_lep, x2 boundaries
  double jac = 1.;
  bool status;
  
  //phase space factors (should be (1/2/pi)^3 for each free particle)
  double wtvj=1./8./M_PI; //pq -> V+j phase space
  double wtv=1./8./M_PI/2./2./M_PI; //V-> ll phase space
  double wt=wtvj*wtv/2./M_PI;
  jac = jac * wt;
  
  double r3[3] = {x[0], x[1], x[2]};
  //two alternatives for phase space generation
  //status = phasespace::gen_mqty(r3, jac, true);
  //phasespace::calcexpy(); //calculate exp(y) and exp(-y), since they are used in genV4p() and genx2()
  status = phasespace::gen_myqt(r3, jac, true);
  phasespace::calcmt(); //calculate mt, since it is used in genV4p() and genx2()
  
  if (!status)
    return 0.;
  
  jac = jac *2.*phasespace::qt;

  //Set factorization and renormalization scales
  if (opts.dynamicscale)
    {
      muf = phasespace::m*opts.kmufac;
      mur = phasespace::m*opts.kmuren;
      double mur2 = mur*mur;
      scaleset_(mur2); //set renormalization and factorization scales, and calculate ason2pi and ason4pi
    }

  //Generate the boson 4-momentum
  phasespace::set_phiV(-M_PI+2.*M_PI*x[6]);
  phasespace::genRFaxes(phasespace::naive);
  phasespace::genV4p();

  //move dx2 integration here
  //Calculate Bjorken x1 and x2 (integration is performed in dx2)
  phasespace::gen_x2(x[5], jac);
  if (jac == 0.)
    {
      //cout << "x2 with jac = 0 " << phasespace::x2 << endl;
      return 0.;
    }

  //Generate incoming partons
  phasespace::genp12();

  //Generate jet from momentum conservation
  phasespace::genp5();

  //reject event if any s(i,j) is too small
  double s15=2.*(phasespace::p1[3]*phasespace::p5[3]-phasespace::p1[0]*phasespace::p5[0]-phasespace::p1[1]*phasespace::p5[1]-phasespace::p1[2]*phasespace::p5[2]);
  double s25=2.*(phasespace::p2[3]*phasespace::p5[3]-phasespace::p2[0]*phasespace::p5[0]-phasespace::p2[1]*phasespace::p5[1]-phasespace::p2[2]*phasespace::p5[2]);
  //if (-s15 < cutoff_.cutoff_ || -s25 < cutoff_.cutoff_)
  if (-s15 < scutoff || -s25 < scutoff)
    {
      //      cout << "failed s cut off " << s15 << "  " << s25 << endl;
      return 0.;
    }

  //perform dOmega integration (start loop on phi_lep, costh)
  double r2[2] = {x[3], x[4]};
  phasespace::gen_costhphi(r2, jac);
  phasespace::genl4p();

  if (!Kinematics::Cuts::KeepThisEvent(phasespace::p3, phasespace::p4))
    return 0.;

  //calculate V+j matrix elements
  double p[4][mxpart];
  fillp(p);
  
  double msq[11][11];
  if(opts.nproc == 3)
    //qqb_z_g_(p,msq);
    qqb_z1jet_(p,msq);
  else
    qqb_w_g_(p,msq);

  // Load central PDF and QCD coupling
  //if (pdferr) then
  //int npdf = 0;
  //dysetpdf_(npdf);
  //double gsqcentral=gsq;
  //     skip PDF loop in the preconditioning phase
  //int maxpdf=0;
  //if (doFill.ne.0) maxpdf = totpdf-1
      
  //     start PDF loop
  //  for (int np=0; np <= maxpdf; np++)
  //    {
  //      dysetpdf_(np);
  //      hists_setpdf_(np);
  //     intitialise xmsq to 0
  //cout << muf << "  " << x1 << "  " << x2 << endl;
  // calculate PDFs
  double fx1[2*MAXNF+1],fx2[2*MAXNF+1];
  fdist_(opts.ih1,phasespace::x1,muf,fx1);
  fdist_(opts.ih2,phasespace::x2,muf,fx2);

  double xmsq = convolute(fx1, fx2, msq);

  double shad = pow(opts.sroot,2);
  xmsq=xmsq*gevfb/(2.*phasespace::x1*phasespace::x2*shad);
  
  //f[np+1]=xmsq;;

  //if (doFill.ne.0)
  //{
  double val=xmsq*jac;

  //} end PDFs loop

  if (isnan_ofast(xmsq))
    {
      //      cout << "xmsq in vjloint is nan" << endl;
      return 0.;
    }

  end_time = clock();
  if (opts.timeprofile)
    cout << setw (3) << "m" << setw(10) << phasespace::m << setw(4) << "qt" << setw(10) <<  phasespace::qt
	 << setw(8) << "result" << setw(10) << xmsq*jac
      	 << setw(10) << "tot time" << setw(10) << float( end_time - begin_time ) /  CLOCKS_PER_SEC
	 << endl;
  
  return xmsq*jac;
}

void vjloint::fillp(double p[4][mxpart])
{
  p[0][0] = phasespace::p1[0];
  p[1][0] = phasespace::p1[1];
  p[2][0] = phasespace::p1[2];
  p[3][0] = phasespace::p1[3];

  p[0][1] = phasespace::p2[0];
  p[1][1] = phasespace::p2[1];
  p[2][1] = phasespace::p2[2];
  p[3][1] = phasespace::p2[3];

  p[0][2] = phasespace::p3[0];
  p[1][2] = phasespace::p3[1];
  p[2][2] = phasespace::p3[2];
  p[3][2] = phasespace::p3[3];

  p[0][3] = phasespace::p4[0];
  p[1][3] = phasespace::p4[1];
  p[2][3] = phasespace::p4[2];
  p[3][3] = phasespace::p4[3];

  p[0][4] = phasespace::p5[0];
  p[1][4] = phasespace::p5[1];
  p[2][4] = phasespace::p5[2];
  p[3][4] = phasespace::p5[3];
}

double vjloint::convolute(double fx1[2*MAXNF+1], double fx2[2*MAXNF+1], double msq[2*MAXNF+1][2*MAXNF+1])
{
  double xmsq;
  if (opts.nproc == 1)
    {
      xmsq =fx1[0]*fx2[5]*msq[5][0];
      xmsq+=fx1[0]*fx2[7]*msq[7][0];
      xmsq+=fx1[0]*fx2[9]*msq[9][0];
      xmsq+=fx1[2]*fx2[5]*msq[5][2];
      xmsq+=fx1[2]*fx2[7]*msq[7][2];
      xmsq+=fx1[2]*fx2[9]*msq[9][2];
      xmsq+=fx1[4]*fx2[5]*msq[5][4];
      xmsq+=fx1[4]*fx2[7]*msq[7][4];
      xmsq+=fx1[4]*fx2[9]*msq[9][4];
      xmsq+=fx1[5]*fx2[0]*msq[0][5];
      xmsq+=fx1[5]*fx2[2]*msq[2][5];
      xmsq+=fx1[5]*fx2[4]*msq[4][5];
      xmsq+=fx1[5]*fx2[7]*msq[7][5];
      xmsq+=fx1[5]*fx2[9]*msq[9][5];
      xmsq+=fx1[7]*fx2[2]*msq[2][7];
      xmsq+=fx1[7]*fx2[4]*msq[4][7];
      xmsq+=fx1[7]*fx2[5]*msq[5][7];
      xmsq+=fx1[9]*fx2[0]*msq[0][9];
      xmsq+=fx1[9]*fx2[2]*msq[2][9];
      xmsq+=fx1[9]*fx2[4]*msq[4][9];
      xmsq+=fx1[9]*fx2[5]*msq[5][9];
    }
  else if (opts.nproc == 2)
    {
      xmsq =fx1[1 ]*fx2[5 ]*msq[5 ][1 ];
      xmsq+=fx1[1 ]*fx2[6 ]*msq[6 ][1 ];
      xmsq+=fx1[1 ]*fx2[8 ]*msq[8 ][1 ];
      xmsq+=fx1[1 ]*fx2[10]*msq[10][1 ];
      xmsq+=fx1[3 ]*fx2[5 ]*msq[5 ][3 ];
      xmsq+=fx1[3 ]*fx2[6 ]*msq[6 ][3 ];
      xmsq+=fx1[3 ]*fx2[8 ]*msq[8 ][3 ];
      xmsq+=fx1[3 ]*fx2[10]*msq[10][3 ];
      xmsq+=fx1[5 ]*fx2[1 ]*msq[1 ][5 ];
      xmsq+=fx1[5 ]*fx2[3 ]*msq[3 ][5 ];
      xmsq+=fx1[5 ]*fx2[6 ]*msq[6 ][5 ];
      xmsq+=fx1[5 ]*fx2[8 ]*msq[8 ][5 ];
      xmsq+=fx1[5 ]*fx2[10]*msq[10][5 ];
      xmsq+=fx1[6 ]*fx2[1 ]*msq[1 ][6 ];
      xmsq+=fx1[6 ]*fx2[3 ]*msq[3 ][6 ];
      xmsq+=fx1[6 ]*fx2[5 ]*msq[5 ][6 ];
      xmsq+=fx1[8 ]*fx2[1 ]*msq[1 ][8 ];
      xmsq+=fx1[8 ]*fx2[3 ]*msq[3 ][8 ];
      xmsq+=fx1[8 ]*fx2[5 ]*msq[5 ][8 ];
      xmsq+=fx1[10]*fx2[1 ]*msq[1 ][10];
      xmsq+=fx1[10]*fx2[3 ]*msq[3 ][10];
      xmsq+=fx1[10]*fx2[5 ]*msq[5 ][10];
    }
  else if (opts.nproc == 3)
    {
      xmsq =fx1[0 ]*fx2[5 ]*msq[5 ][0 ];
      xmsq+=fx1[0 ]*fx2[10]*msq[10][0 ];
      xmsq+=fx1[1 ]*fx2[5 ]*msq[5 ][1 ];
      xmsq+=fx1[1 ]*fx2[9 ]*msq[9 ][1 ];
      xmsq+=fx1[2 ]*fx2[5 ]*msq[5 ][2 ];
      xmsq+=fx1[2 ]*fx2[8 ]*msq[8 ][2 ];
      xmsq+=fx1[3 ]*fx2[5 ]*msq[5 ][3 ];
      xmsq+=fx1[3 ]*fx2[7 ]*msq[7 ][3 ];
      xmsq+=fx1[4 ]*fx2[5 ]*msq[5 ][4 ];
      xmsq+=fx1[4 ]*fx2[6 ]*msq[6 ][4 ];
      xmsq+=fx1[5 ]*fx2[0 ]*msq[0 ][5 ];
      xmsq+=fx1[5 ]*fx2[1 ]*msq[1 ][5 ];
      xmsq+=fx1[5 ]*fx2[2 ]*msq[2 ][5 ];
      xmsq+=fx1[5 ]*fx2[3 ]*msq[3 ][5 ];
      xmsq+=fx1[5 ]*fx2[4 ]*msq[4 ][5 ];
      xmsq+=fx1[5 ]*fx2[6 ]*msq[6 ][5 ];
      xmsq+=fx1[5 ]*fx2[7 ]*msq[7 ][5 ];
      xmsq+=fx1[5 ]*fx2[8 ]*msq[8 ][5 ];
      xmsq+=fx1[5 ]*fx2[9 ]*msq[9 ][5 ];
      xmsq+=fx1[5 ]*fx2[10]*msq[10][5 ];
      xmsq+=fx1[6 ]*fx2[4 ]*msq[4 ][6 ];
      xmsq+=fx1[6 ]*fx2[5 ]*msq[5 ][6 ];
      xmsq+=fx1[7 ]*fx2[3 ]*msq[3 ][7 ];
      xmsq+=fx1[7 ]*fx2[5 ]*msq[5 ][7 ];
      xmsq+=fx1[8 ]*fx2[2 ]*msq[2 ][8 ];
      xmsq+=fx1[8 ]*fx2[5 ]*msq[5 ][8 ];
      xmsq+=fx1[9 ]*fx2[1 ]*msq[1 ][9 ];
      xmsq+=fx1[9 ]*fx2[5 ]*msq[5 ][9 ];
      xmsq+=fx1[10]*fx2[0 ]*msq[0 ][10];
      xmsq+=fx1[10]*fx2[5 ]*msq[5 ][10];
    }
  /*
  double xmsq = 0.;
  for (int j = 0; j < 2*MAXNF+1; j++)
    for (int k = 0; k < 2*MAXNF+1; k++)
      {
	//if (isnan_ofast(msq[k][j]))
	//cout << j << "  " << k << "  " << fx1[j] << "  " << fx2[k] << "  " << msq[k][j] << endl;
	if (msq[k][j] == 0.) continue;
	//     gsq/gsqcentral correct for a possibly different value of alphas in the PDF (at O(alphas))
	xmsq=xmsq+fx1[j]*fx2[k]*msq[k][j]; //*(gsq/gsqcentral)
	//cout << j << "  " << k << "  " << fx1[j] << "  " << fx2[k] << "  " << msq[k][j] << endl;
      }
  */

  return xmsq;
}
  
