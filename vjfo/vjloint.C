#include "vjloint.h"

#include "settings.h"
#include "mcfm_interface.h"
#include "cuts.h"
#include "mesq.h"
#include "pdf.h"
#include "phasespace.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iomanip>

using namespace std;

double vjloint::jac;

void vjloint::init()
{
}

double vjloint::vlint(const double x[6])
{
  genps(x);
  //generate phase space as m, qt, y, costh, phi_lep, x2
  //Jacobian of the change of variables from the unitary hypercube x[6] to the m, qt, y, costh, phi_lep, x2 boundaries
  jac = 1.;

  double r3[3] = {x[0], x[1], x[2]};
  phasespace::gen_mqty(r3, jac);
  if (jac == 0.)
    return 0.;

  //Set factorization scale
  double muf;
  if (opts.dynamicscale)
    muf = phasespace::m*opts.kmufac;
  else
    muf = opts.rmass*opts.kmufac;
  
  //Generate the boson 4-momentum
  phasespace::set_phiV(0.);
  phasespace::genV4p();

  //perform dOmega integration (start loop on phi_lep, costh)
  double r2[2] = {x[3], x[4]};
  phasespace::gen_costhphi(r2, jac);

  //Generate leptons 4-momenta
  //p4 is the lepton and p3 is the antilepton
  phasespace::genl4p();

  //Apply cuts on leptons
  if (!cuts::lep(phasespace::p3,phasespace::p4))
    return 0.;
  //jac = jac*binner_(p3,p4);

  //perform dx2 integration (start loop on x2)
  //Calculate Bjorken x1 and x2 (integration is performed in dx2)
  phasespace::gen_x2(x[5], jac);
  if (jac == 0.)
    return 0.;

  //Generate incoming partons
  phasespace::genp12();

  //Generate jet from momentum conservation
  phasespace::genp5();

  jac = jac *2.*phasespace::qt;

  //phase space factors (should be (1/2/pi)^3 for each free particle)
  double wtvj=1./8./M_PI; //pq -> V+j phase space
  double wtv=1./8./M_PI/2./2./M_PI; //V-> ll phase space
  double wt=wtvj*wtv/2./M_PI;

  jac = jac * wt;

  //reject event if any s(i,j) is too small
  double s15=2.*(phasespace::p1[3]*phasespace::p5[3]-phasespace::p1[0]*phasespace::p5[0]-phasespace::p1[1]*phasespace::p5[1]-phasespace::p1[2]*phasespace::p5[2]);
  double s25=2.*(phasespace::p2[3]*phasespace::p5[3]-phasespace::p2[0]*phasespace::p5[0]-phasespace::p2[1]*phasespace::p5[1]-phasespace::p2[2]*phasespace::p5[2]);
  if (-s15 < cutoff_.cutoff_ || -s25 < cutoff_.cutoff_)
    return 0.;
  
  //calculate V+j matrix elements
  double p[4][12];
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

  double msq[11][11];
  if(opts.nproc == 3)
    qqb_z_g_(p,msq);
  else
    qqb_w_g_(p,msq);

  // Load central PDF and QCD coupling
  //if (pdferr) then
  int npdf = 0;
  dysetpdf_(npdf);

  //double gsqcentral=gsq;

  //     skip PDF loop in the preconditioning phase
  int maxpdf=0;
  //if (doFill.ne.0) maxpdf = totpdf-1
      
  //     start PDF loop
  //  for (int np=0; np <= maxpdf; np++)
  //    {
  //      dysetpdf_(np);
  //      hists_setpdf_(np);
  //     intitialise xmsq to 0
  double xmsq=0.;

  //cout << muf << "  " << x1 << "  " << x2 << endl;
  // calculate PDFs
  double fx1[2*MAXNF+1],fx2[2*MAXNF+1];
  fdist_(opts.ih1,phasespace::x1,muf,fx1);
  fdist_(opts.ih2,phasespace::x2,muf,fx2);

  if (opts.nproc == 3)
    {
      xmsq+=fx1[0 ]*fx2[5 ]*msq[5 ][0 ]; 
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
  else
    for (int j = 0; j < 2*MAXNF+1; j++)
      for (int k = 0; k < 2*MAXNF+1; k++)
	{
	  //cout << j << "  " << k << "  " << fx1[j] << "  " << fx2[k] << "  " << msq[k][j] << endl;
	  if (msq[k][j] == 0.) continue;
	  //     gsq/gsqcentral correct for a possibly different value of alphas in the PDF (at O(alphas))
	  xmsq=xmsq+fx1[j]*fx2[k]*msq[k][j]; //*(gsq/gsqcentral)
	  //cout << j << "  " << k << "  " << fx1[j] << "  " << fx2[k] << "  " << msq[k][j] << endl;
	}
  double shad = pow(opts.sroot,2);
  //xmsq=xmsq*gevfb/(2.*x1*x2*shad);
  double fbGeV2 = 0.389379e12;
  xmsq=xmsq*fbGeV2/(2.*phasespace::x1*phasespace::x2*shad);
  
  //f[np+1]=xmsq;;

  //if (doFill.ne.0) then
  // {
  //   val=xmsq*wgt
  //   hists_fill_(p3,p4,val);
  // }

  //} end PDFs loop

  //end vcth loop
  //end dOmega loop

  if (xmsq != xmsq)
    {
      cout << "xmsq in vjloint is nan" << endl;
      return 0.;
    }
  
  return xmsq*jac;
}

//generate phase space
void vjloint::genps(const double x[6])
{
}
