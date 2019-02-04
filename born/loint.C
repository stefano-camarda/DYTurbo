#include "loint.h"
#include "mesq.h"
#include "parton.h"
#include "settings.h"
#include "omegaintegr.h"
#include "pdf.h"
#include "scales.h"
#include "phasespace.h"
#include "abint.h"
#include "resconst.h"
#include "dyres_interface.h"
#include "dynnlo_interface.h"
#include "isnan.h"

#include <iostream>
#include <string.h>

using namespace std;

//double loint::lint(double costh, double m, double y, int mode, double f[])
void loint::lint(double costh, double m, double y, int mode, double f[2])
{
  f[0] = 0.;
  f[1] = 0.;

  double m2 = m*m;
  double exppy = exp(y);
  double expmy = 1./exppy;
  double tau = sqrt(m2/pow(opts.sroot,2));
  
  mesq::setpropagators(m);

  double cthmom0, cthmom1, cthmom2;
  if (mode == 0)
    {
      //apply here lepton cuts
      cthmom0 = 1.;
      cthmom1 = costh;
      cthmom2 = pow(costh,2);
    }
  else if (mode == 1)
    {
      phasespace::set_mqtyphi(m, 0., y);
      omegaintegr::genV4p();
      omegaintegr::cthmoments(cthmom0,cthmom1,cthmom2);
    }

  mesq::setmesq(cthmom0, cthmom1, cthmom2);
  //cout << cthmom0 << "  " << cthmom1 << "  " << cthmom2 << endl;

  //calculate Bjorken x1 x2
  double x1 = tau*exppy;
  double x2 = tau*expmy;

  if (x1 >= 1 || x2 >= 1)
    return;
  
  /*
  // **************** Check matrix element calculation
  //    phase space generation
  // *****************************************
  //     generate p3 and p4 4-momenta
  //     pick whatever value of phi and phi_lep
  double phi = 2*M_PI*0.25;
  double phi_lep = 0.75*2*M_PI;
  phasespace::set_mqtyphi(m, 0., y, phi);
  phasespace::set_cth(costh);
  phasespace::set_philep(phi_lep);
  phasespace::calcexpy();
  phasespace::calcmt();
  phasespace::genV4p();
  phasespace::genl4p();

  double p3[4], p4[4];
  phasespace::getp3(p3);
  phasespace::getp4(p4);

  double p1[4];
  p1[0] = 0.;
  p1[1] = 0.;
  p1[2] = -x1*opts.sroot/2.; 
  p1[3] = -x1*opts.sroot/2.;

  double p2[4];
  p2[0] = 0.;
  p2[1] = 0.;
  p2[2] = x2*opts.sroot/2.; 
  p2[3] = -x2*opts.sroot/2.;

  double p[4][12];
  p[0][0] = p1[0];
  p[1][0] = p1[1];
  p[2][0] = p1[2];
  p[3][0] = p1[3];

  p[0][1] = p2[0];
  p[1][1] = p2[1];
  p[2][1] = p2[2];
  p[3][1] = p2[3];

  p[0][2] = p3[0];
  p[1][2] = p3[1];
  p[2][2] = p3[2];
  p[3][2] = p3[3];

  p[0][3] = p4[0];
  p[1][3] = p4[1];
  p[2][3] = p4[2];
  p[3][3] = p4[3];
  //     End of phase space generation
  //     *****************************************

  // Compute Born matrix element
  double msqc[11][11];
  if(opts.nproc == 3)
    qqb_z_(p,msqc);
  else
    qqb_w_(p,msqc);

  initsigma_cpp_(m, cthmom0, cthmom1, cthmom2);

  double fbGeV2=0.38937966e12;

  cout << m << "  " << y << "  "  << costh << endl;  
  fcomplex za[12][12];
  fcomplex zb[12][12];
  int N = 5;
  spinoru_(N,p,za,zb);
  cout << real(cx(za[2][1])*cx(zb[0][3]))/m/m << "  " << (1 + costh)/2 << endl;
  cout << real(cx(za[2][0])*cx(zb[1][3]))/m/m << "  " << (1 - costh)/2 <<  endl;
  if (opts.nproc == 3)
    {
      cout << mesq::mesqij[0]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::Ub][parton::U] << "  " << sigmaij_.sigmaij_[parton::ub][parton::u ]* M_PI*6./fbGeV2*(2*m2) << endl;//* M_PI*6./fbGeV2*(2*m2) << endl;
      cout << mesq::mesqij[1]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::U][parton::Ub] << "  " << sigmaij_.sigmaij_[parton::u ][parton::ub ]* M_PI*6./fbGeV2*(2*m2) << endl;//* M_PI*6./fbGeV2*(2*m2) << endl;
      cout << mesq::mesqij[2]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::Db][parton::D] << "  " << sigmaij_.sigmaij_[parton::db][parton::d ]* M_PI*6./fbGeV2*(2*m2) << endl;//* M_PI*6./fbGeV2*(2*m2) << endl;
      cout << mesq::mesqij[3]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::D][parton::Db] << "  " << sigmaij_.sigmaij_[parton::d ][parton::db ]* M_PI*6./fbGeV2*(2*m2) << endl;//* M_PI*6./fbGeV2*(2*m2) << endl;
      cout << mesq::mesqij[4]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::Sb][parton::S] << "  " << sigmaij_.sigmaij_[parton::sb][parton::s ]* M_PI*6./fbGeV2*(2*m2) << endl;//* M_PI*6./fbGeV2*(2*m2) << endl;
      cout << mesq::mesqij[5]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::S][parton::Sb] << "  " << sigmaij_.sigmaij_[parton::s ][parton::sb ]* M_PI*6./fbGeV2*(2*m2) << endl;//* M_PI*6./fbGeV2*(2*m2) << endl;
      cout << mesq::mesqij[6]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::Cb][parton::C] << "  " << sigmaij_.sigmaij_[parton::cb][parton::c ]* M_PI*6./fbGeV2*(2*m2) << endl;//* M_PI*6./fbGeV2*(2*m2) << endl;
      cout << mesq::mesqij[7]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::C][parton::Cb] << "  " << sigmaij_.sigmaij_[parton::c ][parton::cb ]* M_PI*6./fbGeV2*(2*m2) << endl;//* M_PI*6./fbGeV2*(2*m2) << endl;
      cout << mesq::mesqij[8]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::Bb][parton::B] << "  " << sigmaij_.sigmaij_[parton::bb][parton::b ]* M_PI*6./fbGeV2*(2*m2) << endl;//* M_PI*6./fbGeV2*(2*m2) << endl;
      cout << mesq::mesqij[9]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::B][parton::Bb] << "  " << sigmaij_.sigmaij_[parton::b ][parton::bb ]* M_PI*6./fbGeV2*(2*m2) << endl;//* M_PI*6./fbGeV2*(2*m2) << endl;
    }
  else if (opts.nproc == 1)
    {
      cout << mesq::mesqij[0 ]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::Db][parton::U ] << endl;
      cout << mesq::mesqij[1 ]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::U ][parton::Db] << endl;
      cout << mesq::mesqij[2 ]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::Sb][parton::U ] << endl;
      cout << mesq::mesqij[3 ]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::U ][parton::Sb] << endl;
      cout << mesq::mesqij[4 ]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::Bb][parton::U ] << endl;
      cout << mesq::mesqij[5 ]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::U ][parton::Bb] << endl;
      cout << mesq::mesqij[6 ]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::Sb][parton::C ] << endl;
      cout << mesq::mesqij[7 ]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::C ][parton::Sb] << endl;
      cout << mesq::mesqij[8 ]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::Db][parton::C ] << endl;
      cout << mesq::mesqij[9 ]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::C ][parton::Db] << endl;
      cout << mesq::mesqij[10]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::Bb][parton::C ] << endl;
      cout << mesq::mesqij[11]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::C ][parton::Bb] << endl;
    }
  else if (opts.nproc == 2)
    {
      cout << mesq::mesqij[0 ]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::Ub][parton::D ] << endl;
      cout << mesq::mesqij[1 ]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::D ][parton::Ub] << endl;
      cout << mesq::mesqij[2 ]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::Ub][parton::S ] << endl;
      cout << mesq::mesqij[3 ]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::S ][parton::Ub] << endl;
      cout << mesq::mesqij[4 ]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::Ub][parton::B ] << endl;
      cout << mesq::mesqij[5 ]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::B ][parton::Ub] << endl;
      cout << mesq::mesqij[6 ]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::Cb][parton::S ] << endl;
      cout << mesq::mesqij[7 ]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::S ][parton::Cb] << endl;
      cout << mesq::mesqij[8 ]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::Cb][parton::D ] << endl;
      cout << mesq::mesqij[9 ]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::D ][parton::Cb] << endl;
      cout << mesq::mesqij[10]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::Cb][parton::B ] << endl;
      cout << mesq::mesqij[11]* M_PI*6./fbGeV2*(2*m2) << "  " << msqc[parton::B ][parton::Cb] << endl;
    }
  cout << endl;

  
  //     End of ME check
  //*******************************************************
  */
  
  //Set scales
  scales::set(m);
  scales::mcfm();
  double muf = scales::fac;
  double mur = scales::ren;

  /*
  //Set factorization scale
  double muf, mur;
  if (opts.dynamicscale)
    {
      muf = m*opts.kmufac;
      mur = m*opts.kmuren;
      double mur2 = mur*mur;
      scaleset_(mur2);
    }
  else
    {
      muf = opts.rmass*opts.kmufac;
      mur = opts.rmass*opts.kmuren;
    }
  */

  //PDFs
  double fx1[2*MAXNF+1],fx2[2*MAXNF+1];
  fdist_(opts.ih1,x1,muf,fx1);
  fdist_(opts.ih2,x2,muf,fx2);

  /*
  double xmsq = 0.;
  if (opts.nproc == 3)
    {
      xmsq += real(mesq::mesqij[0])*fx1[parton::U ]*fx2[parton::Ub];
      xmsq += real(mesq::mesqij[1])*fx1[parton::Ub]*fx2[parton::U ];
      xmsq += real(mesq::mesqij[2])*fx1[parton::D ]*fx2[parton::Db];
      xmsq += real(mesq::mesqij[3])*fx1[parton::Db]*fx2[parton::D ];
      xmsq += real(mesq::mesqij[4])*fx1[parton::S ]*fx2[parton::Sb];
      xmsq += real(mesq::mesqij[5])*fx1[parton::Sb]*fx2[parton::S ];
      xmsq += real(mesq::mesqij[6])*fx1[parton::C ]*fx2[parton::Cb];
      xmsq += real(mesq::mesqij[7])*fx1[parton::Cb]*fx2[parton::C ];
      xmsq += real(mesq::mesqij[8])*fx1[parton::B ]*fx2[parton::Bb];
      xmsq += real(mesq::mesqij[9])*fx1[parton::Bb]*fx2[parton::B ];
    }
  else if (opts.nproc == 1)
    {
      xmsq += real(mesq::mesqij[0])*fx1[parton::U ]*fx2[parton::Db];
      xmsq += real(mesq::mesqij[1])*fx1[parton::Db]*fx2[parton::U ];
      xmsq += real(mesq::mesqij[2])*fx1[parton::U ]*fx2[parton::Sb];
      xmsq += real(mesq::mesqij[3])*fx1[parton::Sb]*fx2[parton::U ];
      xmsq += real(mesq::mesqij[4])*fx1[parton::U ]*fx2[parton::Bb];
      xmsq += real(mesq::mesqij[5])*fx1[parton::Bb]*fx2[parton::U ];
      xmsq += real(mesq::mesqij[6])*fx1[parton::C ]*fx2[parton::Sb];
      xmsq += real(mesq::mesqij[7])*fx1[parton::Sb]*fx2[parton::C ];
      xmsq += real(mesq::mesqij[8])*fx1[parton::C ]*fx2[parton::Db];
      xmsq += real(mesq::mesqij[9])*fx1[parton::Db]*fx2[parton::C ];
      xmsq += real(mesq::mesqij[10])*fx1[parton::C ]*fx2[parton::Bb];
      xmsq += real(mesq::mesqij[11])*fx1[parton::Bb]*fx2[parton::C ];
    }
  else if (opts.nproc == 2)
    {
      //!!! remove parenthesis!!!
      xmsq += real(mesq::mesqij[0]*fx1[parton::D ]*fx2[parton::Ub]);
      xmsq += real(mesq::mesqij[1]*fx1[parton::Ub]*fx2[parton::D ]);
      xmsq += real(mesq::mesqij[2]*fx1[parton::S ]*fx2[parton::Ub]);
      xmsq += real(mesq::mesqij[3]*fx1[parton::Ub]*fx2[parton::S ]);
      xmsq += real(mesq::mesqij[4]*fx1[parton::B ]*fx2[parton::Ub]);
      xmsq += real(mesq::mesqij[5]*fx1[parton::Ub]*fx2[parton::B ]);
      xmsq += real(mesq::mesqij[6]*fx1[parton::S ]*fx2[parton::Cb]);
      xmsq += real(mesq::mesqij[7]*fx1[parton::Cb]*fx2[parton::S ]);
      xmsq += real(mesq::mesqij[8]*fx1[parton::D ]*fx2[parton::Cb]);
      xmsq += real(mesq::mesqij[9]*fx1[parton::Cb]*fx2[parton::D ]);
      xmsq += real(mesq::mesqij[10]*fx1[parton::B ]*fx2[parton::Cb]);
      xmsq += real(mesq::mesqij[11]*fx1[parton::Cb]*fx2[parton::B ]);
    }
  */


  // start QCD
  /*****************************************************/
  double LF, LR;
  if (opts.order >= 1)
    LF = log(m2/pow(muf,2));
  if (opts.order >= 2)
    LR = log(m2/pow(mur,2));

  double lx1, lx2;
  double pqqintx1, pqqintx2;
  double d0intx1, d0intx2;
  double d1intx1, d1intx2;
  if (opts.order >= 1)
    {
      lx1 = log(x1);
      lx2 = log(x2);
      pqqintx1 = pqqint_(x1);
      pqqintx2 = pqqint_(x2);
    }
  if (opts.order >= 2)
    {
      d0intx1 = d0int_(x1);
      d0intx2 = d0int_(x2);
      d1intx1 = d1int_(x1);
      d1intx2 = d1int_(x2);
    }
 
  //Start the fast alfa beta integration

  //Preliminary loop for caching
  double cz1[abint::abdim];
  double cz2[abint::abdim];
  // --> cache the logs, and pass them as parameters to the CxxPxx, Pxxxx functions
  //double lz1[abint::abdim];
  //double lz2[abint::abdim];
  //double l1z1[abint::abdim];
  //double l1z2[abint::abdim];
  double lzoz1[abint::abdim];
  double lzoz2[abint::abdim];
  if (opts.order >= 1)
    for (int ab = 0; ab < abint::abdim; ab++)
      {
	cz1[ab] = pow(x1,abint::abx[ab]);
	cz2[ab] = pow(x2,abint::abx[ab]);
	//lz1[ab] = log(cz1[ab]);
	//lz2[ab] = log(cz2[ab]);
	//l1z1[ab] = log(1.-cz1[ab]);
	//l1z2[ab] = log(1.-cz2[ab]);
	lzoz1[ab] = log(1.-cz1[ab])/(1.-cz1[ab]);
	lzoz2[ab] = log(1.-cz2[ab])/(1.-cz2[ab]);
      }

  //cache scaled PDFs
  double fx1p[abint::abdim][2*MAXNF+1],fx2p[abint::abdim][2*MAXNF+1];
  if (opts.order >= 1)
    for (int ab = 0; ab < abint::abdim; ab++)
      {
	double fx1temp[2*MAXNF+1],fx2temp[2*MAXNF+1];
	double xx1 = pow(x1,(1-abint::abx[ab]));
	double xx2 = pow(x2,(1-abint::abx[ab]));
	fdist_(opts.ih1,xx1,muf,fx1temp);
	fdist_(opts.ih2,xx2,muf,fx2temp);
	memcpy(fx1p[ab], fx1temp, (2*MAXNF+1)*sizeof(double));
	memcpy(fx2p[ab], fx2temp, (2*MAXNF+1)*sizeof(double));
      }

  double sumfx1p[abint::abdim];
  double sumfx2p[abint::abdim];
  if (opts.order >= 2)
    for (int ab = 0; ab < abint::abdim; ab++)
      {
	sumfx1p[ab] = 0.;
	sumfx2p[ab] = 0.;
	for (int f = 0; f < MAXNF; f++)
	  {
	    sumfx1p[ab] += fx1p[ab][f]+fx1p[ab][parton::charge_conj(parton::pdgid(f))];
	    sumfx2p[ab] += fx2p[ab][f]+fx2p[ab][parton::charge_conj(parton::pdgid(f))];
	  }
      }
  

  /*
  // Start calculation
  double tdelta = 0; //LO
  double th1st = 0;  //NLO
  double th1stF = 0; //NLO
  double th2st = 0;  //NNLO
  double tcga = 0;   //NNLO
  double tgamma2 = 0;//NNLO
  double tgaga = 0;  //NNLO

  //loop on born subprocesses, i.e. born incoming partons ij
  for (int sp = 0; sp < mesq::totpch; sp++)
    {
      //simplify notation
      double bornmesqij = real(mesq::mesqij[sp]); //born level amplitudes
      parton::pdgid i = mesq::pid1[sp];         //parton 1
      parton::pdgid j = mesq::pid2[sp];         //parton 2
      parton::pdgid g = parton::G;              //gluon
      parton::pdgid im = parton::charge_conj(i);
      parton::pdgid jm = parton::charge_conj(j);

      //LO
      //Simplest term without convolutions
      tdelta += fx1[i]*fx2[j]*bornmesqij;
      if (opts.order == 0) continue;

      //NLO
      //H1st delta term
      th1st += 2*resconst::C1qqdelta*fx1[i]*fx2[j]*bornmesqij;
      
      //alfa loop (first leg)
      for (int a = 0; a < abint::abdim; a++)
	{
	  if (cz1[a] >= 1) continue;
	  
	  //H1st non delta terms
	  th1st += (fx1p[a][i]*cqq_(cz1[a])+fx1p[a][g]*cqg_(cz1[a]))*(-lx1)*fx2[j]*bornmesqij * abint::abw[a];
	  
	  //H1st muf dependence, gammaqq and gammaqg:
	  th1stF += (-lx1*((fx1p[a][i]-fx1[i]*cz1[a])*pqq_(cz1[a])+fx1p[a][g]*dypqg_(cz1[a])))*fx2[j]*bornmesqij * abint::abw[a];
	}
      th1stF += -pqqintx1*fx1[i]*fx2[j]*bornmesqij;
      
      //beta loop (second leg)
      for (int b = 0; b < abint::abdim; b++)
	{
	  if (cz2[b] >= 1) continue;

	  //H1st non delta terms
	  th1st += (fx2p[b][j]*cqq_(cz2[b])+fx2p[b][g]*cqg_(cz2[b]))*(-lx2)*fx1[i]*bornmesqij * abint::abw[b];
	  
	  //H1st muf dependence, gammaqq and gammaqg:
	  th1stF += (-lx2*((fx2p[b][j]-fx2[j]*cz2[b])*pqq_(cz2[b])+fx2p[b][g]*dypqg_(cz2[b])))*fx1[i]*bornmesqij * abint::abw[b];
	}
      th1stF += -pqqintx2*fx1[i]*fx2[j]*bornmesqij;
      if (opts.order == 1) continue;

      //NNLO
      //alfa loop
      double diffg1f = 0;
      double diffg10 = 0;
      double diff1 = 0;
      double diffc1f = 0;
      double diffc10 = 0;
      for (int a = 0; a < abint::abdim; a++)
	{
	  if (cz1[a] >= 1) continue;

	  //h2st qqbar contribution from c1*c1 (without delta term) -> regular-delta
	  th2st += fx1p[a][i]*cqq_(cz1[a])*(-lx1)*fx2[j]*resconst::C1qqdelta*bornmesqij * abint::abw[a];

	  //h2st qg contribution from c1*c1 ->regular-delta
	  th2st += fx1p[a][g]*cqg_(cz1[a])*(-lx1)*fx2[j]*resconst::C1qqdelta*bornmesqij * abint::abw[a];

	  //H2st qqbar channel: D0(z), first leg
	  th2st += 0.5*(-lx1*(fx1p[a][i]-fx1[i]*cz1[a])*resconst::H2qqD0/(1.-cz1[a]))*fx2[j]*bornmesqij * abint::abw[a];
	  th2st += -0.5*resconst::H2qqD0*d0intx1*fx1[i]*fx2[j]*bornmesqij * abint::abw[a];

	  //C2qq, regular part, first leg
	  th2st += fx1p[a][i]*c2qqreg_(cz1[a])*(-lx1)*fx2[j]*bornmesqij * abint::abw[a];

	  //C2qg, first leg
	  th2st += fx1p[a][g]*c2qg_(cz1[a])*(-lx1)*fx2[j]*bornmesqij * abint::abw[a];

	  //Cqqbar contribution: first leg
	  th2st += fx1p[a][im]*c2qqb_(cz1[a])*(-lx1)*fx2[j]*bornmesqij * abint::abw[a];

	  //Cqqp contribution: first leg
	  th2st += (sumfx1p[a]-(fx1p[a][i]+fx1p[a][parton::charge_conj(parton::pdgid(i))]))*c2qqp_(cz1[a])*(-lx1)*fx2[j]*bornmesqij * abint::abw[a];

	  //Terms needed for NNLO scale dependence
	  //(gamma+gamma)*(gamma+gamma) term

	  //First part: one gamma for each leg
	  diffg1f += (-lx1*(fx1p[a][i]-fx1[i]*cz1[a])*pqq_(cz1[a]) - pqqintx1*fx1[i]) * abint::abw[a];
	  diffg10 += -lx1*fx1p[a][g]*dypqg_(cz1[a]) * abint::abw[a];

	  //Second part: gamma*gamma terms
	  //Pij * Pjk = D1ijjk (log(1-z)/(1-z))_+ + D0ijjk/(1-z)_+ 
	  //          + Pijjk(z) + Deltaijjk delta(1-z)
	  //First leg
	  diff1 += (-lx1*((fx1p[a][i]-fx1[i]*cz1[a])
			  *(resconst::D0qqqq/(1-cz1[a])+resconst::D1qqqq*lzoz1[a])
			  +fx1p[a][i]*pqqqq_(cz1[a])+fx1p[a][g]*(pqqqg_(cz1[a])+pqggg_(cz1[a])))
		    +(resconst::Deltaqqqq-resconst::D0qqqq*d0intx1-resconst::D1qqqq*d1intx1)
		    *fx1[i]) * abint::abw[a];

	  //Include Pqggq
	  diff1 += -lx1*sumfx1p[a]*pqggq_(cz1[a]) * abint::abw[a];
	  //End of (gamma+gamma)*(gamma+gamma) term

	  //Start  (C+C)*(gamma+gamma) term
	  //C first leg, gamma second leg
	  diffc1f += (-lx1*fx1p[a][i]*cqq_(cz1[a])+resconst::C1qqdelta*fx1[i]) * abint::abw[a];
	  diffc10 += -lx1*fx1p[a][g]*cqg_(cz1[a]) * abint::abw[a];
      
	  //C*gamma: first leg (ignore delta term in Cqq: taken into account with th1stF)
	  tcga += (fx1p[a][i]*cqqpqq_(cz1[a])+fx1p[a][g]*(cqqpqg_(cz1[a])+cqgpgg_(cz1[a])))*(-lx1)*fx2[j]*bornmesqij * abint::abw[a];

	  //Add Cqg*Pgq contribution
	  tcga += sumfx1p[a]*cqgpgq_(cz1[a])*(-lx1)*fx2[j]*bornmesqij * abint::abw[a];

	  //Start 2-loop AP
	  // Gluon + pure singlet
	  //f == gluon piece
	  tgamma2 += fx1p[a][g]*p2qg_(cz1[a])*(-lx1)*fx2[j]*bornmesqij * abint::abw[a];

	  //f != gluon piece
	  tgamma2 += sumfx1p[a]*p2qqs_(cz1[a])*(-lx1)*fx2[j]*bornmesqij * abint::abw[a];

	  //P2qq non-singlet: regular part
	  tgamma2 += fx1p[a][i]*p2qqv_(cz1[a])*(-lx1)*fx2[j]*bornmesqij * abint::abw[a];

	  //P2qq non-singlet: 1/(1-z)_+
	  tgamma2 += 2./3.*resconst::Kappa*(-lx1*(fx1p[a][i]-fx1[i]*cz1[a])/(1-cz1[a])-d0intx1*fx1[i])*fx2[j]*bornmesqij * abint::abw[a];

	  //P2qqb non singlet
	  tgamma2 += fx1p[a][im]*p2qqbv_(cz1[a])*(-lx1)*fx2[j]*bornmesqij * abint::abw[a];
	}
	  
      //beta loop
      double diffg2f = 0;
      double diffg20 = 0;
      double diff2 = 0;
      double diffc2f = 0;
      double diffc20 = 0;
      for (int b = 0; b < abint::abdim; b++)
	{
	  if (cz2[b] >= 1) continue;
	  
	  //h2st qqbar contribution from c1*c1 (without delta term) -> regular-delta
	  th2st += fx2p[b][j]*cqq_(cz2[b])*(-lx2)*fx1[i]*resconst::C1qqdelta*bornmesqij * abint::abw[b];

	  //h2st qg contribution from c1*c1 ->regular-delta
	  th2st += fx2p[b][g]*cqg_(cz2[b])*(-lx2)*fx1[i]*resconst::C1qqdelta*bornmesqij * abint::abw[b];

	  //H2st, qqbar channel: D0(z), second leg
	  th2st += 0.5*(-lx2*(fx2p[b][j]-fx2[j]*cz2[b])*resconst::H2qqD0/(1.-cz2[b]))*fx1[i]*bornmesqij * abint::abw[b];
	  th2st += -0.5*resconst::H2qqD0*d0intx2*fx1[i]*fx2[j]*bornmesqij * abint::abw[b];

	  //C2qq, regular part, second leg
	  th2st += fx2p[b][j]*c2qqreg_(cz2[b])*(-lx2)*fx1[i]*bornmesqij * abint::abw[b];

	  //C2qg, second leg
	  th2st += fx2p[b][g]*c2qg_(cz2[b])*(-lx2)*fx1[i]*bornmesqij * abint::abw[b];

	  //Cqqbar contribution: second leg
	  th2st += fx2p[b][jm]*c2qqb_(cz2[b])*(-lx2)*fx1[i]*bornmesqij * abint::abw[b];

	  //Cqqp contribution: second leg
	  th2st += (sumfx2p[b]-(fx2p[b][j]+fx2p[b][parton::charge_conj(parton::pdgid(j))]))*c2qqp_(cz2[b])*(-lx2)*fx1[i]*bornmesqij * abint::abw[b];

	  //Terms needed for NNLO scale dependence
	  //(gamma+gamma)*(gamma+gamma) term

	  //First part: one gamma for each leg
	  diffg2f += (-lx2*(fx2p[b][j]-fx2[j]*cz2[b])*pqq_(cz2[b]) - pqqintx2*fx2[j]) * abint::abw[b];
	  diffg20 += -lx2*fx2p[b][g]*dypqg_(cz2[b]) * abint::abw[b];

	  //Second part: gamma*gamma terms
	  //Pij * Pjk = D1ijjk (log(1-z)/(1-z))_+ + D0ijjk/(1-z)_+ 
	  //          + Pijjk(z) + Deltaijjk delta(1-z)
	  //Second leg
	  diff2 += (-lx2*((fx2p[b][j]-fx2[j]*cz2[b])
			  *(resconst::D0qqqq/(1-cz2[b])+resconst::D1qqqq*lzoz2[b])
			  +fx2p[b][j]*pqqqq_(cz2[b])+fx2p[b][g]*(pqqqg_(cz2[b])+pqggg_(cz2[b])))
		    +(resconst::Deltaqqqq-resconst::D0qqqq*d0intx2-resconst::D1qqqq*d1intx2)
		    *fx2[j])* abint::abw[b];
	  //Include Pqggq
	  diff2 += -lx2*sumfx2p[b]*pqggq_(cz2[b]) * abint::abw[b];
	  //End of (gamma+gamma)*(gamma+gamma) term

	  //Start  (C+C)*(gamma+gamma) term
	  //gamma first leg, C second leg
	  diffc2f += (-lx2*fx2p[b][j]*cqq_(cz2[b])+resconst::C1qqdelta*fx2[j]) * abint::abw[b];
	  diffc20 += -lx2*fx2p[b][g]*cqg_(cz2[b]) * abint::abw[b];

	  //C*gamma: second leg (ignore delta term in Cqq: taken into account with th1stF)
	  tcga += (fx2p[b][j]*cqqpqq_(cz2[b])+fx2p[b][g]*(cqqpqg_(cz2[b])+cqgpgg_(cz2[b])))*(-lx2)*fx1[i]*bornmesqij * abint::abw[b];

	  //Add Cqg*Pgq contribution
	  tcga += sumfx2p[b]*cqgpgq_(cz2[b])*(-lx2)*fx1[i]*bornmesqij * abint::abw[b];

	  //Start 2-loop AP
	  // Gluon + pure singlet
	  //f == gluon piece
	  tgamma2 += fx2p[b][g]*p2qg_(cz2[b])*(-lx2)*fx1[i]*bornmesqij * abint::abw[b];

	  //f != gluon piece
	  tgamma2 += sumfx2p[b]*p2qqs_(cz2[b])*(-lx2)*fx1[i]*bornmesqij * abint::abw[b];

	  //P2qq non-singlet: regular part
	  tgamma2 += fx2p[b][j]*p2qqv_(cz2[b])*(-lx2)*fx1[i]*bornmesqij * abint::abw[b];

	  //P2qq non-singlet: 1/(1-z)_+
	  tgamma2 += 2./3.*resconst::Kappa*(-lx2*(fx2p[b][j]-fx2[j]*cz2[b])/(1-cz2[b])-d0intx2*fx2[j])*fx1[i]*bornmesqij * abint::abw[b];

	  //P2qqb non singlet
	  tgamma2 += fx2p[b][jm]*p2qqbv_(cz2[b])*(-lx2)*fx1[i]*bornmesqij * abint::abw[b];
	}

      tgaga=tgaga+2*(diffg10*diffg20+diffg1f*diffg2f+diffg10*diffg2f+diffg1f*diffg20)*bornmesqij;
      tgaga += diff1*fx2[j]*bornmesqij;
      tgaga += diff2*fx1[i]*bornmesqij;

      tcga += bornmesqij*(diffc10*diffg20+diffc1f*diffg2f+diffc10*diffg2f+diffc1f*diffg20);
      tcga += bornmesqij*(diffg10*diffc20+diffg1f*diffc2f+diffg10*diffc2f+diffg1f*diffc20);

      //alfa-beta loop
      for (int a = 0; a < abint::abdim; a++)
	for (int b = 0; b < abint::abdim; b++)
	  {
	    //h2st gg contribution
	    th2st += fx1p[a][g]*cqg_(cz1[a])*(-lx1)*fx2p[b][g]*cqg_(cz2[b])*(-lx2)*bornmesqij * abint::abw[a]*abint::abw[b];

	    //h2st qqbar contribution from c1*c1 (without delta term) -> regular*regular
	    th2st += fx1p[a][i]*cqq_(cz1[a])*(-lx1)*fx2p[b][j]*cqq_(cz2[b])*(-lx2)*bornmesqij * abint::abw[a]*abint::abw[b];

	    //h2st qg contribution from c1*c1 -> regular*regular
	    th2st += fx1p[a][g]*cqg_(cz1[a])*(-lx1)*fx2p[b][j]*cqq_(cz2[b])*(-lx2)*bornmesqij * abint::abw[a]*abint::abw[b];
	    th2st += fx1p[a][i]*cqq_(cz1[a])*(-lx1)*fx2p[b][g]*cqg_(cz2[b])*(-lx2)*bornmesqij * abint::abw[a]*abint::abw[b];
	  }
    }

  double asopi = qcdcouple_.ason2pi_*2.;
  
  //LO term
  xmsq=tdelta;

  //NLO terms
  if (opts.order >= 1)
    xmsq += asopi*(th1st+LF*th1stF);

  //NNLO terms
  if (opts.order >= 2)
    {
      xmsq += pow(asopi,2)*(tdelta*resconst::H2qqdelta+th2st);

      //add scale dependence at NNLO
      xmsq += pow(asopi,2)*(0.5*resconst::beta0*pow(LF,2)*th1stF
			    +tgamma2*LF
			    -resconst::beta0*LR*(th1st+LF*th1stF)
			    +LF*tcga+0.5*pow(LF,2)*tgaga);


      //Include missing delta term from C*gamma (no factor 2 here !)
      xmsq += pow(asopi,2)*(LF*resconst::C1qqdelta*th1stF);

      //Include missing term from contact term in 2 loop AP
      xmsq += pow(asopi,2)*(2*resconst::Delta2qq*tdelta)*LF;
    }
  */
  /****************************/

  // Start calculation
  double asopi = qcdcouple_.ason2pi_*2.;
  
  //loop on born subprocesses, i.e. born incoming partons ij
  double lumi[mesq::totpch];
  for (int sp = 0; sp < mesq::totpch; sp++)
    {
      //simplify notation
      //double bornmesqij = real(mesq::mesqij[sp]); //born level amplitudes
      parton::pdgid i = mesq::pid1[sp];         //parton 1
      parton::pdgid j = mesq::pid2[sp];         //parton 2
      parton::pdgid g = parton::G;              //gluon
      parton::pdgid im = parton::charge_conj(i);
      parton::pdgid jm = parton::charge_conj(j);

      //LO term
      //Simplest term without convolutions
      double tdelta = fx1[i]*fx2[j];

      lumi[sp] = tdelta;

      if (opts.order == 0) continue;
      //end LO
      
      //NLO terms
      double th1st = 0;
      double th1stF = 0;
      //H1st delta term
      th1st += 2*resconst::C1qqdelta*fx1[i]*fx2[j];
      
      //alfa loop (first leg)
      for (int a = 0; a < abint::abdim; a++)
	{
	  if (cz1[a] >= 1) continue;
	  
	  //H1st non delta terms
	  th1st += (fx1p[a][i]*cqq_(cz1[a])+fx1p[a][g]*cqg_(cz1[a]))*(-lx1)*fx2[j] * abint::abw[a];
	  
	  //H1st muf dependence, gammaqq and gammaqg:
	  th1stF += (-lx1*((fx1p[a][i]-fx1[i]*cz1[a])*pqq_(cz1[a])+fx1p[a][g]*dypqg_(cz1[a])))*fx2[j] * abint::abw[a];
	}
      th1stF += -pqqintx1*fx1[i]*fx2[j];
      
      //beta loop (second leg)
      for (int b = 0; b < abint::abdim; b++)
	{
	  if (cz2[b] >= 1) continue;

	  //H1st non delta terms
	  th1st += (fx2p[b][j]*cqq_(cz2[b])+fx2p[b][g]*cqg_(cz2[b]))*(-lx2)*fx1[i] * abint::abw[b];
	  
	  //H1st muf dependence, gammaqq and gammaqg:
	  th1stF += (-lx2*((fx2p[b][j]-fx2[j]*cz2[b])*pqq_(cz2[b])+fx2p[b][g]*dypqg_(cz2[b])))*fx1[i] * abint::abw[b];
	}
      th1stF += -pqqintx2*fx1[i]*fx2[j];

      lumi[sp] += asopi*(th1st+LF*th1stF);

      if (opts.order == 1) continue;
      //end NLO

      //NNLO terms
      double th2st = 0;
      double tcga = 0;
      double tgamma2 = 0;
      double tgaga = 0;
      
      //alfa loop
      double diffg1f = 0;
      double diffg10 = 0;
      double diff1 = 0;
      double diffc1f = 0;
      double diffc10 = 0;
      for (int a = 0; a < abint::abdim; a++)
	{
	  if (cz1[a] >= 1) continue;

	  //h2st qqbar contribution from c1*c1 (without delta term) -> regular-delta
	  th2st += fx1p[a][i]*cqq_(cz1[a])*(-lx1)*fx2[j]*resconst::C1qqdelta * abint::abw[a];

	  //h2st qg contribution from c1*c1 ->regular-delta
	  th2st += fx1p[a][g]*cqg_(cz1[a])*(-lx1)*fx2[j]*resconst::C1qqdelta * abint::abw[a];

	  //H2st qqbar channel: D0(z), first leg
	  th2st += 0.5*(-lx1*(fx1p[a][i]-fx1[i]*cz1[a])*resconst::H2qqD0/(1.-cz1[a]))*fx2[j] * abint::abw[a];
	  th2st += -0.5*resconst::H2qqD0*d0intx1*fx1[i]*fx2[j] * abint::abw[a];

	  //C2qq, regular part, first leg
	  th2st += fx1p[a][i]*c2qqreg_(cz1[a])*(-lx1)*fx2[j] * abint::abw[a];

	  //C2qg, first leg
	  th2st += fx1p[a][g]*c2qg_(cz1[a])*(-lx1)*fx2[j] * abint::abw[a];

	  //Cqqbar contribution: first leg
	  th2st += fx1p[a][im]*c2qqb_(cz1[a])*(-lx1)*fx2[j] * abint::abw[a];

	  //Cqqp contribution: first leg
	  th2st += (sumfx1p[a]-(fx1p[a][i]+fx1p[a][parton::charge_conj(parton::pdgid(i))]))*c2qqp_(cz1[a])*(-lx1)*fx2[j] * abint::abw[a];

	  //Terms needed for NNLO scale dependence
	  //(gamma+gamma)*(gamma+gamma) term

	  //First part: one gamma for each leg
	  diffg1f += (-lx1*(fx1p[a][i]-fx1[i]*cz1[a])*pqq_(cz1[a]) - pqqintx1*fx1[i]) * abint::abw[a];
	  diffg10 += -lx1*fx1p[a][g]*dypqg_(cz1[a]) * abint::abw[a];

	  //Second part: gamma*gamma terms
	  //Pij * Pjk = D1ijjk (log(1-z)/(1-z))_+ + D0ijjk/(1-z)_+ 
	  //          + Pijjk(z) + Deltaijjk delta(1-z)
	  //First leg
	  diff1 += (-lx1*((fx1p[a][i]-fx1[i]*cz1[a])
			  *(resconst::D0qqqq/(1-cz1[a])+resconst::D1qqqq*lzoz1[a])
			  +fx1p[a][i]*pqqqq_(cz1[a])+fx1p[a][g]*(pqqqg_(cz1[a])+pqggg_(cz1[a])))
		    +(resconst::Deltaqqqq-resconst::D0qqqq*d0intx1-resconst::D1qqqq*d1intx1)
		    *fx1[i]) * abint::abw[a];

	  //Include Pqggq
	  diff1 += -lx1*sumfx1p[a]*pqggq_(cz1[a]) * abint::abw[a];
	  //End of (gamma+gamma)*(gamma+gamma) term

	  //Start  (C+C)*(gamma+gamma) term
	  //C first leg, gamma second leg
	  diffc1f += (-lx1*fx1p[a][i]*cqq_(cz1[a])+resconst::C1qqdelta*fx1[i]) * abint::abw[a];
	  diffc10 += -lx1*fx1p[a][g]*cqg_(cz1[a]) * abint::abw[a];
      
	  //C*gamma: first leg (ignore delta term in Cqq: taken into account with th1stF)
	  tcga += (fx1p[a][i]*cqqpqq_(cz1[a])+fx1p[a][g]*(cqqpqg_(cz1[a])+cqgpgg_(cz1[a])))*(-lx1)*fx2[j] * abint::abw[a];

	  //Add Cqg*Pgq contribution
	  tcga += sumfx1p[a]*cqgpgq_(cz1[a])*(-lx1)*fx2[j] * abint::abw[a];

	  //Start 2-loop AP
	  // Gluon + pure singlet
	  //f == gluon piece
	  tgamma2 += fx1p[a][g]*p2qg_(cz1[a])*(-lx1)*fx2[j] * abint::abw[a];

	  //f != gluon piece
	  tgamma2 += sumfx1p[a]*p2qqs_(cz1[a])*(-lx1)*fx2[j] * abint::abw[a];

	  //P2qq non-singlet: regular part
	  tgamma2 += fx1p[a][i]*p2qqv_(cz1[a])*(-lx1)*fx2[j] * abint::abw[a];

	  //P2qq non-singlet: 1/(1-z)_+
	  tgamma2 += 2./3.*resconst::Kappa*(-lx1*(fx1p[a][i]-fx1[i]*cz1[a])/(1-cz1[a])-d0intx1*fx1[i])*fx2[j] * abint::abw[a];

	  //P2qqb non singlet
	  tgamma2 += fx1p[a][im]*p2qqbv_(cz1[a])*(-lx1)*fx2[j] * abint::abw[a];
	}
	  
      //beta loop
      double diffg2f = 0;
      double diffg20 = 0;
      double diff2 = 0;
      double diffc2f = 0;
      double diffc20 = 0;
      for (int b = 0; b < abint::abdim; b++)
	{
	  if (cz2[b] >= 1) continue;
	  
	  //h2st qqbar contribution from c1*c1 (without delta term) -> regular-delta
	  th2st += fx2p[b][j]*cqq_(cz2[b])*(-lx2)*fx1[i]*resconst::C1qqdelta * abint::abw[b];

	  //h2st qg contribution from c1*c1 ->regular-delta
	  th2st += fx2p[b][g]*cqg_(cz2[b])*(-lx2)*fx1[i]*resconst::C1qqdelta * abint::abw[b];

	  //H2st, qqbar channel: D0(z), second leg
	  th2st += 0.5*(-lx2*(fx2p[b][j]-fx2[j]*cz2[b])*resconst::H2qqD0/(1.-cz2[b]))*fx1[i] * abint::abw[b];
	  th2st += -0.5*resconst::H2qqD0*d0intx2*fx1[i]*fx2[j] * abint::abw[b];

	  //C2qq, regular part, second leg
	  th2st += fx2p[b][j]*c2qqreg_(cz2[b])*(-lx2)*fx1[i] * abint::abw[b];

	  //C2qg, second leg
	  th2st += fx2p[b][g]*c2qg_(cz2[b])*(-lx2)*fx1[i] * abint::abw[b];

	  //Cqqbar contribution: second leg
	  th2st += fx2p[b][jm]*c2qqb_(cz2[b])*(-lx2)*fx1[i] * abint::abw[b];

	  //Cqqp contribution: second leg
	  th2st += (sumfx2p[b]-(fx2p[b][j]+fx2p[b][parton::charge_conj(parton::pdgid(j))]))*c2qqp_(cz2[b])*(-lx2)*fx1[i] * abint::abw[b];

	  //Terms needed for NNLO scale dependence
	  //(gamma+gamma)*(gamma+gamma) term

	  //First part: one gamma for each leg
	  diffg2f += (-lx2*(fx2p[b][j]-fx2[j]*cz2[b])*pqq_(cz2[b]) - pqqintx2*fx2[j]) * abint::abw[b];
	  diffg20 += -lx2*fx2p[b][g]*dypqg_(cz2[b]) * abint::abw[b];

	  //Second part: gamma*gamma terms
	  //Pij * Pjk = D1ijjk (log(1-z)/(1-z))_+ + D0ijjk/(1-z)_+ 
	  //          + Pijjk(z) + Deltaijjk delta(1-z)
	  //Second leg
	  diff2 += (-lx2*((fx2p[b][j]-fx2[j]*cz2[b])
			  *(resconst::D0qqqq/(1-cz2[b])+resconst::D1qqqq*lzoz2[b])
			  +fx2p[b][j]*pqqqq_(cz2[b])+fx2p[b][g]*(pqqqg_(cz2[b])+pqggg_(cz2[b])))
		    +(resconst::Deltaqqqq-resconst::D0qqqq*d0intx2-resconst::D1qqqq*d1intx2)
		    *fx2[j])* abint::abw[b];
	  //Include Pqggq
	  diff2 += -lx2*sumfx2p[b]*pqggq_(cz2[b]) * abint::abw[b];
	  //End of (gamma+gamma)*(gamma+gamma) term

	  //Start  (C+C)*(gamma+gamma) term
	  //gamma first leg, C second leg
	  diffc2f += (-lx2*fx2p[b][j]*cqq_(cz2[b])+resconst::C1qqdelta*fx2[j]) * abint::abw[b];
	  diffc20 += -lx2*fx2p[b][g]*cqg_(cz2[b]) * abint::abw[b];

	  //C*gamma: second leg (ignore delta term in Cqq: taken into account with th1stF)
	  tcga += (fx2p[b][j]*cqqpqq_(cz2[b])+fx2p[b][g]*(cqqpqg_(cz2[b])+cqgpgg_(cz2[b])))*(-lx2)*fx1[i] * abint::abw[b];

	  //Add Cqg*Pgq contribution
	  tcga += sumfx2p[b]*cqgpgq_(cz2[b])*(-lx2)*fx1[i] * abint::abw[b];

	  //Start 2-loop AP
	  // Gluon + pure singlet
	  //f == gluon piece
	  tgamma2 += fx2p[b][g]*p2qg_(cz2[b])*(-lx2)*fx1[i] * abint::abw[b];

	  //f != gluon piece
	  tgamma2 += sumfx2p[b]*p2qqs_(cz2[b])*(-lx2)*fx1[i] * abint::abw[b];

	  //P2qq non-singlet: regular part
	  tgamma2 += fx2p[b][j]*p2qqv_(cz2[b])*(-lx2)*fx1[i] * abint::abw[b];

	  //P2qq non-singlet: 1/(1-z)_+
	  tgamma2 += 2./3.*resconst::Kappa*(-lx2*(fx2p[b][j]-fx2[j]*cz2[b])/(1-cz2[b])-d0intx2*fx2[j])*fx1[i] * abint::abw[b];

	  //P2qqb non singlet
	  tgamma2 += fx2p[b][jm]*p2qqbv_(cz2[b])*(-lx2)*fx1[i] * abint::abw[b];
	}

      tgaga=tgaga+2*(diffg10*diffg20+diffg1f*diffg2f+diffg10*diffg2f+diffg1f*diffg20);
      tgaga += diff1*fx2[j];
      tgaga += diff2*fx1[i];

      tcga += (diffc10*diffg20+diffc1f*diffg2f+diffc10*diffg2f+diffc1f*diffg20);
      tcga += (diffg10*diffc20+diffg1f*diffc2f+diffg10*diffc2f+diffg1f*diffc20);

      //alfa-beta loop
      for (int a = 0; a < abint::abdim; a++)
	for (int b = 0; b < abint::abdim; b++)
	  {
	    //h2st gg contribution
	    th2st += fx1p[a][g]*cqg_(cz1[a])*(-lx1)*fx2p[b][g]*cqg_(cz2[b])*(-lx2) * abint::abw[a]*abint::abw[b];

	    //h2st qqbar contribution from c1*c1 (without delta term) -> regular*regular
	    th2st += fx1p[a][i]*cqq_(cz1[a])*(-lx1)*fx2p[b][j]*cqq_(cz2[b])*(-lx2) * abint::abw[a]*abint::abw[b];

	    //h2st qg contribution from c1*c1 -> regular*regular
	    th2st += fx1p[a][g]*cqg_(cz1[a])*(-lx1)*fx2p[b][j]*cqq_(cz2[b])*(-lx2) * abint::abw[a]*abint::abw[b];
	    th2st += fx1p[a][i]*cqq_(cz1[a])*(-lx1)*fx2p[b][g]*cqg_(cz2[b])*(-lx2) * abint::abw[a]*abint::abw[b];
	  }

      //NNLO terms
      lumi[sp] += pow(asopi,2)*(tdelta*resconst::H2qqdelta+th2st);
      
      //add scale dependence at NNLO
      lumi[sp] += pow(asopi,2)*(0.5*resconst::beta0*pow(LF,2)*th1stF
				+tgamma2*LF
				-resconst::beta0*LR*(th1st+LF*th1stF)
				+LF*tcga+0.5*pow(LF,2)*tgaga);
      
      //Include missing delta term from C*gamma (no factor 2 here !)
      lumi[sp] += pow(asopi,2)*(LF*resconst::C1qqdelta*th1stF);
      
      //Include missing term from contact term in 2 loop AP
      lumi[sp] += pow(asopi,2)*(2*resconst::Delta2qq*tdelta)*LF;
    }
  /****************************/


  double xmsq = 0.;
  for (int sp = 0; sp < mesq::totpch; sp++)
    xmsq += lumi[sp]*real(mesq::mesqij[sp]);
  
  double xmsq_m = 0.;
  if (opts.helicity >= 0)
    {
      if (opts.helicity == 0)
	{
	  cthmom0 = 4./3.; //20./3. * (0.5*2. - 1.5*2./3.  ) +2./3.*2.;
	  cthmom1 = 0.;
	  cthmom2 = -4./3.; //20./3. * (0.5*2./3.-1.5*2./5.  ) +2./3.*2./3.;
	}
      else if (opts.helicity == 4)
	{
	  cthmom0=0.;
	  cthmom1=4*cthmom2*(y < 0 ? -1. : 1.);
	  cthmom2=0.;
	}
      else
	{
	  cthmom0=0.;
	  cthmom1=0.;
	  cthmom2=0.;
	}
      mesq::setmesq(cthmom0, cthmom1, cthmom2);
      for (int sp = 0; sp < mesq::totpch; sp++)
	xmsq_m += lumi[sp]*real(mesq::mesqij[sp]);
      //xmsq_m = convolute(fx1,fx2,mesq::mesqij);
    }
  
  double shad = pow(opts.sroot,2);

  //double fbGeV2=0.38937966e12;
  //double fac = M_PI*6./fbGeV2*(2*m2);
  //double flux = fbGeV2/(2.*x1*x2*shad);
  //double ps = 1./shad;
  //double norm = 1./16./M_PI / 2. / M_PI;
  //xmsq = xmsq * ps * fac * flux * norm;

  xmsq = xmsq/shad * 3./8. /2./M_PI;
  xmsq_m = xmsq_m/shad * 3./8. /2./M_PI;

  if (isnan_ofast(xmsq))
    cout << m << " " << y << " " << costh << "  " << xmsq << endl;

  f[0] = xmsq;
  f[1] = xmsq_m;
  
  return;
}

double loint::convolute(double fx1[2*MAXNF+1], double fx2[2*MAXNF+1], complex <double> mesq[12])
{
  double xmsq = 0.;

  for (int sp = 0; sp < mesq::totpch; sp++)
    xmsq += fx1[mesq::pid1[sp]]*fx2[mesq::pid2[sp]]*real(mesq[sp]);

  return xmsq;
}
