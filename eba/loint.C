#include "loint.h"
#include "mesq.h"
#include "parton.h"
#include "settings.h"
#include "omegaintegr.h"
#include "pdf.h"
#include "phasespace.h"

#include <iostream>

using namespace std;

double loint::lint(double costh, double m, double y, int mode)
{
  double m2 = m*m;
  double exppy = exp(y);
  double expmy = 1./exppy;

  if (y >= 0)
    phasespace::setcthbounds(opts.costhmin,opts.costhmax);
  else
    phasespace::setcthbounds(-opts.costhmax,-opts.costhmin);
  
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

  cthmom0=cthmom0;
  cthmom1=cthmom1;
  cthmom2=cthmom2;
  mesq::setmesq(cthmom0, cthmom1, cthmom2);

  //calculate Bjorken x1 x2
  double x1 = sqrt(m2/pow(opts.sroot,2))*exppy;
  double x2 = sqrt(m2/pow(opts.sroot,2))*expmy;

  /*
  // **************** Check matrix element calculation
  //    phase space generation
  // *****************************************
  //     generate p3 and p4 4-momenta
  //     pick whatever value of phi and phi_lep
  double phi = 2*M_PI*0.25;
  double phi_lep = 0.75*2*M_PI;
  phasespace::set_mqtycth(m, 0., y, costh);
  genV4p(m, 0., y, phi);
  genl4p(costh, phi_lep);

  double p3[4], p4[4];
  getp3(p3);
  getp4(p4);

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
  */
  /*
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
  cout << endl;
  */
//     End of ME check
//*******************************************************

  
  
  
  //Set factorization scale
  double muf;
  if (opts.dynamicscale)
    muf = m*opts.kmufac;
  else
    muf = opts.rmass*opts.kmufac;
  
  double fx1[2*MAXNF+1],fx2[2*MAXNF+1];
  fdist_(opts.ih1,x1,muf,fx1);
  fdist_(opts.ih2,x2,muf,fx2);
  
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
      xmsq += real(mesq::mesqij[0]*fx1[parton::D ]*fx2[parton::Ub]);
      xmsq += real(mesq::mesqij[1]*fx1[parton::Ub]*fx2[parton::D ]);
      xmsq += real(mesq::mesqij[2]*fx1[parton::S ]*fx2[parton::Ub]);
      xmsq += real(mesq::mesqij[3]*fx1[parton::Ub]*fx2[parton::S ]);
      xmsq += real(mesq::mesqij[4]*fx1[parton::B ]*fx2[parton::Ub]);
      xmsq += real(mesq::mesqij[5]*fx1[parton::Ub]*fx2[parton::B ]);
      xmsq += real(mesq::mesqij[6]*fx1[parton::C ]*fx2[parton::Sb]);
      xmsq += real(mesq::mesqij[7]*fx1[parton::Sb]*fx2[parton::C ]);
      xmsq += real(mesq::mesqij[8]*fx1[parton::D ]*fx2[parton::Cb]);
      xmsq += real(mesq::mesqij[9]*fx1[parton::Cb]*fx2[parton::D ]);
      xmsq += real(mesq::mesqij[10]*fx1[parton::B ]*fx2[parton::Cb]);
      xmsq += real(mesq::mesqij[11]*fx1[parton::Cb]*fx2[parton::B ]);
    }

  double shad = pow(opts.sroot,2);
  

  //double fbGeV2=0.38937966e12;
  //double fac = M_PI*6./fbGeV2*(2*m2);
  //double flux = fbGeV2/(2.*x1*x2*shad);
  //double ps = 1./shad;
  //double norm = 1./16./M_PI / 2. / M_PI;
  //xmsq = xmsq * ps * fac * flux * norm;

  xmsq = xmsq/shad * 3./8. /2./M_PI;
  
  //  cout << m << " " << y << " " << costh << "  " << xmsq << endl;

  return xmsq;
}
