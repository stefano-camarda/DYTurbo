#include "vjint.h"
#include "settings.h"
#include "interface.h"
#include "resconst.h"
#include <iomanip>

#include "LHAPDF/LHAPDF.h"

#include <math.h>

double vjint::brz;
double vjint::brw;

void vjint::init()
{
  //gevpb_.gevpb_ = 3.8937966e8; //MCFM 6.8 value
  gevpb_.gevpb_ = 0.389379e9; //dyres value

  flagch_.flagch_ = 0;

  vegint_.iter_ = 20;
  vegint_.locall_ = 10000;
  vegint_.nlocall_ = 20000;

  para_.ss_ =  opts.sroot;
  para_.s_ = pow(para_.ss_,2);

  pdf_.ih1_ = opts.ih1;
  pdf_.ih2_ = opts.ih2;

  vjorder_.iord_ = opts.order - 1;

  if (opts.nproc == 3)
    {
      if (opts.useGamma && !opts.zerowidth)
	prodflag_.prodflag_ = 5; 
      else
	prodflag_.prodflag_ = 3;
    }
  else if (opts.nproc == 1)
    prodflag_.prodflag_ = 21;
  else if (opts.nproc == 2)
    prodflag_.prodflag_ = 22;

  ckm_.vud_ = cabib_.Vud_;
  ckm_.vus_ = cabib_.Vus_;
  ckm_.vub_ = cabib_.Vub_;
  ckm_.vcd_ = cabib_.Vcd_;
  ckm_.vcs_ = cabib_.Vcs_;
  ckm_.vcb_ = cabib_.Vcb_;
  ckm_.vtd_ = 0.;
  ckm_.vts_ = 0.;
  ckm_.vtb_ = 0.;

  double cw2 = 1.- ewcouple_.xw_;
  em_.aemmz_ = sqrt(2.)* ewcouple_.Gf_ *pow(dymasses_.wmass_,2)* ewcouple_.xw_ /M_PI;

  brz = 1./dymasses_.zwidth_*em_.aemmz_/24.*dymasses_.zmass_ *(pow(-1.+2.*ewcouple_.xw_,2)+pow(2.*ewcouple_.xw_,2))/(ewcouple_.xw_*cw2); //=0.033638
  brw = 1./dymasses_.wwidth_*em_.aemmz_*dymasses_.wmass_/(12.*ewcouple_.xw_); //=0.10906 

  if (opts.zerowidth)
    {
      sigs_.sigz_ = brz; // /(16.*cw2)
      sigs_.sigw_ = brw; // /4.
    }
  
  //quarks are ordered according to mass:
  //1,2,3,4,5,6
  //u,d,s,c,b,t
  
  //double equ = 2./3.;  //up-quarks electric charge
  //double eqd = -1./3.; //down-quarks electric charge
  //quarks_.eq_[0]=equ;                 
  //quarks_.eq_[1]=eqd;
  //quarks_.eq_[2]=eqd;
  //quarks_.eq_[3]=equ;
  //quarks_.eq_[4]=eqd;
  quarks_.eq_[0] = ewcharge_.Q_[MAXNF+2];
  quarks_.eq_[1] = ewcharge_.Q_[MAXNF+1];
  quarks_.eq_[2] = ewcharge_.Q_[MAXNF+3];
  quarks_.eq_[3] = ewcharge_.Q_[MAXNF+4];
  quarks_.eq_[4] = ewcharge_.Q_[MAXNF+5];

  //definition of 'generalized' ckm matrix:
  //    (uu ud us uc ub ut)
  //    (du dd ds dc db dt)
  //ckm=(su sd ss sc sb st)
  //    (cu cd cs cc cb ct)
  //    (bu bd bs bc bb bt)
  //    (tu td ts tc tb tt)
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
      quarks_.ckm_[j][i] = 0.;
  
  quarks_.ckm_[1][0]=ckm_.vud_;
  quarks_.ckm_[2][0]=ckm_.vus_;
  quarks_.ckm_[4][0]=ckm_.vub_;
  quarks_.ckm_[0][1]=ckm_.vud_;
  quarks_.ckm_[3][1]=ckm_.vcd_;
  quarks_.ckm_[5][1]=ckm_.vtd_;
  quarks_.ckm_[0][2]=ckm_.vus_;
  quarks_.ckm_[3][2]=ckm_.vcs_;
  quarks_.ckm_[5][2]=ckm_.vts_;
  quarks_.ckm_[1][3]=ckm_.vcd_;
  quarks_.ckm_[2][3]=ckm_.vcs_;
  quarks_.ckm_[4][3]=ckm_.vcb_;
  quarks_.ckm_[0][4]=ckm_.vub_;
  quarks_.ckm_[3][4]=ckm_.vcb_;
  quarks_.ckm_[5][4]=ckm_.vtb_;
  quarks_.ckm_[1][5]=ckm_.vtd_;
  quarks_.ckm_[2][5]=ckm_.vts_;
  quarks_.ckm_[4][5]=ckm_.vtb_;
      
  //definition of 'delta' matrix
  for (int i = 0; i < MAXNF; i++)
    for (int j = 0; j < MAXNF; j++)
      quarks_.delta_[j][i] = 0.;

  for (int i = 0; i < MAXNF; i++)
    quarks_.delta_[i][i] = 1.;

  //definition of tau3's Pauli matrix
  for (int i = 0; i < MAXNF; i++)
    for (int j = 0; j < MAXNF; j++)
      quarks_.tau3_[j][i] = 0.;

  //quarks_.tau3_[0][0]=1.;
  //quarks_.tau3_[1][1]=-1.;
  //quarks_.tau3_[2][2]=-1.;
  //quarks_.tau3_[3][3]=1.;
  //quarks_.tau3_[4][4]=-1.;
  quarks_.tau3_[0][0] = ewcharge_.tau_[MAXNF+2];
  quarks_.tau3_[1][1] = ewcharge_.tau_[MAXNF+1];
  quarks_.tau3_[2][2] = ewcharge_.tau_[MAXNF+3];
  quarks_.tau3_[3][3] = ewcharge_.tau_[MAXNF+4];
  quarks_.tau3_[4][4] = ewcharge_.tau_[MAXNF+5];

  //definition of constants and couplings
  const2_.ca_ = 3.;
  const2_.xnc_ = 3.;
  const2_.cf_ = 4./3.;
  const2_.tr_ = resconst::NF/2.;
  const2_.pi_ = M_PI;
  couplings_.alpha0_ = em_.aemmz_;
  couplings_.xw_ = ewcouple_.xw_;
  couplings_.sw_ = sqrt(couplings_.xw_); // sin_w
  couplings_.cw_ = sqrt(1.-couplings_.xw_); // cos_w

//c.....read the input file
//      if(ic.eq.1) then
//         ih1=1                      !beam 1 @ LHC!
//         ih2=1                      !beam 2 @ LHC!
//      elseif(ic.eq.-1) then
//         ih1=1                      !beam 1 @ Tevatron!
//         ih2=-1                     !beam 2 @ Tevatron!
//      endif
}


double vjint::vint(double m, double pt, double y)
{
  //set scales and alpha strong
  if (opts.dynamicscale)
    {
      scales2_.xmur_ = m*opts.kmuren;
      scales2_.xmuf_ = m*opts.kmufac;
    }
  else
    {
      scales2_.xmur_ = opts.rmass*opts.kmuren;
      scales2_.xmuf_ = opts.rmass*opts.kmufac;
    }
  
  scales2_.xmur2_ = pow(scales2_.xmur_,2);
  scales2_.xmuf2_ = pow(scales2_.xmuf_,2);
  asnew_.as_ = LHAPDF::alphasPDF(scales2_.xmur_)/M_PI;

  //calculate propagators
  double q2 = m*m;
  double cw2 = 1.- ewcouple_.xw_;
  sigs_.sigz_ = brz*dymasses_.zwidth_*q2/(M_PI*dymasses_.zmass_) / (pow(q2-pow(dymasses_.zmass_,2),2)+pow(dymasses_.zmass_,2)*pow(dymasses_.zwidth_,2)); // /(16.*cw2)
  sigs_.sigw_ = brw*dymasses_.wwidth_*q2/(M_PI*dymasses_.wmass_) / (pow(q2-pow(dymasses_.wmass_,2),2)+pow(dymasses_.wmass_,2)*pow(dymasses_.wwidth_,2)); // !/4.
  sigs_.siggamma_ = em_.aemmz_/(3.*M_PI*q2);
  sigs_.sigint_ = -em_.aemmz_/(6.*M_PI)*(q2-pow(dymasses_.zmass_,2)) /(pow(q2-pow(dymasses_.zmass_,2),2)+pow(dymasses_.zmass_,2)*pow(dymasses_.zwidth_,2)) * (-1.+4.*ewcouple_.xw_)/(2.*sqrt(ewcouple_.xw_*cw2));  // !sqrt(sw2/cw2)/2.

  //set phase space variables
  internal_.q_ = m;
  internal_.q2_ = pow(m,2);
  internal_.qt_ = pt;
  yv_.yv_ = y;
  yv_.expyp_ = exp(y);
  yv_.expym_ = exp(-y);

  if (opts.zerowidth)
    internal_.q_ = opts.rmass;
  
  //call cross section calculation
  int ord = opts.order - 1;
  double res, err, chi2a;
  qtdy_(res,err,chi2a,y,y,ord);
  //cout << setprecision(16) << "C++ " << m << "  " << pt << "  " << "  " << y << "  " << res << endl;

  //Apply conversion factors:
  //ds/dqt = ds/dqt2 * dqt2/dqt
  //fb = 1000 * pb
  return 2*pt*res*1000.;
}
