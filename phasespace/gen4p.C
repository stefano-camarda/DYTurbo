#include "phasespace.h"
#include "settings.h"

#include <iostream>
#include <math.h>

//rest frame axes
double phasespace::kap1[4];
double phasespace::xax[3];
double phasespace::yax[3];
double phasespace::zax[3];

phasespace::restframeid CurrentRF;
#pragma omp threadprivate(CurrentRF)

//the m, qt, y, and phi variables are taken from the phasespace:: namelist
void phasespace::genV4p()
{
  //Generate the boson 4-momentum

  //vector boson momentum: pV[3]^2-pV[0]^2-pV[1]^2-pV[2]^2=m2
  double cosphiV, sinphiV;
  if (phiV == 0.)
    {
      cosphiV = 1.;
      sinphiV = 0.;
    }
  else
    {
      cosphiV = cos(phiV);
      sinphiV = sqrt(max(0.,1.-pow(cosphiV,2)))*(phiV>0 ? 1 : -1);
    }
  pV[0]=qt*cosphiV;           //px
  pV[1]=qt*sinphiV;           //py
  pV[2]=0.5*mt*(exppy-expmy); //pz
  pV[3]=0.5*mt*(exppy+expmy); //E

  //Calculate the boost 4-vector from the boson rest frame to the laboratory frame (-pV)
  boostv_(m,pV,gam,beta);
}

void phasespace::genl4p()
{
  //simple phase space generation in the naive dilepton rest frame
  //with axes x = (1,0,0); y = (0,1,0); z=(0,0,1)
  //p3 is the lepton and p4 is the antilepton

  //Generate the 4-momentum of the first lepton in the boson rest frame, using costh and phi_lep
  double p3cm[4];
  if (qt == 0. || CurrentRF == naive)
    genp_(costh, phi_lep, m, p3cm);
  else
    {
      //Start from the z axis, and rotate by angle theta with respect to y axis
      double rot1[3];
      double c = phasespace::costh;
      double s = sqrt(max(0.,1.-pow(c,2)));
      phasespace::rotate(zax, c, s, yax, rot1);

      //rotate by angle phi_lep with respect to z axis
      double rot2[3];
      c = cos(phasespace::phi_lep-M_PI);      //Use the same M_PI rotation convention as in kinematic
      s = sqrt(max(0.,1.-pow(c,2)))*((phasespace::phi_lep-M_PI)>0 ? 1 : -1);
      phasespace::rotate(rot1, c, s, zax, rot2);

      p3cm[3]=phasespace::m/2.; //E
      p3cm[0]=p3cm[3]*rot2[0];  //px
      p3cm[1]=p3cm[3]*rot2[1];  //py
      p3cm[2]=p3cm[3]*rot2[2];  //pz
    }
  //Boost to go in the lab frame
  dyboost_(gam, beta, p3cm, p3);
  
  //momentum of the second lepton
  p4[0]=pV[0]-p3[0];      //px
  p4[1]=pV[1]-p3[1];      //py
  p4[2]=pV[2]-p3[2];      //pz
  //p4[3]=pV[3]-p3[3];      //E
  //ensure energy is calculated with enough precision
  p4[3]=sqrt(p4[0]*p4[0]+p4[1]*p4[1]+p4[2]*p4[2]);
}

void phasespace::genp12()
{
  //Generate the incoming partons
  p1[3]=-x1*opts.sroot*0.5;
  p1[0]=0.;
  p1[1]=0.;
  p1[2]=-x1*opts.sroot*0.5;

  p2[3]=-x2*opts.sroot*0.5;
  p2[0]=0.;
  p2[1]=0.;
  p2[2]=+x2*opts.sroot*0.5;
}

void phasespace::genp5()
{
  //Generate jet 4-momentum from momentum conservation
  p5[0]=-p1[0]-p2[0]-pV[0];
  p5[1]=-p1[1]-p2[1]-pV[1];
  p5[2]=-p1[2]-p2[2]-pV[2];
  //p5[3]=-p1[3]-p2[3]-pV[3];
  //ensure energy is calculates with enough precision
  p5[3]=sqrt(p5[0]*p5[0]+p5[1]*p5[1]+p5[2]*p5[2]);
}

void phasespace::genRFaxes(restframeid RF)
{
  CurrentRF = RF;
  //In some cases the generation of the first lepton 4-momentum is trivial,
  //i.e. it is done with respect to the reference system x = (1,0,0); y = (0,1,0); z=(0,0,1)
  //Do not need to calculate axes in these cases:
  //a) When qt=0, as in fixed order predictions, because all rest frames are equivalent
  //b) In the naive prescription, the rest frame axes are the native x,y,z axes of the laboratory frame
  if (qt == 0. || RF == naive)
    return;

  double m2 = pow(phasespace::m,2);
  double qt2 = pow(phasespace::qt,2);
  double mt2 = m2+qt2;

  //recoil prescriptions are defined by the values of kt1 and kt2
  double kt1, kt2;

  //CS frame prescription
  if (RF == CS)
    {
      kt1 = phasespace::pV[0]/2.;
      kt2 = phasespace::pV[1]/2.;
    }

  //naive prescription (also called MY prescription in DYRES)
  if (RF == naive)
    {
      kt1=(1.+phasespace::pV[2]/(sqrt(m2)+phasespace::pV[3]))*phasespace::pV[0]/2; //this prescription should be still symmetric, because for kt_1(y)=kt_2(-y) and kt_2(y)=kt_1(-y)  
      kt2=(1.+phasespace::pV[2]/(sqrt(m2)+phasespace::pV[3]))*phasespace::pV[1]/2;
    }  

  //alternative k1t = 0 prescription
  if (RF == kt0)
    {
      kt1 = 0; //this prescription is unphysical. It should be at least: kt1 = pV[0] * (phasespace::y < 0 ? 1 : 0); kt2 = pV[1] * (phasespace::y < 0 ? 1 : 0)
      kt2 = 0;
    }
  
  //zeta1 as in Eq.(26) of arXiv:1507.06937
  double ktdotpV = phasespace::pV[0]*kt1+phasespace::pV[1]*kt2;
  double zeta1 = 1./m2/2.*(m2+2.*(ktdotpV)+sqrt(pow((m2+2.*(ktdotpV)),2)-4.*mt2*(pow(kt1,2)+pow(kt2,2))));
  double qP1 = (phasespace::pV[3]-phasespace::pV[2])*opts.sroot/2.;
  
  //kap1 is the colliding parton a1 after the lorentz transformation from the boson rest frame to the laboratory frame
  kap1[3] = opts.sroot/2.*(zeta1*m2/2./qP1+(pow(kt1,2)+pow(kt2,2))/zeta1*qP1/m2/pow(opts.sroot,2)*2.);
  kap1[0] = kt1;
  kap1[1] = kt2;
  kap1[2] = opts.sroot/2.*(zeta1*m2/2./qP1-(pow(kt1,2)+pow(kt2,2))/zeta1*qP1/m2/pow(opts.sroot,2)*2.);

  //Determine the x,y,z axes of the boson rest frame for any general kt1 prescription
  //Calculate z axis by boosting back kap1 from the laboratory frame to the boson rest frame
  double bt[3];
  bt[0]=-phasespace::beta[0];
  bt[1]=-phasespace::beta[1];
  bt[2]=-phasespace::beta[2];
  
  double bdotk1=kap1[0]*bt[0]+kap1[1]*bt[1]+kap1[2]*bt[2];
  zax[0]=kap1[0]+phasespace::gam*bt[0]*(phasespace::gam/(phasespace::gam+1)*bdotk1-kap1[3]);
  zax[1]=kap1[1]+phasespace::gam*bt[1]*(phasespace::gam/(phasespace::gam+1)*bdotk1-kap1[3]);
  zax[2]=kap1[2]+phasespace::gam*bt[2]*(phasespace::gam/(phasespace::gam+1)*bdotk1-kap1[3]);

  //normalise z axis to unity
  double zaxabs = sqrt(pow(zax[0],2)+pow(zax[1],2)+pow(zax[2],2));
  zax[0] = zax[0] / zaxabs;
  zax[1] = zax[1] / zaxabs;
  zax[2] = zax[2] / zaxabs;

  //define y axis ortogonal to z and V direction
  yax[0] = zax[1]*phasespace::pV[2] - zax[2]*phasespace::pV[1];
  yax[1] = zax[2]*phasespace::pV[0] - zax[0]*phasespace::pV[2];
  yax[2] = zax[0]*phasespace::pV[1] - zax[1]*phasespace::pV[0];
  
  //normalise y axis to unity
  double yaxabs = sqrt(pow(yax[0],2)+pow(yax[1],2)+pow(yax[2],2));
  yax[0] = yax[0] / yaxabs;
  yax[1] = yax[1] / yaxabs;
  yax[2] = yax[2] / yaxabs;

  //define x axis ortogonal to y and z direction
  xax[0] = yax[1]*zax[2] - yax[2]*zax[1];
  xax[1] = yax[2]*zax[0] - yax[0]*zax[2];
  xax[2] = yax[0]*zax[1] - yax[1]*zax[0];
  
  //normalise x axis to unity
  double xaxabs = sqrt(pow(yax[0],2)+pow(yax[1],2)+pow(yax[2],2));
  xax[0] = xax[0] / xaxabs;
  xax[1] = xax[1] / xaxabs;
  xax[2] = xax[2] / xaxabs;

  /*
  cout << endl;
  cout << "rest frame axes x,y,z" << endl;
  cout << xax[0] << " " << xax[1] << "  " << xax[2] << endl;
  cout << yax[0] << " " << yax[1] << "  " << yax[2] << endl;
  cout << zax[0] << " " << zax[1] << "  " << zax[2] << endl;
  cout << "boson momentum " << endl;
  cout << phasespace::pV[0] << " " << phasespace::pV[1] << "  " << phasespace::pV[2] << endl;
  */
}
