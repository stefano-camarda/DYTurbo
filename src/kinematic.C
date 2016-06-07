#include "kinematic.h"

#include <iostream>
#include <algorithm>

using namespace std;

double kinematic::lp[4], kinematic::lm[4];
double kinematic::v[4];
double kinematic::m, kinematic::m2, kinematic::qt, kinematic::qt2, kinematic::y;
double kinematic::phiZ;
double kinematic::costh, kinematic::phi_lep;
double kinematic::mtrans;

void kinematic::set(double p3[4], double p4[4])
{
  //lepton
  lm[0] = p3[0];
  lm[1] = p3[1];
  lm[2] = p3[2];
  lm[3] = p3[3];

  //anti-lepton
  lp[0] = p4[0];
  lp[1] = p4[1];
  lp[2] = p4[2];
  lp[3] = p4[3];

  //boson
  v[0] = lm[0]+lp[0];
  v[1] = lm[1]+lp[1];
  v[2] = lm[2]+lp[2];
  v[3] = lm[3]+lp[3];
}

void kinematic::calc_vb()
{
  //calculate boson variables
  calcM();
  calcQt();
  calcY();
  calcPhiZ();
  calcMt();
}

void kinematic::calc_angles()
{
  //calculate letpon angles in CS rest frame
  calcCosThCS();
  calcPhiCS();
}

void kinematic::calcM()
{
  m2 = v[3]*v[3] -v[0]*v[0] -v[1]*v[1] -v[2]*v[2];
  m = sqrt(m2);
}
void kinematic::calcQt()
{
  qt2 = v[0]*v[0]+v[1]*v[1];
  qt = sqrt(qt2);
}
void kinematic::calcY()
{
  y = 0.5*log((v[3]+v[2])/(v[3]-v[2]));
}
void kinematic::calcPhiZ()
{
  phiZ = atan2(v[1],v[0]);
}
void kinematic::calcMt()
{
  //mt = sqrt(2*ptl*ptnu*(1-cos(phi))
  mtrans = 2*(sqrt(lm[0]*lm[0]+lm[1]*lm[1])*sqrt(lp[0]*lp[0]+lp[1]*lp[1])-(lm[0]*lp[0]+lm[1]*lp[1]));
  mtrans = sqrt(max(mtrans,0.));
}
double kinematic::Vplus(double p[4])
{
  return (p[3]+p[2]);
}
double kinematic::Vminus(double p[4])
{
  return (p[3]-p[2]);
}

//Original formulas of costh and phi in CS frame are in Eq.(2.1) of Phys.Rev. D16 (1977) 2219
void kinematic::calcCosThCS()
{
  costh = (Vplus(lm)*Vminus(lp) - Vplus(lp)*Vminus(lm));
  costh /= sqrt(m2*(m2+qt2));
  costh *=  v[2] < 0. ? -1 : 1; //sign flip according to boson rapidity
}

void kinematic::calcPhiCS()
{
  double plxCS, plyCS;
  if (qt == 0) //if qt = 0, use the original x and y axis to determine phiCS in the boson rest frame
    {
      plxCS = lm[0];
      plyCS = lm[1];
    }
  else // if qt > 0, use the boson direction in the transverse plane as x axis in the boson rest frame
    {
      /*******************************************************************/
      //Mirkes definition as in Nucl.Phys. B387 (1992) 385 Eq.(22)
      //first needs to rotate the lab frame so that the x axis lies in the event plane defined by the boson and proton directions
      double c = v[0]/qt;//cos(-phiZ);
      double s = -v[1]/qt;//sin(-phiZ);
      double plx = c*lm[0] - s*lm[1];
      double ply = s*lm[0] + c*lm[1];

      //Now apply formulas (22) of Nucl.Phys. B387 (1992) 385
      plxCS = 0.5 * m / sqrt(m2+qt2) * (2.*plx - qt);
      plyCS = ply;
  
      /*******************************************************************/
      //Original formula, as in Eq.(2.1) of Phys.Rev. D16 (1977) 2219 
      double delta[2];
      delta[0] = lm[0]-lp[0];
      delta[1] = lm[1]-lp[1];

      //unit vector in the transverse plane from the projection of the boson momentum
      double qht[2];
      qht[0] = v[0]/qt;
      qht[1] = v[1]/qt;

      //unit vector in the transverse plane perpendicular to the boson momentum and beam 1
      double rht[2];
      rht[0] = -v[1]/qt;
      rht[1] = v[0]/qt;
  
      plxCS = 0.5 * m/sqrt(m2+qt2) * (delta[0]*qht[0]+delta[1]*qht[1]);
      plyCS = 0.5 * (delta[0]*rht[0]+delta[1]*rht[1]);
      /*******************************************************************/
    }
  double sign = v[2] > 0. ? 1 : -1; //sign flip according to boson rapidity
  phi_lep = atan2(-sign*plyCS,-plxCS); //rotate by M_PI (just a convention, it means that the x axis is opposite to the boson direction in the transverse plane)
}

//this function is not yet used, needs to be checked
void kinematic::rotate(double vin[3], double phi, double ax[3], double vout[])
{
  double c = cos(phi);
  double s = sin(phi);
  vout[0]=(c+ax[0]*ax[0]*(1-c))      *vin[0] + (ax[0]*ax[1]*(1-c)-ax[2]*s)*vin[1] + (ax[0]*ax[2]*(1-c)+ax[1]*s)*vin[2];
  vout[1]=(ax[1]*ax[0]*(1-c)+ax[2]*s)*vin[0] + (c+ax[1]*ax[1]*(1-c))      *vin[1] + (ax[1]*ax[2]*(1-c)-ax[0]*s)*vin[2];
  vout[2]=(ax[2]*ax[0]*(1-c)-ax[1]*s)*vin[0] + (ax[2]*ax[1]*(1-c)+ax[0]*s)*vin[1] + (c+ax[2]*ax[2]*(1-c))      *vin[2];
}
