#include "phasespace.h"

#include <iostream>
#include <math.h>

//integration boundaries
double phasespace::mmin;
double phasespace::mmax;
double phasespace::qtmin;
double phasespace::qtmax;
double phasespace::ymin;
double phasespace::ymax;

double phasespace::cthmin;
double phasespace::cthmax;

//global variables for the generation of the phase space
double phasespace::m;
double phasespace::qt;
double phasespace::y;
double phasespace::phiV;

double phasespace::costh;

using namespace std;

void phasespace::setbounds(double m1, double m2, double qt1, double qt2, double y1, double y2)
{
  mmin = m1;
  mmax = m2;
  qtmin = qt1;
  qtmax = qt2;
  ymin = y1;
  ymax = y2;
    
  //Check ordering
  if (mmin > mmax || qtmin > qtmax || ymin > ymax)
    {
      cout << "Error on integration boundaries" << endl;
      exit(-1);
    }
}

void phasespace::setcthbounds(double cth1, double cth2)
{
  cthmin = cth1;
  cthmax = cth2;
    
  //Check ordering
  if (cthmin > cthmax)
    {
      cout << "Error on costh integration boundaries" << endl;
      exit(-1);
    }
}

void phasespace::set_mqtyphi(double M, double Qt, double Y, double PhiV)
{
  m = M;
  qt = Qt;
  y = Y;
  phiV = PhiV;
}
void phasespace::set_m(double M) {m = M;}
void phasespace::set_qt(double Qt) {qt = Qt;}
void phasespace::set_y(double Y) {y = Y;}
void phasespace::set_phiV(double PhiV) {phiV = PhiV;}

void phasespace::set_cth(double Costh) {costh = Costh;}

//fortran functions
void setqt_(double &qtt) {phasespace::qt = qtt;}
void sety_(double &yy) {phasespace::y = yy;}

//This is a boson variables binner: add mass here (not yet needed, beacause the mass is always in the phase space boundaries, will be needed when doing also mass bins)
//the boson binner function is used only with the MCFM integrands, it is needed because the generation of the phase space is not done with respect to m, pt, y
int binner_(double p3[4], double p4[4])
{
  double qt = sqrt((float)pow(p3[0]+p4[0],2)+pow(p3[1]+p4[1],2));
  if (qt < phasespace::qtmin || qt > phasespace::qtmax)
    return false;

  double y = 0.5 *log((p3[3] + p4[3] + p3[2] + p4[2]) / (p3[3] + p4[3] - p3[2] - p4[2]));
  if (y < phasespace::ymin || y > phasespace::ymax)
    return false;

  double m = sqrt(pow(p3[3]+p4[3],2) - pow(p3[2]+p4[2],2) - qt*qt);
  if (m < phasespace::mmin || m > phasespace::mmax)
    return false;
  
  //cout << "qt " << qt << " y " << y << endl;
  return true;
}

//Make also a lepton variables binner, which is needed for producing distributions as a function of lepton variables with quadrature integration
