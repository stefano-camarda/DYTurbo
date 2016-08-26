#include "phasespace.h"

#include <iostream>

//integration boundaries
double phasespace::mmin;
double phasespace::mmax;
double phasespace::qtmin;
double phasespace::qtmax;
double phasespace::ymin;
double phasespace::ymax;

double phasespace::cthmin = -1;
double phasespace::cthmax = +1;

//global variables for the generation of the phase space
double phasespace::m;
double phasespace::qt;
double phasespace::y;
double phasespace::phiV;

//decay angles
double phasespace::costh;
double phasespace::phi_lep;

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
void phasespace::set_philep(double Phi_lep) {phi_lep = Phi_lep;}

//fortran functions
void setqt_(double &qtt) {phasespace::qt = qtt;}
void sety_(double &yy) {phasespace::y = yy;}
