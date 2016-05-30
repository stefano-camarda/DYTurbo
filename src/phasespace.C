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

//global variables for phase space
double phasespace::m;
double phasespace::qt;
double phasespace::y;
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

void phasespace::set_mqtycth(double M, double Qt, double Y, double Costh)
{
  m = M;
  qt = Qt;
  y = Y;
  costh = Costh;
}
void phasespace::set_mqty(double M, double Qt, double Y)
{
  m = M;
  qt = Qt;
  y = Y;
}
void phasespace::set_m(double M) {m = M;}
void phasespace::set_qt(double Qt) {qt = Qt;}
void phasespace::set_y(double Y) {y = Y;}
void phasespace::set_cth(double Costh) {costh = Costh;}

//fortran functions
void setqt_(double &qtt) {phasespace::qt = qtt;}
void sety_(double &yy) {phasespace::y = yy;}

//this is a boson variables binner: add mass here (not really needed, beacause the mass is always in the phase space boundaries)
//and make a lepton variables binner
int binner_(double p3[4], double p4[4])
{
  double qt = sqrt((float)pow(p3[0]+p4[0],2)+pow(p3[1]+p4[1],2));
  if (qt < phasespace::qtmin || qt > phasespace::qtmax)
    return false;

  double y = 0.5 *log((p3[3] + p4[3] + p3[2] + p4[2]) / (p3[3] + p4[3] - p3[2] - p4[2]));
  if (y < phasespace::ymin || y > phasespace::ymax)
    return false;

  //cout << "qt " << qt << " y " << y << endl;
  return true;
}

