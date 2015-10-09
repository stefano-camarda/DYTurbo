#include "switch.h"

#include <math.h>
#include <TMath.h>

double switching::k = 3./4.;
double switching::cutoff = 0.01;
int switching::mode = 1;

void switching::init(int md, double kk)
{
  switching::k = kk;
  switching::mode = md;
}

double switching::swtch(double qt, double m)
{
  double swtch=1.;

  if (mode == 1)
    if (qt >= m*k) swtch=exp(-pow((m*k-qt),2)/pow((m/2.),2)); // GAUSS SWITCH
  
  /*
  if (qt >= m)       swtch=exp(-pow((m-qt),2)/pow((m/2),2));          // GAUSS SWITCH
  if (qt >= m/2)     swtch=exp(-pow((m/2-qt),2)/pow((m/2),2));        // GAUSS SWITCH
  if (qt >= m/2)     swtch=exp(-pow(2*(m/2-qt),2)/pow((m/2),2));      // GAUSS SWITCH FASTER
  if (qt >= m/2)     swtch=exp(((m*m)/4-(qt*qt)) / (m*m) );           // EXP SWITCH
  if (qt >= m/2)     swtch=exp(((m*m)/4-(qt*qt))/(m*m)*4);            // EXP SWITCH FASTER
  if (qt >= m/2)     swtch=(cos(TMath::Pi()/50*(qt-45.6))+1)/2;// COS SWITCH
  if (qt >= 95.6)    swtch=0;   
  */

  if (swtch <= cutoff)
    swtch = 0;

  return swtch;
}

double switching::qtlimit(double m)
{
  double limit=10000;
  if (mode == 1)
    limit = m*(sqrt(log(1./cutoff))+k);
  return limit;
}

double switching_(double &qt, double &m)
{
  switching::swtch(qt, m);
}
