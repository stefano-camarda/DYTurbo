#include "phasespace.h"

#include <iostream>
#include <math.h>

//This is a boson variables binner: add mass here (not yet needed, because the mass is always in the phase space boundaries, will be needed when doing also mass bins)
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

  //move the costh cut to the lepton variables binner
  if (phasespace::cthmin != -1 && phasespace::cthmax != 1)
    {
      double costh = ((p3[3]+p3[2])*(p4[3]-p4[2])-(p4[3]+p4[2])*(p3[3]-p3[2]))/sqrt(m*m*(m*m+qt*qt));
      costh *= (y < 0. ? -1 : 1); //sign flip according to boson rapidity
      if (costh < phasespace::cthmin || costh > phasespace::cthmax)
	return false;
    }
  
  //cout << "qt " << qt << " y " << y << endl;
  return true;
}

//Make also a lepton variables binner, which is needed for producing distributions as a function of lepton variables with quadrature integration
