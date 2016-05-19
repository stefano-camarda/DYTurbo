#include "cuts.h"
#include "settings.h"
#include "interface.h"

#include <iostream>
#include <math.h>

double cuts::getEtMiss(double p3[4], double p4[4]){
  double EtMiss = 0.;
  if (opts.nproc == 1)
    EtMiss = getPt(p3);
  else if (opts.nproc == 2)
    EtMiss = getPt(p4);
  return EtMiss;
}

double cuts::getLPt(double p3[4], double p4[4]){
  double LPt = 0.;
  if (opts.nproc == 1)
    LPt = getPt(p4);
  else if (opts.nproc == 2)
    LPt = getPt(p3);
  return LPt;
}

double cuts::getLY(double p3[4], double p4[4]){
  double LY = 0.;
  if (opts.nproc == 1)
    LY = getY(p4);
  else if (opts.nproc == 2)
    LY = getY(p3);
  return LY;
}

double cuts::getPt(double p[4]){
    return sqrt(pow(p[0],2)+pow(p[1],2));
}
double cuts::getY(double p[4]){
    return 0.5 *log((p[3] + p[2]) / (p[3] - p[2]));
}
double cuts::getEta(double p[4]){
    return atanh(p[2]/sqrt((float) p[0]*p[0] + p[1]*p[1] + p[2]*p[2]));
}
double cuts::getMt(double p3[4], double p4[4])
{
  //double p[4];
  //p[2]=p3[2]+p4[2];
  //p[3]=p3[3]+p4[3];
  //return sqrt((float)pow(p[3],2)-pow(p[2],2));

  double cosphi = (p3[0]*p4[0]+p3[1]*p4[1])/sqrt((pow(p3[0],2)+pow(p3[1],2))*(pow(p4[0],2)+pow(p4[1],2)));
  double mtrans = 2.*sqrt(pow(p3[0],2)+pow(p3[1],2))*sqrt(pow(p4[0],2)+pow(p4[1],2))*(1.-cosphi);
  mtrans = sqrt(max(mtrans,0.));

  /*
  TLorentzVector el, nu;
  if (opts.nproc == 2)
    {
      el.SetPxPyPzE(p3[0],p3[1],p3[2],p3[3]);
      nu.SetPxPyPzE(p4[0],p4[1],p4[2],p4[3]);
    }
  else if (opts.nproc == 1)
    {
      el.SetPxPyPzE(p4[0],p4[1],p4[2],p4[3]);
      nu.SetPxPyPzE(p3[0],p3[1],p3[2],p3[3]);
    }
  double Mt = sqrt(2 * el.Pt() * nu.Pt() * (1 - TMath::Cos(el.DeltaPhi(nu))));
  */

  return mtrans;
}
double cuts::getM(double p3[4], double p4[4]){
    double p[4];
    p[0] = p3[0] + p4[0];
    p[1] = p3[1] + p4[1];
    p[2] = p3[2] + p4[2];
    p[3] = p3[3] + p4[3];
    return sqrt(float(p[3]*p[3] - p[0]*p[0] - p[1]*p[1] - p[2]*p[2]));
}

//fortran interface
int cuts_(double p[4][12], int &njet){
    double p3[4];
    double p4[4];
    for (int i=0; i<4; i++){
        p3[i] = p[i][2];
        p4[i] = p[i][3];
    }
    //MCFM expects opposite logic false=accept event
    return !cuts::lep(p3,p4);
}

//C native function
bool cuts::lep(double p3[4], double p4[4])
{
  if (!opts.makelepcuts) return true;
  return decide_fiducial(p3,p4);
}

bool cuts::decide_fiducial( double p3[4], double p4[4] ){
  switch (opts.fiducial)
    {
    case cuts::D0    : return fiducial_D0    (p3,p4); break;
    case cuts::CDF   : return fiducial_CDF   (p3,p4); break;
    case cuts::ATLAS : return fiducial_ATLAS (p3,p4); break;
    case cuts::CMS7  : return fiducial_CMS7  (p3,p4); break;
    case cuts::CMS8  : return fiducial_CMS8  (p3,p4); break;
    case cuts::GENEXP :
      {
	/*****************************/
	// Can we clean up this or is it still needed?

	// Fiducial type is zero so in normal run there is no cut. But,
	// because for wwidth I need fiducial / full I decided to run one job
	// with full integral (makelepcuts=false) and plots fill only fiducial
	// (fiducial!=0). In case I want to run standard full-integral
	// full-plotting I need to sef fiducial=0 and makecuts=false
	if (!opts.makelepcuts) return true; //< so this is here for plotting.
	/*****************************/

	if (opts.lptcut > 0)
	  if (opts.nproc == 3)
	    {
	      float pt3 = sqrt((float)pow(p3[0],2)+(float)pow(p3[1],2));
	      if (pt3 < opts.lptcut)
		return false;
	      float pt4 = sqrt((float)pow(p4[0],2)+(float)pow(p4[1],2));
	      if (pt4 < opts.lptcut)
		return false;
	    }
	  else
	    if (getLPt(p3, p4) < opts.lptcut)
	      return false;
	
	if (opts.lycut < 100)
	  if (opts.nproc == 3)
	    {
	      float y3 = 0.5 *log(((float)p3[3] + (float)p3[2]) / ((float)p3[3] - (float)p3[2]));
	      if (fabs(y3) > opts.lycut)
		return false;
	      float y4 = 0.5 *log(((float)p4[3] + (float)p4[2]) / ((float)p4[3] - (float)p4[2]));
	      if (fabs(y4) > opts.lycut)
		return false;
	    }
	  else
	    if (fabs(getLY(p3, p4)) > opts.lycut)
	      return false;

	if (opts.mtcut > 0)
	    if (getMt(p3, p4) < opts.mtcut)
	      return false;

	if (opts.etmisscut > 0)
	    if (getEtMiss(p3, p4) < opts.etmisscut)
	      return false;
	
	if (opts.l1ptcut > 0 || opts.l2ptcut > 0
	    || opts.l1ycut < 100 || opts.l2ycut < 100)
	  {
	    float pt3 = sqrt((float)pow(p3[0],2)+(float)pow(p3[1],2));
	    float pt4 = sqrt((float)pow(p4[0],2)+(float)pow(p4[1],2));
	    float y3 = 0.5 *log(((float)p3[3] + (float)p3[2]) / ((float)p3[3] - (float)p3[2]));
	    float y4 = 0.5 *log(((float)p4[3] + (float)p4[2]) / ((float)p4[3] - (float)p4[2]));

	    float pt1, pt2, y1, y2;
	    if (pt3 > pt4)
	      {
		pt1 = pt3; y1 = y3;
		pt2 = pt4; y2 = y4;
	      }
	    else
	      {
		pt1 = pt4; y1 = y4;
		pt2 = pt3; y2 = y3;
	      }

	    if (pt1 < opts.l1ptcut)
	      return false;
	    if (pt2 < opts.l2ptcut)
	      return false;
	    if (fabs(y1) > opts.l1ycut)
	      return false;
	    if (fabs(y2) > opts.l2ycut)
	      return false;
	  }

	
	break;
      }
    case cuts::CUSTOM : return user_cuts(p3,p4); break;
    default:
      cout << "not recognised cuts" << endl;
    }
  return true;
}

/// fiducial cuts
bool cuts::fiducial_D0(double p3[4], double p4[4]){
    // electrons only
    double pt3 = getPt(p3);
    if (pt3<25) return false;
    double aeta3 = fabs(getEta(p3));
    if (aeta3>2.5 || (aeta3>1.1 && aeta3<1.5) ) return false;
    double pt4 = getPt(p4);
    if (opts.nproc==3){ // z
        if (pt4<25) return false;
        double aeta4 = fabs(getEta(p4));
        double m34 = getM(p3,p4);
        if (aeta4>2.5 || (aeta4>1.1 && aeta4<1.5) ) return false;
        if (m34 < 75) return false;
        if (m34 > 105) return false;
    } else { // W
        if (pt4<25) return false;
    }
    return true; // keep
}

bool cuts::fiducial_CDF(double p3[4], double p4[4]){
    //  muons
    double PTCUT=20; // 20:muon, 25:electrons
    double pt3 = getPt(p3);
    if (pt3<PTCUT) return false;
    double aeta3 = fabs(getEta(p3));
    if (aeta3>1.0 ) return false;
    //if (aeta3>2.0 || (aeta3>1.0 && aeta3<1.2) ) return false;
    //if (aeta3<1.0 && pt3<25) return false;
    double pt4 = getPt(p4);
    if (opts.nproc==3){ // z
        if (pt4<PTCUT) return false;
        double aeta4 = fabs(getEta(p4));
        if (aeta4>1.0 ) return false;
        double m34 = getM(p3,p4);
        if (m34 <  66.) return false;
        if (m34 > 116.) return false;
    } else { // W
        if (pt4<PTCUT) return false;
    }
    return true; // keep
    // electrons have pt3 25 pt4 25
}

bool cuts::fiducial_ATLAS(double p3[4], double p4[4]){
    double pt3 = getPt(p3);
    if (pt3<20) return false;
    double aeta3 = fabs(getEta(p3));
    if (aeta3>2.4 ) return false;
    /* electron
    if (aeta3>2.47 ) return false;
    if (aeta3>1.37 && aeta3 < 1.52 ) return false;
     */
    double pt4 = getPt(p4);
    if (opts.nproc==3){ // z
        if (pt4<20) return false;
        double aeta4 = fabs(getEta(p4));
        if (aeta4>2.4 ) return false;
        double m34 = getM(p3,p4);
        if (m34 <  66.) return false;
        if (m34 > 116.) return false;
    } else { // W
        if (pt4<25) return false;
        double mt = getMt(p3,p4);
        if (mt<40) return false;
    }
    return true; // keep
}

bool cuts::fiducial_CMS7(double p3[4], double p4[4]){
    // electron pt 25
    // electron eta 0 < aeta < 1.44
    // electron eta 1.57 < aeta < 2.5
    double pt3 = getPt(p3);
    if (pt3<25) return false;
    double aeta3 = fabs(getEta(p3));
    if (aeta3>2.1 ) return false;
    double pt4 = getPt(p4);
    if (opts.nproc==3){ // z
        if (pt4<25) return false;
        double aeta4 = fabs(getEta(p4));
        if (aeta4>2.1 ) return false;
        double m34 = getM(p3,p4);
        if (m34 <  60.) return false;
        if (m34 > 120.) return false;
    } else { // W
        //if (pt4<25) return false;
        //double mt = getMt(p3,p4);
        //if (mt<40) return false;
        return true;
    }
    return true; // keep
}

bool cuts::fiducial_CMS8(double p3[4], double p4[4]){
    // electron pt 25
    // electron eta 0 < aeta < 1.44
    // electron eta 1.57 < aeta < 2.5
    double pt3 = getPt(p3);
    if (pt3<25) return false;
    double aeta3 = fabs(getEta(p3));
    if (aeta3>2.1 ) return false;
    double pt4 = getPt(p4);
    if (opts.nproc==3){ // z
        if (pt4<25) return false;
        double aeta4 = fabs(getEta(p4));
        if (aeta4>2.1 ) return false;
        double m34 = getM(p3,p4);
        if (m34 <  60.) return false;
        if (m34 > 120.) return false;
    } else { // W
        //if (pt4<25) return false;
        //double mt = getMt(p3,p4);
        //if (mt<40) return false;
        return true;
    }
    return true; // keep
}

