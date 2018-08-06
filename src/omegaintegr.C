// --------- TO DO -----------------
//make this collection of functions a self-consistent library
//by removing dependences on settings.h and interface.h,
//introducing a namespace, and adding a init function for initialisation of settings
// ---------------------------------

#include "omegaintegr.h"
#include "settings.h"
#include "interface.h"
#include "phasespace.h"
#include "gaussrules.h"
#include "KinematicCuts.h"
#include "cuba.h"

#include <ctime>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <vector>


//rest frame axes
double omegaintegr::kap1[4];
double omegaintegr::xax[3];
double omegaintegr::yax[3];
double omegaintegr::zax[3];

//double omegaintegr::cthmin;
//double omegaintegr::cthmax;

//void omegaintegr::setcosthbounds(double min, double max)
//{
//  cthmin = min;
//  cthmax = max;
//}

//fortran interface
void genv4p_()
{
  omegaintegr::genV4p();
}

void cthmoments_(double &cthmom0, double &cthmom1, double &cthmom2)
{
  omegaintegr::cthmoments(cthmom0, cthmom1, cthmom2);
}

void genrfaxes_()
{
  omegaintegr::genRFaxes();
}


//the m, qt, y, and phi variables are taken from the phasespace:: namelist
void omegaintegr::genV4p()
{
  //Generate the boson 4-momentum
  phasespace::calcmt();
  phasespace::calcexpy();
  phasespace::genV4p();

  //Generate axes of the rest frame
  genRFaxes();
  return;
}


void omegaintegr::genRFaxes()
{
  phasespace::restframeid RF;
  if (opts.qtrec_naive) RF = phasespace::naive;
  else if (opts.qtrec_cs) RF = phasespace::CS;
  else if (opts.qtrec_kt0) RF = phasespace::kt0;

  phasespace::genRFaxes(RF);
  return;
  
  //In some cases the generation of the first lepton 4-momentum is trivial,
  //i.e. it is done with respect to the reference system x = (1,0,0); y = (0,1,0); z=(0,0,1)
  //Do not need to calculate axes in these cases:
  //a) When qt=0, as in fixed order predictions, because all rest frames are equivalent
  //b) In the naive prescription, the rest frame axes are the native x,y,z axes of the laboratory frame
  //c) With the cuba integration option, the phase space is generated uniformly
  //and costh is calculated using Catani's formulas arXiv:1507.0693 (25-32)
  if (phasespace::qt == 0. || opts.qtrec_naive || opts.cubaint)
    return;

  double m2 = pow(phasespace::m,2);
  double qt2 = pow(phasespace::qt,2);
  double mt2 = m2+qt2;

  //recoil prescriptions are defined by the values of kt1 and kt2
  double kt1, kt2;

  //CS frame prescription
  if (opts.qtrec_cs)
    {
      kt1 = phasespace::pV[0]/2.;
      kt2 = phasespace::pV[1]/2.;
    }

  //naive prescription (also called MY prescription in DYRES)
  if (opts.qtrec_naive)
    {
      kt1=(1.+phasespace::pV[2]/(sqrt(m2)+phasespace::pV[3]))*phasespace::pV[0]/2; //this prescription should be still symmetric, because for kt_1(y)=kt_2(-y) and kt_2(y)=kt_1(-y)  
      kt2=(1.+phasespace::pV[2]/(sqrt(m2)+phasespace::pV[3]))*phasespace::pV[1]/2;
    }  

  //alternative k1t = 0 prescription
  if (opts.qtrec_kt0)
    {
      kt1 = 0; //this prescription is unphysical. It should be at least: kt1 = pV[0] * (phasespace::y < 0 ? 1 : 0); kt2 = pV[1] * (phasespace::y < 0 ? 1 : 0)
      kt2 = 0;
    }
  
  //zeta1 as in Eq.(26) of arXiv:1507.06937
  double ktdotpV = phasespace::pV[0]*kt1+phasespace::pV[1]*kt2;
  double zeta1 = 1./m2/2.*(m2+2.*(ktdotpV)+sqrt(pow((m2+2.*(ktdotpV)),2)-4.*mt2*(pow(kt1,2)+pow(kt2,2))));
  double qP1 = (phasespace::pV[3]-phasespace::pV[2])*opts.sroot/2.;
  //double qP2 = (phasespace::pV[3]+phasespace::pV[2])*opts.sroot/2.; //qP2 is not used
  
  //kap1 is the colliding parton a1 after the lorentz transformation from the boson rest frame to the laboratory frame
  kap1[3] = opts.sroot/2.*(zeta1*m2/2./qP1+(pow(kt1,2)+pow(kt2,2))/zeta1*qP1/m2/pow(opts.sroot,2)*2.);
  kap1[0] = kt1;
  kap1[1] = kt2;
  kap1[2] = opts.sroot/2.*(zeta1*m2/2./qP1-(pow(kt1,2)+pow(kt2,2))/zeta1*qP1/m2/pow(opts.sroot,2)*2.);
  //  kap1[2] *= phasespace::y < 0 ? -1 : 1;

  //Determine the x,y,z axes of the boson rest frame for any general kt1 prescription
  //Calculate z axis by boosting back kap1 from the laboratory frame to the boson rest frame
  double bt[3];
  bt[0]=-phasespace::beta[0];
  bt[1]=-phasespace::beta[1];
  bt[2]=-phasespace::beta[2];
  
  double bdotk1=kap1[0]*bt[0]+kap1[1]*bt[1]+kap1[2]*bt[2];
  //zax[3]=gam*(kap1[3]-bdotk1); //bug bug bug!!! zax has only 3 elements!!!
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

//This function is not used, kept for cross checks
void omegaintegr::genzaxisCS()
{
  //CS frame prescription
  //take incoming partons, and boost them in the boson rest frame
  double bt[3];
  bt[0]=-phasespace::beta[0];
  bt[1]=-phasespace::beta[1];
  bt[2]=-phasespace::beta[2];

  double p1[3];
  p1[0] = 0;
  p1[1] = 0;
  p1[2] = 1;

  double p2[3];
  p2[0] = 0;
  p2[1] = 0;
  p2[2] = -1;

  double p1cm[3];
  double bdotp1=p1[0]*bt[0]+p1[1]*bt[1]+p1[2]*bt[2];
  p1cm[0]=p1[0]+phasespace::gam*bt[0]*(phasespace::gam/(phasespace::gam+1)*bdotp1-p1[3]);
  p1cm[1]=p1[1]+phasespace::gam*bt[1]*(phasespace::gam/(phasespace::gam+1)*bdotp1-p1[3]);
  p1cm[2]=p1[2]+phasespace::gam*bt[2]*(phasespace::gam/(phasespace::gam+1)*bdotp1-p1[3]);
  double p1cmabs = sqrt(pow(p1cm[0],2)+pow(p1cm[1],2)+pow(p1cm[2],2));

  double p2cm[3];
  double bdotp2=p2[0]*bt[0]+p2[1]*bt[1]+p2[2]*bt[2];
  p2cm[0]=p2[0]+phasespace::gam*bt[0]*(phasespace::gam/(phasespace::gam+1)*bdotp2-p2[3]);
  p2cm[1]=p2[1]+phasespace::gam*bt[1]*(phasespace::gam/(phasespace::gam+1)*bdotp2-p2[3]);
  p2cm[2]=p2[2]+phasespace::gam*bt[2]*(phasespace::gam/(phasespace::gam+1)*bdotp2-p2[3]);
  double p2cmabs = sqrt(pow(p2cm[0],2)+pow(p2cm[1],2)+pow(p2cm[2],2));

  //define the polar axis in CS frame, which is the axis bisecting the angle between p1 and p2 incoming partons
  zax[0] = p1cm[0]/p1cmabs - p2cm[0]/p2cmabs;
  zax[1] = p1cm[1]/p1cmabs - p2cm[1]/p2cmabs;
  zax[2] = p1cm[2]/p1cmabs - p2cm[2]/p2cmabs;

  //normalise z axis to unity
  double zaxabs = sqrt(pow(zax[0],2)+pow(zax[1],2)+pow(zax[2],2));
  zax[0] = zax[0] / zaxabs;
  zax[1] = zax[1] / zaxabs;
  zax[2] = zax[2] / zaxabs;
}

void omegaintegr::genl4p(double costh, double phi_lep)
//void omegaintegr::genl4p(float costh, float phi_lep) //quite significant speed up with float numbers for the sine cosine functions
{
  //Generate the 4-momenta of the leptons in the boson rest frame, using costh and phi_lep
  //p4 is the lepton and p3 is the antilepton
  
  phasespace::set_cth(costh);
  phasespace::set_philep(phi_lep);
  phasespace::genl4p();
  return;
    
  //In some cases the generation of the first lepton 4-momentum is trivial,
  //i.e. it is done with respect to the reference system x = (1,0,0); y = (0,1,0); z=(0,0,1)
  //a) With the cuba integration option, the phase space is generated uniformly
  //and costh is calculated using Catani's formulas arXiv:1507.0693 (25-32)
  //b) With the naive prescription, don't need to bother much since the rest frame axes
  //are the native x,y,z axes of the laboratory frame
  //c) When qt=0, because all reference frames are equivalent
  if (opts.cubaint || phasespace::qt == 0. || opts.qtrec_naive)
    phasespace::genl4p();
  else
    {
      //With all other general prescriptions, calculate the lepton 4-momenta p3 and p4
      //from the x,y,z axes of the rest frame as coordinate system
      
      //Start from the z axis, and rotate by angle theta with respect to y axis
      double rot1[3];
      double c = phasespace::costh;
      double s = sqrt(max(0.,1.-pow(c,2)));
      //rot1[0]=(c+yax[0]*yax[0]*(1-c))       *zax[0] + (yax[0]*yax[1]*(1-c)-yax[2]*s)*zax[1] + (yax[0]*yax[2]*(1-c)+yax[1]*s)*zax[2];
      //rot1[1]=(yax[1]*yax[0]*(1-c)+yax[2]*s)*zax[0] + (c+yax[1]*yax[1]*(1-c))       *zax[1] + (yax[1]*yax[2]*(1-c)-yax[0]*s)*zax[2];
      //rot1[2]=(yax[2]*yax[0]*(1-c)-yax[1]*s)*zax[0] + (yax[2]*yax[1]*(1-c)+yax[0]*s)*zax[1] + (c+yax[2]*yax[2]*(1-c))       *zax[2];
      phasespace::rotate(zax, c, s, yax, rot1);

      //rotate by angle phi_lep with respect to z axis
      double rot2[3];
      //c = cos(M_PI/2.-phi_lep); //why I was previously using this pi/2 rotation?
      //s = sin(M_PI/2.-phi_lep);
      c = cos(phasespace::phi_lep-M_PI);      //Use the same M_PI rotation convention as in kinematic
      s = sqrt(max(0.,1.-pow(c,2)))*((phasespace::phi_lep-M_PI)>0 ? 1 : -1);
      //c = cos(phi_lep-M_PI);
      //s = sin(phi_lep-M_PI);
      //rot2[0]=(c+zax[0]*zax[0]*(1-c))       *rot1[0] + (zax[0]*zax[1]*(1-c)-zax[2]*s)*rot1[1] + (zax[0]*zax[2]*(1-c)+zax[1]*s)*rot1[2];
      //rot2[1]=(zax[1]*zax[0]*(1-c)+zax[2]*s)*rot1[0] + (c+zax[1]*zax[1]*(1-c))       *rot1[1] + (zax[1]*zax[2]*(1-c)-zax[0]*s)*rot1[2];
      //rot2[2]=(zax[2]*zax[0]*(1-c)-zax[1]*s)*rot1[0] + (zax[2]*zax[1]*(1-c)+zax[0]*s)*rot1[1] + (c+zax[2]*zax[2]*(1-c))       *rot1[2];
      phasespace::rotate(rot1, c, s, zax, rot2);

      double p3cm[4];
      p3cm[3]=phasespace::m/2.;           //E
      p3cm[0]=p3cm[3]*rot2[0];            //px
      p3cm[1]=p3cm[3]*rot2[1];            //py
      p3cm[2]=p3cm[3]*rot2[2];            //pz

      /****************************
      //check that the generated p3cm is at costh from zax
      cout << "p4 generated in boson rest frame " << p3cm[0] << " " << p3cm[1] << " " << p3cm[2] << endl;
      cout << "angle between p3cm and polar axis (should be costh)" << (zax[0]*p3cm[0]+zax[1]*p3cm[1]+zax[2]*p3cm[2])
	/ sqrt(pow(zax[0],2)+pow(zax[1],2)+pow(zax[2],2))
	/ sqrt(pow(p3cm[0],2)+pow(p3cm[1],2)+pow(p3cm[2],2)) << "  " << costh << endl;
      ****************************/

      //Boost to go in the lab frame
      dyboost_(phasespace::gam, phasespace::beta, p3cm, phasespace::p3);

      //  momentum of the second lepton
      phasespace::p4[3]=phasespace::pV[3]-phasespace::p3[3];      //E
      phasespace::p4[0]=phasespace::pV[0]-phasespace::p3[0];      //px
      phasespace::p4[1]=phasespace::pV[1]-phasespace::p3[1];      //py
      phasespace::p4[2]=phasespace::pV[2]-phasespace::p3[2];      //pz
      
      /****************************
      //check:
      //generate a lepton p4 according to the naive prescription and
      //check that costh of Catani's formula is equal to
      //the angle between the polar axis of the rest frame and the lepton p4
      p3cm[3]=phasespace::m/2.;
      p3cm[0]=p3cm[3]*sin(acos(costh))*sin(phi_lep);
      p3cm[1]=p3cm[3]*sin(acos(costh))*cos(phi_lep);
      p3cm[2]=p3cm[3]*costh;
      double bdotp=p3cm[0]*beta[0]+p3cm[1]*beta[1]+p3cm[2]*beta[2];
      p4[3]=gam*(p3cm[3]-bdotp);
      p4[0]=p3cm[0]+gam*beta[0]*(gam/(gam+1)*bdotp-p3cm[3]);
      p4[1]=p3cm[1]+gam*beta[1]*(gam/(gam+1)*bdotp-p3cm[3]);
      p4[2]=p3cm[2]+gam*beta[2]*(gam/(gam+1)*bdotp-p3cm[3]);
      double costhcs = (zax[0]*p3cm[0]+zax[1]*p3cm[1]+zax[2]*p3cm[2])
	/ sqrt(pow(zax[0],2)+pow(zax[1],2)+pow(zax[2],2))
	/ sqrt(pow(p3cm[0],2)+pow(p3cm[1],2)+pow(p3cm[2],2));
      cout << endl;
      cout << "rapidity " << phasespace::y << " costh  " <<  costh 
	   << " costhcs calculated " << costhcs << " costhcs catani " << (1.-4.*(kap1[3]*p4[3]-kap1[2]*p4[2]-kap1[1]*p4[1]-kap1[0]*p4[0])/(phasespace::m*phasespace::m)) << endl;
      ****************************/
    }
  /****************************
  //check with kinematic calculation, valid for CS frame
  kinematic::set(p4, p3);
  kinematic::calc_vb();
  kinematic::calc_angles();
  cout << " generated  " << costh << "  " << phi_lep << endl;
  cout << " calculated " << kinematic::costh << "  " << kinematic::phi_lep << endl;
  cout << endl;
  ****************************/
}

//costh for any general kt prescription
double omegaintegr::costhCS()
{
  //change this to p3  !!! should use p3 here instead of p4 !!!
  return (1.-4.*(kap1[3]*phasespace::p4[3]-kap1[2]*phasespace::p4[2]-kap1[1]*phasespace::p4[1]-kap1[0]*phasespace::p4[0])/(phasespace::m*phasespace::m));
}

//Algorithm to determine all the subintervals of costh between -1 and +1 in which the lepton cuts are satisfied
void omegaintegr::costhbound(double phi_lep, vector<double> &min, vector<double> &max)
{
  if (!opts.makecuts)
    {
      min.push_back(phasespace::getcthmin());
      max.push_back(phasespace::getcthmax());
      return;
    }

  bool status;
  double c1 = phasespace::getcthmin();
  double c2 = phasespace::getcthmax();
  genl4p(c1, phi_lep);
  if (Kinematics::Cuts::KeepThisEvent(phasespace::p3, phasespace::p4))
    {
      min.push_back(c1);
      status = true;
    }
  else
    status = false;

  while (true)
    {
      if (status) //search for maximum
	{
	  double tempmax = phasespace::getcthmax()*1.1+0.1;
	  int nc = opts.ncstart;
	  for (int precision = 0; precision < 6; precision++)
	    {
	      double hc=(c2-c1)/nc;
	      for(int i=0;i<=nc;i++)
		{
		  double costh = i*hc+c1;
		  genl4p(costh, phi_lep);
		  if (!Kinematics::Cuts::KeepThisEvent(phasespace::p3, phasespace::p4))
		    {
		      tempmax = costh; //tempmax = costh_CS;
		      c2 = i*hc+c1;
		      c1 = (i-1)*hc+c1;
		      nc = 10;
		      break;
		    }
		}
	      if (tempmax > phasespace::getcthmax())
		{
		  max.push_back(phasespace::getcthmax());
		  break;
		}
	    }
	  if (tempmax > phasespace::getcthmax())
	    break;
	  max.push_back(tempmax);
	  status = false;
	  c1 = tempmax;
	  c2 = phasespace::getcthmax();
	}
      else //search for minimum
	{
	  double tempmin = phasespace::getcthmax()*1.1+0.1;
	  int nc = opts.ncstart;
	  for (int precision = 0; precision < 6; precision++)
	    {
	      double hc=(c2-c1)/nc;
	      for(int i=0;i<=nc;i++)
		{
		  double costh = i*hc+c1;
		  genl4p(costh, phi_lep);
		  if (Kinematics::Cuts::KeepThisEvent(phasespace::p3, phasespace::p4))
		    {
		      tempmin = costh; //tempmin = costh_CS;
		      c2 = i*hc+c1;
		      c1 = (i-1)*hc+c1;
		      nc = 10;
		      break;
		    }
		}
	      if (tempmin > phasespace::getcthmax())
		break;
	    }
	  if (tempmin > phasespace::getcthmax())
	    break;
	  min.push_back(tempmin);
	  status = true;
	  c1 = tempmin;
	  c2 = phasespace::getcthmax();
	}
    }
}

//Perform integration of costheta moments
void omegaintegr::cthmoments(double &cthmom0, double &cthmom1, double &cthmom2)
{
  clock_t begin_time, end_time;
  if (!opts.makecuts)
    {
      cthmom0 = phasespace::getcthmax()-phasespace::getcthmin();
      cthmom1 = (pow(phasespace::getcthmax(),2)-pow(phasespace::getcthmin(),2))/2.;
      cthmom2 = (pow(phasespace::getcthmax(),3)-pow(phasespace::getcthmin(),3))/3.;
      return;
    }

  //cuba integration (this allows to use CS prescription of DYRES)
  if (opts.cubaint)
    {
      begin_time = clock();
      const int nd = 2;   //dimensions of the integral
      const int nc = 3;  //components of the integrand
      void *userdata;
      const int nvec = 1;
      const int mineval = opts.suavepoints;
      const int maxeval = opts.suavepoints;
      const int nbatch = 10000;
      const double epsrel = 0.0;
      const double epsabs = 0.0;
      void *spin=NULL;
      int neval;
      int fail;
      double integral[3];
      double error[3];
      double prob[3];
      const char *statefile = "";

      //suave
      const int flags = 8;
      const int seed = 0;
      const int nnew = int(opts.suavepoints/10);
      const int nmin = int(opts.suavepoints/10);
      const double flatness = 100;
      int nregions;
      Suave(nd, nc,
	    (integrand_t) thphiintegrand, userdata, nvec,
	    epsrel, epsabs,
	    flags, seed,
	    mineval, maxeval,
	    nnew, nmin,
	    flatness, statefile, spin,
	    &nregions, &neval, &fail,
	    integral, error, prob);

      cubawait(&spin);
      end_time = clock();

      if (opts.verbose)
	cout << "cuba point " << setw(10) << phasespace::m  << setw(10) << phasespace::qt  << setw(10) << phasespace::y
	     << setw(13) << integral[0] << setw(13)  << error[0]
	     << setw(13) << integral[1] << setw(13)  << error[1]
	     << setw(13) << integral[2] << setw(13)  << error[2]
	     << "   time " << float( end_time - begin_time ) /  CLOCKS_PER_SEC << endl;

      cthmom0 = integral[0];
      cthmom1 = integral[1];
      cthmom2 = integral[2];
    }

  //trapezoidal rule and costh boundaries
  if (opts.trapezint)
    {
      begin_time = clock();
      double phi = 0.;
      double phi1 = 0.;
      double phi2 = 2. * M_PI;
      double hphi=(phi2-phi1)/opts.nphitrape;
      cthmom0 = 0.;
      cthmom1 = 0.;
      cthmom2 = 0.;
      vector <double> cthmin;
      vector <double> cthmax;
      vector<double>::iterator itmn;
      vector<double>::iterator itmx;
      for(int i=0;i<=opts.nphitrape;i++)
	{
	  double phi_lep = i*hphi+phi1;
	  cthmin.clear();
	  cthmax.clear();
	  costhbound(phi_lep, cthmin, cthmax);
	  if (i == 0 || i == opts.nphitrape)
	    {
	      itmn = cthmin.begin(); itmx = cthmax.begin();
	      for (; itmn != cthmin.end(); itmn++, itmx++)
		{
		  cthmom0 += (*itmx - *itmn);
		  cthmom1 += (pow(*itmx,2) - pow(*itmn,2)) / 2.;
		  cthmom2 += (pow(*itmx,3) - pow(*itmn,3)) / 3.;
		}
	    }
	  else
	    {
	      itmn = cthmin.begin(); itmx = cthmax.begin();
	      for (; itmn != cthmin.end(); itmn++, itmx++)
		{
		  cthmom0 += 2.*(*itmx - *itmn);
		  cthmom1 += 2.*(pow(*itmx,2) - pow(*itmn,2)) / 2.;
		  cthmom2 += 2.*(pow(*itmx,3) - pow(*itmn,3)) / 3.;
		}
	    }
	}
      cthmom0 = cthmom0 / opts.nphitrape / 2.;
      cthmom1 = cthmom1 / opts.nphitrape / 2.;
      cthmom2 = cthmom2 / opts.nphitrape / 2.;
      end_time = clock();

      if (opts.verbose)
	cout << "trap point " << setw(10) << phasespace::m  << setw(10) << phasespace::qt  << setw(10) << phasespace::y
	     << setw(13) << cthmom0 << setw(13)  << 0
	     << setw(13) << cthmom1 << setw(13)  << 0
	     << setw(13) << cthmom2 << setw(13)  << 0
	     << "   time " << float( end_time - begin_time ) /  CLOCKS_PER_SEC << endl;
    }


  //quadrature rule and costh boundaries
  if (opts.quadint)
    {
      begin_time = clock();
      double phi1 = 0.;
      double phi2 = 2. * M_PI;
      double hphi=(phi2-phi1)/opts.phiintervals;
      cthmom0 = 0.;
      cthmom1 = 0.;
      cthmom2 = 0.;
      vector <double> cthmin;
      vector <double> cthmax;
      vector<double>::iterator itmn;
      vector<double>::iterator itmx;

      //for qt == 0 there is no phi dependence and the integration is reduced by one dimension, to only dcosth
      if (phasespace::qt == 0.)
	{
	  costhbound(0., cthmin, cthmax);
	  itmn = cthmin.begin();
	  itmx = cthmax.begin();
	  for (; itmn != cthmin.end(); itmn++, itmx++)
	    {
	      cthmom0 += (*itmx - *itmn);
	      cthmom1 += (pow(*itmx,2) - pow(*itmn,2)) / 2.;
	      cthmom2 += (pow(*itmx,3) - pow(*itmn,3)) / 3.;
	    }
	}
      else
	{
	  for(int i=0; i < opts.phiintervals; i++)
	    {
	      double phimin = i*hphi+phi1;
	      double phimax = (i+1)*hphi+phi1;
	      
	      double phic=0.5*(phimin+phimax);
	      double phim=0.5*(phimax-phimin);
	      
	      for(int i=0; i < opts.phirule; i++)
		{
		  double phi_lep = phic+phim*gr::xxx[opts.phirule-1][i];
		  cthmin.clear();
		  cthmax.clear();
		  costhbound(phi_lep, cthmin, cthmax);
		  itmn = cthmin.begin(); itmx = cthmax.begin();
		  for (; itmn != cthmin.end(); itmn++, itmx++)
		    {
		      cthmom0 += (*itmx - *itmn)*gr::www[opts.phirule-1][i]*phim;
		      cthmom1 += (pow(*itmx,2) - pow(*itmn,2)) / 2.*gr::www[opts.phirule-1][i]*phim;
		      cthmom2 += (pow(*itmx,3) - pow(*itmn,3)) / 3.*gr::www[opts.phirule-1][i]*phim;
		    }
		}
	    }
	  cthmom0 /= (phi2-phi1);
	  cthmom1 /= (phi2-phi1);
	  cthmom2 /= (phi2-phi1);
	}
      end_time = clock();

      if (opts.verbose)
	cout << "quad point " << setw(10) << phasespace::m  << setw(10) << phasespace::qt  << setw(10) << phasespace::y
	     << setw(13) << cthmom0 << setw(13)  << 0
	     << setw(13) << cthmom1 << setw(13)  << 0
	     << setw(13) << cthmom2 << setw(13)  << 0
	     << "   time " << float( end_time - begin_time ) /  CLOCKS_PER_SEC << endl;
    }
}


integrand_t omegaintegr::thphiintegrand(const int &ndim, const double x[], const int &ncomp, double f[])
{
  double costh=-1.+2.*x[0];
  double jac=2.;

  //Integral in 2 dimensions: costh and phi_lep
  double phi_lep = x[1] * 2. * M_PI;

  /*
  //Integral in 3 dimensions: costh, phiZ and phi_lep
  double phi_lep = x[1] * 2. * M_PI;   // direction of first lepton decay in the center of mass frame 
  double phi = x[2] * 2. * M_PI;   // Vector boson azimuthal angle
  //generate boson 4-momentum, with m, qt, y and phi
  genV4p(m, qt, y, phi);
  */

  genl4p(costh, phi_lep);
  if (Kinematics::Cuts::KeepThisEvent(phasespace::p3, phasespace::p4))
    {
      //evaluate couplings using cos(theta) in CS framework
      double costh_CS = costhCS();
      f[0] = 1.*jac;
      f[1] = costh_CS*jac;
      f[2] = pow(costh_CS, 2)*jac;
      //      f[0] = 1.*jac;
      //      f[1] = costh*jac;
      //      f[2] = pow(costh, 2)*jac;
    }
  else
    {
      f[0] = 0.;
      f[1] = 0.;
      f[2] = 0.;
    }
  return 0;
}


double omegaintegr::costh_qtrec()
{
  //Calculate costh according to a qt-recoil prescription
  phasespace::restframeid RF;
  if (opts.qtrec_naive) RF = phasespace::naive;
  else if (opts.qtrec_cs) RF = phasespace::CS;
  else if (opts.qtrec_kt0) RF = phasespace::kt0;

  double kt1,kt2;
  double m2 = phasespace::m*phasespace::m;
  
  //CS frame prescription
  if (RF == phasespace::CS)
    {
      kt1 = phasespace::pV[0]/2.;
      kt2 = phasespace::pV[1]/2.;
    }

  //MY (DYRES) prescription
  if (RF == phasespace::naive)
    {
      kt1=(1.+phasespace::pV[2]/(sqrt(m2)+phasespace::pV[3]))*phasespace::pV[0]/2;
      kt2=(1.+phasespace::pV[2]/(sqrt(m2)+phasespace::pV[3]))*phasespace::pV[1]/2;
    }  

  //alternative k1t = 0 prescription
  if (RF == phasespace::kt0)
    {
      kt1 = 0;
      kt2 = 0;
    }
  
  double mt2 = phasespace::m*phasespace::m + phasespace::qt*phasespace::qt;
  double zeta1=1./m2/2.*(m2+2.*(phasespace::pV[0]*kt1+phasespace::pV[1]*kt2)+sqrt(pow((m2+2.*(phasespace::pV[0]*kt1+phasespace::pV[1]*kt2)),2)-4.*mt2*(pow(kt1,2)+pow(kt2,2))));
  
  double qP1=(phasespace::pV[3]-phasespace::pV[2])*opts.sroot/2.;
  double qP2=(phasespace::pV[3]+phasespace::pV[2])*opts.sroot/2.;

  double kap1[4];
  kap1[3]=opts.sroot/2.*(zeta1*m2/2./qP1+(pow(kt1,2)+pow(kt2,2))/zeta1*qP1/m2/pow(opts.sroot,2)*2.);
  kap1[0]=kt1;
  kap1[1]=kt2;
  kap1[2]=opts.sroot/2.*(zeta1*m2/2./qP1-(pow(kt1,2)+pow(kt2,2))/zeta1*qP1/m2/pow(opts.sroot,2)*2.);

  return 1.-4.*(kap1[3]*phasespace::p3[3]-kap1[2]*phasespace::p3[2]-kap1[1]*phasespace::p3[1]-kap1[0]*phasespace::p3[0])/m2;
}
