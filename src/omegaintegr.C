// --------- TO DO -----------------
//make this collection of functions a self-consistent library
//by removing dependences on settings.h and interface.h,
//introducing a namespace, and adding a init function for initialisation of settings
// ---------------------------------

#include "omegaintegr.h"
#include "settings.h"
#include "interface.h"
#include "phasespace.h"
#include "cuba.h"

#include <ctime>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <vector>

//Vector boson 4-momentum and boost
double omegaintegr::pV[4];
double omegaintegr::gam;
double omegaintegr::beta[3];

//leptons 4-momenta
double omegaintegr::p4[4];
double omegaintegr::p3[4];

//rest frame axes
double omegaintegr::kap1[4];
double omegaintegr::xax[3];
double omegaintegr::yax[3];
double omegaintegr::zax[3];

//quadrature rules
const double xq4[4]={-0.3399810435848563,0.3399810435848563,-0.8611363115940526,0.8611363115940526};
const double wq4[4]={0.6521451548625461,0.6521451548625461,0.3478548451374538,0.3478548451374538};

const double xq8[8]={-0.1834346424956498,0.1834346424956498,-0.5255324099163290,0.5255324099163290,-0.7966664774136267,0.7966664774136267,-0.9602898564975363,0.9602898564975363};
const double wq8[8]={0.3626837833783620,0.3626837833783620,0.3137066458778873,0.3137066458778873,0.2223810344533745,0.2223810344533745,0.1012285362903763,0.1012285362903763};

//fortran interface
void genv4p_()
{
  omegaintegr::genV4p();
}

void cthmoments_(double &cthmom0, double &cthmom1, double &cthmom2)
{
  omegaintegr::cthmoments(cthmom0, cthmom1, cthmom2);
}


//instead of passing m, qt, y, and phi as variables, take them from the phasespace:: namelist
void omegaintegr::genV4p() //double m, double qt, double y, double phi)
{
  //Generate the boson 4-momentum
  double m2=phasespace::m*phasespace::m;
  double exppy=exp(phasespace::y);
  double expmy=1./exppy;
  double qt2=pow(phasespace::qt,2);
  double mt2=m2+qt2;

  //vector boson momentum: pV[3]^2-pV[0]^2-pV[1]^2-pV[2]^2=m2
  pV[0]=phasespace::qt*cos(phasespace::phiV);                  //px
  pV[1]=phasespace::qt*sin(phasespace::phiV);                  //py
  pV[2]=0.5*sqrt(mt2)*(exppy-expmy);   //pz
  pV[3]=0.5*sqrt(mt2)*(exppy+expmy);   //E

  //Calculte the boost 4-vector from the lab frame to the rest frame
  gam=pV[3]/phasespace::m;
  beta[0]=-pV[0]/pV[3];
  beta[1]=-pV[1]/pV[3];
  beta[2]=-pV[2]/pV[3];

  //if qt=0, as in fixed order predictions, don't bother calculating axes of the rest frame (all rest frames are equivalent at qt=0)
  if (phasespace::qt == 0.)
    return;
  
  //recoil prescriptions are defined by the values of kt1 and kt2
  double kt1, kt2;

  //CS frame prescription
  if (opts.qtrec_cs)
    {
      kt1 = pV[0]/2.;
      kt2 = pV[1]/2.;
    }

  //naive prescription (also called MY prescription in DYRES)
  if (opts.qtrec_naive)
    {
      kt1=(1.+pV[2]/(sqrt(m2)+pV[3]))*pV[0]/2;
      kt2=(1.+pV[2]/(sqrt(m2)+pV[3]))*pV[1]/2;
    }  

  //alternative k1t = 0 prescription
  if (opts.qtrec_kt0)
    {
      kt1 = 0;
      kt2 = 0;
    }
  
  //zeta1 as in Eq.(26) of arXiv:1507.0693
  double zeta1=1./m2/2.*(m2+2.*(pV[0]*kt1+pV[1]*kt2)+sqrt(pow((m2+2.*(pV[0]*kt1+pV[1]*kt2)),2)-4.*mt2*(pow(kt1,2)+pow(kt2,2))));
  double qP1=(pV[3]-pV[2])*energy_.sroot_/2.;
  double qP2=(pV[3]+pV[2])*energy_.sroot_/2.;
  //kap1 is the colliding parton a1 in the boson rest frame
  kap1[3]=energy_.sroot_/2.*(zeta1*m2/2./qP1+(pow(kt1,2)+pow(kt2,2))/zeta1*qP1/m2/pow(energy_.sroot_,2)*2.);
  kap1[0]=kt1;
  kap1[1]=kt2;
  kap1[2]=energy_.sroot_/2.*(zeta1*m2/2./qP1-(pow(kt1,2)+pow(kt2,2))/zeta1*qP1/m2/pow(energy_.sroot_,2)*2.);

  /*
  if (opts.qtrec_cs)   //CS frame prescription
    {
      //take incoming partons, and boost them in the boson rest frame
      double bt[3];
      bt[0]=-beta[0];
      bt[1]=-beta[1];
      bt[2]=-beta[2];

      double p1[4];
      p1[3] = 0.5*m*exppy;
      p1[0] = 0;
      p1[1] = 0;
      p1[2] = 0.5*m*exppy;

      double p2[4];
      p2[3] = 0.5*m*expmy;
      p2[0] = 0;
      p2[1] = 0;
      p2[2] = -0.5*m*expmy;

      double p1cm[4];
      double bdotp1=p1[0]*bt[0]+p1[1]*bt[1]+p1[2]*bt[2];
      p1cm[3]=gam*(p1[3]-bdotp1);
      p1cm[0]=p1[0]+gam*bt[0]*(gam/(gam+1)*bdotp1-p1[3]);
      p1cm[1]=p1[1]+gam*bt[1]*(gam/(gam+1)*bdotp1-p1[3]);
      p1cm[2]=p1[2]+gam*bt[2]*(gam/(gam+1)*bdotp1-p1[3]);

      double p2cm[4];
      double bdotp2=p2[0]*bt[0]+p2[1]*bt[1]+p2[2]*bt[2];
      p2cm[3]=gam*(p2[3]-bdotp2);
      p2cm[0]=p2[0]+gam*bt[0]*(gam/(gam+1)*bdotp2-p2[3]);
      p2cm[1]=p2[1]+gam*bt[1]*(gam/(gam+1)*bdotp2-p2[3]);
      p2cm[2]=p2[2]+gam*bt[2]*(gam/(gam+1)*bdotp2-p2[3]);

      //define the polar axis in CS frame, which is the axis bisecting the angle between p1 and p2 incoming partons
      double p1cmabs = sqrt(pow(p1cm[0],2)+pow(p1cm[1],2)+pow(p1cm[2],2));
      double p2cmabs = sqrt(pow(p2cm[0],2)+pow(p2cm[1],2)+pow(p2cm[2],2));
      zax[0] = p1cm[0]/p1cmabs - p2cm[0]/p2cmabs;
      zax[1] = p1cm[1]/p1cmabs - p2cm[1]/p2cmabs;
      zax[2] = p1cm[2]/p1cmabs - p2cm[2]/p2cmabs;
    }
  */
  
  //Determine the x,y,z axes of the boson rest frame for any general kt1 prescription
  //Calculate z axis by boosting kap1 from the laboratory to the rest frame
  double bt[3];
  bt[0]=-beta[0];
  bt[1]=-beta[1];
  bt[2]=-beta[2];
  
  double bdotk1=kap1[0]*bt[0]+kap1[1]*bt[1]+kap1[2]*bt[2];
  zax[3]=gam*(kap1[3]-bdotk1); //bug bug bug!!! zax has only 3 elements!!!
  zax[0]=kap1[0]+gam*bt[0]*(gam/(gam+1)*bdotk1-kap1[3]);
  zax[1]=kap1[1]+gam*bt[1]*(gam/(gam+1)*bdotk1-kap1[3]);
  zax[2]=kap1[2]+gam*bt[2]*(gam/(gam+1)*bdotk1-kap1[3]);

  //normalise z axis to unity
  double zaxabs = sqrt(pow(zax[0],2)+pow(zax[1],2)+pow(zax[2],2));
  zax[0] = zax[0] / zaxabs;
  zax[1] = zax[1] / zaxabs;
  zax[2] = zax[2] / zaxabs;

  //define y axis ortogonal to z and V direction
  yax[0] = zax[1]*pV[2] - zax[2]*pV[1];
  yax[1] = zax[2]*pV[0] - zax[0]*pV[2];
  yax[2] = zax[0]*pV[1] - zax[1]*pV[0];
  
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
  cout << pV[0] << " " << pV[1] << "  " << pV[2] << endl;
  */
}

void omegaintegr::genl4p(double costh, double phi_lep)
//void omegaintegr::genl4p(float costh, float phi_lep) //quite significant speed up with float numbers for the sine cosine functions
{
  //phi_lep = M_PI/2.;

  //Generate the 4-momentum of the first lepton in the boson rest frame, using costh and phi_lep
  double p4cm[4];

  //In some cases the generation of the first lepton 4-momentum is trivial,
  //i.e. it is done with respect to the reference system x = (1,0,0); y = (0,1,0); z=(0,0,1)
  //a) With the cuba integration option, the phase space is generated uniformly
  //and costh is calculated using Catani's formulas arXiv:1507.0693 (25-32)
  //b) With the naive prescription, don't need to bother much since the rest frame axes
  //are the native x,y,z axes of the laboratory frame
  //c) When qt=0, because all reference frames are equivalent
  if (opts.cubaint || phasespace::qt == 0. || opts.qtrec_naive)
    {
      double sintheta = sqrt(max(0.,1.-pow(costh,2)));
      double cosphi = cos(phi_lep);
      double sinphi = sqrt(max(0.,1.-pow(cosphi,2)));
      p4cm[3]=phasespace::m/2.;         //E
      p4cm[0]=p4cm[3]*sintheta*sinphi;  //px
      p4cm[1]=p4cm[3]*sintheta*cosphi;  //py
      p4cm[2]=p4cm[3]*costh;            //pz
    }
  else
    {
      //with all other general prescriptions, calculate the lepton 4-momenta p3 and p4
      //from the x,y,z axes of the rest frame as coordinate system
      
      /****************************/
      /*
	cout << endl;
	//check:
	//generate a lepton p4 according to the naive prescription and
	//check that costh of Catani's formula is equal to
	//the angle between the polar axis of the rest frame and the lepton p4
	p4cm[3]=phasespace::m/2.;
	p4cm[0]=p4cm[3]*sin(acos(costh))*sin(phi_lep);
	p4cm[1]=p4cm[3]*sin(acos(costh))*cos(phi_lep);
	p4cm[2]=p4cm[3]*costh;
	double bdotp=p4cm[0]*beta[0]+p4cm[1]*beta[1]+p4cm[2]*beta[2];
	p4[3]=gam*(p4cm[3]-bdotp);
	p4[0]=p4cm[0]+gam*beta[0]*(gam/(gam+1)*bdotp-p4cm[3]);
	p4[1]=p4cm[1]+gam*beta[1]*(gam/(gam+1)*bdotp-p4cm[3]);
	p4[2]=p4cm[2]+gam*beta[2]*(gam/(gam+1)*bdotp-p4cm[3]);
	double costhcs = (zax[0]*p4cm[0]+zax[1]*p4cm[1]+zax[2]*p4cm[2])
	/ sqrt(pow(zax[0],2)+pow(zax[1],2)+pow(zax[2],2))
	/ sqrt(pow(p4cm[0],2)+pow(p4cm[1],2)+pow(p4cm[2],2));
	cout << "rapidity " << phasespace::y << " costh  " <<  costh 
	<< " costhcs calculated " << costhcs << " costhcs catani " << (1.-4.*(kap1[3]*p4[3]-kap1[2]*p4[2]-kap1[1]*p4[1]-kap1[0]*p4[0])/(phasespace::m*phasespace::m)) << endl;
      */
      /****************************/


      //Start from the z axis, and rotate by angle theta with respect to y axis
      double rot1[3];
      double c = costh;
      //double s = sin(acos(costh));
      double s = sqrt(max(0.,1.-pow(c,2)));
      rot1[0]=(c+yax[0]*yax[0]*(1-c))       *zax[0] + (yax[0]*yax[1]*(1-c)-yax[2]*s)*zax[1] + (yax[0]*yax[2]*(1-c)+yax[1]*s)*zax[2];
      rot1[1]=(yax[1]*yax[0]*(1-c)+yax[2]*s)*zax[0] + (c+yax[1]*yax[1]*(1-c))       *zax[1] + (yax[1]*yax[2]*(1-c)-yax[0]*s)*zax[2];
      rot1[2]=(yax[2]*yax[0]*(1-c)-yax[1]*s)*zax[0] + (yax[2]*yax[1]*(1-c)+yax[0]*s)*zax[1] + (c+yax[2]*yax[2]*(1-c))       *zax[2];

      //rotate by angle phi_lep with respect to z axis
      double rot2[3];
      c = cos(M_PI/2.-phi_lep);
      //s = sin(M_PI/2.-phi_lep);
      s = sqrt(max(0.,1.-pow(c,2)));
      rot2[0]=(c+zax[0]*zax[0]*(1-c))       *rot1[0] + (zax[0]*zax[1]*(1-c)-zax[2]*s)*rot1[1] + (zax[0]*zax[2]*(1-c)+zax[1]*s)*rot1[2];
      rot2[1]=(zax[1]*zax[0]*(1-c)+zax[2]*s)*rot1[0] + (c+zax[1]*zax[1]*(1-c))       *rot1[1] + (zax[1]*zax[2]*(1-c)-zax[0]*s)*rot1[2];
      rot2[2]=(zax[2]*zax[0]*(1-c)-zax[1]*s)*rot1[0] + (zax[2]*zax[1]*(1-c)+zax[0]*s)*rot1[1] + (c+zax[2]*zax[2]*(1-c))       *rot1[2];

      p4cm[3]=phasespace::m/2.;           //E
      p4cm[0]=p4cm[3]*rot2[0];            //px
      p4cm[1]=p4cm[3]*rot2[1];            //py
      p4cm[2]=p4cm[3]*rot2[2];            //pz
  
      /*
      //check that the generated p4cm is at costh from zax
      cout << "p4 generated in boson rest frame " << p4cm[0] << " " << p4cm[1] << " " << p4cm[2] << endl;
      cout << "angle between p4cm and polar axis (should be costh)" << (zax[0]*p4cm[0]+zax[1]*p4cm[1]+zax[2]*p4cm[2])
      / sqrt(pow(zax[0],2)+pow(zax[1],2)+pow(zax[2],2))
      / sqrt(pow(p4cm[0],2)+pow(p4cm[1],2)+pow(p4cm[2],2)) << "  " << costh << endl;
      */
    }

  // Boost to go in the lab frame
  double bdotp=p4cm[0]*beta[0]+p4cm[1]*beta[1]+p4cm[2]*beta[2];
  p4[3]=gam*(p4cm[3]-bdotp);
  p4[0]=p4cm[0]+gam*beta[0]*(gam/(gam+1)*bdotp-p4cm[3]);
  p4[1]=p4cm[1]+gam*beta[1]*(gam/(gam+1)*bdotp-p4cm[3]);
  p4[2]=p4cm[2]+gam*beta[2]*(gam/(gam+1)*bdotp-p4cm[3]);

  //  momentum of the second lepton
  p3[3]=pV[3]-p4[3];      //E
  p3[0]=pV[0]-p4[0];      //px
  p3[1]=pV[1]-p4[1];      //py
  p3[2]=pV[2]-p4[2];      //pz
}

//CS FRAME PRESCRIPTION
double omegaintegr::costhCS()
{
  return (1.-4.*(kap1[3]*p4[3]-kap1[2]*p4[2]-kap1[1]*p4[1]-kap1[0]*p4[0])/(phasespace::m*phasespace::m));
}

//Algorithm to determine all the subintervals of costh between -1 and +1 in which the lepton cuts are satisfied
void omegaintegr::costhbound(double phi_lep, vector<double> &min, vector<double> &max)
{
  if (!opts.makelepcuts)
    {
      min.push_back(-1);
      max.push_back(+1);
      return;
    }

  bool status;
  double c1 = -1;
  double c2 = +1;
  genl4p(c1, phi_lep);
  if (cuts::lep(p3, p4))
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
	  double tempmax = +1.1;
	  int nc = opts.ncstart;
	  for (int precision = 0; precision < 6; precision++)
	    {
	      double hc=(c2-c1)/nc;
	      for(int i=0;i<=nc;i++)
		{
		  double costh = i*hc+c1;
		  genl4p(costh, phi_lep);
		  if (!cuts::lep(p3, p4))
		    {
		      tempmax = costh; //tempmax = costh_CS;
		      c2 = i*hc+c1;
		      c1 = (i-1)*hc+c1;
		      nc = 10;
		      break;
		    }
		}
	      if (tempmax > +1)
		{
		  max.push_back(1);
		  break;
		}
	    }
	  if (tempmax > +1)
	    break;
	  max.push_back(tempmax);
	  status = false;
	  c1 = tempmax;
	  c2 = +1;
	}
      else //search for minimum
	{
	  double tempmin = +1.1;
	  int nc = opts.ncstart;
	  for (int precision = 0; precision < 6; precision++)
	    {
	      double hc=(c2-c1)/nc;
	      for(int i=0;i<=nc;i++)
		{
		  double costh = i*hc+c1;
		  genl4p(costh, phi_lep);
		  if (cuts::lep(p3, p4))
		    {
		      tempmin = costh; //tempmin = costh_CS;
		      c2 = i*hc+c1;
		      c1 = (i-1)*hc+c1;
		      nc = 10;
		      break;
		    }
		}
	      if (tempmin > +1)
		break;
	    }
	  if (tempmin > +1)
	    break;
	  min.push_back(tempmin);
	  status = true;
	  c1 = tempmin;
	  c2 = +1;
	}
    }
}

//Perform integration of costheta moments
void omegaintegr::cthmoments(double &cthmom0, double &cthmom1, double &cthmom2)
{
  clock_t begin_time, end_time;
  if (!opts.makelepcuts)
    {
      cthmom0 = 2.;
      cthmom1 = 0.;
      cthmom2 = 2./3.;
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
      double phi = 0.;
      double phi1 = 0.;
      double phi2 = 2. * M_PI;
      double hphi=(phi2-phi1)/opts.quadnphi;
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
	  for(int i=0; i< opts.quadnphi; i++)
	    {
	      double phimin = i*hphi+phi1;
	      double phimax = (i+1)*hphi+phi1;
	      
	      double phic=0.5*(phimin+phimax);
	      double phim=0.5*(phimax-phimin);
	      
	      for(int i=0; i<4; i++)
		{
		  double phi_lep = phic+phim*xq4[i];
		  cthmin.clear();
		  cthmax.clear();
		  costhbound(phi_lep, cthmin, cthmax);
		  itmn = cthmin.begin(); itmx = cthmax.begin();
		  for (; itmn != cthmin.end(); itmn++, itmx++)
		    {
		      cthmom0 += (*itmx - *itmn)*wq4[i]*phim;
		      cthmom1 += (pow(*itmx,2) - pow(*itmn,2)) / 2.*wq4[i]*phim;
		      cthmom2 += (pow(*itmx,3) - pow(*itmn,3)) / 3.*wq4[i]*phim;
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
  if (cuts::lep(p3, p4))
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

void omegaintegr::getp3(double p[4])
{
  p[0] = p3[0];
  p[1] = p3[1];
  p[2] = p3[2];
  p[3] = p3[3];
}
void omegaintegr::getp4(double p[4])
{
  p[0] = p4[0];
  p[1] = p4[1];
  p[2] = p4[2];
  p[3] = p4[3];
}
