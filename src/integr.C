#include "integr.h"
#include "settings.h"
#include "interface.h"
#include "plotter.h"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <ctime>
#include <math.h>
#include <iomanip>
#include <vector>
#include <random>

//integration boundaries
double mmin, mmax;
double qtmin, qtmax;
double ymin, ymax;

//global variables for phase space
double _m, _qt, _y, _costh;

//Vector boson 4-momentum and boost
double pV[4];
double gam;
double beta[3];

//leptons 4-momenta
double p4[4];
double p3[4];

//CS framework
double kap1[4];
double xax[3];
double yax[3];
double zax[3];

//quadrature rules
const double xq4[4]={-0.3399810435848563,0.3399810435848563,-0.8611363115940526,0.8611363115940526};
const double wq4[4]={0.6521451548625461,0.6521451548625461,0.3478548451374538,0.3478548451374538};

const double xq8[8]={-0.1834346424956498,0.1834346424956498,-0.5255324099163290,0.5255324099163290,-0.7966664774136267,0.7966664774136267,-0.9602898564975363,0.9602898564975363};
const double wq8[8]={0.3626837833783620,0.3626837833783620,0.3137066458778873,0.3137066458778873,0.2223810344533745,0.2223810344533745,0.1012285362903763,0.1012285362903763};

integrand_t thphiintegrand(const int &ndim, const double x[], const int &ncomp, double f[]);


void setbounds(double m1, double m2, double qt1, double qt2, double y1, double y2)
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

void set(double costh, double mm, double qtt, double yy)
{
  _costh = costh;
  _m = mm;
  _qt = qtt;
  _y = yy;
}

void set(double mm, double qtt, double yy)
{
  _m = mm;
  _qt = qtt;
  _y = yy;
}
void genV4p(double m, double qt, double y, double phi)
{
  //Generate the boson invariant mass
  double m2=m*m;
  double expy=exp(y);
  double expmy=exp(-y);
  double qt2=pow(qt,2);
  double mt2=m2+qt2;

  //vector boson momentum: pV[3]^2-pV[0]^2-pV[1]^2-pV[2]^2=m2
  pV[0]=qt*cos(phi);
  pV[1]=qt*sin(phi);
  pV[2]=0.5*sqrt(mt2)*(expy-expmy);
  pV[3]=0.5*sqrt(mt2)*(expy+expmy);

  //boost
  gam=pV[3]/m;
  beta[0]=-pV[0]/pV[3];
  beta[1]=-pV[1]/pV[3];
  beta[2]=-pV[2]/pV[3];

  //recoil prescriptions
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
  
  double zeta1=1./m2/2.*(m2+2.*(pV[0]*kt1+pV[1]*kt2)+sqrt(pow((m2+2.*(pV[0]*kt1+pV[1]*kt2)),2)-4.*mt2*(pow(kt1,2)+pow(kt2,2))));
  double qP1=(pV[3]-pV[2])*energy_.sroot_/2.;
  double qP2=(pV[3]+pV[2])*energy_.sroot_/2.;
  kap1[3]=energy_.sroot_/2.*(zeta1*m2/2./qP1+(pow(kt1,2)+pow(kt2,2))/zeta1*qP1/m2/pow(energy_.sroot_,2)*2.);
  kap1[0]=kt1;
  kap1[1]=kt2;
  kap1[2]=energy_.sroot_/2.*(zeta1*m2/2./qP1-(pow(kt1,2)+pow(kt2,2))/zeta1*qP1/m2/pow(energy_.sroot_,2)*2.);

  zax[0] = 0.;
  zax[1] = 0.;
  zax[2] = 0.;

  /*
  if (opts.qtrec_cs)   //CS frame prescription
    {
      //take incoming partons, and boost them in the boson rest frame
      double bt[3];
      bt[0]=-beta[0];
      bt[1]=-beta[1];
      bt[2]=-beta[2];

      double p1[4];
      p1[3] = 0.5*m*expy;
      p1[0] = 0;
      p1[1] = 0;
      p1[2] = 0.5*m*expy;

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
  
  /************/
  //General kt1 prescription: Calculate z axis by boosting kap1 from the laboratory to the rest frame
  double bt[3];
  bt[0]=-beta[0];
  bt[1]=-beta[1];
  bt[2]=-beta[2];
  
  double bdotk1=kap1[0]*bt[0]+kap1[1]*bt[1]+kap1[2]*bt[2];
  zax[3]=gam*(kap1[3]-bdotk1);
  zax[0]=kap1[0]+gam*bt[0]*(gam/(gam+1)*bdotk1-kap1[3]);
  zax[1]=kap1[1]+gam*bt[1]*(gam/(gam+1)*bdotk1-kap1[3]);
  zax[2]=kap1[2]+gam*bt[2]*(gam/(gam+1)*bdotk1-kap1[3]);
  /************/

  //normalise to unity
  double zaxabs = sqrt(pow(zax[0],2)+pow(zax[1],2)+pow(zax[2],2));
  zax[0] = zax[0] / zaxabs;
  zax[1] = zax[1] / zaxabs;
  zax[2] = zax[2] / zaxabs;

  //define y axis ortogonal to z and V direction
  yax[0] = zax[1]*pV[2] - zax[2]*pV[1];
  yax[1] = zax[2]*pV[0] - zax[0]*pV[2];
  yax[2] = zax[0]*pV[1] - zax[1]*pV[0];
  //normalise to unity
  double yaxabs = sqrt(pow(yax[0],2)+pow(yax[1],2)+pow(yax[2],2));
  yax[0] = yax[0] / yaxabs;
  yax[1] = yax[1] / yaxabs;
  yax[2] = yax[2] / yaxabs;

  //define x axis ortogonal to y and z direction
  xax[0] = yax[1]*zax[2] - yax[2]*zax[1];
  xax[1] = yax[2]*zax[0] - yax[0]*zax[2];
  xax[2] = yax[0]*zax[1] - yax[1]*zax[0];
  //normalise to unity
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
void genv4p_()
{
  genV4p(_m, _qt, _y, 0);
}

inline void genl4p(float costh, float phi_lep) //quite significant speed up with float numbers for the sine cosine functions
{
  //phi_lep = M_PI/2.;
  // momentum of the first lepton in the boson rest frame
  double p4cm[4];
  if (opts.cubaint)
    {
      //with the cuba integration option, the phase space is generated uniformly
      //and costh is calculated using Catani's formulas arXiv:1507.0693 (25-32)
      p4cm[3]=_m/2.;
      p4cm[0]=p4cm[3]*sin(acos(costh))*sin(phi_lep);
      p4cm[1]=p4cm[3]*sin(acos(costh))*cos(phi_lep);
      p4cm[2]=p4cm[3]*costh;
    }
  else
    if (opts.qtrec_naive)
      {
	//with the naive prescription, don't need to bother much since the rest frame axes
	//are the native x,y,z axes of the laboratory frame
	p4cm[3]=_m/2.;
	p4cm[0]=p4cm[3]*sin(acos(costh))*sin(phi_lep);
	p4cm[1]=p4cm[3]*sin(acos(costh))*cos(phi_lep);
	p4cm[2]=p4cm[3]*costh;
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
	  p4cm[3]=_m/2.;
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
	  cout << "rapidity " << _y << " costh  " <<  costh 
	  << " costhcs calculated " << costhcs << " costhcs catani " << (1.-4.*(kap1[3]*p4[3]-kap1[2]*p4[2]-kap1[1]*p4[1]-kap1[0]*p4[0])/(_m*_m)) << endl;
	*/
	/****************************/


	//Start from the z axis, and rotate by angle theta with respect to y axis
	double rot1[3];
	double c = costh;
	double s = sin(acos(costh));
	rot1[0]=(c+yax[0]*yax[0]*(1-c))       *zax[0] + (yax[0]*yax[1]*(1-c)-yax[2]*s)*zax[1] + (yax[0]*yax[2]*(1-c)+yax[1]*s)*zax[2];
	rot1[1]=(yax[1]*yax[0]*(1-c)+yax[2]*s)*zax[0] + (c+yax[1]*yax[1]*(1-c))       *zax[1] + (yax[1]*yax[2]*(1-c)-yax[0]*s)*zax[2];
	rot1[2]=(yax[2]*yax[0]*(1-c)-yax[1]*s)*zax[0] + (yax[2]*yax[1]*(1-c)+yax[0]*s)*zax[1] + (c+yax[2]*yax[2]*(1-c))       *zax[2];
  
	//rotate by angle phi_lep with respect to z axis
	double rot2[3];
	c = cos(M_PI/2.-phi_lep);
	s = sin(M_PI/2.-phi_lep);
	rot2[0]=(c+zax[0]*zax[0]*(1-c))       *rot1[0] + (zax[0]*zax[1]*(1-c)-zax[2]*s)*rot1[1] + (zax[0]*zax[2]*(1-c)+zax[1]*s)*rot1[2];
	rot2[1]=(zax[1]*zax[0]*(1-c)+zax[2]*s)*rot1[0] + (c+zax[1]*zax[1]*(1-c))       *rot1[1] + (zax[1]*zax[2]*(1-c)-zax[0]*s)*rot1[2];
	rot2[2]=(zax[2]*zax[0]*(1-c)-zax[1]*s)*rot1[0] + (zax[2]*zax[1]*(1-c)+zax[0]*s)*rot1[1] + (c+zax[2]*zax[2]*(1-c))       *rot1[2];

	p4cm[3]=_m/2.;
	p4cm[0]=p4cm[3]*rot2[0];
	p4cm[1]=p4cm[3]*rot2[1];
	p4cm[2]=p4cm[3]*rot2[2];
  
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
  p3[3]=pV[3]-p4[3];
  p3[0]=pV[0]-p4[0];
  p3[1]=pV[1]-p4[1];
  p3[2]=pV[2]-p4[2];
}

//CS FRAME PRESCRIPTION
double costhCS()
{
  return (1.-4.*(kap1[3]*p4[3]-kap1[2]*p4[2]-kap1[1]*p4[1]-kap1[0]*p4[0])/(_m*_m));
}

void setcosth(double costh) {_costh = costh;}
void setm(double mm) {_m = mm;}
void sety(double yy) {_y = yy;}
void sety_(double &yy) {_y = yy;}
void setqt(double qtt) {_qt = qtt;}
void setqt_(double &qtt) {_qt = qtt;}

double resy(double y, void* param)
{
  sety(y);

  double a = 66;
  double b = 116;
  double result = 0;
  double  err = 0;

  size_t limit = 1000;
  gsl_integration_workspace * gslw = gsl_integration_workspace_alloc (limit);
  gsl_function gslfm;
  gslfm.function = &resm;
  clock_t begin_time = clock();
  gsl_integration_qag(&gslfm, a, b, 0, 1E-2, limit, 1, gslw, &result, &err);
  //size_t neval;
  //gsl_integration_qng(&gslfm, a, b, 0, 1E-2, &result, &err, &neval);
  clock_t end_time = clock();
  cout << "Call to function resy, y: " << y 
  	 << " 66-116 mass integral: " << result 
  	 << " error: " << err 
  	 << " time: " << float( end_time - begin_time ) /  CLOCKS_PER_SEC <<endl;
  return result;
}

double resth(double costh, void* param)
{
  setcosth(costh);

  double a = 0;
  double b = 5;
  double result = 0;
  double  err = 0;

  size_t limit = 1000;
  gsl_integration_workspace * gslw = gsl_integration_workspace_alloc (limit);
  gsl_function gslfy;
  gslfy.function = &resy;
  clock_t begin_time = clock();
  gsl_integration_qag(&gslfy, a, b, 0, 1E-2, limit, 1, gslw, &result, &err);
  //size_t neval;
  //gsl_integration_qng(&gslfy, a, b, 0, 1E-2, &result, &err, &neval);
  clock_t end_time = clock();
  cout << "Call to function resth, costh: " << costh
         << " 0-5 rapidity integral: " << result 
         << " error: " << err
         << " time: " << float( end_time - begin_time ) /  CLOCKS_PER_SEC <<endl;
  return result*2;
}

double resm(double mm, void* param)
{
  int mode = 0.;
  double result = resumm_(_costh, mm, _qt, _y,mode);
  //  cout << "Call to function resm, mass: " << mm << " integral: " << result << endl;
  return result;
}

void costhbound(double phi_lep, vector<double> &min, vector<double> &max)
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
  if (cuts(p3, p4))
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
		  if (!cuts(p3, p4))
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
		  if (cuts(p3, p4))
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

//convert this trapezoidal rule to gaussian quadratures in fixed segments of phi
void cthmoments_(double &cthmom0, double &cthmom1, double &cthmom2)
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
	cout << "cuba point " << setw(10) << _m  << setw(10) << _qt  << setw(10) << _y
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
	cout << "trap point " << setw(10) << _m  << setw(10) << _qt  << setw(10) << _y
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
      end_time = clock();

      if (opts.verbose)
	cout << "quad point " << setw(10) << _m  << setw(10) << _qt  << setw(10) << _y
	     << setw(13) << cthmom0 << setw(13)  << 0
	     << setw(13) << cthmom1 << setw(13)  << 0
	     << setw(13) << cthmom2 << setw(13)  << 0
	     << "   time " << float( end_time - begin_time ) /  CLOCKS_PER_SEC << endl;
    }
}


integrand_t thphiintegrand(const int &ndim, const double x[], const int &ncomp, double f[])
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
  if (cuts(p3, p4))
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

integrand_t resintegrand2d(const int &ndim, const double x[], const int &ncomp, double f[])
{
  clock_t begin_time, end_time;

  begin_time = clock();

  //Jacobian of the change of variables from the unitary hypercube x[3] to the m, y, qt boundaries
  double jac = 1.;

  // Generate the boson invariant mass between the integration boundaries
  double wsqmin = pow(mmin,2);
  double wsqmax = pow(mmax,2);
  double x1=x[0];
  double m2,wt;
  breitw_(x1,wsqmin,wsqmax,opts.rmass,opts.rwidth,m2,wt);
  double m=sqrt(m2);
  jac=jac*wt;
  //Dynamic scale (not implemented yet)
  //if(dynamicscale) call scaleset(m2)

  //integrate up to the qT kinematical limit
  /*
  double dexpy=exp(y);
  double dexpmy=exp(-y);
  double cosh2y=pow(((dexpy+dexpmy)*0.5),2);
  double qtmax=sqrt(pow((pow(energy_.sroot_,2)+m2),2)/(4.*pow(energy_.sroot_,2)*cosh2y) - m2);
  double qtmin=0.1;
  double qt=qtmin+qtmax*x[2];
  jac=jac*(qtmax);
  */

  //integrate between qtmin and qtmax
  double qtmn = max(0.02, qtmin);
  double qt=qtmn+(qtmax-qtmn)*x[1];
  jac=jac*(qtmax-qtmn);

  //set global variables to m, qt
  set(m, qt, 0);

  //generate boson 4-momentum, with m, qt, y and phi=0
  genV4p(m, qt, 0, 0.);

  //Perform quadrature rule integration in rapidity and nested costh integration
  int nocuts = !opts.makelepcuts;

  //Limit y boundaries to the kinematic limit in y
  double ylim = 0.5*log(pow(energy_.sroot_,2)/m2);
  double ymn = min(max(-ylim, ymin),ylim);
  double ymx = max(min(ylim, ymax),-ylim);

  clock_t ybt, yet;
  ybt = clock();
  //  rapintegrals_(ymin,ymax,m,nocuts);
  rapintegrals_(ymn,ymx,m,nocuts);
  yet = clock();
  
  set(m, qt, 0);
  genV4p(m, qt, 0, 0.);

  //  SWITCHING FUNCTIONS
  double swtch=1.;
  if (qt >= m*3/4.)  swtch=exp(-pow((m*3/4.-qt),2)/pow((m/2.),2)); // GAUSS SWITCH
  if (swtch <= 0.01) swtch = 0;

  //Call to the resummation part
  double costh = 0;
  double y = 0.5*(ymx+ymn); //needed to get the correct IFIT index of the approximate PDF
  int mode = 2;

  clock_t rbt = clock();
  if (swtch < 0.01)
    f[0]=0.;
  else
    f[0]=resumm_(costh,m,qt,y,mode)/(8./3.);
  clock_t ret = clock();

  if (f[0] != f[0])
    f[0]=0.;  //avoid nans
	   
  f[0] = f[0]*jac*swtch;
  end_time = clock();
  if (opts.timeprofile)
    cout << setw (3) << "m" << setw(10) << m << setw(4) << "qt" << setw(10) <<  qt
	 << setw(8) << "result" << setw(10) << f[0]
	 << setw(10) << "tot time" << setw(10) << float( end_time - begin_time ) /  CLOCKS_PER_SEC
	 << setw(10) << "rapint"  << setw(10) << float( yet - ybt ) /  CLOCKS_PER_SEC
	 << setw(10) << "resumm"  << setw(10) << float( ret - rbt ) /  CLOCKS_PER_SEC
	 << endl;
  return 0;
}


integrand_t resintegrand3d(const int &ndim, const double x[], const int &ncomp, double f[])
{
  clock_t begin_time, end_time;

  begin_time = clock();

  //Jacobian of the change of variables from the unitary hypercube x[3] to the m, y, qt boundaries
  double jac = 1.;

  // Generate the boson invariant mass between the integration boundaries
  double wsqmin = pow(mmin,2);
  double wsqmax = pow(mmax,2);
  double x1=x[0];
  double m2,wt;
  breitw_(x1,wsqmin,wsqmax,opts.rmass,opts.rwidth,m2,wt);
  double m=sqrt(m2);
  jac=jac*wt;

  /*
  //split the mass in 3 pieces,
  //and unweight the breit wigner only from -5 width to +5 width
  //no sense, I have to split the mass before entering the integration
  //and set breit wigner unweighting or not
  double bwmn = opts.rmass -5*opts.rwidth;
  double bwmx = opts.rmass +5*opts.rwidth;
  double bwmnsq = pow(bwmn,2);
  double bwmxsq = pow(bwmx,2);
  double m,m2,wt;
  double xl = 0.25;
  double xu = 0.75;
  if (x[0] < xl)
    {
      double x1=x[0]/(xl-0.);
      m=mmin+(bwmn-mmin)*x[0];
      jac=jac*(bwmn-mmin)/(xl-0.);
    }
  else if (x[0] > xu)
    {
      double x1=(x[0]-xu)/(1.-xu);
      m=bwmx+(mmax-bwmx)*x1;
      jac=jac*(mmax-bwmx)/(1.-xu);
    }
  else
    {
      double xbw=(x[0]-xl)/(xu-xl);
      breitw_(xbw,bwmnsq,bwmxsq,opts.rmass,opts.rwidth,m2,wt);
      m=sqrt(m2);
      wt=wt/(xu-xl);
      jac=jac*wt;
    }
  */

  //  double m=mmin+(mmax-mmin)*x[0];
  //  jac=jac*(mmax-mmin);

  //Dynamic scale (not implemented yet)
  //if(dynamicscale) call scaleset(m2)

  //Limit y boundaries to the kinematic limit in y
  double ylim = 0.5*log(pow(energy_.sroot_,2)/m2);
  double ymn = min(max(-ylim, ymin),ylim);
  double ymx = max(min(ylim, ymax),-ylim);

  //integrate between ymin and ymax
  double y=ymn+(ymx-ymn)*x[1];
  jac=jac*(ymx-ymn);

  /*
  //Integrate the full phase space
  double y=-ylim+2.*ylim*x[1];
  jac=jac*2.*ylim;
  */


  //integrate up to the qT kinematical limit
  /*
  double dexpy=exp(y);
  double dexpmy=exp(-y);
  double cosh2y=pow(((dexpy+dexpmy)*0.5),2);
  double qtmax=sqrt(pow((pow(energy_.sroot_,2)+m2),2)/(4.*pow(energy_.sroot_,2)*cosh2y) - m2);
  double qtmin=0.1;
  double qt=qtmin+qtmax*x[2];
  jac=jac*(qtmax);
  */

  //integrate between qtmin and qtmax
  double qtmn = max(0.02, qtmin);
  double qt=qtmn+(qtmax-qtmn)*x[2];
  jac=jac*(qtmax-qtmn);

  //set global variables to m, qt, y
  set(m, qt, y);

  //generate boson 4-momentum, with m, qt, y and phi=0
  genV4p(m, qt, y, 0.);

  //  SWITCHING FUNCTIONS
  double swtch=1.;
  if (qt >= m*3/4.)  swtch=exp(-pow((m*3/4.-qt),2)/pow((m/2.),2)); // GAUSS SWITCH
  if (swtch <= 0.01) swtch = 0;

  //In this point of phase space (m, qt, y) the costh integration is performed by 
  //calculating the 0, 1 and 2 moments of costh
  //that are the integrals int(dcosth dphi1 dphi2), int(costh dcosth dphi1 dphi2), and int(costh^2 dcosth dphi1 dphi2) convoluted with cuts
  //Then the epxressions 1, costh and costh^2 in sigmaij are substituted by these costh moments
  double costh = 0;
  int mode = 1;
  if (swtch < 0.01)
    f[0]=0.;
  else
    //evaluate the resummed cross section
    f[0]=resumm_(costh,m,qt,y,mode)/(8./3.);

  if (f[0] != f[0])
    f[0]=0.;  //avoid nans
	   
  f[0] = f[0]*jac*swtch;

  end_time = clock();
  if (opts.timeprofile)
    cout << setw (3) << "m" << setw(10) << m << setw(4) << "qt" << setw(10) <<  qt
	 << setw(8) << "result" << setw(10) << f[0]
	 << setw(10) << "tot time" << setw(10) << float( end_time - begin_time ) /  CLOCKS_PER_SEC
	 << endl;
  return 0;
}

//Perform the integration as in the original dyres code
integrand_t resintegrand4d(const int &ndim, const double x[], const int &ncomp, double f[],
			   void* userdata, const int &nvec, const int &core,
			   double &weight, const int &iter)
{
  clock_t begin_time, end_time;

  double r[ndim];
  for (int i = 0; i < ndim; i++)
    r[i]=x[i];

  double resumm,tempp,ran2;
  int ii,j,ij;
  double m,m2,qt2,s,wgt;
  double p[5][5],pV[4],p4cm[4],ptemp[4],kap1[5],kap1b[5];
  double pjet[4][12];
  double phi,cos_th,phi_lep;
  int n2,n3;
  double msq,wt,swtch,kk,zeta1,zeta1b,kt1,kt2,mt2,costh_CS,qP1,qP2,ritmx;

  double wsqmin = pow(mmin,2);
  double wsqmax = pow(mmax,2);
  
  double lowintHst=0;
  double lowintHst0=0;

  int jstop = 0;

  double azloopmax = 500;
  double azloop = 0;

  ptemp[0]=0.;
  ptemp[1]=0.;
  ptemp[2]=0.;
  ptemp[3]=0.;
  for (int i=0; i < 5; i++)
    {
      p[i][1]=0.;
      p[i][2]=0.;
      p[i][3]=0.;
      p[i][4]=0.;
    }

  // Generate the boson invariant mass

  double jac=1.;
  double x1=r[0];
  breitw_(x1,wsqmin,wsqmax,opts.rmass,opts.rwidth,msq,wt);
  m2=msq;
  m=sqrt(m2);
  jac=jac*wt;
  //     Dynamic scale
  //      if(dynamicscale) call scaleset(m2)


  //Limit y boundaries to the kinematic limit in y
  double ylim = 0.5*log(pow(energy_.sroot_,2)/m2);
  double ymn = min(max(-ylim, ymin),ylim);
  double ymx = max(min(ylim, ymax),-ylim);

  //integrate between ymin and ymax
  double y=ymn+(ymx-ymn)*r[1];
  jac=jac*(ymx-ymn);

  double dexpy=exp(y);
  double dexpmy=exp(-y);
  double cosh2y=pow(((dexpy+dexpmy)*0.5),2);

  //   qT kinematical limit
  //  qtmax=sqrt(pow((pow(energy_.sroot_,2)+m2),2)/(4.*pow(energy_.sroot_,2)*cosh2y) - m2);
  //  qtmin=0.1;
  //  qt=qtmin+qtmax*r[2];
  //  jac=jac*(qtmax);

  double qtmn = max(0.02, qtmin);
  double qt=qtmn+(qtmax-qtmn)*r[2];
  jac=jac*(qtmax-qtmn);

  qt2=pow(qt,2);

  double xx1=sqrt(m2/pow(energy_.sroot_,2))*dexpy;
  double xx2=sqrt(m2/pow(energy_.sroot_,2))*dexpmy;


  // incoming quarks
  p[1][1]=0.;
  p[1][2]=0.;
  p[1][3]=xx1*0.5*energy_.sroot_;
  p[1][4]=xx1*0.5*energy_.sroot_;
  p[2][1]=0.;
  p[2][2]=0.;
  p[2][3]=-xx2*0.5*energy_.sroot_;
  p[2][4]=xx2*0.5*energy_.sroot_;


  // First lepton direction: Cos of the polar angle
  cos_th=-1.+2.*r[3];
  jac=jac*2.;
  
  mt2=m2+qt2;

  int mode = 0;
  // LOOP over (vector boson and lepton) azimuthal angles 
  for (int j=0; j < azloopmax; j++)
    // Vector boson azimuthal angle
    {
      phi = 2.*M_PI*((double)rand()/(double)RAND_MAX);

      // Lepton decay in the center of mass frame 
      // First lepton direction:  azimuthal angle
      phi_lep = 2.*M_PI*((double)rand()/(double)RAND_MAX);

      //  vector boson momentum: pV(4)^2-pV(1)^2-pV(2)^2-pV(3)^2=m2
      pV[0]=qt*cos(phi);
      pV[1]=qt*sin(phi);
      pV[2]=0.5*sqrt(mt2)*(exp(y)-exp(-y));
      pV[3]=0.5*sqrt(mt2)*(exp(y)+exp(-y));

      // momentum of the first lepton 
      p4cm[3]=m/2.;
      p4cm[0]=p4cm[3]*sin(acos(cos_th))*sin(phi_lep);
      p4cm[1]=p4cm[3]*sin(acos(cos_th))*cos(phi_lep);
      p4cm[2]=p4cm[3]*cos_th;

      // Boost to go in the Lab frame
      boost_(m,pV,p4cm,ptemp);
      p[4][1]=ptemp[0];
      p[4][2]=ptemp[1];
      p[4][3]=ptemp[2];
      p[4][4]=ptemp[3];

      //  momentum of the second lepton
      p[3][4]=pV[3]-p[4][4];
      p[3][1]=pV[0]-p[4][1];
      p[3][2]=pV[1]-p[4][2];
      p[3][3]=pV[2]-p[4][3];

      //CS frame prescription
      if (opts.qtrec_cs)
	{
	  kt1 = pV[0]/2.;
	  kt2 = pV[1]/2.;
	}

      //MY (DYRES) prescription
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

      zeta1=1./m2/2.*(m2+2.*(pV[0]*kt1+pV[1]*kt2)+sqrt(pow((m2+2.*(pV[0]*kt1+pV[1]*kt2)),2)-4.*mt2*(pow(kt1,2)+pow(kt2,2))));

      qP1=(pV[3]-pV[2])*energy_.sroot_/2.;
      qP2=(pV[3]+pV[2])*energy_.sroot_/2.;

      kap1[4]=energy_.sroot_/2.*(zeta1*m2/2./qP1+(pow(kt1,2)+pow(kt2,2))/zeta1*qP1/m2/pow(energy_.sroot_,2)*2.);
      kap1[1]=kt1;
      kap1[2]=kt2;
      kap1[3]=energy_.sroot_/2.*(zeta1*m2/2./qP1-(pow(kt1,2)+pow(kt2,2))/zeta1*qP1/m2/pow(energy_.sroot_,2)*2.);

      costh_CS=1.-4.*(kap1[4]*p[4][4]-kap1[3]*p[4][3]-kap1[2]*p[4][2]-kap1[1]*p[4][1])/m2;

       // see whether this point will pass cuts - if it will not, do not
       // bother calculating the matrix elements for it, instead bail out
      for (int i=0; i < 4; i++)
	{
	  pjet[i][0] = p[1][i+1];
	  pjet[i][1] = p[2][i+1];
	  pjet[i][2] = p[3][i+1];
	  pjet[i][3] = p[4][i+1];
	}
      int njet = 0;
       //       cout << qt << "  " <<  y << "  " <<  m << "  " << sqrt(p[3][1]*p[3][1] + p[3][2] * p[3][2]) << endl;
      if (opts.makelepcuts)
	if (cuts_(pjet,njet)) continue;
			    
      //  SWITCHING FUNCTIONS
      swtch=1.;
      if (qt >= m*3/4.)  swtch=exp(-pow((m*3/4.-qt),2)/pow((m/2.),2)); // GAUSS SWITCH

      if (swtch <= 0.01) break;

      if (jstop != 1)
	{
	  jstop=1;

	  //Call to the resummation part
	  if (swtch < 0.01)
	    tempp=0.;
	  else
	    tempp=resumm_(costh_CS,m,qt,y,mode)/(8./3.);  // CS PRESCRIPTION
	    //tempp=resumm_(cos_th,m,qt,y,mode)/(8./3.);
	   if (tempp != tempp)
	     tempp=0.;    //  avoid nans
	   
	   lowintHst0=jac*tempp;
	   lowintHst0=lowintHst0*swtch; // SWITCHING
	 }
      hists.FillEvent(p[3],p[4],weight*lowintHst0/azloopmax);
      azloop=azloop+1;
    }

  lowintHst=lowintHst0*float(azloop)/float(azloopmax);
  
  f[0] = lowintHst;

  return 0;
}
