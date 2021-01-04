#include "icoeff.h"
#include "ifunc.h"
#include "constants.h"
#include "ccoeff.h"

#include "anomalous.h"
#include "mellinint.h"
#include "mesq.h"
#include "settings.h"
#include "gaussrules.h"
#include "resconst.h"
#include "psi.h"
#include "clock_real.h"
//#include <ginac/ginac.h>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace constants;
using namespace resconst;

double icoeff::tlog[rule];
double icoeff::tlin[rule];
double icoeff::faclog[rule];
double icoeff::faclin[rule];

complex <double> *icoeff::kernlog;
complex <double> *icoeff::kernlin;
complex <double> *icoeff::kernlog_1;
complex <double> *icoeff::kernlin_1;
complex <double> *icoeff::kernlog_2;
complex <double> *icoeff::kernlin_2;
  
double icoeff::c3qqpzlog[rule];
double icoeff::c3qqbzlog[rule];
double icoeff::c3qqbpzlog[rule];
double icoeff::c3qgzlin[rule];
double icoeff::c3qqzlin[rule];

double icoeff::C1qq_delta;
double icoeff::C2qq_delta;
double icoeff::C3qq_delta;

double icoeff::C2qq_plus;
double icoeff::C3qq_plus;

complex <double> *icoeff::c1qg;
complex <double> *icoeff::c1qq;
complex <double> *icoeff::c1qqb;
complex <double> *icoeff::c1qqp;
complex <double> *icoeff::c1qqbp;

complex <double> *icoeff::c2qg;
complex <double> *icoeff::c2qq;
complex <double> *icoeff::c2qqb;
complex <double> *icoeff::c2qqp;
complex <double> *icoeff::c2qqbp;

complex <double> *icoeff::c3qg;
complex <double> *icoeff::c3qq;
complex <double> *icoeff::c3qqb;
complex <double> *icoeff::c3qqp;
complex <double> *icoeff::c3qqbp;

complex <double> *icoeff::c3qg_1;
complex <double> *icoeff::c3qq_1;
complex <double> *icoeff::c3qqb_1;
complex <double> *icoeff::c3qqp_1;
complex <double> *icoeff::c3qqbp_1;
complex <double> *icoeff::c3qg_2;
complex <double> *icoeff::c3qq_2;
complex <double> *icoeff::c3qqb_2;
complex <double> *icoeff::c3qqp_2;
complex <double> *icoeff::c3qqbp_2;

//C functions in z space
double C1qg(double z)   {return ifunc::I1qg(z);};
double C1qq(double z)   {return ifunc::I1qq(z);};

double C2qg(double z)   {return ifunc::I2qg(z) + icoeff::C1qq_delta*C1qg(z);};
double C2qq(double z)   {return ifunc::I2qq(z) + icoeff::C1qq_delta*C1qq(z);};
double C2qqb(double z)  {return ifunc::I2qqb(z);};
double C2qqp(double z)  {return ifunc::I2qqp(z);};
double C2qqbp(double z) {return ifunc::I2qqp(z);};

//double C3qg(double z)   {return ifunc::I3qg(z)   + icoeff::C1qq_delta*C2qg(z) + (icoeff::C2qq_delta - pow(icoeff::C1qq_delta,2))*C1qg(z);};
//double C3qq(double z)   {return ifunc::I3qq(z)   + icoeff::C1qq_delta*C2qq(z) + (icoeff::C2qq_delta - pow(icoeff::C1qq_delta,2))*C1qq(z);};
//double C3qqp(double z)  {return ifunc::I3qqp(z)  + icoeff::C1qq_delta*C2qqp(z);};
//double C3qqb(double z)  {return ifunc::I3qqb(z)  + icoeff::C1qq_delta*C2qqb(z);};
//double C3qqbp(double z) {return ifunc::I3qqbp(z) + icoeff::C1qq_delta*C2qqbp(z);};

//Use cached values
double C3qg(double z)   {return ifunc::i3qg   + icoeff::C1qq_delta*C2qg(z) + (icoeff::C2qq_delta - pow(icoeff::C1qq_delta,2))*C1qg(z);};
double C3qq(double z)   {return ifunc::i3qq   + icoeff::C1qq_delta*C2qq(z) + (icoeff::C2qq_delta - pow(icoeff::C1qq_delta,2))*C1qq(z);};
double C3qqp(double z)  {return ifunc::i3qqp  + icoeff::C1qq_delta*C2qqp(z);};
double C3qqb(double z)  {return ifunc::i3qqb  + icoeff::C1qq_delta*C2qqb(z);};
double C3qqbp(double z) {return ifunc::i3qqbp + icoeff::C1qq_delta*C2qqbp(z);};

void icoeff::allocate()
{
  if (opts.mellin1d)
    {
      c1qg = new complex <double>[mellinint::mdim*2];
      c1qq = new complex <double>[mellinint::mdim*2];
      c1qqb = new complex <double>[mellinint::mdim*2];
      c1qqp = new complex <double>[mellinint::mdim*2];
      c1qqbp = new complex <double>[mellinint::mdim*2];

      c2qg = new complex <double>[mellinint::mdim*2];
      c2qq = new complex <double>[mellinint::mdim*2];
      c2qqb = new complex <double>[mellinint::mdim*2];
      c2qqp = new complex <double>[mellinint::mdim*2];
      c2qqbp = new complex <double>[mellinint::mdim*2];

      c3qg = new complex <double>[mellinint::mdim*2];
      c3qq = new complex <double>[mellinint::mdim*2];
      c3qqb = new complex <double>[mellinint::mdim*2];
      c3qqp = new complex <double>[mellinint::mdim*2];
      c3qqbp = new complex <double>[mellinint::mdim*2];

      kernlog = new complex <double> [rule*mellinint::mdim*2];
      kernlin = new complex <double> [rule*mellinint::mdim*2];
    }
  else
    {
      c3qg_1 = new complex <double>[mellinint::mdim*2];
      c3qq_1 = new complex <double>[mellinint::mdim*2];
      c3qqb_1 = new complex <double>[mellinint::mdim*2];
      c3qqp_1 = new complex <double>[mellinint::mdim*2];
      c3qqbp_1 = new complex <double>[mellinint::mdim*2];

      c3qg_2 = new complex <double>[mellinint::mdim*2];
      c3qq_2 = new complex <double>[mellinint::mdim*2];
      c3qqb_2 = new complex <double>[mellinint::mdim*2];
      c3qqp_2 = new complex <double>[mellinint::mdim*2];
      c3qqbp_2 = new complex <double>[mellinint::mdim*2];
      
      kernlog_1 = new complex <double> [rule*mellinint::mdim*2];
      kernlin_1 = new complex <double> [rule*mellinint::mdim*2];

      kernlog_2 = new complex <double> [rule*mellinint::mdim*2];
      kernlin_2 = new complex <double> [rule*mellinint::mdim*2];
    }
}

//Need to be called after ccoeff::delta()
void icoeff::delta()
{
  //Convert from as/pi to as/4/pi Normalisation

  //Delta coefficients
  C1qq_delta = ccoeff::C1qq_delta *4.;
  C2qq_delta = ccoeff::C2qq_delta *16.;
  C3qq_delta = ccoeff::C3qq_delta *64.;
  
  //Plus-distribution coefficients
  C2qq_plus = ((28*(zeta3)-(808./27.))*(CA)*(CF)+(224./27.)*(NF)*(CF)*(TF));
  C3qq_plus = ((CA)*(CF)*(NF)*(TF)*(-((1648*(zeta2))/81.)-((1808*(zeta3))/27.)+((40*(zeta4))/3.)+(125252./729.))+(pow(CA,2))*(CF)*(-(176./3.)*(zeta3)*(zeta2)+((6392*(zeta2))/81.)+((12328*(zeta3))/27.)+((154*(zeta4))/3.)-192*(zeta5)-(297029./729.))+(CF)*(pow(NF,2))*(pow(TF,2))*(-((128*(zeta3))/9.)-(7424./729.))+(pow(CF,2))*(NF)*(TF)*(-((608*(zeta3))/9.)-32*(zeta4)+(3422./27.)));
}

void icoeff::init()
{
  // boundaries of integration
  double xmax = 1;
  //double xminlog = pow(bins.mbins.front()/opts.sroot,2);
  //double xminlin = pow(bins.mbins.front()/opts.sroot,2);
  double xminlog = 1e-14;
  double xminlin = 0;
  double ll = log(xmax/xminlog);
  double cc = 0.5;
  double mm = 0.5;
  for (int i = 0; i < rule; i++)
    {
      double x = cc+mm*gr::xxx[rule-1][i];
      tlog[i] = xminlog*exp(ll*x);
      tlin[i] = xminlin+(xmax-xminlin)*x;
      double jaclog = mm * tlog[i] * ll;
      double jaclin = mm*(xmax-xminlin);
      faclog[i] = jaclog * gr::www[rule-1][i];
      faclin[i] = jaclin * gr::www[rule-1][i];
    }

  //Write out
  fstream testfile(string(SHAREDIR)+"/i3.bin", ios::in | ios::binary);
  if (testfile)
    testfile.close(); //If the i3.bin file exists, skip the expensive C3 calculation
  else
    {
      //GiNaC::Digits = 12;
      clock_t begin_time, end_time;
      begin_time = clock();  
      cout << "First time initialisation of C3 coefficients, please wait... " << flush;
      for (int i = 0; i < rule; i++)
	{
	  double zlog = tlog[i];
	  double zlin = tlin[i];

	  //begin_time = clock();
	  ifunc::I3(zlog);
	  c3qqpzlog[i]  = C3qqp(zlog);
	  c3qqbzlog[i]  = C3qqb(zlog);
	  c3qqbpzlog[i] = C3qqbp(zlog);

	  ifunc::I3(zlin);
	  c3qgzlin[i]   = C3qg(zlin);
	  c3qqzlin[i]   = C3qq(zlin);
	  //end_time = clock();
	  //cout << "z = " << z << " C Coefficients evaluated in "  << float(end_time - begin_time)/CLOCKS_PER_SEC*1000 << " ms" << endl;
	}

      fstream off(string(SHAREDIR)+"/i3.bin", ios::out | ios::binary);
      off.write((char*)&c3qgzlin,rule*sizeof(double));
      off.write((char*)&c3qqzlin,rule*sizeof(double));
      off.write((char*)&c3qqbzlog,rule*sizeof(double));
      off.write((char*)&c3qqpzlog,rule*sizeof(double));
      off.write((char*)&c3qqbpzlog,rule*sizeof(double));
      off.close();

      end_time = clock();
      cout << " C Coefficients evaluated in "  << float(end_time - begin_time)/CLOCKS_PER_SEC << " s" << endl;
      cout << endl;
    }

  //Read in
  fstream iff(string(SHAREDIR)+"/i3.bin", ios::in | ios::binary);
  if (iff)
    {
      iff.read((char*)&c3qgzlin,rule*sizeof(double));
      iff.read((char*)&c3qqzlin,rule*sizeof(double));
      iff.read((char*)&c3qqbzlog,rule*sizeof(double));
      iff.read((char*)&c3qqpzlog,rule*sizeof(double));
      iff.read((char*)&c3qqbpzlog,rule*sizeof(double));
      iff.close();
      return;
    }
}

void icoeff::calc1d()
{
  for (int i = 0; i < rule; i++)
    for (int m = 0; m < mellinint::mdim; m++)
      {
	//cout << i << "  " << m << "  " << i*mellinint::mdim+m << endl;
	kernlog[i*mellinint::mdim+m] = pow(tlog[i], mellinint::Np[m]-1.);
	kernlin[i*mellinint::mdim+m] = pow(tlin[i], mellinint::Np[m]-1.);
      }

  fill(c1qg,   c1qg+2*mellinint::mdim, 0.);
  fill(c1qq,   c1qq+2*mellinint::mdim, 0.);
  fill(c1qqb,  c1qqb+2*mellinint::mdim, 0.);
  fill(c1qqp,  c1qqp+2*mellinint::mdim, 0.);
  fill(c1qqbp, c1qqbp+2*mellinint::mdim, 0.);

  fill(c2qg,   c2qg+2*mellinint::mdim, 0.);
  fill(c2qq,   c2qq+2*mellinint::mdim, 0.);
  fill(c2qqb,  c2qqb+2*mellinint::mdim, 0.);
  fill(c2qqp,  c2qqp+2*mellinint::mdim, 0.);
  fill(c2qqbp, c2qqbp+2*mellinint::mdim, 0.);

  fill(c3qg,   c3qg+2*mellinint::mdim, 0.);
  fill(c3qq,   c3qq+2*mellinint::mdim, 0.);
  fill(c3qqb,  c3qqb+2*mellinint::mdim, 0.);
  fill(c3qqp,  c3qqp+2*mellinint::mdim, 0.);
  fill(c3qqbp, c3qqbp+2*mellinint::mdim, 0.);

  //integral_0^1{ x^(N-1) fx dx}
  for (int i = 0; i < rule; i++)
    for (int m = 0; m < mellinint::mdim; m++)
      {
	int idx = anomalous::index(m,mesq::positive);
	/*
	c1qg[idx]   += faclin[i]*kernlin[i*mellinint::mdim+m] * C1qg(tlin[i]);
	c1qq[idx]   += faclin[i]*kernlin[i*mellinint::mdim+m] * C1qq(tlin[i]);
	
	c2qg[idx]   += faclin[i]*kernlin[i*mellinint::mdim+m] * C2qg(tlin[i]);
	c2qq[idx]   += faclin[i]*kernlin[i*mellinint::mdim+m] * C2qq(tlin[i]);
	c2qqb[idx]  += faclog[i]*kernlog[i*mellinint::mdim+m] * C2qqb(tlog[i]);
	c2qqp[idx]  += faclog[i]*kernlog[i*mellinint::mdim+m] * C2qqp(tlog[i]);
	c2qqbp[idx] += faclog[i]*kernlog[i*mellinint::mdim+m] * C2qqbp(tlog[i]); // (= C2qqp(z))
	*/
	  
	c3qg[idx]   += faclin[i]*kernlin[i*mellinint::mdim+m] * c3qgzlin[i];
	c3qq[idx]   += faclin[i]*kernlin[i*mellinint::mdim+m] * c3qqzlin[i];
	c3qqp[idx]  += faclog[i]*kernlog[i*mellinint::mdim+m] * c3qqpzlog[i];
	c3qqb[idx]  += faclog[i]*kernlog[i*mellinint::mdim+m] * c3qqbzlog[i];
	c3qqbp[idx] += faclog[i]*kernlog[i*mellinint::mdim+m] * c3qqbpzlog[i];
	
	//cout << i << "  " << kern[i*mellinint::mdim+m] << endl;
      }
 
  //Add delta pieces
  for (int m = 0; m < mellinint::mdim; m++)
    {
      int idx = anomalous::index(m,mesq::positive);
      c1qq[idx] += C1qq_delta;
      c2qq[idx] += C2qq_delta;
      c3qq[idx] += C3qq_delta;
    }

  //Add plus-distribution pieces
  for (int m = 0; m < mellinint::mdim; m++)
    {
      //fcomplex xn = fcx(mellinint::Np[m]);
      //fcomplex ps0n;
      //psi0_(xn,ps0n);
      //complex <double> s1nm1 = cx(ps0n)+euler;

      //The Mellin transform of the plus distribution is Mel[1/(1-x)+] = -S1(N-1) = -(psi(N) + gE)
      complex <double> s1nm1 = cpsi0(mellinint::Np[m])+euler;
      int idx = anomalous::index(m,mesq::positive);
      c2qq[idx] += (-s1nm1)*C2qq_plus;
      c3qq[idx] += (-s1nm1)*(C3qq_plus+C1qq_delta*C2qq_plus);
    }

  ////Normalise from as/4/pi to as/2/pi (as C1 and C2 coefficients in DYRes)
  //Normalise from as/4/pi to as/pi
  for (int m = 0; m < mellinint::mdim; m++)
    {
      int idx = anomalous::index(m,mesq::positive);

      c1qg[idx]   /= 4.; //2.;
      c1qq[idx]   /= 4.; //2.;
	  
      c2qg[idx]   /= 16.; //4.;
      c2qq[idx]   /= 16.; //4.;
      c2qqb[idx]  /= 16.; //4.;
      c2qqp[idx]  /= 16.; //4.;
      c2qqbp[idx] /= 16.; //4.;
	  
      c3qg[idx]   /= 64.; //8.;
      c3qq[idx]   /= 64.; //8.;
      c3qqp[idx]  /= 64.; //8.;
      c3qqb[idx]  /= 64.; //8.;
      c3qqbp[idx] /= 64.; //8.;
    }
  
  //Compute negative branch
  for (int m = 0; m < mellinint::mdim; m++)
   {
     int idxp = anomalous::index(m,mesq::positive);
     int idxm = anomalous::index(m,mesq::negative);

     c1qg[idxm]   = conj(c1qg[idxp]);
     c1qq[idxm]   = conj(c1qq[idxp]);

     c2qg[idxm]   = conj(c2qg[idxp]);
     c2qq[idxm]   = conj(c2qq[idxp]);
     c2qqb[idxm]  = conj(c2qqb[idxp]);
     c2qqp[idxm]  = conj(c2qqp[idxp]);
     c2qqbp[idxm] = conj(c2qqbp[idxp]);

     c3qg[idxm]   = conj(c3qg[idxp]);
     c3qq[idxm]   = conj(c3qq[idxp]);
     c3qqb[idxm]  = conj(c3qqb[idxp]);
     c3qqp[idxm]  = conj(c3qqp[idxp]);
     c3qqbp[idxm] = conj(c3qqbp[idxp]);
   }

  /*
  //Check
  for (int m = 0; m < mellinint::mdim; m++)
    {
      int idx = anomalous::index(m,mesq::positive);

      cout << setprecision(16);
      cout << "C1QG " << m << " N " << mellinint::Np[m] << " anomalous " << anomalous::C1QG[idx] << " icoeff " << c1qg[idx] << endl;
      cout << "C1QQ " << m << " N " << mellinint::Np[m] << " anomalous " << anomalous::C1QQ[idx] << " icoeff " << c1qq[idx] << endl;

      cout << "C2QG  "  << m << " N " << mellinint::Np[m] << " anomalous " << anomalous::C2qgM[idx]                              << " icoeff " << c2qg[idx] << endl;
      cout << "C2QQ  "  << m << " N " << mellinint::Np[m] << " anomalous " << (anomalous::C2NSqqM[idx] +anomalous::C2SqqbM[idx]) << " icoeff " << c2qq[idx] << endl;
      cout << "C2QQB "  << m << " N " << mellinint::Np[m] << " anomalous " << (anomalous::C2NSqqbM[idx]+anomalous::C2SqqbM[idx]) << " icoeff " << c2qqb[idx] << endl;
      cout << "C2QQP "  << m << " N " << mellinint::Np[m] << " anomalous " << anomalous::C2SqqbM[idx]                            << " icoeff " << c2qqp[idx] << endl;
      cout << "C2QQBP " << m << " N " << mellinint::Np[m] << " anomalous " << anomalous::C2SqqbM[idx]                            << " icoeff " << c2qqbp[idx] << endl;

      cout << "C3QG  "  << m << " N " << mellinint::Np[m] << " icoeff " << c3qg[idx] << endl;
      cout << "C3QQ  "  << m << " N " << mellinint::Np[m] << " icoeff " << c3qq[idx] << endl;
      cout << "C3QQB "  << m << " N " << mellinint::Np[m] << " icoeff " << c3qqb[idx] << endl;
      cout << "C3QQP "  << m << " N " << mellinint::Np[m] << " icoeff " << c3qqp[idx] << endl;
      cout << "C3QQBP " << m << " N " << mellinint::Np[m] << " icoeff " << c3qqbp[idx] << endl;
      cout << endl << endl;


//      cout << setprecision(16);
//      cout << "C1QG " << m << " N " << mellinint::Np[m] << " ccoeff " << ccoeff::C1qg[idx] << " icoeff " << c1qg[idx] << endl;
//      cout << "C1QQ " << m << " N " << mellinint::Np[m] << " ccoeff " << ccoeff::C1qq[idx] << " icoeff " << c1qq[idx] << endl;
//
//      cout << "C2QG  "  << m << " N " << mellinint::Np[m] << " ccoeff " << ccoeff::C2qg[idx]   << " icoeff " << c2qg[idx] << endl;
//      cout << "C2QQ  "  << m << " N " << mellinint::Np[m] << " ccoeff " << ccoeff::C2qq[idx]   << " icoeff " << c2qq[idx] << endl;
//      cout << "C2QQB "  << m << " N " << mellinint::Np[m] << " ccoeff " << ccoeff::C2qqb[idx]  << " icoeff " << c2qqb[idx] << endl;
//      cout << "C2QQP "  << m << " N " << mellinint::Np[m] << " ccoeff " << ccoeff::C2qqp[idx]  << " icoeff " << c2qqp[idx] << endl;
//      cout << "C2QQBP " << m << " N " << mellinint::Np[m] << " ccoeff " << ccoeff::C2qqbp[idx] << " icoeff " << c2qqbp[idx] << endl;
    }
  */
}

void icoeff::calc2d()
{
  for (int i = 0; i < rule; i++)
    for (int m = 0; m < mellinint::mdim; m++)
      {
	kernlog_1[i*mellinint::mdim+m] = pow(tlog[i], mellinint::Np_1[m]-1.);
	kernlin_1[i*mellinint::mdim+m] = pow(tlin[i], mellinint::Np_1[m]-1.);

	kernlog_2[i*mellinint::mdim+m] = pow(tlog[i], mellinint::Np_2[m]-1.);
	kernlin_2[i*mellinint::mdim+m] = pow(tlin[i], mellinint::Np_2[m]-1.);
      }
  
  fill(c3qg_1,   c3qg_1+2*mellinint::mdim, 0.);
  fill(c3qq_1,   c3qq_1+2*mellinint::mdim, C3qq_delta); //Add delta piece
  fill(c3qqb_1,  c3qqb_1+2*mellinint::mdim, 0.);
  fill(c3qqp_1,  c3qqp_1+2*mellinint::mdim, 0.);
  fill(c3qqbp_1, c3qqbp_1+2*mellinint::mdim, 0.);

  fill(c3qg_2,   c3qg_2+2*mellinint::mdim, 0.);
  fill(c3qq_2,   c3qq_2+2*mellinint::mdim, C3qq_delta); //Add delta piece
  fill(c3qqb_2,  c3qqb_2+2*mellinint::mdim, 0.);
  fill(c3qqp_2,  c3qqp_2+2*mellinint::mdim, 0.);
  fill(c3qqbp_2, c3qqbp_2+2*mellinint::mdim, 0.);

  for (int i = 0; i < rule; i++)
    for (int m = 0; m < mellinint::mdim; m++)
      {
	int idx = anomalous::index(m,mesq::positive);
	c3qg_1[idx]   += faclin[i]*kernlin_1[i*mellinint::mdim+m] * c3qgzlin[i];
	c3qq_1[idx]   += faclin[i]*kernlin_1[i*mellinint::mdim+m] * c3qqzlin[i];
	c3qqp_1[idx]  += faclog[i]*kernlog_1[i*mellinint::mdim+m] * c3qqpzlog[i];
	c3qqb_1[idx]  += faclog[i]*kernlog_1[i*mellinint::mdim+m] * c3qqbzlog[i];
	c3qqbp_1[idx] += faclog[i]*kernlog_1[i*mellinint::mdim+m] * c3qqbpzlog[i];

	c3qg_2[idx]   += faclin[i]*kernlin_2[i*mellinint::mdim+m] * c3qgzlin[i];
	c3qq_2[idx]   += faclin[i]*kernlin_2[i*mellinint::mdim+m] * c3qqzlin[i];
	c3qqp_2[idx]  += faclog[i]*kernlog_2[i*mellinint::mdim+m] * c3qqpzlog[i];
	c3qqb_2[idx]  += faclog[i]*kernlog_2[i*mellinint::mdim+m] * c3qqbzlog[i];
	c3qqbp_2[idx] += faclog[i]*kernlog_2[i*mellinint::mdim+m] * c3qqbpzlog[i];
      }
 
  //Add plus-distribution pieces
  for (int m = 0; m < mellinint::mdim; m++)
    {
      int idx = anomalous::index(m,mesq::positive);

      complex <double> s1nm1_1 = cpsi0(mellinint::Np_1[m])+euler;
      complex <double> s1nm1_2 = cpsi0(mellinint::Np_2[m])+euler;
      c3qq_1[idx] += (-s1nm1_1)*(C3qq_plus+C1qq_delta*C2qq_plus);
      c3qq_2[idx] += (-s1nm1_2)*(C3qq_plus+C1qq_delta*C2qq_plus);
    }

  //Normalise from as/4/pi to as/pi
  for (int m = 0; m < mellinint::mdim; m++)
    {
      int idx = anomalous::index(m,mesq::positive);
      c3qg_1[idx]   /= 64.;
      c3qq_1[idx]   /= 64.;
      c3qqp_1[idx]  /= 64.;
      c3qqb_1[idx]  /= 64.;
      c3qqbp_1[idx] /= 64.;

      c3qg_2[idx]   /= 64.;
      c3qq_2[idx]   /= 64.;
      c3qqp_2[idx]  /= 64.;
      c3qqb_2[idx]  /= 64.;
      c3qqbp_2[idx] /= 64.;
    }
  
  //Compute negative branch
  for (int m = 0; m < mellinint::mdim; m++)
   {
     int idxp = anomalous::index(m,mesq::positive);
     int idxm = anomalous::index(m,mesq::negative);
     c3qg_1[idxm]   = conj(c3qg_1[idxp]);
     c3qq_1[idxm]   = conj(c3qq_1[idxp]);
     c3qqb_1[idxm]  = conj(c3qqb_1[idxp]);
     c3qqp_1[idxm]  = conj(c3qqp_1[idxp]);
     c3qqbp_1[idxm] = conj(c3qqbp_1[idxp]);

     c3qg_2[idxm]   = conj(c3qg_2[idxp]);
     c3qq_2[idxm]   = conj(c3qq_2[idxp]);
     c3qqb_2[idxm]  = conj(c3qqb_2[idxp]);
     c3qqp_2[idxm]  = conj(c3qqp_2[idxp]);
     c3qqbp_2[idxm] = conj(c3qqbp_2[idxp]);
   }

  /*
  //Check
  for (int m = 0; m < mellinint::mdim; m++)
    {
      int idx = anomalous::index(m,mesq::positive);

      cout << setprecision(16);
      cout << "C3QG  "  << m << " N " << mellinint::Np_1[m] << " icoeff " << c3qg_1[idx] << endl;
      cout << "C3QQ  "  << m << " N " << mellinint::Np_1[m] << " icoeff " << c3qq_1[idx] << endl;
      cout << "C3QQB "  << m << " N " << mellinint::Np_1[m] << " icoeff " << c3qqb_1[idx] << endl;
      cout << "C3QQP "  << m << " N " << mellinint::Np_1[m] << " icoeff " << c3qqp_1[idx] << endl;
      cout << "C3QQBP " << m << " N " << mellinint::Np_1[m] << " icoeff " << c3qqbp_1[idx] << endl;
      cout << "C3QG  "  << m << " N " << mellinint::Np_2[m] << " icoeff " << c3qg_2[idx] << endl;
      cout << "C3QQ  "  << m << " N " << mellinint::Np_2[m] << " icoeff " << c3qq_2[idx] << endl;
      cout << "C3QQB "  << m << " N " << mellinint::Np_2[m] << " icoeff " << c3qqb_2[idx] << endl;
      cout << "C3QQP "  << m << " N " << mellinint::Np_2[m] << " icoeff " << c3qqp_2[idx] << endl;
      cout << "C3QQBP " << m << " N " << mellinint::Np_2[m] << " icoeff " << c3qqbp_2[idx] << endl;
      cout << endl << endl;
    }
  */
}

void icoeff::release()
{
  if (opts.mellin1d)
    {
      delete[] c1qg;
      delete[] c1qq;
      delete[] c1qqb;
      delete[] c1qqp;
      delete[] c1qqbp;
      delete[] c2qg;
      delete[] c2qq;
      delete[] c2qqb;
      delete[] c2qqp;
      delete[] c2qqbp;
      delete[] c3qg;
      delete[] c3qq;
      delete[] c3qqb;
      delete[] c3qqp;
      delete[] c3qqbp;

      delete [] kernlog;
      delete [] kernlin;
    }
  else
    {
      delete[] c3qg_1;
      delete[] c3qq_1;
      delete[] c3qqb_1;
      delete[] c3qqp_1;
      delete[] c3qqbp_1;
      delete[] c3qg_2;
      delete[] c3qq_2;
      delete[] c3qqb_2;
      delete[] c3qqp_2;
      delete[] c3qqbp_2;

      delete [] kernlog_1;
      delete [] kernlin_1;
      delete [] kernlog_2;
      delete [] kernlin_2;
    }
}
