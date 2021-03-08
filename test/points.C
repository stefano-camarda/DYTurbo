#include "config.h"

#include "dyturbo.h"
#include "omegaintegr.h"
#include "settings.h"
#include "interface.h"
#include "phasespace.h"
#include "resintegr.h"
#include "ctintegr.h"
#include "finintegr.h"
#include "bornintegr.h"
#include "finitemapping.h"
#include "resint.h"
#include "rapint.h"
#include "ctint.h"
#include "qtint.h"
#include "vjint.h"
#include "loint.h"
#include "scales.h"
#include "sudakovff.h"
#include "besselint.h"

#include <iostream>
#include <iomanip>

using namespace std;

#include "lines.C"

void test_resum_speed(double costh,double m,double qt,double y,int mode);
void test_ct_speed(double costh,double m,double qt,double y,int mode);

int main( int argc , char * argv[])
{
  double begin_time, end_time;

  /***********************************/
  //Initialization
  try {
      DYTurbo::Init(argc,argv);
  } catch (QuitProgram &e) {
      // print help and die
      printf("%s \n",e.what());
      return 0;
  }
  /***********************************/

  /***********************************/
  //Initialization
  ///@todo: print out EW parameters and other settings
  // just a check
  if (!opts.silent) opts.dumpAll();
  if (!opts.silent) DYTurbo::PrintTable::Settings();
  /***********************************/

  double costh, m, qt, y;
  //  std::cout << std::setprecision(15);
  int mode = 0;
  /*****************************************/
  //If using the DYRES approximation for PDFs, make sure that the PDF fit is initialised in the same way
  //Need to throw a random point according to a breit wigner, which is used to determine xtauf in the PDF fit
  if (opts_.approxpdf_ == 1)
    {
      srand(opts.rseed);
      double wsqmin = pow(phasespace::mmin,2);
      double wsqmax = pow(phasespace::mmax,2);
      double x1=((double)rand()/(double)RAND_MAX);
      double m2,wt;
      breitw_(x1,wsqmin,wsqmax,opts.rmass,opts.rwidth,m2,wt);
      cout << "Initialise PDF fit with mass = " << sqrt(m2) << " xtauf = " << m2 / opts.sroot << endl;
      costh = 0.; m = sqrt(m2); qt = 1; y = 0;
      test_resum_speed(costh,m,qt,y,mode);
    }
  /****************************************/
  
  /**************************************/
  //Checks for resummed cross section
  if (!opts.silent)
    {
      cout << "To match the numbers between fortran and C++ set:" << endl;
      cout << "mellinrule = 64     #number of nodes" << endl;
      cout << "zmax = 27.          #upper" << endl;
      cout << "cpoint = 1          " << endl;
      cout << "mellin1d = false          " << endl;
      cout << "bintaccuracy = 1.0e-2  #accuracy" << endl;
      cout << endl;
    }
  //  std::cout << std::setprecision(15);

  //  costh = 0.; m = opts.rmass; qt = 1; y = 0;
  //  test_resum_speed(costh,m,qt,y,mode);

//  costh = 0.1; m = 91; qt = 5; y = 0.2;
//  test_resum_speed(costh,m,qt,y,mode);
//
//  costh = 0.5; m = 70; qt = 10; y = 1.5;
//  test_resum_speed(costh,m,qt,y,mode);
//
//  costh = -1.0; m = 110; qt = 20; y = -2.5;
//  test_resum_speed(costh,m,qt,y,mode);
//
//  costh = 0.1; m = 91; qt = 5; y = 3.5;
//  test_resum_speed(costh,m,qt,y,mode);
//
//  costh = 0.1; m = 91; qt = 5; y = 4.0;
//  test_resum_speed(costh,m,qt,y,mode);
//
//  costh = 0.1; m = 91; qt = 5; y = 0.2;
//  test_resum_speed(costh,m,qt,y,mode);

//  costh = 0.1; m = 91; qt = 5; y = 0.2;
//  test_ct_speed(costh,m,qt,y,mode);

  double f[2];
  /*
  costh = 0.; m = opts.rmass; qt = 0; y = 0.; mode = 1;
  loint::lint(costh, m, y, mode, f);
  cout << "born cross section at m = " << m << " y = " << y << ": " << setprecision(12) << f[0]*2*M_PI*2*m/1000. << " pb" << endl;

  costh = 0.; m = opts.rmass; qt = 0; y = 2.4; mode = 1;
  loint::lint(costh, m, y, mode, f);
  cout << "born cross section at m = " << m << " y = " << y << ": " << setprecision(12) << f[0]*2*M_PI*2*m/1000. << " pb" << endl;


  costh = 0.; m = opts.rmass; qt = 0; y = 0; mode = 3;
  //opts.mellin1d = true;
  bool modlog = opts.modlog;
  int evmode = opts.evolmode;
  opts.modlog = true;
  opts.evolmode = 0;
  resint::rint(costh, m, qt, y, mode, f);
  cout << "born cross section at m = " << m << ": " << f[0]/(8./3.)*2*m/1000. << endl;
  opts.modlog = modlog;
  opts.evolmode = evmode;
  */


  double qtarr[52] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,
	       45,50,55,60,65,70,75,80,85,90,95,100};

  cout << endl << endl;
  cout << "#qt" << "\t" << "dsigma / dqT" << endl;
  std::cout << std::setprecision(15);
  //opts.mellin1d = false;
  if (opts.mellin1d)
    mode = 2;
  else
    mode = 1;
  for (int i = 0; i < 52; i++)
    {
      costh = 0.; m = opts.rmass; qt = qtarr[i]; y = 0;
      phasespace::set_m(m);
      phasespace::set_y(y);
      phasespace::set_qt(qt);
      scales::set(phasespace::m);
      scales::mcfm();
      phasespace::set_phiV(0.);
      omegaintegr::genV4p();
      resint::rint(costh,m,qt,y,mode,f);
      double value = f[0]/(8./3.)*2*m/1000;
      cout << qt << "\t" << value << endl;
    }

  /*
  cout << endl << endl;
  std::cout << std::setprecision(12);
  cout << "#qt" << "\t" << "dsigma / dqT" << endl;
  for (int i = 0; i < 52; i++)
    {
      costh = 0.; m = opts.rmass; qt = qtarr[i]; y = 0;
      phasespace::set_mqtyphi(m, qt, y);//set global variables to costh, m, qt, y
      phasespace::set_cth(costh);//set global variables to costh, m, qt, y
      phasespace::calcexpy();
      phasespace::calcmt();
      omegaintegr::genV4p();//generate boson 4-momentum, with m, qt, y and phi=0
      double vj = vjint::vint(m,qt,y);
      cout << qt << "\t" << vj*2*m/1000 << endl;
    }
  */
  
//  cout << endl;
//  cout << "#b" << "\t" << "S(b)" << endl;
//  double value = resint::rint(costh,m,qt,y,mode)/(8./3.)*2*m/1000;
//  for (int i = 1; i <= 200; i++)
//    {
//      complex <double> bb = 2*double(i)/200.;
//      complex <double> sudak=sudakov::sff(bb);
//      cout << real(bb) << "\t" << real(sudak) << endl;
//    }

//  //plot the b-integrand
//  costh = 0.; m = opts.rmass; qt = 0.5; y = 0; mode = 3;
//  //double f[2];
//  resint::rint(costh,m,qt,y,mode,f);
//  cout << "{" << endl;
//  cout << "TGraph *g = new TGraph();" << endl;
//  double bmax = resint::bc/opts.bcf;
//  for (int i = 0; i < 1000; i++)
//    {
//      double x = i/1000. * bmax * 2;
//      complex <double> b;
//      //      double b = b_mb*0.95+i*(b_mb*1.05 -b_mb*0.95)/100.;
//      //double b = 0.+i*(b_mb*5 -0)/1000.;
//      if (x < resint::bc)
//	b = x;
//      else
//	{
//	  complex <double> jacu = complex <double> (1.,tan(M_PI/opts.phibr));
//	  b = resint::bc + jacu*(x-resint::bc);
//	}
//      //      cout << b << "  " << b_mb << "  " << besselint::bint(b) << endl;;
//      cout << "g->SetPoint(g->GetN(), " << b << ", " << besselint::bint(b) << ");" << endl;
//    }
//  cout << "g->Draw();" << endl;
//  cout << "}" << endl;
  
  
  //costhline();
  //ptline();
  //ptomline();
  //yline();
  //mline();
  //mlinebw();
  //xline();
  //ptavar();
  //ptgvar();

  return 0;
  
  /**************************************/
  //Checks for finite order cross section
  //born level variables (6 dimensions)
  //double m, qt, y, costh;
  double phicm, phiZ;
  costh = 0.0;  m = 91.1876;  qt = 5.;  y = 0.5;  phicm = 0.0;  phiZ = 0.;

  //variables to be integrated (4 dimensions in total)
  //1 collinear PDF dimension, common to real, virtual and lowint
  double zcth;
  zcth = 0.5;   //polar angle of Z in the CM frame
  //3 dimensions of real radiation
  double mjj, phijj, costhjj;
  mjj = 10.;    //invariant mass of the dijet (set to 0 to have the correct virtual phase space mapping)
  phijj = 0.3;  //phi of the dijet in dijet rest frame
  costhjj = 0.1;//costh of the dijet in dijet rest frame
  //1 dimension of virtual
  double vz;
  vz = sqrt(0.95);
  //2 dimensions for the counterterm
  double alpha,beta;
  beta = 0.1;
  alpha = 0.1;

  return 0;
  //call function wrappers, which map the variables into the unity hypercube of the vegas integration
  cout << " check phase space mapping " << endl;
  dyreal(m, y, qt, phicm, phiZ, costh, zcth, mjj, costhjj, phijj);
  dyvirt(m, y, qt, phicm, phiZ, costh, zcth, vz);
  dyct(m, y, qt, phicm, phiZ, costh, alpha, beta);
  /**************************************/


  return 0;
}

void test_resum_speed(double costh,double m,double qt,double y,int mode){
    double begin_time, end_time;
    double value;
    begin_time = clock_real();
    double mcosth=-costh;
    value = resumm_(mcosth,m,qt,y,mode);
    end_time = clock_real();
    cout << setw(10) << "Result" << setw(10) << "(fortran)" << setw(15) << value
         << setw(10) << "time "  << setw(15) << float(end_time - begin_time) << "s" << endl;
    begin_time = clock_real();
    double f[opts.totpdf];
    resint::rint(costh,m,qt,y,mode,f);
    value = f[0];
    end_time = clock_real();
    cout << setw(10) << "Result" << setw(10) << "(C++)" << setw(15) << value
	 << setw(10) << "time "  << setw(15) << float(end_time - begin_time) << "s" << endl;
    return;
}

void test_ct_speed(double costh,double m,double qt,double y,int mode){
    double begin_time, end_time;
    double value;
    begin_time = clock_real();
    double f[opts.totpdf];
    double weight=1.;
    value = ctint_(costh,m,qt,y,mode,f);
    end_time = clock_real();
    cout << setw(10) << "Result" << setw(15) << value
         << setw(10) << "time "  << setw(15) << float(end_time - begin_time) << "s" << endl;
    return;
}

