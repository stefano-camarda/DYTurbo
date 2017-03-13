#include "config.h"

#include "src/dyturbo.h"
#include "src/omegaintegr.h"
#include "src/settings.h"
#include "src/interface.h"
#include "phasespace.h"
#include "resintegr.h"
#include "ctintegr.h"
#include "finintegr.h"
#include "finitemapping.h"
#include "resint.h"
#include "rapint.h"
#include "ctint.h"
#include "qtint.h"
#include "vjint.h"

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
  opts.dumpAll();
  DYTurbo::PrintTable::Settings();
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
  cout << "To match the numbers between fortran and C++ set:" << endl;
  cout << "mellinrule = 64     #number of nodes" << endl;
  cout << "zmax = 27.          #upper" << endl;
  cout << "bintaccuracy = 1.0e-2  #accuracy" << endl;
  cout << endl;
  //  std::cout << std::setprecision(15);

  costh = 0.; m = opts.rmass; qt = 1; y = 0;
  test_resum_speed(costh,m,qt,y,mode);

  costh = 0.1; m = 91; qt = 5; y = 0.2;
  test_resum_speed(costh,m,qt,y,mode);

  costh = 0.5; m = 70; qt = 10; y = 1.5;
  test_resum_speed(costh,m,qt,y,mode);

  costh = -1.0; m = 110; qt = 20; y = -2.5;
  test_resum_speed(costh,m,qt,y,mode);

  costh = 0.1; m = 91; qt = 5; y = 3.5;
  test_resum_speed(costh,m,qt,y,mode);

  costh = 0.1; m = 91; qt = 5; y = 4.0;
  test_resum_speed(costh,m,qt,y,mode);

  costh = 0.1; m = 91; qt = 5; y = 0.2;
  test_resum_speed(costh,m,qt,y,mode);

  costh = 0.1; m = 91; qt = 5; y = 0.2;
  test_ct_speed(costh,m,qt,y,mode);

  //costhline();
  //ptline();
  //yline();
  //mline();
  //mlinebw();
  xline();
  //ptavar();
  //ptgvar();

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
    value = resint::rint(costh,m,qt,y,mode);
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

