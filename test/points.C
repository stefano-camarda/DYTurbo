#include "config.h"
#include <iostream>
#include <iomanip>

#include "init.h"
#include "integr.h"
#include "settings.h"
#include "interface.h"
#include "resintegr.h"
#include "ctintegr.h"
#include "finintegr.h"
#include "finitemapping.h"
#include "printsettings.h"
#include "resint.h"

using namespace std;

#include "lines.C"

void test_resum_speed(double costh,double m,double qt,double y,int mode);
void test_ct_speed(double costh,double m,double qt,double y,int mode);

int main( int argc , const char * argv[])
{
  double begin_time, end_time;

  /***********************************/
  //Initialization
  string conf_file;
  if (argc>1) {
      conf_file = argv[1];
  }
  SMparameters();
  opts.readfromfile(conf_file.c_str());
  opts.initDyresSettings();
  gaussinit_();
  iniflavreduce_();
  dyturboinit();
  rescinit_();
  //bins.init();
  bins.readfromfile(conf_file.c_str());
  /***********************************/

  /***********************************/
  //Initialization
  ///@todo: print out EW parameters and other settings
  // just a check
  opts.dumpAll();
  printsettings();
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
      double wsqmin = pow(opts.mlow,2);
      double wsqmax = pow(opts.mhigh,2);
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
  mline();
  //mlinebw();
  //xline();
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
    value = resumm_(costh,m,qt,y,mode);
    end_time = clock_real();
    cout << setw(10) << "Result" << setw(15) << value
         << setw(10) << "time "  << setw(15) << float(end_time - begin_time) << "s" << endl;
    begin_time = clock_real();
    value = resint::rint(costh,m,qt,y,mode);
    end_time = clock_real();
    cout << setw(10) << "Result" << setw(15) << value
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

