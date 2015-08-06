#include <iostream>
#include <LHAPDF/LHAPDF.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <ctime>
#include <cuba.h>
#include <iomanip>
#include "config.h"

#include "integr.h"
#include "settings.h"
#include "interface.h"
#include "finintegr.h"
#include "finitemapping.h"
#include "cubacall.h"




using namespace std;

// void yline();
// void mline();
// void costhline();
// void ptline();
// void ptavar();
// void ptgvar();

void test_resum_speed(double costh,double m,double qt,double y,int mode);

int main( int argc , const char * argv[])
{

  cout << endl << endl;
  cout << " (        )              (           )   " << endl;
  cout << " )\\ )  ( /(   *   )      )\\ )  (  ( /(   " << endl;
  cout << "(()/(  )\\())` )  /(   ( (()/(( )\\ )\\())  " << endl;
  cout << " /(_))((_)\\  ( )(_))  )\\ /(_))((_|(_)\\   " << endl;
  cout << "(_))___ ((_)(_(_())_ ((_|_))((_)_  ((_)  " << endl;
  cout << " |   \\ \\ / /|_   _| | | | _ \\| _ )/ _ \\  " << endl;
  cout << " | |) \\ V /   | | | |_| |   /| _ \\ (_) | " << endl;
  cout << " |___/ |_|    |_|  \\___/|_|_\\|___/\\___/  " << endl;
  cout << "                                   v " << VERSION << endl;
  cout << endl;
  cout << endl;

  clock_t begin_time, end_time;

  /***********************************/
  //initialise settings
  //opts.init();
  string conf_file = "input/CT10nlo_settings.in";
  if (argc>1) {
      conf_file = argv[1];
  }
  opts.readfromfile(conf_file.c_str());
  opts.initDyresSettings();
  dyinit_();
  //  setup_();
  //bins.init();
  bins.readfromfile(conf_file.c_str());
  //force number of cores to 0 (no parallelization)
  cubacores(4,1000000);
  //To do: print out EW parameters and other settings
  // just a check
  //opts.dumpAll();
  //return 0;
  /***********************************/

  double costh, m, qt, y;
  double value, error;

  /**************************************/
  //Checks for resummed cross section
  //  std::cout << std::setprecision(15);
  int mode = 0;
  costh = 0.3; m = 91; qt = 1; y = 0;
  test_resum_speed(costh,m,qt,y,mode);

  costh = 0.1; m = 91; qt = 5; y = 0.2;
  test_resum_speed(costh,m,qt,y,mode);

  costh = 0.5; m = 70; qt = 10; y = 1.5;
  test_resum_speed(costh,m,qt,y,mode);

  costh = -1.0; m = 110; qt = 20; y = -2.5;
  test_resum_speed(costh,m,qt,y,mode);

  costh = 0.1; m = 91; qt = 5; y = 1.5;
  test_resum_speed(costh,m,qt,y,mode);

  costh = 0.1; m = 91; qt = 5; y = 0.2;
  test_resum_speed(costh,m,qt,y,mode);

  //ptline();
  //yline();
  //mline();
  //ptavar();
  //ptgvar();
  /**************************************/


  /**************************************/
  //Checks for finite order cross section
  //born level variables (6 dimensions)
  // m, qt, y, costh;
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

  //work in progress for rewriting the counterterm with the same integration structure as the resummed part
  /*
  double cthmom0 = 0;
  double cthmom1 = 0;
  double cthmom2 = 0;
  cout << countterm_(costh,m,qt,y,alpha,beta,cthmom0,cthmom1,cthmom2) << endl;
  */
  /**************************************/

  // return 0;

  // Cuba integration
  cout << endl;
  cout << "Start integration" << endl;
  begin_time = clock();
  for (vector<double>::iterator qit = bins.qtbins.begin(); qit != bins.qtbins.end()-1; qit++)
    {
      //Set integration boundaries
      setbounds(opts.mlow, opts.mhigh, *qit, *(qit+1), opts.ylow, opts.yhigh);
      clock_t b_time = clock();
      if (opts.int2d) integr2d(value, error);
      if (opts.int3d) integr3d(value, error);
      if (opts.int4d) integr4d(value, error);
      clock_t e_time = clock();
      value = value / (*(qit+1) - *qit);
      error = error / (*(qit+1) - *qit);
      cout << setw(3) << "bin" << setw(5) << *qit << setw(2) << "-" << setw(5) << *(qit+1)
	   << setw(10) << value << setw(5) << "+/-" << setw(10) << error
	   << setw(6) << "time" << setw(10) <<  float(e_time - b_time) / CLOCKS_PER_SEC << endl;
    }
  end_time = clock();
  cout << endl;
  cout << setw(10) << "time "  << setw(15) << float(end_time - begin_time) / CLOCKS_PER_SEC << endl;

  cout << endl;
  cout << "Start integration of counterterm" << endl;
  begin_time = clock();
  for (vector<double>::iterator qit = bins.qtbins.begin(); qit != bins.qtbins.end()-1; qit++)
    {
      //Set integration boundaries
      setbounds(opts.mlow, opts.mhigh, *qit, *(qit+1), opts.ylow, opts.yhigh);
      clock_t b_time = clock();
      ctintegr(value, error);
      clock_t e_time = clock();
      value = value / (*(qit+1) - *qit);
      error = error / (*(qit+1) - *qit);
      cout << setw(3) << "bin" << setw(5) << *qit << setw(2) << "-" << setw(5) << *(qit+1)
	   << setw(10) << value << setw(5) << "+/-" << setw(10) << error
	   << setw(6) << "time" << setw(10) <<  float(e_time - b_time) / CLOCKS_PER_SEC << endl;
    }
  end_time = clock();
  cout << endl;
  cout << setw(10) << "time "  << setw(15) << float(end_time - begin_time) / CLOCKS_PER_SEC << endl;

  if (opts.order == 1)
    {
      cout << "Start integration of Z+j LO" << endl;
      begin_time = clock();
      for (vector<double>::iterator qit = bins.qtbins.begin(); qit != bins.qtbins.end()-1; qit++)
	{
	  //Set integration boundaries
	  setbounds(opts.mlow, opts.mhigh, *qit, *(qit+1), opts.ylow, opts.yhigh);
	  clock_t b_time = clock();
	  lowintegr(value, error);
	  clock_t e_time = clock();
	  value = value / (*(qit+1) - *qit);
	  error = error / (*(qit+1) - *qit);
	  cout << setw(3) << "bin" << setw(5) << *qit << setw(2) << "-" << setw(5) << *(qit+1)
	       << setw(10) << value << setw(5) << "+/-" << setw(10) << error
	       << setw(6) << "time" << setw(10) <<  float(e_time - b_time) / CLOCKS_PER_SEC << endl;
	}
      end_time = clock();
      cout << endl;
      cout << setw(10) << "time "  << setw(15) << float(end_time - begin_time) / CLOCKS_PER_SEC << endl;
    }
  if (opts.order == 2)
    {
      cout << endl;
      cout << "Start integration of real" << endl;
      begin_time = clock();
      for (vector<double>::iterator qit = bins.qtbins.begin(); qit != bins.qtbins.end()-1; qit++)
	{
	  //Set integration boundaries
	  setbounds(opts.mlow, opts.mhigh, *qit, *(qit+1), opts.ylow, opts.yhigh);
	  clock_t b_time = clock();
	  realintegr(value, error);
	  clock_t e_time = clock();
	  value = value / (*(qit+1) - *qit);
	  error = error / (*(qit+1) - *qit);
	  cout << setw(3) << "bin" << setw(5) << *qit << setw(2) << "-" << setw(5) << *(qit+1)
	       << setw(10) << value << setw(5) << "+/-" << setw(10) << error
	       << setw(6) << "time" << setw(10) <<  float(e_time - b_time) / CLOCKS_PER_SEC << endl;
	}
      end_time = clock();
      cout << endl;
      cout << setw(10) << "time "  << setw(15) << float(end_time - begin_time) / CLOCKS_PER_SEC << endl;

      cout << endl;
      cout << "Start integration of virtual" << endl;
      begin_time = clock();
      for (vector<double>::iterator qit = bins.qtbins.begin(); qit != bins.qtbins.end()-1; qit++)
	{
	  //Set integration boundaries
	  setbounds(opts.mlow, opts.mhigh, *qit, *(qit+1), opts.ylow, opts.yhigh);
	  clock_t b_time = clock();
	  virtintegr(value, error);
	  clock_t e_time = clock();
	  value = value / (*(qit+1) - *qit);
	  error = error / (*(qit+1) - *qit);
	  cout << setw(3) << "bin" << setw(5) << *qit << setw(2) << "-" << setw(5) << *(qit+1)
	       << setw(10) << value << setw(5) << "+/-" << setw(10) << error
	       << setw(6) << "time" << setw(10) <<  float(e_time - b_time) / CLOCKS_PER_SEC << endl;
	}
      end_time = clock();
      cout << endl;
      cout << setw(10) << "time "  << setw(15) << float(end_time - begin_time) / CLOCKS_PER_SEC << endl;
    }

  return 0;
}

void test_resum_speed(double costh,double m,double qt,double y,int mode){
    clock_t begin_time, end_time;
    double value;
    begin_time = clock();
    value = resumm_(costh,m,qt,y,mode);
    end_time = clock();
    cout << setw(10) << "Result" << setw(15) << value
         << setw(10) << "time "  << setw(15) << float(end_time - begin_time) / CLOCKS_PER_SEC << endl;
    return;
}
