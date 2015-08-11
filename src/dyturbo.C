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

#include "lines.C"

void test_resum_speed(double costh,double m,double qt,double y,int mode);
void print_head();
void print_line();
void print_qtbin(vector<double>::iterator it_qt);
void print_result(double val, double err, clock_t btime , clock_t etime);

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
  string conf_file;
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
  cubacores(opts.cubacores,1000000); // < move this to cubainit
  //To do: print out EW parameters and other settings
  // just a check
  opts.dumpAll();
  //return 0;
  /***********************************/

  double costh, m, qt, y;
  double value, error, totval, toterror2;

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

  //costhline();
  //ptline();
  //yline();
  //mline();
  //mlinebw();
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

  // decide what terms to calculate
  opts.doLO   = (opts.doLO   && opts.order == 1);
  opts.doREAL = (opts.doREAL && opts.order == 2);
  opts.doVIRT = (opts.doVIRT && opts.order == 2);
  // Cuba integration
  cout << endl << "Start integration of";
  if (opts.doRES  ) cout << " resummation";
  if (opts.doCT   ) cout << " counterterm";
  if (opts.doLO   ) cout << " finite order";
  if (opts.doREAL ) cout << " real part";
  if (opts.doVIRT ) cout << " virt part";
  cout << endl;

  print_head();
  begin_time = clock();
  for (vector<double>::iterator qit = bins.qtbins.begin(); qit != bins.qtbins.end()-1; qit++)
    {
      //Set integration boundaries
      totval=toterror2=0;
      setbounds(opts.mlow, opts.mhigh, *qit, *(qit+1), opts.ylow, opts.yhigh);
      print_qtbin(qit);
      clock_t bb_time = clock();

      // resummation
      if (opts.doRES) {
          clock_t b_time = clock();
          if (opts.int2d) integr2d(value, error);
          if (opts.int3d) integr3d(value, error);
          if (opts.int4d) integr4d(value, error);
          clock_t e_time = clock();
          value = value / (*(qit+1) - *qit);
          error = error / (*(qit+1) - *qit);
          print_result(value,error,b_time,e_time);
          totval += value;
          toterror2 += error*error;
      }
      // counter term
      if (opts.doCT) {
          clock_t b_time = clock();
          ctintegr(value, error);
          clock_t e_time = clock();
          value = value / (*(qit+1) - *qit);
          error = error / (*(qit+1) - *qit);
          print_result(value,error,b_time,e_time);
          totval += value;
          toterror2 += error*error;
      }
      // leading order
      if (opts.doLO) {
          clock_t b_time = clock();
          lowintegr(value, error);
          clock_t e_time = clock();
          value = value / (*(qit+1) - *qit);
          error = error / (*(qit+1) - *qit);
          print_result(value,error,b_time,e_time);
          totval += value;
          toterror2 += error*error;
      }
      // real part
      if (opts.doREAL) {
          clock_t b_time = clock();
          realintegr(value, error);
          clock_t e_time = clock();
          value = value / (*(qit+1) - *qit);
          error = error / (*(qit+1) - *qit);
          print_result(value,error,b_time,e_time);
          totval += value;
          toterror2 += error*error;
      }
      // virt part
      if (opts.doVIRT) {
          clock_t b_time = clock();
          virtintegr(value, error);
          clock_t e_time = clock();
          value = value / (*(qit+1) - *qit);
          error = error / (*(qit+1) - *qit);
          print_result(value,error,b_time,e_time);
          totval += value;
          toterror2 += error*error;
      }
      // total
      clock_t ee_time = clock();
      print_result (totval, sqrt(toterror2), bb_time, ee_time);
      cout << endl;
    }
  print_line();
  end_time = clock();
  cout << endl;
  cout << setw(10) << "time "  << setw(15) << float(end_time - begin_time) / CLOCKS_PER_SEC << endl;


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

void print_head(){
    print_line();
    cout << "| " << setw(13) << "qt bin" << " | ";
    if (opts.doRES ) cout << setw(36) << "resummed "      << " | ";
    if (opts.doCT  ) cout << setw(36) << "counter term "  << " | ";
    if (opts.doREAL) cout << setw(36) << "real part "     << " | ";
    if (opts.doVIRT) cout << setw(36) << "virtual part "  << " | ";
    if (opts.doLO  ) cout << setw(36) << "Z+j LO "        << " | ";
    if (true       ) cout << setw(36) << "TOTAL "          << " | ";
    cout << endl;
    print_line();

};
void print_line(){
    int N = 18;
    if (opts.doRES ) N += 39;
    if (opts.doCT  ) N += 39;
    if (opts.doREAL) N += 39;
    if (opts.doVIRT) N += 39;
    if (opts.doLO  ) N += 39;
    if (true       ) N += 39;
    cout<<string(N,'-').c_str() <<endl;
}

void print_qtbin(vector<double>::iterator it_qt){
    //      2 + 5 + 3 + 5 + 3 = 18
   cout << "| " << setw(5) << *it_qt << " - " << setw(5) << *(it_qt+1) << " | " <<  flush; 
}

void print_result(double val, double err, clock_t btime , clock_t etime){
    // 10 + 4 + 10 + 1 + 10 + 1 + 3 = 39
    cout << setw(10) << val
         << setw(4)  << "+/-"
         << setw(10) << err
         << "("
         << setw(10) << float(etime - btime) / CLOCKS_PER_SEC
         << ")"
         << " | "
         << flush;
}
