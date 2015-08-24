#include "config.h"
#include <iostream>
#include <LHAPDF/LHAPDF.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <sys/time.h>
#include <cuba.h>
#include <iomanip>

#include "integr.h"
#include "settings.h"
#include "interface.h"
#include "finintegr.h"
#include "finitemapping.h"
#include "cubacall.h"
#include "plotter.h"




using namespace std;

#include "lines.C"

void test_resum_speed(double costh,double m,double qt,double y,int mode);
void test_ct_speed(double costh,double m,double qt,double y,int mode);
void print_head();
void print_line();
void print_qtbin();
void print_ybin();
void print_result(double val, double err, double btime , double etime);
void normalise_result(double &value, double &error);

double clock_real();

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

  double begin_time, end_time;

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
  ///@todo: print out EW parameters and other settings
  // just a check
  opts.dumpAll();
  //return 0;
  // histogram output
  hists.Init();
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
  /**************************************/

  // Cuba integration
  cout << endl << "Start integration of";
  if (opts.doRES  ) cout << " resummation";
  if (opts.doCT   ) cout << " counterterm";
  if (opts.doLO   ) cout << " finite order";
  if (opts.doREAL ) cout << " real part";
  if (opts.doVIRT ) cout << " virt part";
  cout << endl;

  print_head();
  begin_time = clock_real();
  for (vector<double>::iterator yit = bins.ybins.begin(); yit != bins.ybins.end()-1; yit++)
  {
      for (vector<double>::iterator qit = bins.qtbins.begin(); qit != bins.qtbins.end()-1; qit++)
      {

      //Set integration boundaries
      totval=toterror2=0;
      setbounds(opts.mlow, opts.mhigh, *qit, *(qit+1), *yit, *(yit+1) );//  opts.ylow, opts.yhigh);
      print_qtbin();
      print_ybin();
      double  bb_time = clock_real();

      // resummation
      if (opts.doRES) {
          double b_time = clock_real();
          if (opts.int2d) {
              cacheyrapint_(ymin, ymax);
              integr2d(value, error);
          }
          if (opts.int3d) integr3d(value, error);
          if (opts.int4d) integr4d(value, error);
          double e_time = clock_real();
          normalise_result(value,error);
          print_result(value,error,b_time,e_time);
          hists.FillResult( plotter::Resum , value, error, e_time-b_time );
          totval += value;
          toterror2 += error*error;
      }
      // counter term
      if (opts.doCT) {
          double b_time = clock_real();
	  if (opts.ctint2d) ctintegr2d(value, error);
	  if (opts.ctint3d) ctintegr3d(value, error);
          if (opts.ctintvegas) ctintegr(value, error);
          double e_time = clock_real();
          normalise_result(value,error);
          print_result(value,error,b_time,e_time);
          hists.FillResult( plotter::CT , value, error, e_time-b_time );
          totval += value;
          toterror2 += error*error;
      }
      // leading order
      if (opts.doLO) {
          double b_time = clock_real();
          lowintegr(value, error);
          double e_time = clock_real();
          normalise_result(value,error);
          print_result(value,error,b_time,e_time);
          hists.FillResult( plotter::LO , value, error, e_time-b_time );
          totval += value;
          toterror2 += error*error;
      }
      // real part
      if (opts.doREAL) {
          double b_time = clock_real();
          realintegr(value, error);
          double e_time = clock_real();
          normalise_result(value,error);
          print_result(value,error,b_time,e_time);
          hists.FillResult( plotter::Real , value, error, e_time-b_time );
          totval += value;
          toterror2 += error*error;
      }
      // virt part
      if (opts.doVIRT) {
          double b_time = clock_real();
          virtintegr(value, error);
          double e_time = clock_real();
          normalise_result(value,error);
          print_result(value,error,b_time,e_time);
          hists.FillResult( plotter::Virt , value, error, e_time-b_time );
          totval += value;
          toterror2 += error*error;
      }
      // total
      double ee_time = clock_real();
      print_result (totval, sqrt(toterror2), bb_time, ee_time);
      hists.FillResult( plotter::Total ,  totval, sqrt(toterror2), ee_time-bb_time );
      cout << endl;
      }
    }
  print_line();
  end_time = clock_real();
  cout << endl;
  cout << setw(10) << "time "  << setw(15) << float(end_time - begin_time) << endl;

  hists.Finalise();

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
    return;
}


void normalise_result(double &value, double &error){
    value /= qtmax - qtmin;
    error /= qtmax - qtmin;
    //value /= ymax  - ymin;
    //error /= ymax  - ymin;
}

void test_ct_speed(double costh,double m,double qt,double y,int mode){
    double begin_time, end_time;
    double value;
    begin_time = clock_real();
    value = countterm_(costh,m,qt,y,mode);
    end_time = clock_real();
    cout << setw(10) << "Result" << setw(15) << value
         << setw(10) << "time "  << setw(15) << float(end_time - begin_time) << "s" << endl;
    return;
}


void print_head(){
    print_line();
    cout << "| ";
    cout << setw(13) << "qt bin" << " | ";
    cout << setw(13) << "y bin"  << " | ";
    if (opts.doRES ) cout << setw(38) << "resummed "      << " | ";
    if (opts.doCT  ) cout << setw(38) << "counter term "  << " | ";
    if (opts.doREAL) cout << setw(38) << "real part "     << " | ";
    if (opts.doVIRT) cout << setw(38) << "virtual part "  << " | ";
    if (opts.doLO  ) cout << setw(38) << "Z+j LO "        << " | ";
    if (true       ) cout << setw(38) << "TOTAL "         << " | ";
    cout << endl;
    print_line();

};
void print_line(){
    int N = 18+18;
    if (opts.doRES ) N += 41;
    if (opts.doCT  ) N += 41;
    if (opts.doREAL) N += 41;
    if (opts.doVIRT) N += 41;
    if (opts.doLO  ) N += 41;
    if (true       ) N += 41;
    cout<<string(N,'-').c_str() <<endl;
}

void print_qtbin(){
    //      2 + 5 + 3 + 5 + 3 = 18
   cout << "| " << setw(5) << qtmin << " - " << setw(5) << qtmax << " | " <<  flush; 
}
void print_ybin(){
    //      2 + 5 + 3 + 5 + 3 = 18
   cout << setw(5) << ymin << " - " << setw(5) << ymax << " | " <<  flush; 
}


void print_result(double val, double err, double btime , double etime){
    // 10 + 4 + 10 + 2 + 10 + 2 + 3 = 39
    cout << setw(10) << val
         << setw(4)  << "+/-"
         << setw(10) << err
         << " ("
         << setw(10) << float(etime - btime)
         << "s)"
         << " | "
         << flush;
}

double clock_real(){
    struct timeval now;
    gettimeofday(&now, NULL);
    return now.tv_sec+(now.tv_usec/1000000.0); // in sec with micro second precission
}
