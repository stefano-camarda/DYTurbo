#include "config.h"
#include <iostream>
#include <LHAPDF/LHAPDF.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <cuba.h>
#include <iomanip>

#include "init.h"
#include "integr.h"
#include "resintegr.h"
#include "settings.h"
#include "interface.h"
#include "finintegr.h"
#include "finitemapping.h"
#include "cubacall.h"
#include "plotter.h"
#include "printsettings.h"

using namespace std;

#include "lines.C"

void test_resum_speed(double costh,double m,double qt,double y,int mode);
void test_ct_speed(double costh,double m,double qt,double y,int mode);
void print_head();
void print_line();
void print_qtbin();
void print_ybin();
void print_result(double val, double err, double btime , double etime);

ofstream outfile;
void open_file();
void save_qtbin();
void save_ybin();
void save_result(vector <double> vals, double err);
void close_file();

void vadd(vector <double> &totvals, vector <double> vals);

void normalise_result(double &value, double &error);

double TotXSec ;

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
  cubacores(opts.cubacores,1000000);   //< set number of cores (move this to cubainit)
  cubaexit((void (*)()) exitfun,NULL); //< merge at the end of the run
  // histogram output
  hists.Init();
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
  //mline();
  //mlinebw();
  //xline();
  //ptavar();
  //ptgvar();
  //return 0;
  /**************************************/


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

  // Cuba integration
  double value, error, totval, toterror2;
  vector <double> vals;
  vector <double> totvals;

  cout << endl << "Start integration of";
  if (opts.doRES  ) cout << " resummation";
  if (opts.doVV   ) cout << " double virtual";
  if (opts.doCT   ) cout << " counterterm";
  if (opts.doLO   ) cout << " finite order";
  if (opts.doREAL ) cout << " real part";
  if (opts.doVIRT ) cout << " virt part";
  cout << endl;

  print_head();
  open_file();
  begin_time = clock_real();
  TotXSec = 0.;
  for (vector<double>::iterator yit = bins.ybins.begin(); yit != bins.ybins.end()-1; yit++)
  {
      for (vector<double>::iterator qit = bins.qtbins.begin(); qit != bins.qtbins.end()-1; qit++)
      {

      totval=toterror2=0;
      totvals.clear();
      for (int i = 0; i < opts.totpdf; i++)
	totvals.push_back(0);
      //Set integration boundaries
      setbounds(opts.mlow, opts.mhigh, *qit, *(qit+1), *yit, *(yit+1) );//  opts.ylow, opts.yhigh);
      print_qtbin();
      print_ybin();
      save_qtbin();
      save_ybin();
      double  bb_time = clock_real();

      // resummation
      if (opts.doRES) {
          double b_time = clock_real();
          if (opts.resint2d) {
	    cacheyrapint_(ymin, ymax);
	    //C++ resum
	    /*
	    rapint::cache(ymin, ymax);
	    */
	    //end C++ resum
	    resintegr2d(value, error);
          }
          if (opts.resint3d) resintegr3d(value, error);
          if (opts.resintvegas) resintegrMC(value, error);
          double e_time = clock_real();
          normalise_result(value,error);
          print_result(value,error,b_time,e_time);
          hists.FillResult( plotter::Resum , value, error, e_time-b_time );
          totval += value;
          toterror2 += error*error;
      }
      // double virtual
      if (opts.doVV) {
          double b_time = clock_real();
          doublevirtintegr(value, error);
          double e_time = clock_real();
          normalise_result(value,error);
          print_result(value,error,b_time,e_time);
	  //          hists.FillResult( plotter::DoubleV , value, error, e_time-b_time );
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
          realintegr(vals, error);
	  value = vals[0];
          double e_time = clock_real();
          normalise_result(value,error);
          print_result(value,error,b_time,e_time);
          hists.FillResult( plotter::Real , value, error, e_time-b_time );
          totval += value;
          toterror2 += error*error;
	  vadd(totvals,vals);
      }
      // virt part
      if (opts.doVIRT) {
          double b_time = clock_real();
          virtintegr(vals, error);
	  value = vals[0];
	  double e_time = clock_real();
          normalise_result(value,error);
          print_result(value,error,b_time,e_time);
          hists.FillResult( plotter::Virt , value, error, e_time-b_time );
          totval += value;
          toterror2 += error*error;
	  vadd(totvals,vals);
      }
      // total
      double ee_time = clock_real();
      print_result (totval, sqrt(toterror2), bb_time, ee_time);
      save_result(totvals,sqrt(toterror2));
      hists.FillResult( plotter::Total ,  totval, sqrt(toterror2), ee_time-bb_time );
      cout << endl;
      }
    }
  print_line();
  close_file();
  end_time = clock_real();
  cout << endl;
  cout << setw(10) << "time "  << setw(15) << float(end_time - begin_time) << endl;

  hists.Finalise(TotXSec);

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
    TotXSec+=value;
    //value /= qtmax - qtmin;
    //error /= qtmax - qtmin;
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
    if (opts.doVV  ) cout << setw(38) << "double virtual "<< " | ";
    if (opts.doCT  ) cout << setw(38) << "counter term "  << " | ";
    if (opts.doREAL) cout << setw(38) << "real part "     << " | ";
    if (opts.doVIRT) cout << setw(38) << "virtual part "  << " | ";
    if (opts.doLO  ) cout << setw(38) << "Z+j LO "        << " | ";
    if (true       ) cout << setw(38) << "TOTAL "         << " | ";
    cout << endl;
    print_line();
}
void print_line(){
    int N = 18+18;
    if (opts.doRES ) N += 41;
    if (opts.doVV )  N += 41;
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

void open_file()
{
  outfile.open("results.txt");
  outfile << "qt1" << " ";
  outfile << "qt2" << " ";
  outfile << "y1"  << " ";
  outfile << "y2"  << " ";
  outfile << "xs(fb)" << " ";
  for (int i = 1; i < opts.totpdf; i++)
    {
      char pdf[10];
      sprintf(pdf, "pdf%d", i);
      outfile << pdf << " ";
    }
  outfile << "staterr" << endl;
  outfile << flush;
}
void save_qtbin()
{
  outfile << qtmin << " " << qtmax << " " << flush; 
}
void save_ybin()
{
  outfile << ymin << " " << ymax << " " << flush; 
}
void save_result(vector <double> vals, double err)
{
  double central = vals[0];
  for (vector<double>::iterator it = vals.begin(); it != vals.end(); it++)
    {
      if (it - vals.begin() == 0)
	outfile << *(vals.begin()) << " ";
      else
	outfile << *it-*vals.begin() << " ";
    }
  outfile << " " << err << endl;
  outfile << flush;
}
void close_file()
{
  outfile.close();
  system("mv results.txt temp.txt && column -t temp.txt > results.txt && rm temp.txt");
}

void vadd(vector <double> &totvals, vector <double> vals)
{
  vector<double>::iterator it = totvals.begin();
  vector<double>::iterator i = vals.begin();
  for (; it != totvals.end(); it++,i++)
    *it += *i;
}
