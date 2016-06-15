#include "init.h"
#include "settings.h"
#include "interface.h"
//#include "integr.h"
#include "phasespace.h"
#include "resintegr.h"
#include "ctintegr.h"
#include "cubacall.h"
#include "plotter.h"
#include "rapint.h"
#include "clock_real.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

using namespace std;

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

int main( int argc , char * argv[])
{

  double begin_time, end_time;

  /***********************************/
  //Initialization
  try {
      dyturboinit(argc,argv);
  } catch (QuitProgram &e) {
      // print help and die
      printf("%s \n",e.what());
      return 0;
  }
  /***********************************/

  /*****************************************/
  //If using the DYRES approximation for PDFs, make sure that the PDF fit is initialised in the same way
  //Need to throw a random point according to a breit wigner, which is used to determine xtauf in the PDF fit
  double costh, m, qt, y;
  int mode = 0;
  if (opts.doRES || opts.doCT){
      if (opts_.approxpdf_ == 1) {
          srand(opts.rseed);
          double wsqmin = pow(opts.mlow,2);
          double wsqmax = pow(opts.mhigh,2);
          double x1=((double)rand()/(double)RAND_MAX);
          double m2,wt;
          breitw_(x1,wsqmin,wsqmax,opts.rmass,opts.rwidth,m2,wt);
          cout << "Initialise PDF fit with mass = " << sqrt(m2) << " xtauf = " << m2 / opts.sroot << endl;
          costh = 0.; m = sqrt(m2); qt = 1; y = 0;
          for (int ipdf=0; ipdf<opts.totpdf; ipdf++){
              //setpdf_(&ipdf);
              //setmellinpdf_(&ipdf);
              resumm_(costh,m,qt,y,mode);
          }
      }
      else {
          costh = 0.; m = opts.rmass; qt = 1; y = 0;
          for (int ipdf=0; ipdf<opts.totpdf; ipdf++){
              //setpdf_(&ipdf);
              //setmellinpdf_(&ipdf);
              resumm_(costh,m,qt,y,mode);
          }
      }
  }
  double f[opts.totpdf];
  ctint_(costh,m,qt,y,mode,f);
  /****************************************/
  

  // Cuba integration
  double value, error, totval, toterror2;
  vector <double> vals;
  vector <double> totvals;

  cout << endl << "Start integration of";
  if (opts.doRES  ) cout << " resummation";
  if (opts.doVV   ) cout << " double virtual";
  if (opts.doCT   ) cout << " counterterm";
  if (opts.doVJ   ) cout << " V+j finite order";
  if (opts.doLO   ) cout << " V+j LO";
  if (opts.doREAL ) cout << " V+j NLO real";
  if (opts.doVIRT ) cout << " V+J NLO virt";
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
      phasespace::setbounds(opts.mlow, opts.mhigh, *qit, *(qit+1), *yit, *(yit+1) );//  opts.ylow, opts.yhigh);
      print_qtbin();
      print_ybin();
      save_qtbin();
      save_ybin();
      double  bb_time = clock_real();

      // resummation
      if (opts.doRES) {
          double b_time = clock_real();
          if (opts.resint2d) {
	    if (opts.resumcpp)
	      //C++ resum
	      rapint::cache(phasespace::ymin, phasespace::ymax);
	      //end C++ resum
	    else	    
	      cacheyrapint_(phasespace::ymin, phasespace::ymax);
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
          doublevirtintegr(vals, error);
          double e_time = clock_real();
	  value = vals[0];
          normalise_result(value,error);
          print_result(value,error,b_time,e_time);
	  hists.FillResult( plotter::VV , value, error, e_time-b_time );
          totval += value;
          toterror2 += error*error;
	  vadd(totvals,vals);
      }
      // counter term
      if (opts.doCT) {
          double b_time = clock_real();
	  if (opts.ctint2d) ctintegr2d(vals, error);
	  if (opts.ctint3d) ctintegr3d(vals, error);
	  if (opts.ctintvegas6d) ctintegrMC(vals, error);
	  if (opts.ctintvegas8d) ctintegr(vals, error);
          double e_time = clock_real();
	  value = vals[0];
          normalise_result(value,error);
          print_result(value,error,b_time,e_time);
          hists.FillResult( plotter::CT , value, error, e_time-b_time );
          totval += value;
          toterror2 += error*error;
	  vadd(totvals,vals);
      }
      //analytical
      if (opts.doVJ) {
          double b_time = clock_real();
          vjintegr3d(value, error);
          double e_time = clock_real();
	  vals.clear();
	  vals.push_back(value);
          normalise_result(value,error);
          print_result(value,error,b_time,e_time);
          hists.FillResult( plotter::VJ , value, error, e_time-b_time );
          totval += value;
          toterror2 += error*error;
	  vadd(totvals,vals);
      }
      // leading order
      if (opts.doLO) {
          double b_time = clock_real();
          lowintegr(vals, error);
          double e_time = clock_real();
	  value = vals[0];
          normalise_result(value,error);
          print_result(value,error,b_time,e_time);
          hists.FillResult( plotter::LO , value, error, e_time-b_time );
          totval += value;
          toterror2 += error*error;
	  vadd(totvals,vals);
      }
      // real part
      if (opts.doREAL) {
          double b_time = clock_real();
          realintegr(vals, error);
          double e_time = clock_real();
	  value = vals[0];
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
	  double e_time = clock_real();
	  value = vals[0];
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

void normalise_result(double &value, double &error){
    TotXSec+=value;
    if (opts.ptbinwidth)
      {
	value /= phasespace::qtmax - phasespace::qtmin;
	error /= phasespace::qtmax - phasespace::qtmin;
      }
    if (opts.ybinwidth)
      {
	value /= phasespace::ymax  - phasespace::ymin;
	error /= phasespace::ymax  - phasespace::ymin;
      }
}

void print_head(){
    print_line();
    cout << "| ";
    cout << setw(13) << "qt bin" << " | ";
    cout << setw(13) << "y bin"  << " | ";
    if (opts.doRES ) cout << setw(38) << "resummed "      << " | ";
    if (opts.doVV  ) cout << setw(38) << "double virtual "<< " | ";
    if (opts.doCT  ) cout << setw(38) << "counter term "  << " | ";
    if (opts.doVJ  ) cout << setw(38) << "V+j fixed ord "           << " | ";
    if (opts.doLO  ) cout << setw(38) << "V+j LO "        << " | ";
    if (opts.doREAL) cout << setw(38) << "V+j NLO real "     << " | ";
    if (opts.doVIRT) cout << setw(38) << "V+j NLO virt "  << " | ";
    if (true       ) cout << setw(38) << "TOTAL "         << " | ";
    cout << endl;
    print_line();
}
void print_line(){
    int N = 18+18;
    if (opts.doRES ) N += 41;
    if (opts.doVV )  N += 41;
    if (opts.doCT  ) N += 41;
    if (opts.doVJ  ) N += 41;
    if (opts.doLO  ) N += 41;
    if (opts.doREAL) N += 41;
    if (opts.doVIRT) N += 41;
    if (true       ) N += 41;
    cout<<string(N,'-').c_str() <<endl;
}

void print_qtbin(){
    //      2 + 5 + 3 + 5 + 3 = 18
   cout << "| " << setw(5) << phasespace::qtmin << " - " << setw(5) << phasespace::qtmax << " | " <<  flush; 
}
void print_ybin(){
    //      2 + 5 + 3 + 5 + 3 = 18
   cout << setw(5) << phasespace::ymin << " - " << setw(5) << phasespace::ymax << " | " <<  flush; 
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
  outfile << phasespace::qtmin << " " << phasespace::qtmax << " " << flush; 
}
void save_ybin()
{
  outfile << phasespace::ymin << " " << phasespace::ymax << " " << flush; 
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
  int rc = system("mv results.txt temp.txt && column -t temp.txt > results.txt && rm temp.txt");
}

void vadd(vector <double> &totvals, vector <double> vals)
{
  vector<double>::iterator it = totvals.begin();
  vector<double>::iterator i = vals.begin();
  for (; it != totvals.end(); it++,i++)
    *it += *i;
}
