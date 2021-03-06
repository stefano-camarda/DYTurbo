#include "config.h"
#include <iostream>
#include <iomanip>

#include <ctime>
#include <math.h>

using namespace std;

//#include "lines.C"

void test_resum_speed(double costh,double m,double qt,double y);
void test_resum_int(double costh,double m,double qt,double y);

extern "C"
{
  double resumm_(double &costh, double &mm, double &qtt, double &yy);
  double realint_(double r[22], double &wgt) {};
  double virtint_(double r[22], double &wgt) {};
  double countint_(double r[22], double &wgt) {};
  double lowinthst_(double r[22], double &wgt);
  double lowint_(double r[22], double &wgt) {};
  void setup_();
  void dyinit_();
  extern struct {
    double facscale_;
  } facscale_;
  extern struct {
    double sroot_;
  } energy_;
  extern struct {
    double wsqmin_;
    double wsqmax_;
  } limits_;
  extern struct {
    int n2_;
    int n3_;
    double mass2_;
    double width2_;
    double mass3_;
    double width3_;
  } breit_;
}

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

  /***********************************/
  //Initialization
  dyinit_();
  /**************************************/

  //Checks for resummed cross section
  double costh, m, qt, y;
  //  std::cout << std::setprecision(15);
  costh = 0.; m = facscale_.facscale_; qt = 1; y = 0;
  test_resum_speed(costh,m,qt,y);

  for (double pt = 1.; pt <= 30; pt += 1.)
    {
      costh = 0.; m = facscale_.facscale_; qt = pt; y = 2.;
      //test_resum_speed(costh,m,qt,y);
      test_resum_int(costh,m,qt,y);
    }
  return 0;
  
  costh = 0.; m = facscale_.facscale_; qt = 80; y = 0;
  test_resum_speed(costh,m,qt,y);

  costh = 0.1; m = 91; qt = 5; y = 0.2;
  test_resum_speed(costh,m,qt,y);

  costh = 0.5; m = 70; qt = 10; y = 1.5;
  test_resum_speed(costh,m,qt,y);

  costh = -1.0; m = 110; qt = 20; y = -2.5;
  test_resum_speed(costh,m,qt,y);

  costh = 0.1; m = 91; qt = 5; y = 3.5;
  test_resum_speed(costh,m,qt,y);

  costh = 0.1; m = 91; qt = 5; y = 4.0;
  test_resum_speed(costh,m,qt,y);

  costh = 0.1; m = 91; qt = 5; y = 0.2;
  test_resum_speed(costh,m,qt,y);

  return 0;
  //costhline();
  //ptline();
  //  yline();
  //  mline();
  //mlinebw();
  //xline();
  //ptavar();
  //ptgvar();
  return 0;
}

void test_resum_speed(double costh,double m,double qt,double y){
  clock_t begin_time, end_time;
  double value;
  begin_time = clock();
  value = resumm_(costh,m,qt,y);
  end_time = clock();
  //  cout << setw(10) << "Result" << setw(15) << value
  cout << setw(10) << qt << setw(15) << value
       << setw(10) << "time "  << setw(15) << float(end_time - begin_time) << "s" << endl;
  return;
}

void test_resum_int(double costh,double m,double qt,double y){
  clock_t begin_time, end_time;
  double value;
  begin_time = clock();
  double r[4];
  double wgt = 1;

  double rmass = breit_.mass3_;
  double rwidth = breit_.width3_;
  double sqrts = energy_.sroot_; //7000.;
  double s = sqrts*sqrts;

  //convert invariant mass to 0-1 for Breit-Wigner weighting
  double almin, almax, tanal, al, x1;
  almin=atan((limits_.wsqmin_-rmass*rmass)/rmass/rwidth);
  almax=atan((limits_.wsqmax_-rmass*rmass)/rmass/rwidth);
  tanal = (m*m - rmass*rmass)/rmass/rwidth;
  al = atan(tanal);
  x1 = (al-almin)/(almax-almin);

  //rapidity
  double ymax=0.5*log(s/m/m);
  double x2 = (y + ymax)/(2.*ymax);

  //qt
  double qtmin=0.1;
  double cosh2y=pow((exp(y)+exp(-y))*0.5,2);
  double qtmax = sqrt(pow(s+m*m,2)/(4.*s*cosh2y) - m*m);
  double x3 = (qt-qtmin)/qtmax;

  //costh
  double x4 = (costh+1.)/2.;
  
  //phase space mapping
  r[0] = x1;                                           //q2
  r[1] = x2;                                          //y of the Z boson
  r[2] = x3;                                         //costh of the system (qt)
  r[3] = x4;                                          //costh of dilepton in Z rest frame
  
  value = lowinthst_(r,wgt);
  end_time = clock();
  //  std::cout << std::setprecision(15);
  //  cout << setw(10) << "Result" << setw(15) << value
  cout << setw(10) << qt << setw(15) << value
       << setw(10) << "time "  << setw(15) << float(end_time - begin_time) << "s" << endl;
  return;
}
