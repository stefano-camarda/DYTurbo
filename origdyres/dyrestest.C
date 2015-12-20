#include "config.h"
#include <iostream>
#include <iomanip>

#include <ctime>

using namespace std;

//#include "lines.C"

void test_resum_speed(double costh,double m,double qt,double y);

extern "C"
{
  double resumm_(double &costh, double &mm, double &qtt, double &yy);
  double realint_(double r[22], double &wgt) {};
  double virtint_(double r[22], double &wgt) {};
  double countint_(double r[22], double &wgt) {};
  double lowinthst_(double r[22], double &wgt) {};
  double lowint_(double r[22], double &wgt) {};
  void setup_();
  void dyinit_();
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
  costh = 0.3; m = 91; qt = 1; y = 0;
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
  cout << setw(10) << "Result" << setw(15) << value
       << setw(10) << "time "  << setw(15) << float(end_time - begin_time) << "s" << endl;
  return;
}
