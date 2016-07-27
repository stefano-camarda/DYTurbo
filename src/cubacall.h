#ifndef cubacall_h
#define cubacall_h

const int last_iter=4;
static int ICALL=0; // number of integrand calls
void tell_to_grid_we_are_alive();

#include <vector>
using namespace std;

//resummation
void resintegr2d(vector <double> &res, double &err); //missing PDF variations
void resintegr3d(vector <double> &res, double &err); //missing PDF variations
void resintegrMC(vector <double> &res, double &err); //missing PDF variations

//fixed order born configuration
void bornintegr2d(vector <double> &res, double &err);  //missing PDF variations
void bornintegr3d(vector <double> &res, double &err);  //missing PDF variations
void bornintegrMC(vector <double> &res, double &err);

//fixed order V+j
void vjintegr3d(vector <double> &res, double &err); //missing PDF variations
void vjlointegr(vector <double> &res, double &err);
void vjrealintegr(vector <double> &res, double &err);
void vjvirtintegr(vector <double> &res, double &err);

void v2jintegr(vector <double> &res, double &err);

//counterterm
void ctintegr(vector <double> &res, double &err);
void ctintegrMC(vector <double> &res, double &err);
void ctintegr3d(vector <double> &res, double &err);
void ctintegr2d(vector <double> &res, double &err);

void exitfun(void * input, const int &core);

#endif
