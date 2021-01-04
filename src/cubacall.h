#ifndef cubacall_h
#define cubacall_h

//--> Make this a namelist, with a different name (f.i. integration)
//--> Wrap into general methods which can handle different integration routines (Cuba, cubature, Smolyak, Sparse grids, etc...)

const int last_iter=4;
static int ICALL=0; // number of integrand calls
void tell_to_grid_we_are_alive();

#include <vector>
using namespace std;

//resummation
void resintegr1d(vector <double> &res, double &err); //missing Ai
void resintegr2d(vector <double> &res, double &err); //missing Ai
void resintegr2d_my(vector <double> &res, double &err); //missing Ai
void resintegr3d(vector <double> &res, double &err); //(missing Ai)
void resintegrMC(vector <double> &res, double &err); //missing PDF variations

//fixed order born configuration
void bornintegr2d(vector <double> &res, double &err);    //Ai as dimensions
void bornintegrMC4d(vector <double> &res, double &err);  //missing PDF variations
void bornintegrMC6d(vector <double> &res, double &err);

//fixed order V+j
void vjintegr3d(vector <double> &res, double &err);   //missing PDF variations
void vjlointegr5d(vector <double> &res, double &err); //Ai as dimensions
void vjlointegr7d(vector <double> &res, double &err); //missing PDF variations
void vjlointegr(vector <double> &res, double &err);
void vjrealintegr(vector <double> &res, double &err);
void vjvirtintegr(vector <double> &res, double &err);

void v2jintegr(vector <double> &res, double &err);

//counterterm
void ctintegr(vector <double> &res, double &err);
void ctintegrMC(vector <double> &res, double &err);
void ctintegr3d(vector <double> &res, double &err);
void ctintegr2d(vector <double> &res, double &err); //PDF variations, swicth to Ai
void ctintegr1d(vector <double> &res, double &err);

//Cuba initialisation
namespace cuba
{
  void init();
  void initfun(void * input, const int &core);
  void exitfun(void * input, const int &core);
}

#endif
