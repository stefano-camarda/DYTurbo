#ifndef cubacall_h
#define cubacall_h

#include <vector>
using namespace std;

void resintegr2d(double &res, double &err);
void resintegr3d(double &res, double &err);
void resintegrMC(double &res, double &err);

void vjintegr3d(double &res, double &err);
void lowintegr(vector <double> &res, double &err);
void realintegr(vector <double> &res, double &err);
void virtintegr(vector <double> &res, double &err);
void doublevirtintegr(vector <double> &res, double &err);

void ctintegr(vector <double> &res, double &err);
void ctintegrMC(vector <double> &res, double &err);
void ctintegr3d(vector <double> &res, double &err);
void ctintegr2d(vector <double> &res, double &err);

void exitfun(void * input, const int &core);

#endif
