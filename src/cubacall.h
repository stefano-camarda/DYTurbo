#ifndef cubacall_h
#define cubacall_h

#include <vector>
using namespace std;

void resintegr2d(double &res, double &err);
void resintegr3d(double &res, double &err);
void resintegrMC(double &res, double &err);

void lowintegr(double &res, double &err);
void realintegr(vector <double> &res, double &err);
void virtintegr(vector <double> &res, double &err);
void ctintegr(double &res, double &err);
void ctintegr3d(double &res, double &err);
void ctintegr2d(double &res, double &err);
void doublevirtintegr(double &res, double &err);

void exitfun(void * input, const int &core);

#endif
