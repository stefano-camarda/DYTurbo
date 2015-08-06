#ifndef cubacall_h
#define cubacall_h

void integr2d(double &res, double &err);
void integr3d(double &res, double &err);
void integr4d(double &res, double &err);

void lowintegr(double &res, double &err);
void realintegr(double &res, double &err);
void virtintegr(double &res, double &err);
void ctintegr(double &res, double &err);

#endif
