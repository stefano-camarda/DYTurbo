#ifndef specialfunctions_interface_h
#define specialfunctions_interface_h

#include "fcomplex.h"

extern "C"
{
  void ch12n_ (int &n, fcomplex &z, int &nm, fcomplex *chf1, fcomplex *chd1, fcomplex *chf2, fcomplex *chd2 );
  void jyzo_ (int &n, int &nt, double *rj0, double *rj1, double *ry0, double *ry1);
}

#endif
