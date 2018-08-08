#ifndef clenshawcurtisrules_h
#define clenshawcurtisrules_h

namespace cc
{
  extern double xxx[65][65];
  extern double www[65][65];
  extern double cosw[65][65];
  extern double sinw[65][65];

  extern const double xxx9[9],www9[9];
  extern const double xxx33[33],www33[33];
  extern const double xxx65[65],www65[65];
  
  extern void init();
  extern void setw(double omega);
}

#endif
