#ifndef clenshawcurtisrules_h
#define clenshawcurtisrules_h

#define CCNMAX 500

namespace cc
{
  extern double xxx[CCNMAX][CCNMAX];
  extern double www[CCNMAX][CCNMAX];
  extern double cosw[65][65];
  extern double sinw[65][65];

  extern const double xxx9[9],www9[9];
  extern const double xxx33[33],www33[33];
  extern const double xxx65[65],www65[65];
  
  extern void init();
  extern void setw(double omega);
}

#endif
