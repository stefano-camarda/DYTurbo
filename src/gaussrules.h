#ifndef gaussrules_h
#define gaussrules_h

#define GRNMAX 500

namespace gr
{
  extern double xxx[GRNMAX][GRNMAX];
  extern double www[GRNMAX][GRNMAX];

  extern const double xxx2[2],www2[2];
  extern const double xxx3[3],www3[3];
  extern const double xxx4[4],www4[4];
  extern const double xxx5[5],www5[5];
  extern const double xxx8[8],www8[8];
  extern const double xxx10[10],www10[10];
  extern const double xxx20[20],www20[20];
  extern const double xxx24[24],www24[24];
  extern const double xxx40[40],www40[40];
  extern const double xxx50[50],www50[50];
  extern const double xxx64[64],www64[64];
  
  extern void init();
}

#endif
