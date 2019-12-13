#ifndef bequad_h
#define bequad_h

#define BQNMAX 40

namespace bq
{
  extern double xxx[BQNMAX][BQNMAX];
  extern double www[BQNMAX][BQNMAX];

  extern const double xxx1[1],www1[1];
  extern const double xxx2[2],www2[2];
  extern const double xxx3[3],www3[3];
  extern const double xxx4[4],www4[4];
  extern const double xxx5[5],www5[5];
  extern const double xxx20[20],www20[20];
  extern const double xxx40[40],www40[40];
  
  extern void init();
}

#endif
