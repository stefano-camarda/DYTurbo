#ifndef parton_h
#define parton_h

namespace parton
{
  //DYRES convention
  // bb cb sb db ub g u d s c b
  // -5 -4 -3 -2 -1 0 1 2 3 4 5
  enum partid {bb=0, cb=1, sb=2, db=3, ub=4, g=5, u=6, d=7, s=8, c=9, b=10};

  //LHAPDF convention
  // bb cb sb ub db g d u s c b
  // -5 -4 -3 -2 -1 0 1 2 3 4 5
  enum pdgid {Bb=0, Cb=1, Sb=2, Ub=3, Db=4, G=5, D=6, U=7, S=8, C=9, B=10};

  inline pdgid charge_conj(pdgid i)  {return pdgid(int((i-G))*(-1)+G);}
}
#endif
