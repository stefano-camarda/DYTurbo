#ifndef ifunc_h
#define ifunc_h

namespace ifunc
{
  extern void test();
  extern double I1qq(double z);
  extern double I1qg(double z);

  extern double I2qqp(double z);
  extern double I2qqb(double z);
  extern double I2qq(double z);
  extern double I2qg(double z);

  extern void I3(double z);

  extern double i3qqp;
  extern double i3qqbp;
  extern double i3qq;
  extern double i3qqb;
  extern double i3qg;
  
  extern double Id33(double z);
  extern double Istqqp(double z);
  extern double Istqqb(double z);
  extern double Istqq(double z);
  extern double I3qqp(double z);
  extern double I3qqbp(double z);
  extern double I3qq(double z);
  extern double I3qqb(double z);
  extern double I3qg(double z);
}

#endif
