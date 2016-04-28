#ifndef switch_h
#define switch_h

namespace switching {
  void init();
  double swtch(double qt, double m);
  double qtlimit(double m);
  extern double k;
  extern double delta;
  extern int mode;
  extern const double cutoff;
  extern const double tolerance;
}

extern "C" {
  double  switching_(double &qt, double &m);
}

#endif
