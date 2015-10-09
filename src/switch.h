#ifndef switch_h
#define switch_h

namespace switching {
  void init(int md, double kk);
  double swtch(double qt, double m);
  double qtlimit(double m);
  extern double k;
  extern double cutoff;
  extern int mode;
}

extern "C" {
  double  switching_(double &qt, double &m);
}

#endif
