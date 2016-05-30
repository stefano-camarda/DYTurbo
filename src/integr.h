#ifndef integr_h
#define integr_h

extern "C" {
  void cthmoments_(double &cthmom0, double &cthmom1, double &cthmom2);
  void genv4p_();
}

const int last_iter=4;
static int ICALL=0; // number of integrand calls
void tell_to_grid_we_are_alive();
extern void genV4p(double m, double qt, double y, double phi);

#endif
