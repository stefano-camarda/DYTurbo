#ifndef hankel_interface_h
#define hankel_interface_h

extern "C"
{
  void dhankl_(double &bmax, int &nb, int &nrel, double &tol, int &ntol, int *nord, double (*fun1)(double&), int *ijrel, double *dwork, double *dans, double *arg, int &nofun1, int &ierr);
}

#endif
