#ifndef dynnlo_interface_h
#define dynnlo_interface_h

extern "C"
{
  //interface to DYNNLO fortran functions and common blocks
  double c2qqreg_(double& z);
  double c2qqp_(double& z);
  double c2qqb_(double& z);
  double c2qg_(double& z);
}

#endif
