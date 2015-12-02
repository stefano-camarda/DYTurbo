#ifndef pdf_h
#define pdf_h

extern void setalphas();

extern "C" {
  void fdist_(int& ih, double& x, double& xmu, double fx[11]);
}

#endif
