#ifndef pdf_h
#define pdf_h

extern void setalphas();
extern void setg();

extern "C" {
  void fdist_(int& ih, double& x, double& xmu, double fx[11]);
  void setpdf_(int& member);
  void setmellinpdf_(int& member);
}



#endif
