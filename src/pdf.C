#include "pdf.h"
#include "interface.h"
#include "settings.h"

#include <LHAPDF/LHAPDF.h>
#include <LHAPDF/LHAGlue.h>
#include <math.h>

void pdfini_(){
  printf(" CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n\n");
  LHAPDF::initPDFSetByName(opts.LHAPDFset);
  LHAPDF::initPDF(opts.LHAPDFmember);
  // initialization of alphas
  couple_.amz_=LHAPDF::alphasPDF(dymasses_.zmass_) ;
  printf("\n");
  //setalphas();
}

void setpdf_(int& member){
    LHAPDF::initPDF(member);
    // initialization of alphas
    setalphas();
}

void setalphas()
{
  couple_.amz_=LHAPDF::alphasPDF(dymasses_.zmass_) ;
  double scale = fabs(scale_.scale_);
  qcdcouple_.as_=dyalphas_(scale,couple_.amz_,nlooprun_.nlooprun_);
  qcdcouple_.ason2pi_=qcdcouple_.as_/(2*M_PI);
  qcdcouple_.ason4pi_=qcdcouple_.as_/(4*M_PI);
  qcdcouple_.gsq_=4*M_PI*qcdcouple_.as_;
}

void fdist_(int& ih, double& x, double& xmu, double fx[11])
{
  double fPDF[13];

  //set to zero if x out of range
  if (x > 1.)
    for (int i = -5; i <=5; i++)
      fx[5+i]=0.;
 
  LHAPDF::xfx(x,xmu,fPDF);
  if (ih == 1) //proton
    for (int i = -5; i <=5; i++)
      fx[5+i]=fPDF[6+i]/x;
  else if (ih == -1) //antiproton
    for (int i = -5; i <=5; i++)
      fx[5+i]=fPDF[6-i]/x;
}
