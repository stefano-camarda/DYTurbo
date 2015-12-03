#include "pdf.h"
#include "interface.h"
#include "settings.h"

#include <LHAPDF/LHAPDF.h>
#include <LHAPDF/LHAGlue.h>
#include <math.h>

void pdfini_()
{
  printf(" ==Initialize PDF set from LHAPDF==\n\n");
  LHAPDF::initPDFSetByName(opts.LHAPDFset);
  LHAPDF::initPDF(opts.LHAPDFmember);
  // initialization of alphas
  couple_.amz_=LHAPDF::alphasPDF(dymasses_.zmass_) ;
  if (opts.PDFerrors && LHAPDF::numberPDF() > 1)
    {
      opts.totpdf = LHAPDF::numberPDF()+1;
      pdferropts_.pdferr_ = true;
      pdferropts_.totpdf_ = LHAPDF::numberPDF()+1;
    }
  else
    {
      opts.totpdf = 1;
      pdferropts_.pdferr_ = false;
      pdferropts_.totpdf_ = 1;
    }
  printf("\n");

  LHAPDF::Info& cfg = LHAPDF::getConfig();
  cfg.set_entry("Verbosity", 0);
  //setalphas();
}

void setpdf_(int& member)
{
  LHAPDF::initPDF(member);
  setalphas();
}

//set value of alphas
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
