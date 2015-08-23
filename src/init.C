#include "interface.h"
#include "LHAPDF/LHAPDF.h"
#include "settings.h"
#include "init.h"


void pdfini_(){
    printf(" CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n\n");
    LHAPDF::initPDFSetByName(opts.LHAPDFset);
    LHAPDF::initPDF(opts.LHAPDFmember);
    // initialization of alphas
    couple_.amz_=LHAPDF::alphasPDF(dymasses_.zmass_) ;
    printf("\n");
}

void initalphas()
{
  double pi=3.14159265358979;
  double twopi=2*pi;
  double fourpi=4*pi;

  couple_.amz_=LHAPDF::alphasPDF(dymasses_.zmass_) ;
  double scale = fabs(scale_.scale_);
  qcdcouple_.as_=dyalphas_(scale,couple_.amz_,nlooprun_.nlooprun_);
  qcdcouple_.ason2pi_=qcdcouple_.as_/twopi;
  qcdcouple_.ason4pi_=qcdcouple_.as_/fourpi;
  qcdcouple_.gsq_=fourpi*qcdcouple_.as_;
}
