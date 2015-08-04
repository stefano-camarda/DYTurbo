#include "interface.h"
#include "LHAPDF/LHAPDF.h"
#include "settings.h"


void pdfini_(){
    printf("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n\n");
    printf("LHAPDF initialization:");
    LHAPDF::initPDFSetByName(opts.LHAPDFset);
    LHAPDF::initPDF(opts.LHAPDFmember);
    // initialization of alphas
    couple_.amz_=LHAPDF::alphasPDF(masses_.zmass_) ;
}
