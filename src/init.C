#include "interface.h"
#include "LHAPDF/LHAPDF.h"
#include "settings.h"


void pdfini_(){
    printf(" CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n\n");
    LHAPDF::initPDFSetByName(opts.LHAPDFset);
    LHAPDF::initPDF(opts.LHAPDFmember);
    // initialization of alphas
    couple_.amz_=LHAPDF::alphasPDF(dymasses_.zmass_) ;
    printf("\n");
}
