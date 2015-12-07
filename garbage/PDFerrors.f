      logical PDFerrors
      integer maxPDFsets,currentPDF
c--- 40 is my choice for the maximum number of PDF error sets
c--- 50 is my choice for the maximum number of dipoles
      double precision PDFxsec(0:40),PDFxsec_nd(0:40,0:50),
     . PDFwgt(0:40)
      common/PDFerrors/PDFerrors,maxPDFsets,PDFxsec,PDFwgt  
