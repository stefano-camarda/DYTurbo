      DOUBLE PRECISION HIST(150,100),XHIS(150,100),HDEL(150),HMIN(150)
     &,HMAX(150),HAVG(150),HINT(150),HSIG(150)
      COMMON/HISTOR/HIST,XHIS,HDEL,HMIN,HMAX,HAVG,HINT,HSIG

      CHARACTER TITLE*100,BOOK*3
      COMMON/HISTOC/BOOK(150),TITLE(150)                         
      INTEGER NBIN(150),IHIS(150,100),IUSCORE(150),IOSCORE(150),
     & IENT(150),NHIST
      COMMON/HISTOI/NBIN,IHIS,IUSCORE,IOSCORE,IENT,NHIST
