**********************************************************************
C    SIMPLE HISTOGRAMMING PACKAGE --  SIMPLIFIED VERSION OF HBOOK
C    BY Michelangelo Mangano    NOVEMBER 1988
C    LAST REVISED NOVEMBER 9, 1988  
c     (minor modifications by I Hinchliffe   1 May, 89)
C**********************************************************************
C
C Fills up to 100 histograms with up to 100 bins. 
C Gives a data file (to be specified in the calling program by assigning 
C a file name to unit 98) and a topdrawer file (to be specified in the 
C calling program by assigning a file name to unit 99).
C
C INITIALIZATION:
C Call once INIHIST; this just resets a few counters and logicals
C Call MBOOK(N,'TITLE',DEL,XMIN,XMAX) for each histogram to be booked.
C N (an integer) is the label of the histogram;
C 'TITLE' is the name of the histogram (no more then 100 characters);
C DEL (real*8) is the bin size;
C XMIN (real*8) is the lower limit of the first bin;
C XMAX (real*8)is the upper limit of the last  bin
C Example:
C      call mbook(2,'pt distribution',1.,10.,70.)
C This call initializes histogram number 2, called 'pt distribution';
C The bin size will be 1. (possibly GeV, if that's what you want), the
C first bin being  10.<x<11. and the last one being 69.<x<70.
C
C FILLING:
C When it's time, call MFILL(N,X,Y); this will add Y (real*8) to the bin 
C in which X (real*8) happens to be, within histogram N. 
C
C PLAYING AROUND:
C At the end of the day you may want to sum, divide, cancel, etc.etc.
C various histograms (bin by bin). Then you call MOPERA(I,'O',J,K,X,Y). 
C The 1-character string O can take the following values:
C +  : sums       X*(hist I) with Y*(hist J) and puts the result in hist K;
C -  : subtracts  X*(hist I) with Y*(hist J) and puts the result in hist K;
C *  : multiplies X*(hist I) with Y*(hist J) and puts the result in hist K;
C /  : divides    X*(hist I) with Y*(hist J) and puts the result in hist K;
C F  : multiplies hist I by the factor X, and puts the result in hist K;
C R  : takes the square root of  hist  I, and puts the result in hist K;if
C      the value at a given bin is less than or equal to 0., puts 0. in K
C S  : takes the square      of  hist  I, and puts the result in hist K;
C L  : takes the log_10 of  hist  I, and puts the result in hist K; if the
C      value at a given bin is less than or equal to 0., puts 0. in K
C M  : statistical analysis; if I contains the weights (let's say WGT),
C      J contains variable times weight (F*WGT) and K contains the
C      variable squared times the weight (F**2*WGT), then, after using 'M',
C      J will contain the average value of the variable <F> and K will 
C      contain the sigma of the average: sigma=sqrt(<F**2>-<F>**2).
C      If WGT=1. for all the entries, then it is enough to put I=J, and
C      it is not necessary to book a hist with the weights.
C V  : estimates errors for vegas evaluation of differential distributions.
C      Fill I with the values of
C      the functions do integrate times the Vegas weight (fun*wgt); fill
C      J with fun**2*wgt; then K will contain an estimate of the error
C      of the integration. Putting X=1/(#of iterations) performs the 
C      average over the iterations, and gives the right normalization to 
C      the differential distribution, I, and to the errors, K. J stays the same.
C
C FINAL ACCOUNTING:
C Now we can finalize our histograms; MFINAL(N) will calculate the integral
C of the histogram N, the mean value of the X variable and its RMS.
C If we now want to renormalize the hist's, we can call MNORM(N,X), which
C will normalize the integral to X  -- CAUTION: do not call MNORM before
C MFINAL, it will blow up.
C
C OUTPUT:
C To get a .dat file containing the values of the histograms, together with
C some information (like integral, mean values, etc.etc.) call MPRINT(N),
C for each hist N that you want in the .dat file. Before the call to MPRINT
C you want to open unit 98 and give it a name:
C     OPEN(UNIT=98,NAME='NAME.DAT',STATUS='NEW')
C If you want a topdrawer file with a plot of the hist values, call 
C MTOP(N,M,'X','Y','SCALE'). The points of the plot will be taken from histogram
C N, the error bars from histogram M. 'SCALE', character*(*), determines
C the scale for y, logarithmic or linear (SCALE=LOG,LIN).
C If you do not want error bars, keep
C a histogram of zeros, or just call a hist that had not been booked.
C X will appear as a 'bottom title', and Y will appear as a 'left title'.
C The top title is by default the name of the histogram itself.
C A little box below the plot will contain some information on the plot
C itself. Before calling MTOP,                     
C     OPEN(UNIT=99,NAME='NAME.TOP',STATUS='NEW')
c Empty histograms are not put out by MTOP.
C--------------------------------------------------------------------------
      BLOCK DATA HISTOSET
      include 'histo.f'
      data book/
     & ' NO','NO',' NO','NO',' NO','NO',' NO','NO',' NO','NO',
     & ' NO','NO',' NO','NO',' NO','NO',' NO','NO',' NO','NO',
     & ' NO','NO',' NO','NO',' NO','NO',' NO','NO',' NO','NO',
     & ' NO','NO',' NO','NO',' NO','NO',' NO','NO',' NO','NO',
     & ' NO','NO',' NO','NO',' NO','NO',' NO','NO',' NO','NO',
     & ' NO','NO',' NO','NO',' NO','NO',' NO','NO',' NO','NO',
     & ' NO','NO',' NO','NO',' NO','NO',' NO','NO',' NO','NO',
     & ' NO','NO',' NO','NO',' NO','NO',' NO','NO',' NO','NO',
     & ' NO','NO',' NO','NO',' NO','NO',' NO','NO',' NO','NO',
     & ' NO','NO',' NO','NO',' NO','NO',' NO','NO',' NO','NO',
     & ' NO','NO',' NO','NO',' NO','NO',' NO','NO',' NO','NO',
     & ' NO','NO',' NO','NO',' NO','NO',' NO','NO',' NO','NO',
     & ' NO','NO',' NO','NO',' NO','NO',' NO','NO',' NO','NO',
     & ' NO','NO',' NO','NO',' NO','NO',' NO','NO',' NO','NO',
     & ' NO','NO',' NO','NO',' NO','NO',' NO','NO',' NO','NO'/
      END

      SUBROUTINE MBOOK(N,TIT,DEL,XMIN,XMAX)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*(*) TIT
      include 'histo.f'
      NHIST=MAX(N,NHIST)
      TITLE(N)=TIT                     
      BOOK(N)='YES'
      HDEL(N)=DEL
      HMIN(N)=XMIN
      HMAX(N)=XMAX
      NNBIN=INT((XMAX-XMIN)/DEL)
      IF (NNBIN .GT. 150) THEN
      WRITE(6,*) XMAX,XMIN,DEL,NNBIN,' BIN SIZE TOO LARGE'
      DEL=(XMAX-XMIN)/99.d0
      NNBIN=INT((XMAX-XMIN)/DEL)
      ENDIF
      NBIN(N)=NNBIN
      IENT(N)=0
      IUSCORE(N)=0
      IOSCORE(N)=0
      HAVG(N)=0.d0
      HINT(N)=0.d0
      DO 1 I=1,NBIN(N)
      XHIS(N,I)=HMIN(N)+HDEL(N)*(DFLOAT(I)-0.5d0)
   1  HIST(N,I)=0.d0
      END

      SUBROUTINE MFILL(N,X,Y)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      include 'histo.f'
      I=INT((X-HMIN(N))/HDEL(N)+1)
      IF(I.GT.0.AND.I.LE.NBIN(N))  THEN
      IENT(N)=IENT(N)+1
      IHIS(N,I)=IHIS(N,I)+1
      HIST(N,I)=HIST(N,I)+Y/hdel(n)
c     we are renormalising the weights by the bin width
      ELSEIF(I.LE.0) THEN
      IUSCORE(N)=IUSCORE(N)+1
      ELSEIF(I.GT.NBIN(N)) THEN
      IOSCORE(N)=IOSCORE(N)+1
      ENDIF
      END

      SUBROUTINE MOPERA(I,OPER,J,K,X,Y)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER OPER*1
      include 'histo.f'
      IF(NBIN(I).NE.NBIN(J).AND.(OPER.EQ.'+'.OR.OPER.EQ.'-'.OR.OPER.EQ.
     &'*'.OR.OPER.EQ.'/'.OR.OPER.EQ.'M')) GO TO 10
      DO L=1,NBIN(I)
      IF(OPER.EQ.'+') THEN
      HIST(K,L)=X*HIST(I,L) + Y*HIST(J,L)
      ELSEIF(OPER.EQ.'-') THEN
      HIST(K,L)=X*HIST(I,L) - Y*HIST(J,L)
      ELSEIF(OPER.EQ.'*') THEN
      HIST(K,L)=X*HIST(I,L) * Y*HIST(J,L)
      ELSEIF(OPER.EQ.'/') THEN
        IF(Y.EQ.0.d0.OR.HIST(J,L).EQ.0.d0) THEN
          HIST(K,L)=0.d0
          ELSE
          HIST(K,L)=X*HIST(I,L) / (Y*HIST(J,L))
        ENDIF
      ELSEIF(OPER.EQ.'F') THEN
      HIST(K,L)=X*HIST(I,L)
      ELSEIF(OPER.EQ.'R') THEN
        IF(HIST(I,L).GT.0.d0) THEN
        HIST(K,L)=X*SQRT(HIST(I,L))
        ELSE
        HIST(K,L)=0.d0
        ENDIF
      ELSEIF(OPER.EQ.'S') THEN
      HIST(K,L)=X*HIST(I,L)**2
      ELSEIF(OPER.EQ.'l') THEN
        IF(HIST(I,L).EQ.0.d0.OR.J.EQ.0) THEN
             HIST(K,L)=0.d0
             ELSE
             HIST(K,L)=X*LOG10(Y*HIST(I,L))
        ENDIF
      ELSEIF(OPER.EQ.'M') THEN
        IF(I.NE.J) XNORM=HIST(I,L)
        IF(I.EQ.J) XNORM=DFLOAT(IHIS(J,L))
        IF(XNORM.NE.0.d0) THEN
        XAVG=HIST(J,L)/XNORM
        HIST(K,L)=SQRT(ABS(-XAVG**2+HIST(K,L)/XNORM)/DFLOAT(IHIS(I,L)))
        HIST(J,L)=XAVG 
        ELSE
        HIST(K,L)=0.d0
        HIST(J,L)=0.d0                           
        ENDIF
      ELSEIF(OPER.EQ.'V') THEN                 
        XAVG=HIST(I,L)*X
        XSQAVG=HIST(J,L)*X
        XNORM=DFLOAT(IHIS(I,L))*X
        IF(XNORM.NE.0.d0) THEN
        HIST(K,L)=SQRT(ABS(XSQAVG-XAVG**2)/XNORM)
        HIST(I,L)=XAVG
        ELSE
        HIST(K,L)=0.d0
        ENDIF
C
          ELSEIF(OPER.EQ.'A') THEN
            HIST(J,L) = HIST(J,L) + HIST(I,L)
            IHIS(J,L) = IHIS(J,L) + IHIS(I,L)
            HIST(K,L) = HIST(K,L) + HIST(I,L)**2
            IHIS(K,L) = IHIS(K,L) + 1
            IENT(K) = IENT(K)+1
            HIST(I,L) = 0
            IHIS(I,L) = 0
         IENT(J) = IENT(J)+IENT(I)
         IUSCORE(J) = IUSCORE(J) + IUSCORE(I)
         IOSCORE(J) = IOSCORE(J) + IOSCORE(I)
         IUSCORE(K) = IUSCORE(K) + 1
         IOSCORE(K) = IOSCORE(K) + 1
         IENT(I) = 0
         IUSCORE(I) = 0
         IOSCORE(I) = 0
        ELSEIF(OPER.EQ.'E') THEN
c If I contains the accumulated weights, J the accumulated squares of the
c weights and IHIS(J,1) the number of accumulated entries, 'E' will add
c the average value of I to K and will put in J the quadratic dispersion.
         IF(IHIS(J,1).NE.0) THEN
            XXX = 1./IHIS(J,1)
         ELSE
            XXX = 0
         ENDIF
            XSUM   = HIST(I,L)
            XSUMSQ = HIST(J,L)
            HIST(K,L)=HIST(K,L) + XXX*XSUM
            IHIS(K,L)=IHIS(K,L) + IHIS(I,L)
            HIST(J,L)=XXX*dSQRT(dABS(XSUMSQ-XSUM**2*XXX))
         IENT(K)=IENT(K)+IENT(I)
         IUSCORE(K) = IUSCORE(K)+IUSCORE(I)
         IOSCORE(K) = IOSCORE(K)+IOSCORE(I)
       ELSEIF(OPER.EQ.'Q') THEN
            HIST(K,L) = dSQRT(HIST(J,L)**2+HIST(I,L)**2)
C
      ELSE
      WRITE(98,5) OPER
   5  FORMAT(' ****** OPERATION ="',A1,'" UNKNOWN ********'/)
      RETURN
      ENDIF
      END DO
      RETURN
  10  WRITE(98,20) I,J
  20  FORMAT(' ****** INCOMPATIBLE OPERATION HIST ',I2,' &',I2,
     &                                                   '*******'/)
      END
     
      SUBROUTINE MZERO(N)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      include 'histo.f'
      BOOK(N)='RES'
      IENT(N)=0
      IUSCORE(N)=0
      IOSCORE(N)=0
      HAVG(N)=0.d0
      HINT(N)=0.d0
      DO 1 I=1,NBIN(N)
   1  HIST(N,I)=0.d0
      END

      SUBROUTINE MRESET(N)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      include 'histo.f'
      BOOK(N)='RES'
      END

      SUBROUTINE MFINAL(N)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      include 'histo.f'
      IF(BOOK(N).NE.'YES') RETURN
      AVG=0.d0
      XIN=0.d0                                
c      SIG=0.d0
      DO 1, J=1,NBIN(N)
      AVG=AVG+HIST(N,J)*XHIS(N,J)
   1  XIN=XIN+HIST(N,J)
      IF(XIN.EQ.0.d0) GO TO 10
      HAVG(N)=AVG/XIN
c      DO 2, J=1,NBIN(N)
c   2  SIG=HIST(N,J)*(XHIS(N,J)-HAVG(N))**2+SIG
c      IF(SIG.GE.0.)HSIG(N)=SQRT(SIG/XIN)
      HINT(N)=XIN*hdel(n)
      RETURN
  10  BOOK(N)=' NO'
      END               

      SUBROUTINE MNORM(N,X)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      include 'histo.f'
      IF(BOOK(N).NE.'YES')RETURN
      DO 1, I=1,NBIN(N)
    1 HIST(N,I)=HIST(N,I)/HINT(N)*X
      HINT(N)=X
      END

      SUBROUTINE MPRINT(N)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      include 'histo.f'
c      DATA INI/0/
c      IF(INI.EQ.0) THEN
c      CALL IDATE(IMON,IDAY,IYEAR)
c      CALL TIME(CTIME)
c      INI=1
c      ENDIF
      IF(BOOK(N).NE.'YES') then
      write(98,21) n
      RETURN
      end if
c      WRITE(98,7) N,IYEAR,IMON,IDAY,CTIME(1:5)
      WRITE(98,8) N
      WRITE(98,*) TITLE(N)
      WRITE(98,10) (XHIS(N,I),HIST(N,I),I=1,NBIN(N))
      WRITE(98,15) HAVG(N),HSIG(N),HINT(N)
      WRITE(98,20) IENT(N),IUSCORE(N),IOSCORE(N)
c    7 FORMAT(4X,'HIST = ',I3,'   19',I2,'-',I2,'-',I2,1X,A5/)
    8 FORMAT(4X,'HIST = ',I3)
   10 FORMAT(4X,2G13.6)
   15 FORMAT(/' AVG =',E10.3,4X,' RMS =',E10.3,' INTEGRAL =',E10.3,/)
   20 FORMAT('ENTRIES=',I10,1X,'U`FLOW=',I10,1X,'O`FLOW=',I10,//)
   21 FORMAT(' HISTOGRAM ',I3,' IS EMPTY')
      END

      SUBROUTINE MTOP(N,M,BTIT,LTIT,SCALE)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*(*) LTIT,BTIT,SCALE
      include 'histo.f'
      
      IF(BOOK(N).NE.'YES') RETURN
      WRITE(99,101) TITLE(N),TITLE(N),TITLE(N),SCALE,HMIN(N),HMAX(N)
  100 FORMAT( /1x,                               
     &' SET WINDOW Y 2.5 TO 7.'/,1X,
     &' SET WINDOW X 2.5 TO 10.'/,1X,
     &' SET SYMBOL 5O SIZE 1.8'/,1X,
     &' TITLE TOP ','"',A50,'"',/1X,
     &' TITLE BOTTOM ','"',A50,'"',/1X,
     &' TITLE LEFT ','"',A50,'"',/1X,
     &' SET SCALE Y ',A5,/1X,
     &' (SET TICKS TOP OFF)   '/1x,     
     &' SET LIMITS X ',F10.5,' ',F10.5,/1X,
     &' SET ORDER X Y DY ')
  101 FORMAT( /1x,                               
     &' SET WINDOW Y 2.5 TO 7.'/,1X,
     &' SET WINDOW X 2.5 TO 10.'/,1X,
     &' SET SYMBOL 5O SIZE 1.8'/,1X,
     &' TITLE TOP ','"',A,' distribution"',/1X,
     &' TITLE BOTTOM ','"',A,'"',/1X,
     &' TITLE LEFT ','"dS/d',A,' [fb]"',/1X,
     &' CASE       ','" G"',/1X,
     &' SET SCALE Y ',A5,/1X,
     &' (SET TICKS TOP OFF)   '/1x,     
     &' SET LIMITS X ',F10.5,' ',F10.5,/1X,
     &' SET ORDER X Y DY ')
      DO 1 J=1,NBIN(N)
      IF(HIST(N,J).EQ.0.) GO TO 1
      WRITE(99,'(3X,G13.6,2(2X,G13.6))')  
     &                            XHIS(N,J),HIST(N,J),HIST(M,J)
    1 CONTINUE
      WRITE(99,200)
  200 FORMAT('   PLOT')
      WRITE(99,300) HINT(N),HAVG(N),HSIG(N),IENT(N),IUSCORE(N)
     &   ,IOSCORE(N)
  300 FORMAT( /1x,                               
     &' BOX 7. 0.75 SIZE 9. 1.5'/,1X,
     &' SET WINDOW Y 0. TO 2.'/,1X,
     &' SET TITLE SIZE -1.5'/1X,
     &' TITLE 2.8 1.2 "INTGRL =',E12.5,'   AVGE =',E12.5,
     &             '   RMS =',E12.5,'"',/1X,
     &' TITLE 2.8 0.8 "Entries =',I9,2x,'U`flow =',I9,2X
     &                                 ,'O`flow =',I9,'"',/1X,
     &' SET TITLE SIZE -2')
      WRITE(99,400)
  400 FORMAT('   NEW PLOT')
      END
C*******************************************************************
C     END OF THE HISTOGRAMMING PACKAGE
C*******************************************************************
