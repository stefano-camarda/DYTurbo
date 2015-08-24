 
       COMPLEX*16 function wlgama(x)
       IMPLICIT COMPLEX*16 (A - Z)
       WLGAMA = (X - cmplx(0.5,0.0)) * LOG (X)  -
     1  X + cmplx(0.91893853,0.0) + 
     2  cmplx(1.0,0.0)/(cmplx(12.0,0.0)* X)
     3  * (cmplx(1.0,0.0)- cmplx(1.0,0.0)/(cmplx(30.,0.)* X*X) 
     4  * (cmplx(1.,0.)- cmplx(1.,0.0)/(cmplx(3.5,0.0) * X*X)
     5  * (cmplx(1.,0.0)- cmplx(4.,0.0)/(cmplx(3.,0.0)* X*X))))
       return
       end
 
       Real*8 function dlgama(x)
       IMPLICIT real*8 (A - Z)
       dlgama = (X - (0.5)) * LOG (X)  -
     1  X + (0.91893853) + 
     2  (1.0)/((12.0)* X)
     3  * ((1.0)- (1.0)/((30.)* X*X) 
     4  * ((1.)- (1.)/((3.5) * X*X)
     5  * ((1.)- (4.)/((3.)* X*X))))
       return
       end
 
 
*
* $Id: dadapt.F,v 1.1.1.1 1996/04/01 15:02:13 mclareni Exp $
*
* $Log: dadapt.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:13  mclareni
* Mathlib gen
*
*
c#include "gen/pilot.h"
      SUBROUTINE DADAPT(F,A,B,NSEG,RELTOL,ABSTOL,RES,ERR)
c#if !defined(CERNLIB_DOUBLE)
c#include "gen/imp128.inc"
c      CHARACTER*6 NAME
c      NAME = 'DADAPT'
c      CALL MTLPRT(NAME,'D102',
c     +'not available on this machine - see documentation')
c      RETURN
c      END
c#endif
c#if defined(CERNLIB_DOUBLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
C     RES = Estimated Integral of F from A to B,
C     ERR = Estimated absolute error on RES.
C     NSEG  specifies how the adaptation is to be done:
C        =0   means use previous binning,
C        =1   means fully automatic, adapt until tolerance attained.
C        =n>1 means first split interval into n equal segments,
C             then adapt as necessary to attain tolerance.
C     The specified tolerances are:
C            relative: RELTOL ;  absolute: ABSTOL.
C        It stops when one OR the other is satisfied, or number of
C        segments exceeds NDIM.  Either TOLA or TOLR (but not both!)
C        can be set to zero, in which case only the other is used.
 
      EXTERNAL F
 
      PARAMETER (NDIM=100)
      PARAMETER (R1 = 1, HF = R1/2)
 
      DIMENSION XLO(NDIM),XHI(NDIM),TVAL(NDIM),TERS(NDIM)
      SAVE XLO,XHI,TVAL,TERS,NTER
      DATA NTER /0/
 
      IF(NSEG .LE. 0)  THEN
       IF(NTER .EQ. 0) THEN
        NSEGD=1
        GO TO 2
       ENDIF
       TVALS=0
       TERSS=0
       DO 1 I = 1,NTER
       CALL DGS56P(F,XLO(I),XHI(I),TVAL(I),TE)
       TERS(I)=TE**2
       TVALS=TVALS+TVAL(I)
       TERSS=TERSS+TERS(I)
    1  CONTINUE
       ROOT= SQRT(2*TERSS)
       GO TO 9
      ENDIF
      NSEGD=MIN(NSEG,NDIM)
    2 XHIB=A
      BIN=(B-A)/NSEGD
      DO 3 I = 1,NSEGD
      XLO(I)=XHIB
      XLOB=XLO(I)
      XHI(I)=XHIB+BIN
      IF(I .EQ. NSEGD) XHI(I)=B
      XHIB=XHI(I)
      CALL DGS56P(F,XLOB,XHIB,TVAL(I),TE)
      TERS(I)=TE**2
    3 CONTINUE
      NTER=NSEGD
      DO 4 ITER = 1,NDIM
      TVALS=TVAL(1)
      TERSS=TERS(1)
      DO 5 I = 2,NTER
      TVALS=TVALS+TVAL(I)
      TERSS=TERSS+TERS(I)
    5 CONTINUE
      ROOT= SQRT(2*TERSS)
      IF(ROOT .LE. ABSTOL .OR. ROOT .LE. RELTOL*ABS(TVALS)) GO TO 9
      IF(NTER .EQ. NDIM) GO TO 9
      BIGE=TERS(1)
      IBIG=1
      DO 6 I = 2,NTER
      IF(TERS(I) .GT. BIGE) THEN
       BIGE=TERS(I)
       IBIG=I
      ENDIF
    6 CONTINUE
      NTER=NTER+1
      XHI(NTER)=XHI(IBIG)
      XNEW=HF*(XLO(IBIG)+XHI(IBIG))
      XHI(IBIG)=XNEW
      XLO(NTER)=XNEW
      CALL DGS56P(F,XLO(IBIG),XHI(IBIG),TVAL(IBIG),TE)
      TERS(IBIG)=TE**2
      CALL DGS56P(F,XLO(NTER),XHI(NTER),TVAL(NTER),TE)
      TERS(NTER)=TE**2
    4 CONTINUE
    9 RES=TVALS
      ERR=ROOT
      RETURN
      END
c#endif
*
* $Id: dgs56p.F,v 1.1.1.1 1996/04/01 15:02:14 mclareni Exp $
*
* $Log: dgs56p.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:14  mclareni
* Mathlib gen
*
*
c#include "gen/pilot.h"
      SUBROUTINE DGS56P(F,A,B,RES,ERR)
c#if !defined(CERNLIB_DOUBLE)
c#include "gen/imp128.inc"
c      CHARACTER*6 NAME
c      NAME = 'DGS56P'
c      CALL MTLPRT(NAME,'D106',
c     +'not available on this machine - see documentation')
c      RETURN
c      END
c#endif
c#if defined(CERNLIB_DOUBLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      PARAMETER (R1 = 1, HF = R1/2)
      DIMENSION X5(5),W5(5),X6(6),W6(6)
 
      DATA (X5(I),W5(I),I=1,5)
     1/4.6910077030668004D-02, 1.1846344252809454D-01,
     2 2.3076534494715846D-01, 2.3931433524968324D-01,
     3 5.0000000000000000D-01, 2.8444444444444444D-01,
     4 7.6923465505284154D-01, 2.3931433524968324D-01,
     5 9.5308992296933200D-01, 1.1846344252809454D-01/
 
      DATA (X6(I),W6(I),I=1,6)
     1/3.3765242898423989D-02, 8.5662246189585178D-02,
     2 1.6939530676686775D-01, 1.8038078652406930D-01,
     3 3.8069040695840155D-01, 2.3395696728634552D-01,
     4 6.1930959304159845D-01, 2.3395696728634552D-01,
     5 8.3060469323313225D-01, 1.8038078652406930D-01,
     6 9.6623475710157601D-01, 8.5662246189585178D-02/
 
      RANG=B-A
      E5=0
      E6=0
      DO 1 I = 1,5
      E5=E5+W5(I)*F(A+RANG*X5(I))
      E6=E6+W6(I)*F(A+RANG*X6(I))
    1 CONTINUE
      E6=E6+W6(6)*F(A+RANG*X6(6))
      RES=HF*(E6+E5)*RANG
      ERR=ABS((E6-E5)*RANG)
      RETURN
      END



      SUBROUTINE DCSPLN(N,X,M,Y,NDIM,MODE,A,B,C,D)
      CHARACTER NAMEN*(*),NAMET*(*)
      CHARACTER*80 ERRTXT

      PARAMETER (NAMEN = 'DCSPLN', NAMET = 'DCSPLT')
      LOGICAL LNT

      DIMENSION X(0:*),Y(0:NDIM,*)
      DIMENSION A(NDIM,*),B(NDIM,*),C(NDIM,*),D(NDIM,*)

      PARAMETER (Z1 = 1, C1 = 3*Z1/2, C2 = Z1/3, C3 = 2*Z1/3)
      PARAMETER (C4 = Z1/2, C5 = Z1/6, C6 = 2*Z1/15, C7 = 7*Z1/60)

      LNT=.FALSE.
      GO TO 50


      ENTRY DCSPNT(N,X,M,Y,NDIM,MODE,A,B,C,D)

      LNT=.TRUE.
   50 IF(N .LT. 2 .OR. M .LT. 1 .OR. NDIM .LT. N .OR.
     1   MODE .NE. 0. AND. MODE .NE. 1) THEN
       IF(N .LT. 2) THEN
        WRITE(ERRTXT,101) N
c        IF(.NOT.LNT) CALL MTLPRT(NAMEN,'E211.1',ERRTXT)
c        IF(     LNT) CALL MTLPRT(NAMET,'E211.1',ERRTXT)
       ENDIF
       IF(M .LT. 1) THEN
        WRITE(ERRTXT,102) M
c        IF(.NOT.LNT) CALL MTLPRT(NAMEN,'E211.2',ERRTXT)
c        IF(     LNT) CALL MTLPRT(NAMET,'E211.2',ERRTXT)
       ENDIF
       IF(NDIM .LT. N) THEN
        WRITE(ERRTXT,103) NDIM,N
c        IF(.NOT.LNT) CALL MTLPRT(NAMEN,'E211.3',ERRTXT)
c        IF(     LNT) CALL MTLPRT(NAMET,'E211.3',ERRTXT)
       ENDIF
       IF(MODE .NE. 0 .AND. MODE .NE. 1) THEN
        WRITE(ERRTXT,104) MODE
c        IF(.NOT.LNT) CALL MTLPRT(NAMEN,'E211.4',ERRTXT)
c        IF(     LNT) CALL MTLPRT(NAMET,'E211.4',ERRTXT)
       ENDIF
       RETURN
      ENDIF
      DO 1 I = 1,N
    1 D(I,1)=X(I)-X(I-1)
      DO 2 K = 1,M
      DO 2 I = 1,N
    2 A(I,K)=(Y(I,K)-Y(I-1,K))/D(I,1)
      IF(MODE .EQ. 1) THEN
       IF(N .EQ. 2) THEN
        T1=1/(D(1,1)+D(2,1))
        DO 3 K = 1,M
    3   C(2,K)=T1*(A(2,K)-A(1,K))
       ELSE
        DO 4 I = 2,N
    4   B(I,1)=1/(D(I,1)+D(I-1,1))
        DO 5 K = 1,M
    5   C(1,K)=0
        B(1,1)=1

        DO 6 I = 2,N-1
        T1=3*B(I,1)
        T2=B(I,1)*D(I-1,1)
        T3=1/(2+T2*B(I-1,1))
        B(I,1)=(T2-1)*T3
        DO 6 K = 1,M
    6   C(I,K)=(T1*(A(I,K)-A(I-1,K))-T2*C(I-1,K))*T3

        T1=3*B(N,1)
        T2=B(N,1)*D(N-1,1)
        T3=1/(3-T2*(1-B(N-1,1)))
        DO 7 K = 1,M
    7   C(N,K)=(T1*(A(N,K)-A(N-1,K))-T2*C(N-1,K))*T3
       END IF

       DO 8 I = N-1,2,-1
       T1=B(I,1)
       DO 8 K = 1,M
    8  C(I,K)=T1*C(I+1,K)+C(I,K)
       DO 9 K = 1,M
    9  C(1,K)=C(2,K)
       IF(.NOT.LNT) THEN
        DO 10 K = M,1,-1
        B(1,K)=A(1,K)-C(2,K)*D(1,1)
        D(1,K)=0
        DO 11 I = 2,N-1
        B(I,K)=A(I,K)-C2*(C(I+1,K)+2*C(I,K))*D(I,1)
   11   D(I,K)=(C(I+1,K)-C(I,K))/(3*D(I,1))
        B(N,K)=A(N,K)-C(N,K)*D(N,1)
   10   D(N,K)=0
        DO 12 K = 1,M
        DO 12 I = 1,N
   12   A(I,K)=Y(I-1,K)
       ENDIF
      ELSE
       IF(N .EQ. 2) THEN
        T1=C1/(D(1,1)+D(2,1))
        DO 23 K = 1,M
   23   C(2,K)=T1*(A(2,K)-A(1,K))
       ELSE
        DO 24 I = 2,N
   24   B(I,1)=1/(D(I,1)+D(I-1,1))
        DO 25 K = 1,M
   25   C(1,K)=0
        B(1,1)=0

        DO 26 I = 2,N-1
        T1=3*B(I,1)
        T2=B(I,1)*D(I-1,1)
        T3=1/(2+T2*B(I-1,1))
        B(I,1)=(T2-1)*T3
        DO 26 K = 1,M
   26   C(I,K)=(T1*(A(I,K)-A(I-1,K))-T2*C(I-1,K))*T3

        T1=3*B(N,1)
        T2=B(N,1)*D(N-1,1)
        T3=1/(2+T2*B(N-1,1))
        DO 27 K = 1,M
   27   C(N,K)=(T1*(A(N,K)-A(N-1,K))-T2*C(N-1,K))*T3
       END IF

       DO 28 I = N-1,2,-1
       T1=B(I,1)
       DO 28 K = 1,M
   28  C(I,K)=T1*C(I+1,K)+C(I,K)
       DO 29 K = 1,M
   29  C(1,K)=0
       IF(.NOT.LNT) THEN
        T1=C2*D(1,1)
        T2=C2/D(1,1)
        T3=C3*D(N,1)
        T4=C2/D(N,1)
        DO 30 K = M,1,-1
        B(1,K)=A(1,K)-T1*C(2,K)
        D(1,K)=T2*C(2,K)
        DO 31 I = 2,N-1
        B(I,K)=A(I,K)-C2*(C(I+1,K)+2*C(I,K))*D(I,1)
   31   D(I,K)=(C(I+1,K)-C(I,K))/(3*D(I,1))
        B(N,K)=A(N,K)-T3*C(N,K)
   30   D(N,K)=-T4*C(N,K)
        DO 32 K = 1,M
        DO 32 I = 1,N
   32   A(I,K)=Y(I-1,K)
       ENDIF
      ENDIF

      IF(LNT) THEN
       DO 41 K = 1,M
       T1=D(1,1)**2
       A(1,K)=C4*(Y(0,K)+Y(1,K)-C5*(C(1,K)+C(2,K))*T1)*D(1,1)
       B(1,K)=C2*(Y(0,K)+C4*Y(1,K)-(C6*C(1,K)+C7*C(2,K))*T1)*T1
       DO 42 I = 2,N-1
       T1=D(I,1)**2
       A(I,K)=A(I-1,K)+
     1        C4*(Y(I-1,K)+Y(I,K)-C5*(C(I,K)+C(I+1,K))*T1)*D(I,1)
   42  B(I,K)=B(I-1,K)+C2*(Y(I-1,K)+C4*Y(I,K)-
     1                (C6*C(I,K)+C7*C(I+1,K))*T1)*T1+A(I-1,K)*D(I,1)
       T1=D(N,1)**2
       A(N,K)=A(N-1,K)+
     1        C4*(Y(N-1,K)+Y(N,K)-C(N,K)*C5*(1+MODE)*T1)*D(N,1)
   41  B(N,K)=B(N-1,K)+C2*(Y(N-1,K)+C4*Y(N,K)-
     1        C(N,K)*(C6+MODE*C7)*T1)*T1+A(N-1,K)*D(N,1)
      ENDIF
      RETURN
  101 FORMAT('ILLEGAL N = ',I5,' < 2')
  102 FORMAT('ILLEGAL M = ',I5,' < 1')
  103 FORMAT('ILLEGAL NDIM = ',I5,' < ',I5,' = N')
  104 FORMAT('ILLEGAL MODE = ',I5)
      END


      FUNCTION DGAUSS(F,A,B,EPS)
C.----------------------------------------------------------------------
C.
C.    GAUSS INTEGRAL OF THE FUNCTION F IN INTERVAL A,B
C.    LAST UPDATE: 10/04/88
C.
C.----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION W(12),X(12)
      EXTERNAL F
      DATA CONST/1.E-12/
      DATA W
     &/0.101228536290376, 0.222381034453374, 0.313706645877887,
     & 0.362683783378362, 0.027152459411754, 0.062253523938648,
     & 0.095158511682493, 0.124628971255534, 0.149595988816577,
     & 0.169156519395003, 0.182603415044924, 0.189450610455069/
      DATA X
     &/0.960289856497536, 0.796666477413627, 0.525532409916329,
     & 0.183434642495650, 0.989400934991650, 0.944575023073233,
     & 0.865631202387832, 0.755404408355003, 0.617876244402644,
     & 0.458016777657227, 0.281603550779259, 0.095012509837637/
C--
C--   INITIALISE
      DELTA=CONST*ABS(A-B)
      DGAUSS=0.
      AA=A
C--
C--   ITERATION LOOP
   10 Y=B-AA
C--
C--   EPSILON REACHED ??
      IF (ABS(Y).LE.DELTA) RETURN
   20 BB=AA+Y
      C1=0.5*(AA+BB)
      C2=C1-AA
      S8=0.
      S16=0.
      DO 30 I=1,4
        U=X(I)*C2
   30 S8=S8+W(I)*(F(C1+U)+F(C1-U))
      DO 40 I=5,12
        U=X(I)*C2
   40 S16=S16+W(I)*(F(C1+U)+F(C1-U))
      S8=S8*C2
      S16=S16*C2
      IF (ABS(S16-S8).GT.EPS*(1.0+ABS(S16))) GOTO 50
      DGAUSS=DGAUSS+S16
      AA=BB
      GOTO 10
   50 Y=0.5*Y
c      IF (ABS(Y).LE.DELTA) CALL ERRORD(34,' ',0.)
      GOTO 20
      END
