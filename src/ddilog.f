c     dilog function called from vjet/utils.f
      function Li2(x)
      implicit none
      double precision x,li2,ddilog,fli2
      external ddilog,fli2

      print *,x
c      Li2 = ddilog(x)
      Li2 = fli2(x)
      
      return
      end
      
c     dilog function from ancont of J. Bluemlein
      
      DOUBLE PRECISION FUNCTION FR1(X)
C     --------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  SUBSIDIARTY ROUTINE FOR  LI2(X)
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      Y=LOG(1.0D0-X)
      Y2=Y*Y
C
      T = (-1.D0+(-1.D0/4.D0+(-1.D0/36.D0+(1.D0/3600.D0+(
     &-1.D0/211680.D0+(1.D0/10886400.D0+(-1.D0/526901760.D0
     &+(691.D0/16999766784000.D0+(-1.D0/1120863744000.D0
     &+3617.D0/0.18140058832896D18*Y2)*Y2)*Y2)*Y2)*Y2)*Y2)
     &*Y2)*Y)*Y)*Y
C
      FR1=T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FLI2(X)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  LI2(X)
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      EXTERNAL FR1
C
      PI = 3.141592653589793238462643D0
      ZETA2 = PI**2/6.0D0       
C
      IF(X.LT.-1.0D0.OR.X.GT.1.0D0) GOTO 100
      IF(X.EQ.1.0D0)  GOTO 5
      IF(X.EQ.-1.0D0) GOTO 6
      IF(X.LT.0.0D0)  GOTO 1
      IF(X.EQ.0.0D0)  GOTO 2
      IF(X.GT.0.5D0)  GOTO 3
      T=FR1(X)
      GOTO 200
100   WRITE(6,*) 'FLI2 -> NOT ALLOWED,X=',X,'STOP ***'
      STOP
1     Y=-X
      IF(Y.GT.0.5D0) GOTO 4
      Y2=Y**2
      T=FR1(Y2)/2.0D0-FR1(Y)
      GOTO 200
2     T=0.0D0
      GOTO 200
3     XM=1.0D0-X
      T=-FR1(XM)-LOG(X)*LOG(XM)+ZETA2
      GOTO 200
4     YM=1.0D0-Y
      T1=-FR1(YM)-LOG(Y)*LOG(YM)+ZETA2
      Y2=Y**2
      IF(Y2.LT.0.5) THEN
      T2=FR1(Y2)
      ELSE
      T2=-FR1(1.0D0-Y2)-LOG(Y2)*LOG(1.0D0-Y2)+ZETA2
      ENDIF
      T=T2/2.0D0-T1
      GOTO 200
5     T=ZETA2
      GOTO 200
6     T=-ZETA2/2.0D0
      GOTO 200
C 
200   FLI2=T
C
      RETURN
      END


c     dilog function from algorithm C332 of CERNLIB
      
      DOUBLE PRECISION FUNCTION DDILOG(X)
      IMPLICIT NONE

      DOUBLE PRECISION X,Y,T,S,A,PI3,PI6,ZERO,ONE,HALF,MALF,MONE,MTWO
      DOUBLE PRECISION C(0:18),H,ALFA,B0,B1,B2

      INTEGER I
      DATA ZERO /0.0D0/, ONE /1.0D0/
      DATA HALF /0.5D0/, MALF /-0.5D0/, MONE /-1.0D0/, MTWO /-2.0D0/
      DATA PI3 /3.28986 81336 96453D0/, PI6 /1.64493 40668 48226D0/

      DATA C( 0) / 0.42996 69356 08137 0D0/
      DATA C( 1) / 0.40975 98753 30771 1D0/
      DATA C( 2) /-0.01858 84366 50146 0D0/
      DATA C( 3) / 0.00145 75108 40622 7D0/
      DATA C( 4) /-0.00014 30418 44423 4D0/
      DATA C( 5) / 0.00001 58841 55418 8D0/
      DATA C( 6) /-0.00000 19078 49593 9D0/
      DATA C( 7) / 0.00000 02419 51808 5D0/
      DATA C( 8) /-0.00000 00319 33412 7D0/
      DATA C( 9) / 0.00000 00043 45450 6D0/
      DATA C(10) /-0.00000 00006 05784 8D0/
      DATA C(11) / 0.00000 00000 86121 0D0/
      DATA C(12) /-0.00000 00000 12443 3D0/
      DATA C(13) / 0.00000 00000 01822 6D0/
      DATA C(14) /-0.00000 00000 00270 1D0/
      DATA C(15) / 0.00000 00000 00040 4D0/
      DATA C(16) /-0.00000 00000 00006 1D0/
      DATA C(17) / 0.00000 00000 00000 9D0/
      DATA C(18) /-0.00000 00000 00000 1D0/

      double precision fli2
      external fli2

      
      IF(X .EQ. ONE) THEN
       DDILOG=PI6
       RETURN
      ELSE IF(X .EQ. MONE) THEN
       DDILOG=MALF*PI6
       RETURN
      END IF

      if (x.gt.-1.0d0.and.x.lt.1.0d0) then
         ddilog = fli2(x)
         return
      endif
      
      T=-X
      IF(T .LE. MTWO) THEN
       Y=MONE/(ONE+T)
       S=ONE
       A=-PI3+HALF*(LOG(-T)**2-LOG(ONE+ONE/T)**2)
      ELSE IF(T .LT. MONE) THEN
       Y=MONE-T
       S=MONE
       A=LOG(-T)
       A=-PI6+A*(A+LOG(ONE+ONE/T))
      ELSE IF(T .LE. MALF) THEN
       Y=(MONE-T)/T
       S=ONE
       A=LOG(-T)
       A=-PI6+A*(MALF*A+LOG(ONE+T))
      ELSE IF(T .LT. ZERO) THEN
       Y=-T/(ONE+T)
       S=MONE
       A=HALF*LOG(ONE+T)**2
      ELSE IF(T .LE. ONE) THEN
       Y=T
       S=ONE
       A=ZERO
      ELSE
       Y=ONE/T
       S=MONE
       A=PI6+HALF*LOG(T)**2
      END IF

      H=Y+Y-ONE
      ALFA=H+H
      B1=ZERO
      B2=ZERO
      DO I = 18,0,-1
         B0=C(I)+ALFA*B1-B2
         B2=B1
         B1=B0
      ENDDO
      DDILOG=-(S*(B0-H*B2)+A)
      RETURN
      END
