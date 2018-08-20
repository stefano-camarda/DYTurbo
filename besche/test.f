      PROGRAM MAIN
C     DRIVER PROGRAM FOR SUBROUTINE BESCHE
C     COMPUTATION OF THE INTEGRAL OF EXP(-2*X)*J(ALFA(I)*X,NU)
C     OVER (0,30)
C     ALFA(I)=1,10,100,1000,10000,100000
C     NU=0,1,...,10
      IMPLICIT NONE
      DOUBLE PRECISION BESINT
      DOUBLE PRECISION AB,ALFA,DIF,EXACT,FUN,RESULT
      DOUBLE PRECISION DSQRT
      INTEGER I,INFO,K,N,NU,NUM,NUU
      DIMENSION ALFA(6),INFO(6),RESULT(6)
      EXTERNAL FUN

c      print*,BESINT(0.1D0),BESINT(0.2D0),BESINT(0.3D0),BESINT(0.9D0)
c      return
      NUM = 6
      N = 50
      ALFA(1) = 1.0D+00
      DO 10 K=2,NUM
         ALFA(K) = 1.0D+01*ALFA(K-1)
 10   CONTINUE
      WRITE(6,900)
      DO 30 NUU=1,11
         NU = NUU-1
C     COMPUTATION OF THE INTEGRALS
         CALL BESCHE(FUN,3.0D+01,ALFA,NUM,NU,N,RESULT,INFO)
C     COMPUTATION OF THE EXACT VALUE OF THE INTEGRALS AND OF
C     THE ERROR
         DO 20 I=1,NUM
            AB = DSQRT(4.0D+00+ALFA(I)**2)
            EXACT = ((AB-2.0D+00)/ALFA(I))**NU/AB
            DIF = EXACT-RESULT(I)
            WRITE(6,901) NU,ALFA(I),EXACT,RESULT(I),DIF,INFO(I),N
c            WRITE(6,901) NU,ALFA(I),RESULT(I),DIF,IER,N
 20      CONTINUE
 30   CONTINUE
 900  FORMAT(2HNU,6X,7HALFA(I),6X,14HEXACT INTEGRAL,8X,
     .     10HABS. ERROR,3X,4HINFO,4X,1HN/1X)
 901  FORMAT(I2,3X,F10.1,3X,D23.16,2X,D23.16,2X,D9.2,2X,I5,2X,I5)
      STOP
      END
      
      DOUBLE PRECISION FUNCTION FUN(X)
      DOUBLE PRECISION DEXP,X
      FUN=DEXP(-2.0D+00*X)
      RETURN
      END
