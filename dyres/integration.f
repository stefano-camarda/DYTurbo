      FUNCTION ADPINT (F, A, B, AERR, RERR, ERREST, IER)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      integer oldint
     
c.....Integral of F(X) from A to B, with error
c.....less than ABS(AERR) + ABS(RERR*INTEGRAL)
c.....Best estimate of error returned in ERREST.
c.....Error code is IER: zero if OK, non-zero if in trouble.
     
      EXTERNAL F
      PARAMETER (MAXINT = 500)
     
c     Work space:

      COMMON / ADPWRK / U(MAXINT), V(MAXINT), FU(MAXINT),
     >     FV(MAXINT), FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), NUMINT
      SAVE / ADPWRK /
      
      IER = 0
      NUMINT = 5
      DX = (B-A)/ NUMINT
      DO 1  I = 1, NUMINT
         IF (I .EQ. 1)  THEN
            U(I) = A
            FU(I) = F(U(I))
         ELSE
            U(I) = V(I-1)
            FU(I) = FV(I-1)
         ENDIF
         IF (I .EQ. NUMINT) THEN
            V(I) = B
         ELSE
            V(I) = A + DX * I
         ENDIF
         FV(I) = F(V(I))
         CALL ADPCAL(F,I)
 1    CONTINUE
      
 2    CONTINUE
      
c.....Error estimate:
      
      ADPINT = 0.
      ERREST = 0.
      DO 3  I = 1, NUMINT
         ADPINT = ADPINT + RESULT(I)
         ERREST = ERREST + ERR(I)
 3    CONTINUE
      TARGET = ABS(AERR) + ABS(RERR * ADPINT)
      IF (ERREST .GT. TARGET)  THEN
         OLDINT = NUMINT
         DO 4 I = 1, OLDINT
            IF (ERR(I)*2*OLDINT .GT. TARGET)  CALL ADPSPL(F,I,IER)
 4       CONTINUE
         IF (IER .EQ. 0)  GOTO 2
      ENDIF
      RETURN
      END
      
      FUNCTION INTUSE ()
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

c.....Return number of intervals used last call to ADPINT

      PARAMETER (MAXINT = 500)
      COMMON / ADPWRK / U(MAXINT), V(MAXINT), FU(MAXINT),
     >     FV(MAXINT), FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), NUMINT
      INTUSE = NUMINT
      RETURN
      END
      
      SUBROUTINE ADPCAL (F,I)

c.....Fill in details of interval I given endpoints
      
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL F
      PARAMETER (MAXINT = 500)
      COMMON / ADPWRK / U(MAXINT), V(MAXINT), FU(MAXINT),
     >     FV(MAXINT), FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), NUMINT
      
      FW(I) = F( (U(I) + V(I)) /2.)
      DX = V(I) - U(I)
      RESULT(I) = DX * (FU(I) + 4. * FW(I) + FV(I)) / 6.
      ERR(I) = ABS(DX * (FU(I) - 2. * FW(I) + FV(I)) / 12.)
      RETURN
      END
      
      SUBROUTINE ADPSPL (F, I, IER)

c.....Split interval I

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL F
      PARAMETER (MAXINT = 500)
      COMMON / ADPWRK / U(MAXINT), V(MAXINT), FU(MAXINT),
     >     FV(MAXINT), FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), NUMINT
      
      IF (NUMINT .GE. MAXINT)  THEN
         IER = 1
         RETURN
      ENDIF
      NUMINT = NUMINT + 1
      V(NUMINT) = V(I)
      U(NUMINT) = (U(I) + V(I)) / 2.
      V(I) = U(NUMINT)
      FV(NUMINT) = FV(I)
      FU(NUMINT) = FW(I)
      FV(I) = FW(I)
      CALL ADPCAL (F, I)
      CALL ADPCAL (F, NUMINT)
      RETURN
      END
