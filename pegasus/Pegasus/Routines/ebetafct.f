*
* ..File ebetafct.f
*
*
* ..The complex Euler Beta function,  EBETA (Z1, Z2),  calculated from 
*    the asymtotic expansion of the logarithm of the Gamma function at 
*    Re(Z1) > 5 and Re(Z2) > 5  together with the functional equation. 
*    For  Re(Z1) < - 10  the argument  Z1  is first inverted, i.e., the
*    asymtotic expansion is invoked for  Z1' = 1 - Z1 - Z2.
*
* =====================================================================
*
*
       FUNCTION EBETA (Z1, Z2)
*
       IMPLICIT DOUBLE COMPLEX (A - Z)
       INTEGER K
       PARAMETER ( ONE = (1.D0, 0.D0), HLF = (0.5D0, 0.D0), 
     1             IPI = (0.D0, 3.14159 26535 89793 D0) )
*
* ..The asymtotic expansion for ln Gamma(x) (with X2 = 1/x^2)
*
       DATA C1 / (9.18938 53320 46727 D-1, 0.D0) /
       DATA C2 / (8.33333 33333 33333 D-2, 0.D0) /
       DATA C3 / (3.33333 33333 33333 D-2, 0.D0) /
       DATA C4 / (2.85714 28514 28514 D-1, 0.D0) /
       DATA C5 / (0.75 D0, 0.D0) /
*
       LNGAM (X,X2) = (X - HLF) * LOG (X) - X + C1 + C2 / X
     1              * (ONE - C3* X2 * (ONE - C4* X2 * (ONE - C5* X2)))
*
* ---------------------------------------------------------------------
*
* ..Inversion of the first argument
*
       ZZ1 = Z1
       IF (DBLE (Z1) .LT. -10.D0) THEN
         ZZ1 = ONE - Z1 - Z2
       END IF
*
* ---------------------------------------------------------------------
*
* ..Shift of the arguments using the functional equation
*
       NUM = ONE
       DEN = ONE
* 
       IF ( DBLE (ZZ1) .LT. 5.D0) THEN
       DO 1 K = 1, 5-INT(ZZ1)
         NUM = NUM * (ZZ1+Z2)
         DEN = DEN * ZZ1
         ZZ1 = ZZ1 + ONE
  1    CONTINUE
       END IF
*
       ZZ2 = Z2
       IF ( DBLE (Z2) .LT. 5.D0) THEN
       DO 2 K = 1, 5-INT(ZZ2)
         NUM = NUM * (ZZ1+ZZ2)
         DEN = DEN * ZZ2
         ZZ2 = ZZ2 + ONE 
  2    CONTINUE
       END IF
*
* ---------------------------------------------------------------------
*
* ..Use of the asymtotic expansions at the shifted arguments
*
       LG1  = LNGAM ( ZZ1, ONE/(ZZ1*ZZ1) )
       LG2  = LNGAM ( ZZ2, ONE/(ZZ2*ZZ2) )
       LG12 = LNGAM ( ZZ1+ZZ2, ONE/((ZZ1+ZZ2)*(ZZ1+ZZ2)) )
*
       EBETA = NUM / DEN * EXP ( LG1 + LG2 - LG12)
*
* ---------------------------------------------------------------------
*
* ..The additional factor arising from the inversion of Z1
*
       IF (DBLE (Z1) .LT. -10.D0) THEN
*
         IF (DIMAG(Z1) .GT. 0.D0) THEN
           AUX = EXP (IPI*Z1)
           EBETA = (AUX - EXP(IPI*Z2)) / (AUX - ONE) * EBETA
         ELSE
           AUX = EXP (-IPI*Z1)
           EBETA = (EXP(-IPI*Z2) - AUX) / (ONE - AUX) * EBETA
         END IF
*
       END IF
*
* ---------------------------------------------------------------------
*
       RETURN
       END
*
* =================================================================av==
