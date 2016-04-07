*
* ..File: ans2pmom.f     
*
*
* ..The polarized counterpart to  ANS2MOM  in  ans2mom.f.  
*   So far a dummy routine, only initializing a bunch of zeros. 
*
* =====================================================================
*
*
       SUBROUTINE ANS2PMOM 
*
       IMPLICIT DOUBLE COMPLEX (A - Z)
       INTEGER NMAX, NDIM, KN
       PARAMETER ( ZERO = (0.D0, 0.D0) )
       PARAMETER ( NDIM = 144 )
*
* ---------------------------------------------------------------------
*
* ..Input common-block
*
       COMMON / NNUSED / NMAX
*
* ..Output common-block 
*
       COMMON / ANS2   / A2NS (NDIM)
*
* ..Begin of the Mellin-N loop
*
       DO 1 KN = 1, NMAX
*
* ..Output to the array 
*
       A2NS(KN) = ZERO
*
* ---------------------------------------------------------------------
*
  1    CONTINUE
*
       RETURN
       END
*
* =================================================================av==
