*
* ..File: asg2pmom.f     
*
*
* ..The polarized counterpart to  ASG2MOM  in  asg2mom.f.  
*   So far a dummy routine, only initializing a bunch of zeros. 
*
* =====================================================================
*
*
       SUBROUTINE ASG2PMOM 
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
       COMMON / ASG2   / A2SG (NDIM, 2, 2)
*
* ..Begin of the Mellin-N loop
*
       DO 1 KN = 1, NMAX
*
* ..Output to the array 
*
       A2SG(KN,1,1) = ZERO
       A2SG(KN,1,2) = ZERO
       A2SG(KN,2,1) = ZERO
       A2SG(KN,2,2) = ZERO
*
* ---------------------------------------------------------------------
*
  1    CONTINUE
*
       RETURN
       END
*
* =================================================================av==
