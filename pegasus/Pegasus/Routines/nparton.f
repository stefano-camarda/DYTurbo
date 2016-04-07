*
* ..File nparton.f
*
*
* ..The subroutine  NPARTON  for one (in general complex) Mellin moment 
*    N  of the parton densities evolved to the scale  M2 = mu_f^2  by 
*    EVNFFN  (for IVFNS = 0)  or  EVNVFN.  
*
* ..The results are returned as  f_i(N,mu_f^2)  by the array  PDFN.
*    See the header of evnffn.f for the notation for f_i depending
*    on the input parameter  IPSTD.  For  IPOL = 0  we are dealing
*    with the unpolarized, otherwise with the polarized case.
*
* ..Note that, unlike its x-space counterpart  XPARTON  the present
*    routine includes calls of  INITPOL (.., EPAR)  and  INITINP (IPAR) 
*    or  INITPINP (IPAR),  i.e., for the N-space evolution this routine 
*    is the only one the user needs to call.  Input and initialization
*    parameters are specified as usual depending on  EPAR  and  IPAR.
* 
* =====================================================================
*
*
       SUBROUTINE NPARTON (PDFN, ASOUT, N, M2, IPSTD, IPOL, EPAR, IPAR)

       IMPLICIT DOUBLE PRECISION (A - Z)
       INTEGER NDIM, NF, NFF, IPSTD, IVFNS, IPOL, EPAR, IPAR, K
       PARAMETER ( NDIM = 144 )
       DOUBLE COMPLEX N, NPDF(NDIM,-6:6), PDFN(-6:6)
*
* ---------------------------------------------------------------------
*
* ..Input common blocks 
* 
       COMMON / NFFIX  / NFF
       COMMON / VARFLV / IVFNS 
       COMMON / FRRAT  / LOGFR
       COMMON / ASINP  / AS0, M20
       COMMON / ASFTHR / ASC, M2C, ASB, M2B, AST, M2T
*
* ---------------------------------------------------------------------
*
* ..Compute the splitting functions etc for the specific value of N
*
       CALL INITMOM (N, IPOL, EPAR)
*
* ..Initialize the input parameters and initial moments
*
       IF ( IPOL .EQ. 0 ) THEN
         CALL INITINP (IPAR)
       ELSE
         CALL INITPINP (IPAR)
       END IF
*
* ---------------------------------------------------------------------
*
* ..Values of a_s and N_f for the call of EVNFFN or EVNVFN below
*
       R2  = M2 * EXP(-LOGFR)
       IF (IVFNS .EQ. 0) THEN
*
*   Fixed number of flavours
*
         NF  = NFF
         R20 = M20 * R2/M2
         ASI = AS0
         ASF = AS (R2, R20, AS0, NF)
*
       ELSE
*
* ..Variable number of flavours
*
         IF (M2 .GT. M2T) THEN
           NF = 6
           R2T = M2T * R2/M2
           ASI = AST
           ASF = AS (R2, R2T, AST, NF)
*
         ELSE IF (M2 .GT. M2B) THEN
           NF = 5
           R2B = M2B * R2/M2
           ASI = ASB
           ASF = AS (R2, R2B, ASB, NF)
*
         ELSE IF (M2 .GT. M2C) THEN
           NF = 4
           R2C = M2C * R2/M2
           ASI = ASC
           ASF = AS (R2, R2C, ASC, NF)
*
         ELSE
           NF = 3
           R20 = M20 * R2/M2
           ASI = AS0
           ASF = AS (R2, R20, AS0, NF)
*       
         END IF
*
       END IF
*
* ..Output assignment of a_s
*
       ASOUT = ASF
*
* ---------------------------------------------------------------------
*
* ..Calculation of the moments of the parton densities and output
*
       IF (IVFNS .EQ. 0) THEN
*
          CALL EVNFFN (PDFN, ASI, ASF, NF, 1, 1, IPSTD)
       ELSE
*
          CALL EVNVFN (NPDF, ASI, ASF, NF, 1, 1, IPSTD)
       END IF
*
       DO 1 K = -6, 6
         PDFN(K) = NPDF(1,K)
  1    CONTINUE
*
* ---------------------------------------------------------------------
*
       RETURN                                                           
       END                                                             
*
* =================================================================av==
