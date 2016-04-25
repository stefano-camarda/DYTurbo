*
* ..File e1nsn.f      (requires previous calls of BETAFCT, PNS0MOM and
*                      PNS1MOM)
*                  
*
* ..The subroutine  ENS1N  for the NLO hadronic non-singlet evolution 
*    kernels in N-space at a fixed numer of flavours  NF.  The routine 
*    returns the kernels  ENS(K)  for the evolution of the `+' (K=1),
*    and `-' (K=2,3) combinations of quark densities for a moment given 
*    via the counter  KN. The initial and final scales are specified by 
*    the respective values  ALPI  and  ALPF  for  a_s = alpha_s/(4 pi),
*    assuming a fixed ratio of mu_f and mu_r.
*
* ..Depending on  IMODE  given in the common-block  EVMOD  the kernels 
*    are determined in one of four schemes for solving the evolution 
*    equations at NLO.  Unlike the singlet and the NNLO non-singlet 
*    cases, also the iterated solutions can be written in a closed form
*    in the present case.  LOGFR = ln (mu_f^2/mu_r^2)  is taken from 
*    the common-block  FRRAT.
*
* ..The non-singlet splitting functions and the coefficients of the
*    beta function are taken from the common-blocks  PNS0,  PNS1  and 
*    BETA,  respectively.
*
* =====================================================================
*
*
       SUBROUTINE ENS1N (ENS, ALPI, ALPF, S, KN, NF)
*
       IMPLICIT DOUBLE COMPLEX (A - Z)
       INTEGER NDIM, NFMIN, NFMAX, KN, NF, IMODE
       PARAMETER (NDIM = 144, NFMIN = 3, NFMAX = 6)
       DOUBLE PRECISION BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX),
     1                  BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)
       DOUBLE PRECISION B0I, B10, ALPI, ALPF, S, LOGFR
       DIMENSION ENS(3)
       PARAMETER ( ONE = (1.D0, 0.D0) )
*
* ---------------------------------------------------------------------
*
* ..Input common-blocks 
*
       COMMON / EVMOD  / IMODE
       COMMON / FRRAT  / LOGFR
       COMMON / BETA   / BETA0, BETA1, BETA2, BETA3
       COMMON / PNS0   / P0NS (NDIM, NFMIN:NFMAX)
       COMMON / PNS1   / P1NS (NDIM, NFMIN:NFMAX, 3)
*
* ---------------------------------------------------------------------
*
* ..Some abbreviations 
*
       B0I  = 1./ BETA0(NF)
       B10  = BETA1(NF) * B0I 
*
* ..The LO evolution operator
*
       LNS  = EXP (S * P0NS(KN,NF) * B0I)
*
* ..The non-singlet NLO splitting-function combinations U1 = - R1
*   (including the contribution from mu_r unequal mu_f) 
*
       AUX = P0NS(KN,NF) * (LOGFR + B10* B0I)
       U1PLS = - P1NS(KN,NF,1) * B0I + AUX
       U1MIN = - P1NS(KN,NF,2) * B0I + AUX
*
* ---------------------------------------------------------------------
*
* ..The NLO evolution operator in four approximations
*
       IF ( IMODE .EQ. 1 ) THEN
*
* ..1) Iterated solution for the equation truncated for dq/dln(mu^2)
*
         AUX  = LOG( (1.+ B10* ALPF) / (1.+ B10* ALPI) ) / B10
         ENS(1) = LNS * EXP( AUX * U1PLS )
         ENS(2) = LNS * EXP( AUX * U1MIN ) 
*
* ..2) Iterated solution for the equation truncated for dq/da_s
*
       ELSE IF ( IMODE .EQ. 2 ) THEN
*
         ENS(1) = LNS * EXP ( (ALPF - ALPI) * U1PLS )
         ENS(2) = LNS * EXP ( (ALPF - ALPI) * U1MIN )
*
       ELSE IF ( IMODE .EQ. 3 ) THEN
*
* ..3) Truncated solution without expansion of U^(-1) 
*      (only possible in the non-singlet sector)
*
         ENS(1) = LNS * (ONE + ALPF * U1PLS) / (ONE + ALPI * U1PLS)
         ENS(2) = LNS * (ONE + ALPF * U1MIN) / (ONE + ALPI * U1MIN)
*
       ELSE
*
* ..4) fully truncated solution -- default.
*
         ENS(1) = LNS * ( ONE + (ALPF - ALPI) * U1PLS )
         ENS(2) = LNS * ( ONE + (ALPF - ALPI) * U1MIN )
*
       END IF
*
* ---------------------------------------------------------------------
*
* ..Third array element (for technical reasons only)
*
       ENS(3) = ENS(2)
*
* ---------------------------------------------------------------------
* 
       RETURN
       END
*
* =================================================================av==
