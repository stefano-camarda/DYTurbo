*
* ..File: evnfthr.f   (requires previous calls of INPLMOM and INITMOM)
*
*
* ..The subroutine  EVNFTHR  for the evolution of  a_s = alpha_s/(4 pi)
*    and the N-space parton densities from a three-flavour initial 
*    scale to the four- to six-flavour thresholds (identified with 
*    the squares of the corresponding quark masses).  The results are 
*    written to the common-blocks  ASFTHR  and  PAxTHR, x = C, B, T.
*    For  MT2  or  MB2 and MT2  above 10^10 (GeV^2) the corresponding 
*    parts of the calculations are skipped.
*
* ..The three-flavour input moments are given by the common-block PAINP
*    for an external array of Mellin moments.  The input scale  M20 = 
*    mu_(f,0)^2  and the corresponding value  AS0  of a_s  are provided
*    by  ASINP.  The fixed scale logarithm  LOGFR = ln (mu_f^2/mu_r^2)
*    is specified in  FRRAT.  The order ('n' in N^nLO) of the evolution
*    is defined by  NPORD  in the common-block  ORDER.   
*
* ..The operator matrix elements  A2NS  and  A2SG  required for the 
*    N_f matching of the parton densities at NNLO are taken from the
*    common-blocks  ANS2  and  ASG2,  respectively.  The corresponding
*    alpha_s matching is done by the function  ASNF1  in asmatch.f
*
* =====================================================================
*
*
       SUBROUTINE EVNASTHR (MC2, MB2, MT2)
*
       IMPLICIT DOUBLE COMPLEX (A - Z)
       INTEGER NPORD, NMAX, KN
       include 'dimensions.f'
       DOUBLE PRECISION MC2, MB2, MT2, M20, M2C, M2B, M2T, R20, R2C, 
     1                  R2B, R2T, AS, ASNF1, AS0, ASC, ASB, AST, ASC2, 
     2                  ASB2, AST2, ASC3, ASB4, AST5, LOGFR, SC, SB, ST
*
* ---------------------------------------------------------------------
* 
* ..Input common blocks
*  
       COMMON / NNUSED / NMAX
       COMMON / ORDER  / NPORD
       COMMON / ASINP  / AS0, M20
       COMMON / FRRAT  / LOGFR
*
* ..Output common blocks
*
       COMMON / ASFTHR / ASC, M2C, ASB, M2B, AST, M2T
*
* ---------------------------------------------------------------------
*
* ..Coupling constants at and evolution distances to/between thresholds
* 
       R20 = M20 * EXP(-LOGFR)
*
* ..Charm
*
       M2C  = MC2
       R2C  = M2C * R20/M20
       ASC3 = AS (R2C, R20, AS0, 3)
       SC   = LOG (AS0 / ASC3)
       ASC  = ASNF1 (ASC3, -LOGFR, 3)
*
* ..Bottom 
*
       M2B  = MB2
       R2B  = M2B * R20/M20
       ASB4 = AS (R2B, R2C, ASC, 4)
       SB   = LOG (ASC / ASB4)
       ASB  = ASNF1 (ASB4, -LOGFR, 4)
*
* ..Top
*
       M2T  = MT2
       R2T  = M2T * R20/M20
       AST5 = AS (R2T, R2B, ASB, 5)
       ST   = LOG (ASB / AST5)
       AST  = ASNF1 (AST5, -LOGFR, 5)
*
* ..No non-trivial terms in the pdf matching at LO and NLO:
*
       IF (NPORD .LT. 2) THEN
         ASC2 = 0.D0
         ASB2 = 0.D0
         AST2 = 0.D0
       ELSE
         ASC2 = ASC * ASC
         ASB2 = ASB * ASB
         AST2 = AST * AST
       END IF
       RETURN
       END
