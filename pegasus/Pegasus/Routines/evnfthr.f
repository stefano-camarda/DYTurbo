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
       SUBROUTINE EVNFTHR (MC2, MB2, MT2)
*
       IMPLICIT DOUBLE COMPLEX (A - Z)
       INTEGER NPORD, NDIM, NMAX, KN
       PARAMETER (NDIM = 144)
       DOUBLE PRECISION MC2, MB2, MT2, M20, M2C, M2B, M2T, R20, R2C, 
     1                  R2B, R2T, AS, ASNF1, AS0, ASC, ASB, AST, ASC2, 
     2                  ASB2, AST2, ASC3, ASB4, AST5, LOGFR, SC, SB, ST
       DIMENSION ESG(2,2), ENS(3)
*
* ---------------------------------------------------------------------
* 
* ..Input common blocks
*  
       COMMON / NNUSED / NMAX
       COMMON / ORDER  / NPORD
       COMMON / ASINP  / AS0, M20
       COMMON / FRRAT  / LOGFR
       COMMON / PAINP  / VAI (NDIM), M3I (NDIM), M8I (NDIM), 
     1                   SGI (NDIM), P3I (NDIM), P8I (NDIM), GLI (NDIM)
       COMMON / ANS2   / A2NS (NDIM)
       COMMON / ASG2   / A2SG (NDIM,2,2)
*
* ..Output common blocks
*
       COMMON / ASFTHR / ASC, M2C, ASB, M2B, AST, M2T
       COMMON / PACTHR / VAC (NDIM), M3C (NDIM), M8C (NDIM), M15C(NDIM),
     1                   SGC (NDIM), P3C (NDIM), P8C (NDIM), P15C(NDIM),
     2                   GLC (NDIM) 
       COMMON / PABTHR / VAB (NDIM), M3B (NDIM), M8B (NDIM), M15B(NDIM),
     1                   SGB (NDIM), P3B (NDIM), P8B (NDIM), P15B(NDIM),
     2                   M24B(NDIM), P24B(NDIM), GLB (NDIM)
       COMMON / PATTHR / VAT (NDIM), M3T (NDIM), M8T (NDIM), M15T(NDIM),
     1                   SGT (NDIM), P3T (NDIM), P8T (NDIM), P15T(NDIM),
     2                   M24T(NDIM), P24T(NDIM), M35T(NDIM), P35T(NDIM),
     3                   GLT (NDIM)
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
*
* ---------------------------------------------------------------------
*
* ..Begin of the Mellin-N loop
*
       DO 1 KN = 1, NMAX
*
* ..The kernels for the evolution from M20 to M2C
*
       IF (NPORD .EQ. 0) THEN
         CALL ENSG0N (ENS, ESG, AS0, ASC3, SC, KN, 3)
       ELSE IF (NPORD .EQ. 1) THEN
         CALL ENS1N (ENS, AS0, ASC3, SC, KN, 3)
         CALL ESG1N (ESG, AS0, ASC3, SC, KN, 3)
       ELSE
         CALL ENS2N (ENS, AS0, ASC3, SC, KN, 3, 1, 3)
         CALL ESG2N (ESG, AS0, ASC3, SC, KN, 3)
       END IF
*
* ..The N_f=3 parton distributions at the four-flavour threshold
*
       M3C3 = ENS(2) * M3I(KN)
       M8C3 = ENS(2) * M8I(KN)
       VAC3 = ENS(3) * VAI(KN)
*
       P3C3 = ENS(1) * P3I(KN)
       P8C3 = ENS(1) * P8I(KN)
       SGC3 = ESG(1,1) * SGI(KN) + ESG(1,2) * GLI(KN) 
       GLC3 = ESG(2,1) * SGI(KN) + ESG(2,2) * GLI(KN)
*
* ---------------------------------------------------------------------
*
* ..The N_f=4 parton distributions at the four-flavour threshold
*
       NSTHR = 1.+ ASC2 * A2NS(KN)
       HQTHR = ASC2 * (A2SG(KN,1,1) * SGC3 + A2SG(KN,1,2) * GLC3)    
*
       M3C(KN)  = M3C3 * NSTHR
       M8C(KN)  = M8C3 * NSTHR
       VAC(KN)  = VAC3 * NSTHR
       M15C(KN) = VAC(KN)
*
       P3C(KN)  = P3C3 * NSTHR 
       P8C(KN)  = P8C3 * NSTHR 
*
       SGC(KN)  = SGC3 * NSTHR + HQTHR 
       P15C(KN) = SGC(KN) - 4.* HQTHR
       GLC(KN) = GLC3 + ASC2* (A2SG(KN,2,1)* SGC3 + A2SG(KN,2,2)* GLC3)
*    
* ---------------------------------------------------------------------
*
* ..The kernels for the evolution from M2C to M2B
*
       IF (MB2 .GT. 1.D10) GO TO 1
       IF (NPORD .EQ. 0) THEN
         CALL ENSG0N (ENS, ESG, ASC, ASB4, SB, KN, 4)
       ELSE IF (NPORD .EQ. 1) THEN
         CALL ENS1N (ENS, ASC, ASB4, SB, KN, 4)
         CALL ESG1N (ESG, ASC, ASB4, SB, KN, 4)
       ELSE
         CALL ENS2N (ENS, ASC, ASB4, SB, KN, 4, 1, 3)
         CALL ESG2N (ESG, ASC, ASB4, SB, KN, 4)
       END IF
*
* ..The N_f=4 parton distributions at the five-flavour threshold
*
       M3B4  = ENS(2) * M3C(KN)
       M8B4  = ENS(2) * M8C(KN)
       M15B4 = ENS(2) * M15C(KN)
       VAB4  = ENS(3) * VAC(KN)

       P3B4  = ENS(1) * P3C(KN)
       P8B4  = ENS(1) * P8C(KN)
       P15B4 = ENS(1) * P15C(KN)
       SGB4  = ESG(1,1) * SGC(KN) + ESG(1,2) * GLC(KN) 
       GLB4  = ESG(2,1) * SGC(KN) + ESG(2,2) * GLC(KN)
*
* ---------------------------------------------------------------------
*
* ..The N_f=5 parton distributions at the five-flavour threshold
*
       NSTHR = 1.+ ASB2 * A2NS(KN)
       HQTHR = ASB2 * (A2SG(KN,1,1) * SGB4 + A2SG(KN,1,2) * GLB4)    
*
       M3B(KN)  = M3B4 * NSTHR
       M8B(KN)  = M8B4 * NSTHR
       M15B(KN) = M15B4 * NSTHR
       VAB(KN)  = VAB4 * NSTHR
       M24B(KN) = VAB(KN)
*
       P3B(KN)  = P3B4 * NSTHR 
       P8B(KN)  = P8B4 * NSTHR 
       P15B(KN) = P15B4 * NSTHR 
*
       SGB(KN)  = SGB4 * NSTHR + HQTHR 
       P24B(KN) = SGB(KN) - 5.* HQTHR
       GLB(KN) = GLB4 + ASB2* (A2SG(KN,2,1)* SGB4 + A2SG(KN,2,2)* GLB4)
*    
* ---------------------------------------------------------------------
*
* ..The kernels for the evolution from M2B to M2T
*
       IF (MT2 .GT. 1.D10) GO TO 1
       IF (NPORD .EQ. 0) THEN
         CALL ENSG0N (ENS, ESG, ASB, AST5, ST, KN, 5)
       ELSE IF (NPORD .EQ. 1) THEN
         CALL ENS1N (ENS, ASB, AST5, ST, KN, 5)
         CALL ESG1N (ESG, ASB, AST5, ST, KN, 5)
       ELSE
         CALL ENS2N (ENS, ASB, AST5, ST, KN, 5, 1, 3)
         CALL ESG2N (ESG, ASB, AST5, ST, KN, 5)
       END IF
*
* ..The N_f=5 parton distributions at the six-flavour threshold
*
       M3T5  = ENS(2) * M3B(KN)
       M8T5  = ENS(2) * M8B(KN)
       M15T5 = ENS(2) * M15B(KN)
       M24T5 = ENS(2) * M24B(KN)
       VAT5  = ENS(3) * VAB(KN)

       P3T5  = ENS(1) * P3B(KN)
       P8T5  = ENS(1) * P8B(KN)
       P15T5 = ENS(1) * P15B(KN)
       P24T5 = ENS(1) * P24B(KN)
       SGT5  = ESG(1,1) * SGB(KN) + ESG(1,2) * GLB(KN) 
       GLT5  = ESG(2,1) * SGB(KN) + ESG(2,2) * GLB(KN)
*
* ---------------------------------------------------------------------
*
* ..The N_f=6 parton distributions at the six-flavour threshold
*
       NSTHR = 1.+ AST2 * A2NS(KN)
       HQTHR = AST2 * (A2SG(KN,1,1) * SGT5 + A2SG(KN,1,2) * GLT5)    
*
       M3T(KN)  = M3T5 * NSTHR
       M8T(KN)  = M8T5 * NSTHR
       M15T(KN) = M15T5 * NSTHR
       M24T(KN) = M24T5 * NSTHR
       VAT(KN)  = VAT5 * NSTHR
       M35T(KN) = VAT(KN)
*
       P3T(KN)  = P3T5 * NSTHR 
       P8T(KN)  = P8T5 * NSTHR 
       P15T(KN) = P15T5 * NSTHR 
       P24T(KN) = P24T5 * NSTHR 
*
       SGT(KN)  = SGT5 * NSTHR + HQTHR 
       P35T(KN) = SGT(KN) - 6.* HQTHR
       GLT(KN) = GLT5 + AST2* (A2SG(KN,2,1)* SGT5 + A2SG(KN,2,2)* GLT5)
*    
* ---------------------------------------------------------------------
*
  1    CONTINUE 
*
       RETURN
       END
*
* =================================================================av==
