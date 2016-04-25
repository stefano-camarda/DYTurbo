*
* ..File: evnvfn.f     (requires previous calls of INPLMOM, INITMOM and
*                       EVNFTHR)
*
*
* ..The subroutine  EVNVFN  for the N-space evolution of the parton
*    distributions for a variable number  NF = 3...6  of effectively 
*    massless flavours.  The results are returned by  PDFN  for the 
*    elements  NLOW -- NHIGH  of an external array of Mellin moments.
*    The notation for the second argument of  PDFN  reads
*
*      PDFN(KN,0) = g,  PDFN(KN,1) = u_v,  PDFN(KN,-1) = u + ubar, ...
*
*    for IPSTD = 0, and otherwise
*
*      PDFN(KN,0) = g,  PDFN(KN,1) = u,    PDFN(KN,-1) = ubar,  ...  .
*      
* ..The input moments are given by the common-blocks  PAINP and PAxTHR,
*    x = C, B, T,  for the initial scale and the thresholds for charm, 
*    bottom, top evolution.  The initial and final scales are specified 
*    by the respective values  ASI  and  ASF  of  a_s = alpha_s/(4 pi), 
*    assuming a fixed ratio  mu_f/mu_r.  The order ('n' in N^nLO) of 
*    the evolution is defined by  NPORD  in the common-block  ORDER. 
*     
* =====================================================================
*
*
       SUBROUTINE EVNVFN (PDFN, ASI, ASF, NF, NLOW, NHIGH, IPSTD)
*
       IMPLICIT DOUBLE COMPLEX (A - Z)
       INTEGER NF, NLOW, NHIGH, NPORD, IPSTD, NDIM, KN
       PARAMETER (NDIM = 144)
       DOUBLE PRECISION ASI, ASF, S, THRD
       PARAMETER ( THRD = 1.D0/3.D0 )
       DIMENSION PDFN (NDIM, -6:6), ESG(2,2), ENS(3)
*
* ---------------------------------------------------------------------
* 
* ..Input common blocks
* 
       COMMON / ORDER  / NPORD
       COMMON / PAINP  / VAI (NDIM), M3I (NDIM), M8I(NDIM),
     1                   SGI (NDIM), P3I (NDIM), P8I(NDIM), GLI (NDIM)
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
* ..Evolution distance
*
       S = LOG (ASI / ASF)
*
* ..Begin of the Mellin-N loop
*
       DO 1 KN = NLOW, NHIGH
*
* ..Non-singlet and singlet evolution kernels up to NNLO
*
       IF (NPORD .EQ. 0) THEN
         CALL ENSG0N (ENS, ESG, ASI, ASF, S, KN, NF)
       ELSE IF (NPORD .EQ. 1) THEN
         CALL ENS1N (ENS, ASI, ASF, S, KN, NF)
         CALL ESG1N (ESG, ASI, ASF, S, KN, NF)
       ELSE
         CALL ENS2N (ENS, ASI, ASF, S, KN, NF, 1, 3)
         CALL ESG2N (ESG, ASI, ASF, S, KN, NF)
       END IF
*
* ---------------------------------------------------------------------
*
       IF (NF .EQ. 3) THEN
*
* ..Three-flavour evolution
*
         M3N = ENS(2) * M3I(KN)
         M8N = ENS(2) * M8I(KN)
         VAN = ENS(3) * VAI(KN)
*
         P3N = ENS(1) * P3I(KN)
         P8N = ENS(1) * P8I(KN)
         SGN = ESG(1,1) * SGI(KN) + ESG(1,2) * GLI(KN) 
         GLN = ESG(2,1) * SGI(KN) + ESG(2,2) * GLI(KN)
*
         M15N = VAN
         P15N = SGN
         M24N = VAN
         P24N = SGN
         M35N = VAN
         P35N = SGN
*
* ---------------------------------------------------------------------
*         
       ELSE IF (NF .EQ. 4) THEN
*
* ..Four-flavour evolution
*
         M3N = ENS(2) * M3C(KN)
         M8N = ENS(2) * M8C(KN)
         VAN = ENS(3) * VAC(KN)
*
         P3N = ENS(1) * P3C(KN)
         P8N = ENS(1) * P8C(KN)
         SGN = ESG(1,1) * SGC(KN) + ESG(1,2) * GLC(KN) 
         GLN = ESG(2,1) * SGC(KN) + ESG(2,2) * GLC(KN)
*
         M15N= ENS(2) * M15C(KN)
         P15N= ENS(1) * P15C(KN)
         M24N = VAN
         P24N = SGN
         M35N = VAN
         P35N = SGN
*
* ---------------------------------------------------------------------
*         
       ELSE IF (NF .EQ. 5) THEN
*
* ..Five-flavour evolution
*
         M3N = ENS(2) * M3B(KN)
         M8N = ENS(2) * M8B(KN)
         VAN = ENS(3) * VAB(KN)
*
         P3N = ENS(1) * P3B(KN)
         P8N = ENS(1) * P8B(KN)
         SGN = ESG(1,1) * SGB(KN) + ESG(1,2) * GLB(KN) 
         GLN = ESG(2,1) * SGB(KN) + ESG(2,2) * GLB(KN)
*
         M15N = ENS(2) * M15B(KN)
         P15N = ENS(1) * P15B(KN)
         M24N = ENS(2) * M24B(KN)
         P24N = ENS(1) * P24B(KN)
         M35N = VAN
         P35N = SGN
*

* ---------------------------------------------------------------------
*
       ELSE 
*
* ..Six-flavour evolution
*
         M3N = ENS(2) * M3T(KN)
         M8N = ENS(2) * M8T(KN)
         VAN = ENS(3) * VAT(KN)
*
         P3N = ENS(1) * P3T(KN)
         P8N = ENS(1) * P8T(KN)
         SGN = ESG(1,1) * SGT(KN) + ESG(1,2) * GLT(KN)
         GLN = ESG(2,1) * SGT(KN) + ESG(2,2) * GLT(KN)
*
         M15N = ENS(2) * M15T(KN)
         P15N = ENS(1) * P15T(KN)
         M24N = ENS(2) * M24T(KN)
         P24N = ENS(1) * P24T(KN)
         M35N = ENS(2) * M35T(KN)
         P35N = ENS(1) * P35T(KN)
*
       END IF
*
* ---------------------------------------------------------------------
*
* ..Flavour decomposition of the `minus' sector: q-qbar distributions
*
       TV = (VAN - M35N) * 0.5 * THRD
       BV = TV + (M35N - M24N) * 0.2  
       CV = BV + (M24N - M15N) * 0.25 
       SV = CV + (M15N - M8N) * THRD 
       DV = SV + (M8N - M3N) * 0.5 
       UV = SV + (M8N + M3N) * 0.5 
*
* ..Flavour decomposition of the `plus' sector:  q+qbar distributions
*
       TP = (SGN - P35N) * 0.5 * THRD
       BP = TP + (P35N - P24N) * 0.2
       CP = BP + (P24N - P15N) * 0.25
       SP = CP + (P15N - P8N) * THRD
       DP = SP + (P8N - P3N) * 0.5
       UP = SP + (P8N + P3N) * 0.5
*
* ---------------------------------------------------------------------
*
* ..Output to the array
*
       PDFN(KN,0)  = GLN
       IF (IPSTD .EQ. 0) THEN
         PDFN(KN,1)  = UV
         PDFN(KN,2)  = DV
         PDFN(KN,3)  = SV
         PDFN(KN,4)  = CV
         PDFN(KN,5)  = BV
         PDFN(KN,6)  = TV
         PDFN(KN,-1) = UP
         PDFN(KN,-2) = DP
         PDFN(KN,-3) = SP
         PDFN(KN,-4) = CP
         PDFN(KN,-5) = BP
         PDFN(KN,-6) = TP
       ELSE
         PDFN(KN,1)  = 0.5 * (UV + UP)
         PDFN(KN,2)  = 0.5 * (DV + DP)
         PDFN(KN,3)  = 0.5 * (SV + SP)
         PDFN(KN,4)  = 0.5 * (CV + CP)
         PDFN(KN,5)  = 0.5 * (BV + BP)
         PDFN(KN,6)  = 0.5 * (TV + TP)
         PDFN(KN,-1) = 0.5 * (UP - UV)
         PDFN(KN,-2) = 0.5 * (DP - DV)
         PDFN(KN,-3) = 0.5 * (SP - SV)
         PDFN(KN,-4) = 0.5 * (CP - CV)
         PDFN(KN,-5) = 0.5 * (BP - BV)
         PDFN(KN,-6) = 0.5 * (TP - TV)
       END IF
*
* ---------------------------------------------------------------------
*
  1    CONTINUE 
*
       RETURN
       END
*
* =================================================================av==
