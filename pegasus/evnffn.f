*
* ..File: evnffn.f    (requires previous calls of INPLMOM and INITMOM)
*
*
* ..The subroutine  EVNFFN  for the N-space evolution of the parton
*    distributions for a fixed number  NF = 3 ... 5  of effectively 
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
* ..The light-parton input moments are given by the common-block  PAINP,
*    for  NF = 4 [5] in addition  c = cbar = 0  [and  b = bbar = 0]  is 
*    assumed at the input scale.  This scale and the output scale are 
*    specified by the values  ASI  and  ASF  of  a_s = alpha_s/(4 pi), 
*    assuming a fixed ratio  mu_f/mu_r.  The order ('n' in N^nLO) of 
*    the evolution is defined by  NPORD  in the common-block  ORDER. 
*     
* =====================================================================
*
*
       SUBROUTINE EVNFFN (PDFN, ASI, ASF, NF, NLOW, NHIGH, IPSTD)
*
       IMPLICIT DOUBLE COMPLEX (A - Z)
       INTEGER NF, NLOW, NHIGH, NPORD, IPSTD, NDIM, KN
       PARAMETER (NDIM = 144)
       DOUBLE PRECISION ASI, ASF, S, THRD, SXTH
       PARAMETER ( THRD = 1.D0/3.D0, SXTH = THRD/2 )
       DIMENSION PDFN (NDIM, -6:6), ESG(2,2), ENS(3)
       PARAMETER ( ZERO = (0.D0, 0.D0) ) 
*
* ---------------------------------------------------------------------
* 
* ..Input common blocks
* 
       COMMON / PAINP  / VAI (NDIM), M3I (NDIM), M8I(NDIM), 
     1                   SGI (NDIM), P3I (NDIM), P8I(NDIM), GLI (NDIM)
       COMMON / ORDER  / NPORD
*
* ---------------------------------------------------------------------
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
         CALL ENS2N (ENS, ASI, ASF, S, KN, NF, 1, 3 )
         CALL ESG2N (ESG, ASI, ASF, S, KN, NF)
       END IF
*
* ---------------------------------------------------------------------
*
* ..Evolution of the singlet and light-flavour non-singlet combinations
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
* ..The additional non-singlet combinations involving charm and bottom
*
       IF (NF .EQ. 3) THEN
         M15N = VAN
         P15N = SGN
         M24N = VAN
         P24N = SGN
       ELSE IF (NF .EQ. 4) THEN
         M15N = ENS(2) * VAI(KN)
         P15N = ENS(1) * SGI(KN)
         M24N = VAN
         P24N = SGN
       ELSE 
         M15N = ENS(2) * VAI(KN)
         P15N = ENS(1) * SGI(KN)
         M24N = M15N
         P24N = P15N
       END IF
*
* ---------------------------------------------------------------------
*
* ..Flavour decomposition of the `minus' sector:  valence distributions
*
       BV = (VAN - M24N) * 0.2  
       CV = BV + (M24N - M15N) * 0.25 
       SV = CV + (M15N - M8N) * THRD 
       DV = SV + (M8N - M3N) * 0.5 
       UV = SV + (M8N + M3N) * 0.5 
*
* ..Flavour decomposition of the `plus' sector: antiquark distributions
*
       BP = (SGN - P24N) * 0.2 
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
       PDFN(KN,6)  = ZERO
       PDFN(KN,-6) = ZERO
*
       IF (IPSTD .EQ. 0) THEN
         PDFN(KN,1)  = UV
         PDFN(KN,2)  = DV
         PDFN(KN,3)  = SV
         PDFN(KN,4)  = CV
         PDFN(KN,5)  = BV
         PDFN(KN,-1) = UP 
         PDFN(KN,-2) = DP 
         PDFN(KN,-3) = SP 
         PDFN(KN,-4) = CP 
         PDFN(KN,-5) = BP 
       ELSE
         PDFN(KN,1)  = 0.5 * (UV + UP)
         PDFN(KN,2)  = 0.5 * (DV + DP)
         PDFN(KN,3)  = 0.5 * (SV + SP)
         PDFN(KN,4)  = 0.5 * (CV + CP)
         PDFN(KN,5)  = 0.5 * (BV + BP)
         PDFN(KN,-1) = 0.5 * (UP - UV) 
         PDFN(KN,-2) = 0.5 * (DP - DV)
         PDFN(KN,-3) = 0.5 * (SP - SV)
         PDFN(KN,-4) = 0.5 * (CP - CV)
         PDFN(KN,-5) = 0.5 * (BP - BV)
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
