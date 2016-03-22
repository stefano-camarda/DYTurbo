*
* ..File ensg0n.f     (requires previous calls of BETAFCT, PNS0MOM and 
*                      LSGMOM) 
*
* ..The subroutine  ENSG0N  for the LO hadronic non-singlet and singlet 
*    evolution operators in N-space at a fixed number of flavours  NF.  
*    The kernels  ENS  and  ESG(K1,K2),  K1, K2 = 1, 2  with  1 = q, 
*    2 = g,  are returned for a moment N specified via the counter  KN. 
*    The initial and final scales are specified by the respective 
*    values  ASI  and  ASF  of  a_s = alpha_s/(4 pi),  assuming a 
*    fixed ratio of mu_r and mu_f.
*
* ..The non-singlet splitting functions, the eigenvalue decomposition 
*    of the singlet splitting-function matrix and the coefficients of 
*    the beta function are taken from the common-blocks   PNS0,  LSG, 
*    and  BETA,  respectively.
*
* =====================================================================
*
*
       SUBROUTINE ENSG0N (ENS, ESG, ASI, ASF, S, KN, NF)
*
       IMPLICIT DOUBLE COMPLEX (A - Z)
       INTEGER NDIM, NFMIN, NFMAX, KN, NF, K1, K2
       PARAMETER (NDIM = 144, NFMIN = 3, NFMAX = 6)
       DOUBLE PRECISION BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX),
     1                  BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)
       DOUBLE PRECISION ASI, ASF, S
       DIMENSION ENS(3), ESG(2,2), ER(2)
*
* ---------------------------------------------------------------------
*
* ..Input common-blocks 
*
       COMMON / BETA   / BETA0, BETA1, BETA2, BETA3
       COMMON / PNS0   / P0NS (NDIM, NFMIN:NFMAX)
       COMMON / LSG    / R(NDIM, NFMIN:NFMAX, 2),
     1                   E(NDIM, NFMIN:NFMAX, 2, 2, 2)
*
* ---------------------------------------------------------------------
*
* ..The non-singlet evolution kernel
*   (second and third entries for technical reasons only) 
*
       ENS(1) = EXP (S * P0NS(KN,NF) / BETA0(NF))
       ENS(2) = ENS(1)
       ENS(3) = ENS(1)
*
* ..The singlet evolution operator
*  
       ER(1) = EXP (S * R(KN,NF,1))
       ER(2) = EXP (S * R(KN,NF,2))
*
       DO 1 K1 = 1, 2
       DO 1 K2 = 1, 2
         ESG(K1,K2) =   E(KN,NF,1,K1,K2) * ER(1) 
     1                + E(KN,NF,2,K1,K2) * ER(2)
   1   CONTINUE
*
* ---------------------------------------------------------------------
* 
       RETURN
       END
*
* =================================================================av==
