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
       include 'pnsg_inc.f'
       INTEGER K1, K2
       DOUBLE PRECISION PGBETA0 (NFMIN:NFMAX), PGBETA1 (NFMIN:NFMAX),
     1                  PGBETA2 (NFMIN:NFMAX), PGBETA3 (NFMIN:NFMAX)
       DOUBLE PRECISION ASI, ASF, S
       DIMENSION ENS(3), ESG(2,2), ER(2)
*
* ---------------------------------------------------------------------
*
* ..Input common-blocks 
*
       COMMON / PGBETA   / PGBETA0, PGBETA1, PGBETA2, PGBETA3
!       COMMON / PNS0   / P0NS (NDIM, NFMIN:NFMAX)
!!$OMP THREADPRIVATE(/PNS0/)
!       COMMON / LSG    / R(NDIM, NFMIN:NFMAX, 2),
!     1                   E(NDIM, NFMIN:NFMAX, 2, 2, 2)
!!$OMP THREADPRIVATE(/LSG/)
*
* ---------------------------------------------------------------------
*
* ..The non-singlet evolution kernel
*   (second and third entries for technical reasons only) 
*
       ENS(1) = EXP (S * P0NS(KN,NF) / PGBETA0(NF))
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
