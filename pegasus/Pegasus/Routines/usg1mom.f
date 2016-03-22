*
* ..File: usg1mom.f      (requires previous calls of PSG0MOM, PSG1MOM,
*                         LNSMOM and BETAFCT) 
*
*
* ..The subroutine USG1MOM for the NLO (a_s^1) singlet evolution matrix 
*    U1 in N-space for  NF = NFLOW... NFHIGH  (maximally NFMIN... NFMAX
*    = 3... 6)  massless quark flavours.  The matrix is written to the 
*    common-block  U1SG  on an NDIM-dim. array of complex moments.
*
* ..The eigenvalue decomposition of the LO splitting-function matrix,
*    the singlet splitting functions up to two loops, and the coeffi-
*    cients of the beta function are taken from the common-blocks  LSG, 
*    PSG,  PSG1,  and  BETA,  respectively.  The constant scale log
*    LOGFR = ln (mu_f^2/mu_r^2)  is taken from the common-block  FRRAT.
*    The n_f range is specified by the common-block  NFUSED.
*
* ..The output common block  R1SG  is used by the routine USG1HMOM.
*
* ..The notation follows  J. Bluemlein and A.V., PRD 58 (1998) 014020
*
* =====================================================================
*
*
       SUBROUTINE USG1MOM 
*
       IMPLICIT DOUBLE COMPLEX (A - Z)
       INTEGER NMAX, NDIM, NFMIN, NFMAX, NFLOW, NFHIGH, KN, NF, 
     1         J1, J2, K1, K2
       PARAMETER (NDIM = 144, NFMIN = 3, NFMAX = 6)
       DOUBLE PRECISION BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX),
     1                  BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)
       DOUBLE PRECISION B0I, B10S, LOGFR
       DIMENSION RT1(2,2), EM(2,2), EP(2,2)
*
* ---------------------------------------------------------------------
*
* ..Input common-blocks
*
       COMMON / NNUSED / NMAX
       COMMON / NFUSED / NFLOW, NFHIGH
       COMMON / BETA   / BETA0, BETA1, BETA2, BETA3
       COMMON / PSG0   / P0SG (NDIM, NFMIN:NFMAX, 2, 2)
       COMMON / PSG1   / P1SG (NDIM, NFMIN:NFMAX, 2, 2)
       COMMON / LSG    / R(NDIM, NFMIN:NFMAX, 2),
     1                   E(NDIM, NFMIN:NFMAX, 2, 2, 2)  
       COMMON / FRRAT  / LOGFR
*
* ..Output common-blocks
*
       COMMON / U1SG   / U1(NDIM, NFMIN:NFMAX, 2, 2)
       COMMON / R1SG   / R1(NDIM, NFMIN:NFMAX, 2, 2)
*
* ---------------------------------------------------------------------
*
* ..Mellin-N and flavour-number loops
*
       DO 1 KN = 1, NMAX
       DO 2 NF = NFLOW, NFHIGH
*
* ..Some abbreviations and the elements of R1
*   (including the contribution from mu_r unequal mu_f)
*
       B0I = 1./ BETA0(NF)
       B10S = BETA1(NF) * B0I * B0I
       RDIFF = R(KN,NF,1) - R(KN,NF,2)
*
       DO 11 J1 = 1, 2
       DO 11 J2 = 1, 2
         EM(J1,J2) = E(KN,NF,1,J1,J2)      
         EP(J1,J2) = E(KN,NF,2,J1,J2)      
         R1(KN,NF,J1,J2) =   P1SG(KN,NF,J1,J2) * B0I 
     1                     - P0SG(KN,NF,J1,J2) * (LOGFR + B10S)
  11   CONTINUE
*
* ---------------------------------------------------------------------
*
* ..The matrix elements of U1
*
       DO 12 K1 = 1, 2
       DO 12 K2 = 1, 2
         U1(KN,NF,K1,K2) = CMPLX (0.D0, 0.D0)
       DO 12 J1 = 1, 2
       DO 12 J2 = 1, 2
         U1(KN,NF,K1,K2) =  U1(KN,NF,K1,K2)  
     1        - EM(K1,J1) * R1(KN,NF,J1,J2) * EM(J2,K2)
     2        - EP(K1,J1) * R1(KN,NF,J1,J2) * EP(J2,K2)
     3        - EP(K1,J1) * R1(KN,NF,J1,J2) * EM(J2,K2) / (1.- RDIFF)
     4        - EM(K1,J1) * R1(KN,NF,J1,J2) * EP(J2,K2) / (1.+ RDIFF) 
  12   CONTINUE
*
* ---------------------------------------------------------------------
*
  2    CONTINUE
  1    CONTINUE
*
       RETURN
       END
*
* =================================================================av==
