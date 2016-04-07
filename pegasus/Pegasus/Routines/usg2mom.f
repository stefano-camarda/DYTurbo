*
* ..File: usg2mom.f      (requires previous calls of BETAFCT, PSG0MOM, 
*                         PSG1MOM, PSG2MOM, LNSMOM and USG1MOM)
*
*
* ..The subroutine USG2MOM for the NNLO (a_s^2) singlet evolution matrix
*    U2  in N-space for  NF = NFLOW ... NFHIGH  massless quark flavours.
*    The matrix is written to the common-block  U2SG  for an external 
*    NDIM-dimensional array of complex moments.
*
* ..The eigenvalue decomposition of the LO splitting-function matrix,
*    the singlet splitting functions up to three loops, the a_s^1
*    evolution matrix, and the coefficients of the beta function are 
*    taken from the common-blocks  LSG,  PSG0, PSG1, PSG2,  U1SG  and 
*    BETA,  resp. The constant scale log  LOGFR = ln (mu_f^2/mu_r^2)  
*    is taken from the common-block  FRRAT.  The n_f range is specified 
*    by the common-block  NFUSED.
*
* ..The output common block  R2SG  is used by the routine USG2HMOM.
*
* =====================================================================
*
*
       SUBROUTINE USG2MOM 
*
       IMPLICIT DOUBLE COMPLEX (A - Z)
       INTEGER NMAX, NDIM, NFMIN, NFMAX, NFLOW, NFHIGH, KN, NF, 
     1         J1, J2, K1, K2
       PARAMETER (NDIM = 144, NFMIN = 3, NFMAX = 6)
       DOUBLE PRECISION BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX),
     1                  BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)
       DOUBLE PRECISION B0, B1, B0I, B10, B20, LOGFR
       DIMENSION  R0(2,2), R1(2,2), RT2(2,2), EM(2,2), EP(2,2)
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
       COMMON / PSG2   / P2SG (NDIM, NFMIN:NFMAX, 2, 2)
       COMMON / LSG    / R(NDIM, NFMIN:NFMAX, 2),
     1                   E(NDIM, NFMIN:NFMAX, 2, 2, 2)  
       COMMON / U1SG   / U1(NDIM, NFMIN:NFMAX, 2, 2)
       COMMON / FRRAT  / LOGFR
*
* ..Output common-blocks
*
       COMMON / U2SG   / U2(NDIM, NFMIN:NFMAX, 2, 2)
       COMMON / R2SG   / R2(NDIM, NFMIN:NFMAX, 2, 2)
*
* ---------------------------------------------------------------------
*
* ..Mellin-N and flavour-number loops
*
       DO 1 KN = 1, NMAX
       DO 2 NF = NFLOW, NFHIGH
*
* ..Some abbreviations and the elements of R1 and R2
*   (including the contributions from mu_r unequal mu_f)
*
       B0  = BETA0(NF)
       B0I = 1./B0 
       B1  = BETA1(NF) 
       B10 = B1 * B0I 
       B20 = BETA2(NF) * B0I
       RDIFF = R(KN,NF,1) - R(KN,NF,2)
*
       DO 11 J1 = 1, 2
       DO 11 J2 = 1, 2
         EM(J1,J2) =  E(KN,NF,1,J1,J2)      
         EP(J1,J2) =  E(KN,NF,2,J1,J2)      
         R0(J1,J2) =  P0SG(KN,NF,J1,J2) * B0I
         R1(J1,J2) =  P1SG(KN,NF,J1,J2) * B0I 
     1              - R0(J1,J2) * (B10 + B0 * LOGFR)
         R2(KN,NF,J1,J2) 
     1             =  P2SG(KN,NF,J1,J2) * B0I
     2              - R1(J1,J2)* (B10 + 2.* B0* LOGFR)
     3              - R0(J1,J2)* (B20 + 3.* B1* LOGFR + B0**2* LOGFR**2)
  11   CONTINUE
*
* ---------------------------------------------------------------------
*
* ..The elements of RT2 = R(tilde)_2
*
       DO 12 K1 = 1, 2
       DO 12 K2 = 1, 2
         RT2(K1,K2) = R2(KN,NF,K1,K2)
       DO 12 J1 = 1, 2
         RT2(K1,K2) = RT2(K1,K2) + R1(K1,J1) * U1(KN,NF,J1,K2)
  12   CONTINUE
*
* ..The matrix elements of U2
*
       DO 13 K1 = 1, 2
       DO 13 K2 = 1, 2
         U2(KN,NF,K1,K2) = CMPLX (0.D0, 0.D0)
       DO 13 J1 = 1, 2
       DO 13 J2 = 1, 2
         U2(KN,NF,K1,K2) =  U2(KN,NF,K1,K2)  
     1        - EM(K1,J1) * RT2(J1,J2) * EM(J2,K2) / 2.
     2        - EP(K1,J1) * RT2(J1,J2) * EP(J2,K2) / 2.
     3        - EP(K1,J1) * RT2(J1,J2) * EM(J2,K2) / (2.- RDIFF)
     4        - EM(K1,J1) * RT2(J1,J2) * EP(J2,K2) / (2.+ RDIFF) 
  13   CONTINUE
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
