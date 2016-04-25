*
* ..File: usg2hmom.f    (requires previous calls of BETAFCT, LSGMOM
*                        USG1MOM and USG2MOM) 
*
*
* ..The subroutine  USG2HMOM  for the higher-order (a_s^3 etc) terms 
*    of the NNLO singlet matrix U in N-space for  NF = NFLOW ... NFHIGH  
*    massless quark flavours. The matrix is written to the common-block 
*    U2HSG  for an external NDIM-dimensional array of complex moments.
*
* ..The eigenvalue decomposition of the LO splitting-function matrix,
*    the first and second-order matrices R_1, U_1 and R_2, U_2, and the 
*    coefficients of the beta function are taken from the common-blocks 
*    LSG,  R1SG, U1SG,  R2SG, U2SG  and  BETA,  respectively. 
*
* ..Depending on  IMODE  given in  EVMOD,  the evolution equations are 
*    truncated for  dq/dln Q^2  (IMODE = 1)  or  dq/da_s  (else).
*
* =====================================================================
*
*
       SUBROUTINE USG2HMOM 
*
       IMPLICIT DOUBLE COMPLEX (A - Z)
       INTEGER NMAX, NDIM, NFMIN, NFMAX, NFLOW, NFHIGH, NUMAX, IMODE, 
     1         KN, NF, J1, J2, KO, K1, K2
       PARAMETER (NDIM = 144, NFMIN = 3, NFMAX = 6, NUMAX = 20)
       DOUBLE PRECISION BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX),
     1                  BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)
       DOUBLE PRECISION B0I, B10, B20
       DIMENSION EM(2,2), EP(2,2), RH(NUMAX,2,2), RT(NUMAX,2,2)
*
* ---------------------------------------------------------------------
*
* ..Input common-blocks
*
       COMMON / NNUSED / NMAX
       COMMON / NFUSED / NFLOW, NFHIGH
       COMMON / EVMOD  / IMODE
       COMMON / BETA   / BETA0, BETA1, BETA2, BETA3
       COMMON / LSG    / R(NDIM, NFMIN:NFMAX, 2),
     1                   E(NDIM, NFMIN:NFMAX, 2, 2, 2)  
       COMMON / U1SG   / U1(NDIM, NFMIN:NFMAX, 2, 2)
       COMMON / R1SG   / R1(NDIM, NFMIN:NFMAX, 2, 2)
       COMMON / U2SG   / U2(NDIM, NFMIN:NFMAX, 2, 2)
       COMMON / R2SG   / R2(NDIM, NFMIN:NFMAX, 2, 2)
*
* ..Output common-block
*
       COMMON / U2HSG  / U2H(NUMAX, NDIM, NFMIN:NFMAX, 2, 2)
*
* ---------------------------------------------------------------------
*
* ..Mellin-N and flavour-number loops
*
       DO 1 KN = 1, NMAX
       DO 2 NF = NFLOW, NFHIGH
*
* ..Some abbreviations for LO, NLO and NNLO quantities computed before
*
       B0I = 1./ BETA0(NF)
       B10 = BETA1(NF) * B0I
       B20 = BETA2(NF) * B0I
       RDIFF = R(KN,NF,1) - R(KN,NF,2)
*
       DO 11 J1 = 1, 2
       DO 11 J2 = 1, 2
         EM(J1,J2) = E(KN,NF,1,J1,J2)      
         EP(J1,J2) = E(KN,NF,2,J1,J2)      
         RH(1,J1,J2) = R1(KN,NF,J1,J2)
         RH(2,J1,J2) = R2(KN,NF,J1,J2)
         U2H(1,KN,NF,J1,J2) = U1(KN,NF,J1,J2)
         U2H(2,KN,NF,J1,J2) = U2(KN,NF,J1,J2)
  11   CONTINUE
*
* ---------------------------------------------------------------------
*
* ..Loop over the order in alpha_s
*
       DO 3 KO = 3, NUMAX
*
* ..The combinations  RH(KO) = R_KO  of the splitting functions 
*
       DO 12 K1 = 1, 2
       DO 12 K2 = 1, 2
       IF (IMODE .EQ. 1) THEN
         RH(KO,K1,K2) = - B10 * RH(KO-1,K1,K2) - B20 * RH(KO-2,K1,K2)
       ELSE
         RH(KO,K1,K2) = DCMPLX (0.D0, 0.D0)
       END IF
  12   CONTINUE
*                          
* ..The intermediate matrices RT(KO) = R(tilde)_KO 
*
       DO 13 K1 = 1, 2
       DO 13 K2 = 1, 2
         RT(KO,K1,K2) =  RH(KO,K1,K2)
       DO 13 J1 = 1, KO-1
       DO 13 J2 = 1, 2
         RT(KO,K1,K2) =  RT(KO,K1,K2) 
     1                 + RH(J1,K1,J2) * U2H(KO-J1,KN,NF,J2,K2)
  13     CONTINUE
*
* ---------------------------------------------------------------------
*
* ..The NNLO evolution matrices U2H(KO) = U_KO for KO > 2
*
       DO 14 K1 = 1, 2
       DO 14 K2 = 1, 2
         U2H(KO,KN,NF,K1,K2) = DCMPLX (0.D0, 0.D0)
       DO 14 J1 = 1, 2
       DO 14 J2 = 1, 2
         U2H(KO,KN,NF,K1,K2) =  U2H(KO,KN,NF,K1,K2)  
     1        - EM(K1,J1) * RT(KO,J1,J2) * EM(J2,K2) / KO
     2        - EP(K1,J1) * RT(KO,J1,J2) * EP(J2,K2) / KO
     3        - EP(K1,J1) * RT(KO,J1,J2) * EM(J2,K2) / (KO - RDIFF)
     4        - EM(K1,J1) * RT(KO,J1,J2) * EP(J2,K2) / (KO + RDIFF) 
  14   CONTINUE
*
* ---------------------------------------------------------------------
*
  3    CONTINUE
  2    CONTINUE
  1    CONTINUE
*
       RETURN
       END
*
* =================================================================av==
