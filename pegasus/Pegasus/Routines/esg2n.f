*
* ..File esg2n.f       (requires previous calls of LSGMOM, USG1MOM,
*                       USG2MOM, and of USG2HMOM for IMODE = 1 or 2)
*                  
*
* ..The subroutine  ESG2N  for the NNLO hadronic singlet evolution 
*    operators in N-space at a fixed number of flavours  NF.  
*    The kernels  ESG(K1,K2),  K1, K2 = 1, 2  with  1 = q  and  2 = g,  
*    are returned for a Mellin moment N specified via the counter  KN. 
*    The initial and final scales are specified by the respective 
*    values  ASI  and  ASF  of  a_s = alpha_s/(4 pi),  assuming a 
*    fixed ratio of mu_r and mu_f.
*
* ..Depending on  IMODE  given in the common-block  EVMOD,  the kernels 
*    are determined in one of three schemes for solving the evolution 
*    equations at NNLO.  The order  NUORD (=< NUMAX = 20)  in a_s for 
*    the two iterative schemes (IMODE = 1, 2) is specified in  ITORD.
*      
* ..The eigenvalue decomposition of the LO splitting-function matrix,
*    U_1 and U_2, and (IMODE = 1, 2) the higher-order U-matrices are 
*    taken from the common-blocks  LSG,  U1SG/U2SG  and  U2HSG,  resp.
*    The 2-dimensional Konecker symbol is provided by  KRON2D.  
*
* =====================================================================
*
*
       SUBROUTINE ESG2N (ESG, ASI, ASF, S, KN, NF)
*
       IMPLICIT DOUBLE COMPLEX (A - Z)
       INTEGER NDIM, NFMIN, NFMAX, NUMAX, NUORD, KN, NF, IMODE, 
     1         KO, K1, K2, J1, J2
       PARAMETER (NDIM = 144, NFMIN = 3, NFMAX = 6, NUMAX = 20)
       DOUBLE PRECISION ASI, ASF, S, ASFO, ASIO
       DIMENSION ESG(2,2), ER(2), L(2,2), UF(2,2), UI(2,2), UM(2,2), 
     1           UM2(2,2)
       PARAMETER ( ZERO = (0.D0, 0.D0) )
*
* ---------------------------------------------------------------------
*
* ..Input common-blocks 
*
       COMMON / KRON2D / D(2,2)
       COMMON / EVMOD  / IMODE
       COMMON / ITORD  / NUORD
       COMMON / LSG    / R(NDIM, NFMIN:NFMAX, 2),
     1                   E(NDIM, NFMIN:NFMAX, 2, 2, 2)
       COMMON / U1SG   / U1(NDIM, NFMIN:NFMAX, 2, 2)
       COMMON / U2SG   / U2(NDIM, NFMIN:NFMAX, 2, 2)
       COMMON / U2HSG  / U2H(NUMAX, NDIM, NFMIN:NFMAX, 2, 2)
*
* ---------------------------------------------------------------------
*
* ..The LO evolution operator 
*
       ER(1) = EXP (S * R(KN,NF,1))
       ER(2) = EXP (S * R(KN,NF,2))
*
       DO 11 K1 = 1, 2
       DO 11 K2 = 1, 2
         L(K1,K2) =   E(KN,NF,1,K1,K2) * ER(1) 
     1              + E(KN,NF,2,K1,K2) * ER(2)
  11   CONTINUE
*
* ---------------------------------------------------------------------
*
       IF ( (IMODE .NE. 1) .AND. (IMODE .NE. 2) ) THEN
*
* ..The second order term UM2 of the expanded inverse of U
*
       DO 12 K1 = 1, 2
       DO 12 K2 = 1, 2
         UM2(K1,K2) = - U2(KN,NF,K1,K2)
       DO 12 J1 = 1, 2
         UM2(K1,K2) = UM2(K1,K2) + U1(KN,NF,K1,J1) * U1(KN,NF,J1,K2)
  12   CONTINUE
*
* ..The truncated NNLO evolution operator -- default
*
       DO 13 K1 = 1, 2
       DO 13 K2 = 1, 2
         ESG(K1,K2) = L(K1,K2)
       DO 13 J1 = 1, 2
         ESG(K1,K2) = ESG(K1,K2)  +  ASF * ( U1(KN,NF,K1,J1) + ASF * 
     1                U2(KN,NF,K1,J1) ) * L(J1,K2)  -  L(K1,J1) * ASI  
     2                * ( U1(KN,NF,J1,K2) - ASI * UM2(J1,K2) )
       DO 13 J2 = 1, 2
         ESG(K1,K2) = ESG(K1,K2) - ASF * ASI *
     3                U1(KN,NF,K1,J1) * L(J1,J2) * U1(KN,NF,J2,K2)
  13   CONTINUE
       RETURN
*
       END IF      
*
* ---------------------------------------------------------------------
*
* ..Iterative solutions: the U-matrices at the final and initial scales 
*
       DO 14 K1 = 1, 2
       DO 14 K2 = 1, 2
         UF(K1,K2) = D(K1,K2) 
         UI(K1,K2) = D(K1,K2) 
  14   CONTINUE
*
       ASFO = 1.D0
       ASIO = 1.D0
       DO 15 KO = 1, NUORD
         ASFO = ASFO * ASF
         ASIO = ASIO * ASI
       DO 15 K1 = 1, 2
       DO 15 K2 = 1, 2
         UF(K1,K2) = UF(K1,K2) + ASFO * U2H(KO,KN,NF,K1,K2)
         UI(K1,K2) = UI(K1,K2) + ASIO * U2H(KO,KN,NF,K1,K2)
  15   CONTINUE 
*
* ..The full inverse UM of UI
*
       DETINV = 1./ ( UI(1,1)*UI(2,2) - UI(1,2)*UI(2,1) ) 
       UM(1,1) =  UI(2,2) * DETINV
       UM(1,2) = -UI(1,2) * DETINV
       UM(2,1) = -UI(2,1) * DETINV
       UM(2,2) =  UI(1,1) * DETINV
*
* ---------------------------------------------------------------------
* 
* ..The NNLO evolution operators for IMODE = 1 or 2
*
       DO 16 K1 = 1, 2
       DO 16 K2 = 1, 2
         ESG(K1,K2) = ZERO
       DO 16 J1 = 1, 2
       DO 16 J2 = 1, 2
         ESG(K1,K2) = ESG(K1,K2) + UF(K1,J1) * L(J1,J2) * UM(J2,K2) 
  16   CONTINUE
*
* ---------------------------------------------------------------------
* 
       RETURN
       END
*
* =================================================================av==
