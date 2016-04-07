*
* ..File: lsgmom.f    (requires previous calls of PSG0MOM and BETAFCT)
*
*
* ..The subroutine  LSGMOM  for the eigenvalue decomposition of the LO
*    singlet splitting-function matrix P0SG (provided by the commom-
*    block  PSG0  on an NDIM-dimensional array of complex moments N), 
*    divided by the lowest coefficient BETA0 of the beta function of 
*    QCD (given in the common-block  BETA).  
*
* ..The eigenvalues  R(L),  L = 1, 2,  and the projectors  E(L,I,J), 
*    with  I,J = 1/2 = q/g  are stored in the common-block  LSG  for 
*    NF = NFLOW... NFHIGH  (maximally  NFMIN... NFMAX = 3...6) massless
*    quark flavours as specified by the common-block  NFUSED. 
*
* =====================================================================
*
*
       SUBROUTINE LSGMOM 
*
       IMPLICIT DOUBLE COMPLEX (A - Z)
       INTEGER NMAX, NDIM, NFMIN, NFMAX, NFLOW, NFHIGH, KN, NF
       PARAMETER (NDIM = 144, NFMIN = 3, NFMAX = 6)
       DOUBLE PRECISION BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX),
     1                  BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)
*
* ---------------------------------------------------------------------
*
* ..Input common-blocks
*
       COMMON / NNUSED / NMAX
       COMMON / NFUSED / NFLOW, NFHIGH
       COMMON / PSG0   / P0SG (NDIM, NFMIN:NFMAX, 2, 2)
       COMMON / BETA   / BETA0, BETA1, BETA2, BETA3
*
* ..Output common-block
*
       COMMON / LSG    / R(NDIM, NFMIN:NFMAX, 2),
     1                   E(NDIM, NFMIN:NFMAX, 2, 2, 2)                  
*
* ---------------------------------------------------------------------
*
* ..Mellin-N and flavour-number loops
*
       DO 1 KN = 1, NMAX
       DO 2 NF = NFLOW, NFHIGH
*
* ..The elements of R0
*
       B0I = 1./ BETA0(NF)
*
       RQQ = P0SG(KN,NF,1,1) * B0I
       RQG = P0SG(KN,NF,1,2) * B0I   
       RGQ = P0SG(KN,NF,2,1) * B0I
       RGG = P0SG(KN,NF,2,2) * B0I
*
* ..Eigenvalues and projection matrices 
*
       RX = SQRT ((RGG - RQQ) * (RGG - RQQ) + 4. * RQG * RGQ)
       R(KN,NF,2) = 0.5 * (RQQ + RGG + RX)
       R(KN,NF,1) = 0.5 * (RQQ + RGG - RX)
       RDFI = 1./ (R(KN,NF,1) - R(KN,NF,2)) 
*
       E(KN,NF,1,1,1) = RDFI * (RQQ - R(KN,NF,2)) 
       E(KN,NF,1,1,2) = RDFI *  RQG 
       E(KN,NF,1,2,1) = RDFI *  RGQ
       E(KN,NF,1,2,2) = RDFI * (RGG - R(KN,NF,2)) 
*
       E(KN,NF,2,1,1) = 1.- E(KN,NF,1,1,1)
       E(KN,NF,2,1,2) =   - E(KN,NF,1,1,2)
       E(KN,NF,2,2,1) =   - E(KN,NF,1,2,1)
       E(KN,NF,2,2,2) = 1.- E(KN,NF,1,2,2)
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
