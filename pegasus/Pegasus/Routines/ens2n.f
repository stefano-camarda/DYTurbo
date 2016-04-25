*
* ..File ens2n.f       (requires previous calls of BETAFCT, PNS0MOM, 
*                       and UNS2MOM)
*                  
*
* ..The subroutine  ENS2N  for the NNLO hadronic non-singlet evolution 
*    kernels in N-space at a fixed numer of flavours  NF.  The routine 
*    returns, for  K = MSMIN ... NSMAX,  the kernels  ENS(K)  for the 
*    evolution of the `+' (K=1), `-' (K=2) and 'V' (K=3)  quark combin-
*    ations for a moment N given via the counter  KN.  The initial and 
*    final scales are specified by the respective values  ASI  and  ASF 
*    for a_s = alpha_s/(4 pi), assuming a fixed ratio of mu_r and mu_f.
*
* ..Depending on  IMODE  given in the common-block  EVMOD  the kernels
*    are determined in one of four schemes for solving the non-singlet
*    evolution equations at NNLO. 
*
* ..The lowest-order splitting function and the coefficients of the
*    kernels U and of the beta function are given by the common-blocks 
*    PNS0,  U2NS  and  BETA,  respectively.
*
* =====================================================================
*
*
       SUBROUTINE ENS2N (ENS, ASI, ASF, S, KN, NF, NSMIN, NSMAX)
*
       IMPLICIT DOUBLE COMPLEX (A - Z)
       INTEGER KN, NF, NSMIN, NSMAX, NDIM, NFMIN, NFMAX, NUMAX, NUORD, 
     1         IMODE, KO, K1
       PARAMETER (NDIM = 144, NFMIN = 3, NFMAX = 6, NUMAX = 20)
       DOUBLE PRECISION BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX),
     1                  BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)
       DOUBLE PRECISION ASI, ASF, S, LOGFR, AFI1, AFI2, ASIO, ASFO
       DIMENSION ENS(3), UNSF(NUMAX), UNSI(NUMAX)
       PARAMETER ( ONE = (1.D0, 0.D0) )
*
* ---------------------------------------------------------------------
*
* ..Input common-blocks 
*
       COMMON / BETA   / BETA0, BETA1, BETA2, BETA3
       COMMON / PNS0   / P0NS (NDIM, NFMIN:NFMAX)
       COMMON / U2NS   / UNS2(NUMAX, NDIM, NFMIN:NFMAX, 3)
       COMMON / EVMOD  / IMODE
       COMMON / ITORD  / NUORD
*
* ---------------------------------------------------------------------
*
* ..The LO evolution operator 
*
       LNS = EXP (S * P0NS(KN,NF) / BETA0(NF))
*
* ---------------------------------------------------------------------
*
       IF ( (IMODE .EQ. 1) .OR. (IMODE .EQ. 2) ) THEN
*
* ..1,2) the two iterated solutions (difference: UNS2 for KO > 2)
*
       DO 21 K1 = NSMIN, NSMAX
         UNSF(K1) = ONE
         UNSI(K1) = ONE
  21   CONTINUE
*
       ASFO = 1.D0
       ASIO = 1.D0
       DO 22 KO = 1, NUORD
         ASFO = ASFO * ASF
         ASIO = ASIO * ASI
       DO 22 K1 = NSMIN, NSMAX
         UNSF(K1) = UNSF(K1) + ASFO * UNS2(KO,KN,NF,K1) 
         UNSI(K1) = UNSI(K1) + ASIO * UNS2(KO,KN,NF,K1) 
  22   CONTINUE
*
       DO 11 K1 = NSMIN, NSMAX
         ENS(K1) = UNSF(K1) * LNS / UNSI(K1)
  11   CONTINUE
*
* ---------------------------------------------------------------------
*
       ELSE IF ( IMODE .EQ. 3 ) THEN
*
* ..3) the truncated solution without expansion of U^(-1)
*
       DO 12 K1 = NSMIN, NSMAX
         ENS(K1) =  ( ONE + ASF * (UNS2(1,KN,NF,K1) 
     1                    + ASF *  UNS2(2,KN,NF,K1)) ) / 
     2              ( ONE + ASI * (UNS2(1,KN,NF,K1) 
     3                    + ASI *  UNS2(2,KN,NF,K1)) ) * LNS
  12   CONTINUE
*
       ELSE
*
* ..4) the fully truncated solution -- default
*
       AFI1 = (ASF - ASI) * (1.- ASI)
       AFI2 = ASF*ASF - ASI*ASI
*
       DO 13 K1 = NSMIN, NSMAX
         ENS(K1) = LNS * ( ONE + AFI1 * UNS2(1,KN,NF,K1) 
     1                         + AFI2 * UNS2(2,KN,NF,K1) )
  13   CONTINUE
*
       END IF
*
* ---------------------------------------------------------------------
* 
       RETURN
       END
*
* =================================================================av==
