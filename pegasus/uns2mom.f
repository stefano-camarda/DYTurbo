*
* ..File: uns2mom.f     (requires previous calls of BETAFCT, PNS0MOM,
*                        PNS1MOM and PNS2MOM) 
*
*
* ..The subroutine UNS2MOM for the NNLO non-singlet evolution operators 
*    U_a, a = +, - and V  (including the higher-order pieces for the 
*    iterative solutions) in Mellin space for  NF = NFLOW ... NFHIGH  
*    (maximally  NFMIN... NFMAX = 3... 6)  massless quark flavours. 
*    The results  UNS2  are written to the common-block  U2NS  for an 
*    external NDIM-dimensional array of complex moments. The notation
*    for the last array argument reads  1 = +, 2 = -, 3 = V.
*
* ..The splitting functions up to three loops and the coefficients of 
*    the beta function are taken from the common-blocks  PNS0,  PNS1, 
*    PNS2  and  BETA,  respectively.  The constant scale logarithm
*    LOGFR = ln (mu_f^2/mu_r^2)  is taken from the common-block  FRRAT.
*    The n_f range is specified by the common-block  NFUSED.
* 
* ..Depending on  IMODE  provided by  EVMOD,  the evolution equations 
*    are truncated for  dq/dln Q^2  (IMODE = 1)  or  dq/da_s  (else).
*
* =====================================================================
*
*
       SUBROUTINE UNS2MOM 
*
       IMPLICIT DOUBLE COMPLEX (A - Z)
       INTEGER NMAX, NDIM, NFMIN, NFMAX, NFLOW, NFHIGH, NUMAX, IMODE, 
     1         KN, NF, KO, K1, K2
       PARAMETER (NDIM = 144, NFMIN = 3, NFMAX = 6, NUMAX = 20)
       DOUBLE PRECISION BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX),
     1                  BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)
       DOUBLE PRECISION B0, B0I, B1, B10, B20, LOGFR
       DIMENSION RNS(NUMAX,3)
*
* ---------------------------------------------------------------------
*
* ..Input common-blocks
*
       COMMON / NNUSED / NMAX
       COMMON / NFUSED / NFLOW, NFHIGH
       COMMON / EVMOD  / IMODE
       COMMON / FRRAT  / LOGFR
       COMMON / BETA   / BETA0, BETA1, BETA2, BETA3
       COMMON / PNS0   / P0NS (NDIM, NFMIN:NFMAX)
       COMMON / PNS1   / P1NS (NDIM, NFMIN:NFMAX, 3)
       COMMON / PNS2   / P2NS (NDIM, NFMIN:NFMAX, 3)
*
* ..Output common-block
*
       COMMON / U2NS   / UNS2(NUMAX, NDIM, NFMIN:NFMAX, 3)
*
* ---------------------------------------------------------------------
*
* ..Mellin-N and flavour-number loops
*
       DO 1 KN = 1, NMAX
       DO 2 NF = NFLOW, NFHIGH
*
* ..Some abbreviations
*
       B0  = BETA0(NF)
       B0I = 1./B0
       B1  = BETA1(NF)
       B10 = B1 * B0I
       B20 = BETA2(NF) * B0I
*
       AUX1 = B10 + B0 * LOGFR
       AUX2 = B10 + 2.* B0 * LOGFR
       AUX3 = B20 + B0**2 * LOGFR**2 + 3.* B1 * LOGFR
*
* ..The first- and second-order quantities
*
       R0NS  = P0NS(KN,NF) * B0I
       DO 11 K1 = 1, 3
         RNS(1,K1) =  P1NS(KN,NF,K1) * B0I - R0NS * AUX1
         RNS(2,K1) =  P2NS(KN,NF,K1) * B0I - RNS(1,K1) * AUX2 
     1              - R0NS * AUX3
         UNS2(1,KN,NF,K1) = - RNS(1,K1)
         UNS2(2,KN,NF,K1) = - 0.5 * ( RNS(2,K1) - RNS(1,K1)*RNS(1,K1) )
  11   CONTINUE
*
* ---------------------------------------------------------------------
*
* ..Loop over the higher orders in alpha_s 
*
       DO 3 KO = 3, NUMAX
*
* ..The combinations  RH(KO) = R_KO  of the splitting functions 
*
       DO 12 K1 = 1, 3
       IF (IMODE .EQ. 1) THEN
         RNS(KO,K1) = - B10 * RNS(KO-1,K1) - B20 * RNS(KO-2,K1)
       ELSE
         RNS(KO,K1) = DCMPLX (0.D0, 0.D0)
       END IF
  12   CONTINUE
*
* ..The KO > 2 coefficients U2NS(KO) = U_KO of the evolution operators
*
       DO 13 K1 = 1, 3
         UNS2(KO,KN,NF,K1) = - RNS(KO,K1) / KO
       DO 13 K2 = 1, KO-1
         UNS2(KO,KN,NF,K1) = UNS2(KO,KN,NF,K1)
     1                       - (RNS(K2,K1) * UNS2(KO-K2,KN,NF,K1)) / KO
  13   CONTINUE
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
