*
* ..File: psg1mom.f      (requires previous call of PNS1MOM)
*
*
* ..The subroutine  PSG1MOM  for the NLO (alpha_s^2) flavour-singlet 
*    QCD splitting functions  P1SG  in N-space for  NF = NFMIN... NFMAX 
*    = 3... 6  massless quark flavours in the MS(bar) scheme.
*    The coupling constant is normalized as  a_s = alpha_s/(4*pi). 
*
* ..These quantities are determined on an external NDIM-dimensional 
*    array  NA  of complex Mellin moments provided by the common-block 
*    MOMS.  The results are written to the common-block  PSG1. The last 
*    two array arguments refer to the quark-gluon mixing matrix with 
*    1 = q  and  2 = g.
*
* ..The colour factors  CF, CA and TF  are taken from the common-block 
*    COLOUR.  The simple harmonic sums S_i(N) are provided by the 
*    common-block  HSUMS. 
*
* =====================================================================
*
*
       SUBROUTINE PSG1MOM 
*
       IMPLICIT DOUBLE COMPLEX (A - Z)
       INTEGER NMAX, NDIM, NFMIN, NFMAX, KN, NF
       PARAMETER (NDIM = 144, NFMIN = 3, NFMAX = 6)
       DOUBLE PRECISION CF, CA, TR
*
* ---------------------------------------------------------------------
*
* ..Input common-blocks 
*
       COMMON / MOMS   / NA (NDIM)
       COMMON / NNUSED / NMAX
       COMMON / HSUMS  / S(NDIM,6)
       COMMON / COLOUR / CF, CA, TR
*
       COMMON / PNS1   / P1NS (NDIM, NFMIN:NFMAX, 3)
       COMMON / SPSUMS / SSCHLP(NDIM), SSTR2P(NDIM), SSTR3P(NDIM)
*
* ..Output common-block
*
       COMMON / PSG1   / P1SG (NDIM, NFMIN:NFMAX, 2, 2)
*
* ---------------------------------------------------------------------
*
* ..Mellin-N loop
*
       DO 1 KN = 1, NMAX
*
* ..Some abbreviations 
*
       N  = NA(KN) 
       S1 = S(KN,1)
       S2 = S(KN,2)
*
       NS = N * N
       NT = NS * N
       NFO = NT * N
       NFI = NFO * N
       NSI = NFI * N
       NSE = NSI * N
       NE = NSE * N
       NN = NE * N
*
       NM = N - 1.
       N1 = N + 1.
       N2 = N + 2.
       NMS = NM * NM
       N1S = N1 * N1
       N1T = N1S * N1
       N2S = N2 * N2
       N2T = N2S * N2
*
* ---------------------------------------------------------------------
*
* ..The contributions to P1SG as given in Floratos et al. (1981) 
*   (SSTR2P, SSTR3P, SSCHLP  are taken from the common block  SPSUMS)
*
* ..Pure singlet (PS) and QG
* 
       PPSA = (5.* NFI + 32.* NFO + 49.* NT + 38.* NS + 28.* N + 8.) 
     1        / (NM * NT * N1T * N2S) * 2.
*
       PQGA = (-2.* S1 * S1 + 2.* S2 - 2.* SSTR2P(KN)) 
     1          * (NS + N + 2.) / (N * N1 * N2) 
     2        + (8.* S1 * (2.* N + 3.)) / (N1S * N2S)
     3        + 2.* (NN + 6.* NE + 15. * NSE + 25.* NSI + 36.* NFI
     4          + 85.* NFO + 128.* NT + 104.* NS + 64.* N + 16.)
     5          / (NM * NT * N1T * N2T)
       PQGB = (2.* S1 * S1 - 2.* S2 + 5.) * (NS + N + 2.)
     1          / (N * N1 * N2)   -   4.* S1 / NS
     2        + (11.* NFO + 26.* NT + 15.* NS + 8.* N + 4.)
     3          / (NT * N1T * N2) 
*
* ---------------------------------------------------------------------
*
* ..GQ and GG
*
       PGQA = (- S1 * S1 + 5.* S1 - S2) * (NS + N + 2.) 
     1          / (NM * N * N1)  -  2.* S1 / N1S
     2        - (12.* NSI + 30.* NFI + 43.* NFO + 28.* NT - NS
     3          - 12.* N - 4.) / (2.* NM * NT * N1T) 
       PGQB = (S1*S1 + S2 - SSTR2P(KN)) * (NS + N + 2.) / (NM * N * N1)
     1        - S1 * (17.* NFO + 41.* NS - 22.* N - 12.) 
     2          / (3.* NMS * NS * N1)
     3        + (109.* NN + 621.* NE + 1400.* NSE + 1678.* NSI
     4          + 695.* NFI - 1031.* NFO - 1304.* NT - 152.* NS
     5          + 432.* N + 144.) / (9.* NMS * NT * N1T * N2S)
       PGQC = (S1 - 8./3.) * (NS + N + 2.) / (NM * N * N1)  +  1./ N1S
       PGQC = 4./3.* PGQC
*  
       PGGA = - (2.* NFI + 5.* NFO + 8.* NT + 7.* NS - 2.* N - 2.)
     1          * 8.* S1 / (NMS * NS * N1S * N2S) -  67./9.* S1 + 8./3.
     2        - 4.* SSTR2P(KN) * (NS + N + 1.) / (NM * N * N1 * N2)
     3        + 2.* S1 * SSTR2P(KN) - 4.* SSCHLP(KN) + 0.5 * SSTR3P(KN)
     4        + (457.* NN + 2742.* NE + 6040.* NSE + 6098.* NSI
     5          + 1567.* NFI - 2344.* NFO - 1632.* NT + 560.* NS
     6          + 1488.* N + 576.) / (18.* NMS * NT * N1T * N2T)
       PGGB = (38.* NFO + 76.* NT + 94.* NS + 56.* N + 12.) * (-2.)
     1          / (9.* NM * NS * N1S * N2)  +  20./9.* S1  -  4./3.
       PGGC = (2.* NSI + 4.* NFI + NFO - 10.* NT - 5.* NS - 4.* N
     1          - 4.) * (-2.) / (NM * NT * N1T * N2)  -  1.
*
* ---------------------------------------------------------------------
*
* ..Flavour-number loop and output to the array 
*   (P1NS is provided by the common block PNS1)
*
       DO 2 NF = NFMIN, NFMAX
*
       P1SG(KN,NF,1,1) = P1NS (KN,NF,1) + TR*NF * CF * PPSA * 4.
       P1SG(KN,NF,1,2) = TR*NF * (CA * PQGA + CF * PQGB) * 4.
       P1SG(KN,NF,2,1) = ( CF*CF * PGQA + CF*CA * PGQB 
     1                   + TR*NF * CF * PGQC) * 4.
       P1SG(KN,NF,2,2) = ( CA*CA * PGGA 
     1                   + TR*NF * (CA * PGGB + CF * PGGC) ) * 4.
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
