*
* ..File: psg1pmom.f      (requires previous call of PNSP1MOM)
*
*
* ..The subroutine  PSGP1MOM  for the NLO (alpha_s^2) polarized
*    flavour-singlet QCD splitting functions  P1SG  in N-space for 
*    NF = NFMIN ... NFMAX = 3 ... 6  massless quark flavours in the 
*    MS(bar) scheme.  The coupling constant is  a_s = alpha_s/(4*pi). 
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
       SUBROUTINE PSG1PMOM 
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
       COMMON / SPSUMS / SSCHLM(NDIM), SSTR2M(NDIM), SSTR3M(NDIM)
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
* ..The contributions to P1SG of Mertig and van Neerven (1995) using 
*    the rewriting by GRSV (1995). In this case the odd-N continuations 
*    SSTR2M, SSTR3M, SSCHLM  are taken from the common block  SPSUMS
*
* ..Pure singlet (PS) and QG
* 
       PPSA = - 2.* N2 * (1.+ 2.* N + NT) / (NT* N1T) 
*
       PQGA = ( (S1 * S1 - S2 + SSTR2M(KN)) * NM / (N * N1) 
     1        - 4.* S1 / (N * N1S)  -  (- 2. - 7.* N + 3.* NS - 4.* NT 
     2           + NFO + NFI) / (NT * N1T) ) * (-2.0)
       PQGB = ( (- S1*S1 + S2 + 2.* S1 / N) * NM / (N * N1)  
     1        -  NM * (1. + 3.5 * N + 4.* NS + 5.* NT + 2.5 * NFO) 
     2          / (NT * N1T)  +  4.* NM / (NS * N1S) ) * (-2.0)
*
* ..GQ and GG
*
       PGQA = ( 2.* (S1*S1 + S2) * N2 / (N * N1)  -  2.* S1 * N2 
     1          * (1.+ 3.* N) / (N * N1S)  -  N2 * (2.+ 15.* N 
     2          + 8.* NS - 12.* NT - 9.* NFO) / (NT * N1T) 
     3        + 8.* N2 / (NS * N1S) ) * (-0.5)
       PGQB = ( (- S1*S1 - S2 + SSTR2M(KN))* N2 / (N * N1) 
     1        +  S1 * (12.+ 22.* N + 11.* NS) / (3.* NS * N1)  
     2        -  (36.+ 72.* N + 41.* NS + 254.* NT + 271.* NFO 
     3          + 76.* NFI) / (9.* NT* N1T) ) * (-1.0)
       PGQC = (- S1 * N2 / (3.* N * N1) 
     1        +  N2 * (2.+ 5.* N) / (9.* N * N1S) ) * (-4.0)
*  
       PGGA = ( - 4.* S1 * SSTR2M(KN) - SSTR3M(KN) + 8.* SSCHLM(KN)
     1        +  8.* SSTR2M(KN) / (N * N1)  +  2.* S1 * (72.+ 144.* N
     2          + 67.* NS + 134.* NT + 67.* NFO) / (9.* NS* N1S)
     3        - (144.+ 258.* N + 7.* NS + 698.* NT + 469.* NFO
     4          + 144.* NFI + 48.* NSI) / (9.* NT * N1T) ) * (-0.5)
       PGGB = ( - 5.* S1 / 9. + (- 3.+ 13.* N + 16.* NS + 6.* NT 
     1          + 3.* NFO) / (9.* NS* N1S) ) * (-4.0)
       PGGC = ( 4.+ 2.* N - 8.* NS + NT + 5.* NFO + 3.* NFI + NSI) 
     1          / (NT * N1T) * (-1.0)
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
