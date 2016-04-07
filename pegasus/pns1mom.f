*
* ..File: pns1mom.f  
*
*
* ..The subroutine  PNS1MOM  for the NLO (alpha_s^2) non-singlet QCD 
*    splitting functions  P1NS  in N-space for  NF = NFMIN ... NFMAX 
*    = 3...6  massless quark flavours in the MS(bar) scheme.
*    The coupling constant is normalized as  a_s = alpha_s/(4*pi). 
*
* ..These quantities are determined on an external NDIM-dimensional 
*    array  NA  of complex Mellin moments provided by the common-block 
*    MOMS.  The results are written to the common-block  PNS1.  The 
*    last array argument refers to the `+', `-' and `V' NS combinations
*    with  1 = +,  2 = -,  and  3 = V = -  (at this order). 
*
* ..The colour factors  CF, CA and TF  are taken from the common-block
*    COLOUR.  The simple harmonic sums S_i(N) are provided by the
*    common-block  HSUMS.  The lowest integer values of the Riemann
*    Zeta-function are provided by the common-block  RZETA.
*
* ..The output common block  SPSUMS  is used by the routine PSG1MOM.
*
* =====================================================================
*
*
       SUBROUTINE PNS1MOM 
*
       IMPLICIT DOUBLE COMPLEX (A - Z)
       INTEGER NMAX, NDIM, NFMIN, NFMAX, KN, NF
       PARAMETER (NDIM = 144, NFMIN = 3, NFMAX = 6)
       DOUBLE PRECISION ZETA(6), CF, CA, TR
*
* ---------------------------------------------------------------------
*
* ..Input common-blocks 
*
       COMMON / MOMS   / NA (NDIM)
       COMMON / NNUSED / NMAX
       COMMON / HSUMS  / S(NDIM,6)
       COMMON / COLOUR / CF, CA, TR
       COMMON / RZETA  / ZETA
*
* ..Output common-blocks 
*
       COMMON / PNS1   / P1NS (NDIM, NFMIN:NFMAX, 3)
       COMMON / SPSUMS / SSCHLP(NDIM), SSTR2P(NDIM), SSTR3P(NDIM)
*
* ---------------------------------------------------------------------
*
* ..Begin of the Mellin-N loop
*
       DO 1 KN = 1, NMAX
*
* ..Some abbreviations 
*
       N  = NA(KN) 
       S1 = S(KN,1)
       S2 = S(KN,2)
*
       N1 = N + 1.
       N2 = N + 2.
       NS = N * N
       NT = NS * N
       NFO = NT * N
       N1S = N1 * N1
       N1T = N1S * N1
*
* ---------------------------------------------------------------------
*
* ..Analytic continuations of the occuring sums as given in GRV (1990) 
*   (with an improved parametrization of the moments of  Sp(x)/(1+x).)
*
       N3 = N + 3.
       N4 = N + 4.
       N5 = N + 5.
       N6 = N + 6.
       S11 = S1  + 1./N1
       S12 = S11 + 1./N2
       S13 = S12 + 1./N3
       S14 = S13 + 1./N4
       S15 = S14 + 1./N5
       S16 = S15 + 1./N6
       SPMOM = 1.0000D0 * (ZETA(2) - S1 / N ) / N  -
     1         0.9992D0 * (ZETA(2) - S11/ N1) / N1 +
     2         0.9851D0 * (ZETA(2) - S12/ N2) / N2 -
     3         0.9005D0 * (ZETA(2) - S13/ N3) / N3 +
     4         0.6621D0 * (ZETA(2) - S14/ N4) / N4 -
     5         0.3174D0 * (ZETA(2) - S15/ N5) / N5 +
     6         0.0699D0 * (ZETA(2) - S16/ N6) / N6  
*
       SLC = - 5./8.D0 * ZETA(3)
       SLV = - ZETA(2)/2.* (PSI(N1/2.) - PSI(N/2.)) + S1/NS + SPMOM
       SSCHLM = SLC - SLV
       SSTR2M = ZETA(2) - DPSI (N1/2.,1)
       SSTR3M = 0.5 * DPSI (N1/2.,2) + ZETA(3)
       SSCHLP(KN) = SLC + SLV
       SSTR2P(KN) = ZETA(2) - DPSI (N2/2.,1)
       SSTR3P(KN) = 0.5 * DPSI (N2/2.,2) + ZETA(3)
*
* ---------------------------------------------------------------------
*
* ..The contributions to P1NS as given in Gonzalez-Arroyo et al. (1979) 
*   (Note that the anomalous dimensions in the literature often differ 
*    from these moments of the splitting functions by factors -1 or -2,
*    in addition to possible different normalizations of the coupling)
*
       PNMA = ( 16.* S1 * (2.* N + 1.) / (NS * N1S) +
     1          16.* (2.* S1 - 1./(N * N1)) * ( S2 - SSTR2M ) +
     2          64.* SSCHLM + 24.* S2 - 3. - 8.* SSTR3M -
     3          8.* (3.* NT + NS -1.) / (NT * N1T) +
     4          16.* (2.* NS + 2.* N +1.) / (NT * N1T) ) * (-0.5)
       PNPA = ( 16.* S1 * (2.* N + 1.) / (NS * N1S) +
     1          16.* (2.* S1 - 1./(N * N1)) * ( S2 - SSTR2P(KN) ) +
     2          64.* SSCHLP(KN) + 24.* S2 - 3. - 8.* SSTR3P(KN) -
     3          8.* (3.* NT + NS -1.) / (NT * N1T) -
     4          16.* (2.* NS + 2.* N +1.) / (NT * N1T) ) * (-0.5)
       PNSB = ( S1 * (536./9. + 8.* (2.* N + 1.) / (NS * N1S)) -
     1          (16.* S1 + 52./3.- 8./(N * N1)) * S2 - 43./6. -
     2          (151.* NFO + 263.* NT + 97.* NS + 3.* N + 9.) *
     3          4./ (9.* NT * N1T) ) * (-0.5)
       PNSC = ( -160./9.* S1 + 32./3.* S2 + 4./3. +
     1          16.* (11.* NS + 5.* N - 3.) / (9.* NS * N1S) ) * (-0.5)
*
* ---------------------------------------------------------------------
*
* ..Flavour-number loop and output to the array
*
       DO 2 NF = NFMIN, NFMAX
*
       P1NS(KN,NF,1) = CF * ((CF-CA/2.)* PNPA + CA* PNSB + TR*NF* PNSC)
       P1NS(KN,NF,2) = CF * ((CF-CA/2.)* PNMA + CA* PNSB + TR*NF* PNSC)
       P1NS(KN,NF,3) = P1NS(KN,NF,2)
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
