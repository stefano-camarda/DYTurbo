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
!       INTEGER NMAX, NDIM, NFMIN, NFMAX, KN, NF
!       include 'dimensions.f'
!       DOUBLE PRECISION ZETA(6), CF, CA, TR
!*
!* ---------------------------------------------------------------------
!*
!* ..Input common-blocks 
!*
!       COMMON / MOMS   / NA (NDIM)
!       COMMON / NNUSED / NMAX
!       COMMON / HSUMS  / S(NDIM,6)
!       COMMON / COLOUR / CF, CA, TR
!       COMMON / RZETA  / ZETA
!*
!* ..Output common-blocks 
!*
!       COMMON / PNS1   / P1NS (NDIM, NFMIN:NFMAX, 3)
!       COMMON / SPSUMS / SSCHLP(NDIM), SSTR2P(NDIM), SSTR3P(NDIM)
!*
!* ---------------------------------------------------------------------
      include 'pnsg_inc.f'
      include 'const_inc.f'
      include 'moms_inc.f'
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
       N1 = N + 1D0
       N2 = N + 2D0
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
!       N3 = N + 3D0
!       N4 = N + 4D0
!       N5 = N + 5D0
!       N6 = N + 6D0
!       S11 = S1  + 1D0/N1
!       S12 = S11 + 1D0/N2
!       S13 = S12 + 1D0/N3
!       S14 = S13 + 1D0/N4
!       S15 = S14 + 1D0/N5
!       S16 = S15 + 1D0/N6
!       SPMOM = 1.0000D0 * (ZETA(2) - S1 / N ) / N  -
!     1         0.9992D0 * (ZETA(2) - S11/ N1) / N1 +
!     2         0.9851D0 * (ZETA(2) - S12/ N2) / N2 -
!     3         0.9005D0 * (ZETA(2) - S13/ N3) / N3 +
!     4         0.6621D0 * (ZETA(2) - S14/ N4) / N4 -
!     5         0.3174D0 * (ZETA(2) - S15/ N5) / N5 +
!     6         0.0699D0 * (ZETA(2) - S16/ N6) / N6  
       SPMOM=fun6(N-1D0) ! --> Improved approximation for the Mellin transform of: Li2(x)/(1+x)
*
       SLC = - 5D0/8D0 * ZETA(3)
       SLV = - ZETA(2)/2D0* (PSI(N1/2D0) - PSI(N/2D0)) + S1/NS + SPMOM
       SSCHLM = SLC - SLV
       SSTR2M = ZETA(2) - DPSI (N1/2D0,1)
       SSTR3M = 0.5D0 * DPSI (N1/2D0,2) + ZETA(3)
       SSCHLP(KN) = SLC + SLV
       SSTR2P(KN) = ZETA(2) - DPSI (N2/2D0,1)
       SSTR3P(KN) = 0.5D0 * DPSI (N2/2D0,2) + ZETA(3)
*
* ---------------------------------------------------------------------
*
* ..The contributions to P1NS as given in Gonzalez-Arroyo et al. (1979) 
*   (Note that the anomalous dimensions in the literature often differ 
*    from these moments of the splitting functions by factors -1 or -2,
*    in addition to possible different normalizations of the coupling)
*
       PNMA = ( 16D0* S1 * (2D0* N + 1D0) / (NS * N1S) +
     1          16D0* (2D0* S1 - 1D0/(N * N1)) * ( S2 - SSTR2M ) +
     2          64D0* SSCHLM + 24D0* S2 - 3D0 - 8D0* SSTR3M -
     3          8D0* (3D0* NT + NS -1D0) / (NT * N1T) +
     4          16D0* (2D0* NS + 2D0* N +1D0) / (NT * N1T) ) * (-0.5D0)
       PNPA = ( 16D0* S1 * (2D0* N + 1D0) / (NS * N1S) +
     1          16D0* (2D0* S1 - 1D0/(N * N1)) * ( S2 - SSTR2P(KN) ) +
     2          64D0* SSCHLP(KN) + 24D0* S2 - 3D0 - 8D0* SSTR3P(KN) -
     3          8D0* (3D0* NT + NS -1D0) / (NT * N1T) -
     4          16D0* (2D0* NS + 2D0* N +1D0) / (NT * N1T) ) * (-0.5D0)
       PNSB = ( S1 * (536D0/9D0 + 8D0* (2D0* N + 1D0) / (NS * N1S)) -
     1          (16D0* S1 + 52D0/3D0- 8D0/(N * N1)) * S2 - 43D0/6D0 -
     2          (151D0* NFO + 263D0* NT + 97D0* NS + 3D0* N + 9D0) *
     3          4D0/ (9D0* NT * N1T) ) * (-0.5D0)
       PNSC = ( -160D0/9D0* S1 + 32D0/3D0* S2 + 4D0/3D0 +
     1      16D0* (11D0* NS + 5D0* N - 3D0) / (9D0* NS * N1S) )*(-0.5D0)
*
* ---------------------------------------------------------------------
*
* ..Flavour-number loop and output to the array
*
       DO 2 NF = NFMIN, NFMAX
*
       P1NS(KN,NF,1) = CF * ((CF-CA/2D0)* PNPA + CA* PNSB + TR*NF* PNSC)
       P1NS(KN,NF,2) = CF * ((CF-CA/2D0)* PNMA + CA* PNSB + TR*NF* PNSC)
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
