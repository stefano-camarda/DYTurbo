*
* ..File: psg2old.f     (requires previous call of PNS2MOM)
*
*
* ..The routine  PSG2MOM  for the approx. NNLO (alpha_s^3) singlet QCD 
*    splitting functions  P2NS  in N-space for  NF = NFMIN ... NFMAX 
*    = 3...6  massless quark flavours in the MS(bar) scheme.
*    The coupling constant is normalized as  a_s = alpha_s/(4*pi). 
*
* ..These quantities are determined on an external NDIM-dimensional 
*    array  NA  of complex Mellin moments provided by the common-block 
*    MOMS.  The results are written to the common-block  PSG2. The last
*    two array arguments refer to the quark-gluon mixing matrix with
*    1 = q  and  2 = g.
*
* ..The two approximations spanning the estimated residual uncertainty 
*    are invoked by  IAPP2 = 1  or  2,  for any other value of IAPP2
*    (input common-block  P2APPR) the central results are returned.  
*
* ..The QCD colour factors have been hard-wired in the approximations.
*    The harmonic sums S_i(N) are provided by the common-block  HSUMS.
*
* ..Reference: W.L. van Neerven and A. Vogt, hep-ph/0007362
*
* =====================================================================
*
*
       SUBROUTINE PSG2MOM 
*
       IMPLICIT DOUBLE COMPLEX (A - Z)
       INTEGER NMAX, NDIM, NFMIN, NFMAX, IAPP2, KN, NF
       PARAMETER (NDIM = 144, NFMIN = 3, NFMAX = 6)
*
* ---------------------------------------------------------------------
*
* ..Input common-blocks 
*
       COMMON / MOMS   / NA (NDIM)
       COMMON / NNUSED / NMAX
       COMMON / HSUMS  / S(NDIM,6)
       COMMON / P2APPR / IAPP2
*
       COMMON / PNS2   / P2NS (NDIM, NFMIN:NFMAX, 3)
*
* ..Output common-block 
*
       COMMON / PSG2   / P2SG (NDIM, NFMIN:NFMAX, 2, 2)
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
       S3 = S(KN,3)
       S4 = S(KN,4)
*
       NM = N - 1.
       N1 = N + 1.
       N2 = N + 2.
       NI = 1./N
       NMI = 1./NM
       N1I = 1./N1
       N2I = 1./N2
*
       S1M = S1 - NI
       S11 = S1 + N1I
       S21 = S2 + N1I*N1I
*
* ---------------------------------------------------------------------
*
*  ..Moments of the basic x-space functions 
*
       A0  = - S1M 
       B1  = - S1 * NI
       B11 = - S11 * N1I
       B2  = (S1**2 + S2) * NI
       B21 = (S11**2 + S21) * N1I
       B3  = - (S1**3 + 3.*S1*S2 + 2.*S3) * NI
       B4  = (S1**4 + 6.*S1**2*S2 + 8.*S1*S3 + 3.*S2**2 + 6.*S4) * NI
*
       C0  = NI
       CM  = NMI
       C1  = N1I
       C2  = N2I
       C3  = 1./(N+3.)
*
       D1  = - NI*NI
       D1M = - NMI*NMI
       D11 = - N1I*N1I
       D2  = 2.* NI**3
       D21 = 2.* N1I**3
       D3  = - 6.* NI**4
       D4  = 24.* NI**5
*
* ---------------------------------------------------------------------
*
* ..The approximate contributions to P2SG 
*
* ..QQ (pure singlet) and QG
*
       PPS1A = - 229.497 * (B1 - B11) - 722.99  * (C2 - C3)
     1         + 2678.77 * (C0 - C1)  - 560.20  * (CM - C0)
     2         + 2008.61 * D1  + 998.15 * D2 - 3584./27. * D1M
       PPS1B = + 73.845  * (B2 - B21) + 305.988 * (B1 - B11)
     1         + 2063.19 * (C1 - C2)  - 387.95  * (CM - C0)
     2         + 1999.35 * D11 - 732.68 * D1 - 3584./27. * D1M
       PPS2  = - 7.282   * (B1 - B11) - 38.779  * (C2 - C3)
     1         + 32.022  * (C1 - C2)  - 6.252   * (C0 - C1)
     2         + 1.767   * (CM - C0)  + 7.453   * D2
*
       PQG1A = - 31.830  * B3 + 1252.267 * B1 + 1999.89 * C1
     1         + 1722.47 * C0 +  1223.43 * D2 - 1334.61 * CM
     2         - 896./3. * D1M
       PQG1B = + 19.428  * B4 + 159.833  * B3 + 309.384 * B2
     1         + 2631.00 * (C0 - C1) -  67.25 * D2 - 776.793 * CM
     2         - 896./3. * D1M
       PQG2  = - 0.9085  * B2 -  35.803  * B1 - 128.023 * C0
     1         + 200.929 * (C0 - C1) + 40.542 * D1 +  3.284  * CM
*
* ---------------------------------------------------------------------
*
* ..GQ and GG
*
       PGQ0A = + 13.1212 * B4 + 126.665 * B3 + 308.536 * B2
     1         + 361.21  * C0 - 2113.45 * D1 - 17.965  * D1M
       PGQ0B = - 4.5108  * B4 - 66.618  * B3 - 231.535 * B2
     1         - 1224.22 * (C0 - C1) + 240.08 * D2
     2         + 379.60  * (D1M + 4.*CM)
       PGQ1A = + 2.4428  * B4 + 27.763  * B3 + 80.549  * B2
     1         - 227.14  * C0 - 151.04  * D2 + 65.91   * D1M
       PGQ1B = - 1.4028  * B4 - 11.638  * B3 + 164.964 * B1
     1         - 1066.78 * (C0 - C1) - 182.08 * D2
     2         + 138.54  * (D1M +2.*CM)
       PGQ2  = + 1.9361  * B2 + 11.178  * B1 + 11.632  * C0
     1         - 15.145  * (C0 - C1) + 3.354 * D1 - 2.133 * CM
*
       PGG0A = + 2626.38 * A0 + 4424.168
     1         - 732.715 * B2 - 20640.07 * C1 - 15428.58 * (C0 - C2)
     2         - 15213.6 * D2 + 16700.88 * CM + 2675.85 * D1M
       PGG0B = + 2678.22 * A0 + 4590.570
     1         + 3748.934* B1 - 35974.45 * (C0 + C2) + 60879.62 * C1
     2         + 2002.96 * D2 + 9762.09 * CM + 2675.85 * D1M
       PGG1A = - 415.71  * A0 - 548.569
     1         - 425.708 * B1 + 914.548 * C2 - 1122.86 * C0
     2         - 444.21  * D2 + 376.98  * CM + 157.18  * D1M
       PGG1B = - 412.00  * A0 - 534.951
     1         + 62.630  * B2 + 801.90 * C0 + 1891.40 * D1
     2         + 813.78  * D2 + 1.360  * CM + 157.18  * D1M
       PGG2  = - 16./9.D0* A0 + 6.4882
     1         + 37.6417 * C2 - 72.926  * C1 + 32.349  * C0
     2         - 0.991   * D2 + 2.818   * CM
*
* ---------------------------------------------------------------------
*
* ..Flavour-number loop and output to the array (depending on IAPP2)
*
       DO 2 NF = NFMIN, NFMAX
*
       IF (IAPP2 .EQ. 1) THEN
         P2SG(KN,NF,1,1) = P2NS(KN,NF,1) + NF * PPS1A + NF**2 * PPS2 
         P2SG(KN,NF,1,2) =         NF * PQG1A + NF**2 * PQG2 
         P2SG(KN,NF,2,1) = PGQ0A + NF * PGQ1A + NF**2 * PGQ2
         P2SG(KN,NF,2,2) = PGG0A + NF * PGG1A + NF**2 * PGG2
*
       ELSE IF (IAPP2 .EQ. 2) THEN
         P2SG(KN,NF,1,1) = P2NS(KN,NF,1) + NF * PPS1B + NF**2 * PPS2 
         P2SG(KN,NF,1,2) =         NF * PQG1B + NF**2 * PQG2 
         P2SG(KN,NF,2,1) = PGQ0B + NF * PGQ1B + NF**2 * PGQ2
         P2SG(KN,NF,2,2) = PGG0B + NF * PGG1B + NF**2 * PGG2
*
       ELSE
         P2SG(KN,NF,1,1) = P2NS(KN,NF,1) + NF * 0.5 * (PPS1A + PPS1B)
     1                     + NF**2 * PPS2
         P2SG(KN,NF,1,2) = NF * 0.5 * (PQG1A + PQG1B) + NF**2 * PQG2
         P2SG(KN,NF,2,1) = 0.5 * (PGQ0A + PGQ0B + NF * (PGQ1A + PGQ1B))
     1                     + NF**2 * PGQ2
         P2SG(KN,NF,2,2) = 0.5 * (PGG0A + PGG0B + NF * (PGG1A + PGG1B))
     1                     + NF**2 * PGG2
       END IF
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
