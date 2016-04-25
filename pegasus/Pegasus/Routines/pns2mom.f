*
* ..File: pns2mom.f
*
*
* ..The subroutine  PNS2MOM  for the parametrized NNLO non-singlet QCD 
*    splitting functions  P2NS  in N-space for  NF = NFMIN ... NFMAX 
*    = 3...6  massless quark flavours in the MS(bar) scheme.
*    The coupling constant is normalized as  a_s = alpha_s/(4*pi). 
*
* ..These quantities are determined on an external NDIM-dimensional 
*    array  NA  of complex Mellin moments provided by the common-block 
*    MOMS.  The results are written to the common-block  PNS2.  The 
*    last array argument refers to the `plus', `minus', and `valence'
*    NS combinations with  1 = +,  2 = -,  and  3 = V.
*
* ..The QCD colour factors have been hard-wired in the parametrizations.
*    The harmonic sums S_i(N) and the lowest integer values of the Zeta 
*    function are provided by the common-blocks  HSUMS  and  RZETA. 
*
* =====================================================================
*
*
       SUBROUTINE PNS2MOM 
*
       IMPLICIT DOUBLE COMPLEX (A - Z)
       INTEGER NMAX, NDIM, NFMIN, NFMAX, IAPP2, KN, NF
       PARAMETER (NDIM = 144, NFMIN = 3, NFMAX = 6)
       DOUBLE PRECISION ZETA(6)
*
* ---------------------------------------------------------------------
*
* ..Input common-blocks
*
       COMMON / MOMS   / NA (NDIM)
       COMMON / NNUSED / NMAX
       COMMON / HSUMS  / S(NDIM,6)
       COMMON / RZETA  / ZETA
*
* ..Output common-block
*
       COMMON / PNS2   / P2NS (NDIM, NFMIN:NFMAX, 3)
*
* ---------------------------------------------------------------------
*
* ..Begin of the Mellin-N loop
*
       DO 1 KN = 1, NMAX
*
* ..Some abbreviations and the simple harmonic sums for complex N
*
       N  = NA(KN) 
       S1 = S(KN,1)
       S2 = S(KN,2)
       S3 = S(KN,3)
*
       NI = 1./N
       NI2 = NI*NI
       NI3 = NI*NI2
       NM = N - 1.
       NMI = 1./NM
*
       N1 = N + 1.
       N1I = 1./N1
       N1I2 = N1I*N1I
       N1I3 = N1I*N1I2
       N2 = N + 2.
       N2I = 1./N2
*
       S1M = S1 - NI
       S11 = S1 + N1I
       S12 = S11 + N2I
*
* ---------------------------------------------------------------------
*
* ..The moments of the functions employed in the parametrizations:
*
* ...1/(1-x)_+ [A0]  and  x^a ln^b (1-x) [B`b(a)' with M for a = -1]
*
       A0  = - S1M
       B1  = - S1 * NI
* ...with special care for the first moment of x^-1 ln(1-x)
       IF ( ( DABS (DIMAG(N)) .LT. 1.D-5 ) .AND.
     ,      ( DABS ( DBLE(N) - 1.D0 ) .LT. 1.D-5 ) ) THEN
         B1M = - ZETA(2)
       ELSE
         B1M = - S1M * NMI
       ENDIF
       B11 = - S11 * N1I
       B12 = - S12 * N2I
*
* ...x^a [C`a']
*
       C0 = NI
       C1 = N1I
       C2 = N2I
       C3 = 1./(N+3.)
       C4 = 1./(N+4.)
*
* ...x^a ln^b x [D`b(a)']
*
       D1  = - NI2
       D2  = 2.* NI3
       D3  = - 6.* NI2*NI2
       D31 = - 6.* N1I2*N1I2
       D4  = 24.* NI2*NI3
*
* ...x^a ln^b x ln(1-x) [E`b(a)']
*
       E1  = S1*NI2 + (S2-ZETA(2))*NI
       E2  = 2.* ( - S1*NI3 + (ZETA(2)-S2)*NI2 - (S3-ZETA(3))*NI )
*
* ---------------------------------------------------------------------
*
* ..The parametrized n_f^{0,1} components of P_ns^(2)i, i = +, -, v
*    as given in S. Moch, J. Vermaseren and A. Vogt, hep-ph/0403192
*
       PP20 = + 1174.898 * A0 + 1295.384 + 714.1 * B1 - 522.1 * C3
     1        + 243.6 * C2 - 3135.* C1 + 1641.1 * C0 + 1258.* D1
     2        + 294.9 * D2 + 800/27.D0 * D3 + 128/81.D0 * D4
     3        + 563.9 * E1 + 256.8 * E2
       PP21 = - 183.187 * A0 - 173.927 - 5120/81.D0 * B1
     1        + 44.79 * C3 + 72.94 * C2 + 381.1 * C1 - 197.0 * C0
     2        - 152.6 * D1 - 2608./81.D0 * D2 - 192./81.D0 * D3
     3        - 56.66 * E1 - 1.497 * D31 
*
       PM20 = + 1174.898 * A0 + 1295.470 + 714.1 * B1 - 433.2 * C3
     1        + 297.0 * C2 - 3505.*C1 + 1860.2 * C0 + 1465.2 * D1
     2        + 399.2 * D2 + 320./9.D0 * D3 + 116./81.D0 * D4
     3        + 684.0 * E1 + 251.2 * E2
       PM21 = - 183.187 * A0 - 173.933 - 5120/81.D0 * B1
     1        + 34.76 * C3 + 77.89 * C2 + 406.5 * C1 - 216.62 * C0
     2        - 172.69 * D1 - 3216./81D0 * D2 - 256./81.D0 * D3
     3        - 65.43 * E1 - 1.136 * D31
*
       PS2  = - 163.9 * (B1M-B1)-7.208 * (B11-B12) + 4.82 * (C3-C4)
     1        - 43.12 * (C2-C3) + 44.51 * (C1-C2) + 151.49 * (C0-C1)
     2        + 178.04 * D1 + 6.892 * D2 - 40./27.D0 * (2.*D3 - D4)
     2        - 173.1 * E1 + 46.18 * E2
*
* ..The exact n_f^2 contribution first determined by J.A. Gracey in
*    hep-ph/9401214
*
       PF2  = - ( 17./72.D0 - 2./27.D0 * S1 - 10./27.D0 * S2
     1        + 2./9.D0 * S3 - (12.* N**4 + 2.* N**3 - 12.* N**2
     2        - 2.* N + 3.)/(27.* N**3 * N1**3) ) * 32./3.D0
*
* ---------------------------------------------------------------------
*
* ..Flavour-number loop and output to the array 
*
       DO 2 NF = NFMIN, NFMAX
*
         P2NS(KN,NF,1) = PP20 + NF * (PP21 + NF * PF2)
         P2NS(KN,NF,2) = PM20 + NF * (PM21 + NF * PF2)
         P2NS(KN,NF,3) = P2NS(KN,NF,2) + NF * PS2 
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
