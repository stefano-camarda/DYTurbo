*
* ..File: pns2old.f
*
*
* ..The subroutine  PNS2MOM  for the approximate NNLO non-singlet QCD 
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
* ..The two approximations spanning the estimated residual uncertainty 
*    are invoked by  IAPP2 = 1  or  2,  for any other value of IAPP2
*    (input common-block  P2APPR)  the central results are returned.  
*
* ..The QCD colour factors have been hard-wired in the approximations.
*    The harmonic sums S_i(N) are provided by the common-block  HSUMS. 
*
*
* ..Reference: W.L. van Neerven and A. Vogt, hep-ph/0007362
*
* =====================================================================
*
*
       SUBROUTINE PNS2MOM 
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
       N1 = N + 1.
       NI = 1./N
       N1I = 1./N1
       N2I = 1./(N + 2.)
       N3I = 1./(N + 3.)
       N4I = 1./(N + 4.)
*
       S1M = S1 - NI
       S11 = S1 + N1I
       S21 = S1 + N1I*N1I
*
* ---------------------------------------------------------------------
*
*  ..Moments of the basic x-space functions 
*
       A0  = - S1M 
       B1  = - S1 * NI
       B11 = - S11 * N1I
       B2  = (S1**2 + S2) * NI
*
       C0  = NI
       C1  = N1I
       C2  = N2I
       C3  = N3I
       C4  = N4I
*
       D1  = - NI*NI
       D12 = - N2I * N2I
       D2  = 2.* NI**3
       D21 = 2.* N1I**3
       D3  = - 6.* NI**4
       D4  = 24.* NI**5
*
* ---------------------------------------------------------------------
*
* ..The approximate contributions to P2NS 
*
* ..`+' case
*
       PPA = + 1183.762 * A0 + 1347.032 + 1047.590 * B1 - 843.884 * C2
     1       - 98.65 * (C0 - C1) - 33.71 * D2 + 1.580 * (D4 + 4.* D3)
       PPB = + 1182.774 * A0 + 1351.088 - 147.692 * B2 - 2602.738 * C2
     1       - 170.11 * C0 + 148.47 * D1 + 1.580 * (D4 - 4.* D3)
*
       PP1A = - 183.148 * A0 - 174.402 + 9.649 * B2 + 406.171 * C2
     1        + 32.218 * (C0 - C1) + 5.976 * D2 + 1.60 * D3 
       PP1B = - 183.931 * A0 - 178.208 - 89.941 * B1 + 218.482 * C2
     1        + 9.623 * C0 + 0.910 * D2 - 1.60 * D3 
*
* ..`-' case
*
       PMA = + 1185.229 * A0 + 1365.458 - 157.387 * B2 - 2741.42 * C2
     1       - 490.43 * (C0 - C1) + 67.00 * D2 + 10.005 * D3 
     2       + 1.432 * D4
       PMB = + 1174.348 * A0 + 1286.799 + 115.099 * B2 + 1581.05 * B1
     1       + 267.33 * (C0 - C1) - 127.65 * D2 - 25.22 * D3
     2       + 1.432 * D4
*
       PM1A = - 184.765 * A0 - 184.289 + 17.989 * B2 + 355.636 * C2
     1        - 73.407 * (B1 - B11) + 11.491 * D2 + 1.928 * D3 
       PM1B = - 183.718 * A0 - 177.762 + 11.999 * B2 + 397.546 * C2
     1        + 41.949 * (C0 - C1) - 1.477 * D2 - 0.538 * D3 
*
* ---------------------------------------------------------------------
*
* ..Difference of `V' and `-'
*
       PSA = - 1441.57 * (C2 - C3) + 12603.59 * (C1 - C2)
     1       - 15450.01 * (C0 - C1) + 7876.93 * D21 - 4260.29 * D1
     2       - 229.27 * D2 + 4.4075 * D3 
       PSB = - 704.67 * (C3 - C4) + 3310.32 * (C2 - C3)
     1       + 2144.81 * (C1 - C2) - 244.68 * (C0 - C1) + 4490.81 * D12
     2       + 42.875 * D1 - 11.0165 * D3 
*
* ..The exact NF^2 contribution of Gracey (1994) (CF = 4/3 inserted)
*
       PNSF2 = ( 17./72.D0 - 2./27.D0 * S1 - 10./27.D0 * S2
     1         + 2./9.D0 * S3 - (12.* N**4 + 2.* N**3 - 12.* N**2
     2         - 2.* N + 3.) / (3.* N * N1)**3 ) * (-32./3.D0)
*
* ---------------------------------------------------------------------
*
* ..Flavour-number loop and output to the array (depending on IAPP2)
*
       DO 2 NF = NFMIN, NFMAX
*
       IF (IAPP2 .EQ. 1) THEN
         P2NS(KN,NF,1) = PPA + NF * PP1A + NF**2 * PNSF2 
         P2NS(KN,NF,2) = PMA + NF * PM1A + NF**2 * PNSF2  
         P2NS(KN,NF,3) = PMA + NF * (PM1A+PSA) + NF**2 * PNSF2 
 
       ELSE IF (IAPP2 .EQ. 2) THEN
         P2NS(KN,NF,1) = PPB + NF * PP1B + NF**2 * PNSF2 
         P2NS(KN,NF,2) = PMB + NF * PM1B + NF**2 * PNSF2  
         P2NS(KN,NF,3) = PMB + NF * (PM1B+PSB) + NF**2 * PNSF2 
*
       ELSE
         P2NS(KN,NF,1) = 0.5 * (PPA+PPB + NF * (PP1A+PP1B)) 
     1                   + NF**2 * PNSF2
         P2NS(KN,NF,2) = 0.5 * (PMA+PMB + NF * (PM1A+PM1B)) 
     1                   + NF**2 * PNSF2
         P2NS(KN,NF,3) = 0.5 * (PMA+PMB + NF * (PM1A+PM1B+PSA+PSB)) 
     1                   + NF**2 * PNSF2
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
