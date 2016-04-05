*
* ..File: psg2mom.f     (requires previous call of PNS2MOM)
*
*
* ..The routine  PSG2MOM  for the parametrized NNLO (alpha_s^3) singlet 
*    QCD splitting functions  P2NS  in N-space for  NF = NFMIN .. NFMAX
*    = 3...6  massless quark flavours in the MS(bar) scheme.
*    The coupling constant is normalized as  a_s = alpha_s/(4*pi). 
*
* ..These quantities are determined on an external NDIM-dimensional 
*    array  NA  of complex Mellin moments provided by the common-block 
*    MOMS.  The results are written to the common-block  PSG2. The last
*    two array arguments refer to the quark-gluon mixing matrix with
*    1 = q  and  2 = g.
*
* ..The QCD colour factors have been hard-wired in the approximations.
*    The harmonic sums S_i(N) and the lowest integer values of the Zeta
*    function are provided by the common-blocks  HSUMS  and  RZETA.
*
* =====================================================================
*
*
       SUBROUTINE PSG2MOM 
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
       NI = 1./N
       NI2 = NI*NI
       NI3 = NI*NI2
       NM = N - 1.
       NMI = 1./NM
       NMI2 = NMI*NMI
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
       S2M = S2 - NI2
       S21 = S2 + N1I2
       S31 = S3 + N1I3
*
* ---------------------------------------------------------------------
*
* ..The moments of the functions employed in the parametrizations:
*
* ...1/(1-x)_+ [A0]  and  x^a ln^b (1-x) [B`b(a)' with M for a = -1]
*
       A0  = - S1M
       B1  = - S1 * NI
       B1M = - S1M * NMI
       B11 = - S11 * N1I
       B2  = (S1**2 + S2) * NI
       B2M = (S1M**2 + S2M) * NMI
       B21 = (S11**2 + S21) * N1I
       B3  = - (S1**3 + 3.*S1*S2 + 2.*S3) * NI
       B31 = - (S11**3 + 3.*S11*S21 + 2.*S31) * N1I
       B4  = (S1**4 + 6.*S1**2*S2 + 8.*S1*S3 + 3.*S2**2 + 6.*S4) * NI
*
* ...x^a [C`a']
*
       C0 = NI
       CM = NMI
       C1 = N1I
       C2 = N2I
       C3 = 1./(N+3.)
       C4 = 1./(N+4.)
*
* ...x^a ln^b x [D`b(a)']
*
       D1  = - NI2
       D1M = - NMI2
       D11 = - N1I2
       D2  = 2.* NI3
       D21 = 2.* N1I3
       D3  = - 6.* NI2*NI2
       D31 = - 6.* N1I2*N1I2
       D4  = 24.* NI2*NI3
       D41 = 24.* N1I2*N1I3
*
* ...x^a ln^b x ln(1-x) [E`b(a)']
*
       E1  = S1*NI2 + (S2-ZETA(2))*NI
       E11 = S11*N1I2 + (S21-ZETA(2))*N1I
       E2  = 2.* ( - S1*NI3 + (ZETA(2)-S2)*NI2 - (S3-ZETA(3))*NI )
*
* ---------------------------------------------------------------------
*
* ..The parametrized n_f-components of P_ps^(2) and P_qg^(2) ...
*   [ P_qq^(2) is obtained below by adding the non-singlet quantity
*     P_ns^(2)+ provided by the subroutine  P2NS2MOM  to P_ps^(2) ]
*
       PS1 = - 3584./27.* (D1M-D1) - 506.* (CM-C0) + 160./27.* (D4-D41)
     ,       - 400./9.* (D3-D31) + 131.4 * (D2-D21) - 661.6 * (D1-D11)
     ,       - 5.926 * (B3-B31) - 9.751 * (B2-B21) - 72.11 * (B1-B11)
     ,       + 177.4 * (C0-C1) + 392.9 * (C1-C2) - 101.4 * (C2-C3)
     ,       - 57.04 * (E1-E11)
       PS2 = 256./81.* (CM-C0) + 32./27.* (D3-D31) + 17.89 * (D2-D21)
     ,       + 61.75 * (D1-D11) + 1.778 * (B2-B21) + 5.944 * (B1-B11)
     ,       + 100.1 * (C0-C1) - 125.2 * (C1-C2) + 49.26 * (C2-C3)
     ,       - 12.59 * (C3-C4) - 1.889 * (E1-E11)
*
       QG1 = - 896./3.* D1M - 1268.3 * CM + 536./27.* D4 - 44./3.* D3
     ,       + 881.5 * D2 + 424.9 * D1 + 100./27.* B4 - 70./9.* B3
     ,       - 120.5 * B2 + 104.42 * B1 + + 2522.* C0 - 3316. * C1
     ,       + 2126. * C2 + 1823. * E1 - 25.22 * E2 - 252.5 * D31
       QG2 =   1112./243.* CM - 16./9.* D4 - 376./27.* D3 - 90.8 * D2
     ,       - 254.0 * D1 + 20./27.* B3 + 200./27.* B2 - 5.496 * B1
     ,       - 252.0 * C0 + 158.0 * C1 + 145.4 * C2 - 139.28 * C3
     ,       - 53.09 * E1 - 80.616 * E2 - 98.07 * D21 + 11.70 * D31
*
* ...and of P^(2)_gq and P^(2)_gg  [GQ2 is exact], all as given by 
*     A. Vogt, S. Moch and J. Vermaseren in hep-ph/0404111
*
       GQ0 = 1189.3 * D1M + 6163.1 * CM - 4288./81. * D4 + 1568./9.* D3
     ,       - 1794.* D2 + 4033.* D1 + 400./81.* B4 + 2200./27.* B3
     ,       + 606.3 * B2 + 2193.* B1 - 4307.* C0 + 489.3 * C1
     ,       + 1452.* C2 + 146.* C3 - 447.3 * E2 - 972.9 * D21
       GQ1 =   71.082 * D1M - 46.41 * CM + 128./27.* D4 + 704/81.* D3
     ,       + 20.39 * D2 + 174.8 * D1 - 400./81.* B3 - 68.069 * B2
     ,       - 296.7 * B1 - 183.8 * C0 + 33.35 * C1 - 277.9 * C2
     ,       + 108.6 * D21 - 49.68 * E1
       GQ2 = ( 64.* (- CM + C0 + 2.* C1) + 320.* (B1M - B1 + 0.8 * B11)
     ,       + 96.* (B2M - B2 + 0.5 * B21) ) / 27.
*
       GG0 = 2675.8 * D1M + 14214.* CM - 144.D0 * D4 + 72.D0 * D3
     ,       - 7471.* D2 + 274.4 * D1 - 20852.* C0 + 3968.* C1
     ,       - 3363.* C2 + 4848.* C3 + 7305.* E1 + 8757.* E2
     ,       + 3589.* B1 + 4425.894 + 2643.521 * A0
       GG1 =   157.27 * D1M + 182.96 * CM + 512./27.D0 * D4
     ,       + 832./9.D0 * D3 + 491.3 * D2 + 1541.* D1 - 350.2 * C0
     ,       + 755.7 * C1 - 713.8 * C2 + 559.3 * C3 + 26.15 * E1
     ,       - 808.7 * E2 - 320.D0 * B1 - 528.723 - 412.172 * A0
       GG2 = - 680./243.D0 * CM - 32./27.D0 * D3 + 9.680 * D2
     ,       - 3.422 * D1 - 13.878 * C0 + 153.4 * C1 - 187.7 * C2
     ,       + 52.75 * C3 - 115.6 * E1 + 85.25 * E11 - 63.23 * E2
     ,       + 6.4630 - 16./9.D0 * A0
*
* ---------------------------------------------------------------------
*
* ..Flavour-number loop and output to the array 
*
       DO 2 NF = NFMIN, NFMAX
*
         P2SG(KN,NF,1,1) = P2NS(KN,NF,1) + NF * ( PS1 + NF * PS2 )
         P2SG(KN,NF,1,2) =       NF * ( QG1 + NF * QG2 )
         P2SG(KN,NF,2,1) = GQ0 + NF * ( GQ1 + NF * GQ2 )
         P2SG(KN,NF,2,2) = GG0 + NF * ( GG1 + NF * GG2 )
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
