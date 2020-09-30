C...ANOMALOUS DIMENSIONS FOR LEADING AND NEXT TO LEADING ORDER
C...EVOLUTION OF PARTON DENSITIES AND WILSON COEFFICIENTS FOR
C...NLO STRUCTURE FUNCTIONS :
       SUBROUTINE ANOM (ANS, AM, AP, AL, BE, AB, RMIN, RPLUS, RQQ, RQG,
     1                  RGQ, RGG, C2Q, C2G, CDYQ, CDYG, XN, FR,
     2                  QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F,
     3                  QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, C2QI, C2GF,
     4                  CDYQI, CDYGI)
       IMPLICIT DOUBLE COMPLEX (A - Z)
       DOUBLE PRECISION B0, B1, B02, B10
       INTEGER N, F, FR, nnf
       COMMON /ANOMSCHEME/ N
C...ANOMALOUS DIMENSIONS AND RELATED QUANTITIES IN LEADING ORDER :
      COMMON/NFLAVORS/nnF
       F = FR
       if (nnf.eq.3) F = 3
       N = 1
       B0 = 11.- 2./3.* FR
       B02 = 2.* B0
       QQ = QQI
       QG = F * QGF
       GQ = GQI
       GG = GGI + F * GGF
       SQ = SQRT ((GG - QQ) * (GG - QQ) + 4.* QG * GQ)
       GP = 0.5 * (QQ + GG + SQ)
       GM = 0.5 * (QQ + GG - SQ)
       ANS = QQ / B02
       AM = GM / B02
       AP = GP / B02
       AL = (QQ - GP) / (GM - GP)
       BE = QG / (GM - GP)
       AB = GQ / (GM - GP)
C...NEXT TO LEADING ORDER : ANOMALOUS DIMENSIONS AND WILSON COEFFICIENTS
C...IN THE MS-BAR FACTORIZATION SCHEME OF BARDEEN ET AL. (1981) :
       NS1M = NS1MI + F * NS1F
       NS1P = NS1PI + F * NS1F
       QQ1 = NS1P + F * QQ1F
       QG1 = F * QG1F
       GQ1 = GQ1I + F * GQ1F
       GG1 = GG1I + F * GG1F
       C2Q = C2QI
       C2G = F * C2GF
       CDYQ = CDYQI
       CDYG = CDYGI
C...CHANGE TO THE SCHEME OF ALTARELLI ET AL. (1979) FOR N = 2 OR TO THE
C...DIS SCHEME FOR N = 3 :
       IF (N .EQ. 2) THEN
          DEL = -0.5 * QG
          C2G = C2G + DEL
          QQ1 = QQ1 - DEL * GQ
          QG1 = QG1 + DEL * (QQ - GG - B02 - QG)
          GQ1 = GQ1 + DEL * GQ
          GG1 = GG1 + DEL * (GQ + B02)
       ELSE IF (N .EQ. 3) THEN
          NS1P = NS1P + B02 * C2Q
          NS1M = NS1M + B02 * C2Q
          QQ1 = QQ1 + C2Q * (QG + B02) + C2G * GQ
          QG1 = QG1 + C2Q * QG + C2G * (GG - QQ + QG + B02)
          GQ1 = GQ1 - C2Q * (QQ - GG + GQ + B02) - C2G * GQ
          GG1 = GG1 - C2Q * QG - C2G * (GQ + B02)
          C2Q = 0.
          C2G = 0.
       END IF
       XN1 = XN + 1.
C       C3Q = C2Q - 8./3.* (1./ XN + 1./ XN1)
C       CLQ = 16./ (3.* XN1)
C       CLG = 8.* F / (XN1 * (XN + 2.))
C...COMBINATIONS OF ANOMALOUS DIMENSIONS FOR NLO - SINGLET EVOLUTION :
       B1 = 102 - 38./3.* FR
       B10 = B1 / B0
       RMIN = (NS1M - QQ * B10) / B02
       RPLUS = (NS1P - QQ * B10) / B02
       RQQ = (QQ1 - QQ * B10) / B02
       RQG = (QG1 - QG * B10) / B02
       RGQ = (GQ1 - GQ * B10) / B02
       RGG = (GG1 - GG * B10) / B02
       RETURN
       END
C
C
C...CALCULATION OF ANOMALOUS DIMENSIONS AND WILSON COEFFICIENTS
C...UP TO THEIR DEPENDENCE OF THE NUMBER OF ACTIVE FLAVOURS F :
       SUBROUTINE ANCALC (QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F,
     1                    QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, XN)
       IMPLICIT DOUBLE COMPLEX (A - Z)
       DOUBLE PRECISION ZETA2, ZETA3,  Q2, Q2MUR, Q2MUF
c       COMMON / SCALES / Q2, Q2MUR, Q2MUF
C
c       LMQ = LOG(Q2/Q2MUF)
       XNS = XN * XN
       XN1 = XN + 1.
       XN2 = XN + 2.
       XNM = XN - 1.
C...LEADING ORDER :
       CPSI = PSIFN (XN1) + 0.577216
       QQI = (8./3.) * (-3.- 2./(XN * XN1) + 4.* CPSI)
       QGF = -4.* (XNS + XN +2.) / (XN * XN1 * XN2)
       GQI = -(16./3.) * (XNS + XN + 2.) / (XN * XN1 * XNM)
       GGI = -22.- 24./(XN * XNM) - 24./(XN1 * XN2) + 24.* CPSI
       GGF = 4./3.
C...NEXT OT LEADING ORDER :
       XNT = XNS * XN
       XNFO = XNT * XN
       XN1S = XN1 * XN1
       XN1T = XN1S * XN1
C...ANALYTIC CONTINUATIONS OF N-SUMS AS GIVEN IN GLUECK ET AL. (1990) :
       ZETA2 = 1.644934
       ZETA3 = 1.202057
       CPSI1 = ZETA2 - PSIFN1 (XN1)
*       CPSI2 = 0.5 * PSIFN2 (XN1) - ZETA2
*       SPMOM = 1.01/XN1 - 0.846/XN2 + 1.155/(XN+3.) - 1.074/(XN+4.) +
*     1         0.55/(XN+5.)
c       SPMOM = 1.004D0 / XN1 - 0.846D0 / XN2 + 1.342D0 / (XN+3.) -
c     1         1.532D0 / (XN+4.) + 0.839D0 / (XN+5.)
       SPMOM=ACG3(XN-1D0) ! --> Improved approximation for the Mellin transform of: Li2(x)/(1+x)
*
       SLC = -5./8.* ZETA3
       SLV = - ZETA2/2.* (PSIFN (XN1/2.) - PSIFN (XN/2.))
     1       + CPSI/XNS + SPMOM
       SSCHLM = SLC - SLV
       SSTR2M = ZETA2 - PSIFN1 (XN1/2.)
       SSTR3M = 0.5 * PSIFN2 (XN1/2.) + ZETA3
       SSCHLP = SLC + SLV
       SSTR2P = ZETA2 - PSIFN1 (XN2/2.)
       SSTR3P = 0.5 * PSIFN2 (XN2/2.) + ZETA3
C...NON-SINGLET PIECES AS GIVEN IN CURCI ET AL. (1980) :
       NS1MA = 16.* CPSI * (2.* XN + 1.) / (XNS * XN1S) +
     1         16.* (2.* CPSI - 1./(XN * XN1)) * ( CPSI1 - SSTR2M ) +
     2         64.* SSCHLM + 24.* CPSI1 - 3. - 8.* SSTR3M -
     3         8.* (3.* XNT + XNS -1.) / (XNT * XN1T) +
     4         16.* (2.* XNS + 2.* XN +1.) / (XNT * XN1T)
       NS1PA = 16.* CPSI * (2.* XN + 1.) / (XNS * XN1S) +
     1         16.* (2.* CPSI - 1./(XN * XN1)) * ( CPSI1 - SSTR2P ) +
     2         64.* SSCHLP + 24.* CPSI1 - 3. - 8.* SSTR3P -
     3         8.* (3.* XNT + XNS -1.) / (XNT * XN1T) -
     4         16.* (2.* XNS + 2.* XN +1.) / (XNT * XN1T)
       NS1B = CPSI * (536./9. + 8.* (2.* XN + 1.) / (XNS * XN1S)) -
     1        (16.* CPSI + 52./3.- 8./(XN * XN1)) * CPSI1 - 43./6. -
     2        (151.* XNFO + 263.* XNT + 97.* XNS + 3.* XN + 9.) *
     3        4./ (9.* XNT * XN1T)
       NS1C = -160./9.* CPSI + 32./3.* CPSI1 + 4./3. +
     1        16.* (11.* XNS + 5.* XN - 3.) / (9.* XNS * XN1S)
       NS1MI = -2./9.* NS1MA + 4.* NS1B
       NS1PI = -2./9.* NS1PA + 4.* NS1B
       NS1F = 2./3. * NS1C
C...SINGLET PIECES AS GIVEN IN FLORATOS ET AL. (1981) :
       XNFI = XNFO * XN
       XNSI = XNFI * XN
       XNSE = XNSI * XN
       XNE = XNSE * XN
       XNN = XNE * XN
       XNMS = XNM * XNM
       XN2S = XN2 * XN2
       XN2T = XN2S * XN2
       QQ1F = (5.* XNFI + 32.* XNFO + 49.* XNT + 38.* XNS + 28.* XN
     1          + 8.) / (XNM * XNT * XN1T * XN2S) * (-32./3.)
       QG1A = (-2.* CPSI * CPSI + 2.* CPSI1 - 2.* SSTR2P)
     1          * (XNS + XN + 2.) / (XN * XN1 * XN2)
     2        + (8.* CPSI * (2.* XN + 3.)) / (XN1S * XN2S)
     3        + 2.* (XNN + 6.* XNE + 15. * XNSE + 25.* XNSI + 36.* XNFI
     4          + 85.* XNFO + 128.* XNT + 104.* XNS + 64.* XN + 16.)
     5          / (XNM * XNT * XN1T * XN2T)
       QG1B = (2.* CPSI * CPSI - 2.* CPSI1 + 5.) * (XNS + XN + 2.)
     1          / (XN * XN1 * XN2)   -   4.* CPSI / XNS
     2        + (11.* XNFO + 26.* XNT + 15.* XNS + 8.* XN + 4.)
     3          / (XNT * XN1T * XN2)
       QG1F = - 12.* QG1A - 16./3.* QG1B
       GQ1A = (-2.* CPSI * CPSI + 10.* CPSI - 2.* CPSI1)
     1          * (XNS + XN + 2.) / (XNM * XN * XN1)  -  4.* CPSI / XN1S
     2        - (12.* XNSI + 30.* XNFI + 43.* XNFO + 28.* XNT - XNS
     3          - 12.* XN - 4.) / (XNM * XNT * XN1T)
       GQ1B = (CPSI * CPSI + CPSI1 - SSTR2P) * (XNS + XN + 2.)
     1          / (XNM * XN * XN1)
     2        - CPSI * (17.* XNFO + 41.* XNS - 22.* XN - 12.)
     3          / (3.* XNMS * XNS * XN1)
     4        + (109.* XNN + 621.* XNE + 1400.* XNSE + 1678.* XNSI
     5          + 695.* XNFI - 1031.* XNFO - 1304.* XNT - 152.* XNS
     6          + 432.* XN + 144.) / (9.* XNMS * XNT * XN1T * XN2S)
       GQ1C = (CPSI - 8./3.) * (XNS + XN + 2.) / (XNM * XN * XN1)
     1        + 1./ XN1S
       GQ1I = - 64./9.* GQ1A - 32.* GQ1B
       GQ1F = - 64./9.* GQ1C
       GG1A = 16./9.* (38.* XNFO + 76.* XNT + 94.* XNS + 56.* XN + 12.)
     1          / (XNM * XNS * XN1S * XN2)   -   160./9.* CPSI + 32./3.
       GG1B = (2.* XNSI + 4.* XNFI + XNFO - 10.* XNT - 5.* XNS - 4.* XN
     1          - 4.) * 16. / (XNM * XNT * XN1T * XN2)   +   8.
       GG1C = (2.* XNFI + 5.* XNFO + 8.* XNT + 7.* XNS - 2.* XN - 2.)
     1          * 64.* CPSI / (XNMS * XNS * XN1S * XN2S)
     2        + 536./9.* CPSI - 64./3.
     3        + 32.* SSTR2P * (XNS + XN + 1.) / (XNM * XN * XN1 * XN2)
     4        - 16.* CPSI * SSTR2P + 32.* SSCHLP - 4.* SSTR3P
     5        - 4.* (457.* XNN + 2742.* XNE + 6040.* XNSE + 6098.* XNSI
     6          + 1567.* XNFI - 2344.* XNFO - 1632.* XNT + 560.* XNS
     7          + 1488.* XN + 576.) / (9.* XNMS * XNT * XN1T * XN2T)
       GG1I = 9.* GG1C
       GG1F = 3./2.* GG1A + 2./3.* GG1B

C...WILSON COEFFICIENTS :
c       C2QI = 4./3.* (2.* CPSI * CPSI - 2.* CPSI1 + 3.* CPSI - 9.
c     1       - 2.* CPSI / (XN * XN1) + 3./ XN + 4./ XN1 + 2./ XNS)
c       C2QI = C2QI - LMQ * QQI/2.      
c       C2GF = - 2.* (CPSI * (XNS + XN + 2.) / (XN * XN1 * XN2)
c     1       + 1./ XN - 1./ XNS - 6./ XN1 + 6./ XN2)
c       C2GF = C2GF - LMQ * QGF/2.    
C... DRELL-YAN COEFFICIENTS :
c       CDYQI = 4./3. * (-8. + 8.*ZETA2 + 2./XNS + 2./XN1S - 
c     1                   4.*CPSI/XN/XN1 + 4.*CPSI*CPSI +
c     2                   (3. + 2./XN/XN1 - 4.*CPSI)*LMQ )
c       CDYGI = 1./2. * ( (4. + 14.*XN + 22.*XNS + 11.*XNT + XNFO)/
c     1                   (XNS*XN1S*XN2S) - 
c     2                    2.*(2. + XN + XNS)*CPSI/(XN*XN1*XN2) +
c     3                   (1./XN - 2./XN1 + 2./XN2)*LMQ )       
       RETURN
       END
C
C

C
C
C...PSI - FUNCTION FOR COMPLEX ARGUMENT
       DOUBLE COMPLEX FUNCTION PSIFN (Z)
       DOUBLE COMPLEX Z, ZZ, RZ, DZ, SUB
       SUB = CMPLX (0.,0.)
       ZZ = Z
  1    CONTINUE
       IF (DREAL (ZZ) .LT. 10.) THEN
         SUB = SUB - 1./ ZZ
         ZZ = ZZ + 1.
         GOTO 1
       END IF
       RZ = 1./ ZZ
       DZ = RZ * RZ
       PSIFN = SUB + LOG(ZZ) - RZ/2.- DZ/2520. * ( 210.+ DZ * (-21.+
     1         10.*DZ ))
       RETURN
       END
C
C
C...FIRST DERIVATIVE OF THE PSI - FUNCTION FOR COMPLEX ARGUMENT :
       DOUBLE COMPLEX FUNCTION PSIFN1 (Z)
       DOUBLE COMPLEX Z, ZZ, RZ, DZ, SUB
       SUB = CMPLX (0.,0.)
       ZZ = Z
  1    CONTINUE
       IF (DREAL (ZZ) .LT. 10.) THEN
         SUB = SUB + 1./ (ZZ * ZZ)
         ZZ = ZZ + 1.
         GOTO 1
       END IF
       RZ = 1./ ZZ
       DZ = RZ * RZ
       PSIFN1 = SUB + RZ + DZ/2. * ( 1 + RZ/630. * ( 210.- DZ * ( 42.-
     1         DZ * ( 30.- 42.*DZ ))))
       RETURN
       END
C
C
C...SECOND DERIVATIVE OF THE PSI - FUNCTION FOR COMPLEX ARGUMENT :
       DOUBLE COMPLEX FUNCTION PSIFN2 (Z)
       DOUBLE COMPLEX Z, ZZ, RZ, DZ, SUB
       SUB = CMPLX (0.,0.)
       ZZ = Z
  1    CONTINUE
       IF (DREAL (ZZ) .LT. 10.) THEN
         SUB = SUB - 2./ (ZZ * ZZ * ZZ)
         ZZ = ZZ + 1.
         GOTO 1
       END IF
       RZ = 1./ ZZ
       DZ = RZ * RZ
       PSIFN2 = SUB - DZ/60. * ( 60.+ RZ * ( 60.+ RZ * ( 30.- DZ *
     1         ( 10.- DZ * ( 10.- DZ * ( 18.- 50.* DZ ))))))
       RETURN
       END
