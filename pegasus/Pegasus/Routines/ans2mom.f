*
* ..File: ans2mom.f     
*
*
* ..The subroutine  ANS2MOM  for the NNLO (alpha_s^2) heavy quark 
*    contribution  A2NS  to the non-singlet operator matrix element 
*    (OME) in N-space in the MS(bar) scheme for mu_f^2 = m_H^2.
*    The coupling constant is normalized as  a_s = alpha_s/(4*pi).
*
* ..This quantity, presented in Appendix B of Buza, Matiounine, Smith 
*    and van Neerven, Eur. Phys. J. C1 (1998) 301 (BSMN), is required 
*    for the N_f matching of the NNLO parton densities.
*
* ..The results (written to the common-block  ANS2)  are computed on an 
*    external NDIM-dimensional array  NA  of complex Mellin moments 
*    provided by the common-block  MOMS. 
*
* ..The SU(N_colours=3) colour factors  CF, CA and TF  are taken from
*    the common-block  COLOUR.  The simple harmonic sums S_i(N) are
*    provided by the common-block  HSUMS,  and the lowest integer
*    values of the Riemann Zeta-function are provided by  RZETA.
*
* =====================================================================
*
*
       SUBROUTINE ANS2MOM 
*
       IMPLICIT DOUBLE COMPLEX (A - Z)
       INTEGER NMAX, NDIM, KN
       PARAMETER (NDIM = 144)
       DOUBLE PRECISION ZETA(6), CF, CA, TR
*
* ---------------------------------------------------------------------
*
* ..Input common-blocks 
*
       COMMON / MOMS   / NA (NDIM)
       COMMON / NNUSED / NMAX
       COMMON / HSUMS  / S(NDIM,6)
       COMMON / RZETA  / ZETA
       COMMON / COLOUR / CF, CA, TR
*
* ..Output common-block 
*
       COMMON / ANS2   / A2NS (NDIM)
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
*
       N1 = N + 1.
       NI = 1./N
       N1I = 1./N1
*
       S1M = S1 - NI
       S2M = S2 - NI*NI
       S3M = S3 - NI**3
       S21 = S2 + N1I*N1I
       S31 = S3 + N1I**3
*
* ---------------------------------------------------------------------
*
*  ..Moments of the basic x-space functions 
*
       A0 = - S1M
*
       C0 = NI
       C1 = N1I
*
       D1  = - NI*NI
       D11 = - N1I*N1I
*
       G1  = S2M - ZETA(2)  
       G12 = S21 - ZETA(2)  
       G2  = - 2.* ( S3M - ZETA(3) )  
       G22 = - 2.* ( S31 - ZETA(3) )  
*
* ---------------------------------------------------------------------
*
* ..The moments of the OME A_{qq,H}^{NS,(2)} given in Eq. (B.4) of BMSN 
*
       A2QQ = 224./27.D0 * A0 - 8./3.D0 * ZETA(3) + 40/9.D0 * ZETA(2) 
     1        + 73./18.D0 + 44./27.D0 * C0 - 268./27.D0 * C1 
     2        + 8./3.D0 * (D1 - D11) + 20./9.D0 * (G1 + G12) 
     3        + 2./3.D0 * (G2 + G22)
*
* ..Output to the array 
*
       A2NS(KN) = CF*TR * A2QQ
*
* ---------------------------------------------------------------------
*
  1    CONTINUE
*
       RETURN
       END
*
* =================================================================av==
