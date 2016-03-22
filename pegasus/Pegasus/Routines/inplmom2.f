*
* ..File: inplmom2.f
*
*
* ..The subroutine  INPLMOM2  (version 2) for the initial light-parton 
*    distributions of a non-strange hadron in N-space.  The moments are 
*    stored in the common-block  PAINP for an external NDIM-dimensional
*    array  NA  of complex moments specified in the common-block  MOMS.
*
* ..The functional form of the distributions for this version reads
*
*        xf = Nf * x^Pf(2) * (1 - x)^Pf(3)  
*             * ( 1 + x^0.5 * Pf(4) + x * Pf(5) + x^1.5 * Pf(6))    (1)
*
*    for  f = UV = u - ubar,  DV = d - dbar,  DL = dbar - ubar, 
*    LS = 2 * (dbar + ubar),  SM = s - sbar,  SS = s + sbar,   GL = g.
*
* ..The normalization factors  NUV and NDV  are fixed by the respective 
*    quark numbers given by  PUV(1)  and  PDV(1)  (i.e.  PUV(1) = 2 and
*    PDV(1) = 1  for the proton).
*    NDL = PDL(1);  for  IMOMIN = 0  also  NLS = PLS(1), NSS = PSS(1).
*    For non-zero  IMOMIN,  PLS(1) and PSS(1)  represent the respective
*    momentum fractions.  NGL  is fixed by the momentum sum  PGL(1)  of
*    all partons.  For non-zero  ISSIMP,  SS = PSS(1) * LS  is used for 
*    the strange sea instead of (1), and  SM  is set to zero. Otherwise
*    the coefficient of x^PSM(4) is such that the first moment is zero.
*
* ..The normalization factors  Nf  and the momentum fractions  Af  are 
*    written to the common-block PANORM.  The flag  IINNEW  in  INPNEW
*    is set to '1' at the end of this routine.
*
* =====================================================================
*
*
       SUBROUTINE INPLMOM2 (PUV,PDV,PDL,PLS,PSM,PSS,PGL,IMOMIN,ISSIMP)
*
       IMPLICIT DOUBLE COMPLEX (A - Z)            
       INTEGER NMAX, NDIM, IMOMIN, ISSIMP, IINNEW, KN
       PARAMETER (NDIM = 144)
       PARAMETER ( E = (1.D0, 0.D0), ZERO = (0.D0, 0.D0) )
       DOUBLE PRECISION PUV(6), PDV(6), PDL(6), PLS(6), PSM(6), PSS(6),
     1         PGL(6), AUV, ADV, ALS, ASS, AGL, NUV, NDV, NLS, NSS, NGL
*
* ---------------------------------------------------------------------
*
* ..Input common-blocks
*
       COMMON / MOMS   / NA(NDIM)
       COMMON / NNUSED / NMAX
*
* ..Output common-blocks
*
       COMMON / PAINP  / VAI (NDIM), M3I (NDIM), M8I(NDIM), 
     1                   SGI (NDIM), P3I (NDIM), P8I(NDIM), GLI (NDIM)
       COMMON / PANORM / AUV, ADV, ALS, ASS, AGL, 
     1                   NUV, NDV, NLS, NSS, NGL
       COMMON / INPNEW / IINNEW
*
* ---------------------------------------------------------------------
*
* ..Normalizations of the valence distributions
*
       NUV = PUV(1)/ ( EBETA (PUV(2)*E, PUV(3)+E) *
     1               (1.+ PUV(5) * PUV(2) / (PUV(2)+ PUV(3)+ 1.))
     2               + EBETA (0.5+PUV(2), PUV(3)+E) * (PUV(4) + 
     3               PUV(6) * (PUV(2)+ 0.5) / (PUV(2)+ PUV(3)+ 1.5)) )
       NDV = PDV(1)/ ( EBETA (PDV(2)*E, PDV(3)+E) *
     1               (1.+ PDV(5) * PDV(2) / (PDV(2)+ PDV(3)+ 1.))
     2               + EBETA (0.5+PDV(2), PDV(3)+E) * (PDV(4) + 
     3               PDV(6) * (PDV(2)+ 0.5) / (PDV(2)+ PDV(3)+ 1.5)) )
*
* ..Second moments of the valence distributions (for the momentum sum) 
*
       AUV = NUV * ( EBETA (PUV(2)+E, PUV(3)+E) *
     1               (1.+ PUV(5) * (PUV(2)+1.) / (PUV(2) + PUV(3)+ 2.))
     2               + EBETA (1.5+PUV(2), PUV(3)+E) * (PUV(4) + 
     3               PUV(6) * (PUV(2)+ 1.5) / (PUV(2)+ PUV(3)+ 2.5)) )
       ADV = NDV * ( EBETA (PDV(2)+E, PDV(3)+E) *
     1               (1.+ PDV(5) * (PDV(2)+1.) / (PDV(2) + PDV(3)+ 2.))
     2               + EBETA (1.5+PUV(2), PUV(3)+E) * (PUV(4) + 
     3               PUV(6) * (PUV(2)+ 1.5) / (PUV(2)+ PUV(3)+ 2.5)) )
*
* --------------------------------------------------------------------- 
*
* ..Sea quark normalization factors NXS and momenta AXS, X = L, S
*
       IF (IMOMIN .EQ. 0) THEN
*
         NLS = PLS(1)
         ALS = NLS * ( EBETA (PLS(2)+E, PLS(3)+E) *
     1               (1.+ PLS(5) * (PLS(2)+ 1.) / (PLS(2)+ PLS(3)+ 2.))
     2             + EBETA (1.5+PLS(2), PLS(3)+E) * (PLS(4) + 
     3               PLS(6) * (PLS(2)+ 1.5) / (PLS(2)+ PLS(3)+ 2.5)) )
         IF (ISSIMP .EQ. 0) THEN
           NSS = PSS(1)
           ASS = NSS * ( EBETA (PSS(2)+E, PSS(3)+E) * 
     1               (1.+ PSS(5) * (PSS(2)+ 1.) / (PSS(2)+ PSS(3)+ 2.))
     2               + EBETA (1.5+PSS(2), PSS(3)+E) * (PSS(4) + 
     3               PSS(6) * (PSS(2)+ 1.5) / (PSS(2)+ PSS(3)+ 2.5)) )
         ELSE
           NSS = PSS(1) * NLS
           ASS = PSS(1) * ALS
         END IF
*
       ELSE
*
         ALS = PLS(1)
         NLS = ALS / ( EBETA (PLS(2)+E, PLS(3)+E) *
     1               (1.+ PLS(5) * (PLS(2)+ 1.) / (PLS(2)+ PLS(3)+ 2.))
     2               + EBETA (1.5+PLS(2), PLS(3)+E) * (PLS(4) + 
     3               PLS(6) * (PLS(2)+ 1.5) / (PLS(2)+ PLS(3)+ 2.5)) )
         IF (ISSIMP .EQ. 0) THEN
           ASS = PSS(1)
           NSS = ASS / ( EBETA (PSS(2)+E, PSS(3)+E) *
     1               (1.+ PSS(5) * (PSS(2)+ 1.) / (PSS(2)+ PSS(3)+ 2.))
     2               + EBETA (1.5+PSS(2), PSS(3)+E) * (PSS(4) + 
     3               PSS(6) * (PSS(2)+ 1.5) / (PSS(2)+ PSS(3)+ 2.5)) )
         ELSE
           NSS = PSS(1) * NLS
           ASS = PSS(1) * ALS
         END IF
* 
       END IF
*
* --------------------------------------------------------------------- 
*
* ..The vanishing first moment of s-s implemented by fixing PSM(5)
*
       IF (ISSIMP .EQ. 0) THEN
         PSM(5) = - (PSM(2)+ PSM(3)+ 1.)/PSM(2) 
     1            * (E + EBETA (0.5+PSM(2), PSM(3)+E) * (PSM(4) +
     2               PSM(6) * (PSM(2)+ 0.5) / (PSM(2)+ PSM(3)+ 1.5))
     3               / EBETA (PSM(2)*E, PSM(3)+E) )
       END IF
*
* ..The gluon normalization in terms of the total parton momentum
*      
       AGL = PGL(1) - AUV - ADV - ALS - ASS
       NGL = AGL / ( EBETA (PGL(2)+E, PGL(3)+E) * 
     1               (1.+ PGL(5) * (PGL(2)+ 1.) / (PGL(2)+ PGL(3)+ 2.))
     2             + EBETA (1.5+PGL(2), PGL(3)+E) * (PGL(4) + 
     3               PGL(6) * (PGL(2)+ 1.5) / (PGL(2)+ PGL(3)+ 2.5)) )
* 
* ---------------------------------------------------------------------
*
* ..Begin of the Mellin-N loop 
*
       DO 1 KN = 1, NMAX
       N = NA(KN)
*
* ..Valence quark distributions 

       UVI = NUV * ( EBETA (N-1.+PUV(2), PUV(3)+E) *
     1               (1.+ PUV(5) * (PUV(2)+N-1.) / (PUV(2)+PUV(3)+N))
     2             + EBETA (N-0.5+PUV(2), PUV(3)+E) * (PUV(4) + 
     3               PUV(6) * (PUV(2)+N-0.5) / (PUV(2)+PUV(3)+N+0.5)) )
*
       DVI = NDV * ( EBETA (N-1.+PDV(2), PDV(3)+E) *
     1               (1.+ PDV(5) * (PDV(2)+N-1.) / (PDV(2)+PDV(3)+N))
     2             + EBETA (N-0.5+PDV(2), PDV(3)+E) * (PDV(4) + 
     3               PDV(6) * (PDV(2)+N-0.5) / (PDV(2)+PDV(3)+N+0.5)) )
* 
* ---------------------------------------------------------------------
*
* ..Sea quark and gluon distributions  
*
       DLI = PDL(1) * ( EBETA (N-1.+PDL(2), PDL(3)+E) * 
     1               (1.+ PDL(5) * (PDL(2)+N-1.) / (PDL(2)+PDL(3)+N))
     2             + EBETA (N-0.5+PDL(2), PDL(3)+E) * (PDL(4) + 
     3               PDL(6) * (PDL(2)+N-0.5) / (PDL(2)+PDL(3)+N+0.5)) )
*
       LSI = NLS * ( EBETA (N-1.+PLS(2), PLS(3)+E) *
     1               (1.+ PLS(5) * (PLS(2)+N-1.) / (PLS(2)+PLS(3)+N))
     2             + EBETA (N-0.5+PLS(2), PLS(3)+E) * (PLS(4) + 
     3               PLS(6) * (PLS(2)+N-0.5) / (PLS(2)+PLS(3)+N+0.5)) )
*
       IF (ISSIMP .EQ. 0) THEN
         SMI = PSM(1) * ( EBETA (N-1.+PSM(2), PSM(3)+E) *
     1               (1.+ PSM(5) * (PSM(2)+N-1.) / (PSM(2)+PSM(3)+N))
     2             + EBETA (N-0.5+PSM(2), PSS(3)+E) * (PSM(4) + 
     3               PSM(6) * (PSM(2)+N-0.5) / (PSM(2)+PSM(3)+N+0.5)) )
         SSI = NSS * ( EBETA (N-1.+PSS(2), PSS(3)+E) *
     1               (1.+ PSS(5) * (PSS(2)+N-1.) / (PSS(2)+PSS(3)+N))
     2             + EBETA (N-0.5+PSS(2), PSS(3)+E) * (PSS(4) + 
     3               PSS(6) * (PSS(2)+N-0.5) / (PSS(2)+PSS(3)+N+0.5)) )
       ELSE
         SMI = ZERO
         SSI = LSI * PSS(1)
       END IF
*
       GLI(KN) = NGL * ( EBETA (N-1.+PGL(2), PGL(3)+E) *
     1               (1.+ PGL(5) * (PGL(2)+N-1.) / (PGL(2)+PGL(3)+N))
     2               + EBETA (N-0.5+PGL(2), PGL(3)+E) * (PGL(4) + 
     3               PGL(6) * (PGL(2)+N-0.5) / (PGL(2)+PGL(3)+N+0.5)) )
* 
* ---------------------------------------------------------------------
*
* ..Arrays of non-singlet and singlet quark combinations for N_f = 3
*
       VAI(KN) = UVI + DVI + SMI
       M3I(KN) = UVI - DVI
       M8I(KN) = VAI(KN) - 3.* SMI
*
       SGI(KN) = UVI + DVI + LSI + SSI
       P3I(KN) = M3I(KN) - 2.* DLI
       P8I(KN) = SGI(KN) - 3.* SSI
*
* ---------------------------------------------------------------------
*
  1    CONTINUE 
       IINNEW = 1
*
       RETURN
       END
*
* =================================================================av==
