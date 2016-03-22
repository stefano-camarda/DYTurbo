*
* ..File: inppmom2.f   (polarized) 
*
*
* ..The subroutine  INPPMOM2  (version 2) for the initial polarized
*    light-quark distributions of a hadron in N-space.  The moments are 
*    stored in the common-block  PAINP for an external NDIM-dimensional
*    array  NA  of complex moments specified in the common-block  MOMS.
*
* ..The functional form of the distributions for this version reads
*
*        xf = Nf * Pf(1) * x^Pf(2) * (1 - x)^Pf(3)  
*             * ( 1 + x^0.5 * Pf(4) + x * Pf(5) + x^1.5 * Pf(6))    (1)
*
*    for  f = UV = u - ubar,  DV = d - dbar,  DL = dbar - ubar, 
*    LS = 2 * (dbar + ubar),  SM = s - sbar,  SS = s + sbar,   GL = g. 
*
*    For  IMOMIN = 0,  all  Nf = 1  in (1). Otherwise the  Nf  are 
*    chosen such that  Pf(1)  represents the  respective first moment. 
*
* ..The normalization factors  Nf  and the first moments  Af  are 
*    written to the common-block PANORM.  The flag  IINNEW  in  INPNEW
*    is set to '1' at the end of this routine.
*
* =====================================================================
*
*
       SUBROUTINE INPPMOM2 (PUV,PDV,PDL,PLS,PSM,PSS,PGL,IMOMIN,ISSIMP)
*
       IMPLICIT DOUBLE COMPLEX (A - Z)            
       INTEGER NMAX, NDIM, IMOMIN, ISSIMP, IINNEW, KN
       PARAMETER (NDIM = 144)
       PARAMETER ( E = (1.D0, 0.D0), ZERO = (0.D0, 0.D0) )
       DOUBLE PRECISION PUV(6), PDV(6), PDL(6), PLS(6), PSM(6), PSS(6),
     1                  PGL(6), AUV, ADV, ALS, ASM, ASS, AGL, NUV, NDV, 
     2                  NLS, NSM, NSS, NGL
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
* ..The normalization factors and first moments ...
*
       IF (IMOMIN .EQ. 0) THEN
*
* ...if the former are the input parameters ...
*
         NUV = PUV(1)
         NDV = PDV(1)
         NDL = PDL(1)
         NLS = PLS(1)
         NGL = PGL(1)
*
         AUV = NUV * ( EBETA (PUV(2)*E, PUV(3)+E) *
     1               (1.+ PUV(5) * PUV(2) / (PUV(2)+ PUV(3)+ 1.))
     2               + EBETA (0.5+PUV(2), PUV(3)+E) * (PUV(4) +
     3               PUV(6) * (PUV(2)+ 0.5) / (PUV(2)+ PUV(3)+ 1.5)) )
         ADV = NDV * ( EBETA (PDV(2)*E, PDV(3)+E) *
     1               (1.+ PDV(5) * PDV(2) / (PDV(2)+ PDV(3)+ 1.))
     2               + EBETA (0.5+PDV(2), PDV(3)+E) * (PDV(4) +
     3               PDV(6) * (PDV(2)+ 0.5) / (PDV(2)+ PDV(3)+ 1.5)) )
         ADL = NDL * ( EBETA (PDL(2)*E, PDL(3)+E) *
     1               (1.+ PDL(5) * PDL(2) / (PDL(2)+ PDL(3)+ 1.))
     2               + EBETA (0.5+PDL(2), PDL(3)+E) * (PDL(4) +
     3               PDL(6) * (PDL(2)+ 0.5) / (PDL(2)+ PDL(3)+ 1.5)) )
         ALS = NLS * ( EBETA (PLS(2)*E, PLS(3)+E) *
     1               (1.+ PLS(5) * PLS(2) / (PLS(2)+ PLS(3)+ 1.))
     2               + EBETA (0.5+PLS(2), PLS(3)+E) * (PLS(4) +
     3               PLS(6) * (PLS(2)+ 0.5) / (PLS(2)+ PLS(3)+ 1.5)) )
         AGL = NGL * ( EBETA (PGL(2)*E, PGL(3)+E) *
     1               (1.+ PGL(5) * PGL(2) / (PGL(2)+ PGL(3)+ 1.))
     2               + EBETA (0.5+PGL(2), PGL(3)+E) * (PGL(4) +
     3               PGL(6) * (PGL(2)+ 0.5) / (PGL(2)+ PGL(3)+ 1.5)) )
*
         IF (ISSIMP .EQ. 0) THEN
*
           NSM = PSM(1)
           NSS = PSS(1)
           ASM = NSM * ( EBETA (PSM(2)*E, PSM(3)+E) *
     1               (1.+ PSM(5) * PSM(2) / (PSM(2)+ PSM(3)+ 1.))
     2               + EBETA (0.5+PSM(2), PSM(3)+E) * (PSM(4) +
     3               PSM(6) * (PSM(2)+ 0.5) / (PSM(2)+ PSM(3)+ 1.5)) )
           ASS = NSS * ( EBETA (PSS(2)*E, PSS(3)+E) *
     1               (1.+ PSS(5) * PSS(2) / (PSS(2)+ PSS(3)+ 1.))
     2               + EBETA (0.5+PSS(2), PSS(3)+E) * (PSS(4) +
     3               PSS(6) * (PSS(2)+ 0.5) / (PSS(2)+ PSS(3)+ 1.5)) )
         ELSE
*
           NSM = 0.D0
           ASM = 0.D0
           NSS = PSS(1) * NLS
           ASS = PSS(1) * ALS
*
         END IF
*
* ---------------------------------------------------------------------
*
       ELSE
*
* ...and if the latter are the input parameters  
*
         AUV = PUV(1)
         ADV = PDV(1)
         ADL = PDL(1)
         ALS = PLS(1)
         AGL = PGL(1)
*
         NUV = AUV / ( EBETA (PUV(2)*E, PUV(3)+E) *
     1               (1.+ PUV(5) * PUV(2) / (PUV(2)+ PUV(3)+ 1.))
     2               + EBETA (0.5+PUV(2), PUV(3)+E) * (PUV(4) +
     3               PUV(6) * (PUV(2)+ 0.5) / (PUV(2)+ PUV(3)+ 1.5)) )
         NDV = ADV / ( EBETA (PDV(2)*E, PDV(3)+E) *
     1               (1.+ PDV(5) * PDV(2) / (PDV(2)+ PDV(3)+ 1.))
     2               + EBETA (0.5+PDV(2), PDV(3)+E) * (PDV(4) +
     3               PDV(6) * (PDV(2)+ 0.5) / (PDV(2)+ PDV(3)+ 1.5)) )
         NDL = ADL / ( EBETA (PDL(2)*E, PDL(3)+E) *
     1               (1.+ PDL(5) * PDL(2) / (PDL(2)+ PDL(3)+ 1.))
     2               + EBETA (0.5+PDL(2), PDL(3)+E) * (PDL(4) +
     3               PDL(6) * (PDL(2)+ 0.5) / (PDL(2)+ PDL(3)+ 1.5)) )
         NLS = ALS / ( EBETA (PLS(2)*E, PLS(3)+E) *
     1               (1.+ PLS(5) * PLS(2) / (PLS(2)+ PLS(3)+ 1.))
     2               + EBETA (0.5+PLS(2), PLS(3)+E) * (PLS(4) +
     3               PLS(6) * (PLS(2)+ 0.5) / (PLS(2)+ PLS(3)+ 1.5)) )
         NGL = AGL / ( EBETA (PGL(2)*E, PGL(3)+E) *
     1               (1.+ PGL(5) * PGL(2) / (PGL(2)+ PGL(3)+ 1.))
     2               + EBETA (0.5+PGL(2), PGL(3)+E) * (PGL(4) +
     3               PGL(6) * (PGL(2)+ 0.5) / (PGL(2)+ PGL(3)+ 1.5)) )
*
         IF (ISSIMP .EQ. 0) THEN
*
           ASM = PSM(1)
           ASS = PSS(1)
           NSM = ASM / ( EBETA (PSM(2)*E, PSM(3)+E) *
     1               (1.+ PSM(5) * PSM(2) / (PSM(2)+ PSM(3)+ 1.))
     2               + EBETA (0.5+PSM(2), PSM(3)+E) * (PSM(4) +
     3               PSM(6) * (PSM(2)+ 0.5) / (PSM(2)+ PSM(3)+ 1.5)) )
           NSS = ASS / ( EBETA (PSS(2)*E, PSS(3)+E) *
     1               (1.+ PSS(5) * PSS(2) / (PSS(2)+ PSS(3)+ 1.))
     2               + EBETA (0.5+PSS(2), PSS(3)+E) * (PSS(4) +
     3               PSS(6) * (PSS(2)+ 0.5) / (PSS(2)+ PSS(3)+ 1.5)) )
         ELSE
*
           NSM = 0.D0
           ASM = 0.D0
           NSS = PSS(1) * NLS
           ASS = PSS(1) * ALS
*
         END IF
       END IF
*
* --------------------------------------------------------------------- 
*
* ..Begin of the Mellin-N loop 
*
       DO 1 KN = 1, NMAX
       N = NA(KN)
*
* ..Up and down quark distributions 
*
       UVI = NUV * ( EBETA (N-1.+PUV(2), PUV(3)+E) *
     1               (1.+ PUV(5) * (PUV(2)+N-1.) / (PUV(2)+PUV(3)+N))
     2             + EBETA (N-0.5+PUV(2), PUV(3)+E) * (PUV(4) +
     3               PUV(6) * (PUV(2)+N-0.5) / (PUV(2)+PUV(3)+N+0.5)) )
       DVI = NDV * ( EBETA (N-1.+PDV(2), PDV(3)+E) *
     1               (1.+ PDV(5) * (PDV(2)+N-1.) / (PDV(2)+PDV(3)+N))
     2             + EBETA (N-0.5+PDV(2), PDV(3)+E) * (PDV(4) +
     3               PDV(6) * (PDV(2)+N-0.5) / (PDV(2)+PDV(3)+N+0.5)) )
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
* ---------------------------------------------------------------------
*
* ..Strange quark and gluon distributions  
*
       IF (ISSIMP .EQ. 0) THEN
         SMI = NSM * ( EBETA (N-1.+PSM(2), PSM(3)+E) *
     1               (1.+ PSM(5) * (PSM(2)+N-1.) / (PSM(2)+PSM(3)+N))
     2             + EBETA (N-0.5+PSM(2), PSM(3)+E) * (PSM(4) +
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
* ..Arrays of non-singlet and singlet quark combinations for n_f = 3
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
