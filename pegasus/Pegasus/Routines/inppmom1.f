*
* ..File: inppmom1.f   (polarized)
*
*
* ..The subroutine  INPPMOM1  (version 1) for the initial polarized
*    light-quark distributions of a hadron in N-space.  The moments are 
*    stored in the common-block  PAINP for an external NDIM-dimensional
*    array  NA  of complex moments specified in the common-block  MOMS.
*
* ..The functional form of the distributions for this version reads
*
*        xf = Nf * Pf(1) * x^Pf(2) * (1 - x)^Pf(3)  
*             * ( 1 + x^Pf(4) * Pf(5) + x * Pf(6) )                 (1)
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
       SUBROUTINE INPPMOM1 (PUV,PDV,PDL,PLS,PSM,PSS,PGL,IMOMIN,ISSIMP)
*
       IMPLICIT DOUBLE COMPLEX (A - Z)            
       INTEGER NMAX, NDIM, IMOMIN, ISSIMP, IINNEW, KN
       PARAMETER (NDIM = 144)
       PARAMETER ( E = (1.D0, 0.D0), ZERO = (0.D0, 0.D0) )
       DOUBLE PRECISION PUV(6), PDV(6), PDL(6), PLS(6), PSM(6), PSS(6),
     1                  PGL(6), AUV, ADV, ADL, ALS, ASM, ASS, AGL, NUV,
     2                  NDV, NDL, NLS, NSM, NSS, NGL
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
         AUV = NUV * ( EBETA (PUV(2)*E, PUV(3)+E)
     1               * (1.+ PUV(6) * PUV(2) / (PUV(2) + PUV(3) + 1.))
     2               + PUV(5) * EBETA (E*(PUV(2)+PUV(4)), PUV(3)+E) )
         ADV = NDV * ( EBETA (PDV(2)*E, PDV(3)+E)
     1               * (1.+ PDV(6) * PDV(2) / (PDV(2) + PDV(3) + 1.))
     2               + PDV(5) * EBETA (E*(PDV(2)+PDV(4)), PDV(3)+E) )
         ADL = NDL * ( EBETA (PDL(2)*E, PDL(3)+E)
     1               * (1.+ PDL(6) * PDL(2) / (PDL(2) + PDL(3) + 1.))
     2               + PDL(5) * EBETA (E*(PDL(2)+PDL(4)), PDL(3)+E) )
         ALS = NLS * ( EBETA (PLS(2)*E, PLS(3)+E)
     1               * (1.+ PLS(6) * PLS(2) / (PLS(2) + PLS(3) + 1.))
     2               + PLS(5) * EBETA (E*(PLS(2)+PLS(4)), PLS(3)+E) )
         AGL = NGL * ( EBETA (PGL(2)*E, PGL(3)+E)
     1               * (1.+ PGL(6) * PGL(2) / (PGL(2) + PGL(3) + 1.))
     2               + PGL(5) * EBETA (E*(PGL(2)+PGL(4)), PGL(3)+E) )
*
         IF (ISSIMP .EQ. 0) THEN
*
           NSM = PSM(1)
           NSS = PSS(1)
           ASM = NSM * ( EBETA (PSM(2)*E, PSM(3)+E)
     1               * (1.+ PSM(6) * PSM(2) / (PSM(2) + PSM(3) + 1.))
     2               + PSM(5) * EBETA (E*(PSM(2)+PSM(4)), PSM(3)+E) )
           ASS = NSS * ( EBETA (PSS(2)*E, PSS(3)+E)
     1               * (1.+ PSS(6) * PSS(2) / (PSS(2) + PSS(3) + 1.))
     2               + PSS(5) * EBETA (E*(PSS(2)+PSS(4)), PSS(3)+E) )
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
         NUV = AUV / ( EBETA (PUV(2)*E, PUV(3)+E)
     1               * (1.+ PUV(6) * PUV(2) / (PUV(2) + PUV(3) + 1.))
     2               + PUV(5) * EBETA (E*(PUV(2)+PUV(4)), PUV(3)+E) )
         NDV = ADV / ( EBETA (PDV(2)*E, PDV(3)+E)
     1               * (1.+ PDV(6) * PDV(2) / (PDV(2) + PDV(3) + 1.))
     2               + PDV(5) * EBETA (E*(PDV(2)+PDV(4)), PDV(3)+E) )
         NDL = ADL / ( EBETA (PDL(2)*E, PDL(3)+E)
     1               * (1.+ PDL(6) * PDL(2) / (PDL(2) + PDL(3) + 1.))
     2               + PDL(5) * EBETA (E*(PDL(2)+PDL(4)), PDL(3)+E) )
         NLS = ALS / ( EBETA (PLS(2)*E, PLS(3)+E)
     1               * (1.+ PLS(6) * PLS(2) / (PLS(2) + PLS(3) + 1.))
     2               + PLS(5) * EBETA (E*(PLS(2)+PLS(4)), PLS(3)+E) )
         NGL = AGL / ( EBETA (PGL(2)*E, PGL(3)+E)
     1               * (1.+ PGL(6) * PGL(2) / (PGL(2) + PGL(3) + 1.))
     2               + PGL(5) * EBETA (E*(PGL(2)+PGL(4)), PGL(3)+E) )
*
         IF (ISSIMP .EQ. 0) THEN
*
           ASM = PSM(1)
           ASS = PSS(1)
           NSM = ASM / ( EBETA (PSM(2)*E, PSM(3)+E)
     1               * (1.+ PSM(6) * PSM(2) / (PSM(2) + PSM(3) + 1.))
     2               + PSM(5) * EBETA (E*(PSM(2)+PSM(4)), PSM(3)+E) )
           NSS = ASS / ( EBETA (PSS(2)*E, PSS(3)+E)
     1               * (1.+ PSS(6) * PSS(2) / (PSS(2) + PSS(3) + 1.))
     2               + PSS(5) * EBETA (E*(PSS(2)+PSS(4)), PSS(3)+E) )
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

       UVI = NUV * ( EBETA (N-1.+PUV(2), PUV(3)+E)
     1           * (1.+ PUV(6) * (PUV(2)+N-1.) / (PUV(2) + PUV(3) + N))
     2           + PUV(5) * EBETA (N-1.+PUV(2)+PUV(4), PUV(3)+E) )
       DVI = NDV * ( EBETA (N-1.+PDV(2), PDV(3)+E)
     1           * (1.+ PDV(6) * (PDV(2)+N-1.) / (PDV(2) + PDV(3) + N))
     2           + PDV(5) * EBETA (N-1.+PDV(2)+PDV(4), PDV(3)+E) )
*
       DLI = PDL(1) * ( EBETA (N-1.+PDL(2), PDL(3)+E)
     1           * (1.+ PDL(6) * (PDL(2)+N-1.) / (PDL(2) + PDL(3) + N))
     2           + PDL(5) * EBETA (N-1.+PDL(2)+PDL(4), PDL(3)+E) )
       LSI = NLS * ( EBETA (N-1.+PLS(2), PLS(3)+E)
     1           * (1.+ PLS(6) * (PLS(2)+N-1.) / (PLS(2) + PLS(3) + N))
     2           + PLS(5) * EBETA (N-1.+PLS(2)+PLS(4), PLS(3)+E) )
* 
* ---------------------------------------------------------------------
*
* ..Strange quark and gluon distributions  
*
       IF (ISSIMP .EQ. 0) THEN
         SMI = NSM * ( EBETA (N-1.+PSM(2), PSM(3)+E)
     1           * (1.+ PSM(6) * (PSM(2)+N-1.) / (PSM(2) + PSM(3) + N))
     2           + PSM(5) * EBETA (N-1.+PSM(2)+PSM(4), PSM(3)+E) )
         SSI = NSS * ( EBETA (N-1.+PSS(2), PSS(3)+E)
     1           * (1.+ PSS(6) * (PSS(2)+N-1.) / (PSS(2) + PSS(3) + N))
     2           + PSS(5) * EBETA (N-1.+PSS(2)+PSS(4), PSS(3)+E) )
       ELSE
         SMI = ZERO
         SSI = LSI * PSS(1)
       END IF
*
       GLI(KN) = NGL * ( EBETA (N-1.+PGL(2), PGL(3)+E)
     1          * (1.+ PGL(6) * (PGL(2)+N-1.) / (PGL(2) + PGL(3) + N))
     2          + PGL(5) * EBETA (N-1.+PGL(2)+PGL(4), PGL(3)+E) )
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
