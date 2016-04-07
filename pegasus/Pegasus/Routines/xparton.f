*
* ..File xparton.f
*
*
* ..The subroutine  XPARTON  for the Mellin inversion of the parton 
*    densities evolved to the scale  M2 = mu_f^2  by  EVNFFN  (for 
*    IVFNS = 0)  or  EVNVFN.  
*    The results are returned as  x*f_i(x,mu_f^2)  by the array  PDFX 
*    (see the header of evnffn.f for the notation for f_i).
* 
* .. L1 and L2  are the upper and lower limits for the inversion loop. 
*
* ..The integration parameters (for the complex-N contour and for the 
*    Gauss quadratures) are taken from the common-blocks  NCONT,  MOMS 
*    WEIGHTS  and  INVFST.  The initial scale M20 and the corresponding
*    value  AS0  of  a_s = alpha_s/(4 pi) at  R20  are specified in the
*    common-block  ASINP.  The corresponding numbers for/at the heavy-
*    quark thresholds, used for a non-vanishing  IVFNS in VARFLV,  are
*    given in  ASFTHR.  The fixed scale logaritm  LOGFR = ln (M2/R2)  
*    is provided by  FRRAT.
* 
* ..The routine saves couplings and moments. They are re-used for the
*    next inversion in case  M2  and the inputs have not been changed.
*
* =====================================================================
*
*
       SUBROUTINE XPARTON (PDFX, ASOUT, X, M2, L1, L2, IPSTD)

       IMPLICIT DOUBLE PRECISION (A - Z)
       INTEGER NF, NFF, L1, L2, IPSTD, NDIM, IFAST, IVFNS, IMOLD, 
     1         IINNEW, NM, NMSAVE, I1, I2
       PARAMETER (NDIM = 144)
       DIMENSION WN(NDIM), PDFX(-6:6)
       PARAMETER ( IPI = 1.D0 / 3.1415 92653 58979 D0 )
       DOUBLE COMPLEX CC, CCP, XNM 
       DOUBLE COMPLEX CEX(NDIM), NA(NDIM), PDFN(NDIM,-6:6)
*
* ---------------------------------------------------------------------
*
* ..Input common blocks 
* 
       COMMON / NCONT  / C, CC                              
       COMMON / WEIGHTS/ WN
       COMMON / MOMS   / NA  
       COMMON / INVFST / IFAST
       COMMON / INPNEW / IINNEW
       COMMON / NFFIX  / NFF
       COMMON / VARFLV / IVFNS 
       COMMON / FRRAT  / LOGFR
       COMMON / ASINP  / AS0, M20
       COMMON / ASFTHR / ASC, M2C, ASB, M2B, AST, M2T
*
* ..Variables to be saved for the next call
*
       SAVE PDFN, ASI, ASF, NF, M2SAVE, NMSAVE
*
* ---------------------------------------------------------------------
*
* ..Check for new scale mu_f^2 or new input call 
*
       IMOLD = 0
       IF ( (IINNEW .EQ. 0) .AND.  (ABS(M2 - M2SAVE) .LT. 1.D-3) ) THEN
         IMOLD = 1
         GO TO 20
       END IF
*
* ---------------------------------------------------------------------
*
* ..Values of a_s and N_f for the call of EVNFFN or EVNVFN below
*
       R2  = M2 * EXP(-LOGFR)
       IF (IVFNS .EQ. 0) THEN
*
*   Fixed number of flavours
*
         NF  = NFF
         R20 = M20 * R2/M2
         ASI = AS0
         ASF = AS (R2, R20, AS0, NF)
*
       ELSE
*
* ..Variable number of flavours
*
         IF (M2 .GT. M2T) THEN
           NF = 6
           R2T = M2T * R2/M2
           ASI = AST
           ASF = AS (R2, R2T, AST, NF)
*
         ELSE IF (M2 .GT. M2B) THEN
           NF = 5
           R2B = M2B * R2/M2
           ASI = ASB
           ASF = AS (R2, R2B, ASB, NF)
*
         ELSE IF (M2 .GT. M2C) THEN
           NF = 4
           R2C = M2C * R2/M2
           ASI = ASC
           ASF = AS (R2, R2C, ASC, NF)
*
         ELSE
           NF = 3
           R20 = M20 * R2/M2
           ASI = AS0
           ASF = AS (R2, R20, AS0, NF)
*       
         END IF
*
       END IF
*
* ..Output assignment of a_s
*
       ASOUT = ASF
*
* ---------------------------------------------------------------------
*
  20   IF (IFAST . EQ. 0) THEN
*
* ..Integration lengths for the accurate Mellin inversion      
*    zmax =  5  for        x < 0.01, zmax = 14  for  0.01 < x < 0.3, 
*    zmax = 32  for  0.3 < x < 0.7,  zmax = 80  for  0.7  < x  
* 
         IF ( X .LT. 0.01 ) THEN                          
           NM = 64          
         ELSE IF ( X .LT. 0.3 ) THEN                   
           NM = 88 
         ELSE IF ( X .LT. 0.7 ) THEN                   
           NM = 112  
         ELSE                                                          
           NM = 144
         END IF
*
* ---------------------------------------------------------------------
*
       ELSE
*
* ..Integration lengths for the fast Mellin inversion      
*    zmax =  5  for        x < 0.01, zmax = 17  for  0.01 < x < 0.3, 
*    zmax = 36  for  0.3 < x < 0.7,  zmax = 80  for  0.7  < x 
*          
         IF ( X .LT. 0.01 ) THEN    
           NM = 32
         ELSE IF ( X .LT. 0.3 ) THEN
           NM = 48 
         ELSE IF ( X .LT. 0.7 ) THEN
           NM = 64
         ELSE  
           NM = 80
         END IF
*
       END IF 
*
* ---------------------------------------------------------------------
*
* ..Calculation of the moments of the parton densities
*   (taking into account stored results of previous calls)
*
       IF (IVFNS .EQ. 0) THEN
*
* ..Fixed number of flavours
*
         IF (IMOLD .EQ. 0) THEN
           CALL EVNFFN (PDFN, ASI, ASF, NF, 1, NM, IPSTD)
         ELSE
           IF ( NMSAVE .LT. NM) 
     1     CALL EVNFFN (PDFN, ASI, ASF, NF, NMSAVE+1, NM, IPSTD)
         END IF
*
       ELSE
*
* ..Variable number of flavours
*
         IF (IMOLD .EQ. 0) THEN
           CALL EVNVFN (PDFN, ASI, ASF, NF, 1, NM, IPSTD)
         ELSE
           IF ( NMSAVE .LT. NM) 
     1     CALL EVNVFN (PDFN, ASI, ASF, NF, NMSAVE+1, NM, IPSTD)
         END IF
*
       END IF
*
* ---------------------------------------------------------------------
*
* ..Factors for the Mellin inversion
*
       LOGX = LOG (X)
       CCP = CC * IPI
*
       DO 3 I1 = 1, NM
         XNM = - NA(I1) * LOGX
         CEX(I1) = EXP (XNM) * CCP
  3    CONTINUE
*
* ---------------------------------------------------------------------
*
* ..The Mellin inversion of  PDFN (KN, L1 ... L2)
*
       DO 2 I2 = L1, L2
         FUN = 0.D0
       DO 1 I1 = 1, NM
         FZ  = DIMAG (PDFN(I1,I2)*CEX(I1))
         FUN = FUN + WN(I1) * FZ
c        IF (MOD(I1,8) .EQ. 0) WRITE (6,*) I1, FUN * X
*        (view convergence of the integral)
  1    CONTINUE
         PDFX(I2) = FUN * X
  2    CONTINUE
*
* ---------------------------------------------------------------------
*
* ..Saving of variables for the next call
*
       M2SAVE = M2
       NMSAVE = NM
       IINNEW = 0
*
       RETURN                                                           
       END                                                             
*
* =================================================================av==
