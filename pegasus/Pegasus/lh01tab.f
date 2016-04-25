*
* ..File lh01tab.f
*
*
* ..The example program of section 4.4 of the manual. Should return the
*    upper part of table 4 (except for the sea quarks at x = 0.9) of 
*    the Les-Houches 2001/2 for the default initialization and input.
*
* =====================================================================
*
*
       PROGRAM LH01TAB
*
       IMPLICIT DOUBLE PRECISION (A - Z)
       DIMENSION PDFX(-6:6), XB(11)
       PARAMETER ( PI = 3.1415 92653 58979 D0 )
       INTEGER K1
*
* ..Access the input parton momentum fractions and normalizations
*
       COMMON / PANORM / AUV, ADV, ALS, ASS, AGL,
     ,                   NUV, NDV, NLS, NSS, NGL
*
* ..The values for mu^2 and x
*
       DATA M2 / 1.D4 /
       DATA XB / 1.D-7,  1.D-6,  1.D-5,  1.D-4,  1.D-3,  1.D-2,
     ,           1.D-1,  3.D-1,  5.D-1,  7.D-1,  9.D-1 /
*
* ---------------------------------------------------------------------
*
* ..General initialization (0 = internal default)
*
       CALL INITEVOL (0)
*
* ..Input initialization (0 = internal default)
*
       CALL INITINP (0)
*
* ..Output of the momentum fractions and normalizations
*
       WRITE(6,10) AUV, ADV, ALS, ASS, AGL
  10   FORMAT (1X,'AUV =',F9.6,2X,'ADV =',F9.6,2X,'ALS =',F9.6,
     ,         2X,'ASS =',F9.6,2X,'AGL =',F9.6)
*
       WRITE(6,11) NUV, NDV, NLS, NSS, NGL
  11   FORMAT (1X,'NUV =',F9.6,2X,'NDV =',F9.6,2X,'NLS =',F9.6,
     ,         2X,'NSS =',F9.6,2X,'NGL =',F9.6,/)
*
* ---------------------------------------------------------------------
*
* ..Loop only over x (M2 is fixed in the Les-Houches tables)
*
       DO 1 K1 = 1, 11
         X = XB (K1)
*
* ..Call of the Mellin inversion, reconstruction of L+ and L-
*
         CALL XPARTON (PDFX, AS, X, M2, -5, 2, 0)
*
         LMI = (PDFX(-2) - PDFX(-1) - PDFX(2) + PDFX(1)) * 0.5
         LPL =  PDFX(-1) + PDFX(-2) - PDFX(1) - PDFX(2)
*
* ..Output to be compared to the upper part of Table 4 
*
         IF (K1 .EQ. 1) WRITE (6,12) AS * 4*PI
  12     FORMAT (2X,'ALPHA_S(MR2(M2)) =  ',F9.6,/)
*
         IF (K1 .EQ. 1) WRITE (6,13)
  13     FORMAT (2X,'x',8X,'xu_v',7X,'xd_v',7X,'xL_-',7X,'xL_+',
     ,           7X,'xs_+',7X,'xc_+',7X,'xb_+',7X,'xg',/)
         WRITE (6,14) X, PDFX(1), PDFX(2), LMI, LPL, PDFX(-3),
     ,                   PDFX(-4), PDFX(-5), PDFX(0)
  14     FORMAT (1PE6.0,1X,8(1PE11.4))
*
   1   CONTINUE
*
       STOP
       END
*
* =================================================================av==
