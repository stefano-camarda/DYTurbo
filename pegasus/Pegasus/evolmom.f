*
* ..File lh01mom.f
*
* ..A sample program for using the N-space option described in 
*    section 4.6 and 5.5 of the manual.
*
* =====================================================================
*
*
       PROGRAM LH01MOM
*
       IMPLICIT DOUBLE PRECISION (A - Z)
       DOUBLE COMPLEX PDFN(-6:6), N
       PARAMETER ( PI = 3.1415 92653 58979 D0 )
       DIMENSION MS (20)
       INTEGER K1
*
* ..The values for mu^2
*
       DATA MS / 2.0D0,  2.7D0,  3.6D0,  5.D0,   7.D0,   1.D1,
     1           1.4D1,  2.D1,   3.D1,   5.D1,   7.D1,   1.D2,
     2           2.D2,   5.D2,   1.D3,   3.D3,   1.D4,   4.D4,
     3           2.D5,   1.D6 /
*
* ..Access the input parton momentum fractions and normalizations
*
       COMMON / PANORM / AUV, ADV, ALS, ASS, AGL,
     ,                   NUV, NDV, NLS, NSS, NGL
*
* ..The moment N
*
       N = CMPLX (2.D0, 0.D0)
*
       WRITE (6,12) DBLE (N)
  12   FORMAT (/,30X,'N =  ',F5.1,/)
*
* ---------------------------------------------------------------------
*
* ..Loop over mu^2 
*
       DO 1 K1 = 1, 20
         M2 = MS(K1)
*
* ..Call of the Mellin inversion, reconstruction of L+ and L-
*
         CALL NPARTON (PDFN, ASOUT, N, M2, 0, 0, 0, 1)
*
         LMI = (PDFN(-2) - PDFN(-1) - PDFN(2) + PDFN(1)) * 0.5
         LPL =  PDFN(-1) + PDFN(-2) - PDFN(1) - PDFN(2)
*
* ..Output 
*
*
         IF (K1 .EQ. 1) WRITE (6,13)
  13     FORMAT (2X,'M2',7X,'xu_v',5X,'xd_v',5X,'xL_-',5X,'xL_+',
     ,           5X,'xs_+',5X,'xc_+',5X,'xb_+',5X,'xg',/)
         WRITE (6,14) M2, DBLE(PDFN(1)), DBLE(PDFN(2)),
     ,                DBLE(LMI), DBLE(LPL), DBLE(PDFN(-3)),
     ,                DBLE(PDFN(-4)), DBLE(PDFN(-5)), DBLE(PDFN(0))
  14     FORMAT (1PE8.1,8(0PF9.6))
*
   1   CONTINUE
*
       STOP
       END
*
* =================================================================av==
