*
* ..File: initpinp.f
*
*
* ..The input initialization for the polarized parton evolution.
*
* ..This routine initializes the parton distributions at the initial
*    scale  M20 (in GeV^2)  and, in the VFNS case, at the flavour
*    thresholds,  Mi2, i = C, B, T, together with the corresponding
*    values of the strong coupling constant. The corresponding routines
*    are  INPPMOMj  and  EVNFTHR  where more information can be found.
*    So far  INPPMOM1,2  are available and selected by  NFORM = 1, 2. 
*
* ..The flavour parameters  IVFNS, NFF  and the (fixed) scale log
*    LOGFR = ln(mu_f^2/mu_r^2)  are provided by the common blocks
*    VARFLV, NFFIX  and  FRRAT.  If called with  IPAR  unequal 1 or 2,
*    the initial conditions for the evolution are those of the 2001
*    Les Houches benchmark.  Other inputs (including the choice of
*    NFORM  are specified by reading the file  usrpinp.f (for IPAR = 1)
*    of by calling the routine  USRPINP  (for IPAR = 2).  
*
* =====================================================================
*
*
       SUBROUTINE INITPINP (IPAR)
* 
       IMPLICIT DOUBLE PRECISION (A - Z)
       INTEGER IPAR, IVFNS, NFF, NF, NFORM, IMOMIN, ISSIMP
       DIMENSION PUV(6), PDV(6), PLM(6), PLP(6), PSM(6), PSP(6), PGL(6)
       PARAMETER ( PI = 3.1415 92653 58979 D0 )
*
* ---------------------------------------------------------------------
*
* ..Input common blocks
*
       COMMON / VARFLV / IVFNS
       COMMON / NFFIX  / NFF
       COMMON / FRRAT  / LOGFR
*
* ..Output common blocks
*
       COMMON / ASINP  / AS0, M20
       COMMON / ASFTHR / ASC, M2C, ASB, M2B, AST, M2T
*
* ---------------------------------------------------------------------
*
* ..Some default settings of scales, couplings and input distributions
*
* ...Initial scale  M20 = M_0^2 (in GeV^2)  and  ASI = alpha_s(M_0^2)
*
       M20 = 2.D0
       ASI = 0.35000 
*
* ...The heavy quark masses squared
*
       MC2 = 2.00001D0
       MB2 = 20.25D0
       MT2 = 3.0625D4
*
* ...Flags for the initial parton distributions
*
       NFORM  = 1
       IMOMIN = 0
       ISSIMP = 0
*
* ---------------------------------------------------------------------
*
* ...Parameters for the Les-Houches-2001 PDFs for use with  inplmom1.f
*
* ...u_v   (PUV(1) is the number of up quarks)
*
       PUV(1) =  1.3D0
       PUV(2) =  0.7D0
       PUV(3) =  3.0D0
       PUV(4) =  0.5D0
       PUV(5) =  0.0D0
       PUV(6) =  3.0D0
*
* ...d_v   (PDV(1) is the number of down quarks)
*
       PDV(1) = -0.5D0
       PDV(2) =  0.7D0
       PDV(3) =  4.0D0
       PDV(4) =  0.5D0
       PDV(5) =  0.0D0
       PDV(6) =  4.0D0
*    _   _
* ...d - u
*
       PLM(1) =  0.0D0
       PLM(2) =  0.7D0
       PLM(3) =  7.0D0
       PLM(4) =  0.5D0
       PLM(5) =  0.0D0
       PLM(6) =  0.0D0
*      _   _
* ...2(d + u)   (for IMOMIN = 1,  PLS(1)  is the momentum fraction)
*
       PLP(1) = -0.2D0
       PLP(2) =  0.3D0
       PLP(3) =  7.0D0
       PLP(4) =  0.5D0
       PLP(5) =  0.0D0
       PLP(6) =  0.0D0
*
* ---------------------------------------------------------------------
*         _
* ...(s - s)   (First moment vanishes. PSM(5) is not used)
*
       PSM(1) =  0.0D0
       PSM(2) =  0.3D0
       PSM(3) =  7.0D0
       PSM(4) =  0.5D0
       PSM(5) =  0.0D0
       PSM(6) =  0.0D0
*         _
* ...(s + s)   (for  ISSIMP = 0, the strange sea is  SP = PSS(1) * LP)
*
       PSP(1) = -0.05D0
       PSP(2) =  0.3D0
       PSP(3) =  7.0D0
       PSP(4) =  0.5D0
       PSP(5) =  0.0D0
       PSP(6) =  0.0D0
*
* ...g   (PGL(1) is the total momentum sum of the partons)
*
       PGL(1) =  1.5D0
       PGL(2) =  0.5D0
       PGL(3) =  5.0D0
       PGL(4) =  0.5D0
       PGL(5) =  0.0D0
       PGL(6) =  0.0D0
*
* ---------------------------------------------------------------------
*
* ..Override these values by reading the file  usrinp.dat (for IPAR=1)
*    or by calling the subroutine  USRINP (for IPAR=2)
*
       IF ( IPAR .EQ. 1) THEN
         OPEN (92,FILE='usrpinp.dat',STATUS='old')
         READ (92,*) M20 
         READ (92,*) ASI
         READ (92,*) MC2
         READ (92,*) MB2
         READ (92,*) MT2
         READ (92,*) NFORM
         READ (92,*) IMOMIN
         READ (92,*) ISSIMP
         READ (92,*) PUV(1), PUV(2), PUV(3), PUV(4), PUV(5), PUV(6)
         READ (92,*) PDV(1), PDV(2), PDV(3), PDV(4), PDV(5), PDV(6)
         READ (92,*) PLM(1), PLM(2), PLM(3), PLM(4), PLM(5), PLM(6)
         READ (92,*) PLP(1), PLP(2), PLP(3), PLP(4), PLP(5), PLP(6)
         READ (92,*) PSM(1), PSM(2), PSM(3), PSM(4), PSM(5), PSM(6)
         READ (92,*) PSP(1), PSP(2), PSP(3), PSP(4), PSP(5), PSP(6)
         READ (92,*) PGL(1), PGL(2), PGL(3), PGL(4), PGL(5), PGL(6)
         CLOSE(92)
       ELSE IF ( IPAR .EQ. 2) THEN
         CALL USRPINP (PUV, PDV, PLM, PLP, PSM, PSP, PGL, M20, ASI,
     1                 MC2, MB2, MT2, NFORM, IMOMIN, ISSIMP)
       END IF
*
* ---------------------------------------------------------------------
*
* ..Stop some nonsense
*
       IF ( (IVFNS .EQ. 1) .AND. (M20 .GT. MC2) ) THEN
         WRITE (6,*) 'Too high mu_0 for VFNS evolution. STOP'
         STOP
       END IF
*
       IF ( (ASI .GT. 2.D0) .OR. (ASI .LT. 2.D-2) ) THEN
         WRITE (6,*) 'alpha_s out of range. STOP'
         STOP
       END IF
*
       IF ( (IVFNS .EQ. 1) .AND. (MC2 .GT. MB2) ) THEN
         WRITE (6,*) 'Wrong charm-bottom mass hierarchy. STOP'
         STOP
       END IF
       IF ( (IVFNS .EQ. 1) .AND. (MB2 .GT. MT2) ) THEN
         WRITE (6,*) 'Wrong bottom-top mass hierarchy. STOP'
         STOP
       END IF
*
       IF ( (NFORM .NE. 1) .AND. (NFORM .NE. 2) ) THEN
         WRITE (6,*) 'Inappropriate value of NFORM. STOP'
         STOP
       END IF
*
       IF ( (PUV(2) .LT. 1.D-3) .OR. (PDV(2) .LT. 1.D-3) ) THEN
         WRITE (6,*) 'Wrong small-x power of u_v or d_v. STOP'
         STOP
       END IF
       IF ( (PUV(3) .LT. 1.D-1) .OR. (PDV(3) .LT. 1.D-1) ) THEN
         WRITE (6,*) 'Wrong large-x power of u_v or d_v. STOP'
         STOP
       END IF
*
       IF (  (PLM(2) .LT. 1.D-3) .OR. (PLP(2) .LT. 1.D-3)
     ,  .OR. (PSP(2) .LT. 1.D-3) .OR. (PGL(2) .LT. 1.D-3) ) THEN
         WRITE (6,*) 'Wrong small-x power of sea or gluon. STOP'
         STOP
       END IF
       IF (  (PLM(3) .LT. 1.D-1) .OR. (PLP(3) .LT. 1.D-1)
     ,  .OR. (PSP(3) .LT. 1.D-1) .OR. (PGL(3) .LT. 1.D-1) ) THEN
         WRITE (6,*) 'Wrong large-x power of sea or gluon. STOP'
         STOP
       END IF
*
*
* ---------------------------------------------------------------------
*
* ...For mu_r unequal mu_f  AS0  is different from the input parameter 
*    ASM = a_s(M_0^2). The VFNS evolution starts with n_f = 3 at M_0^2.
*
       ASM = ASI / (4.* PI)
       R20 = M20 * EXP(-LOGFR)
       NF = NFF
       IF (IVFNS .NE. 0)  NF  = 3
       AS0 = AS (R20, M20, ASM, NF)
*
* ..Input initialization (including threshold values for the VFNS case)
*
       IF ( NFORM .EQ. 1 ) THEN
         CALL INPPMOM1 (PUV,PDV, PLM,PLP, PSM,PSP, PGL, IMOMIN,ISSIMP)
       ELSE IF ( NFORM .EQ. 2 ) THEN
         CALL INPPMOM2 (PUV,PDV, PLM,PLP, PSM,PSP, PGL, IMOMIN,ISSIMP)
       END IF
*
       IF (IVFNS .NE. 0) THEN
         CALL EVNFTHR (MC2, MB2, MT2)
       END IF
*
* ---------------------------------------------------------------------
*
       RETURN
       END
*
* =================================================================av==
