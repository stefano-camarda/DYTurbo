*
* ..File: initmom.f    
*
*
* ..The counterpart  INITMOM  of both  INITEVOL  and INITPOL  for the 
*    evolution at one complex value  N   of the Mellin variable.  IPOL
*    switches between the unpolarized  (IPAR = 0)  and polarized case
*    (otherwise).  IPAR  sets the mode in which the initialization
*    parameters are read in. 
*
* =====================================================================
*
*
       SUBROUTINE INITMOM (N, IPOL, IPAR)
* 
       IMPLICIT DOUBLE COMPLEX (A - Z)
       INTEGER IPOL, IPAR, NDIM, NMAX, IFAST, IMODEV, IVFNS, NFF, NFMIN,
     1         NFMAX, NFLOW, NFHIGH, NPORD, NUORD, NAORD, NASTPS, K
       PARAMETER ( NDIM = 144, NFMIN = 3,  NFMAX = 6 )
*
       DOUBLE PRECISION PI, EMC, ZETA, CF, CA, TR, LOGFR, FR2
*
* ---------------------------------------------------------------------
*
* ..Input common-blocks 
*
       COMMON / EVMOD  / IMODEV
       COMMON / ORDER  / NPORD
       COMMON / VARFLV / IVFNS
       COMMON / NFFIX  / NFF
       COMMON / FRRAT  / LOGFR
*
* ..Output common-blocks 
*
       COMMON / ASPAR  / NAORD, NASTPS
       COMMON / NFUSED / NFLOW, NFHIGH
       COMMON / ITORD  / NUORD
*
       COMMON / NNUSED / NMAX 
       COMMON / MOMS   / NA(NDIM)
       COMMON / HSUMS  / S(NDIM,6)
*
       COMMON / RZETA  / ZETA(6)
       COMMON / KRON2D / D(2,2)   
       COMMON / COLOUR / CF, CA, TR
*
* ---------------------------------------------------------------------
*
* ..Some constants and the two-dimensional Kronecker symbol
*
       PI  = 3.1415 92653 58979 D0
       EMC = 0.5772 15664 90153 D0
       ZETA(1) = EMC
       ZETA(2) = 1.64493 40668 48226 D0
       ZETA(3) = 1.20205 69031 59594 D0
       ZETA(4) = 1.08232 32337 11138 D0
       ZETA(5) = 1.03692 77551 43370 D0
       ZETA(6) = 1.01734 30619 84449 D0
*
       I = DCMPLX (0.D0, 1.D0)
       D(1,1) = DCMPLX (1.D0, 0.D0)
       D(1,2) = DCMPLX (0.D0, 0.D0)
       D(2,1) = DCMPLX (0.D0, 0.D0)
       D(2,2) = DCMPLX (1.D0, 0.D0)
*
* ..QCD colour factors
*
       CA = 3.D0
       CF = 4./3.D0
       TR = 0.5 D0
*
* ---------------------------------------------------------------------
*
* ..Some default settings of the external initialization parameters
*   (iterated VFNS NLO evolution at mu_f/mu_r = 1).  IFAST  is not used
*
       IFAST  = 0
       IVFNS  = 1
       NFF    = 4
       IMODEV = 1
       NPORD  = 1
       FR2    = 1.D0
*
* ..Override these values by reading the file  usrinit.dat (for IPAR=1)
*    or by calling the subroutine  USRINIT (for IPAR=2) -- or their
*    polarized counterparts.
*
       IF ( IPOL .EQ. 0 ) THEN
*
         IF ( IPAR .EQ. 1 ) THEN
           OPEN (91,FILE='usrinit.dat',STATUS='old')
           READ (91,*) IFAST
           READ (91,*) IVFNS
           READ (91,*) NFF
           READ (91,*) IMODEV
           READ (91,*) NPORD
           READ (91,*) FR2
           CLOSE(91)
         ELSE IF ( IPAR .EQ. 2 ) THEN
           CALL USRINIT (IFAST, IVFNS, NFF, IMODEV, NPORD, FR2)
         END IF
*
       ELSE
*
         IF ( IPAR .EQ. 1) THEN
         OPEN (91,FILE='usrpinit.dat',STATUS='old')
           READ (91,*) IFAST
           READ (91,*) IVFNS
           READ (91,*) NFF
           READ (91,*) IMODEV
           READ (91,*) NPORD
           READ (91,*) FR2
           CLOSE(91)
         ELSE IF ( IPAR .EQ. 2) THEN
           CALL USRPINIT (IFAST, IVFNS, NFF, IMODEV, NPORD, FR2)
         END IF
*
       END IF
       LOGFR = LOG(FR2)
*
* ..Stop some nonsense
*
       IF ( (IVFNS .EQ. 0) .AND. (NFF .LT. 3) ) THEN
         WRITE (6,*) 'Wrong flavour number for FFNS evolution. STOP'
         STOP
       END IF
       IF ( (IVFNS .EQ. 0) .AND. (NFF .GT. 5) ) THEN
         WRITE (6,*) 'Wrong flavour number for FFNS evolution. STOP'
         STOP
       END IF
*
       IF ( NPORD .GT. 2 ) THEN
         WRITE (6,*) 'Specified order in a_s too high. STOP'
         STOP
       END IF
*
       IF ( (IVFNS .NE. 0) .AND. (FR2 .GT. 4.001D0) ) THEN
         WRITE (6,*) 'Too low mu_r for VFNS evolution. STOP'
         STOP
       END IF
*
* ---------------------------------------------------------------------
*
* ..Set NMAX and put the one value of N into the array as NA(1)
*   (the rest is not used here)
* 
         NMAX  = 1
         NA(1) = N
* 
* ..The lowest simple harmonic sums for this value of N
*
       DO 1 K = 1, NMAX
         N1 = NA(K) + 1.
         S(K,1) = PSI(N1) + EMC
         S(K,2) = ZETA(2) - DPSI (N1,1)
         S(K,3) = ZETA(3) + 0.5 * DPSI (N1,2)
         S(K,4) = ZETA(4) - 1./6.D0 * DPSI (N1,3)
         S(K,5) = ZETA(5) + 1./24.D0 * DPSI (N1,4)
         S(K,6) = ZETA(6) - 1./120.D0 * DPSI (N1,5)
  1    CONTINUE
*
* ---------------------------------------------------------------------
*
* ..Some more (internal) evolution and initialization parameters 
*
       NAORD  = NPORD
       NASTPS = 20 
*
       IF ( (IMODEV .EQ. 1) .OR. (IMODEV .EQ. 2) ) THEN
         NUORD = 15
       ELSE
         NUORD = NPORD
       END IF
*
       IF ( IVFNS .EQ. 0) THEN
         NFLOW  = NFF
         NFHIGH = NFF
       ELSE 
         NFLOW  = NFMIN
         NFHIGH = NFMAX
       END IF
*
* ---------------------------------------------------------------------
*
* ..The unpolarized or polarized N-space splitting functions, evolution
*    operators etc.
*
       CALL BETAFCT
       CALL PNS0MOM 
*
       IF ( IPOL .EQ. 0 ) THEN
*
* ..The unpolarized case up to NNLO (NPORD = 2)
*
         CALL PSG0MOM 
         CALL LSGMOM 
*
         IF (NPORD .GT. 0) THEN
           CALL PNS1MOM 
           CALL PSG1MOM 
           CALL USG1MOM 
           IF (NUORD .GT. 1) CALL USG1HMOM
         END IF
*
         IF (NPORD .GT. 1) THEN
           CALL PNS2MOM 
           CALL PSG2MOM 
           CALL UNS2MOM 
           CALL USG2MOM 
           IF (NUORD .GT. 2) CALL USG2HMOM
         END IF
*
         IF (IVFNS .NE. 0) THEN
           CALL ANS2MOM
           CALL ASG2MOM
         END IF
*
       ELSE
*
* ---------------------------------------------------------------------
*
* ..The polarized case, presently restricted to NLO (NPORD =< 1)
*
         CALL PSG0PMOM
         CALL LSGMOM
*
         IF (NPORD .GT. 0) THEN
           CALL PNS1PMOM
           CALL PSG1PMOM
           CALL USG1MOM
           IF (NUORD .GT. 1) CALL USG1HMOM
         END IF
*
         IF (IVFNS .NE. 0) THEN
           CALL ANS2PMOM
           CALL ASG2PMOM
         END IF
       END IF
*
* ---------------------------------------------------------------------
*
       RETURN
       END 
*
* =================================================================av==
