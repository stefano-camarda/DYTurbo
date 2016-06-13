C-----------------------------------------------------------------------
C     Martin-Roberts-Stirling Parton Distributions
C
      SUBROUTINE  mrs (nset,x,qsq,pden)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION pden(-6:6)
      scale = dsqrt(qsq)
      mode = nset
      IF (mode.LE.3) THEN
          CALL mrs123 (x,scale,mode,upv,dnv,sea,str,chm,bot,glu)
      ELSEIF (mode.EQ.4.OR.mode.EQ.5) THEN
          CALL mrseb  (x,scale,mode,upv,dnv,sea,str,chm,bot,glu)
      ELSE
          PRINT *,'Bad MRS mode ',nset
          STOP
      ENDIF
      pden(0) = glu/x
      pden(1) = (upv+sea)/x
      pden(2) = (dnv+sea)/x
      pden(-1) = sea/x
      pden(-2) = pden(-1)
      pden(3) = str/x
      pden(-3) = pden(3)
      pden(4) = chm/x
      pden(-4) = pden(4)
      pden(5) = bot/x
      pden(-5) = pden(5)
      pden(6) = 0d0
      pden(-6) = pden(6)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE MRS123(X,SCALE,MODE,UPV,DNV,SEA,STR,CHM,BOT,GLU)
C***************************************************************C
C                                                               C
C  MODE 1 CORRESPONDS TO                                        C
C  MARTIN, ROBERTS, STIRLING (SOFT GLUE)  WITH LAMBDA=.107 GEV  C
C                                                               C
C  MODE 2 CORRESPONDS TO                                        C
C  MARTIN, ROBERTS, STIRLING (HARD GLUE)  WITH LAMBDA=.250 GEV  C
C                                                               C
C  MODE 3 CORRESPONDS TO                                        C
C  MARTIN, ROBERTS, STIRLING (1/RTX GLUE) WITH LAMBDA=.178 GEV  C
C                                                               C
C                         -*-                                   C
C                                                               C
C    (NOTE THAT X TIMES THE PARTON DISTRIBUTION FUNCTION        C
C    IS RETURNED I.E. G(X) = GLU/X ETC, AND THAT "SEA"          C
C    IS THE LIGHT QUARK SEA I.E. UBAR(X)=DBAR(X)=  ...          C
C    = SEA/X FOR A PROTON.  IF IN DOUBT, CHECK THE              C
C    MOMENTUM SUM RULE! NOTE ALSO THAT SCALE=Q IN GEV)          C
C                                                               C
C                         -*-                                   C
C                                                               C
C     (THE RANGE OF APPLICABILITY IS CURRENTLY:                 C
C     10**-4 < X < 1  AND  5 < Q**2 < 1.31 * 10**6              C
C     HIGHER Q**2 VALUES CAN BE SUPPLIED ON REQUEST             C
C     - PROBLEMS, COMMENTS ETC TO SRG$T3@GEN                    C
C                                                               C
C                                                               C
C***************************************************************C
      IMPLICIT REAL*8(A-H,O-Z)
      IF(MODE.EQ.1) CALL STRUC1(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GLU)
      IF(MODE.EQ.2) CALL STRUC2(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GLU)
      IF(MODE.EQ.3) CALL STRUC3(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GLU)
      RETURN
      END
C
      SUBROUTINE STRUC1(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GLU)
C :::::::::::: MARTIN ROBERTS STIRLING :::::SOFT GLUE::::107 MEV::::
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(6,42,19),G(6),XX(42),N0(6)
      SAVE
      DATA XX/1.D-4,2.D-4,4.D-4,6.D-4,8.D-4,
     .       1.D-3,2.D-3,4.D-3,6.D-3,8.D-3,
     .       1.D-2,2.D-2,4.D-2,6.D-2,8.D-2,
     .     .1D0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
     .     .5D0,.525D0,.55D0,.575D0,.6D0,.65D0,.7D0,.75D0,
     .     .8D0,.9D0,1.D0/
      DATA XMIN,XMAX,QSQMIN,QSQMAX/1.D-4,1.D0,5.D0,1310720.D0/
      DATA N0/2,5,4,5,0,0/
      DATA INIT/0/
      IF(INIT.NE.0) GOTO 10
      OPEN(UNIT=31,FILE='mrs1.dat'
     &     ,STATUS='OLD') !,READONLY)
      READ(31,*)
      INIT=1
      DO 20 N=1,41
      DO 20 M=1,19
      READ(31,50)F(1,N,M),F(2,N,M),F(3,N,M),F(4,N,M),F(5,N,M),F(6,N,M)
         DO 25 I=1,6
  25     F(I,N,M)=F(I,N,M)/(1.D0-XX(N))**N0(I)
  20  CONTINUE
      DO 30 J=1,15
      XX(J)=DLOG10(XX(J))+1.1D0
      DO 30 I=1,5
      DO 30 K=1,19
  30  F(I,J,K)=DLOG(F(I,J,K))*F(I,16,K)/DLOG(F(I,16,K))
  50  FORMAT(6F10.5)
      DO 40 I=1,6
      DO 40 M=1,19
  40  F(I,42,M)=0.D0
      CLOSE(31)
  10  CONTINUE
      IF(X.LT.XMIN) X=XMIN
      IF(X.GT.XMAX) X=XMAX
      QSQ=SCALE**2
      IF(QSQ.LT.QSQMIN) QSQ=QSQMIN
      IF(QSQ.GT.QSQMAX) QSQ=QSQMAX
      XXX=X
      IF(X.LT.1.D-1) XXX=DLOG10(X)+1.1D0
      N=0
  70  N=N+1
      IF(XXX.GT.XX(N+1)) GOTO 70
      A=(XXX-XX(N))/(XX(N+1)-XX(N))
      RM=DLOG(QSQ/QSQMIN)/DLOG(2.D0)
      B=RM-DINT(RM)
      M=1+IDINT(RM)
      DO 60 I=1,6
      G(I)= (1.D0-A)*(1.D0-B)*F(I,N,M)+(1.D0-A)*B*F(I,N,M+1)
     .    + A*(1.D0-B)*F(I,N+1,M)  + A*B*F(I,N+1,M+1)
      IF(N.GE.16) GOTO 65
      IF(I.EQ.6) GOTO 65
          FAC=(1.D0-B)*F(I,16,M)+B*F(I,16,M+1)
 
          G(I)=FAC**(G(I)/FAC)
  65  CONTINUE
      G(I)=G(I)*(1.D0-X)**N0(I)
  60  CONTINUE
      UPV=G(1)
      DNV=G(2)
      SEA=G(4)
      STR=G(4)
      CHM=G(5)
      GLU=G(3)
      BOT=G(6)
      RETURN
      END
      SUBROUTINE STRUC2(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GLU)
C :::::::::::: MARTIN ROBERTS STIRLING :::::HARD GLUE::::250 MEV::::
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(6,42,19),G(6),XX(42),N0(6)
      SAVE
      DATA XX/1.D-4,2.D-4,4.D-4,6.D-4,8.D-4,
     .       1.D-3,2.D-3,4.D-3,6.D-3,8.D-3,
     .       1.D-2,2.D-2,4.D-2,6.D-2,8.D-2,
     .     .1D0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
     .     .5D0,.525D0,.55D0,.575D0,.6D0,.65D0,.7D0,.75D0,
     .     .8D0,.9D0,1.D0/
      DATA XMIN,XMAX,QSQMIN,QSQMAX/1.D-4,1.D0,5.D0,1310720.D0/
      DATA N0/2,5,4,5,0,0/
      DATA INIT/0/
      IF(INIT.NE.0) GOTO 10
      OPEN(UNIT=32,FILE='mrs2.dat'
     &     ,STATUS='OLD') !,READONLY)
      READ(32,*)
      INIT=1
      DO 20 N=1,41
      DO 20 M=1,19
      READ(32,50)F(1,N,M),F(2,N,M),F(3,N,M),F(4,N,M),F(5,N,M),F(6,N,M)
         DO 25 I=1,6
  25     F(I,N,M)=F(I,N,M)/(1.D0-XX(N))**N0(I)
  20  CONTINUE
      DO 30 J=1,15
      XX(J)=DLOG10(XX(J))+1.1D0
      DO 30 I=1,5
      DO 30 K=1,19
  30  F(I,J,K)=DLOG(F(I,J,K))*F(I,16,K)/DLOG(F(I,16,K))
  50  FORMAT(6F10.5)
      DO 40 I=1,6
      DO 40 M=1,19
  40  F(I,42,M)=0.D0
      CLOSE(32)
  10  CONTINUE
      IF(X.LT.XMIN) X=XMIN
      IF(X.GT.XMAX) X=XMAX
      QSQ=SCALE**2
      IF(QSQ.LT.QSQMIN) QSQ=QSQMIN
      IF(QSQ.GT.QSQMAX) QSQ=QSQMAX
      XXX=X
      IF(X.LT.1.D-1) XXX=DLOG10(X)+1.1D0
      N=0
  70  N=N+1
      IF(XXX.GT.XX(N+1)) GOTO 70
      A=(XXX-XX(N))/(XX(N+1)-XX(N))
      RM=DLOG(QSQ/QSQMIN)/DLOG(2.D0)
      B=RM-DINT(RM)
      M=1+IDINT(RM)
      DO 60 I=1,6
      G(I)= (1.D0-A)*(1.D0-B)*F(I,N,M)+(1.D0-A)*B*F(I,N,M+1)
     .    + A*(1.D0-B)*F(I,N+1,M)  + A*B*F(I,N+1,M+1)
      IF(N.GE.16) GOTO 65
      IF(I.EQ.6) GOTO 65
          FAC=(1.D0-B)*F(I,16,M)+B*F(I,16,M+1)
 
          G(I)=FAC**(G(I)/FAC)
  65  CONTINUE
      G(I)=G(I)*(1.D0-X)**N0(I)
  60  CONTINUE
      UPV=G(1)
      DNV=G(2)
      SEA=G(4)
      STR=G(4)
      CHM=G(5)
      GLU=G(3)
      BOT=G(6)
      RETURN
      END
      SUBROUTINE STRUC3(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GLU)
C :::::::::::: MARTIN ROBERTS STIRLING :::::1/RX GLUE::::178 MEV::::
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(6,42,19),G(6),XX(42),N0(6)
      SAVE
      DATA XX/1.D-4,2.D-4,4.D-4,6.D-4,8.D-4,
     .       1.D-3,2.D-3,4.D-3,6.D-3,8.D-3,
     .       1.D-2,2.D-2,4.D-2,6.D-2,8.D-2,
     .     .1D0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
     .     .5D0,.525D0,.55D0,.575D0,.6D0,.65D0,.7D0,.75D0,
     .     .8D0,.9D0,1.D0/
      DATA XMIN,XMAX,QSQMIN,QSQMAX/1.D-4,1.D0,5.D0,1310720.D0/
      DATA N0/3,4,5,5,0,0/
      DATA INIT/0/
      IF(INIT.NE.0) GOTO 10
      OPEN(UNIT=33,FILE='mrs3.dat'
     &     ,STATUS='OLD') !,READONLY)
      READ(33,*)
      INIT=1
      DO 20 N=1,41
      DO 20 M=1,19
      READ(33,50)F(1,N,M),F(2,N,M),F(3,N,M),F(4,N,M),F(5,N,M),F(6,N,M)
         DO 25 I=1,6
  25     F(I,N,M)=F(I,N,M)/(1.D0-XX(N))**N0(I)
  20  CONTINUE
      DO 30 J=1,15
      XX(J)=DLOG10(XX(J))+1.1D0
      DO 30 I=1,5
      DO 30 K=1,19
  30  F(I,J,K)=DLOG(F(I,J,K))*F(I,16,K)/DLOG(F(I,16,K))
  50  FORMAT(6F10.5)
      DO 40 I=1,6
      DO 40 M=1,19
  40  F(I,42,M)=0.D0
      CLOSE(33)
  10  CONTINUE
      IF(X.LT.XMIN) X=XMIN
      IF(X.GT.XMAX) X=XMAX
      QSQ=SCALE**2
      IF(QSQ.LT.QSQMIN) QSQ=QSQMIN
      IF(QSQ.GT.QSQMAX) QSQ=QSQMAX
      XXX=X
      IF(X.LT.1.D-1) XXX=DLOG10(X)+1.1D0
      N=0
  70  N=N+1
      IF(XXX.GT.XX(N+1)) GOTO 70
      A=(XXX-XX(N))/(XX(N+1)-XX(N))
      RM=DLOG(QSQ/QSQMIN)/DLOG(2.D0)
      B=RM-DINT(RM)
      M=1+IDINT(RM)
      DO 60 I=1,6
      G(I)= (1.D0-A)*(1.D0-B)*F(I,N,M)+(1.D0-A)*B*F(I,N,M+1)
     .    + A*(1.D0-B)*F(I,N+1,M)  + A*B*F(I,N+1,M+1)
      IF(N.GE.16) GOTO 65
      IF(I.EQ.6) GOTO 65
          FAC=(1.D0-B)*F(I,16,M)+B*F(I,16,M+1)
 
          G(I)=FAC**(G(I)/FAC)
  65  CONTINUE
      G(I)=G(I)*(1.D0-X)**N0(I)
  60  CONTINUE
      UPV=G(1)
      DNV=G(2)
      SEA=G(4)
      STR=G(4)
      CHM=G(5)
      GLU=G(3)
      BOT=G(6)
      RETURN
      END
C Received from J. Stirling February 1989
      SUBROUTINE MRSEB(X,SCALE,MODE,UPV,DNV,SEA,STR,CHM,BOT,GLU)
C***************************************************************C
C                                                               C
C                                                               C
C     NEW VERSIONS !!!! JANUARY 1989  (AS DESCRIBED IN          C
C     "IMPROVED PARTON DISTRIBUTIONS ... " A.D. MARTIN,         C
C     R.G. ROBERTS AND W.J. STIRLING PREPRINT RAL-88-113 )      C
C                                                               C
C  MODE 1 CORRESPONDS TO                                        C
C  MARTIN, ROBERTS, STIRLING (EMC FIT)    WITH LAMBDA= 100 MEV  C
C                                                               C
C  MODE 2  CORRESPONDS TO                                       C
C  MARTIN, ROBERTS, STIRLING (BCDMS FIT)  WITH LAMBDA= 200 MEV  C
C                                                               C
C  (  SOFT GLUE :  X G(X,Q0) = A (1-X)**4.4 )                   C
C                                                               C
C                         -*-                                   C
C                                                               C
C    (NOTE THAT X TIMES THE PARTON DISTRIBUTION FUNCTION        C
C    IS RETURNED I.E. G(X) = GLU/X ETC, AND THAT "SEA"          C
C    IS THE LIGHT QUARK SEA I.E. UBAR(X)=DBAR(X)=  ...          C
C    = SEA/X FOR A PROTON.  IF IN DOUBT, CHECK THE              C
C    MOMENTUM SUM RULE! NOTE ALSO THAT SCALE=Q IN GEV)          C
C                                                               C
C                         -*-                                   C
C                                                               C
C     (THE RANGE OF APPLICABILITY IS CURRENTLY:                 C
C     10**-4 < X < 1  AND  5 < Q**2 < 1.31 * 10**6              C
C     HIGHER Q**2 VALUES CAN BE SUPPLIED ON REQUEST             C
C     - PROBLEMS, COMMENTS ETC TO SRG$T3@GEN                    C
C                                                               C
C                                                               C
C***************************************************************C
      IMPLICIT REAL*8(A-H,O-Z)
      IF(MODE.EQ.4) CALL STRUCE(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GLU)
      IF(MODE.EQ.5) CALL STRUCB(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GLU)
      RETURN
      END
C
      SUBROUTINE STRUCE(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GLU)
C :::::::::::: MARTIN ROBERTS STIRLING :::::SOFT GLUE:::: 91 MEV:::
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(6,42,19),G(6),XX(42),N0(6)
      DATA XX/1.D-4,2.D-4,4.D-4,6.D-4,8.D-4,
     .       1.D-3,2.D-3,4.D-3,6.D-3,8.D-3,
     .       1.D-2,2.D-2,4.D-2,6.D-2,8.D-2,
     .     .1D0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
     .     .5D0,.525D0,.55D0,.575D0,.6D0,.65D0,.7D0,.75D0,
     .     .8D0,.9D0,1.D0/
      DATA XMIN,XMAX,QSQMIN,QSQMAX/1.D-4,1.D0,5.D0,1310720.D0/
      DATA N0/2,5,4,5,0,0/
      DATA INIT/0/
      IF(INIT.NE.0) GOTO 10
      OPEN(UNIT=34,FILE='mrse.dat'
     &     ,STATUS='OLD') !,READONLY)
      READ(34,*)
      INIT=1
      DO 20 N=1,41
      DO 20 M=1,19
      READ(34,50)F(1,N,M),F(2,N,M),F(3,N,M),F(4,N,M),F(5,N,M),F(6,N,M)
         DO 25 I=1,6
  25     F(I,N,M)=F(I,N,M)/(1.D0-XX(N))**N0(I)
  20  CONTINUE
      DO 30 J=1,15
      XX(J)=DLOG10(XX(J))+1.1D0
      DO 30 I=1,5
      DO 30 K=1,19
  30  F(I,J,K)=DLOG(F(I,J,K))*F(I,16,K)/DLOG(F(I,16,K))
  50  FORMAT(6F10.5)
      DO 40 I=1,6
      DO 40 M=1,19
  40  F(I,42,M)=0.D0
      CLOSE (34)
  10  CONTINUE
      IF(X.LT.XMIN) X=XMIN
      IF(X.GT.XMAX) X=XMAX
      QSQ=SCALE**2
      IF(QSQ.LT.QSQMIN) QSQ=QSQMIN
      IF(QSQ.GT.QSQMAX) QSQ=QSQMAX
      XXX=X
      IF(X.LT.1.D-1) XXX=DLOG10(X)+1.1D0
      N=0
  70  N=N+1
      IF(XXX.GT.XX(N+1)) GOTO 70
      A=(XXX-XX(N))/(XX(N+1)-XX(N))
      RM=DLOG(QSQ/QSQMIN)/DLOG(2.D0)
      B=RM-DINT(RM)
      M=1+IDINT(RM)
      DO 60 I=1,6
      G(I)= (1.D0-A)*(1.D0-B)*F(I,N,M)+(1.D0-A)*B*F(I,N,M+1)
     .    + A*(1.D0-B)*F(I,N+1,M)  + A*B*F(I,N+1,M+1)
      IF(N.GE.16) GOTO 65
      IF(I.EQ.6) GOTO 65
          FAC=(1.D0-B)*F(I,16,M)+B*F(I,16,M+1)
          G(I)=FAC**(G(I)/FAC)
  65  CONTINUE
      G(I)=G(I)*(1.D0-X)**N0(I)
  60  CONTINUE
      UPV=G(1)
      DNV=G(2)
      SEA=G(4)
      STR=G(4)
      CHM=G(5)
      GLU=G(3)
      BOT=G(6)
      RETURN
      END
      SUBROUTINE STRUCB(X,SCALE,UPV,DNV,SEA,STR,CHM,BOT,GLU)
C :::::::::::: MARTIN ROBERTS STIRLING :::::SOFT GLUE::::228 MEV:::
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(6,42,19),G(6),XX(42),N0(6)
      DATA XX/1.D-4,2.D-4,4.D-4,6.D-4,8.D-4,
     .       1.D-3,2.D-3,4.D-3,6.D-3,8.D-3,
     .       1.D-2,2.D-2,4.D-2,6.D-2,8.D-2,
     .     .1D0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
     .     .5D0,.525D0,.55D0,.575D0,.6D0,.65D0,.7D0,.75D0,
     .     .8D0,.9D0,1.D0/
      DATA XMIN,XMAX,QSQMIN,QSQMAX/1.D-4,1.D0,5.D0,1310720.D0/
      DATA N0/2,5,4,5,0,0/
      DATA INIT/0/
      IF(INIT.NE.0) GOTO 10
      OPEN(UNIT=35,FILE='mrsb.dat'
     &     ,STATUS='OLD') !,READONLY)
      READ(35,*)
      INIT=1
      DO 20 N=1,41
      DO 20 M=1,19
      READ(35,50)F(1,N,M),F(2,N,M),F(3,N,M),F(4,N,M),F(5,N,M),F(6,N,M)
         DO 25 I=1,6
  25     F(I,N,M)=F(I,N,M)/(1.D0-XX(N))**N0(I)
  20  CONTINUE
      DO 30 J=1,15
      XX(J)=DLOG10(XX(J))+1.1D0
      DO 30 I=1,5
      DO 30 K=1,19
  30  F(I,J,K)=DLOG(F(I,J,K))*F(I,16,K)/DLOG(F(I,16,K))
  50  FORMAT(6F10.5)
      DO 40 I=1,6
      DO 40 M=1,19
  40  F(I,42,M)=0.D0
      CLOSE (35)
  10  CONTINUE
      IF(X.LT.XMIN) X=XMIN
      IF(X.GT.XMAX) X=XMAX
      QSQ=SCALE**2
      IF(QSQ.LT.QSQMIN) QSQ=QSQMIN
      IF(QSQ.GT.QSQMAX) QSQ=QSQMAX
      XXX=X
      IF(X.LT.1.D-1) XXX=DLOG10(X)+1.1D0
      N=0
  70  N=N+1
      IF(XXX.GT.XX(N+1)) GOTO 70
      A=(XXX-XX(N))/(XX(N+1)-XX(N))
      RM=DLOG(QSQ/QSQMIN)/DLOG(2.D0)
      B=RM-DINT(RM)
      M=1+IDINT(RM)
      DO 60 I=1,6
      G(I)= (1.D0-A)*(1.D0-B)*F(I,N,M)+(1.D0-A)*B*F(I,N,M+1)
     .    + A*(1.D0-B)*F(I,N+1,M)  + A*B*F(I,N+1,M+1)
      IF(N.GE.16) GOTO 65
      IF(I.EQ.6) GOTO 65
          FAC=(1.D0-B)*F(I,16,M)+B*F(I,16,M+1)
          G(I)=FAC**(G(I)/FAC)
  65  CONTINUE
      G(I)=G(I)*(1.D0-X)**N0(I)
  60  CONTINUE
      UPV=G(1)
      DNV=G(2)
      SEA=G(4)
      STR=G(4)
      CHM=G(5)
      GLU=G(3)
      BOT=G(6)
      RETURN
      END
