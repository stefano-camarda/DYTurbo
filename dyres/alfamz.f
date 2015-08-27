      DOUBLE PRECISION FUNCTION dyALPHAS(Q,AMZ,NLOOP)
c--- this is simply a wrapper to the LHAPDF implementation of the
c--- running coupling, in the style of the native MCFM routine

c--- Note that the inputs AMZ and NLOOP are not used
      IMPLICIT NONE
      DOUBLE PRECISION Q,AMZ,alphasPDF
      INTEGER NLOOP
            
      dyALPHAS=alphasPDF(Q)

      RETURN
      END

c      DOUBLE PRECISION FUNCTION dyALPHAS(Q,AMZ,NLOOP)
cc     Evaluation of strong coupling constant alpha_S
cc     Author: R.K. Ellis
c
cc     q -- scale at which alpha_s is to be evaluated
cc     amz -- value of alpha_s at the mass of the Z-boson
cc     nloop -- the number of loops (1,2, or 3) at which beta 
cc     function is evaluated to determine running.
cc     the values of the cmass and the bmass should be set
cc     in common block qmass.
c
c      IMPLICIT NONE
c      DOUBLE PRECISION Q,T,AMZ,AMZ0,AMB,AMC,ZMASS,BMASS,CMASS,AS_OUT
c      INTEGER NLOOP,NLOOP0,NF3,NF4,NF5
c      PARAMETER(ZMASS=91.188D0)
c      PARAMETER(NF5=5,NF4=4,NF3=3)
c      COMMON/QMASS/CMASS,BMASS
c      SAVE AMZ0,NLOOP0,AMB,AMC
c      DATA AMZ0,NLOOP0/0D0,0/
c
c      IF (Q .LE. 0D0) THEN 
c         WRITE(6,*) 'q .le. 0 in alphas'
c         WRITE(6,*) 'q= ',Q
c         STOP
c      ENDIF
c      IF (AMZ .LE. 0D0) THEN 
c         WRITE(6,*) 'amz .le. 0 in alphas',AMZ
c         STOP
c      ENDIF
c      IF (CMASS .LE. 0.3D0) THEN 
c         WRITE(6,*) 'cmass .le. 0.3GeV in alphas',CMASS
c         WRITE(6,*) 'COMMON/QMASS/CMASS,BMASS'
c         WRITE(6,*) 'continue with cmass=1.5GeV'
c         CMASS=1.5D0
c      ENDIF
c      IF (BMASS .LE. 0D0) THEN 
c         WRITE(6,*) 'bmass .le. 0 in alphas',BMASS
c         WRITE(6,*) 'COMMON/QMASS/CMASS,BMASS'
c         WRITE(6,*) 'continue with bmass=5.0GeV'
c         BMASS=5D0
c      ENDIF
cc--- establish value of coupling at b- and c-mass and save
c      IF ((AMZ .NE. AMZ0) .OR. (NLOOP .NE. NLOOP0)) THEN
c         AMZ0=AMZ
c         NLOOP0=NLOOP
c         T=2D0*DLOG(BMASS/ZMASS)
c         CALL NEWTON1(T,AMZ,AMB,NLOOP,NF5)
c         T=2D0*DLOG(CMASS/BMASS)
c         CALL NEWTON1(T,AMB,AMC,NLOOP,NF4)
c      ENDIF
c
cc--- evaluate strong coupling at scale q
c      IF (Q  .LT. BMASS) THEN
c           IF (Q  .LT. CMASS) THEN
c             T=2D0*DLOG(Q/CMASS)
c             CALL NEWTON1(T,AMC,AS_OUT,NLOOP,NF3)
c           ELSE
c             T=2D0*DLOG(Q/BMASS)
c             CALL NEWTON1(T,AMB,AS_OUT,NLOOP,NF4)
c           ENDIF
c      ELSE
c      T=2D0*DLOG(Q/ZMASS)
c      CALL NEWTON1(T,AMZ,AS_OUT,NLOOP,NF5)
c      ENDIF
c      dyALPHAS=AS_OUT
c      RETURN
c      END
c
c
c      SUBROUTINE DIFF(Q,AMZ,NLOOP)
c      IMPLICIT NONE
c      DOUBLE PRECISION BETA(3:5),B0(3:5),C1(3:5),C2(3:5)
c      INTEGER NLOOP,J
c      DOUBLE PRECISION Q,QP,QM,AMZ,CMASS,BMASS,X1,X2,X3,EP,DIFF1
c     .     ,dyALPHAS
c      COMMON/QMASS/CMASS,BMASS
cC---     B0=(11.-2.*F/3.)/4./PI
c      DATA B0/0.716197243913527D0,0.66314559621623D0,0.61009394851893D0/
cC---     C1=(102.D0-38.D0/3.D0*F)/4.D0/PI/(11.D0-2.D0/3.D0*F)
c      DATA C1/.565884242104515D0,0.49019722472304D0,0.40134724779695D0/
cC---     C2=(2857.D0/2.D0-5033*F/18.D0+325*F**2/54)
cC---     /16.D0/PI**2/(11.D0-2.D0/3.D0*F)
c      DATA C2/0.453013579178645D0,0.30879037953664D0,0.14942733137107D0/
cC---     DEL=SQRT(4*C2-C1**2)
c
c      X1=dyALPHAS(Q,AMZ,1)
c      X2=dyALPHAS(Q,AMZ,2)
c      X3=dyALPHAS(Q,AMZ,3)
c      J=3
c      IF (Q .GT. CMASS) J=4
c      IF (Q .GT. BMASS) J=5
c      EP=.001D0
c      QP=Q*(1D0+EP)
c      QM=Q*(1D0-EP)
c      IF (NLOOP .EQ.1) THEN 
c      BETA(J)=-B0(J)*X1**2
c      DIFF1=(dyALPHAS(QP,AMZ,1)-dyALPHAS(QM,AMZ,1))/4d0/EP/BETA(J)
c      ENDIF
c      IF (NLOOP .EQ.2) THEN 
c      BETA(J)=-B0(J)*X2**2*(1D0+C1(J)*X2)
c      DIFF1=(dyALPHAS(QP,AMZ,2)-dyALPHAS(QM,AMZ,2))/4d0/EP/BETA(J)
c      ENDIF
c      IF (NLOOP .EQ.3) THEN 
c      BETA(J)=-B0(J)*X3**2*(1D0+C1(J)*X3+C2(J)*X3**2)
c      DIFF1=(dyALPHAS(QP,AMZ,3)-dyALPHAS(QM,AMZ,3))/4d0/EP/BETA(J)
c      ENDIF
c      WRITE(6,*) Q,DIFF1,NLOOP
c      RETURN
c      END
c
cC      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
cC      DOUBLE PRECISION B0(3:5),C1(3:5),C2(3:5),DEL(3:5)
cC      PARAMETER(PI=3.1415926535898D0)
cC      NLOOP=2
cC      AMZ=0.113D0
cC      DO N=3,5
cC      F=DFLOAT(N)
cC      B0(N)=(11.D0-2.D0*F/3.D0)/4.D0/PI
cC      C1(N)=(102.D0-38.D0/3.D0*F)/4.D0/PI/(11.D0-2.D0/3.D0*F)
cC      C2(N)=(2857.D0/2.D0-5033*F/18.D0+325D0*F**2/54D0)
cC     &   /16D0/PI**2/(11.D0-2D0/3D0*F)
cC      DEL(N)=SQRT(4D0*C2(N)-C1(N)**2)
cC      ENDDO
cC      OPEN(UNIT=67,FILE='TEMP.DAT')
cC      WRITE(67,*) B0
cC      WRITE(67,*) C1
cC      WRITE(67,*) C2
cC      WRITE(67,*) DEL
cC      DO N=1,100
cC      Q=DFLOAT(N)+0.1
cC      WRITE(6,*)
cC      CALL DIFF(Q,AMZ,1)
cC      CALL DIFF(Q,AMZ,2)
cC      CALL DIFF(Q,AMZ,3)
cC      ENDDO
cC      STOP
cC      END


