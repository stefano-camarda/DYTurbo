*     Extracted from the MRST demonstration program ggh.f
      FUNCTION ALPHAS(SCALE,LAMBDA,IORDER)
C.
C.   The QCD coupling used in the MRS analysis.
C.   Automatically corrects for NF=3,4,5 with thresholds at m_Q.
C.   The following parameters must be supplied:
C.       SCALE   =  the QCD scale in GeV (real)
C.       LAMBDA  =  the 4 flavour MSbar Lambda parameter in GeV (real)
C.       IORDER  =  0 for LO, 1 for NLO and 2 for NNLO
C. 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 LAMBDA
      DATA TOL,QSCT,QSDT/5D-4,74D0,8.18D0/
      PI=4D0*DATAN(1D0)
      IORD=IORDER
      ITH=0
      FLAV=4D0
      AL=LAMBDA
      AL2=LAMBDA*LAMBDA
      QS=SCALE*SCALE
      T=DLOG(QS/AL2)
      TT=T
C.
      qsctt=qsct/4.
      qsdtt=qsdt/4.
      if(qs.lt.0.5d0) then   !!  running stops below 0.5
          qs=0.5d0
          t=dlog(qs/al2)
          tt=t
      endif

      IF(QS.gt.QSCTT) GO	TO 12  
      IF(QS.lt.QSDTT) GO	TO 312  
   11 CONTINUE
      B0=11-2.*FLAV/3. 
      X1=4.*PI/B0
  5   continue    
      IF(IORD.eq.0) then
      ALPHAS=X1/T
      ELSEIF(IORD.eq.1) then
      B1=102.-38.*FLAV/3.
      X2=B1/B0**2
      AS2=X1/T*(1.-X2*dlog(T)/T)
  95    AS=AS2
        F=-T+X1/AS-X2*dlog(X1/AS+X2)
        FP=-X1/AS**2*(1.-X2/(X1/AS+X2))
        AS2=AS-F/FP
        DEL=ABS(F/FP/AS)
        IF((DEL-TOL).GT.0.) go to 95
      ALPHAS=AS2
      ELSEIF(IORD.eq.2) then
      ALPHAS=qwikalf(t,2,flav)
      ENDIF
      IF(ITH.EQ.0) RETURN
      GO TO (13,14,15) ITH
      print *, 'here'
      GO TO 5
   12 ITH=1
      T=dlog(QSCTT/AL2)
      GO TO 11
   13 ALFQC4=ALPHAS
      FLAV=5.   
      ITH=2
      GO TO 11
   14 ALFQC5=ALPHAS
      ITH=3
      T=TT
      GO TO 11
   15 ALFQS5=ALPHAS
      ALFINV=1./ALFQS5+1./ALFQC4-1./ALFQC5
      ALPHAS=1./ALFINV
      RETURN

  311 CONTINUE
      B0=11-2.*FLAV/3. 
      X1=4.*PI/B0
      IF(IORD.eq.0) then
      ALPHAS=X1/T
      ELSEIF(IORD.eq.1) then
      B1=102.-38.*FLAV/3.
      X2=B1/B0**2
      AS2=X1/T*(1.-X2*dlog(T)/T)
   35 AS=AS2
      F=-T+X1/AS-X2*dlog(X1/AS+X2)
      FP=-X1/AS**2*(1.-X2/(X1/AS+X2))
      AS2=AS-F/FP
      DEL=ABS(F/FP/AS)
      IF((DEL-TOL).GT.0.) go to 35
      ELSEIF(IORD.eq.2) then
      ALPHAS=qwikalf(t,2,flav)
      ENDIF
      IF(ITH.EQ.0) RETURN
      GO TO (313,314,315) ITH
  312 ITH=1
      T=dlog(QSDTT/AL2)
      GO TO 311
  313 ALFQC4=ALPHAS
      FLAV=3.   
      ITH=2
      GO TO 311
  314 ALFQC3=ALPHAS
      ITH=3
      T=TT
      GO TO 311
  315 ALFQS3=ALPHAS
      ALFINV=1./ALFQS3+1./ALFQC4-1./ALFQC3
      ALPHAS=1./ALFINV
      RETURN
      END
      function qwikalf(t,iord,flav)
      implicit real*8(a-h,o-z)
      dimension z3(6),z4(6),z5(6),zz3(6),zz4(6),zz5(6)
      data z3/ -.161667E+01,0.954244E+01,
     .0.768623E+01,0.101523E+00,-.360127E-02,0.457867E-04/
      data z4/ -.172239E+01,0.831185E+01,
     .0.721463E+01,0.835531E-01,-.285436E-02,0.349129E-04/
      data z5/ -.872190E+00,0.572816E+01,
     .0.716119E+01,0.195884E-01,-.300199E-03,0.151741E-05/
      data zz3/-.155611E+02,0.168406E+02,
     .0.603014E+01,0.257682E+00,-.970217E-02,0.127628E-03/
      data zz4/-.106762E+02,0.118497E+02,0.664964E+01,
     .0.112996E+00,-.317551E-02,0.302434E-04/
      data zz5/-.531860E+01,0.708503E+01,0.698352E+01,
     .0.274170E-01,-.426894E-03,0.217591E-05/
      data pi/3.14159/
      nfm2=flav-2.
      x=dsqrt(t)
      x2=x*x
      x3=x*x2
      x4=x*x3
      x5=x*x4
      go to (1,2) iord  !!iord change!!
    1 go to (3,4,5) nfm2
    3 y=z3(1)+z3(2)*x+z3(3)*x2+z3(4)*x3+z3(5)*x4+z3(6)*x5
      go to 10
    4 y=z4(1)+z4(2)*x+z4(3)*x2+z4(4)*x3+z4(5)*x4+z4(6)*x5
      go to 10
    5 y=z5(1)+z5(2)*x+z5(3)*x2+z5(4)*x3+z5(5)*x4+z5(6)*x5
      go to 10
    2 go to (6,7,8) nfm2
    6 y=zz3(1)+zz3(2)*x+zz3(3)*x2+zz3(4)*x3+zz3(5)*x4+zz3(6)*x5
      go to 10
    7 y=zz4(1)+zz4(2)*x+zz4(3)*x2+zz4(4)*x3+zz4(5)*x4+zz4(6)*x5
      go to 10
    8 y=zz5(1)+zz5(2)*x+zz5(3)*x2+zz5(4)*x3+zz5(5)*x4+zz5(6)*x5
      go to 10
   10 qwikalf=4.*pi/y
      return
      end


