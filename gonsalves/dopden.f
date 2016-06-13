C-----------------------------------------------------------------------
C     Duke-Owens' Q**2 dependent parton densities
C     Phys. Rev. D 30, 49 (1984)
C
      SUBROUTINE dopden (nset,x,qsq,pden)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION pden(-6:6)
 
C---  Q for start of evolution
      q0=2.
 
C---  The output distributions for the proton are
C     pden(0) = gluon
C     pden(1) = valence up quark + up sea
C     pden(2) = valence down quark + down sea
C     pden(3) = strange sea
C     pden(4) = charm sea
C     pden(5) = bottom sea
C     pden(6) = top sea
C     pden(-1) = anti-up sea
C     pden(-2) = anti-down sea
C     pden(-3) = anti-strange sea
C     pden(-4) = anti-charm sea
C     pden(-5) = anti-bottom sea
C     pden(-6) = anti-top sea

      IF (nset.LT.1.OR.nset.GT.2) THEN
          PRINT *,'Bad Duke Owens Set ',nset
          STOP
      ENDIF 
      GOTO (1,2) nset
1     CONTINUE
C---  nset = 1  :  set 1  Lambda = 0.2 GeV/c
 
      xlambd=0.2
      s=dlog(qsq/xlambd**2)
      s=s/(2.*dlog(q0/xlambd))
      s=dlog(s)
 
C     valence quark parameters
 
      eta1=0.419+0.004*s-0.007*s**2
      eta2=3.46+0.724*s-0.066*s**2
      gamud=4.40-4.86*s+1.33*s**2
      eta3=0.763-0.237*s+0.026*s**2
      eta4=4.00+0.627*s-0.019*s**2
      gamd=-0.421*s+0.033*s**2
 
C     up-down-strange sea quark parameters
 
      saa=1.265-1.132*s+0.293*s**2
      sa=-0.372*s-0.029*s**2
      sb=8.05+1.59*s-0.153*s**2
      salp=6.31*s-0.273*s**2
      sbet=-10.5*s-3.17*s**2
      sgam=14.7*s+9.80*s**2
 
C     charmed quark parameters
 
      chaa=0.135*s-0.075*s**2
      cha=-0.036-0.222*s-0.058*s**2
      chb=6.35+3.26*s-0.909*s**2
      chalp=-3.03*s+1.50*s**2
      chbet=17.4*s-11.3*s**2
      chgam=-17.9*s+15.6*s**2
 
C     gluon parameters
 
      glaa=1.56-1.71*s+0.638*s**2
      gla=-0.949*s+0.325*s**2
      glb=6.0+1.44*s-1.05*s**2
      glalp=9.0-7.19*s+0.255*s**2
      glbet=-16.5*s+10.9*s**2
      glgam=15.3*s-10.1*s**2
 
      GOTO 3
 
2     CONTINUE
C---  nset = 2  :  set 2  Lambda = 0.4 GeV/c
 
      xlambd=0.4
      s=dlog(qsq/xlambd**2)
      s=s/(2.*dlog(q0/xlambd))
      s=dlog(s)
 
C     valence quark parameters
 
      eta1=0.374+0.014*s
      eta2=3.33+0.753*s-0.076*s**2
      gamud=6.03-6.22*s+1.56*s**2
      eta3=0.761-0.232*s+0.023*s**2
      eta4=3.83+0.627*s-0.019*s**2
      gamd=-0.418*s+0.036*s**2
 
C     up-down-strange sea quark parameters
 
      saa=1.67-1.92*s+0.582*s**2
      sa=-0.273*s-0.164*s**2
      sb=9.15+0.530*s-0.763*s**2
      salp=15.7*s-2.83*s**2
      sbet=-101.*s+44.7*s**2
      sgam=223.*s-117.*s**2
 
C     charmed quark parameters
 
      chaa=0.067*s-0.031*s**2
      cha=-0.120-0.233*s-0.023*s**2
      chb=3.51+3.66*s-0.453*s**2
      chalp=-0.474*s+0.358*s**2
      chbet=9.50*s-5.43*s**2
      chgam=-16.6*s+15.5*s**2
 
C     gluon parameters
 
      glaa=0.879-0.971*s+0.434*s**2
      gla=-1.16*s+0.476*s**2
      glb=4.0+1.23*s-0.254*s**2
      glalp=9.0-5.64*s-0.817*s**2
      glbet=-7.54*s+5.50*s**2
      glgam=-0.596*s+0.126*s**2
 
3     CONTINUE
 
C     gluon density
 
      xglu=glaa*x**gla*(1.-x)**glb*(1.+glalp*x+glbet*x**2+glgam*x**3)
      pden(0)=xglu/x
 
C     valence quark densities
 
      xnud=3./beta(eta1,eta2+1.)
      xnud=xnud/(1.+gamud*eta1/(eta1+eta2+1.))
      xnd=1./beta(eta3,eta4+1.)
      xnd=xnd/(1.+gamd*eta3/(eta3+eta4+1.))
 
      xuvdv=xnud*x**eta1*(1.-x)**eta2*(1.+gamud*x)
      xdv=xnd*x**eta3*(1.-x)**eta4*(1.+gamd*x)
      pden(1)=(xuvdv-xdv)/x
      pden(2)=xdv/x
 
C     up-down-strange sea quark densities
 
      xsea=saa*x**sa*(1.-x)**sb*(1.+salp*x+sbet*x**2+sgam*x**3)
      pden(-1)=xsea/x/6.0
      pden(-2)=pden(3)
      pden(-3)=pden(3)
 
C     charmed quark density
 
      xchm=chaa*x**cha*(1.-x)**chb*(1.+chalp*x+chbet*x**2+chgam*x**3)
      pden(-4)=xchm/x
 
C---  top and bottom are not implemented - return 0
      pden(-5)=0.
      pden(-6)=0.

C---  rearrange
      pden(1)=pden(1)+pden(-1)
      pden(2)=pden(2)+pden(-2)
      do i=3,6
          pden(i)=pden(-i)
      enddo
 
      RETURN
      END
C-----------------------------------------------------------------------
C     Logarithm of the gamma function - numerical recipes
C
      DOUBLE PRECISION FUNCTION gammln (xx)
C---  Routine supposed to be accurate for xx > 1.
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION cof(6)
      DATA cof,stp /76.18009173d0,-86.50532033d0,24.01409822d0,
     &       -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
      DATA half,one,fpf /0.5d0,1.0d0,5.5d0/
      x=xx-one
      tmp=x+fpf
      tmp=(x+half)*dlog(tmp)-tmp
      ser=one
      DO 11 j=1,6
          x=x+one
          ser=ser+cof(j)/x
11    CONTINUE
      gammln=tmp+dlog(stp*ser)
      RETURN
      END
C-----------------------------------------------------------------------
C     The beta function : both arguments positive
C
      DOUBLE PRECISION FUNCTION beta (xx,yy)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DATA one /1.0d0/
      x=xx
      y=yy
      IF (xx.LT.one) x=x+one
      IF (yy.LT.one) y=y+one
      beta=gammln(x)+gammln(y)-gammln(x+y)
      beta=dexp(beta)
      IF (xx.LT.one) THEN
          beta=beta*(one+yy/xx)
          IF (yy.LT.one) beta=beta*(one+x/yy)
      ELSE
          IF (yy.LT.one) beta=beta*(one+xx/yy)
      ENDIF
      RETURN
      END
