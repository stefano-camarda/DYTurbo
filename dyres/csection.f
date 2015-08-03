      function xsection(ax,eta,xmuf,ippbar1,ippbar2)
!!!      function xsection(ax,eta,xmuf2,isetproton,ippbar1,ippbar2)
      IMPLICIT NONE
      dimension sigmaij(-5:5,-5:5)
      DOUBLE PRECISION ax,ax1,ax2,eta,xmuf,xmuf2, FX1(-5:5),FX2(-5:5),
     . sigmaij,x1,x2,QQBN,xsection,FX1NEW(-2:2),FX2NEW(-2:2)
      integer isetproton,ippbar1,ippbar2,nf,si,sj
      common/sigmaij/sigmaij
      COMMON/NFLAVORS/nF
        AX1 = (AX+2*ETA)/2d0 
        AX2 = (AX-2*ETA)/2d0 
        x1=dExp(ax1)
        x2=dExp(ax2)
!      call partons(xmuf2,x1,fx1,5,isetproton,ippbar1)
!      call partons(xmuf2,x2,fx2,5,isetproton,ippbar2)
      call fdist(ippbar1,x1,xmuf,fx1)
      call fdist(ippbar2,x2,xmuf,fx2)

c!!!!   changed definition U <-> D 

      FX1NEW(1)=FX1(2)
      FX1NEW(-1)=FX1(-2)
      FX1NEW(2)=FX1(1)
      FX1NEW(-2)=FX1(-1)
      FX2NEW(1)=FX2(2)
      FX2NEW(-1)=FX2(-2)
      FX2NEW(2)=FX2(1)
      FX2NEW(-2)=FX2(-1)

      FX1(1)=FX1NEW(1)
      FX1(-1)=FX1NEW(-1)
      FX1(2)=FX1NEW(2)
      FX1(-2)=FX1NEW(-2)
      FX2(1)=FX2NEW(1)
      FX2(-1)=FX2NEW(-1)
      FX2(2)=FX2NEW(2)
      FX2(-2)=FX2NEW(-2)

c!!!

      QQBN=0
      
      do si=-nf,nf
      do sj=-nf,nf
C Q Qb + Qb Q
       QQBN = QQBN +  FX1(si)*FX2(sj)*sigmaij(si,sj)
!       write(*,*) si,sj,FX1(si),FX2(sj),sigmaij(si,sj)
      enddo
      enddo       
       
       xsection=QQBN
!       write(*,*) ax,eta,xmuf,ippbar1,ippbar2,xsection
       return
       end
      
!      function xsection2(ax,eta,xmuf2,isetproton,ippbar1,ippbar2)
      function xsection2(ax,eta,xmuf,ippbar1,ippbar2)
      IMPLICIT NONE
      dimension sigmaij(-5:5,-5:5)
      DOUBLE PRECISION ax,ax1,ax2,eta,xmuf2,xmuf, FX1(-5:5),FX2(-5:5),
     .  sigmaij,x1,x2,QQBN,xsection2
      integer isetproton,ippbar1,ippbar2,nf,si,sj
      common/sigmaij/sigmaij
      COMMON/NFLAVORS/nF
        AX1 = (AX+2*ETA)/2d0 
        AX2 = (AX-2*ETA)/2d0 
        x1=dExp(ax1)
        x2=dExp(ax2)
      call DISTRIBUTIONS1(X1,FX1(1),FX1(2),FX1(-1),FX1(-2),
     .     FX1(3),FX1(0),FX1(4),FX1(5))
        FX1(-3)=FX1(3)
        FX1(-4)=FX1(4)
        FX1(-5)=FX1(5)
      call DISTRIBUTIONS2(X2,FX2(1),FX2(2),FX2(-1),FX2(-2),
     .     FX2(3),FX2(0),FX2(4),FX2(5))
        FX2(-3)=FX2(3)
        FX2(-4)=FX2(4)
        FX2(-5)=FX2(5)
                
      QQBN=0      
      do si=-nf,nf
      do sj=-nf,nf
C Q Qb + Qb Q
       QQBN = QQBN +  FX1(si)*FX2(sj)*sigmaij(si,sj)/x1/x2
       enddo
       enddo              
       xsection2=QQBN       
       return
       end
      
c computes pdfs for beam 1 in the approximated result at muf
      SUBROUTINE DISTRIBUTIONS1(X,U,D,US,DS,SS,GL,CH,BO)
      implicit NONE
      double precision A1UV,A2UV,A3UV,A4UV,A5UV,A6UV,A7UV,A8UV
      double precision A1DV,A2DV,A3DV,A4DV,A5DV,A6DV,A7DV,A8DV
      double precision A1US,A2US,A3US,A4US,A5US,A6US,A7US,A8US 
      double precision A1DS,A2DS,A3DS,A4DS,A5DS,A6DS,A7DS,A8DS
      double precision A1SS,A2SS,A3SS,A4SS,A5SS,A6SS,A7SS,A8SS
      double precision A1GL,A2GL,A3GL,A4GL,A5GL,A6GL,A7GL,A8GL
      double precision A1CH,A2CH,A3CH,A4CH,A5CH,A6CH,A7CH,A8CH
      double precision A1BO,A2BO,A3BO,A4BO,A5BO,A6BO,A7BO,A8BO
c
       dimension A1UV(30),A2UV(30),A3UV(30),A4UV(30),
     .       A5UV(30),A6UV(30),A7UV(30),A8UV(30)
       dimension A1DV(30),A2DV(30),A3DV(30),A4DV(30),
     .       A5DV(30),A6DV(30),A7DV(30),A8DV(30)
       dimension A1US(30),A2US(30),A3US(30),A4US(30),
     .       A5US(30),A6US(30),A7US(30),A8US(30)
       dimension A1DS(30),A2DS(30),A3DS(30),A4DS(30),
     .       A5DS(30),A6DS(30),A7DS(30),A8DS(30)
       dimension A1SS(30),A2SS(30),A3SS(30),A4SS(30),
     .       A5SS(30),A6SS(30),A7SS(30),A8SS(30)
       dimension A1GL(30),A2GL(30),A3GL(30),A4GL(30),
     .       A5GL(30),A6GL(30),A7GL(30),A8GL(30)
       dimension A1CH(30),A2CH(30),A3CH(30),A4CH(30),
     .       A5CH(30),A6CH(30),A7CH(30),A8CH(30)
       dimension A1BO(30),A2BO(30),A3BO(30),A4BO(30),
     .       A5BO(30),A6BO(30),A7BO(30),A8BO(30)
c     
      double precision aa,UV,DV,US,DS,SS,GL,CH,BO,U,D,X 
      double precision UTEMP,DTEMP
      INTEGER I,ih1,ih2
      COMMON/ CUV1/ A1UV,A2UV,A3UV,A4UV,A5UV,A6UV,A7UV,A8UV
      COMMON/ CDV1/ A1DV,A2DV,A3DV,A4DV,A5DV,A6DV,A7DV,A8DV
      COMMON/ CUS1/ A1US,A2US,A3US,A4US,A5US,A6US,A7US,A8US
      COMMON/ CDS1/ A1DS,A2DS,A3DS,A4DS,A5DS,A6DS,A7DS,A8DS
      COMMON/ CSS1/ A1SS,A2SS,A3SS,A4SS,A5SS,A6SS,A7SS,A8SS
      COMMON/ CGL1/ A1GL,A2GL,A3GL,A4GL,A5GL,A6GL,A7GL,A8GL
      COMMON/ CCH1/ A1CH,A2CH,A3CH,A4CH,A5CH,A6CH,A7CH,A8CH
      COMMON/ CBO1/ A1BO,A2BO,A3BO,A4BO,A5BO,A6BO,A7BO,A8BO
      COMMON/expp/aa
      COMMON/ IFIT/ I
      COMMON/collider/ih1,ih2

       
      UV=A1UV(I)*x**A2UV(I)*(1-x)**A3UV(I)*(1+A4UV(I)*x+A5UV(I)*x**(0.5)
     .     +A6UV(I)*x**(1.5)+A7UV(I)*X**2 + A8UV(I)*X**(aa))
      DV=A1DV(I)*x**A2DV(I)*(1-x)**A3DV(I)*(1+A4DV(I)*x+A5DV(I)*x**(0.5)
     .     +A6DV(I)*x**(1.5)+A7DV(I)*X**2 + A8DV(I)*X**(aa))

      US=A1US(I)*x**A2US(I)*(1-x)**A3US(I)*(1+A4US(I)*x+A5US(I)*x**(0.5)
     .     +A6US(I)*x**(1.5)+A7US(I)*X**2 + A8US(I)*X**(aa))
      DS=A1DS(I)*x**A2DS(I)*(1-x)**A3DS(I)*(1+A4DS(I)*x+A5DS(I)*x**(0.5)
     .     +A6DS(I)*x**(1.5)+A7DS(I)*X**2 + A8DS(I)*X**(aa))

       U=UV+US
       D=DV+DS

      SS=A1SS(I)*x**A2SS(I)*(1-x)**A3SS(I)*(1+A4SS(I)*x+A5SS(I)*x**(0.5)
     .     +A6SS(I)*x**(1.5)+A7SS(I)*X**2 + A8SS(I)*X**(aa))

      GL=A1GL(I)*x**A2GL(I)*(1-x)**A3GL(I)*(1+A4GL(I)*x+A5GL(I)*x**(0.5)
     .     +A6GL(I)*x**(1.5)+A7GL(I)*X**2 + A8GL(I)*X**(aa))

      CH=A1CH(I)*x**A2CH(I)*(1-x)**A3CH(I)*(1+A4CH(I)*x+A5CH(I)*x**(0.5)
     .     +A6CH(I)*x**(1.5)+A7CH(I)*X**2 + A8CH(I)*X**(aa))

      BO=A1BO(I)*x**A2BO(I)*(1-x)**A3BO(I)*(1+A4BO(I)*x+A5BO(I)*x**(0.5)
     .     +A6BO(I)*x**(1.5)+A7BO(I)*X**2 + A8BO(I)*X**(aa))


       if (ih1.eq.-1) then
c pbar  u<->ubar  d<->dbar       
       UTEMP=U
       U=US
       US=UTEMP

       DTEMP=D
       D=DS
       DS=DTEMP
         
       endif

       RETURN
       END


c computes pdfs for beam 2 in the approximated result at muf
      SUBROUTINE DISTRIBUTIONS2(X,U,D,US,DS,SS,GL,CH,BO)
      implicit NONE
      double precision A1UV,A2UV,A3UV,A4UV,A5UV,A6UV,A7UV,A8UV
      double precision A1DV,A2DV,A3DV,A4DV,A5DV,A6DV,A7DV,A8DV
      double precision A1US,A2US,A3US,A4US,A5US,A6US,A7US,A8US 
      double precision A1DS,A2DS,A3DS,A4DS,A5DS,A6DS,A7DS,A8DS
      double precision A1SS,A2SS,A3SS,A4SS,A5SS,A6SS,A7SS,A8SS
      double precision A1GL,A2GL,A3GL,A4GL,A5GL,A6GL,A7GL,A8GL
      double precision A1CH,A2CH,A3CH,A4CH,A5CH,A6CH,A7CH,A8CH
      double precision A1BO,A2BO,A3BO,A4BO,A5BO,A6BO,A7BO,A8BO
c
       dimension A1UV(30),A2UV(30),A3UV(30),A4UV(30),
     .       A5UV(30),A6UV(30),A7UV(30),A8UV(30)
       dimension A1DV(30),A2DV(30),A3DV(30),A4DV(30),
     .       A5DV(30),A6DV(30),A7DV(30),A8DV(30)
       dimension A1US(30),A2US(30),A3US(30),A4US(30),
     .       A5US(30),A6US(30),A7US(30),A8US(30)
       dimension A1DS(30),A2DS(30),A3DS(30),A4DS(30),
     .       A5DS(30),A6DS(30),A7DS(30),A8DS(30)
       dimension A1SS(30),A2SS(30),A3SS(30),A4SS(30),
     .       A5SS(30),A6SS(30),A7SS(30),A8SS(30)
       dimension A1GL(30),A2GL(30),A3GL(30),A4GL(30),
     .       A5GL(30),A6GL(30),A7GL(30),A8GL(30)
       dimension A1CH(30),A2CH(30),A3CH(30),A4CH(30),
     .       A5CH(30),A6CH(30),A7CH(30),A8CH(30)
       dimension A1BO(30),A2BO(30),A3BO(30),A4BO(30),
     .       A5BO(30),A6BO(30),A7BO(30),A8BO(30)
c     
      double precision aa,UV,DV,US,DS,SS,GL,CH,BO,U,D,X 
      double precision UTEMP,DTEMP
      INTEGER I,ih1,ih2
      COMMON/ CUV2/ A1UV,A2UV,A3UV,A4UV,A5UV,A6UV,A7UV,A8UV
      COMMON/ CDV2/ A1DV,A2DV,A3DV,A4DV,A5DV,A6DV,A7DV,A8DV
      COMMON/ CUS2/ A1US,A2US,A3US,A4US,A5US,A6US,A7US,A8US
      COMMON/ CDS2/ A1DS,A2DS,A3DS,A4DS,A5DS,A6DS,A7DS,A8DS
      COMMON/ CSS2/ A1SS,A2SS,A3SS,A4SS,A5SS,A6SS,A7SS,A8SS
      COMMON/ CGL2/ A1GL,A2GL,A3GL,A4GL,A5GL,A6GL,A7GL,A8GL
      COMMON/ CCH2/ A1CH,A2CH,A3CH,A4CH,A5CH,A6CH,A7CH,A8CH
      COMMON/ CBO2/ A1BO,A2BO,A3BO,A4BO,A5BO,A6BO,A7BO,A8BO
      COMMON/expp/aa
      COMMON/ IFIT/ I
      COMMON/collider/ih1,ih2
             
      UV=A1UV(I)*x**A2UV(I)*(1-x)**A3UV(I)*(1+A4UV(I)*x+A5UV(I)*x**(0.5)
     .     +A6UV(I)*x**(1.5)+A7UV(I)*X**2 + A8UV(I)*X**(aa))
      DV=A1DV(I)*x**A2DV(I)*(1-x)**A3DV(I)*(1+A4DV(I)*x+A5DV(I)*x**(0.5)
     .     +A6DV(I)*x**(1.5)+A7DV(I)*X**2 + A8DV(I)*X**(aa))

      US=A1US(I)*x**A2US(I)*(1-x)**A3US(I)*(1+A4US(I)*x+A5US(I)*x**(0.5)
     .     +A6US(I)*x**(1.5)+A7US(I)*X**2 + A8US(I)*X**(aa))
      DS=A1DS(I)*x**A2DS(I)*(1-x)**A3DS(I)*(1+A4DS(I)*x+A5DS(I)*x**(0.5)
     .     +A6DS(I)*x**(1.5)+A7DS(I)*X**2 + A8DS(I)*X**(aa))

       U=UV+US
       D=DV+DS

      SS=A1SS(I)*x**A2SS(I)*(1-x)**A3SS(I)*(1+A4SS(I)*x+A5SS(I)*x**(0.5)
     .     +A6SS(I)*x**(1.5)+A7SS(I)*X**2 + A8SS(I)*X**(aa))

      GL=A1GL(I)*x**A2GL(I)*(1-x)**A3GL(I)*(1+A4GL(I)*x+A5GL(I)*x**(0.5)
     .     +A6GL(I)*x**(1.5)+A7GL(I)*X**2 + A8GL(I)*X**(aa))

      CH=A1CH(I)*x**A2CH(I)*(1-x)**A3CH(I)*(1+A4CH(I)*x+A5CH(I)*x**(0.5)
     .     +A6CH(I)*x**(1.5)+A7CH(I)*X**2 + A8CH(I)*X**(aa))

      BO=A1BO(I)*x**A2BO(I)*(1-x)**A3BO(I)*(1+A4BO(I)*x+A5BO(I)*x**(0.5)
     .     +A6BO(I)*x**(1.5)+A7BO(I)*X**2 + A8BO(I)*X**(aa))

       if (ih2.eq.-1) then
c pbar  u<->ubar  d<->dbar       
       UTEMP=U 
       U=US
       US=UTEMP

       DTEMP=D
       D=DS
       DS=DTEMP         
       endif
       RETURN
       END



      function etaintegrate(etalim) 
C Integrates over full rapidity range using 20 points gaussian quadratures
       IMPLICIT NONE
       Double precision xxx20(20),www20(20),etalim,eta,xm,xc,x,xnorm,
     .                    xsection,ax,xborn,xborn2,xnormal,resu,
     .                    xsection2,resummed,ss,etaintegrate
       integer jjj,isetproton,ih1,ih2
       include 'scales.h' 
!       COMMON/isetproton/isetproton
       COMMON/collider/ih1,ih2
       COMMON / ETA / ETA
       data xxx20/-0.993128599185d0,-0.963971927278d0,-0.912234428251d0,
     .  -0.839116971822d0,-0.74633190646d0,-0.636053680727d0,
     .  -0.510867001951d0,-0.373706088715d0,-0.227785851142d0,
     .   -0.0765265211335d0,0.0765265211335d0,0.227785851142d0,
     .   0.373706088715d0,0.510867001951d0,0.636053680727d0,
     .   0.74633190646d0,0.839116971822d0,0.912234428251d0,
     .   0.963971927278d0,0.993128599185d0/

      data www20/0.0176140070678d0,0.0406014298819d0,0.0626720482976d0,
     .   0.0832767415506d0,0.101930119826d0,0.118194531969d0,
     .   0.131688638458d0,0.142096109327d0,0.149172986482d0,
     .   0.15275338714d0,0.15275338714d0,0.149172986482d0,
     .   0.142096109327d0,0.131688638458d0,0.118194531969d0,
     .   0.101930119826d0,0.0832767415506d0,0.0626720482976d0,
     .   0.0406014298819d0,0.0176140070678d0/

       ss=0d0
c integral between 1 and etalim  20 points
       xm=0.5d0*(etalim+etalim)
       xc=0.5d0*(etalim-etalim)   
        do  jjj=1,20
       eta=xc+xm*xxx20(jjj)
C Evaluate kinematics and logs       
       call evaluatekin(eta,qt,q2)
C Normalization
      x=(Q2/Shad)
      xnorm=x*qt/2/Q2
C Normalize PDF's to born
       ax=dlog(x)
!       xborn=xsection(ax,eta,muf2,isetproton,ih1,ih2)
!       xborn2=xsection2(ax,eta,muf2,isetproton,ih1,ih2)
       muf=dsqrt(muf2)
       xborn=xsection(ax,eta,muf,ih1,ih2)
       xborn2=xsection2(ax,eta,muf,ih1,ih2)
       xnormal=xborn/xborn2
       if ((xnormal.gt.1.15).or.(xnormal.lt.0.85)) xnormal=1d0
C COMPUTE RESUMMED CONTRIBUTION (b integral 0->infinity)
      call resummation(resu)
      resummed=resu*xnorm *xnormal
      ss=ss+www20(jjj)*(resummed)*xm
          
c       write(6,*)qt,eta,resummed,xnormal
c       if ((resummed.lt.-0.2d0).or.(resummed.gt.15)) then
c        write(6,*)'problem at qt,eta, resu=',qt,eta,resummed
c        write(12,*)'problem at qt,eta, resu=',qt,eta,resummed
c        endif
      enddo 
       etaintegrate=ss
       return
       end
