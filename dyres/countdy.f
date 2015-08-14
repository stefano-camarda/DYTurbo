CC    Counterterm to be subtracted from real+virt to get a finite
CC    cross section at qt->0

C     Version that allows to separate also qg channel

C     Scale dependence included up to NNLO

      double precision function countterm(costh,mm,qtt,yy,mode)
      implicit none
      double precision costh,mm,qtt,yy
      integer mode
      double precision cthmom0,cthmom1,cthmom2
      include 'constants.f'
      include 'realonly.f'
      include 'virtonly.f'
      include 'noglue.f'
      include 'vegas_common.f'
      include 'ptilde.f'
      include 'npart.f'
      include 'scale.f'
      include 'facscale.f'
      include 'zerowidth.f'
      include 'efficiency.f'
      include 'masses.f'
      include 'limits.f'
C
      include 'jetlabel.f'
      include 'qcdcouple.f'
      include 'phasemin.f'
      include 'rescoeff.f'
      include 'dynamicscale.f'
C
      integer ih1,ih2,j,k,l,nd,nmax,nmin,nvec,order
      integer nproc
      common/nproc/nproc
      double precision vector(mxdim),W,val,xint
      double precision sqrts,qtmax
      double precision p(mxpart,4),pjet(mxpart,4),p1ext(4),p2ext(4)
      double precision pswt,rscalestart,fscalestart
      double precision s(mxpart,mxpart),wgt,msq(-nf:nf,-nf:nf)
      double precision msqc(-nf:nf,-nf:nf),xmsq(0:maxd)
      double precision flux,BrnRat,xreal,xreal2
      double precision xx1,xx2,q(mxpart,4)
      double precision m3,m4,m5,qtcut,xqtcut,switch,qt,m
CC
      logical cuts
      double precision x1,x2,dot,ptrans(mxpart,4)
      double precision q2,qt2,shat,Itilde
      double precision fx10(-nf:nf),fx20(-nf:nf)
      double precision fx1p(-nf:nf),fx2p(-nf:nf)
      double precision alfa,beta,diff,Pqq,Pqg,Pqqint,Cqq,Cqg
      double precision xjacq2,xjacqt2,xth,x3,almin,almax
      double precision xmio,fluxborn,pswt0,xmioOLD
      double precision shad,yq,zmax,tauh,Vol,y3
      double precision xx0(2),xx10,xx20,y34,cosh2y34
      double precision logxx10,logxx20
      double precision sig1,sig2,LR,LF,LQ
      double precision sig11,sig12
      double precision sig21,sig22,sig23,sig24
      double precision tdelta,tH1st,tH1stF,tH1stQ,tgaga,tcga,tgamma2
      double precision LL1,LL2,LL3,LL4
      double precision z1,z2,diff1,diff2,cut
      double precision D0int,D1int
      double precision Pqqqq,Pqqqg,Pqggq,Pqggg
      double precision CqqPqq,CqqPqg,CqgPgq,CqgPgg
      double precision P2qg,P2qqV,P2qqbV,P2qqS
      double precision diffg10,diffg20,diffc10,diffc20
      double precision diffg1f,diffg2f,diffc1f,diffc2f
      external Itilde,Pqq,Pqg,Cqq,Cqg,Pqqint,D0int,D1int
      external Pqqqq,Pqqqg,Pqggq,Pqggg,CqqPqq,CqqPqg,CqgPgq,CqgPgg
      external P2qqV,P2qqbV,P2qg,P2qqS

      common/xmio/xmio
      common/xx0/xx0
      common/qtcut/xqtcut
      common/nnlo/order 

CC
CC    Variables passed from virtint or lowint
CC
      common/count/qt2,q2,shat
CC
       COMMON/a_param/a_param,b0p
       double precision a_param,b0p
C
      integer n2,n3,sgnj,sgnk,flgq
      double precision mass2,width2,mass3,width3
      common/breit/n2,n3,mass2,width2,mass3,width3
      common/xreal/xreal,xreal2
      logical bin,first,failed
      logical incldip(0:maxd),includedipole,includereal
      logical creatent,dswhisto
      common/density/ih1,ih2
      common/energy/sqrts
      common/bin/bin
      common/Pext/p1ext,p2ext
      common/nmax/nmax
      common/BrnRat/BrnRat
      common/nmin/nmin
      common/incldip/incldip
      common/outputflags/creatent,dswhisto
      data p/48*0d0/
      data first/.true./
      save first,rscalestart,fscalestart

      real *8 pV(4),p4cm(4),ptemp(4)
      double precision phi, phi_lep, mt2

      COMMON/SIGMAIJ/SIGMAIJ
      double precision sigmaij(-5:5,-5:5)

      double precision kt1,kt2,qP1,qP2,zeta1

      integer ai,aj,bi,bj
      integer alfaintervals,betaintervals
      double precision alfaa,alfab,alfac,alfam
      double precision betaa,betab,betac,betam
      include 'gauss.inc'

      
      countterm=0d0 
      do nd=0,1
         xmsq(nd)=0d0
      enddo     

      m=mm
      q2=mm*mm

      qt=qtt
      qt2=qtt*qtt
      
C     C   Set qtcut limit
      qtcut=xqtcut*dsqrt(q2)
      if (qt.lt.qtcut)  goto 999

C     C   Set qtmax (kinematical limit)
      if(qt2.gt.((sqrts**2+q2)**2/(4d0*sqrts**2)-q2)) goto 999

!     SWITCHING FUNCTIONS
      switch=1d0
      if(qt.ge.m*3/4d0)  switch=dexp(-(m*3/4d0-qt)**2/(m/2d0)**2) ! GAUSS SWITCH
!!!!!!!!!!!!
!      if(qt.ge.m)  switch=dexp(-(m-qt)**2/(m/2d0)**2)           ! GAUSS SWITCH
!      if(qt.ge.m/2d0)  switch=dexp(-(m/2d0-qt)**2/(m/2d0)**2)              ! GAUSS SWITCH
!      if(qt.ge.m/2d0)        switch=dexp(-2d0*(m/2d0-qt)**2/(m/2d0)**2)    ! GAUSS SWITCH FASTER
!      if(qt.ge.m/2d0)        switch=dexp((m2/4d0-qt2)/m2)                  ! EXP SWITCH
!      if(qt.ge.m/2d0)        switch=dexp((m2/4d0-qt2)/m2*4d0)              ! EXP SWITCH FASTER
!      if(qt.ge.m/2d0)        switch=(dcos(pi/50d0*(qt-45.6d0))+1d0)/2d0    ! COS SWITCH
!      if(qt.ge.95.6)         switch=0d0                                    ! COS SWITCH

      if(switch.le.0.01d0)    goto 999    ! 

      shad=sqrts**2

c used in besselkfast
      xmio=dsqrt(qt2/(q2/a_param**2))
      
      Vol=1d0

CC   Dynamic scale
      if(dynamicscale) call scaleset(q2)

CC   LR,LF,LQ

      LR=dlog(q2/scale**2)
      LF=dlog(q2/facscale**2) 
      LQ=2d0*dlog(a_param)


CC   LL1,LL2,LL3,LL4: large log (squared) corresponding to eq. (136) 
CC   In this way normalization is fixed to dsigma/dqt2

      LL1=Itilde(1)/q2**2*a_param**2
      LL2=Itilde(2)/q2**2*a_param**2
      LL3=Itilde(3)/q2**2*a_param**2
      LL4=Itilde(4)/q2**2*a_param**2

c     *****************************************
c     Rewritten phase space generation
      mt=dsqrt(q2+qt2)
      xx0(1)=dsqrt(q2/sqrts**2)*dexp(+yy)
      xx0(2)=dsqrt(q2/sqrts**2)*dexp(-yy)

c---check if x is out of normal range
      if   ((xx0(1) .gt. 1d0) 
     & .or. (xx0(2) .gt. 1d0)
     & .or. (xx0(1) .lt. xmin)
     & .or. (xx0(2) .lt. xmin)) return

c     generate p3 and p4 4-momenta
      phi = twopi*0.25d0
      phi_lep = 0.5d0*twopi
       mt2=q2+qt2

       pV(1)=qt*dcos(phi)
       pV(2)=qt*dsin(phi)
       pV(3)=0.5d0*dsqrt(mt2)*(dexp(yy)-dexp(-yy))
       pV(4)=0.5d0*dsqrt(mt2)*(dexp(yy)+dexp(-yy))

c     CS FRAME PRESCRIPTION
c      kt1=pV(1)/2d0
c      kt2=pV(2)/2d0

c     naive prescription  
       kt1=(1d0+pV(3)/(sqrt(q2)+pV(4)))*pV(1)/2d0;
       kt2=(1d0+pV(3)/(sqrt(q2)+pV(4)))*pV(2)/2d0;

       
       zeta1=1d0/q2/2d0*(q2+2d0*(pV(1)*kt1+pV(2)*kt2)
     &      +dsqrt((q2+2d0*(pV(1)*kt1+pV(2)*kt2))**2
     &      -4d0*mt**2*(kt1**2+kt2**2)))

       qP1=(pV(4)-pV(3))*sqrts/2d0
       qP2=(pV(4)+pV(3))*sqrts/2d0

! incoming quarks
       p(1,4)=sqrts/2d0*(zeta1*q2/2d0/qP1
     &      +(kt1**2+kt2**2)/zeta1*qP1/q2/sqrts**2*2d0)
       p(1,1)=kt1
       p(1,2)=kt2
       p(1,3)=sqrts/2d0*(zeta1*q2/2d0/qP1
     &      -(kt1**2+kt2**2)/zeta1*qP1/q2/sqrts**2*2d0)

       p(2,4)=pV(4)-p(1,4)
       p(2,1)=pV(1)-p(1,1)
       p(2,2)=pV(2)-p(1,2)
       p(2,3)=pV(3)-p(1,3)

       p4cm(4)=m/2d0
       p4cm(1)=p4cm(4)*dsin(dacos(costh))*dsin(phi_lep)
       p4cm(2)=p4cm(4)*dsin(dacos(costh))*dcos(phi_lep)
       p4cm(3)=p4cm(4)*costh

!     Boost to go in the Lab frame
       call boost(m,pV,p4cm,ptemp)
       p(3,1)=ptemp(1)
       p(3,2)=ptemp(2)
       p(3,3)=ptemp(3)
       p(3,4)=ptemp(4)

!  momentum of the second lepton
       p(4,4)=pV(4)-p(3,4)
       p(4,1)=pV(1)-p(3,1)
       p(4,2)=pV(2)-p(3,2)
       p(4,3)=pV(3)-p(3,3)

c      print*,'phase space in counterterm'
c      print*,p(1,1),p(1,2),p(1,3),p(1,4)
c      print*,p(2,1),p(2,2),p(2,3),p(2,4)
c      print*,p(3,1),p(3,2),p(3,3),p(3,4)
c      print*,p(4,1),p(4,2),p(4,3),p(4,4)
c      print*,'mass',sqrt((p(4,4)+p(3,4))**2
c     +     - (p(4,1)+p(3,1))**2
c     +     - (p(4,2)+p(3,2))**2
c     +     - (p(4,3)+p(3,3))**2)
c      print*,'y',0.5d0*log((p(4,4)+p(3,4) + (p(4,3)+p(3,3)))/
c     +     (p(4,4)+p(3,4) - (p(4,3)+p(3,3))))
c      print*,'pt',sqrt((p(4,1)+p(3,1))**2 + (p(4,2)+p(3,2))**2)
c      print*
c     End of phase space generation
c     *****************************************

           
CC    Generate event to be binned
      y34=yy
      cosh2y34=((dexp(y34)+dexp(-y34))*0.5d0)**2

CC   Set qtmax (kinematical limit)
      if(qt2.gt.((sqrts**2+q2)**2/(4d0*sqrts**2*(cosh2y34)-q2)))goto 999


      
C     C Compute Born matrix element
      if(nproc.eq.3)then
         call qqb_z(p,msqc)
      else
         call qqb_w(p,msqc)
      endif

c **************** Check matrix element calculation (works only with naive prescription)
c      call initsigma(m,-costh)
c      print *,1,-1,sigmaij(2,-2)
c      print *,-1,1,sigmaij(-2,2)
c      print *,2,-2,sigmaij(1,-1)
c      print *,-2,2,sigmaij(-1,1)
c
c      print *,1,-1,(msqc(1,-1))/(sigmaij(2,-2)*pi*6d0/fbGeV2*(2*q2))
c      print *,-1,1,(msqc(-1,1))/(sigmaij(-2,2)*pi*6d0/fbGeV2*(2*q2))
c      print *,2,-2,(msqc(2,-2))/(sigmaij(1,-1)*pi*6d0/fbGeV2*(2*q2))
c      print *,-2,2,(msqc(-2,2))/(sigmaij(-1,1)*pi*6d0/fbGeV2*(2*q2))
c      print *,3,-3,(msqc(3,-3))/(sigmaij(3,-3)*pi*6d0/fbGeV2*(2*q2))
c      print *,-3,3,(msqc(-3,3))/(sigmaij(-3,3)*pi*6d0/fbGeV2*(2*q2))
c      print *,4,-4,(msqc(4,-4))/(sigmaij(4,-4)*pi*6d0/fbGeV2*(2*q2))
c      print *,-4,4,(msqc(-4,4))/(sigmaij(-4,4)*pi*6d0/fbGeV2*(2*q2))
c      print *,5,-5,(msqc(5,-5))/(sigmaij(5,-5)*pi*6d0/fbGeV2*(2*q2))
c      print *,-5,5,(msqc(-5,5))/(sigmaij(-5,5)*pi*6d0/fbGeV2*(2*q2))

      if (mode.eq.0) then
         call initsigma(mm,-costh)
      elseif (mode.eq.1) then
         call cthmoments(cthmom0,cthmom1,cthmom2)
         cthmom0=cthmom0*twopi*twopi
         cthmom1=cthmom1*twopi*twopi
         cthmom2=cthmom2*twopi*twopi
         call initsigmacth(mm,cthmom0,cthmom1,cthmom2)
      endif

      do j=-nf,nf
         do k=-nf,nf
            msqc(j,k) = 0d0
         enddo
      enddo
      if (nproc.eq.3) then
         msqc(1,-1)=sigmaij(2,-2)*pi*6d0/fbGeV2*(2*q2)
         msqc(-1,1)=sigmaij(-2,2)*pi*6d0/fbGeV2*(2*q2)
         msqc(2,-2)=sigmaij(1,-1)*pi*6d0/fbGeV2*(2*q2)
         msqc(-2,2)=sigmaij(-1,1)*pi*6d0/fbGeV2*(2*q2)
         msqc(3,-3)=sigmaij(3,-3)*pi*6d0/fbGeV2*(2*q2)
         msqc(-3,3)=sigmaij(-3,3)*pi*6d0/fbGeV2*(2*q2)
         msqc(4,-4)=sigmaij(4,-4)*pi*6d0/fbGeV2*(2*q2)
         msqc(-4,4)=sigmaij(-4,4)*pi*6d0/fbGeV2*(2*q2)
         msqc(5,-5)=sigmaij(5,-5)*pi*6d0/fbGeV2*(2*q2)
         msqc(-5,5)=sigmaij(-5,5)*pi*6d0/fbGeV2*(2*q2)
      endif
CC  

c---  calculate PDF's  
      xx10=xx0(1)
      xx20=xx0(2)
      logxx10=dlog(xx10)
      logxx20=dlog(xx20)
      
!      if(xx10.lt.1d-5)write(*,*) "q2,xx1",q2,xx10
!      if(xx20.lt.1d-5)write(*,*) "q2,xx1",q2,xx20!

      call fdist(ih1,xx10,facscale,fx10)
      call fdist(ih2,xx20,facscale,fx20)
      
C Scaled momentum fractions
      cut=1d-7

c     start the alfa beta integration
      alfaintervals = 20
      betaintervals = 20
      do ai=1,alfaintervals
         alfaa = 0d0+(1d0-0d0)*(ai-1)/alfaintervals
         alfab = 0d0+(1d0-0d0)*ai/alfaintervals
         alfac=0.5d0*(alfaa+alfab)
         alfam=0.5d0*(alfab-alfaa)
         do aj=1,4
            alfa=alfac+alfam*xxx4(aj)
            alfa=cut+(1-cut)*alfa
      z2=xx20**alfa
      call fdist(ih2,xx20**(1-alfa),facscale,fx2p)



      do bi=1,betaintervals
         betaa = 0d0+(1d0-0d0)*(bi-1)/betaintervals
         betab = 0d0+(1d0-0d0)*bi/betaintervals
         betac=0.5d0*(betaa+betab)
         betam=0.5d0*(betab-betaa)
         do bj=1,4
            beta=betac+betam*xxx4(bj)
            beta=cut+(1-cut)*beta
            
      z1=xx10**beta
      call fdist(ih1,xx10**(1-beta),facscale,fx1p)

CC Switch off gluon !!

      if(noglue) then
        fx10(0)=0d0
        fx20(0)=0d0
        fx1p(0)=0d0
        fx2p(0)=0d0
      endif

CC Gluon only !

      if(ggonly) then
       do j=1,5
       fx10(j)=0d0
       fx10(-j)=0d0
       fx1p(j)=0d0
       fx1p(-j)=0d0
       fx20(j)=0d0
       fx20(-j)=0d0
       fx2p(j)=0d0
       fx2p(-j)=0d0
       enddo
      endif

       flgq=1
       if(gqonly)flgq=0




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CC Start construction of the counterterm

        tdelta=0d0
        tH1st=0d0
        tH1stF=0d0
        tH1stQ=0d0
        tgaga=0d0
        tcga=0d0
        tgamma2=0d0

        diffc10=0d0
        diffc1f=0d0
        diffc20=0d0
        diffc2f=0d0

        diffg10=0d0
        diffg1f=0d0
        diffg20=0d0
        diffg2f=0d0

        sig1=0d0
        sig2=0d0

        sig11=0d0
        sig12=0d0
        sig21=0d0      
        sig22=0d0
        sig23=0d0
        sig24=0d0

      
      do j=-nf,nf
      do k=-nf,nf

      if(msqc(j,k).eq.0d0) goto 75


C     Simplest term without convolutions
  
      tdelta=tdelta+fx10(j)*fx20(k)*msqc(j,k)*flgq

C     Start H1st: to be used later

C     H1st delta term

      tH1st=tH1st+2*C1qqdelta*fx10(j)*fx20(k)*msqc(j,k)*flgq

!add resummation scale dependence

      tH1stQ=tH1stQ-(B1q+A1q/2d0*LQ)*LQ*fx10(j)*fx20(k)*msqc(j,k)*flgq

C     H1st: non delta terms, first leg


      tH1st=tH1st+(fx1p(j)*Cqq(z1)*flgq+fx1p(0)*Cqg(z1))
     & *(-logxx10)*fx20(k)*msqc(j,k)


C     H1st: non delta terms, second leg


      tH1st=tH1st+(fx2p(k)*Cqq(z2)*flgq+fx2p(0)*Cqg(z2))         
     & *(-logxx20)*fx10(j)*msqc(j,k)
      

C     H1st: muf dependence (LF factor to be added at the end)


c     gammaqq and gammaqg: first leg      


      diff=-logxx10
     &  *((fx1p(j)-fx10(j)*xx10**beta)*Pqq(z1)*flgq+fx1p(0)*Pqg(z1))
      tH1stF=tH1stF+diff*fx20(k)*msqc(j,k)
      tH1stF=tH1stF-Pqqint(xx10)*fx10(j)*fx20(k)*msqc(j,k)*flgq

c     gammaqq and gammaqg: second leg   


      diff=-logxx20
     &  *((fx2p(k)-fx20(k)*xx20**alfa)*Pqq(z2)*flgq+fx2p(0)*Pqg(z2))
      tH1stF=tH1stF+diff*fx10(j)*msqc(j,k)
      tH1stF=tH1stF-Pqqint(xx20)*fx10(j)*fx20(k)*msqc(j,k)*flgq

CC    End of H1st

      if(order.eq.1) goto 75

CC    Now (gamma+gamma)*(gamma+gamma) term: to be used later

C     First part: one gamma for each leg: FLGQ here is non trivial ! DONE


      diffg1f=-logxx10*(fx1p(j)-fx10(j)*xx10**beta)*Pqq(z1)
     &  - Pqqint(xx10)*fx10(j)


      diffg10=-logxx10*fx1p(0)*Pqg(z1)

      diffg2f=-logxx20*(fx2p(k)-fx20(k)*xx20**alfa)*Pqq(z2)
     &  - Pqqint(xx20)*fx20(k)


      diffg20=-logxx20*fx2p(0)*Pqg(z2)


      tgaga=tgaga+2*
     #   (flgq*diffg10*diffg20+flgq*diffg1f*diffg2f
     #   +diffg10*diffg2f+diffg1f*diffg20)*msqc(j,k)


CC     Second part: gamma*gamma terms

c     Pij * Pjk = D1ijjk (log(1-z)/(1-z))_+ + D0ijjk/(1-z)_+ 
c              + Pijjk(z) + Deltaijjk delta(1-z)

C     First leg

      
      diff1=-logxx10*(flgq*(fx1p(j)-fx10(j)*xx10**beta)
     &    *(D0qqqq/(1-z1)+D1qqqq*dlog(1-z1)/(1-z1))
     &    +fx1p(j)*Pqqqq(z1)*flgq+fx1p(0)*(Pqqqg(z1)+Pqggg(z1)))
     &    +(Deltaqqqq-D0qqqq*D0int(xx10)-D1qqqq*D1int(xx10))
     &    *fx10(j)*flgq


C    Second leg

      
      diff2=-logxx20*(flgq*(fx2p(k)-fx20(k)*xx20**alfa)
     &    *(D0qqqq/(1-z2)+D1qqqq*dlog(1-z2)/(1-z2))
     &    +fx2p(k)*Pqqqq(z2)*flgq+fx2p(0)*(Pqqqg(z2)+Pqggg(z2)))
     &    +(Deltaqqqq-D0qqqq*D0int(xx20)-D1qqqq*D1int(xx20))
     &    *fx20(k)*flgq


C     Include Pqggq

      do l=1,nf
      diff1=diff1-logxx10*(fx1p(l)+fx1p(-l))*Pqggq(z1)*flgq
      diff2=diff2-logxx20*(fx2p(l)+fx2p(-l))*Pqggq(z2)*flgq
      enddo

      tgaga=tgaga+diff1*fx20(k)*msqc(j,k)
      tgaga=tgaga+diff2*fx10(j)*msqc(j,k)



C    End of (gamma+gamma)*(gamma+gamma) term: FLGQ non trivial here ! DONE

C    Start  (C+C)*(gamma+gamma) term

c    gamma first leg, C second leg


      diffc2f=-logxx20*fx2p(k)*Cqq(z2)+C1qqdelta*fx20(k)

      diffc20=-logxx20*fx2p(0)*Cqg(z2)


      tcga=tcga+msqc(j,k)*
     # (flgq*diffg10*diffc20+flgq*diffg1f*diffc2f
     #          +diffg10*diffc2f+diffg1f*diffc20)


c    C first leg, gamma second leg

      diffc1f=-logxx10*fx1p(j)*Cqq(z1)+C1qqdelta*fx10(j)

      diffc10=-logxx10*fx1p(0)*Cqg(z1)

      tcga=tcga+msqc(j,k)*
     # (flgq*diffc10*diffg20+flgq*diffc1f*diffg2f
     #          +diffc10*diffg2f+diffc1f*diffg20)
    

c    C*gamma: first leg (ignore delta term in Cqq: taken into account with tH1stF)

      tcga=tcga
     &     +(fx1p(j)*CqqPqq(z1)*flgq+fx1p(0)*(CqqPqg(z1)+CqgPgg(z1)))
     &     *(-logxx10)*fx20(k)*msqc(j,k) 

c    C*gamma: second leg (ignore delta term in Cqq: taken into account with tH1stF)

      tcga=tcga
     &     +(fx2p(k)*CqqPqq(z2)*flgq+fx2p(0)*(CqqPqg(z2)+CqgPgg(z2)))
     &     *(-logxx20)*fx10(j)*msqc(j,k) 

c    Add Cqg*Pgq contribution

      do l=1,nf
      tcga=tcga+(fx1p(l)+fx1p(-l))*CqgPgq(z1)
     &           *(-logxx10)*fx20(k)*msqc(j,k)*flgq 
      tcga=tcga+(fx2p(l)+fx2p(-l))*CqgPgq(z2)
     &           *(-logxx20)*fx10(j)*msqc(j,k)*flgq 
      enddo

CC  Start 2-loop AP

C   Gluon + pure singlet


      do l=-nf,nf
      if(l.eq.0) then
      tgamma2=tgamma2+fx1p(0)*P2qg(z1)
     & *(-logxx10)*fx20(k)*msqc(j,k)
      tgamma2=tgamma2+fx2p(0)*P2qg(z2)
     & *(-logxx20)*fx10(j)*msqc(j,k)
      else
      tgamma2=tgamma2+fx1p(l)*P2qqS(z1)
     & *(-logxx10)*fx20(k)*msqc(j,k)*flgq
      tgamma2=tgamma2+fx2p(l)*P2qqS(z2)
     & *(-logxx20)*fx10(j)*msqc(j,k)*flgq
      endif
      enddo


C   P2qq non-singlet: regular part

      tgamma2=tgamma2+fx1p(j)*P2qqV(z1)
     & *(-logxx10)*fx20(k)*msqc(j,k)*flgq
      tgamma2=tgamma2+fx2p(k)*P2qqV(z2)
     & *(-logxx20)*fx10(j)*msqc(j,k)*flgq


C   P2qq non-singlet: 1/(1-z)_+


      diff=-logxx10
     &  *(fx1p(j)-fx10(j)*xx10**beta)/(1-z1)
     &  - D0int(xx10)*fx10(j)      
  
      tgamma2=tgamma2+2d0/3*Kappa*diff*fx20(k)*msqc(j,k)*flgq


      diff=-logxx20
     &  *(fx2p(k)-fx20(k)*xx20**alfa)/(1-z2)
     &  - D0int(xx20)*fx20(k)      
  
      tgamma2=tgamma2+2d0/3*Kappa*diff*fx10(j)*msqc(j,k)*flgq

      

C   P2qqb non singlet

      tgamma2=tgamma2+fx1p(-j)*P2qqbV(z1)
     & *(-logxx10)*fx20(k)*msqc(j,k)*flgq

      tgamma2=tgamma2+fx2p(-k)*P2qqbV(z2)
     & *(-logxx20)*fx10(j)*msqc(j,k)*flgq

 75   continue

      enddo
      enddo

CM 7/11  Resummation scale dependence added 


CC   First order

      sig12=-0.5d0*A1q*tdelta
      sig11=-(B1q+A1q*LQ)*tdelta-tH1stF

CC   Second order

      sig24=(A1q)**2/8*tdelta
       
      sig23=-beta0*A1q/3*tdelta-0.5d0*A1q*sig11

      sig22=0.5d0*(beta0*A1q*(LR-LQ)-A2q)*tdelta!
     &     -0.5d0*A1q*(tH1st+tH1stQ+(LF-LQ)*tH1stF)!
     &     -0.5d0*(B1q+A1q*LQ-beta0)*sig11!
     &     +0.5d0*(B1q+A1q*LQ)*tH1stF!
     &     +0.5d0*tgaga!


      sig21=-beta0*(LR-LQ)*sig11
     &     -(B1q+A1q*LQ)*(tH1st+tH1stQ+(LF-LQ)*tH1stF)
     &     -(LF-LQ)*tgaga-(B2q+A2q*LQ)*tdelta+beta0*tH1st-tcga-tgamma2
     &     +(B1q+0.5d0*A1q*LQ)*LQ*tH1stF  


c     Include missing delta term from C*gamma (no factor 2 here !)

      sig21=sig21-C1qqdelta*tH1stF

C     Include missing term from contact term in 2 loop AP

      sig21=sig21-2*Delta2qq*tdelta

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


CC Include as/pi factors and sum O(as) and O(as^2) contributions

      sig1=sig12*LL2+sig11*LL1
      sig2=sig24*LL4+sig23*LL3+sig22*LL2+sig21*LL1


        sig1=sig1*ason2pi*2
        sig2=sig2*(ason2pi*2)**2

        if(order.eq.1)then
           xmsq(1)=xmsq(1)-sig1*www4(aj)*www4(bj)*alfam*betam
        else
           xmsq(1)=xmsq(1)-(sig1+sig2)*www4(aj)*www4(bj)*alfam*betam
        endif
c        print *, alfa,beta, (sig1+sig2)
      enddo
      enddo
c      print *,alfa, (sig1+sig2)
      enddo
      enddo
c     end of alfa beta integration
      
CC Include iacobians (do not include jacobians)

      xmsq(1)=xmsq(1)*q2/shad/Vol


      countterm=0d0
      xint=0d0


C Flux for Born cross section
       fluxborn=fbGeV2/(2*q2)
C Multiply by BORN phase space weight
        xmsq(1)=xmsq(1)*fluxborn/BrnRat*((one/4d0/pi)**3)

! SWITCHING
        xmsq(1)=xmsq(1)*switch


 77    continue


c---Add to total

        xint=xmsq(1)
        val=xmsq(1)*wgt
        
      countterm=xint
c      print *, qtt, mm, yy, countterm
     
      xreal=xreal+xint*wgt/dfloat(itmx)
      xreal2=xreal2+(xint*wgt)**2/dfloat(itmx)
      

      return

 999  countterm=0d0
      ntotzero=ntotzero+1
 
      return
      end

cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c
cCC qq splitting function (with asopi normalization)
c
c      function Pqq(z)
c      implicit none
c      real *8 Pqq,z
c      Pqq=2d0/3*(1+z**2)/(1-z)
c      return
c      end
c
cCC qg splitting function (with asopi normalization)
c
c      function Pqg(z)
c      implicit none
c      real *8 Pqg,z
c      Pqg=0.25d0*(1-2*z*(1-z))
c      return
c      end
c
cCC Non delta term in Cqq coefficient (with asopi normalization)
c
c      function Cqq(z)
c      implicit none
c      real *8 Cqq,z
c      Cqq=2d0/3*(1-z)
c      return
c      end
c
c
cCC Cqg coefficient (with asopi normalization)
c
c      function Cqg(z)
c      implicit none
c      real *8 Cqg,z
c      Cqg=0.5d0*z*(1-z)
c      return
c      end
c
c
cCC Integral of Pqq=1/2 CF (1+x^2)/(1-x) from 0 to z
c
c      function Pqqint(z)
c      implicit none
c      real *8 Pqqint,z
c      Pqqint=-2d0/3*(z+z**2/2+2*dlog(1-z))
c      return
c      end
c
cCC Integral of 1/(1-x) from 0 to z
c
c      function D0int(z)
c      implicit none
c      real *8 D0int,z
c      D0int=-dlog(1-z)
c      return
c      end
c
cCC Integral of log(1-x)/(1-x) from 0 to z
c
c      function D1int(z)
c      implicit none
c      real *8 D1int,z
c      D1int=-0.5d0*dlog(1-z)**2
c      return
c      end
c
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
cC
cC                P*P convolutions
cC
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
cCC Regular part of Pqq*Pqq (checked !)
c
c      function Pqqqq(z)
c      implicit none
c      real *8 Pqqqq,z
c      Pqqqq=4d0/9*(-4*dlog(z)/(1-z)-2*(1-z)
c     &  +(1+z)*(3*dlog(z)-4*dlog(1-z)-3))
c      return
c      end
c
c
cCC Pqq*Pqg (checked !)
c
c      function Pqqqg(z) 
c      implicit none
c      real *8 Pqqqg,z
c      Pqqqg=1d0/3*((z**2+(1-z)**2)*dlog((1-z)/z)
c     &  -(z-0.5d0)*dlog(z)+z-0.25d0)
c      return
c      end
c
cCC Pqg*Pgq (checked !)
c
c      function Pqggq(z)
c      implicit none
c      real *8 Pqggq,z
c      Pqggq=1d0/3*(2d0/3/z+(1+z)*dlog(z)-2d0/3*z**2-0.5d0*(z-1))
c      return
c      end
c
c
cCC Full Pqg*Pgg (checked !)
c
c      function Pqggg(z)
c      implicit none
c      real *8 Pqggg,z,beta0,Pqg
c      integer nf
c      external Pqg
c      nf=5
c      beta0=(33-2*nf)/12d0
c      Pqggg=1.5d0*(1/3d0/z+(z**2-z+0.5d0)*dlog(1-z)
c     &     +(2*z+0.5d0)*dlog(z)+0.25d0+2*z-31d0/12*z**2)
c
c      Pqggg=Pqggg+beta0*Pqg(z)
c      return
c      end
c
c
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
cC
cC                C*P convolutions
cC
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
cCC Cqq*Pqq (without delta term in Cqq) (checked !)
c
c      function CqqPqq(z)
c      implicit none
c      real *8 CqqPqq,z
c      CqqPqq=2d0/9*(1-z)*(4*dlog(1-z)-2*dlog(z)-1)
c      return
c      end
c
cCC Cqq*Pqg (without delta term in Cqq) (checked !)
c
c      function CqqPqg(z)
c      implicit none
c      real *8 CqqPqg,z
c      CqqPqg=(-2+z+z**2-(1+2*z)*dlog(z))/6d0
c      return
c      end
c
cCC Cqg*Pgq (checked !)
c
c      function CqgPgq(z) 
c      implicit none
c      real *8 CqgPgq,z
c      CqgPgq=(1d0/3/z-1+2*z**2/3-z*dlog(z))/3d0
c      return
c      end
c
cCC Cqg*Pgg (checked !)
c
c      function CqgPgg(z)
c      implicit none
c      real *8 CqgPgg,z,beta0
c      integer nf
c      nf=5
c      beta0=(33-2*nf)/12d0
c      CqgPgg=3d0/4*(2*z*(1-z)*dlog(1-z)-4*z*dlog(z)
c     &      +1d0/3/z-1-5*z+17d0*z**2/3)+beta0/2*z*(1-z)
c      return
c      end
c
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
cC
cC           Two loop AP:  pqq of ESW is my 3/2 Pqq
cC
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
cC     Pqq NS: Eq. (4.107) ESW (no 1/(1-x)_+ and delta term)
c
c      function P2qqV(x)
c      implicit none
c      real *8 x,P2qqV,Pqq,pi
c      integer nf
c      external Pqq
c
c      pi=3.14159265358979d0
c      nf=5
c
c      P2qqV=16d0/9*(-(2*dlog(x)*dlog(1-x)+1.5d0*dlog(x))*3d0/2*Pqq(x)
c     &     -(1.5d0+3.5d0*x)*dlog(x)-0.5d0*(1+x)*dlog(x)**2-5*(1-x))
c     &     +4*((0.5d0*dlog(x)**2+11d0/6*dlog(x))*3d0/2*Pqq(x)
c     &     -(67d0/18-pi**2/6)*(1+x)
c     &     +(1+x)*dlog(x)+20d0/3*(1-x))
c     &     +2d0/3d0*nf*(-dlog(x)*Pqq(x)+10d0/9*(1+x)-4d0/3*(1-x))
c
cc     Change to as/pi normalization
c
c      P2qqV=P2qqV/4  
c
c      return
c      end
c
c
cC    Pqqb NS: Eq. (4.108) ESW
c
c      function P2qqbV(x)
c      implicit none
c      real *8 x,P2qqbV,Pqq,S2
c      external Pqq,S2
c
c      P2qqbV=-2d0/9*(3d0*Pqq(-x)*S2(x)+2*(1+x)*dlog(x)+4*(1-x))
c      
cc     Change to as/pi normalization
c
c      P2qqbV=P2qqbV/4 
c
c      return
c      end
c
c
c
cC    Pqg Singlet: Eq. (4.110) ESW (ESW Pqg is 4 times my Pqg)
c
c      function P2qg(x)
c      implicit none
c      real *8 x,P2qg,Pqg,pi,S2,logx,logomxsx
c      external Pqg,S2
c
c      pi=3.14159265358979d0
c      logx=dlog(x)
c      logomxsx=dlog((1-x)/x)
c
c      P2qg=2d0/3*(4-9*x-(1-4*x)*logx-(1-2*x)*logx**2+4*dlog(1-x)
c     &    +(2*logomxsx**2-4*logomxsx-2d0/3*pi**2+10d0)*4*Pqg(x))
c     &    +1.5d0*(182d0/9+14d0/9*x+40d0/9/x+(136d0/3*x-38d0/3)*logx
c     &    -4*dlog(1-x)-(2+8*x)*logx**2+8*Pqg(-x)*S2(x)
c     &    +(-logx**2+44d0/3*logx-2*dlog(1-x)**2+4*dlog(1-x)+pi**2/3
c     &    -218d0/9)*4*Pqg(x))
c
cc     Change to as/pi normalization
c
c      P2qg=P2qg/4d0
c  
cc     Divide by 2 to eliminate 2nf factor
c
c      P2qg=P2qg/2d0
c
c      return
c      end
c
cC     Pqq Pure Singlet appearing in ESW Eq. (4.95)
cC     PSqq=PSqqb
cC     Obtained through Eq.(4.101)
cC     PSqq=1/2/nf (P2qq-P2qqbV-P2qqV) (contains only CF TR=2/3)
c
c      function P2qqS(x)
c      implicit none
c      real *8 P2qqS,x
c
c      P2qqS=2d0/3*(20 - 18*x + 54*x**2 - 56*x**3
c     &    +3*x*(3 + 15*x + 8*x**2)*dlog(x) 
c     &    - 9*x*(1 + x)*dlog(x)**2)/(9*x)
c      
c      P2qqS=P2qqS/4
c
c      return
c      end
c
c
cC    S2: Eq. (4.114) ESW
c
c      function S2(x)
c      implicit none
c      real *8 x,pi,S2,myli2
c      external myli2      
c      pi=3.14159265358979d0
c
c      S2=-2*myli2(-x)+0.5d0*dlog(x)**2-2*dlog(x)*dlog(1+x)-pi**2/6
c      return
c      end
