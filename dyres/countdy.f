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
      integer j,k,l,nd,order
      integer nproc
      common/nproc/nproc
      double precision W,xint
      double precision p(mxpart,4)
      double precision rscalestart,fscalestart
      double precision s(mxpart,mxpart)
      double precision msqc(-nf:nf,-nf:nf),xmsq(0:maxd)
      double precision flux,BrnRat
      double precision xx1,xx2,q(mxpart,4)
      double precision m3,m4,m5,qtcut,xqtcut,switch,qt,m
CC
      double precision x1,x2
      double precision q2,qt2,shat,Itilde
      double precision fx10(-nf:nf),fx20(-nf:nf)
      double precision fx1p(-nf:nf),fx2p(-nf:nf)
      double precision log1x,logz1,logz2
      double precision alfa,beta,diff,Pqq,Pqg,Pqqint,Cqq,Cqg
      double precision xmio,fluxborn,xmioOLD
      double precision shad,yq,zmax,tauh,Vol,y3
      double precision xx0(2),xx10,xx20,y34,cosh2y34
      double precision logxx10,logxx20
      double precision sig1,sig2,LR,LF,LQ
      double precision sig11,sig12
      double precision sig21,sig22,sig23,sig24
      double precision tdelta,tH1st,tH1stF,tH1stQ,tgaga,tcga,tgamma2
      double precision tdeltajk(-nf:nf,-nf:nf)
      double precision tH1stjk(-nf:nf,-nf:nf)
      double precision tH1stQjk(-nf:nf,-nf:nf)
      double precision tH1stFjk(-nf:nf,-nf:nf)
      double precision tgagajk(-nf:nf,-nf:nf)
      double precision tcgajk(-nf:nf,-nf:nf)
      double precision tgamma2jk(-nf:nf,-nf:nf)
      double precision LL1,LL2,LL3,LL4
      double precision z1,z2,diff1,diff2
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
      integer n2,n3,flgq
      logical bin,first
      integer ih1,ih2
      common/density/ih1,ih2
      double precision sqrts
      common/energy/sqrts
      common/bin/bin
      common/BrnRat/BrnRat
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

c common block form ctquadinit
      integer ctintervals,ctrule,ctdim
      parameter (ctintervals=20)
      parameter (ctrule=10)
      parameter (ctdim=ctrule*ctintervals)
      double precision ctx(ctdim)
      double precision ctw(ctdim)
      integer quadpoints
      common/ctweights/ctx,ctw,quadpoints

c     cached variables for fast integration
      integer ii,jj
      double precision cz1(ctdim),cz2(ctdim)
      double precision cfx1p(-nf:nf,ctdim)
      double precision cfx2p(-nf:nf,ctdim)
      double precision oz1(ctdim),oz2(ctdim)
      double precision Pqqint1,Pqqint2
      double precision D0intx1,D0intx2
      double precision D1intx1,D1intx2
      double precision log1z1(ctdim),log1z2(ctdim)
      double precision Cqqz1(ctdim),Cqqz2(ctdim)
      double precision Pqqz1(ctdim),Pqqz2(ctdim)
      double precision Pqggqz1(ctdim),Pqggqz2(ctdim)
      double precision Pqqqqz1(ctdim),Pqqqqz2(ctdim)
      double precision Pqqqgz1(ctdim),Pqqqgz2(ctdim)
      double precision Pqgggz1(ctdim),Pqgggz2(ctdim)
      double precision CqgPgqz1(ctdim),CqgPgqz2(ctdim)
      double precision CqqPqqz1(ctdim),CqqPqqz2(ctdim)
      double precision CqgPggz1(ctdim),CqgPggz2(ctdim)
      double precision CqqPqgz1(ctdim),CqqPqgz2(ctdim)
      double precision P2qqVz1(ctdim),P2qqVz2(ctdim)
      double precision P2qqbVz1(ctdim),P2qqbVz2(ctdim)
      double precision P2qqSz1(ctdim),P2qqSz2(ctdim)
      double precision P2qgz1(ctdim),P2qgz2(ctdim)

      double precision sumfx1p,sumfx2p
      
c     Input large logs from ctqtint
      double precision LL1jk(-nf:nf,-nf:nf),LL2jk(-nf:nf,-nf:nf)
      double precision LL3jk(-nf:nf,-nf:nf),LL4jk(-nf:nf,-nf:nf)
      common/largelogs/LL1jk,LL2jk,LL3jk,LL4jk


      if(first)  then           !! ONE TIME INITIALIZATION
         print *, 'first call to countterm'
         call ctquadinit
         first = .false.
      end if

      countterm=0d0 
      do nd=0,1
         xmsq(nd)=0d0
      enddo     

      m=mm
      q2=mm*mm

CC   Dynamic scale
      if(dynamicscale) call scaleset(q2)

CC   LR,LF,LQ

      LR=dlog(q2/scale**2)
      LF=dlog(q2/facscale**2) 
      LQ=2d0*dlog(a_param)

      if (mode.eq.0.or.mode.eq.1) then
         qt=qtt
         qt2=qtt*qtt
!     SWITCHING FUNCTIONS
         switch=1d0
         if(qt.ge.m*3/4d0)  switch=dexp(-(m*3/4d0-qt)**2/(m/2d0)**2) ! GAUSS SWITCH
!!!!!!!!!!!!
!     if(qt.ge.m)  switch=dexp(-(m-qt)**2/(m/2d0)**2)           ! GAUSS SWITCH
!     if(qt.ge.m/2d0)  switch=dexp(-(m/2d0-qt)**2/(m/2d0)**2)              ! GAUSS SWITCH
!     if(qt.ge.m/2d0)        switch=dexp(-2d0*(m/2d0-qt)**2/(m/2d0)**2)    ! GAUSS SWITCH FASTER
!     if(qt.ge.m/2d0)        switch=dexp((m2/4d0-qt2)/m2)                  ! EXP SWITCH
!     if(qt.ge.m/2d0)        switch=dexp((m2/4d0-qt2)/m2*4d0)              ! EXP SWITCH FASTER
!     if(qt.ge.m/2d0)        switch=(dcos(pi/50d0*(qt-45.6d0))+1d0)/2d0    ! COS SWITCH
!     if(qt.ge.95.6)         switch=0d0                                    ! COS SWITCH

         if(switch.le.0.01d0) return

c     used in besselkfast for Itilde
         xmio=dsqrt(qt2/(q2/a_param**2))
         
C     C   LL1,LL2,LL3,LL4: large log (squared) corresponding to eq. (136) 
C     C   In this way normalization is fixed to dsigma/dqt2

c     depends on xmio, which depends on qt2 and q2      
         LL1=Itilde(1)/q2**2*a_param**2
         LL2=Itilde(2)/q2**2*a_param**2
         LL3=Itilde(3)/q2**2*a_param**2
         LL4=Itilde(4)/q2**2*a_param**2
      elseif (mode.eq.2) then
c qt integration already performed
         switch=1d0
      endif

cc **************** Check matrix element calculation (works only with naive prescription)
cc     Rewritten phase space generation
cc     *****************************************
cc     generate p3 and p4 4-momenta
c      mt=dsqrt(q2+qt2)
c     
cc     pick whatever value of phi and phi_lep
c      phi = twopi*0.25d0
c      phi_lep = 0.5d0*twopi
c      mt2=q2+qt2
c
c      pV(1)=qt*dcos(phi)
c      pV(2)=qt*dsin(phi)
c      pV(3)=0.5d0*dsqrt(mt2)*(dexp(yy)-dexp(-yy))
c      pV(4)=0.5d0*dsqrt(mt2)*(dexp(yy)+dexp(-yy))
c
cc     CS FRAME PRESCRIPTION
cc      kt1=pV(1)/2d0
cc      kt2=pV(2)/2d0
c
cc     naive prescription  
c      kt1=(1d0+pV(3)/(sqrt(q2)+pV(4)))*pV(1)/2d0;
c      kt2=(1d0+pV(3)/(sqrt(q2)+pV(4)))*pV(2)/2d0;
c
c      zeta1=1d0/q2/2d0*(q2+2d0*(pV(1)*kt1+pV(2)*kt2)
c     &     +dsqrt((q2+2d0*(pV(1)*kt1+pV(2)*kt2))**2
c     &     -4d0*mt**2*(kt1**2+kt2**2)))
c      
c      qP1=(pV(4)-pV(3))*sqrts/2d0
c      qP2=(pV(4)+pV(3))*sqrts/2d0
c
c! incoming quarks
c      p(1,4)=sqrts/2d0*(zeta1*q2/2d0/qP1
c     &     +(kt1**2+kt2**2)/zeta1*qP1/q2/sqrts**2*2d0)
c      p(1,1)=kt1
c      p(1,2)=kt2
c      p(1,3)=sqrts/2d0*(zeta1*q2/2d0/qP1
c     &     -(kt1**2+kt2**2)/zeta1*qP1/q2/sqrts**2*2d0)
c      
c      p(2,4)=pV(4)-p(1,4)
c      p(2,1)=pV(1)-p(1,1)
c      p(2,2)=pV(2)-p(1,2)
c      p(2,3)=pV(3)-p(1,3)
c
c      p4cm(4)=m/2d0
c      p4cm(1)=p4cm(4)*dsin(dacos(costh))*dsin(phi_lep)
c      p4cm(2)=p4cm(4)*dsin(dacos(costh))*dcos(phi_lep)
c      p4cm(3)=p4cm(4)*costh
c
c!     Boost to go in the Lab frame
c      call boost(m,pV,p4cm,ptemp)
c      p(3,1)=ptemp(1)
c      p(3,2)=ptemp(2)
c      p(3,3)=ptemp(3)
c      p(3,4)=ptemp(4)
c
c!  momentum of the second lepton
c      p(4,4)=pV(4)-p(3,4)
c      p(4,1)=pV(1)-p(3,1)
c      p(4,2)=pV(2)-p(3,2)
c      p(4,3)=pV(3)-p(3,3)
c
cc      print*,'phase space in counterterm'
cc      print*,p(1,1),p(1,2),p(1,3),p(1,4)
cc      print*,p(2,1),p(2,2),p(2,3),p(2,4)
cc      print*,p(3,1),p(3,2),p(3,3),p(3,4)
cc      print*,p(4,1),p(4,2),p(4,3),p(4,4)
cc      print*,'mass',sqrt((p(4,4)+p(3,4))**2
cc     +     - (p(4,1)+p(3,1))**2
cc     +     - (p(4,2)+p(3,2))**2
cc     +     - (p(4,3)+p(3,3))**2)
cc      print*,'y',0.5d0*log((p(4,4)+p(3,4) + (p(4,3)+p(3,3)))/
cc     +     (p(4,4)+p(3,4) - (p(4,3)+p(3,3))))
cc      print*,'pt',sqrt((p(4,1)+p(3,1))**2 + (p(4,2)+p(3,2))**2)
cc      print*
c
cc     End of phase space generation
cc     *****************************************
c      
cC     C Compute Born matrix element
c      if(nproc.eq.3)then
c         call qqb_z(p,msqc)
c      else
c         call qqb_w(p,msqc)
c      endif
c
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
c
cc     End of ME check
cc*******************************************************
      
      if (mode.eq.0) then
         call initsigma(mm,-costh)
      elseif (mode.eq.1) then
         call cthmoments(cthmom0,cthmom1,cthmom2)
         cthmom0=cthmom0*twopi*twopi
         cthmom1=cthmom1*twopi*twopi
         cthmom2=cthmom2*twopi*twopi
         call initsigmacth(mm,cthmom0,cthmom1,cthmom2)
      endif

      if (mode.eq.0.or.mode.eq.1) then
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
      endif
CC  

c The following alfa beta integral depends only on x1 and x2, which depend on y and m
c---  calculate PDF's  
      xx0(1)=dsqrt(q2/sqrts**2)*dexp(+yy)
      xx0(2)=dsqrt(q2/sqrts**2)*dexp(-yy)
c---check if x is out of normal range
      if   ((xx0(1) .gt. 1d0) 
     & .or. (xx0(2) .gt. 1d0)
     & .or. (xx0(1) .lt. xmin)
     & .or. (xx0(2) .lt. xmin)) return
      xx10=xx0(1)
      xx20=xx0(2)
      logxx10=dlog(xx10)
      logxx20=dlog(xx20)
      
      call fdist(ih1,xx10,facscale,fx10)
      call fdist(ih2,xx20,facscale,fx20)
      
C Scaled momentum fractions
      Pqqint1 = Pqqint(xx10)
      D0intx1 = D0int(xx10)
      D1intx1 = D1int(xx10)
      D0intx2 = D0int(xx20)
      Pqqint2 = Pqqint(xx20)
      D1intx2 = D1int(xx20)

c     Preliminary loop for caching
      do ii=1,quadpoints
         z1 = xx10**ctx(ii)
         z2 = xx20**ctx(ii)
         oz1(ii) = 1d0/(1d0-z1)
         oz2(ii) = 1d0/(1d0-z2)
         z2 = xx20**ctx(ii)
         cz1(ii) = z1
         cz2(ii) = z2
         logz1 = dlog(z1)
         logz2 = dlog(z2)
         log1z1(ii) = dlog(1-z1)
         log1z2(ii) = dlog(1-z2)
         call fdist(ih1,xx10**(1-ctx(ii)),facscale,fx1p)
         call fdist(ih2,xx20**(1-ctx(ii)),facscale,fx2p)
         cfx1p(:,ii)=fx1p
         cfx2p(:,ii)=fx2p
         Cqqz1(ii) = Cqq(z1)
         Cqqz2(ii) = Cqq(z2)
c         Cqgz1(ii) = Cqg(z1)
c         Cqgz2(ii) = Cqg(z2)
         Pqqz1(ii) = Pqq(z1)
         Pqqz2(ii) = Pqq(z2)
c         Pqgz1(ii) = Pqg(z1)
c         Pqgz2(ii) = Pqg(z2)

         Pqggqz1(ii) = Pqggq(z1)
         Pqggqz2(ii) = Pqggq(z2)
         Pqqqqz1(ii) = Pqqqq(z1)
         Pqqqqz2(ii) = Pqqqq(z2)
         Pqqqgz1(ii) = Pqqqg(z1)
         Pqqqgz2(ii) = Pqqqg(z2)
         Pqgggz1(ii) = Pqggg(z1)
         Pqgggz2(ii) = Pqggg(z2)
         CqgPgqz1(ii) = CqgPgq(z1)
         CqgPgqz2(ii) = CqgPgq(z2)
         CqqPqqz1(ii) = CqqPqq(z1)
         CqqPqqz2(ii) = CqqPqq(z2)
         CqgPggz1(ii) = CqgPgg(z1)
         CqgPggz2(ii) = CqgPgg(z2)
         CqqPqgz1(ii) = CqqPqg(z1)
         CqqPqgz2(ii) = CqqPqg(z2)
         P2qqVz1(ii) = P2qqV(z1)
         P2qqVz2(ii) = P2qqV(z2)
         P2qqbVz1(ii) = P2qqbV(z1)
         P2qqbVz2(ii) = P2qqbV(z2)
         P2qqSz1(ii) = P2qqS(z1)
         P2qqSz2(ii) = P2qqS(z2)
         P2qgz1(ii) = P2qg(z1)
         P2qgz2(ii) = P2qg(z2)
      enddo
      flgq=1
c     start the fast alfa beta integration
      do ii=1,quadpoints
         z2=cz2(ii)
         fx2p=cfx2p(:,ii)
         do jj=1,quadpoints
            z1=cz1(jj)
            fx1p=cfx1p(:,jj)

            sumfx1p=0
            sumfx2p=0
            do l=1,nf
               sumfx1p=sumfx1p+fx1p(l)+fx1p(-l)
               sumfx2p=sumfx2p+fx2p(l)+fx2p(-l)
            enddo
            
            if (mode.eq.0.or.mode.eq.1) then
C     Start construction of the counterterm
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
c         k = -j
         if(msqc(j,k).eq.0d0) cycle
         
C     Simplest term without convolutions
      tdelta=tdelta+fx10(j)*fx20(k)*msqc(j,k)*flgq
C     Start H1st: to be used later
C     H1st delta term
      tH1st=tH1st+2*C1qqdelta*fx10(j)*fx20(k)*msqc(j,k)*flgq
!add resummation scale dependence
      tH1stQ=tH1stQ-(B1q+A1q/2d0*LQ)*LQ*fx10(j)*fx20(k)*msqc(j,k)*flgq
C     H1st: non delta terms, first leg
      tH1st=tH1st+(fx1p(j)*Cqqz1(jj)*flgq+fx1p(0)*Cqg(z1))
     & *(-logxx10)*fx20(k)*msqc(j,k)
C     H1st: non delta terms, second leg
      tH1st=tH1st+(fx2p(k)*Cqqz2(ii)*flgq+fx2p(0)*Cqg(z2))         
     & *(-logxx20)*fx10(j)*msqc(j,k)
C     H1st: muf dependence (LF factor to be added at the end)
c     gammaqq and gammaqg: first leg      
      diff=-logxx10
     &  *((fx1p(j)-fx10(j)*z1)*Pqqz1(jj)*flgq+fx1p(0)*Pqg(z1))
      tH1stF=tH1stF+diff*fx20(k)*msqc(j,k)
      tH1stF=tH1stF-Pqqint1*fx10(j)*fx20(k)*msqc(j,k)*flgq
c     gammaqq and gammaqg: second leg   
      diff=-logxx20
     &  *((fx2p(k)-fx20(k)*z2)*Pqqz2(ii)*flgq+fx2p(0)*Pqg(z2))
      tH1stF=tH1stF+diff*fx10(j)*msqc(j,k)
      tH1stF=tH1stF-Pqqint2*fx10(j)*fx20(k)*msqc(j,k)*flgq
CC    End of H1st
      if(order.eq.1) cycle
CC    Now (gamma+gamma)*(gamma+gamma) term: to be used later
C     First part: one gamma for each leg: FLGQ here is non trivial ! DONE
      diffg1f=-logxx10*(fx1p(j)-fx10(j)*z1)*Pqqz1(jj)
     &  - Pqqint1*fx10(j)
      diffg10=-logxx10*fx1p(0)*Pqg(z1)
      diffg2f=-logxx20*(fx2p(k)-fx20(k)*z2)*Pqqz2(ii)
     &  - Pqqint2*fx20(k)
      diffg20=-logxx20*fx2p(0)*Pqg(z2)
      tgaga=tgaga+2*
     &   (flgq*diffg10*diffg20+flgq*diffg1f*diffg2f
     &   +diffg10*diffg2f+diffg1f*diffg20)*msqc(j,k)
CC     Second part: gamma*gamma terms
c     Pij * Pjk = D1ijjk (log(1-z)/(1-z))_+ + D0ijjk/(1-z)_+ 
c              + Pijjk(z) + Deltaijjk delta(1-z)
C     First leg
      diff1=-logxx10*(flgq*(fx1p(j)-fx10(j)*z1)
     &    *(D0qqqq*oz1(jj)+D1qqqq*log1z1(jj)*oz1(jj))
     &    +fx1p(j)*Pqqqqz1(jj)*flgq+fx1p(0)*(Pqqqgz1(jj)+Pqgggz1(jj)))
     &    +(Deltaqqqq-D0qqqq*D0intx1-D1qqqq*D1intx1)
     &    *fx10(j)*flgq
C    Second leg
      diff2=-logxx20*(flgq*(fx2p(k)-fx20(k)*z2)
     &    *(D0qqqq*oz2(ii)+D1qqqq*log1z2(ii)*oz2(ii))
     &    +fx2p(k)*Pqqqqz2(ii)*flgq+fx2p(0)*(Pqqqgz2(ii)+Pqgggz2(ii)))
     &    +(Deltaqqqq-D0qqqq*D0intx2-D1qqqq*D1intx2)
     &    *fx20(k)*flgq
C     Include Pqggq
      diff1=diff1-logxx10*(sumfx1p)*Pqggqz1(jj)*flgq
      diff2=diff2-logxx20*(sumfx2p)*Pqggqz2(ii)*flgq
      tgaga=tgaga+diff1*fx20(k)*msqc(j,k)
      tgaga=tgaga+diff2*fx10(j)*msqc(j,k)
C    End of (gamma+gamma)*(gamma+gamma) term: FLGQ non trivial here ! DONE
C    Start  (C+C)*(gamma+gamma) term
c    gamma first leg, C second leg
      diffc2f=-logxx20*fx2p(k)*Cqqz2(ii)+C1qqdelta*fx20(k)
      diffc20=-logxx20*fx2p(0)*Cqg(z2)
      tcga=tcga+msqc(j,k)*
     # (flgq*diffg10*diffc20+flgq*diffg1f*diffc2f
     #          +diffg10*diffc2f+diffg1f*diffc20)
c    C first leg, gamma second leg
      diffc1f=-logxx10*fx1p(j)*Cqqz1(jj)+C1qqdelta*fx10(j)
      diffc10=-logxx10*fx1p(0)*Cqg(z1)
      tcga=tcga+msqc(j,k)*
     # (flgq*diffc10*diffg20+flgq*diffc1f*diffg2f
     #          +diffc10*diffg2f+diffc1f*diffg20)
c    C*gamma: first leg (ignore delta term in Cqq: taken into account with tH1stF)
      tcga=tcga
     &  +(fx1p(j)*CqqPqqz1(jj)*flgq+fx1p(0)*(CqqPqgz1(jj)+CqgPggz1(jj)))
     &  *(-logxx10)*fx20(k)*msqc(j,k) 
c    C*gamma: second leg (ignore delta term in Cqq: taken into account with tH1stF)
      tcga=tcga
     &  +(fx2p(k)*CqqPqqz2(ii)*flgq+fx2p(0)*(CqqPqgz2(ii)+CqgPggz2(ii)))
     &  *(-logxx20)*fx10(j)*msqc(j,k) 
c    Add Cqg*Pgq contribution
      tcga=tcga+(sumfx1p)*CqgPgqz1(jj)
     &           *(-logxx10)*fx20(k)*msqc(j,k)*flgq 
      tcga=tcga+(sumfx2p)*CqgPgqz2(ii)
     &           *(-logxx20)*fx10(j)*msqc(j,k)*flgq 
      
CC  Start 2-loop AP
C   Gluon + pure singlet
      tgamma2=tgamma2+fx1p(0)*P2qgz1(jj)
     & *(-logxx10)*fx20(k)*msqc(j,k)
      tgamma2=tgamma2+fx2p(0)*P2qgz2(ii)
     & *(-logxx20)*fx10(j)*msqc(j,k)
      tgamma2=tgamma2+sumfx1p*P2qqSz1(jj)
     & *(-logxx10)*fx20(k)*msqc(j,k)*flgq
      tgamma2=tgamma2+sumfx2p*P2qqSz2(ii)
     & *(-logxx20)*fx10(j)*msqc(j,k)*flgq
C   P2qq non-singlet: regular part
      tgamma2=tgamma2+fx1p(j)*P2qqVz1(jj)
     & *(-logxx10)*fx20(k)*msqc(j,k)*flgq
      tgamma2=tgamma2+fx2p(k)*P2qqVz2(ii)
     & *(-logxx20)*fx10(j)*msqc(j,k)*flgq
C   P2qq non-singlet: 1/(1-z)_+
      diff=-logxx10
     &  *(fx1p(j)-fx10(j)*z1)*oz1(jj)
     &  - D0intx1*fx10(j)      
      tgamma2=tgamma2+2d0/3*Kappa*diff*fx20(k)*msqc(j,k)*flgq
      diff=-logxx20
     &  *(fx2p(k)-fx20(k)*z2)*oz2(ii)
     &  - D0intx2*fx20(k)      
      tgamma2=tgamma2+2d0/3*Kappa*diff*fx10(j)*msqc(j,k)*flgq
C   P2qqb non singlet
      tgamma2=tgamma2+fx1p(-j)*P2qqbVz1(jj)
     & *(-logxx10)*fx20(k)*msqc(j,k)*flgq
      tgamma2=tgamma2+fx2p(-k)*P2qqbVz2(ii)
     & *(-logxx20)*fx10(j)*msqc(j,k)*flgq

      enddo
      enddo

CM 7/11  Resummation scale dependence added 
CC   First order
      sig12=-0.5d0*A1q*tdelta
      sig11=-(B1q+A1q*LQ)*tdelta-tH1stF
CC   Second order
      sig24=(A1q)**2/8*tdelta
      sig23=-beta0*A1q/3*tdelta-0.5d0*A1q*sig11
      sig22=0.5d0*(beta0*A1q*(LR-LQ)-A2q)*tdelta
     &     -0.5d0*A1q*(tH1st+tH1stQ+(LF-LQ)*tH1stF)
     &     -0.5d0*(B1q+A1q*LQ-beta0)*sig11
     &     +0.5d0*(B1q+A1q*LQ)*tH1stF
     &     +0.5d0*tgaga
      sig21=-beta0*(LR-LQ)*sig11
     &     -(B1q+A1q*LQ)*(tH1st+tH1stQ+(LF-LQ)*tH1stF)
     &     -(LF-LQ)*tgaga-(B2q+A2q*LQ)*tdelta+beta0*tH1st-tcga-tgamma2
     &     +(B1q+0.5d0*A1q*LQ)*LQ*tH1stF  
c     Include missing delta term from C*gamma (no factor 2 here !)
      sig21=sig21-C1qqdelta*tH1stF
C     Include missing term from contact term in 2 loop AP
      sig21=sig21-2*Delta2qq*tdelta
CC Include as/pi factors and sum O(as) and O(as^2) contributions
      sig1=sig12*LL2+sig11*LL1
      sig2=sig24*LL4+sig23*LL3+sig22*LL2+sig21*LL1
      sig1=sig1*ason2pi*2
      sig2=sig2*(ason2pi*2)**2
      if(order.eq.1)then
         xmsq(1)=xmsq(1)-sig1*ctw(ii)*ctw(jj)
      else
         xmsq(1)=xmsq(1)-(sig1+sig2)*ctw(ii)*ctw(jj)
      endif

      else if (mode.eq.2) then
         do j=-nf,nf
            do k=-nf,nf
c     k = -j
               if(LL1jk(j,k).eq.0d0) cycle
C     Simplest term without convolutions
      tdeltajk(j,k)=fx10(j)*fx20(k)*flgq
C     Start H1st: to be used later
C     H1st delta term
      tH1stjk(j,k)=2*C1qqdelta*fx10(j)*fx20(k)*flgq
!add resummation scale dependence
      tH1stQjk(j,k)=-(B1q+A1q/2d0*LQ)*LQ*fx10(j)*fx20(k)*flgq
C     H1st: non delta terms, first leg
      tH1stjk(j,k)=tH1stjk(j,k)+(fx1p(j)*Cqqz1(jj)*flgq+fx1p(0)*Cqg(z1))
     & *(-logxx10)*fx20(k)
C     H1st: non delta terms, second leg
      tH1stjk(j,k)=tH1stjk(j,k)+(fx2p(k)*Cqqz2(ii)*flgq+fx2p(0)*Cqg(z2))         
     & *(-logxx20)*fx10(j)
C     H1st: muf dependence (LF factor to be added at the end)
c     gammaqq and gammaqg: first leg      
      diff=-logxx10
     &  *((fx1p(j)-fx10(j)*z1)*Pqqz1(jj)*flgq+fx1p(0)*Pqg(z1))
      tH1stFjk(j,k)=diff*fx20(k)
      tH1stFjk(j,k)=tH1stFjk(j,k)-Pqqint1*fx10(j)*fx20(k)*flgq
c     gammaqq and gammaqg: second leg   
      diff=-logxx20
     &  *((fx2p(k)-fx20(k)*z2)*Pqqz2(ii)*flgq+fx2p(0)*Pqg(z2))
      tH1stFjk(j,k)=tH1stFjk(j,k)+diff*fx10(j)
      tH1stFjk(j,k)=tH1stFjk(j,k)-Pqqint2*fx10(j)*fx20(k)*flgq
CC    End of H1st
      if(order.eq.1) cycle
CC    Now (gamma+gamma)*(gamma+gamma) term: to be used later
C     First part: one gamma for each leg: FLGQ here is non trivial ! DONE
      diffg1f=-logxx10*(fx1p(j)-fx10(j)*z1)*Pqqz1(jj)
     &  - Pqqint1*fx10(j)
      diffg10=-logxx10*fx1p(0)*Pqg(z1)
      diffg2f=-logxx20*(fx2p(k)-fx20(k)*z2)*Pqqz2(ii)
     &  - Pqqint2*fx20(k)
      diffg20=-logxx20*fx2p(0)*Pqg(z2)
      tgagajk(j,k)=2*
     &   (flgq*diffg10*diffg20+flgq*diffg1f*diffg2f
     &   +diffg10*diffg2f+diffg1f*diffg20)
CC     Second part: gamma*gamma terms
c     Pij * Pjk = D1ijjk (log(1-z)/(1-z))_+ + D0ijjk/(1-z)_+ 
c              + Pijjk(z) + Deltaijjk delta(1-z)
C     First leg
      diff1=-logxx10*(flgq*(fx1p(j)-fx10(j)*z1)
     &    *(D0qqqq*oz1(jj)+D1qqqq*log1z1(jj)*oz1(jj))
     &    +fx1p(j)*Pqqqqz1(jj)*flgq+fx1p(0)*(Pqqqgz1(jj)+Pqgggz1(jj)))
     &    +(Deltaqqqq-D0qqqq*D0intx1-D1qqqq*D1intx1)
     &    *fx10(j)*flgq
C    Second leg
      diff2=-logxx20*(flgq*(fx2p(k)-fx20(k)*z2)
     &    *(D0qqqq*oz2(ii)+D1qqqq*log1z2(ii)*oz2(ii))
     &    +fx2p(k)*Pqqqqz2(ii)*flgq+fx2p(0)*(Pqqqgz2(ii)+Pqgggz2(ii)))
     &    +(Deltaqqqq-D0qqqq*D0intx2-D1qqqq*D1intx2)
     &    *fx20(k)*flgq
C     Include Pqggq
      diff1=diff1-logxx10*(sumfx1p)*Pqggqz1(jj)*flgq
      diff2=diff2-logxx20*(sumfx2p)*Pqggqz2(ii)*flgq
      
      tgagajk(j,k)=tgagajk(j,k)+diff1*fx20(k)
      tgagajk(j,k)=tgagajk(j,k)+diff2*fx10(j)
C    End of (gamma+gamma)*(gamma+gamma) term: FLGQ non trivial here ! DONE
C    Start  (C+C)*(gamma+gamma) term
c    gamma first leg, C second leg
      diffc2f=-logxx20*fx2p(k)*Cqqz2(ii)+C1qqdelta*fx20(k)
      diffc20=-logxx20*fx2p(0)*Cqg(z2)
      tcgajk(j,k)=
     # (flgq*diffg10*diffc20+flgq*diffg1f*diffc2f
     #          +diffg10*diffc2f+diffg1f*diffc20)
c    C first leg, gamma second leg
      diffc1f=-logxx10*fx1p(j)*Cqqz1(jj)+C1qqdelta*fx10(j)
      diffc10=-logxx10*fx1p(0)*Cqg(z1)
      tcgajk(j,k)=tcgajk(j,k)+
     # (flgq*diffc10*diffg20+flgq*diffc1f*diffg2f
     #          +diffc10*diffg2f+diffc1f*diffg20)
c    C*gamma: first leg (ignore delta term in Cqq: taken into account with tH1stF)
      tcgajk(j,k)=tcgajk(j,k)
     &  +(fx1p(j)*CqqPqqz1(jj)*flgq+fx1p(0)*(CqqPqgz1(jj)+CqgPggz1(jj)))
     &  *(-logxx10)*fx20(k) 
c    C*gamma: second leg (ignore delta term in Cqq: taken into account with tH1stF)
      tcgajk(j,k)=tcgajk(j,k)
     &  +(fx2p(k)*CqqPqqz2(ii)*flgq+fx2p(0)*(CqqPqgz2(ii)+CqgPggz2(ii)))
     &  *(-logxx20)*fx10(j) 
c    Add Cqg*Pgq contribution
      tcgajk(j,k)=tcgajk(j,k)+(sumfx1p)*CqgPgqz1(jj)
     &           *(-logxx10)*fx20(k)*flgq 
      tcgajk(j,k)=tcgajk(j,k)+(sumfx2p)*CqgPgqz2(ii)
     &           *(-logxx20)*fx10(j)*flgq 
CC  Start 2-loop AP
C     Gluon + pure singlet
      tgamma2jk(j,k)=0
      tgamma2jk(j,k)=tgamma2jk(j,k)+fx1p(0)*P2qgz1(jj)
     & *(-logxx10)*fx20(k)
      tgamma2jk(j,k)=tgamma2jk(j,k)+fx2p(0)*P2qgz2(ii)
     & *(-logxx20)*fx10(j)
      tgamma2jk(j,k)=tgamma2jk(j,k)+(sumfx1p)*P2qqSz1(jj)
     & *(-logxx10)*fx20(k)*flgq
      tgamma2jk(j,k)=tgamma2jk(j,k)+(sumfx2p)*P2qqSz2(ii)
     & *(-logxx20)*fx10(j)*flgq
C   P2qq non-singlet: regular part
      tgamma2jk(j,k)=tgamma2jk(j,k)+fx1p(j)*P2qqVz1(jj)
     & *(-logxx10)*fx20(k)*flgq
      tgamma2jk(j,k)=tgamma2jk(j,k)+fx2p(k)*P2qqVz2(ii)
     & *(-logxx20)*fx10(j)*flgq
C   P2qq non-singlet: 1/(1-z)_+
      diff=-logxx10
     &  *(fx1p(j)-fx10(j)*z1)*oz1(jj)
     &  - D0intx1*fx10(j)      
      tgamma2jk(j,k)=tgamma2jk(j,k)+2d0/3*Kappa*diff*fx20(k)*flgq
      diff=-logxx20
     &  *(fx2p(k)-fx20(k)*z2)*oz2(ii)
     &  - D0intx2*fx20(k)      
      tgamma2jk(j,k)=tgamma2jk(j,k)+2d0/3*Kappa*diff*fx10(j)*flgq
C   P2qqb non singlet
      tgamma2jk(j,k)=tgamma2jk(j,k)+fx1p(-j)*P2qqbVz1(jj)
     & *(-logxx10)*fx20(k)*flgq
      tgamma2jk(j,k)=tgamma2jk(j,k)+fx2p(-k)*P2qqbVz2(ii)
     & *(-logxx20)*fx10(j)*flgq
      enddo
      enddo

      sig1 = 0
      sig2 = 0
      do j=-nf,nf
         do k=-nf,nf
            if(LL1jk(j,k).eq.0d0) cycle

CM 7/11  Resummation scale dependence added 
CC   First order
      sig12=-0.5d0*A1q*tdeltajk(j,k)
      sig11=-(B1q+A1q*LQ)*tdeltajk(j,k)
     .     -tH1stFjk(j,k)
CC   Second order
      sig24=(A1q)**2/8*tdeltajk(j,k)
      sig23=-beta0*A1q/3*tdeltajk(j,k)-0.5d0*A1q*sig11
      sig22=0.5d0*(beta0*A1q*(LR-LQ)-A2q)*tdeltajk(j,k)
     &     -0.5d0*A1q*(tH1stjk(j,k)
     &     +tH1stQjk(j,k)+(LF-LQ)*tH1stFjk(j,k))
     &     -0.5d0*(B1q+A1q*LQ-beta0)*sig11
     &     +0.5d0*(B1q+A1q*LQ)*tH1stFjk(j,k)
     &     +0.5d0*tgagajk(j,k)
      sig21=-beta0*(LR-LQ)*sig11
     &     -(B1q+A1q*LQ)*(tH1stjk(j,k)
     &     +tH1stQjk(j,k)+(LF-LQ)*tH1stFjk(j,k))
     &     -(LF-LQ)*tgagajk(j,k)-(B2q+A2q*LQ)*tdeltajk(j,k)
     &     +beta0*tH1stjk(j,k)
     &     -tcgajk(j,k)-tgamma2jk(j,k)
     &     +(B1q+0.5d0*A1q*LQ)*LQ*tH1stFjk(j,k)
c     Include missing delta term from C*gamma (no factor 2 here !)
      sig21=sig21-C1qqdelta*tH1stFjk(j,k)
C     Include missing term from contact term in 2 loop AP
      sig21=sig21-2*Delta2qq*tdeltajk(j,k)
CC Include as/pi factors and sum O(as) and O(as^2) contributions
      sig1=sig1+sig12*LL2jk(j,k)+sig11*LL1jk(j,k)
      sig2=sig2+sig24*LL4jk(j,k)+sig23*LL3jk(j,k)
     &     +sig22*LL2jk(j,k)+sig21*LL1jk(j,k)
      enddo
      enddo
      sig1=sig1*ason2pi*2
      sig2=sig2*(ason2pi*2)**2

      if(order.eq.1)then
         xmsq(1)=xmsq(1)-sig1*ctw(ii)*ctw(jj)
      else
         xmsq(1)=xmsq(1)-(sig1+sig2)*ctw(ii)*ctw(jj)
      endif


c      print *,tdelta,tH1st,tH1stF,tH1stQ,tgaga,tcga,tgamma2

         
      endif
      enddo
      enddo
c     end of alfa beta loops
      
CC Include iacobians (do not include jacobians)
      shad=sqrts**2
      Vol=1d0
      xmsq(1)=xmsq(1)*q2/shad/Vol

      countterm=0d0
      xint=0d0

C Flux for Born cross section
      fluxborn=fbGeV2/(2*q2)
C Multiply by BORN phase space weight
      xmsq(1)=xmsq(1)*fluxborn/BrnRat*((one/4d0/pi)**3)

! SWITCHING
      xmsq(1)=xmsq(1)*switch

c---Add to total
      xint=xmsq(1)
        
      countterm=xint
c      print *, qtt, mm, yy, countterm

      return
      end

c initialize the points of the gaussian quadrature for the alfa and beta integration
      subroutine ctquadinit
      implicit none
      double precision min,max
      double precision a,b,c,m,x,t,jac
      integer i,j
      include 'gauss.inc'

c     Common block (output)
      integer ctintervals,ctrule,ctdim
      parameter (ctintervals=20)
      parameter (ctrule=10)
      parameter (ctdim=ctrule*ctintervals)
      double precision ctx(ctdim)
      double precision ctw(ctdim)
      integer quadpoints
      common/ctweights/ctx,ctw,quadpoints

      min = 1d-7
      max = 1d0
      quadpoints = ctrule*ctintervals
c      do i=1,ctintervals
c         a = min*((max/min)**(real(i-1, 8)/real(ctintervals,8)))
c         b = min*((max/min)**(real(i,8)/real(ctintervals,8)))
c         c=0.5d0*(a+b)
c         m=0.5d0*(b-a)
c         do j=1,ctrule
c            ctx(j+(i-1)*ctrule)=c+m*xxx10(j)
c            ctw(j+(i-1)*ctrule)=www10(j)*m
c         enddo
c      enddo

      do i=1,ctintervals
         a = min+(max-min)*(i-1)/ctintervals
         b = min+(max-min)*i/ctintervals
         c=0.5d0*(a+b)
         m=0.5d0*(b-a)
         do j=1,ctrule
            x=c+m*xxx10(j)
            t=min*(max/min)**x
            jac=t*log(max/min)
            ctx(j+(i-1)*ctrule)=t
            ctw(j+(i-1)*ctrule)=www10(j)*m*jac
         enddo
      enddo
      
      return
      end


c perform qt integration
      subroutine ctqtint(m, y, qtmin, qtmax)
      implicit none
      double precision qtmin,qtmax
      double precision y,m
      double precision xa,xb,xc,xm
      double precision min,max
      double precision q2
      double precision qt,qt2,qta,qtb,qtx,jac
      double precision qtmax2,qtmin2,tiny
      double precision x,w
      double precision switch,xmio
      common/xmio/xmio
      integer qtintervals,qtrule
      parameter (qtintervals=2)
      parameter (qtrule=4)
      integer i,j

      double precision Itilde
      external Itilde
      
      common/a_param/a_param,b0p
      double precision a_param,b0p
      
      include 'constants.f'
      double precision cthmom0,cthmom1,cthmom2
      double precision msqc(-nf:nf,-nf:nf)
      integer jj,kk
      COMMON/SIGMAIJ/SIGMAIJ
      double precision sigmaij(-5:5,-5:5)
      integer nproc
      common/nproc/nproc

      include 'gauss.inc'
c     Common block (output)
      double precision LL1,LL2,LL3,LL4
      double precision LL1jk(-nf:nf,-nf:nf),LL2jk(-nf:nf,-nf:nf)
      double precision LL3jk(-nf:nf,-nf:nf),LL4jk(-nf:nf,-nf:nf)
      common/largelogs/LL1jk,LL2jk,LL3jk,LL4jk

      do jj=-nf,nf
         do kk=-nf,nf
            LL1jk(jj,kk)=0
            LL2jk(jj,kk)=0
            LL3jk(jj,kk)=0
            LL4jk(jj,kk)=0
         enddo
      enddo
      q2=m*m
      qtmin2 = qtmin**2
      qtmax2 = qtmax**2
      tiny = 1D-5;
      qta = 1d0/(1d0+log(qtmax2/tiny))
      qtb = 1d0/(1d0+log(qtmin2/tiny))
      min = 0
      max = 1
      do i=1,qtintervals
         xa = min+(max-min)*(i-1)/qtintervals
         xb = min+(max-min)*i/qtintervals
         xc=0.5d0*(xa+xb)
         xm=0.5d0*(xb-xa)
         do j=1,qtrule
            x=xc+xm*xxx4(j)
            qtx = qta + (qtb-qta) * x
            qt2=tiny*exp(1d0/qtx - 1d0)
            jac=(qtb-qta)*qt2/qtx**2

            w=www4(j)*xm*jac

            qt = sqrt(qt2)

!      SWITCHING FUNCTIONS
            switch=1d0
            if(qt.ge.m*3/4d0)  switch=dexp(-(m*3/4d0-qt)**2/(m/2d0)**2) ! GAUSS SWITCH
!!!!!!!!!!!!
!     if(qt.ge.m)  switch=dexp(-(m-qt)**2/(m/2d0)**2)           ! GAUSS SWITCH
!     if(qt.ge.m/2d0)  switch=dexp(-(m/2d0-qt)**2/(m/2d0)**2)              ! GAUSS SWITCH
!     if(qt.ge.m/2d0)        switch=dexp(-2d0*(m/2d0-qt)**2/(m/2d0)**2)    ! GAUSS SWITCH FASTER
!     if(qt.ge.m/2d0)        switch=dexp((m2/4d0-qt2)/m2)                  ! EXP SWITCH
!     if(qt.ge.m/2d0)        switch=dexp((m2/4d0-qt2)/m2*4d0)              ! EXP SWITCH FASTER
!     if(qt.ge.m/2d0)        switch=(dcos(pi/50d0*(qt-45.6d0))+1d0)/2d0    ! COS SWITCH
!     if(qt.ge.95.6)         switch=0d0                                    ! COS SWITCH

            if(switch.le.0.01d0) cycle

c     xmio is used in besselkfast for Itilde
            xmio=dsqrt(qt2/(q2/a_param**2))

            LL1=Itilde(1)/q2**2*a_param**2*w*switch
            LL2=Itilde(2)/q2**2*a_param**2*w*switch
            LL3=Itilde(3)/q2**2*a_param**2*w*switch
            LL4=Itilde(4)/q2**2*a_param**2*w*switch

            call setqt(qt)
            call genV4p()
            call cthmoments(cthmom0,cthmom1,cthmom2)
            cthmom0=cthmom0*twopi*twopi
            cthmom1=cthmom1*twopi*twopi
            cthmom2=cthmom2*twopi*twopi
            call initsigmacth(m,cthmom0,cthmom1,cthmom2)

            do jj=-nf,nf
               do kk=-nf,nf
                  msqc(jj,kk) = 0d0
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
            
            do jj=-nf,nf
               do kk=-nf,nf
                  if (msqc(jj,kk).ne.0d0) then
                     LL1jk(jj,kk)=LL1jk(jj,kk)+LL1*msqc(jj,kk)
                     LL2jk(jj,kk)=LL2jk(jj,kk)+LL2*msqc(jj,kk)
                     LL3jk(jj,kk)=LL3jk(jj,kk)+LL3*msqc(jj,kk)
                     LL4jk(jj,kk)=LL4jk(jj,kk)+LL4*msqc(jj,kk)
                  endif
               enddo
            enddo
         enddo
      enddo
      return
      end
