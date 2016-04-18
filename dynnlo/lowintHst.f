C     Version that allows to separate the channels
C     Scale dependence included

C     March 2015: Bug in muf dependence corrected

      double precision function lowintHst_dynnlo(r,wgt,f)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'limits.f'
      include 'npart.f'
      include 'vegas_common.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'facscale.f'
      include 'noglue.f'
      include 'process.f'
      include 'phasemin.f'

C
      include 'qcdcouple.f'
      include 'rescoeff.f'
      include 'dynamicscale.f'
      include 'options.f'

      double precision f(*)
      integer ih1,ih2,j,k,l,nvec,flgq
      double precision r(mxdim),W,sqrts,xmsq,val,
     . fx10(-nf:nf),fx20(-nf:nf),p(mxpart,4),pjet(mxpart,4),
     . pswt
      double precision wgt,msqc(-nf:nf,-nf:nf)
      double precision xx(2),flux,BrnRat
      logical includedipole
CC
      logical cuts
      external cuts
      double precision x1p,x2p,fx1p(-nf:nf),fx2p(-nf:nf)
      double precision asopi,z1,z2,alfa,beta,cut,diff
      double precision tdelta,tH1st,tH1stF,xx10,xx20,tH2st
      double precision tgaga,tcga,tgamma2
      double precision diff10,diff20,diffc10,diffc20,diffg10,diffg20
      double precision diff1f,diff2f,diffg1f,diffg2f,diffc1f,diffc2f
      double precision Pggreg,D0int,D1int,Cgq,Pgq,LF,LR
      double precision dot,q2,Ggq,Ggg,spgq1,spgq2,tH2sp
      double precision Pqqint,Cqq,Cqg,Pqq,dyPqg
      double precision C2qqreg,C2qqp,C2qqb,C2qg
C
      double precision diff1,diff2     
      double precision Pqqqq,Pqqqg,Pqggq,Pqggg
      double precision CqqPqq,CqqPqg,CqgPgq,CqgPgg
      double precision P2qg,P2qqV,P2qqbV,P2qqS
C
      double precision beta1,H2qqdelta,H2qqD0
      common/Hstcoeff/beta1,H2qqdelta,H2qqD0

      integer order,a,b
      common/nnlo/order
CC
      integer ndec,nproc
      common/nproc/nproc
CC
      common/density/ih1,ih2
      common/energy/sqrts
      common/x1x2/xx
      common/BrnRat/BrnRat
     

      data p/48*0d0/

      logical binner
      external binner
      external hists_fill

      integer npdf,maxpdf
      
      lowintHst_dynnlo=0d0
      do npdf=0,totpdf-1
         f(npdf+1)=0d0
      enddo

      W=sqrts**2

      npart=2
      call gen2(r,p,pswt,*999)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      nvec=npart+2
      call dotem(nvec,p,s)

      call masscuts(s,*999)
      
c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead bail out
      if (cuts(p,0)) goto 999
      if (.not.binner(p(3,:),p(4,:))) goto 999

      
      xx(1)=-2d0*p(1,4)/sqrts
      xx(2)=-2d0*p(2,4)/sqrts

c     Load central PDF and QCD coupling (not needed)
      if (pdferr) then
         call dysetpdf(0)
      endif

c--- Calculate the required matrix element      
      if(nproc.eq.3) then
         call qqb_z(p,msqc)
      else
         call qqb_w(p,msqc)
      endif
            
      flux=fbGeV2/(2d0*xx(1)*xx(2)*W)
      
C     Compute Q2

      q2=2*dot(p,3,4)

C     Dynamic scale

      if(dynamicscale) call scaleset(q2)

      LF=dlog(q2/facscale**2)
      LR=dlog(q2/scale**2)


C Scaled momentum fractions

      cut=1d-8
   
C ndim here is 6 as for H->2gamma


!      beta=cut+(1-cut)*r(ndim-1)
!      alfa=cut+(1-cut)*r(ndim)
      beta=cut+(1-cut)*r(5)
      alfa=cut+(1-cut)*r(6)


      xx10=xx(1)
      xx20=xx(2)

      z1=xx10**beta
      z2=xx20**alfa

c     skip PDF loop in the preconditioning phase
      maxpdf=0
      if (doFill.ne.0) maxpdf = totpdf-1
      
c     start PDF loop
      do npdf=0,maxpdf
         call dysetpdf(npdf)
         asopi=ason2pi*2
         call hists_setpdf(npdf)
c     intitialise xmsq to 0
         xmsq=0d0
      
c--- calculate PDF's  

      call fdist(ih1,xx10,facscale,fx10)
      call fdist(ih2,xx20,facscale,fx20)

      call fdist(ih1,xx10**(1-beta),facscale,fx1p)
      call fdist(ih2,xx20**(1-alfa),facscale,fx2p)

       if(noglue) then
        fx10(0)=0d0
        fx20(0)=0d0
        fx1p(0)=0d0
        fx2p(0)=0d0
       endif

       if(ggonly) then
        do j=1,nf
        fx10(j)=0d0
        fx10(-j)=0d0
        fx20(j)=0d0
        fx20(-j)=0d0
        fx1p(j)=0d0
        fx1p(-j)=0d0
        fx2p(j)=0d0
        fx2p(-j)=0d0   
        enddo
       endif

CC     TIENI SOLO uubar
c        do j=-nf,1
c        fx10(j)=0d0
c        fx1p(j)=0d0
c        enddo
c        do j=3,nf
c        fx10(j)=0d0
c        fx1p(j)=0d0
c        enddo
c        do j=-nf,-3
c        fx20(j)=0d0
c        fx2p(j)=0d0
c        enddo
c        do j=-1,nf
c        fx20(j)=0d0
c        fx2p(j)=0d0
c        enddo
CC


        flgq=1
        if(gqonly)flgq=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CC Start calculation

      tdelta=0d0
      tH1st=0d0
      tH1stF=0d0
      tH2st=0d0
      tcga=0d0
      tgamma2=0d0
      tgaga=0d0



      do j=-nf,nf
      do k=-nf,nf

      if(msqc(j,k).eq.0d0) goto 75


C     Simplest term without convolutions
  
      tdelta=tdelta+fx10(j)*fx20(k)*msqc(j,k)*flgq

      if(order.eq.0) goto 75

C     Start H1st: to be used later

C     H1st delta term

      tH1st=tH1st+2*C1qqdelta*fx10(j)*fx20(k)*msqc(j,k)*flgq

C     H1st: non delta terms, first leg


      tH1st=tH1st+(fx1p(j)*Cqq(z1)*flgq+fx1p(0)*Cqg(z1))
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)


C     H1st: non delta terms, second leg


      tH1st=tH1st+(fx2p(k)*Cqq(z2)*flgq+fx2p(0)*Cqg(z2))         
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)
      

C     H1st: muf dependence (LF factor to be added at the end)


c     gammaqq and gammaqg: first leg      


      diff=-dlog(xx10)
     &  *((fx1p(j)-fx10(j)*xx10**beta)*Pqq(z1)*flgq+fx1p(0)*dyPqg(z1))
      tH1stF=tH1stF+diff*fx20(k)*msqc(j,k)
      tH1stF=tH1stF-Pqqint(xx10)*fx10(j)*fx20(k)*msqc(j,k)*flgq

c     gammaqq and gammaqg: second leg   


      diff=-dlog(xx20)
     &  *((fx2p(k)-fx20(k)*xx20**alfa)*Pqq(z2)*flgq+fx2p(0)*dyPqg(z2))
      tH1stF=tH1stF+diff*fx10(j)*msqc(j,k)
      tH1stF=tH1stF-Pqqint(xx20)*fx10(j)*fx20(k)*msqc(j,k)*flgq

      if(order.eq.1) goto 75

CC    End of H1st

CC    Start H2 contribution

CC    H2st gg contribution

      tH2st=tH2st+fx1p(0)*Cqg(z1)*(-dlog(xx10))*
     &            fx2p(0)*Cqg(z2)*(-dlog(xx20))*msqc(j,k)*flgq

CC    H2st qqbar contribution from C1*C1 (without delta term)

C     regular*regular

      tH2st=tH2st+fx1p(j)*Cqq(z1)*(-dlog(xx10))*
     &            fx2p(k)*Cqq(z2)*(-dlog(xx20))*msqc(j,k)*flgq

C     regular-delta

      tH2st=tH2st+fx1p(j)*Cqq(z1)*(-dlog(xx10))*
     &            fx20(k)*C1qqdelta*msqc(j,k)*flgq       

      tH2st=tH2st+fx2p(k)*Cqq(z2)*(-dlog(xx20))*
     &            fx10(j)*C1qqdelta*msqc(j,k)*flgq       


CC    H2st qg contribution from C1*C1

C     regular*regular

      tH2st=tH2st+fx1p(0)*Cqg(z1)*(-dlog(xx10))*
     &            fx2p(k)*Cqq(z2)*(-dlog(xx20))*msqc(j,k)

      tH2st=tH2st+fx1p(j)*Cqq(z1)*(-dlog(xx10))*
     &            fx2p(0)*Cqg(z2)*(-dlog(xx20))*msqc(j,k)


C     regular-delta

      tH2st=tH2st+fx1p(0)*Cqg(z1)*(-dlog(xx10))*
     &            fx20(k)*C1qqdelta*msqc(j,k)       

      tH2st=tH2st+fx2p(0)*Cqg(z2)*(-dlog(xx20))*
     &            fx10(j)*C1qqdelta*msqc(j,k)    

CC    H2st qqbar channel: D0(z), first leg

      diff=-dlog(xx10)*(fx1p(j)-fx10(j)*xx10**beta)*H2qqD0/(1-z1)

      tH2st=tH2st+0.5d0*diff*fx20(k)*msqc(j,k)*flgq
      tH2st=tH2st-0.5d0*H2qqD0*D0int(xx10)
     &            *fx10(j)*fx20(k)*msqc(j,k)*flgq

CC    H2st, qqbar channel: D0(z), second leg
      
      diff=-dlog(xx20)*(fx2p(k)-fx20(k)*xx20**alfa)*H2qqD0/(1-z2)

      tH2st=tH2st+0.5d0*diff*fx10(j)*msqc(j,k)*flgq
      tH2st=tH2st-0.5d0*H2qqD0*D0int(xx20)
     &            *fx10(j)*fx20(k)*msqc(j,k)*flgq

CC    C2qq, regular part, first leg

      tH2st=tH2st+fx1p(j)*C2qqreg(z1)
     &                   *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq

CC    C2qq, regular part, second leg

      tH2st=tH2st+fx2p(k)*C2qqreg(z2)
     &                   *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq

CC    C2qg, first leg

      tH2st=tH2st+fx1p(0)*C2qg(z1)*(-dlog(xx10))*fx20(k)*msqc(j,k)

CC    C2qg, second leg

      tH2st=tH2st+fx2p(0)*C2qg(z2)*(-dlog(xx20))*fx10(j)*msqc(j,k)

CC    Cqqbar contribution: first leg

      tH2st=tH2st+fx1p(-j)*C2qqb(z1)
     &                    *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq

CC    Cqqbar contribution: second leg
  
      tH2st=tH2st+fx2p(-k)*C2qqb(z2)
     &                    *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq  

      do a=1,nf

CC    Cqqp contribution: first leg

      if(a.ne.abs(j)) then      
       tH2st=tH2st+(fx1p(a)+fx1p(-a))*
     &       C2qqp(z1)*(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq
      endif

CC    Cqqp contribution: second leg

      if(a.ne.abs(k)) then      
       tH2st=tH2st+(fx2p(a)+fx2p(-a))*
     &       C2qqp(z2)*(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq
      endif

      enddo

CCCC Terms needed for NNLO scale dependence  CCCCCC


CC    (gamma+gamma)*(gamma+gamma) term

C     First part: one gamma for each leg


      diffg1f=-dlog(xx10)*(fx1p(j)-fx10(j)*xx10**beta)*Pqq(z1)
     &  - Pqqint(xx10)*fx10(j)


      diffg10=-dlog(xx10)*fx1p(0)*dyPqg(z1)

      diffg2f=-dlog(xx20)*(fx2p(k)-fx20(k)*xx20**alfa)*Pqq(z2)
     &  - Pqqint(xx20)*fx20(k)


      diffg20=-dlog(xx20)*fx2p(0)*dyPqg(z2)


      tgaga=tgaga+2*
     #   (flgq*diffg10*diffg20+flgq*diffg1f*diffg2f
     #   +diffg10*diffg2f+diffg1f*diffg20)*msqc(j,k)


CC     Second part: gamma*gamma terms

c     Pij * Pjk = D1ijjk (log(1-z)/(1-z))_+ + D0ijjk/(1-z)_+ 
c              + Pijjk(z) + Deltaijjk delta(1-z)

C     First leg

      
      diff1=-dlog(xx10)*(flgq*(fx1p(j)-fx10(j)*xx10**beta)
     &    *(D0qqqq/(1-z1)+D1qqqq*dlog(1-z1)/(1-z1))
     &    +fx1p(j)*Pqqqq(z1)*flgq+fx1p(0)*(Pqqqg(z1)+Pqggg(z1)))
     &    +(Deltaqqqq-D0qqqq*D0int(xx10)-D1qqqq*D1int(xx10))
     &    *fx10(j)*flgq


C    Second leg

      
      diff2=-dlog(xx20)*(flgq*(fx2p(k)-fx20(k)*xx20**alfa)
     &    *(D0qqqq/(1-z2)+D1qqqq*dlog(1-z2)/(1-z2))
     &    +fx2p(k)*Pqqqq(z2)*flgq+fx2p(0)*(Pqqqg(z2)+Pqggg(z2)))
     &    +(Deltaqqqq-D0qqqq*D0int(xx20)-D1qqqq*D1int(xx20))
     &    *fx20(k)*flgq


C     Include Pqggq

      do l=1,nf
      diff1=diff1-dlog(xx10)*(fx1p(l)+fx1p(-l))*Pqggq(z1)*flgq
      diff2=diff2-dlog(xx20)*(fx2p(l)+fx2p(-l))*Pqggq(z2)*flgq
      enddo

      tgaga=tgaga+diff1*fx20(k)*msqc(j,k)
      tgaga=tgaga+diff2*fx10(j)*msqc(j,k)



C    End of (gamma+gamma)*(gamma+gamma) term

C    Start  (C+C)*(gamma+gamma) term

c    gamma first leg, C second leg


      diffc2f=-dlog(xx20)*fx2p(k)*Cqq(z2)+C1qqdelta*fx20(k)

      diffc20=-dlog(xx20)*fx2p(0)*Cqg(z2)


      tcga=tcga+msqc(j,k)*
     # (flgq*diffg10*diffc20+flgq*diffg1f*diffc2f
     #          +diffg10*diffc2f+diffg1f*diffc20)


c    C first leg, gamma second leg

      diffc1f=-dlog(xx10)*fx1p(j)*Cqq(z1)+C1qqdelta*fx10(j)

      diffc10=-dlog(xx10)*fx1p(0)*Cqg(z1)

      tcga=tcga+msqc(j,k)*
     # (flgq*diffc10*diffg20+flgq*diffc1f*diffg2f
     #          +diffc10*diffg2f+diffc1f*diffg20)
    

c    C*gamma: first leg (ignore delta term in Cqq: taken into account with tH1stF)

      tcga=tcga
     &     +(fx1p(j)*CqqPqq(z1)*flgq+fx1p(0)*(CqqPqg(z1)+CqgPgg(z1)))
     &     *(-dlog(xx10))*fx20(k)*msqc(j,k) 

c    C*gamma: second leg (ignore delta term in Cqq: taken into account with tH1stF)

      tcga=tcga
     &     +(fx2p(k)*CqqPqq(z2)*flgq+fx2p(0)*(CqqPqg(z2)+CqgPgg(z2)))
     &     *(-dlog(xx20))*fx10(j)*msqc(j,k) 

c    Add Cqg*Pgq contribution

      do l=1,nf
      tcga=tcga+(fx1p(l)+fx1p(-l))*CqgPgq(z1)
     &           *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq 
      tcga=tcga+(fx2p(l)+fx2p(-l))*CqgPgq(z2)
     &           *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq 
      enddo

CC  Start 2-loop AP

C   Gluon + pure singlet


      do l=-nf,nf
      if(l.eq.0) then
      tgamma2=tgamma2+fx1p(0)*P2qg(z1)
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)
      tgamma2=tgamma2+fx2p(0)*P2qg(z2)
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)
      else
      tgamma2=tgamma2+fx1p(l)*P2qqS(z1)
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq
      tgamma2=tgamma2+fx2p(l)*P2qqS(z2)
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq
      endif
      enddo


C   P2qq non-singlet: regular part

      tgamma2=tgamma2+fx1p(j)*P2qqV(z1)
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq
      tgamma2=tgamma2+fx2p(k)*P2qqV(z2)
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq


C   P2qq non-singlet: 1/(1-z)_+


      diff=-dlog(xx10)
     &  *(fx1p(j)-fx10(j)*xx10**beta)/(1-z1)
     &  - D0int(xx10)*fx10(j)      
  
      tgamma2=tgamma2+2d0/3*Kappa*diff*fx20(k)*msqc(j,k)*flgq


      diff=-dlog(xx20)
     &  *(fx2p(k)-fx20(k)*xx20**alfa)/(1-z2)
     &  - D0int(xx20)*fx20(k)      
  
      tgamma2=tgamma2+2d0/3*Kappa*diff*fx10(j)*msqc(j,k)*flgq

      

C   P2qqb non singlet

      tgamma2=tgamma2+fx1p(-j)*P2qqbV(z1)
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq

      tgamma2=tgamma2+fx2p(-k)*P2qqbV(z2)
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq


CCCCCCCCCCCC   End of NNLO scale dependence CCCCCCCCCCCCCCCCC


 75   continue

      enddo
      enddo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 72   xmsq=tdelta


      if(order.eq.1)then
       xmsq=xmsq+asopi*(tH1st+LF*tH1stF)
      elseif(order.eq.2)then
       xmsq=xmsq+asopi*(tH1st+LF*tH1stF)
     &          +asopi**2*(tdelta*H2qqdelta+tH2st)

CC     add scale dependence at NNLO

       xmsq=xmsq+asopi**2*(0.5d0*beta0*LF**2*tH1stF
     &                   +tgamma2*LF
     &                   -beta0*LR*(tH1st+LF*tH1stF)
     &                   +LF*tcga+0.5d0*LF**2*tgaga)


C     Include missing delta term from C*gamma (no factor 2 here !)

      xmsq=xmsq+asopi**2*(LF*C1qqdelta*tH1stF)

C     Include missing term from contact term in 2 loop AP

      xmsq=xmsq+asopi**2*(2*Delta2qq*tdelta)*LF

      endif     


      xmsq=flux*pswt*xmsq/BrnRat
      f(npdf+1)=xmsq
c      print *,npdf,xmsq
      call getptildejet(0,pjet)
      
c      call dotem(nvec,pjet,s)

      if (doFill.ne.0) then
         val=xmsq*wgt
          call hists_fill(p(3,:),p(4,:),val)
      endif

      enddo                     ! end PDF loop

      lowintHst_dynnlo=f(1)

      return

 999  continue
      
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function C2qqreg(z)
      implicit none
      real *8 Pi,Z2,Z3,myli2,myli3,z,CA,CF
      integer nf

      external myli2,myli3

      Pi=3.14159265358979d0
      z2=Pi**2/6

      Z3=1.20205690316d0

      CF=4d0/3
      CA=3d0
      nf=5

      C2qqreg=   
     & (CF*(-344+24*Pi**2+974*z-1600*CA*z+2052*CF*z+148*nf*z-60*Pi**2*z+
     & 108*CF*Pi**2*z-1188*z**2+1584*CA*z**2-4104*CF*z**2-72*nf*z**2+
     & 72*Pi**2*z**2-216*CF*Pi**2*z**2+830*z**3+16*CA*z**3+2052*CF*z**3-
     & 76*nf*z**3-60*Pi**2*z**3+108*CF*Pi**2*z**3-
     & 272*z**4+24*Pi**2*z**4+
     & 324*CA*z*z2-1728*CF*z*z2 - 648*CA*z**2*z2 + 3456*CF*z**2*z2+
     & 324*CA*z**3*z2-1728*CF*z**3*z2 + 1188*CA*z*Z3 + 864*CF*z*Z3-
     & 324*CA*z**3*Z3-864*CF*z**3*Z3 - 108*CA*z**2*dlog(1-z) +
     & 108*CF*z**2*dlog(1-z)+108*CA*z**3*dlog(1-z)-
     & 108*CF*z**3*dlog(1-z)-
     & 216*CA*z*z2*dlog(1-z) + 216*CF*z*z2*dlog(1-z) -
     & 216*CA*z**3*z2*dlog(1-z)+216*CF*z**3*z2*dlog(1-z)-252*z*dlog(z)+
     & 348*CA*z*dlog(z)-540*CF*z*dlog(z)-
     & 60*nf*z*dlog(z)+612*z**2*dlog(z)-
     & 432*CA*z**2*dlog(z)+1404*CF*z**2*dlog(z)-744*z**3*dlog(z)+
     & 996*CA*z**3*dlog(z)-1728*CF*z**3*dlog(z)-60*nf*z**3*dlog(z)+
     & 384*z**4*dlog(z)-144*dlog(1-z)*dlog(z)+360*z*dlog(1-z)*dlog(z)-
     & 216*CA*z*dlog(1-z)*dlog(z) + 648*CF*z*dlog(1-z)*dlog(z) -
     & 432*z**2*dlog(1-z)*dlog(z) + 432*CA*z**2*dlog(1-z)*dlog(z) -
     & 1296*CF*z**2*dlog(1-z)*dlog(z) + 360*z**3*dlog(1-z)*dlog(z) -
     & 216*CA*z**3*dlog(1-z)*dlog(z) + 648*CF*z**3*dlog(1-z)*dlog(z) -
     & 144*z**4*dlog(1-z)*dlog(z)+216*CA*z*dlog(1-z)**2*dlog(z) -
     & 324*CF*z*dlog(1-z)**2*dlog(z)+216*CA*z**3*dlog(1-z)**2*dlog(z) -
     & 324*CF*z**3*dlog(1-z)**2*dlog(z)+27*z*dlog(z)**2+
     & 99*CA*z*dlog(z)**2-
     & 162*CF*z*dlog(z)**2-18*nf*z*dlog(z)**2+108*CA*z**2*dlog(z)**2 -
     & 108*CF*z**2*dlog(z)**2+45*z**3*dlog(z)**2-9*CA*z**3*dlog(z)**2+
     & 108*CF*z**3*dlog(z)**2-18*nf*z**3*dlog(z)**2-72*z**4*dlog(z)**2-
     & 108*CF*z*dlog(1-z)*dlog(z)**2-108*CF*z**3*dlog(1-z)*dlog(z)**2-
     & 18*z*dlog(z)**3 + 18*CA*z*dlog(z)**3-
     & 18*CF*z*dlog(z)**3+18*z**3*dlog(z)**3 +
     & 18*CA*z**3*dlog(z)**3+18*CF*z**3*dlog(z)**3 -
     & 72*((-1 + z)**2*(2 + (-1+3*CA-6*CF)*z + 2*z**2) -
     & 3*(CA-CF)*z*(1 + z**2)*dlog(1-z)-3*(CA-3*CF)*z*(1+z**2)*dlog(z))*
     & myli2(z) + 216*CA*z*myli3(1-z)-216*CF*z*myli3(1-z) +
     & 216*CA*z**3*myli3(1-z)-216*CF*z**3*myli3(1-z) -
     & 432*CA*z*myli3(z) + 1080*CF*z*myli3(z) -
     & 432*CA*z**3*myli3(z)+1080*CF*z**3*myli3(z)-1944*CF*z*Z3-
     & 216*CF*z**3*Z3))/(864*(-1 + z)*z)




        return
        end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function C2qqp(z)
      implicit none
      real *8 Pi,myli2,z,CF
      integer nf

      external myli2

      Pi=3.14159265358979d0


      CF=4d0/3
     

      C2qqp=(CF*(2*(-1+z)*(-172+143*z-136*z**2+6*Pi**2*(2-z+2*z**2))-
     &   12*(z*(-21 + 30*z - 32*z**2)+
     &  6*(-2+3*z-3*z**2+2*z**3)*dlog(1-z))*
     &  dlog(z)-9*z*(3+3*z+8*z**2)*dlog(z)**2+18*z*(1 + z)*dlog(z)**3-
     &  72*(-2 + 3*z - 3*z**2 + 2*z**3)*myli2(z)))/(864d0*z)

 
      return
      end   

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function C2qqb(z)
      implicit none
      real *8 Pi,Z3,myli2,myli3,z,CF,CA,C2qqp
    

      external myli2,myli3,C2qqp

      Pi=3.14159265358979d0
      Z3=1.20205690316d0

      CF=4d0/3   
      CA=3d0


       C2qqb=C2qqp(z)+
     &  (CF*(-CA+2*CF)*(45-3*Pi**2-2*Pi**2*z-45*z**2+Pi**2*z**2+9*
     &  dlog(z)+42*z*dlog(z)+33*z**2*dlog(z)+12*dlog(1-z)*dlog(z)-
     &  12*z**2*dlog(1-z)*dlog(z)-dlog(z)**3-z**2*
     &  dlog(z)**3+2*Pi**2*dlog(1+z) +
     &  2*Pi**2*z**2*dlog(1+z)-12*dlog(z)*
     &  dlog(1+z)-24*z*dlog(z)*dlog(1+z)-
     &  12*z**2*dlog(z)*dlog(1+z)+6*dlog(z)**2*dlog(1+z)+
     &  6*z**2*dlog(z)**2*dlog(1+z)-4*dlog(1+z)**3-4*z**2*dlog(1+z)**3-
     &  12*((1+z)**2+(1+z**2)*dlog(z))*myli2(-z)-
     &  12*(-1+z**2+dlog(z)+z**2*dlog(z))*myli2(z)+36*myli3(-z)+
     &  36*z**2*myli3(-z)+24*myli3(z)+24*z**2*myli3(z)+
     &  24*myli3(1d0/(1+z))+24*z**2*myli3(1d0/(1+z))-
     &  18*Z3-18*z**2*Z3))/(48*(1+z))


      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     NEW version August 2013
C

      double precision function C2qg(z)
      implicit none
      real *8 z,CF,CA,Pi,Z3,myli2,myli3
      external myli2,myli3


      Pi=3.14159265358979d0
      Z3=1.20205690316d0
      CF=4d0/3
      CA=3d0


      C2qg=(688*CA-1260*CA*z-702*CF*z+1548*CA*z**2+2322*CF*z**2+ 
     -    72*CA*Pi**2*z**2+144*CF*Pi**2*z**2-1192*CA*z**3-2160*CF*z**3- 
     -    144*CF*Pi**2*z**3+324*CA*z**2*Log(1-z)-324*CF*z**2*Log(1-z)- 
     -    432*CA*z**3*Log(1-z)+432*CF*z**3*Log(1-z)- 
     -    216*CA*z**2*Log(1-z)**2+216*CF*z**2*Log(1-z)**2+ 
     -    216*CA*z**3*Log(1-z)**2-216*CF*z**3*Log(1-z)**2 + 
     -    36*CA*z*Log(1 - z)**3 - 36*CF*z*Log(1 - z)**3 - 
     -    72*CA*z**2*Log(1-z)**3 + 72*CF*z**2*Log(1-z)**3 + 
     -    72*CA*z**3*Log(1-z)**3-72*CF*z**3*Log(1-z)**3 + 
     -    504*CA*z*Log(z)+432*CF*z*Log(z)-720*CA*z**2*Log(z) + 
     -    810*CF*z**2*Log(z)+144*CA*Pi**2*z**2*Log(z)
     -    +1632*CA*z**3*Log(z)- 
     -    432*CF*z**3*Log(z) - 432*CF*z**2*Log(1-z)*Log(z) + 
     -    432*CF*z**3*Log(1-z)*Log(z)+108*CF*z*Log(1-z)**2*Log(z)- 
     -    216*CF*z**2*Log(1-z)**2*Log(z)+
     -    216*CF*z**3*Log(1-z)**2*Log(z)- 
     -    54*CA*z*Log(z)**2+27*CF*z*Log(z)**2+216*CA*z**2*Log(z)**2+ 
     -    324*CF*z**2*Log(z)**2 - 792*CA*z**3*Log(z)**2 - 
     -    216*CF*z**3*Log(z)**2 + 108*CF*z*Log(1 - z)*Log(z)**2- 
     -    864*CA*z**2*Log(1-z)*Log(z)**2-
     -    216*CF*z**2*Log(1-z)*Log(z)**2+ 
     -    216*CF*z**3*Log(1-z)*Log(z)**2+36*CA*z*Log(z)**3 - 
     -    18*CF*z*Log(z)**3+72*CA*z**2*Log(z)**3+36*CF*z**2*Log(z)**3- 
     -    72*CF*z**3*Log(z)**3 + 36*CA*Pi**2*z*Log(1 + z) + 
     -    72*CA*Pi**2*z**2*Log(1 + z) + 72*CA*Pi**2*z**3*Log(1+z)+ 
     -    432*CA*z**2*Log(z)*Log(1 + z) + 432*CA*z**3*Log(z)*Log(1+z)+ 
     -    108*CA*z*Log(z)**2*Log(1+z)+216*CA*z**2*Log(z)**2*Log(1+z)+ 
     -    216*CA*z**3*Log(z)**2*Log(1+z)-72*CA*z*Log(1+z)**3 - 
     -    144*CA*z**2*Log(1 + z)**3 - 144*CA*z**3*Log(1 + z)**3 - 
     -    72*(3*(CA - CF)*z*(1 - 2*z + 2*z**2)*Log(1 - z) + 
     -       2*CA*(2 - 3*z + 12*z**2 - 11*z**3 + 6*z**2*Log(z)))*
     -     myLi2(1-z) - 216*CA*z*
     -     (-2*z*(1 + z) + (1 + 2*z + 2*z**2)*Log(z))*myli2(-z)+ 
     -    216*CF*z*Log(z)*myli2(z)-
     -    1728*CA*z**2*Log(z)*myli2(z)- 
     -    432*CF*z**2*Log(z)*myli2(z)+432*CF*z**3*Log(z)*myli2(z)+ 
     -    216*CA*z*myli3(1-z)-216*CF*z*myli3(1-z)- 
     -    432*CA*z**2*myli3(1-z)+432*CF*z**2*myli3(1-z) + 
     -    432*CA*z**3*myli3(1-z) - 432*CF*z**3*myli3(1-z) + 
     -    648*CA*z*myli3(-z) + 1296*CA*z**2*myli3(-z)+ 
     -    1296*CA*z**3*myli3(-z) - 216*CF*z*myli3(z) + 
     -    1728*CA*z**2*myli3(z) + 432*CF*z**2*myli3(z)- 
     -    432*CF*z**3*myli3(z) + 432*CA*z*myli3(1/(1+z))+ 
     -    864*CA*z**2*myli3(1/(1+z)) +
     -    864*CA*z**3*myli3(1/(1+z)) - 
     -    648*CA*z*Z3 + 1728*CF*z*Z3-3456*CF*z**2*Z3 - 
     -    1296*CA*z**3*Z3 + 3456*CF*z**3*Z3)/(1728d0*z)


      return
      end


    
