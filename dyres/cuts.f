      logical function cutsOLD(pjet,njets)
      logical cuts
      implicit none
      include 'constants.f'
      include 'masses.f'
      integer i,j,k,njets
      double precision pjet(mxpart,4),etvec(4)
      double precision pt,etarap
CC
      double precision pt3,pt4,eta3,eta4,pt34,y34,m34,yraptwo
      double precision pt5,eta5,eta6,ptjetmax,m45
      double precision cosphi34,deltaphi,ptmin,pt3dpt4,tmass
      double precision pte,ptmiss,etae
      double precision ppt(1:4),pto(1:4),mz(1:4),dmz(1:4),tmp(1:2)

      integer l(4),m(4),i1,i2,i3,i4

c      logical isol
c      common/isol/isol

      integer nproc
      common/nproc/nproc

CC
      cuts=.false.



CC    Insert here cuts

      pt3=dsqrt(pjet(3,1)**2+pjet(3,2)**2)
      pt4=dsqrt(pjet(4,1)**2+pjet(4,2)**2)      

      if(pt3.ne.0d0) then
       eta3=etarap(3,pjet)
       else
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      eta3=100d0
      endif

      if(pt4.ne.0d0) then
       eta4=etarap(4,pjet)
       else
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      eta4=100d0
      endif

      eta3=etarap(3,pjet)
      eta4=etarap(4,pjet)

      ptmin=min(pt3,pt4)

      pt34=dsqrt((pjet(3,1)+pjet(4,1))**2+(pjet(3,2)+pjet(4,2))**2)

      m34=dsqrt((pjet(3,4)+pjet(4,4))**2-(pjet(3,1)+pjet(4,1))**2
     &     -(pjet(3,2)+pjet(4,2))**2-(pjet(3,3)+pjet(4,3))**2)



      if(pt4.ne.0d0.and.pt3.ne.0d0) then
      y34=yraptwo(3,4,pjet)
       else
C-- set to 100 if this is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
      y34=100d0
      endif

C     Transverse mass

      pt3dpt4=pjet(3,1)*pjet(4,1)+pjet(3,2)*pjet(4,2)

      if(2d0*(pt3*pt4-pt3dpt4).le.0d0) then
!CM       write(*,*) 'WARNING: tmass^2=',2d0*(pt3*pt4-pt3dpt4)
       tmass=0d0
      else
      tmass=dsqrt(2d0*(pt3*pt4-pt3dpt4))
      endif
  
C     Cuts for Z production

      if(nproc.eq.3) then
       if(pt3.lt.20d0) cuts=.true.
       if(pt4.lt.20d0) cuts=.true.
       if(dabs(eta3).gt.2.4d0) cuts=.true.
       if(dabs(eta4).gt.2.4d0) cuts=.true.
       if(m34.lt.66d0.or.m34.gt.116d0) cuts=.true.
!       if(pt34.gt.600d0) cuts=.true.
      endif
C     Cuts for W production

      if(nproc.eq.1) then
       pte=pt4
       etae=eta4
       ptmiss=pt3
      elseif(nproc.eq.2) then
       pte=pt3
       etae=eta3
       ptmiss=pt4
      endif

      if((nproc.eq.1).or.(nproc.eq.2)) then
!       if(pte.lt.30d0) cuts=.true.
!       if(ptmiss.lt.30d0) cuts=.true.
!       if(dabs(etae).gt.2.4d0) cuts=.true.
!       if(tmass.lt.60d0) cuts=.true.
!       if(pt34.gt.30d0) cuts=.true.
       cuts=.false.
      endif


!!!!!!!! NO CUTS !!!!!!!!
!      cuts=.false.
!!!!!!!!

      return     
      end
 
 
 
 
