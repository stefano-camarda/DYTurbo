C     Lepton isolation for W and Z production

      logical function isolation(p)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),Ris,pt,eps,Ehad
      double precision Et5,Et6,R
      integer i,i1,i2
      common/isolabel/i1,i2

      isolation=.true.


c     Leptons are isolated if total transverse energy in a cone
c     of radius Ris is smaller than eps*pt

      Ris=0.4d0
      eps=0.1d0

      Et5=pt(5,p)
      Et6=pt(6,p)

      
C     Loop over leptons

      do i=i1,i2
      
      Ehad=0d0

C     Energy around lepton i

      if(Et6.eq.0d0) then
       if(r(p,i,5).lt.Ris)Ehad=Et5
      else
       if((r(p,i,5).lt.Ris).and.(r(p,i,6).lt.Ris)) then
        Ehad=Et5+Et6
       elseif ((r(p,i,5).lt.Ris).and.(r(p,i,6).gt.Ris)) then
        Ehad=Et5
       elseif ((r(p,i,5).gt.Ris).and.(r(p,i,6).lt.Ris)) then
        Ehad=Et6
       endif
      endif


      if(Ehad.gt.eps*pt(i,p))then
       isolation=.false.
       return
      endif

      enddo


      return
      end
