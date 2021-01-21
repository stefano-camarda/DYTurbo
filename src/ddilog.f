!c     li2 function called from vjet/utils.f
!      double precision function Li2(x)
!      implicit none
!      double precision x,dli2
!      external dli2
!      Li2 = dli2(x)
!      return
!      end
!      
!c     ddilog function called from mcfm/i3m.f, mcfm/lfunctions.f
!      double precision function ddilog(x)
!      implicit none
!      double precision x,dli2
!      external dli2
!      ddilog = dli2(x)
!      return
!      end
!
