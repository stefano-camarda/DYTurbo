      subroutine gen3(r,p,wt3,*)
c----generate 3 dimensional phase space weight and vectors p(7,4)
c----and x1 and x2 given seven random numbers
c----p(5,i) and p(4,i) are set equal to zero
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'phasemin.f'
      integer nu

      double precision r(mxdim),sqrts,wt3,
     . p(mxpart,4),p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4)
      double precision pswt,xjac,xx(2),tau,y

      common/energy/sqrts
      common/x1x2/xx
      data p6/0d0,0d0,0d0,0d0/
      data p7/0d0,0d0,0d0,0d0/

      wt3=0d0
      tau=dexp(dlog(taumin)*r(6))
      y=0.5d0*dlog(tau)*(1d0-2d0*r(7))
      xjac=dlog(taumin)*tau*dlog(tau)

      xx(1)=dsqrt(tau)*dexp(+y)
      xx(2)=dsqrt(tau)*dexp(-y)
 
c---if x's out of normal range alternative return
      if   ((xx(1) .gt. 1d0) 
     & .or. (xx(2) .gt. 1d0)
     & .or. (xx(1) .lt. xmin)
     & .or. (xx(2) .lt. xmin)) return 1

      p1(4)=-xx(1)*sqrts*half
      p1(1)=zip
      p1(2)=zip
      p1(3)=-xx(1)*sqrts*half

      p2(4)=-xx(2)*sqrts*half
      p2(1)=zip
      p2(2)=zip
      p2(3)=+xx(2)*sqrts*half


      call phase3(r,p1,p2,p3,p4,p5,p6,p7,pswt)

      do nu=1,4
      p(1,nu)=p1(nu)
      p(2,nu)=p2(nu)
      p(3,nu)=p3(nu)
      p(4,nu)=p4(nu)
      p(5,nu)=p5(nu)
      p(6,nu)=p6(nu)
      p(7,nu)=p7(nu)
      enddo 
      wt3=xjac*pswt


      if(wt3 .eq. 0d0) return 1

      return
      end
