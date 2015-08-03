C     My phase space for W(Z)+2partons
C     To be used for double real contribution

      subroutine gen4(r,p,wt4,*)
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'process.f'
      include 'phasemin.f'
      integer nu
      double precision r(mxdim)
      double precision wt4,p1(4),p2(4),p3(4),p4(4),p5(4),p6(4)
      double precision p(mxpart,4)
      double precision pswt,xjac,p1ext(4),p2ext(4)
      double precision xx(2),tau,y,sqrts
      common/pext/p1ext,p2ext
      common/x1x2/xx
      common/energy/sqrts


      wt4=0d0

      tau=dexp(dlog(taumin)*r(9))
      y=0.5d0*dlog(tau)*(1d0-2d0*r(10))
      xjac=dlog(taumin)*tau*dlog(tau)

      xx(1)=dsqrt(tau)*dexp(+y)
      xx(2)=dsqrt(tau)*dexp(-y)    


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
   
      call phase4(r,p1,p2,p3,p4,p5,p6,pswt,*999) 


      do nu=1,4
      p(1,nu)=p1(nu)
      p(2,nu)=p2(nu)
      p(3,nu)=p3(nu)
      p(4,nu)=p4(nu)
      p(5,nu)=p5(nu)
      p(6,nu)=p6(nu)
      p(7,nu)=0d0

      enddo 

      wt4=xjac*pswt
      
      return

 999  return 1
      end

