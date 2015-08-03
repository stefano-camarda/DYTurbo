      subroutine genBORN2qT(q2,qt2,shat,r,p,wt,*)
c----generate phase space weight and vectors p(i,4) for i=1,2,3,4
c----and x1 and x2 given seven random numbers and q2
c----all other four momenta must be zero
C GF  07/12: Included various prescriptions to absorb (in the counterterm/asymptotic) the recoil due to soft gluon emissions
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'phasemin.f'
      integer nu

      double precision r(mxdim),sqrts,wt,wt34,
     . p(mxpart,4),p1(4),p2(4),p3(4),p4(4),q(4)
      double precision pswt,xjac,xx(2),tau,tau0,y,q2,yq,shat
      double precision kt1,kt2,qP1,qP2,zeta1,mt,qt2,qt,phi

      common/energy/sqrts
      common/xx0/xx

!
      mt=dsqrt(q2+qt2)
      qt=dsqrt(qt2)
!
      wt=0d0

      tau=shat/sqrts**2
      tau0=q2/sqrts**2
      y=0.5d0*dlog(tau0)*(1d0-2d0*r(7))
      xjac=dabs(dlog(tau0))

      xx(1)=dsqrt(tau0)*dexp(+y)
      xx(2)=dsqrt(tau0)*dexp(-y)

c---if x's out of normal range alternative return
      if   ((xx(1) .gt. 1d0) 
     & .or. (xx(2) .gt. 1d0)
     & .or. (xx(1) .lt. xmin)
     & .or. (xx(2) .lt. xmin)) return 99

!      p1(4)=-xx(1)*sqrts*half
!      p1(1)=zip
!      p1(2)=zip
!      p1(3)=-xx(1)*sqrts*half

!      p2(4)=-xx(2)*sqrts*half
!      p2(1)=zip
!      p2(2)=zip
!      p2(3)=+xx(2)*sqrts*half

!      do nu=1,4
!      q(nu)=-p1(nu)-p2(nu)
!      enddo

      phi=twopi*r(1)

      mt=dsqrt(q2+qt2)
      qt=dsqrt(qt2)

      q(1)=qt*dcos(phi)
      q(2)=qt*dsin(phi)
      q(3)=0.5d0*mt*(dexp(y)-dexp(-y))
      q(4)=0.5d0*mt*(dexp(y)+dexp(-y))

       kt1=q(1)/2d0   !
       kt2=q(2)/2d0   !  CS FRAME PRESCRIPTION

       zeta1=1d0/q2/2d0*(q2+2d0*(q(1)*kt1+q(2)*kt2)
     & +dsqrt((q2+2d0*(q(1)*kt1+q(2)*kt2))**2
     & -4d0*mt**2*(kt1**2+kt2**2)))

       qP1=(q(4)-q(3))*sqrts/2d0
       qP2=(q(4)+q(3))*sqrts/2d0


       p(1,4)=sqrts/2d0*(zeta1*q2/2d0/qP1
     & +(kt1**2+kt2**2)/zeta1*qP1/q2/sqrts**2*2d0)
       p(1,1)=kt1
       p(1,2)=kt2
       p(1,3)=sqrts/2d0*(zeta1*q2/2d0/qP1
     & -(kt1**2+kt2**2)/zeta1*qP1/q2/sqrts**2*2d0)

       p(2,4)=q(4)-p(1,4)
       p(2,1)=q(1)-p(1,1)
       p(2,2)=q(2)-p(1,2)
       p(2,3)=q(3)-p(1,3)

      call phi3m0(r(4),r(5),q,p3,p4,wt34,*99)      

      do nu=1,4
!      p(1,nu)=p1(nu)
!      p(2,nu)=p2(nu)
      p(3,nu)=p3(nu)
      p(4,nu)=p4(nu)
      enddo 
      wt=xjac*wt34
      return
 99   continue
      wt=0d0

      return
      end
