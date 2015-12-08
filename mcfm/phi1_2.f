      subroutine phi1_2(x1,x2,x3,x4,p1,p2,p3,wt,*)
c     massive particle p1 decaying into p2 mass m2 and p3 mass m3.
c     with invariant mass 
c     of particle two s2 and particle three s3 integrated over.
c     vectors returned p2 and p3 are in the same frame as p1 is supplied
c     Expression evaluate is 
c     ds2 ds3 d^4 p2 d^4 p3 (2 pi)^4 delta(p1-p2-p3)/(2 pi)^6 
c     delta(p2^2-s2) delta(p3^2-s3)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'process.f'
      include 'zerowidth.f'
      double precision p1(4),p2(4),p3(4),p3cm(4)
      double precision x1,x2,x3,x4,costh,sinth,phi,cphi,sphi
      double precision wt,wt0,w2,w3
      double precision s2max,s2min,s3max,s3min
      double precision m1,m2,s1,s2,s3,lambda,mass2,width2,mass3,width3
      integer j,n2,n3
      common/breit/n2,n3,mass2,width2,mass3,width3
      common/lambda/lambda,s1,s2,s3
      parameter(wt0=one/8d0/pi)
    

      wt=0d0
      s1=p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2  
      if (s1 .lt. 0d0) return 1
      m1=dsqrt(s1)
      if (
     . zerowidth 
     . .and. (m1 .lt. mass2*dfloat(n2)+mass3*dfloat(n3))
     .    ) return 1
c      s2min=bbsqmin
c      s2max=min(s1,bbsqmax)
      s2min=0d0
      s2max=s1
      if (s2min .gt. s2max) return 1
      if (n2 .eq. 0) then
         w2=s2max-s2min
         s2=s2max*x1+s2min*(1d0-x1)
      elseif (n2 .eq. 1) then
         call breitw(x1,s2min,s2max,mass2,width2,s2,w2)       
      endif

      m2=dsqrt(s2)
      s3min=1d-15
      s3max=(m2-m1)**2
      if (s3max-s3min .lt. 1d-12) return 1
      if (n3 .eq. 0) then
         w3=s3max-s3min
         s3=s3max*x2+s3min*(1d0-x2)
      elseif (n3 .eq. 1) then
         call breitw(x2,s3min,s3max,mass3,width3,s3,w3)       
      endif


      costh=two*x3-one      
      phi=twopi*x4
      sinth=dsqrt(one-costh**2)
      cphi=dcos(phi)
      sphi=dsin(phi)
      lambda=((s1-s2-s3)**2-4d0*s2*s3)

      if (lambda .lt. 0d0) then
c      write(6,*) '(lambda .lt. 0) in phi1_2.f',lambda
c      write(6,*) 'sqrt(s1)',sqrt(s1)
c      write(6,*) 'sqrt(s2)',sqrt(s2)
c      write(6,*) 'sqrt(s3)',sqrt(s3)
c      write(6,*) s3min,s3,s3max,m1,m2,sqrt(s1),sqrt(s2)
      return 1
      endif
      lambda=dsqrt(lambda)
      wt=wt0*w2*w3*lambda/s1


      p3cm(4)=m1/two*(s1+s3-s2)/s1
      p3cm(1)=m1/two*lambda/s1*sinth*sphi
      p3cm(2)=m1/two*lambda/s1*sinth*cphi
      p3cm(3)=m1/two*lambda/s1*costh
      call boost(m1,p1,p3cm,p3)
      do j=1,4
      p2(j)=p1(j)-p3(j)
      enddo
      if (  (p1(4) .lt. 0d0) 
     & .or. (p2(4) .lt. 0d0) 
     & .or. (p3(4) .lt. 0d0)) then 
       if (case(1:5) .ne. 'vlchk') then 
        write(6,*) 'p1',p1(4),p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2,s1
        write(6,*) 'p2',p2(4),p2(4)**2-p2(1)**2-p2(2)**2-p2(3)**2,s2
        write(6,*) 'p3',p3(4),p3(4)**2-p3(1)**2-p3(2)**2-p3(3)**2,s3
        write(6,*) n2,n3
       endif
       return 1
      endif

      return
      end



