C     Phase space generation for V + 2 partons
C     To be used for V+jet real contribution

      subroutine dygen4(r,p,wt4,*)
      implicit none
      include 'constants.f'
      include 'mxdim.f'
c      include 'process.f'
c      include 'masses.f'
      include 'phasemin.f'
      integer nu
      double precision r(mxdim)
      double precision wt4,p1(4),p2(4),p3(4),p4(4),p5(4),p6(4)
      double precision p12(4),p34(4),p56(4)
      double precision p(mxpart,4)
      double precision pswt,xjac
      double precision xx(2),tau,logtau,y,sqrts
      double precision expy
      common/x1x2/xx
      common/energy/sqrts

      double precision wt,wt3456,wt34,wt56,wt0
      integer j
      parameter(wt0=1d0/twopi**2)

      wt4=0d0

      logtau=logtaumin*r(9)
      tau=dexp(logtau)
      y=0.5d0*logtau*(1d0-2d0*r(10))
      xjac=logtaumin*tau*logtau

c      xx(1)=dsqrt(tau)*dexp(+y)
c      xx(2)=dsqrt(tau)*dexp(-y)    
      expy=exp(+y)
      xx(1)=sqrt(tau)*expy
      xx(2)=sqrt(tau)/expy


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
   
c---- generate phase space for 2-->4 process
c---- r(mxdim),p1(4),p2(4) are inputs reversed in sign 
c---- from physical values 
c---- phase space for -p1-p2 --> p3+p4+p5+p6
c---- with all 2 pi's (ie 1/(2*pi)^8)

      do j=1,4
      p12(j)=-p1(j)-p2(j)
      enddo

c     generate the dilepton (p34) and dijet (p56) systems
      call dyphi1_2(r(1),r(2),r(3),r(4),p12,p56,p34,wt3456,*999)
c      call phi1_2(r(1),r(2),r(3),r(4),p12,p56,p34,wt3456,*999)

c     generate the dilepton phase space
      call phi3m0(r(7),r(8),p34,p3,p4,wt34,*999)

c     generate the dijet phase space
      call phi3m0(r(5),r(6),p56,p5,p6,wt56,*999)
      
      pswt=wt0*wt3456*wt34*wt56
      
      do nu=1,4
      p(1,nu)=p1(nu)
      p(2,nu)=p2(nu)
      p(3,nu)=p3(nu)
      p(4,nu)=p4(nu)
      p(5,nu)=p5(nu)
      p(6,nu)=p6(nu)
      enddo 

      wt4=xjac*pswt
      
      return

 999  return 1
      end


      subroutine dyphi1_2(x1,x2,x3,x4,p1,p2,p3,wt,*)
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
      include 'zerowidth.f'
      include 'limits.f'
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
      
c     s1 is the partonic centre-of-mass energy (mass^2 of p12)
      s1=p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2  
      if (s1 .lt. 0d0) return 1
      m1=sqrt(s1)
      if (zerowidth.and.(m1.lt.mass2*dfloat(n2)+mass3*dfloat(n3))) then
         return 1
      endif
      
c     s2 is the mass^2 of the dijet system (p56)
      s2min=0d0
c     s2max=s1                   ! Allow sqrt(s2) up to the mass of p12
      s2max=(m1-sqrt(wsqmin))**2 ! Allow sqrt(s2) up to the mass of p12 minus the minimum mass of p34
      if (s2min .gt. s2max) return 1
      w2=s2max-s2min
      s2=s2max*x1+s2min*(1d0-x1)
      m2=sqrt(s2)
      
c     s3 is the mass^2 of the dilepton system (p34)
      s3min=wsqmin
      s3max=min(wsqmax,(m1-m2)**2)
      if (s3max-s3min .lt. 1d-12) return 1
      call breitw(x2,s3min,s3max,mass3,width3,s3,w3)       

c     angular variables in the p34+p56 rest frame
      costh=two*x3-one      
      phi=twopi*x4
      sinth=sqrt(one-costh**2)
      cphi=cos(phi)
      sphi=sin(phi)
c      sphi=sqrt(max(0,1-cphi**2))
c      if (phi.lt.0) sphi=-sphi
      
      lambda=((s1-s2-s3)**2-4d0*s2*s3)

      if (lambda .lt. 0d0) then
c      write(6,*) '(lambda .lt. 0) in phi1_2.f',lambda
c      write(6,*) 'sqrt(s1)',sqrt(s1)
c      write(6,*) 'sqrt(s2)',sqrt(s2)
c      write(6,*) 'sqrt(s3)',sqrt(s3)
c      write(6,*) s3min,s3,s3max,m1,m2,sqrt(s1),sqrt(s2)
      return 1
      endif
      lambda=sqrt(lambda)
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
         return 1
      endif

      return
      end
      
