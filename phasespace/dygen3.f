      subroutine dygen3(r,p,wt3,*)
c----generate 3 dimensional phase space weight and vectors p(7,4)
c----and x1 and x2 given seven random numbers
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'phasemin.f'

      integer j,nu
      double precision p12(4),p34(4),smin
      double precision wt125,wt34,wt0,m5
      parameter(wt0=1d0/twopi)

      double precision r(mxdim),sqrts,wt3,
     . p(mxpart,4),p1(4),p2(4),p3(4),p4(4),p5(4)
      double precision pswt,xjac,xx(2),tau,logtau,y,expy

      common/energy/sqrts
      common/x1x2/xx

      wt3=0d0
      logtau=logtaumin*r(6)
      tau=exp(logtau)
      y=0.5d0*logtau*(1d0-2d0*r(7))
      xjac=logtaumin*tau*logtau

c     xx(1)=dsqrt(tau)*dexp(+y)
c     xx(2)=dsqrt(tau)*dexp(-y)
      expy=exp(+y)
      xx(1)=sqrt(tau)*expy
      xx(2)=sqrt(tau)/expy
 
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


c----generate phase space for 2-->3 process
c----r(mxdim),p1(4),p2(4) are inputs 
c----incoming p1 and p2 reversed in sign from physical values 
c----i.e. phase space for -p1-p2 --> p3+p4+p5
c----with all 2 pi's (ie 1/(2*pi)^5)
c----(p4,p5) are dummies

c---m5 is the mass of p5
      m5=0d0

      do j=1,4
      p12(j)=-p1(j)-p2(j)
      enddo

c---generate p5 and p34, 
      call dyphi1_2m(m5,r(1),r(2),r(3),p12,p5,p34,wt125,*99)
c---decay 34-system
      call phi3m0(r(4),r(5),p34,p3,p4,wt34,*99)

      pswt=wt0*wt125*wt34
      
      do nu=1,4
      p(1,nu)=p1(nu)
      p(2,nu)=p2(nu)
      p(3,nu)=p3(nu)
      p(4,nu)=p4(nu)
      p(5,nu)=p5(nu)
      enddo 
      wt3=xjac*pswt


      if(wt3 .eq. 0d0) return 1

      return

 99   continue
      wt3=0d0
      return
      end

      subroutine dyphi1_2m(m2,x3,xth,xphi,p1,p2,p3,wt,*)
c     massive particle p1 decaying into p2 mass m2 and p3 mass-squared s3.
c     with invariant mass of particle three s3 integrated over.
c     s3min is the minimum value of s3.
c     Vectors returned p2 and p3 are in the same frame as p1 is supplied.
c     Expression evaluated is 
c     ds3 d^4 p2 d^4 p3 (2 pi)^4 delta(p1-p2-p3)/(2 pi)^6
c     delta(p2^2-m2) delta(p3^2-s3)
      implicit none
      include 'constants.f'
      include 'limits.f'
      double precision p1(4),p2(4),p3(4),p3cm(4)
      double precision x3,xth,xphi,costh,sinth,phi,cphi,sphi
      double precision wt,wt0,w3
      double precision s3max,s3min,xexp
      double precision m1,m2,m3,s1,s2,s3,lambda,xjac,rtxth,
     . mass2,width2,mass3,width3

      integer j,n2,n3
      integer jbranch

      common/breit/n2,n3,mass2,width2,mass3,width3
      parameter(wt0=one/8d0/pi)
      data jbranch/1/
      data xjac,rtxth/1d0,1d0/

      save jbranch,xexp,rtxth,xjac

      xexp=1d0
      wt=0d0
      s1=p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2  
      if (s1 .lt. 0d0) return 1
      m1=dsqrt(s1)
      s2=m2**2
      
c     s3 is the mass^2 of the dilepton system (p34)
      s3min=wsqmin
      s3max=min(wsqmax,(m1-m2)**2)
      if (s3min .gt. s3max) return 1
      call breitw(x3,s3min,s3max,mass3,width3,s3,w3) 

      m3=dsqrt(s3)
      if (m1-m2-m3.lt. 0d0) return 1


      if (jbranch .eq. 1) then
c S.C. to swicth off jbranch switching comment the following line
c        jbranch=2
        rtxth=xth**xexp 
        xjac=1d0/(xexp*xth**(xexp-1d0))
      elseif (jbranch .eq. 2) then 
        jbranch=1
        rtxth=1d0-xth**xexp
        xjac=1d0/(xexp*xth**(xexp-1d0))
      endif


      costh=two*rtxth-one      
      phi=twopi*xphi
      sinth=dsqrt(one-costh**2)
      cphi=dcos(phi)
      sphi=dsin(phi)
      lambda=((s1-s2-s3)**2-4d0*s2*s3)
      if (m1-m2-m3.lt. 0d0) then
      write(6,*) 'lambda in phi1_2m',lambda
      write(6,*) 's1 in phi1_2m',s1
      write(6,*) 's2 in phi1_2m',s2
      write(6,*) 's3 in phi1_2m',s3
      write(6,*) 'm1 in phi1_2m',m1
      write(6,*) 'm2 in phi1_2m',m2
      write(6,*) 'm3 in phi1_2m',m3
      write(6,*) 'm1-m2-m3 in phi1_2m',m1-m2-m3
      write(6,*) 'x3 in phi1_2m',x3
      write(6,*) 'n3 in phi1_2m',n3
      write(6,*) 'mass3 in phi1_2m',mass3
      return 1
      endif
      lambda=dsqrt(lambda)

      wt=wt0*w3*lambda/s1/xjac

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
