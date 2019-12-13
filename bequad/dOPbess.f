c
c
      program dOPbess
c
c
      external dquad
      double precision dnu,x0,a1(100),b1(100),a2(100),
     *b2(100),a3(100),b3(100),a4(100),b4(100),e(100),
     *d1mach,depsma,dendl(4),dendr(4),deps,dxp(1),dyp(1),
     *dxfer(100),dwfer(100),dalpha(40),dbeta(40),dbe(40),
     *dx(100),dw(100),dxm(400),dwm(400),dp0(400),dp1(400),dp2(400)
      logical finld,finrd
      common/d/dnu,x0,a1,b1,a2,b2,a3,b3,a4,b4,e
c
c This routine generates in double precision the first 40
c recurrence coefficients of the orthogonal polynomials for 
c the K-Bessel weight function. Printed are the double-precision
c values of the alpha- and beta-coefficients. The routine
c requires the ORTHPOL routines d1mach, drecur, dgauss, dmcdis.
c dqgp, dsti, and dlancz (see W. Gautschi, "Algorithm 726: 
c ORTHPOL -- a package of routines for generating orthogonal 
c polynomials and Gauss-type quadrature rules", ACM Trans. Math. 
c Software, v. 20 (1994), pp. 21-62).
c
      print *,'hello'
      stop
      nmax=40
      dnu=1.d0/3.d0
      if(dnu.lt..05d0 .or. dnu.gt..95d0) then
        write(*,1)
    1   format(1x,'nu out of range')
        stop
      end if
      x0=1.d0
      call drecur(100,6,0.d0,-dnu,a1,b1,ier1)
      call drecur(100,6,0.d0,dnu,a2,b2,ier2)
      call drecur(100,1,0.d0,0.d0,a3,b3,ier3)
      call drecur(100,7,0.d0,0.d0,a4,b4,ier4)
      if(ier1.ne.0 .or. ier2.ne.0 .or. ier3.ne.0 .or.
     *  ier4.ne.0) then
        write(*,2) ier1,ier2,ier3,ier4
    2     format(1x,'ier1 in drecur=',i5,'  ier2=',i5,
     *    'ier3=',i5,'  ier4=',i5)
        stop
      end if
      finld=.true.
      finrd=.false.
      depsma=d1mach(3)
c
c depsma is the machine double precision.
c
      iq=1
      idelta=1
      irout=1
      mcd=4
      mp=0
      ncpmd=100
c
c Set up the partition for the discretization of the inner product.
c
      dendl(1)=0.d0
      dendl(2)=1.d0
      dendl(3)=2.d0
      dendl(4)=3.d0
      dendr(1)=1.d0
      dendr(2)=2.d0
      dendr(3)=3.d0
      deps=1.d3*depsma
c
c Compute the desired recursion coefficients by the multiple-component
c discretization procedure. 
c
      call dmcdis(nmax,ncpmd,mcd,mp,dxp,dyp,dquad,deps,iq,idelta,irout,
     *  finld,finrd,dendl,dendr,dxfer,dwfer,dalpha,dbeta,ncapd,kountd,
     *  ierrd,ied,dbe,dx,dw,dxm,dwm,dp0,dp1,dp2)
c
c Print the results.
c
      write(*,3)
    3 format(/5x,'k',6x,'dalpha(k)',13x,'dbeta(k)'/)
      do 20 k=1,nmax
        km1=k-1
        write(*,4) km1,dalpha(k),dbeta(k)
    4   format(1x,i5,2d21.13)
   20 continue
      stop
      end 

      subroutine dquad(n,dx,dw,i,ierr)
      double precision dx(n),dw(n),dnu,x0,a1(100),b1(100),a2(100),
     *b2(100),a3(100),b3(100),a4(100),b4(100),e(100),z(23),w(23),
     *pi,c,eps,d1mach,x1,c1,dgamma,x,x2,sum1,t1,dm,sum,t,dbessK
      common/d/dnu,x0,a1,b1,a2,b2,a3,b3,a4,b4,e
      common/da/z,w
      pi=4.d0*datan(1.d0)
      c=2.d0*dcos(.5d0*dnu*pi)/pi
      eps=d1mach(3)
      x1=10.d0
      if(i.eq.1) then
        c1=pi*((.25d0*x0)**(1.d0-dnu))/(dsin(dnu*pi)*dgamma(1.d0-
     *    dnu,ie))
        call dgauss(n,a1,b1,eps,dx,dw,ierg,e)
        do 20 k=1,n
          dx(k)=.5d0*x0*(1.d0+dx(k))
          x=dx(k)
          m=0
          x2=x**2
          sum1=1.d0
          t1=1.d0
   10     m=m+1
          dm=dble(m)
          sum=sum1
          t=t1
          t1=x2*t/(4.d0*dm*(dm-dnu))
          sum1=sum+t1
          if(sum1.ne.sum) goto 10
          dw(k)=c*c1*sum1*dw(k)
   20   continue
      else if(i.eq.2) then
        c1=pi*((.25d0*x0)**(1.d0+dnu))/(dsin(dnu*pi)*dgamma(1.d0+
     *    dnu,ie))
        call dgauss(n,a2,b2,eps,dx,dw,ierg,e)
        do 40 k=1,n
          dx(k)=.5d0*x0*(1.d0+dx(k))
          x=dx(k)
          m=0
          x2=x**2
          sum1=1.d0
          t1=1.d0
   30     m=m+1
          dm=dble(m)
          sum=sum1
          t=t1
          t1=x2*t/(4.d0*dm*(dm+dnu))
          sum1=sum+t1
          if(sum1.ne.sum) goto 30
          dw(k)=-c*c1*sum1*dw(k)
   40   continue
      else if(i.eq.3) then
        call dgauss(n,a3,b3,eps,dx,dw,ierg,e)
        do 50 k=1,n
          dx(k)=.5d0*((x1-x0)*dx(k)+x1+x0)
          dw(k)=.5d0*c*(x1-x0)*dexp(-dx(k))*dbessK(dx(k),dnu)*dw(k)
   50   continue
      else 
        c1=dexp(-x1)
        call dgauss(n,a4,b4,eps,dx,dw,ierg,e)
        do 60 k=1,n
          dx(k)=x1+dx(k)
          dw(k)=c*c1*dbessK(dx(k),dnu)*dw(k)
   60   continue
      end if
      return
      end

      double precision function dbessK(x,dnu)
c
c This computes exp(x)K_nu(x) for positive x.
c
      double precision x,dnu,z(23),w(23),eps,d1mach,al,a(23),
     *b(23),e(23),pi,c1,dser,dgamma,c2,sum
      common/da/z,w
      if(x.lt.0.d0) then
        write(*,1)
    1   format(1x,'invalid x-argument in dbessK')
        stop
      end if
      eps=d1mach(3)
      al=dnu-.5d0
      call drecur(23,7,al,0.d0,a,b,ierr)
      call dgauss(23,a,b,eps,z,w,ierg,e)
      pi=4.d0*datan(1.d0)
      c1=pi/(2.d0*dsin(dnu*pi))
      if(x.le.2.d0) then
        dbessK=c1*dexp(x)*(((.5d0*x)**(-dnu))*dser(x,-dnu)/
     *    dgamma(1.d0-dnu,ie)-((.5d0*x)**dnu)*dser(x,dnu)/
     *    dgamma(1.d0+dnu,ie))
      else
        n=23
        c2=dsqrt(pi)/((2.d0**dnu)*dgamma(dnu+.5d0,ie))
        sum=0.d0
        do 10 k=1,n
          sum=sum+w(k)*(2.d0+z(k)/x)**al
   10   continue
        dbessK=c2*sum/dsqrt(x)
      end if
      return
      end

      double precision function dser(x,dnu)
      double precision x,dnu,x2,sum1,t1,dk,sum,t
      k=0
      x2=x**2
      sum1=1.d0
      t1=1.d0
   10 k=k+1
      dk=dble(k)
      sum=sum1
      t=t1
      t1=x2*t/(4.d0*dk*(dnu+dk))
      sum1=sum+t1
      if(sum1.ne.sum) goto 10
      dser=sum1
      return
      end
  
      double precision function dwf(x,i)
      double precision x
      return
      end

