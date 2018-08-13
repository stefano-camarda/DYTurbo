c     Defines the functions h1(z,v) and h2(z,v) according to
c     the definition in eq.(14) of hep-ph/0002078

      function h1(z)
      implicit none
      complex *16 h1,z,zz
      real *8 t,v,acc,tmp1,tmp2
      real *8 err,aa
      integer i,j
      common/ij/i,j      
!$OMP THREADPRIVATE(/ij/)
      common/v/v
c !$OMP THREADPRIVATE(/v/)
      common/c/zz
!$OMP THREADPRIVATE(/c/)
      external aa
      double precision pi
      parameter(pi=3.14159265358979d0)
      complex *16 ii
      parameter(ii=(0d0,1d0))
      zz=z
      acc=1d-3
c     Real part
      h1=(0d0,0d0)
      j=1
      do i=1,3     
      call dadapt(aa,0d0,1d0,5,1d-15,1d-8,tmp1,err)
      h1=h1+tmp1
      enddo
c     Imaginary part
      j=2
      do i=1,3
      call dadapt(aa,0d0,1d0,5,1d-15,1d-8,tmp2,err)
      h1=h1+II*tmp2
      enddo
      h1=-h1/pi
      return
      end

      function h2(z)
      implicit none
      complex *16 h2,z,zz
      common/c/zz
!$OMP THREADPRIVATE(/c/)
      real *8 t,v,acc,tmp1,tmp2
      real *8 err,aa
      integer i,j
      common/ij/i,j      
!$OMP THREADPRIVATE(/ij/)
      common/v/v
c !$OMP THREADPRIVATE(/v/)
      double precision pi
      parameter(pi=3.14159265358979d0)
      complex *16 ii
      parameter(ii=(0d0,1d0))
      external aa
      zz=z
      acc=1d-3
c     Real part
      h2=(0d0,0d0)
      j=1
      do i=4,6     
      call dadapt(aa,0d0,1d0,5,1d-15,1d-8,tmp1,err)
      h2=h2+tmp1
      enddo
c     Imaginary part
      j=2
      do i=4,6
      call dadapt(aa,0d0,1d0,5,1d-15,1d-8,tmp2,err)
      h2=h2+II*tmp2
      enddo
      h2=-h2/pi
      return
      end
 

c     The function aa(t) give the parametrization of the real (j=1) 
c     or imaginary (j=2) part of the function Exp[-i z sin(teta)]
c     and includes the jacobian
c     i=1,2,3 needed to compute h1
c     i=4,5,6 needed to compute h2

      function aa(t)
      implicit none
      complex *16 zz,a,teta,jac
      real *8 v,aa,t
      integer i,j
      common/c/zz
!$OMP THREADPRIVATE(/c/)
      common/ij/i,j
!$OMP THREADPRIVATE(/ij/)
      common/v/v      
c !$OMP THREADPRIVATE(/v/)
      double precision pi
      parameter(pi=3.14159265358979d0)
      complex *16 ii
      parameter(ii=(0d0,1d0))
      if(i.eq.1) then
      teta=-II*pi*v*t
      jac=II*pi*v
      elseif(i.eq.2) then
      teta=-pi*t
      jac=-pi
      elseif(i.eq.3) then
      teta=-Pi+II*v*pi*t
      jac=II*v*pi
      elseif(i.eq.4) then
      teta=pi+II*pi*v*t
      jac=-II*pi*v
      elseif(i.eq.5) then  
      teta=pi*t
      jac=-pi
      elseif(i.eq.6) then
      teta=-II*v*pi*t
      jac=-II*v*pi
      endif   
      a=jac*exp(-II*zz*sin(teta))
      if(j.eq.1) then
      aa=dreal(a)
      elseif(j.eq.2) then
      aa=dimag(a)
      endif
      return
      end
