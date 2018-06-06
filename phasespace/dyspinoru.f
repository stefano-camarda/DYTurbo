      subroutine spinoru(N,p,za,zb)
c---Calculate spinor products
c---extended to deal with negative energies ie with all momenta outgoing
c---Arbitrary conventions of Bern, Dixon, Kosower, Weinzierl, 
c---za(i,j)*zb(j,i)=s(i,j)
c---  26/08/2016 Fixed numerical instabilities for rt(j) close to 0 or negative
      implicit none
      double complex im,czip,cone
      parameter(im=(0d0,1d0),czip=(0d0,0d0),cone=(1d0,0d0))
      double precision two
      parameter(two=2d0)
      integer mxpart
      parameter(mxpart=6)
      double complex za(mxpart,mxpart),zb(mxpart,mxpart)
      double precision s(mxpart,mxpart)
      common/sprods/s
!$OMP THREADPRIVATE(/sprods/)
      double precision p(mxpart,4),rt(mxpart)
      double complex c23(mxpart),f(mxpart)
      double precision k1(mxpart),k2(mxpart),kp(mxpart)
      double precision cphi,sphi
      integer i,j,N
      integer sig
      integer x,y,z
      integer sigx,sigy,sigz
      
c     See arXiv:hep-ph/9601359 Sec 2.2 for definitions
c     za(i,j) = <i,j>, zb(i,j) = [i,j]
c     rt(j) = k0 + k3 = k+
c     c23(j) = k1 + i*k2
c     axis permutations lead to the same amplitudes (but not to the same <i,j>, [i,j])
c     If a 4-vector is aligned to the k3 axis, the calculation fails because k+ (rt) is zero (division by zero)
c     An axis permutation (or a rotation) should be implemented for this case to prevent failure
c     The current choice is k1=ky; k2=-kz; k3=kx; (since the incoming particles are aligned along z, k3 should never be kz)
c     The calculation also fails if one of the vectors happens to be zero

cc     ======================================================
cc     Calculation based on Eq.(15) and (16) of arXiv:hep-ph/9601359
c      do j=1,N
c         za(j,j)=czip
c         zb(j,j)=czip
c
cC-----positive energy case
c         if (p(j,4) .gt. 0d0) then
c            rt(j)=sqrt(max(0d0,p(j,4)+p(j,1)))
c            c23(j)=dcmplx(p(j,3),-p(j,2))
c            f(j)=cone
c         else
cC-----negative energy case
c            rt(j)=sqrt(max(0d0,-p(j,4)-p(j,1)))
c            c23(j)=dcmplx(-p(j,3),p(j,2))
c            f(j)=im
c         endif
c
c      enddo
c      
c      do i=2,N
c         do j=1,i-1
c            s(i,j)=two*(p(i,4)*p(j,4)-p(i,1)*p(j,1)
c     &           -p(i,2)*p(j,2)-p(i,3)*p(j,3))
cc     In case one rt is small this empirical solution works
c            if (rt(j).lt.1d-4.or.rt(i).lt.1d-4) then
c               za(i,j)=
c     &              conjg(sqrt(dcmplx(s(i,j))))
c               zb(i,j)=-s(i,j)/za(i,j)
c            else
c               za(i,j)=f(i)*f(j)
c     &         *(c23(i)*rt(j)/rt(i)-c23(j)*rt(i)/rt(j))
c
c               if (abs(s(i,j)).lt.1d-9) then
c                  zb(i,j)=-(f(i)*f(j))**2*conjg(za(i,j))
c               else
c                  zb(i,j)=-s(i,j)/za(i,j)
c               endif
c            endif
c            za(j,i)=-za(i,j)
c            zb(j,i)=-zb(i,j)
c            s(j,i)=s(i,j)
c         enddo
c      enddo

cc     ======================================================
cc     Calculation based on Eq.(15) and (16) of arXiv:hep-ph/9601359
cc     do not calculate sij, avoid if conditions
c      do j=1,N
c         za(j,j)=czip
c         zb(j,j)=czip
c         
cc     Compact notation for positive and negative energies
c         sig=int(sign(1d0,p(j,4)))
c         f(j)=sqrt(dcmplx(sig))
c         rt(j)=sqrt(max(0d0,sig*p(j,4)+sig*p(j,1)))
c         c23(j)=dcmplx(sig*p(j,3),-sig*p(j,2))
c      enddo
c
c      do i=2,N
c         do j=1,i-1
c            za(i,j)=f(i)*f(j)*(c23(i)*rt(j)/rt(i)-c23(j)*rt(i)/rt(j))
c            zb(i,j)=-(f(i)*f(j))**2*conjg(za(i,j))
c            za(j,i)=-za(i,j)
c            zb(j,i)=-zb(i,j)
c         enddo
c      enddo
      
c     ======================================================
c     Calculation based on Eq.(16) and (17) of arXiv:hep-ph/9601359
c     do not calculate sij, avoid if conditions, do not set diagonal elements to zero
c     This is the fastest method
      
      do j=1,N
         sig=int(sign(1d0,p(j,4)))
         f(j)=sqrt(dcmplx(sig))
         
         kp(j)=sig*p(j,4)+sig*p(j,1)
         rt(j)=sqrt(kp(j))
         k1(j)=sig*p(j,3)
         k2(j)=-sig*p(j,2)
      enddo
      
      do i=2,N
         do j=1,i-1

            cphi=(k1(i)*kp(j)-k1(j)*kp(i))/(rt(i)*rt(j))
            sphi=(k2(i)*kp(j)-k2(j)*kp(i))/(rt(i)*rt(j))
            za(i,j)=f(i)*f(j)*dcmplx(cphi,sphi)
            zb(i,j)=f(i)*f(j)*dcmplx(-cphi,sphi)
            
            za(j,i)=-za(i,j)
            zb(j,i)=-zb(i,j)
         enddo
      enddo


cc     ======================================================
cc     Compare the two methods
cC     choice of axis
c!     y -> x  
c!     z -> -y 
c!     x -> z  
c      x = 3 
c      y = 2 
c      z = 1 
c      sigx = 1
c      sigy = -1
c      sigz = 1
c
c      do j=1,N
c         sig=int(sign(1d0,p(j,4)))
c         f(j)=sqrt(dcmplx(sig))
c         rt(j)=sqrt(max(0d0,sig*p(j,4)+sig*sigz*p(j,z)))
c         c23(j)=dcmplx(sig*sigx*p(j,x),sig*sigy*p(j,y))
c
c         kp(j)=sig*p(j,4)+sig*sigz*p(j,z)
c         k1(j)=sig*sigx*p(j,x)
c         k2(j)=sig*sigy*p(j,y)
c      enddo
c
c      print *,
c      do i=2,N
c         do j=1,i-1
c            za(i,j)=f(i)*f(j)*(c23(i)*rt(j)/rt(i)-c23(j)*rt(i)/rt(j))
c            zb(i,j)=-(f(i)*f(j))**2*conjg(za(i,j))
c            print *,'Eq.(15-16)',i,j,za(i,j),zb(i,j)
c
c            cphi=(k1(i)*kp(j)-k1(j)*kp(i))/(rt(i)*rt(j))
c            sphi=(k2(i)*kp(j)-k2(j)*kp(i))/(rt(i)*rt(j))
c            za(i,j)=f(i)*f(j)*dcmplx(cphi,sphi)
c            zb(i,j)=f(i)*f(j)*dcmplx(-cphi,sphi)
c            print *,'Eq.(16-17)',i,j,za(i,j),zb(i,j)
c
c            za(j,i)=-za(i,j)
c            zb(j,i)=-zb(i,j)
c         enddo
c      enddo
      
      
      return
      end
