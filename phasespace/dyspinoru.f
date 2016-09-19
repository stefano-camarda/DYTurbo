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
      double complex za(12,12),zb(12,12)
      double precision s(12,12)
      common/sprods/s
      double precision p(12,4),rt(12)
      double complex c23(12),f(12)
      integer i,j,N

c     See arXiv:hep-ph/9601359 Sec 2.2 for definitions
c     za(i,j) = <i,j>, zb(i,j) = [i,j]
c     rt(j) = k0 + k3 = k+
c     c23(j) = k1 + i*k2
      
c---if one of the vectors happens to be zero this routine fails.
      do j=1,N
         za(j,j)=czip
         zb(j,j)=za(j,j)

C-----positive energy case
         if (p(j,4) .gt. 0d0) then
            rt(j)=dsqrt(max(0d0,p(j,4)+p(j,1)))
            c23(j)=dcmplx(p(j,3),-p(j,2))
            f(j)=cone
         else
C-----negative energy case
            rt(j)=dsqrt(max(0d0,-p(j,4)-p(j,1)))
            c23(j)=dcmplx(-p(j,3),p(j,2))
            f(j)=im
         endif
      enddo
      do i=2,N
         do j=1,i-1
            s(i,j)=two*(p(i,4)*p(j,4)-p(i,1)*p(j,1)
     &           -p(i,2)*p(j,2)-p(i,3)*p(j,3))
            if (rt(j).lt.1d-4.or.rt(i).lt.1d-4) then
               za(i,j)=
     &              dconjg(sqrt(dcmplx(s(i,j))))
               zb(i,j)=-dcmplx(s(i,j))/za(i,j)
            else
               za(i,j)=f(i)*f(j)
     &         *(c23(i)*dcmplx(rt(j)/rt(i))-c23(j)*dcmplx(rt(i)/rt(j)))

               if (abs(s(i,j)).lt.1d-9) then
                  zb(i,j)=-(f(i)*f(j))**2*dconjg(za(i,j))
               else
                  zb(i,j)=-dcmplx(s(i,j))/za(i,j)
               endif
            endif
            za(j,i)=-za(i,j)
            zb(j,i)=-zb(i,j)
            s(j,i)=s(i,j)
         enddo
      enddo

      return
      end
