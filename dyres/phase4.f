      subroutine phase4(r,p1,p2,p3,p4,p5,p6,wt,*)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'process.f'
c---- generate phase space for 2-->4 process
c---- r(mxdim),p1(4),p2(4) are inputs reversed in sign 
c---- from physical values 
c---- phase space for -p1-p2 --> p3+p4+p5+p6
c---- with all 2 pi's (ie 1/(2*pi)^8)
      double precision r(mxdim)
      double precision p1(4),p2(4),p3(4),p4(4),p5(4),p6(4)
      double precision p12(4),p34(4),p56(4)
      double precision wt,wt3456,wt34,wt56,wt0
      integer j
      parameter(wt0=1d0/twopi**2)

      do j=1,4
      p12(j)=-p1(j)-p2(j)
      enddo
c p56 is the b-bbar system
      call phi1_2(r(1),r(2),r(3),r(4),p12,p56,p34,wt3456,*99)
      call phi3m0(r(7),r(8),p34,p3,p4,wt34,*99)

   
      call phi3m0(r(5),r(6),p56,p5,p6,wt56,*99)
      
      wt=wt0*wt3456*wt34*wt56
      return
 
 99   wt=0d0
      return 1
      end
