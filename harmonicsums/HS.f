
      function HS(k1,k2,k3,k4,k5,n,eta)

      implicit none
      integer k1,k2,k3,k4,k5,r,n1,n2,n3,n4,n5,nsum
      complex*16 HS,n,sk1k2k3k4k5npr,sk2k3k4k5npr,sk3k4k5npr,sk4k5npr,
     .     sk5npr,HSexpand,
     .     npr,sum1,sum2,sum3,sum4,sum5,factork1(1000),factork2(2:1000),
     .     factork3(3:1000),factork4(4:1000),factork5(5:1000)
      real*8 eta,etak,signkn,etanpr

      if(abs(k1)+abs(k2)+abs(k3)+abs(k4)+abs(k5).gt.5
     .     .or.(k1.eq.0).or.(k2.eq.0.and.k3.ne.0)
     .     .or.(k3.eq.0.and.k4.ne.0).or.(k4.eq.0.and.k5.ne.0))then
         write(6,*)'HS not defined for indices',k1,k2,k3,k4,k5
         stop
      endif

      r=0
      npr=n
c     Increase this number (32.d0) for more accuracy
!      do while(abs(npr).lt.32.d0)
!         r=r+1
!         npr=n+dble(r)
!      enddo

      etanpr=eta
      if(mod(r,2).ne.0)etanpr=-eta

      sk1k2k3k4k5npr=HSexpand(k1,k2,k3,k4,k5,npr,etanpr)
      if(r.ge.1)then
         sk2k3k4k5npr=HSexpand(k2,k3,k4,k5,0,npr,etanpr)
         do nsum=1,r
            factork1(nsum)=etak(eta,k1)*signkn(nsum,k1)
     .           *cdexp(-dabs(dble(k1))*cdlog(dble(nsum)+n))
         enddo
      endif
      if(r.ge.2.and.k2.ne.0)then
         sk3k4k5npr=HSexpand(k3,k4,k5,0,0,npr,etanpr)
         do nsum=2,r
            factork2(nsum)=etak(eta,k2)*signkn(nsum,k2)
     .           *cdexp(-dabs(dble(k2))*cdlog(dble(nsum)+n))
         enddo
      endif
      if(r.ge.3.and.k3.ne.0)then
         sk4k5npr=HSexpand(k4,k5,0,0,0,npr,etanpr)
         do nsum=3,r
            factork3(nsum)=etak(eta,k3)*signkn(nsum,k3)
     .           *cdexp(-dabs(dble(k3))*cdlog(dble(nsum)+n))
         enddo
      endif
      if(r.ge.4.and.k4.ne.0)then
         sk5npr=HSexpand(k5,0,0,0,0,npr,etanpr)
         do nsum=4,r
            factork4(nsum)=etak(eta,k4)*signkn(nsum,k4)
     .           *cdexp(-dabs(dble(k4))*cdlog(dble(nsum)+n))
         enddo
      endif
      if(r.ge.5.and.k5.ne.0)then
         do nsum=5,r
            factork5(nsum)=etak(eta,k5)*signkn(nsum,k5)
     .           *cdexp(-dabs(dble(k5))*cdlog(dble(nsum)+n))
         enddo
      endif

      sum1=sk1k2k3k4k5npr
      do n1=1,r
         sum2=sk2k3k4k5npr
         if(k2.ne.0)then
            do n2=1,r-n1
               sum3=sk3k4k5npr
               if(k3.ne.0)then
                  do n3=1,r-n1-n2
                     sum4=sk4k5npr
                     if(k4.ne.0)then
                        do n4=1,r-n1-n2-n3
                           sum5=sk5npr
                           if(k5.ne.0)then
                              do n5=1,r-n1-n2-n3-n4
                                 sum5=sum5-factork5(n1+n2+n3+n4+n5)
                              enddo
                           endif
                           sum4=sum4-factork4(n1+n2+n3+n4)*sum5
                        enddo
                     endif
                     sum3=sum3-factork3(n1+n2+n3)*sum4
                  enddo
               endif
               sum2=sum2-factork2(n1+n2)*sum3
            enddo
         endif
         sum1=sum1-factork1(n1)*sum2
      enddo

      HS=sum1

      end

c-------------------------------------

      function HSexpand(k1,k2,k3,k4,k5,n,eta)

      implicit complex*16 (I,t)
      
      integer k1,k2,k3,k4,k5
      real*8 eta
      complex*16 HSexpand,n,gammaEplnn,invnp1,invnp2

      invnp1=1.d0/n
      invnp2=invnp1*invnp1

      gammaEplnn=cdlog(n)+0.577215664901532860606512090082d0

      HSexpand=(1.d0,0.d0)

      include 'HSpart.f'

      end

c------------------------------------------

      function signkn(n,k)

      implicit none

      real*8 signkn
      integer n,k

      signkn=1.d0
      if(k.lt.0.d0.and.mod(n,2).ne.0)signkn=-1.d0

      end

c------------------------------------------

      function etak(eta,k)

      implicit none

      real*8 etak,eta
      integer k

      etak=1.d0
      if(k.lt.0.d0.and.eta.lt.0.d0)etak=-1.d0

      end

c------------------------------------------
