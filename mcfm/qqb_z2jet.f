      subroutine qqb_z2jet(p,msq)
      implicit none
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) --> Z +g(p5) +g(p6)
c                          |
c                          --> l(p3)+a(p4)
c                            
c--all momenta incoming
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'msq_cs.f'
      include 'flags.f'
      integer i,j,k,pq,pl,nquark,swap(2),swap1(0:2),nup,ndo,
     .     j1,j2,j3,icol
      integer a,b,c
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac,faclo,
     .   qqbZgg2(2,2),qgZqg2(2,2),
c    .   qbqZgg2(2,2),qbgZqbg2(2,2),gqbZqbg2(2,2),
     .   gqZqg2(2,2),ggZqbq2(2,2),
     .   qqbZgg2_cs(0:2,2,2),qbqZgg2_cs(0:2,2,2),
     .   qgZqg2_cs(0:2,2,2),gqZqg2_cs(0:2,2,2),
     .   qbgZqbg2_cs(0:2,2,2),gqbZqbg2_cs(0:2,2,2),
     .   ggZqbq2_cs(0:2,2,2),ggtemp(0:2)
      double precision tup,tdo
      double complex a111,a112,a121,a211,a122,a212,a221,a222
      double complex b111,b112,b121,b211,b122,b212,b221,b222

      double complex a111v(nf),a112v(nf),a121v(nf),a211v(nf),
     .     a122v(nf),a212v(nf),a221v(nf),a222v(nf)

      double complex b111v(nf),b112v(nf),b121v(nf),b211v(nf),
     .     b122v(nf),b212v(nf),b221v(nf),b222v(nf)
      
      double complex qRb_a(2,2,2),qRb_b(2,2,2)
      double complex qqb_a(2,2,2),qqb_b(2,2,2),prop

      double complex qbq_a(2,2,2),qbq_b(2,2,2)
      double complex qbR_a(2,2,2),qbR_b(2,2,2)

      double complex qq_a(2,2,2),qq_b(2,2,2)
      double complex qR_a(2,2,2),qR_b(2,2,2)

      double complex qbRb_a(2,2,2),qbRb_b(2,2,2)
      double complex qbqb_a(2,2,2),qbqb_b(2,2,2)

      double precision mqq(0:2,fn:nf,fn:nf)
      common/mqq/mqq
      double precision mqq0,mqq1,mqq2

      double precision mqq0v(nf),mqq1v(nf),mqq2v(nf),mqqtotv(nf)

      double precision msq_1(nf,nf),msq_2(nf,nf),
     .     msq_3(nf,nf),msq_4(nf,nf)

      double precision msqd_1(nf),msqd_2(nf),msqd_3(nf),msqd_4(nf)
      
      double precision ll(nf),rr(nf),qq(nf)
      double precision msq_pos(nf,nf)
      
      data swap/2,1/
      save swap
      data swap1/0,2,1/
      save swap1

      double complex bosprop
      external bosprop

      msq=0d0
c      do j=-nf,nf
c      do k=-nf,nf
c      msq(j,k)=0d0
c      enddo
c      enddo


      call spinoru(6,p,za,zb)
c      prop=s(3,4)/dcmplx((s(3,4)-zmass**2),zmass*zwidth)
      prop=bosprop(s(3,4))

c--- calculate 2-quark, 2-gluon amplitudes
      if (Gflag) then
        call z2jetsq(1,2,3,4,5,6,za,zb,qqbZgg2)
        call storecsz(qqbZgg2_cs)
        call z2jetsq(1,5,3,4,2,6,za,zb,qgZqg2)
        call storecsz(qgZqg2_cs)
        call z2jetsq(2,5,3,4,1,6,za,zb,gqZqg2)
        call storecsz(gqZqg2_cs)

        do j=1,2
        do k=1,2
c        qbqZgg2(j,k)=qqbZgg2(swap(j),k)
c        qbgZqbg2(j,k)=qgZqg2(swap(j),k)
c        gqbZqbg2(j,k)=gqZqg2(swap(j),k)
        do i=0,2
        qbqZgg2_cs(i,j,k)=qqbZgg2_cs(swap1(i),swap(j),k)
        qbgZqbg2_cs(i,j,k)=qgZqg2_cs(swap1(i),swap(j),k)
        gqbZqbg2_cs(i,j,k)=gqZqg2_cs(swap1(i),swap(j),k)
        enddo
        enddo
        enddo

c        call z2jetsq(2,1,3,4,5,6,za,zb,qbqZgg2)
c        call storecsz(qbqZgg2_cs)
c        call z2jetsq(5,1,3,4,2,6,za,zb,qbgZqbg2)
c        call storecsz(qbgZqbg2_cs)
c        call z2jetsq(5,2,3,4,1,6,za,zb,gqbZqbg2)
c        call storecsz(gqbZqbg2_cs)
        
C --NB this is the matrix element for gg->Z qb(5) q(6)
        call z2jetsq(5,6,3,4,1,2,za,zb,ggZqbq2)
        call storecsz(ggZqbq2_cs)        

        fac=v*xn/four*(esq*gsq)**2
        do pq=1,2
        do pl=1,2
        do i=0,2        
          qqbZgg2_cs(i,pq,pl) = half*aveqq*fac*qqbZgg2_cs(i,pq,pl)
          qbqZgg2_cs(i,pq,pl) = half*aveqq*fac*qbqZgg2_cs(i,pq,pl)
          gqZqg2_cs(i,pq,pl)  = aveqg*fac*gqZqg2_cs(i,pq,pl)
          qgZqg2_cs(i,pq,pl)  = aveqg*fac*qgZqg2_cs(i,pq,pl)
          gqbZqbg2_cs(i,pq,pl)= aveqg*fac*gqbZqbg2_cs(i,pq,pl)
          qbgZqbg2_cs(i,pq,pl)= aveqg*fac*qbgZqbg2_cs(i,pq,pl)
          ggZqbq2_cs(i,pq,pl) = avegg*fac*ggZqbq2_cs(i,pq,pl) 
       enddo

        qqbZgg2(pq,pl) = qqbZgg2_cs(1,pq,pl)+qqbZgg2_cs(2,pq,pl)
     .                  +qqbZgg2_cs(0,pq,pl) 
        gqZqg2(pq,pl)  = gqZqg2_cs(1,pq,pl) +gqZqg2_cs(2,pq,pl)
     .                  +gqZqg2_cs(0,pq,pl)  
        qgZqg2(pq,pl)  = qgZqg2_cs(1,pq,pl)  +qgZqg2_cs(2,pq,pl)
     .                  +qgZqg2_cs(0,pq,pl)  
c        qbqZgg2(pq,pl) = qbqZgg2_cs(1,pq,pl)+qbqZgg2_cs(2,pq,pl)
c     .                  +qbqZgg2_cs(0,pq,pl) 
c        gqbZqbg2(pq,pl)= gqbZqbg2_cs(1,pq,pl)+gqbZqbg2_cs(2,pq,pl)
c     .                  +gqbZqbg2_cs(0,pq,pl)
c        qbgZqbg2(pq,pl)= qbgZqbg2_cs(1,pq,pl)+qbgZqbg2_cs(2,pq,pl)
c     .                  +qbgZqbg2_cs(0,pq,pl)
        ggZqbq2(pq,pl) = ggZqbq2_cs(1,pq,pl) +ggZqbq2_cs(2,pq,pl)
     .                  +ggZqbq2_cs(0,pq,pl) 
        enddo
        enddo
      endif


      if (Qflag) then
c      call spinoru(6,p,za,zb)

c--- qRb->qRb
      call ampqqb_qqb(1,5,2,6,qRb_a,qRb_b)
c--- qR->qR
c instead of calling ampqqb_qqb(1,5,6,2,qR_a,qR_b)
      do j1=1,2
      do j2=1,2
      do j3=1,2
      qR_a(j1,j2,j3)=+qRb_a(j1,swap(j2),j3)
      qR_b(j1,j2,j3)=-qRb_b(j1,swap(j2),j3)
      enddo
      enddo
      enddo
c--- qbR->qbR
      call ampqqb_qqb(6,1,5,2,qbR_a,qbR_b)

c--- qbRb->qbRb
c instead of calling ampqqb_qqb(5,1,2,6,qbRb_a,qbRb_b)
      do j1=1,2
      do j2=1,2
      do j3=1,2
      qbRb_a(j1,j2,j3)=-qRb_a(swap(j1),j2,j3)
      qbRb_b(j1,j2,j3)=+qRb_b(swap(j1),j2,j3)
      enddo
      enddo
      enddo

c--- qqb->qqb
      call ampqqb_qqb(1,2,5,6,qqb_a,qqb_b)
c--- qbq->qqb
c instead of calling ampqqb_qqb(2,1,5,6,qbq_a,qbq_b)
      do j1=1,2
      do j2=1,2
      do j3=1,2
      qbq_a(j1,j2,j3)=-qqb_a(swap(j1),j2,j3)
      qbq_b(j1,j2,j3)=+qqb_b(swap(j1),j2,j3)
      enddo
      enddo
      enddo

c--- qq->qq
      call ampqqb_qqb(1,6,5,2,qq_a,qq_b)
c--- qbqb->qbqb
c instead of calling ampqqb_qqb(6,1,2,5,qbqb_a,qbqb_b)
      do j1=1,2
      do j2=1,2
      do j3=1,2
      qbqb_a(j1,j2,j3)=-qq_a(swap(j1),swap(j2),j3)
      qbqb_b(j1,j2,j3)=-qq_b(swap(j1),swap(j2),j3)
      enddo
      enddo
      enddo

      faclo=4d0*V*gsq**2*esq**2*aveqq 
      endif



c      if (Gflag) then 
c         
c      do j=-nf,nf
c      do k=-nf,nf
c      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 19
c
c      if     ((j .eq. 0) .and. (k .eq. 0)) then
c
c          do icol=0,2
c          ggtemp(icol)=0d0
c          do nquark=1,nf
c           ggtemp(icol)=ggtemp(icol)
c     .      +abs(Q(nquark)*q1+L(nquark)*l1*prop)**2*ggZqbq2_cs(icol,1,1)
c     .      +abs(Q(nquark)*q1+R(nquark)*r1*prop)**2*ggZqbq2_cs(icol,2,2)
c     .      +abs(Q(nquark)*q1+L(nquark)*r1*prop)**2*ggZqbq2_cs(icol,1,2)
c     .      +abs(Q(nquark)*q1+R(nquark)*l1*prop)**2*ggZqbq2_cs(icol,2,1)
c          enddo
c          msq_cs(icol,j,k)=ggtemp(icol)
c          enddo
c      elseif ((j .gt. 0) .and. (k .lt. 0)) then
c          do icol=0,2
c             msq_cs(icol,j,k)=
c     .       +abs(Q(j)*q1+L(j)*l1*prop)**2*qqbZgg2_cs(icol,1,1)
c     .       +abs(Q(j)*q1+R(j)*r1*prop)**2*qqbZgg2_cs(icol,2,2)
c     .       +abs(Q(j)*q1+L(j)*r1*prop)**2*qqbZgg2_cs(icol,1,2)
c     .       +abs(Q(j)*q1+R(j)*l1*prop)**2*qqbZgg2_cs(icol,2,1)
c          enddo
cc---Statistical factor already included above
c      elseif ((j .lt. 0) .and. (k .gt. 0)) then
c          do icol=0,2
c             msq_cs(icol,j,k)=
c     .       +abs(Q(k)*q1+L(k)*l1*prop)**2*qbqZgg2_cs(icol,1,1)
c     .       +abs(Q(k)*q1+R(k)*r1*prop)**2*qbqZgg2_cs(icol,2,2)
c     .       +abs(Q(k)*q1+L(k)*r1*prop)**2*qbqZgg2_cs(icol,1,2)
c     .       +abs(Q(k)*q1+R(k)*l1*prop)**2*qbqZgg2_cs(icol,2,1)
c          enddo
c      elseif ((j .gt. 0) .and. (k .eq. 0)) then
c          do icol=0,2
c             msq_cs(icol,j,k)=
c     .       +abs(Q(j)*q1+L(j)*l1*prop)**2*qgZqg2_cs(icol,1,1)
c     .       +abs(Q(j)*q1+R(j)*r1*prop)**2*qgZqg2_cs(icol,2,2)
c     .       +abs(Q(j)*q1+L(j)*r1*prop)**2*qgZqg2_cs(icol,1,2)
c     .       +abs(Q(j)*q1+R(j)*l1*prop)**2*qgZqg2_cs(icol,2,1)
c          enddo
c      elseif ((j .lt. 0) .and. (k .eq. 0)) then
c          do icol=0,2
c             msq_cs(icol,j,k)=
c     .       +abs(Q(-j)*q1+L(-j)*l1*prop)**2*qbgZqbg2_cs(icol,1,1)
c     .       +abs(Q(-j)*q1+R(-j)*r1*prop)**2*qbgZqbg2_cs(icol,2,2)
c     .       +abs(Q(-j)*q1+L(-j)*r1*prop)**2*qbgZqbg2_cs(icol,1,2)
c     .       +abs(Q(-j)*q1+R(-j)*l1*prop)**2*qbgZqbg2_cs(icol,2,1)
c          enddo
c      elseif ((j .eq. 0) .and. (k .gt. 0)) then
c          do icol=0,2
c             msq_cs(icol,j,k)=
c     .       +abs(Q(k)*q1+L(k)*l1*prop)**2*gqZqg2_cs(icol,1,1)
c     .       +abs(Q(k)*q1+R(k)*r1*prop)**2*gqZqg2_cs(icol,2,2)
c     .       +abs(Q(k)*q1+L(k)*r1*prop)**2*gqZqg2_cs(icol,1,2)
c     .       +abs(Q(k)*q1+R(k)*l1*prop)**2*gqZqg2_cs(icol,2,1)
c          enddo
c      elseif ((j .eq. 0) .and. (k .lt. 0)) then
c          do icol=0,2
c             msq_cs(icol,j,k)=
c     .       +abs(Q(-k)*q1+L(-k)*l1*prop)**2*gqbZqbg2_cs(icol,1,1)
c     .       +abs(Q(-k)*q1+R(-k)*r1*prop)**2*gqbZqbg2_cs(icol,2,2)
c     .       +abs(Q(-k)*q1+L(-k)*r1*prop)**2*gqbZqbg2_cs(icol,1,2)
c     .       +abs(Q(-k)*q1+R(-k)*l1*prop)**2*gqbZqbg2_cs(icol,2,1)
c          enddo
c      endif
c      msq(j,k)=msq_cs(0,j,k)+msq_cs(1,j,k)+msq_cs(2,j,k)
c
c   19 continue
c      enddo
c      enddo
c      endif

c      if (Gflag) then !(two external quark lines and two external gluons)
      j=0
      k=0
      do icol=0,2
         ggtemp(icol)=0d0
         do nquark=1,nf
            ggtemp(icol)=ggtemp(icol)
     .      +abs(Q(nquark)*q1+L(nquark)*l1*prop)**2*ggZqbq2_cs(icol,1,1)
     .      +abs(Q(nquark)*q1+R(nquark)*r1*prop)**2*ggZqbq2_cs(icol,2,2)
     .      +abs(Q(nquark)*q1+L(nquark)*r1*prop)**2*ggZqbq2_cs(icol,1,2)
     .      +abs(Q(nquark)*q1+R(nquark)*l1*prop)**2*ggZqbq2_cs(icol,2,1)
         enddo
         msq_cs(icol,j,k)=ggtemp(icol)
      enddo
      msq(j,k)=msq_cs(0,j,k)+msq_cs(1,j,k)+msq_cs(2,j,k)
      
      do j=1,nf
         k=-j
         do icol=0,2
            msq_cs(icol,j,k)=
     .       +abs(Q(j)*q1+L(j)*l1*prop)**2*qqbZgg2_cs(icol,1,1)
     .       +abs(Q(j)*q1+R(j)*r1*prop)**2*qqbZgg2_cs(icol,2,2)
     .       +abs(Q(j)*q1+L(j)*r1*prop)**2*qqbZgg2_cs(icol,1,2)
     .       +abs(Q(j)*q1+R(j)*l1*prop)**2*qqbZgg2_cs(icol,2,1)
         enddo
         msq(j,k)=msq_cs(0,j,k)+msq_cs(1,j,k)+msq_cs(2,j,k)
      enddo
c---Statistical factor already included above
      do j=-nf,-1
         k=-j
         do icol=0,2
            msq_cs(icol,j,k)=
     .       +abs(Q(k)*q1+L(k)*l1*prop)**2*qbqZgg2_cs(icol,1,1)
     .       +abs(Q(k)*q1+R(k)*r1*prop)**2*qbqZgg2_cs(icol,2,2)
     .       +abs(Q(k)*q1+L(k)*r1*prop)**2*qbqZgg2_cs(icol,1,2)
     .       +abs(Q(k)*q1+R(k)*l1*prop)**2*qbqZgg2_cs(icol,2,1)
         enddo
         msq(j,k)=msq_cs(0,j,k)+msq_cs(1,j,k)+msq_cs(2,j,k)
      enddo
      
      k=0
      do j=1,nf
         do icol=0,2
            msq_cs(icol,j,k)=
     .       +abs(Q(j)*q1+L(j)*l1*prop)**2*qgZqg2_cs(icol,1,1)
     .       +abs(Q(j)*q1+R(j)*r1*prop)**2*qgZqg2_cs(icol,2,2)
     .       +abs(Q(j)*q1+L(j)*r1*prop)**2*qgZqg2_cs(icol,1,2)
     .       +abs(Q(j)*q1+R(j)*l1*prop)**2*qgZqg2_cs(icol,2,1)
         enddo
         msq(j,k)=msq_cs(0,j,k)+msq_cs(1,j,k)+msq_cs(2,j,k)
      enddo
      do j=-nf,-1
         do icol=0,2
            msq_cs(icol,j,k)=
     .       +abs(Q(-j)*q1+L(-j)*l1*prop)**2*qbgZqbg2_cs(icol,1,1)
     .       +abs(Q(-j)*q1+R(-j)*r1*prop)**2*qbgZqbg2_cs(icol,2,2)
     .       +abs(Q(-j)*q1+L(-j)*r1*prop)**2*qbgZqbg2_cs(icol,1,2)
     .       +abs(Q(-j)*q1+R(-j)*l1*prop)**2*qbgZqbg2_cs(icol,2,1)
         enddo
         msq(j,k)=msq_cs(0,j,k)+msq_cs(1,j,k)+msq_cs(2,j,k)
      enddo
      
      j=0
      do k=1,nf
         do icol=0,2
            msq_cs(icol,j,k)=
     .       +abs(Q(k)*q1+L(k)*l1*prop)**2*gqZqg2_cs(icol,1,1)
     .       +abs(Q(k)*q1+R(k)*r1*prop)**2*gqZqg2_cs(icol,2,2)
     .       +abs(Q(k)*q1+L(k)*r1*prop)**2*gqZqg2_cs(icol,1,2)
     .       +abs(Q(k)*q1+R(k)*l1*prop)**2*gqZqg2_cs(icol,2,1)
         enddo
         msq(j,k)=msq_cs(0,j,k)+msq_cs(1,j,k)+msq_cs(2,j,k)
      enddo
      do k=-nf,-1
         do icol=0,2
            msq_cs(icol,j,k)=
     .       +abs(Q(-k)*q1+L(-k)*l1*prop)**2*gqbZqbg2_cs(icol,1,1)
     .       +abs(Q(-k)*q1+R(-k)*r1*prop)**2*gqbZqbg2_cs(icol,2,2)
     .       +abs(Q(-k)*q1+L(-k)*r1*prop)**2*gqbZqbg2_cs(icol,1,2)
     .       +abs(Q(-k)*q1+R(-k)*l1*prop)**2*gqbZqbg2_cs(icol,2,1)
         enddo
         msq(j,k)=msq_cs(0,j,k)+msq_cs(1,j,k)+msq_cs(2,j,k)
      enddo


c      endif


!      if (Qflag) then !(all four external lines are quarks)
!         
!      do j=-nf,nf
!      do k=-nf,nf
!      
!      do icol=0,2
!      mqq(icol,j,k)=zip
!      enddo
!      
!          if ((j .gt. 0) .and. (k .gt. 0)) then
!c----QQ case
!            if (j .ne. k) then
!            a111=(Q(j)*q1+L(j)*l1*prop)*qR_a(1,1,1)
!     .          +(Q(k)*q1+L(k)*l1*prop)*qR_b(1,1,1)
!            a121=(Q(j)*q1+L(j)*l1*prop)*qR_a(1,2,1)
!     .          +(Q(k)*q1+R(k)*l1*prop)*qR_b(1,2,1)
!            a112=(Q(j)*q1+L(j)*r1*prop)*qR_a(1,1,2)
!     .          +(Q(k)*q1+L(k)*r1*prop)*qR_b(1,1,2)
!            a122=(Q(j)*q1+L(j)*r1*prop)*qR_a(1,2,2)
!     .          +(Q(k)*q1+R(k)*r1*prop)*qR_b(1,2,2)
!            a211=(Q(j)*q1+R(j)*l1*prop)*qR_a(2,1,1)
!     .          +(Q(k)*q1+L(k)*l1*prop)*qR_b(2,1,1)
!            a221=(Q(j)*q1+R(j)*l1*prop)*qR_a(2,2,1)
!     .          +(Q(k)*q1+R(k)*l1*prop)*qR_b(2,2,1)
!            a212=(Q(j)*q1+R(j)*r1*prop)*qR_a(2,1,2)
!     .          +(Q(k)*q1+L(k)*r1*prop)*qR_b(2,1,2)
!            a222=(Q(j)*q1+R(j)*r1*prop)*qR_a(2,2,2)
!     .          +(Q(k)*q1+R(k)*r1*prop)*qR_b(2,2,2)
!            mqq(0,j,k)=zip
!            mqq(1,j,k)=
!     .      +faclo*(abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
!     .             +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)
!            mqq(2,j,k)=zip
!            elseif (j .eq. k) then
!            a111=(Q(j)*q1+L(j)*l1*prop)*(qR_a(1,1,1)+qR_b(1,1,1))
!            b111=(Q(j)*q1+L(j)*l1*prop)*(qq_a(1,1,1)+qq_b(1,1,1))
!            a112=(Q(j)*q1+L(j)*r1*prop)*(qR_a(1,1,2)+qR_b(1,1,2))
!            b112=(Q(j)*q1+L(j)*r1*prop)*(qq_a(1,1,2)+qq_b(1,1,2))
!            a221=(Q(j)*q1+R(j)*l1*prop)*(qR_a(2,2,1)+qR_b(2,2,1))
!            b221=(Q(j)*q1+R(j)*l1*prop)*(qq_a(2,2,1)+qq_b(2,2,1))
!            a222=(Q(j)*q1+R(j)*r1*prop)*(qR_a(2,2,2)+qR_b(2,2,2))
!            b222=(Q(j)*q1+R(j)*r1*prop)*(qq_a(2,2,2)+qq_b(2,2,2))
!
!            a121=(Q(j)*q1+L(j)*l1*prop)*qR_a(1,2,1)
!     .          +(Q(k)*q1+R(k)*l1*prop)*qR_b(1,2,1)
!            b121=(Q(j)*q1+L(j)*l1*prop)*qq_a(1,2,1)
!     .          +(Q(k)*q1+R(k)*l1*prop)*qq_b(1,2,1)
!            a122=(Q(j)*q1+L(j)*r1*prop)*qR_a(1,2,2)
!     .          +(Q(k)*q1+R(k)*r1*prop)*qR_b(1,2,2)
!            b122=(Q(j)*q1+L(j)*r1*prop)*qq_a(1,2,2)
!     .          +(Q(k)*q1+R(k)*r1*prop)*qq_b(1,2,2)
!            a211=(Q(j)*q1+R(j)*l1*prop)*qR_a(2,1,1)
!     .          +(Q(k)*q1+L(k)*l1*prop)*qR_b(2,1,1)
!            b211=(Q(j)*q1+R(j)*l1*prop)*qq_a(2,1,1)
!     .          +(Q(k)*q1+L(k)*l1*prop)*qq_b(2,1,1)
!            a212=(Q(j)*q1+R(j)*r1*prop)*qR_a(2,1,2)
!     .          +(Q(k)*q1+L(k)*r1*prop)*qR_b(2,1,2)
!            b212=(Q(j)*q1+R(j)*r1*prop)*qq_a(2,1,2)
!     .          +(Q(k)*q1+L(k)*r1*prop)*qq_b(2,1,2)
!
!            mqq(0,j,k)=half*faclo*(
!     .      +Dble(a111*Dconjg(b111))+Dble(a112*Dconjg(b112))
!     .      +Dble(a221*Dconjg(b221))+Dble(a222*Dconjg(b222)))*two/xn
!            mqq(1,j,k)=half*faclo*
!     .      (abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
!     .      +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)
!            mqq(2,j,k)=half*faclo*(
!     .      +abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
!     .      +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)
!            endif
!          elseif ((j .lt. 0) .and. (k .lt. 0)) then
!C----QbQb case
!            if (j .ne. k) then
!            a111=(Q(-j)*q1+L(-j)*l1*prop)*qbRb_a(1,1,1)
!     .          +(Q(-k)*q1+L(-k)*l1*prop)*qbRb_b(1,1,1)
!            a121=(Q(-j)*q1+L(-j)*l1*prop)*qbRb_a(1,2,1)
!     .          +(Q(-k)*q1+R(-k)*l1*prop)*qbRb_b(1,2,1)
!
!            a112=(Q(-j)*q1+L(-j)*r1*prop)*qbRb_a(1,1,2)
!     .          +(Q(-k)*q1+L(-k)*r1*prop)*qbRb_b(1,1,2)
!            a122=(Q(-j)*q1+L(-j)*r1*prop)*qbRb_a(1,2,2)
!     .          +(Q(-k)*q1+R(-k)*r1*prop)*qbRb_b(1,2,2)
!
!            a211=(Q(-j)*q1+R(-j)*l1*prop)*qbRb_a(2,1,1)
!     .          +(Q(-k)*q1+L(-k)*l1*prop)*qbRb_b(2,1,1)
!            a221=(Q(-j)*q1+R(-j)*l1*prop)*qbRb_a(2,2,1)
!     .          +(Q(-k)*q1+R(-k)*l1*prop)*qbRb_b(2,2,1)
!
!            a212=(Q(-j)*q1+R(-j)*r1*prop)*qbRb_a(2,1,2)
!     .          +(Q(-k)*q1+L(-k)*r1*prop)*qbRb_b(2,1,2)
!            a222=(Q(-j)*q1+R(-j)*r1*prop)*qbRb_a(2,2,2)
!     .          +(Q(-k)*q1+R(-k)*r1*prop)*qbRb_b(2,2,2)
!            mqq(0,j,k)=zip
!            mqq(1,j,k)=
!     .      +faclo*(abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
!     .             +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)
!            mqq(2,j,k)=zip
!            elseif (j .eq. k) then
!
!            a111=(Q(-j)*q1+L(-j)*l1*prop)*(qbRb_a(1,1,1)+qbRb_b(1,1,1))
!            b111=(Q(-j)*q1+L(-j)*l1*prop)*(qbqb_a(1,1,1)+qbqb_b(1,1,1))
!            a112=(Q(-j)*q1+L(-j)*r1*prop)*(qbRb_a(1,1,2)+qbRb_b(1,1,2))
!            b112=(Q(-j)*q1+L(-j)*r1*prop)*(qbqb_a(1,1,2)+qbqb_b(1,1,2))
!            a221=(Q(-j)*q1+R(-j)*l1*prop)*(qbRb_a(2,2,1)+qbRb_b(2,2,1))
!            b221=(Q(-j)*q1+R(-j)*l1*prop)*(qbqb_a(2,2,1)+qbqb_b(2,2,1))
!            a222=(Q(-j)*q1+R(-j)*r1*prop)*(qbRb_a(2,2,2)+qbRb_b(2,2,2))
!            b222=(Q(-j)*q1+R(-j)*r1*prop)*(qbqb_a(2,2,2)+qbqb_b(2,2,2))
!
!
!            a121=(Q(-j)*q1+L(-j)*l1*prop)*qbRb_a(1,2,1)
!     .          +(Q(-k)*q1+R(-k)*l1*prop)*qbRb_b(1,2,1)
!            a122=(Q(-j)*q1+L(-j)*r1*prop)*qbRb_a(1,2,2)
!     .          +(Q(-k)*q1+R(-k)*r1*prop)*qbRb_b(1,2,2)
!            a211=(Q(-j)*q1+R(-j)*l1*prop)*qbRb_a(2,1,1)
!     .          +(Q(-k)*q1+L(-k)*l1*prop)*qbRb_b(2,1,1)
!            a212=(Q(-j)*q1+R(-j)*r1*prop)*qbRb_a(2,1,2)
!     .          +(Q(-k)*q1+L(-k)*r1*prop)*qbRb_b(2,1,2)
!
!            b121=(Q(-j)*q1+L(-j)*l1*prop)*qbqb_a(1,2,1)
!     .          +(Q(-k)*q1+R(-k)*l1*prop)*qbqb_b(1,2,1)
!            b122=(Q(-j)*q1+L(-j)*r1*prop)*qbqb_a(1,2,2)
!     .          +(Q(-k)*q1+R(-k)*r1*prop)*qbqb_b(1,2,2)
!            b211=(Q(-j)*q1+R(-j)*l1*prop)*qbqb_a(2,1,1)
!     .          +(Q(-k)*q1+L(-k)*l1*prop)*qbqb_b(2,1,1)
!            b212=(Q(-j)*q1+R(-j)*r1*prop)*qbqb_a(2,1,2)
!     .          +(Q(-k)*q1+L(-k)*r1*prop)*qbqb_b(2,1,2)
!
!
!            mqq(0,j,k)=half*faclo*(
!     .      +Dble(a111*Dconjg(b111))+Dble(a112*Dconjg(b112))
!     .      +Dble(a221*Dconjg(b221))+Dble(a222*Dconjg(b222)))*two/xn
!            mqq(1,j,k)=half*faclo*
!     .      (abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
!     .      +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)
!            mqq(2,j,k)=half*faclo*(
!     .      +abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
!     .      +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)
!            endif
!cC---q-qb case
!c         elseif ((j .gt. 0) .and. (k .lt. 0)) then
!c             if (j .ne. -k) then 
!c            a111=(Q(+j)*q1+L(+j)*l1*prop)*qRb_a(1,1,1)
!c     .          +(Q(-k)*q1+L(-k)*l1*prop)*qRb_b(1,1,1)
!c            a112=(Q(+j)*q1+L(+j)*r1*prop)*qRb_a(1,1,2)
!c     .          +(Q(-k)*q1+L(-k)*r1*prop)*qRb_b(1,1,2)
!c            a221=(Q(+j)*q1+R(+j)*l1*prop)*qRb_a(2,2,1)
!c     .          +(Q(-k)*q1+R(-k)*l1*prop)*qRb_b(2,2,1)
!c            a222=(Q(+j)*q1+R(+j)*r1*prop)*qRb_a(2,2,2)
!c     .          +(Q(-k)*q1+R(-k)*r1*prop)*qRb_b(2,2,2)
!c
!c            a121=(Q(+j)*q1+L(+j)*l1*prop)*qRb_a(1,2,1)
!c     .          +(Q(-k)*q1+R(-k)*l1*prop)*qRb_b(1,2,1)
!c            a122=(Q(+j)*q1+L(+j)*r1*prop)*qRb_a(1,2,2)
!c     .          +(Q(-k)*q1+R(-k)*r1*prop)*qRb_b(1,2,2)
!c            a211=(Q(+j)*q1+R(+j)*l1*prop)*qRb_a(2,1,1)
!c     .          +(Q(-k)*q1+L(-k)*l1*prop)*qRb_b(2,1,1)
!c            a212=(Q(+j)*q1+R(+j)*r1*prop)*qRb_a(2,1,2)
!c     .          +(Q(-k)*q1+L(-k)*r1*prop)*qRb_b(2,1,2)
!c            mqq(0,j,k)=zip
!c            mqq(1,j,k)=zip
!c            mqq(2,j,k)=
!c     .      +faclo*(abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
!c     .             +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)
!c
!c            elseif (j .eq. -k) then
!cc--case where final state from annihilation diagrams is the same quark
!c            a111=(Q(j)*q1+L(j)*l1*prop)*(qRb_a(1,1,1)+qRb_b(1,1,1))
!c            b111=(Q(j)*q1+L(j)*l1*prop)*(qqb_a(1,1,1)+qqb_b(1,1,1))
!c
!c            a112=(Q(j)*q1+L(j)*r1*prop)*(qRb_a(1,1,2)+qRb_b(1,1,2))
!c            b112=(Q(j)*q1+L(j)*r1*prop)*(qqb_a(1,1,2)+qqb_b(1,1,2))
!c
!c            a221=(Q(j)*q1+R(j)*l1*prop)*(qRb_a(2,2,1)+qRb_b(2,2,1))
!c            b221=(Q(j)*q1+R(j)*l1*prop)*(qqb_a(2,2,1)+qqb_b(2,2,1))
!c
!c            a222=(Q(j)*q1+R(j)*r1*prop)*(qRb_a(2,2,2)+qRb_b(2,2,2))
!c            b222=(Q(j)*q1+R(j)*r1*prop)*(qqb_a(2,2,2)+qqb_b(2,2,2))
!c
!c            a121=(Q(+j)*q1+L(+j)*l1*prop)*qRb_a(1,2,1)
!c     .          +(Q(-k)*q1+R(-k)*l1*prop)*qRb_b(1,2,1)
!c            a122=(Q(+j)*q1+L(+j)*r1*prop)*qRb_a(1,2,2)
!c     .          +(Q(-k)*q1+R(-k)*r1*prop)*qRb_b(1,2,2)
!c            a211=(Q(+j)*q1+R(+j)*l1*prop)*qRb_a(2,1,1)
!c     .          +(Q(-k)*q1+L(-k)*l1*prop)*qRb_b(2,1,1)
!c            a212=(Q(+j)*q1+R(+j)*r1*prop)*qRb_a(2,1,2)
!c     .          +(Q(-k)*q1+L(-k)*r1*prop)*qRb_b(2,1,2)
!c
!c            b121=(Q(+j)*q1+L(+j)*l1*prop)*qqb_a(1,2,1)
!c     .          +(Q(-k)*q1+R(-k)*l1*prop)*qqb_b(1,2,1)
!c            b122=(Q(+j)*q1+L(+j)*r1*prop)*qqb_a(1,2,2)
!c     .          +(Q(-k)*q1+R(-k)*r1*prop)*qqb_b(1,2,2)
!c            b211=(Q(+j)*q1+R(+j)*l1*prop)*qqb_a(2,1,1)
!c     .          +(Q(-k)*q1+L(-k)*l1*prop)*qqb_b(2,1,1)
!c            b212=(Q(+j)*q1+R(+j)*r1*prop)*qqb_a(2,1,2)
!c     .          +(Q(-k)*q1+L(-k)*r1*prop)*qqb_b(2,1,2)
!c
!c            mqq(0,j,k)=faclo*(
!c     .      +Dble(a111*Dconjg(b111))+Dble(a112*Dconjg(b112))
!c     .      +Dble(a221*Dconjg(b221))+Dble(a222*Dconjg(b222)))*two/xn
!c            mqq(1,j,k)=faclo*(
!c     .      +abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
!c     .      +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)
!c            mqq(2,j,k)=faclo*
!c     .      (abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
!c     .      +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)
!c
!c       if ((j.eq.1).or.(j.eq.3).or.(j.eq.5)) then
!c           nup=2
!c           ndo=nf-3
!c       else
!c           nup=1
!c           ndo=nf-2
!c       endif
!c       if (nf.le. 4) ndo=ndo-1
!c       if (nf.le. 3) nup=nup-1
!c            b111=(Q(+j)*q1+L(+j)*l1*prop)*qqb_a(1,1,1)
!c     .          +(Q(+1)*q1+L(+1)*l1*prop)*qqb_b(1,1,1)
!c            b112=(Q(+j)*q1+L(+j)*r1*prop)*qqb_a(1,1,2)
!c     .          +(Q(+1)*q1+L(+1)*r1*prop)*qqb_b(1,1,2)
!c            b221=(Q(+j)*q1+R(+j)*l1*prop)*qqb_a(2,2,1)
!c     .          +(Q(+1)*q1+R(+1)*l1*prop)*qqb_b(2,2,1)
!c            b222=(Q(+j)*q1+R(+j)*r1*prop)*qqb_a(2,2,2)
!c     .          +(Q(+1)*q1+R(+1)*r1*prop)*qqb_b(2,2,2)
!c            b121=(Q(+j)*q1+L(+j)*l1*prop)*qqb_a(1,2,1)
!c     .          +(Q(+1)*q1+R(+1)*l1*prop)*qqb_b(1,2,1)
!c            b122=(Q(+j)*q1+L(+j)*r1*prop)*qqb_a(1,2,2)
!c     .          +(Q(+1)*q1+R(+1)*r1*prop)*qqb_b(1,2,2)
!c            b211=(Q(+j)*q1+R(+j)*l1*prop)*qqb_a(2,1,1)
!c     .          +(Q(+1)*q1+L(+1)*l1*prop)*qqb_b(2,1,1)
!c            b212=(Q(+j)*q1+R(+j)*r1*prop)*qqb_a(2,1,2)
!c     .          +(Q(+1)*q1+L(+1)*r1*prop)*qqb_b(2,1,2)
!c            
!c      tdo=faclo*(abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
!c     .             +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)
!c
!c            b111=(Q(+j)*q1+L(+j)*l1*prop)*qqb_a(1,1,1)
!c     .          +(Q(+2)*q1+L(+2)*l1*prop)*qqb_b(1,1,1)
!c            b112=(Q(+j)*q1+L(+j)*r1*prop)*qqb_a(1,1,2)
!c     .          +(Q(+2)*q1+L(+2)*r1*prop)*qqb_b(1,1,2)
!c            b221=(Q(+j)*q1+R(+j)*l1*prop)*qqb_a(2,2,1)
!c     .          +(Q(+2)*q1+R(+2)*l1*prop)*qqb_b(2,2,1)
!c            b222=(Q(+j)*q1+R(+j)*r1*prop)*qqb_a(2,2,2)
!c     .          +(Q(+2)*q1+R(+2)*r1*prop)*qqb_b(2,2,2)
!c            b121=(Q(+j)*q1+L(+j)*l1*prop)*qqb_a(1,2,1)
!c     .          +(Q(+2)*q1+R(+2)*l1*prop)*qqb_b(1,2,1)
!c            b122=(Q(+j)*q1+L(+j)*r1*prop)*qqb_a(1,2,2)
!c     .          +(Q(+2)*q1+R(+2)*r1*prop)*qqb_b(1,2,2)
!c            b211=(Q(+j)*q1+R(+j)*l1*prop)*qqb_a(2,1,1)
!c     .          +(Q(+2)*q1+L(+2)*l1*prop)*qqb_b(2,1,1)
!c            b212=(Q(+j)*q1+R(+j)*r1*prop)*qqb_a(2,1,2)
!c     .          +(Q(+2)*q1+L(+2)*r1*prop)*qqb_b(2,1,2)
!c            
!c      tup=faclo*(abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
!c     .          +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)
!c
!c      mqq(1,j,k)=mqq(1,j,k)+dfloat(nup)*tup+dfloat(ndo)*tdo
!c      endif
!c      elseif ((j .lt. 0) .and. (k .gt. 0)) then
!cC---Qb-q case
!c            if (j .ne. -k) then
!c            a111=(Q(-j)*q1+L(-j)*l1*prop)*qbR_a(1,1,1)
!c     .          +(Q(+k)*q1+L(+k)*l1*prop)*qbR_b(1,1,1)
!c            a121=(Q(-j)*q1+L(-j)*l1*prop)*qbR_a(1,2,1)
!c     .          +(Q(+k)*q1+R(+k)*l1*prop)*qbR_b(1,2,1)
!c            a112=(Q(-j)*q1+L(-j)*r1*prop)*qbR_a(1,1,2)
!c     .          +(Q(+k)*q1+L(+k)*r1*prop)*qbR_b(1,1,2)
!c            a122=(Q(-j)*q1+L(-j)*r1*prop)*qbR_a(1,2,2)
!c     .          +(Q(+k)*q1+R(+k)*r1*prop)*qbR_b(1,2,2)
!c            a211=(Q(-j)*q1+R(-j)*l1*prop)*qbR_a(2,1,1)
!c     .          +(Q(+k)*q1+L(+k)*l1*prop)*qbR_b(2,1,1)
!c            a221=(Q(-j)*q1+R(-j)*l1*prop)*qbR_a(2,2,1)
!c     .          +(Q(+k)*q1+R(+k)*l1*prop)*qbR_b(2,2,1)
!c            a212=(Q(-j)*q1+R(-j)*r1*prop)*qbR_a(2,1,2)
!c     .          +(Q(+k)*q1+L(+k)*r1*prop)*qbR_b(2,1,2)
!c            a222=(Q(-j)*q1+R(-j)*r1*prop)*qbR_a(2,2,2)
!c     .          +(Q(+k)*q1+R(+k)*r1*prop)*qbR_b(2,2,2)
!c
!c            mqq(0,j,k)=zip
!c            mqq(1,j,k)=zip
!c            mqq(2,j,k)=
!c     .      +faclo*(abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
!c     .             +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)
!c            elseif (j .eq. -k) then
!c
!c            a111=(Q(-j)*q1+L(-j)*l1*prop)*(qbR_a(1,1,1)+qbR_b(1,1,1))
!c            b111=(Q(-j)*q1+L(-j)*l1*prop)*(qbq_a(1,1,1)+qbq_b(1,1,1))
!c            a112=(Q(-j)*q1+L(-j)*r1*prop)*(qbR_a(1,1,2)+qbR_b(1,1,2))
!c            b112=(Q(-j)*q1+L(-j)*r1*prop)*(qbq_a(1,1,2)+qbq_b(1,1,2))
!c            a221=(Q(-j)*q1+R(-j)*l1*prop)*(qbR_a(2,2,1)+qbR_b(2,2,1))
!c            b221=(Q(-j)*q1+R(-j)*l1*prop)*(qbq_a(2,2,1)+qbq_b(2,2,1))
!c            a222=(Q(-j)*q1+R(-j)*r1*prop)*(qbR_a(2,2,2)+qbR_b(2,2,2))
!c            b222=(Q(-j)*q1+R(-j)*r1*prop)*(qbq_a(2,2,2)+qbq_b(2,2,2))
!c
!c            a121=(Q(-j)*q1+L(-j)*l1*prop)*qbR_a(1,2,1)
!c     .          +(Q(+k)*q1+R(+k)*l1*prop)*qbR_b(1,2,1)
!c            a122=(Q(-j)*q1+L(-j)*r1*prop)*qbR_a(1,2,2)
!c     .          +(Q(+k)*q1+R(+k)*r1*prop)*qbR_b(1,2,2)
!c            a211=(Q(-j)*q1+R(-j)*l1*prop)*qbR_a(2,1,1)
!c     .          +(Q(+k)*q1+L(+k)*l1*prop)*qbR_b(2,1,1)
!c            a212=(Q(-j)*q1+R(-j)*r1*prop)*qbR_a(2,1,2)
!c     .          +(Q(+k)*q1+L(+k)*r1*prop)*qbR_b(2,1,2)
!c
!c            b121=(Q(-j)*q1+L(-j)*l1*prop)*qbq_a(1,2,1)
!c     .          +(Q(+k)*q1+R(+k)*l1*prop)*qbq_b(1,2,1)
!c            b122=(Q(-j)*q1+L(-j)*r1*prop)*qbq_a(1,2,2)
!c     .          +(Q(+k)*q1+R(+k)*r1*prop)*qbq_b(1,2,2)
!c            b211=(Q(-j)*q1+R(-j)*l1*prop)*qbq_a(2,1,1)
!c     .          +(Q(+k)*q1+L(+k)*l1*prop)*qbq_b(2,1,1)
!c            b212=(Q(-j)*q1+R(-j)*r1*prop)*qbq_a(2,1,2)
!c     .          +(Q(+k)*q1+L(+k)*r1*prop)*qbq_b(2,1,2)
!c
!c            mqq(0,j,k)=faclo*(
!c     .      +Dble(a111*Dconjg(b111))+Dble(a112*Dconjg(b112))
!c     .      +Dble(a221*Dconjg(b221))+Dble(a222*Dconjg(b222)))*two/xn
!c            mqq(2,j,k)=faclo*
!c     .      (abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
!c     .      +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)
!c            mqq(1,j,k)=faclo*(
!c     .      +abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
!c     .      +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)
!c
!cc--Here we must also add the contribution of other final state quarks
!cc  unequal to initial annihilating quarks
!c       if ((k.eq.1).or.(k.eq.3).or.(k.eq.5)) then
!c           nup=2
!c           ndo=nf-3
!c       else
!c           nup=1
!c           ndo=nf-2
!c       endif
!c       if (nf.le. 4) ndo=ndo-1
!c       if (nf.le. 3) nup=nup-1
!c            b111=(Q(-j)*q1+L(-j)*l1*prop)*qbq_a(1,1,1)
!c     .          +(Q(+3)*q1+L(+3)*l1*prop)*qbq_b(1,1,1)
!c            b112=(Q(-j)*q1+L(-j)*r1*prop)*qbq_a(1,1,2)
!c     .          +(Q(+3)*q1+L(+3)*r1*prop)*qbq_b(1,1,2)
!c            b221=(Q(-j)*q1+R(-j)*l1*prop)*qbq_a(2,2,1)
!c     .          +(Q(+3)*q1+R(+3)*l1*prop)*qbq_b(2,2,1)
!c            b222=(Q(-j)*q1+R(-j)*r1*prop)*qbq_a(2,2,2)
!c     .          +(Q(+3)*q1+R(+3)*r1*prop)*qbq_b(2,2,2)
!c            b121=(Q(-j)*q1+L(-j)*l1*prop)*qbq_a(1,2,1)
!c     .          +(Q(+3)*q1+R(+3)*l1*prop)*qbq_b(1,2,1)
!c            b122=(Q(-j)*q1+L(-j)*r1*prop)*qbq_a(1,2,2)
!c     .          +(Q(+3)*q1+R(+3)*r1*prop)*qbq_b(1,2,2)
!c            b211=(Q(-j)*q1+R(-j)*l1*prop)*qbq_a(2,1,1)
!c     .          +(Q(+3)*q1+L(+3)*l1*prop)*qbq_b(2,1,1)
!c            b212=(Q(-j)*q1+R(-j)*r1*prop)*qbq_a(2,1,2)
!c     .          +(Q(+3)*q1+L(+3)*r1*prop)*qbq_b(2,1,2)
!c      tdo=faclo*(abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
!c     .          +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)
!c
!c            b111=(Q(-j)*q1+L(-j)*l1*prop)*qbq_a(1,1,1)
!c     .          +(Q(+2)*q1+L(+2)*l1*prop)*qbq_b(1,1,1)
!c            b112=(Q(-j)*q1+L(-j)*r1*prop)*qbq_a(1,1,2)
!c     .          +(Q(+2)*q1+L(+2)*r1*prop)*qbq_b(1,1,2)
!c            b221=(Q(-j)*q1+R(-j)*l1*prop)*qbq_a(2,2,1)
!c     .          +(Q(+2)*q1+R(+2)*l1*prop)*qbq_b(2,2,1)
!c            b222=(Q(-j)*q1+R(-j)*r1*prop)*qbq_a(2,2,2)
!c     .          +(Q(+2)*q1+R(+2)*r1*prop)*qbq_b(2,2,2)
!c            b121=(Q(-j)*q1+L(-j)*l1*prop)*qbq_a(1,2,1)
!c     .          +(Q(+2)*q1+R(+2)*l1*prop)*qbq_b(1,2,1)
!c            b122=(Q(-j)*q1+L(-j)*r1*prop)*qbq_a(1,2,2)
!c     .          +(Q(+2)*q1+R(+2)*r1*prop)*qbq_b(1,2,2)
!c            b211=(Q(-j)*q1+R(-j)*l1*prop)*qbq_a(2,1,1)
!c     .          +(Q(+2)*q1+L(+2)*l1*prop)*qbq_b(2,1,1)
!c            b212=(Q(-j)*q1+R(-j)*r1*prop)*qbq_a(2,1,2)
!c     .          +(Q(+2)*q1+L(+2)*r1*prop)*qbq_b(2,1,2)
!c      tup=faclo*(abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
!c     .          +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)
!c
!c      mqq(1,j,k)=mqq(1,j,k)+dfloat(nup)*tup+dfloat(ndo)*tdo
!c
!c          endif
!          endif
!      msq(j,k)=msq(j,k)+mqq(0,j,k)+mqq(1,j,k)+mqq(2,j,k)
!      enddo
!      enddo
!      endif
!      return
      
      if (.not.Qflag) return

      do j=1,nf
         do k=1,nf
            
c----QQ case
         a111=(Q(j)*q1+L(j)*l1*prop)*qR_a(1,1,1)
     .          +(Q(k)*q1+L(k)*l1*prop)*qR_b(1,1,1)
         a121=(Q(j)*q1+L(j)*l1*prop)*qR_a(1,2,1)
     .          +(Q(k)*q1+R(k)*l1*prop)*qR_b(1,2,1)
         a112=(Q(j)*q1+L(j)*r1*prop)*qR_a(1,1,2)
     .          +(Q(k)*q1+L(k)*r1*prop)*qR_b(1,1,2)
         a122=(Q(j)*q1+L(j)*r1*prop)*qR_a(1,2,2)
     .          +(Q(k)*q1+R(k)*r1*prop)*qR_b(1,2,2)
         a211=(Q(j)*q1+R(j)*l1*prop)*qR_a(2,1,1)
     .          +(Q(k)*q1+L(k)*l1*prop)*qR_b(2,1,1)
         a221=(Q(j)*q1+R(j)*l1*prop)*qR_a(2,2,1)
     .          +(Q(k)*q1+R(k)*l1*prop)*qR_b(2,2,1)
         a212=(Q(j)*q1+R(j)*r1*prop)*qR_a(2,1,2)
     .          +(Q(k)*q1+L(k)*r1*prop)*qR_b(2,1,2)
         a222=(Q(j)*q1+R(j)*r1*prop)*qR_a(2,2,2)
     .          +(Q(k)*q1+R(k)*r1*prop)*qR_b(2,2,2)
         msq_1(j,k)=
     .       +faclo*(abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     .        +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)

c----QbQb case
         a111=(Q(j)*q1+L(j)*l1*prop)*qbRb_a(1,1,1)
     .          +(Q(k)*q1+L(k)*l1*prop)*qbRb_b(1,1,1)
         a121=(Q(j)*q1+L(j)*l1*prop)*qbRb_a(1,2,1)
     .          +(Q(k)*q1+R(k)*l1*prop)*qbRb_b(1,2,1)

         a112=(Q(j)*q1+L(j)*r1*prop)*qbRb_a(1,1,2)
     .          +(Q(k)*q1+L(k)*r1*prop)*qbRb_b(1,1,2)
         a122=(Q(j)*q1+L(j)*r1*prop)*qbRb_a(1,2,2)
     .          +(Q(k)*q1+R(k)*r1*prop)*qbRb_b(1,2,2)

         a211=(Q(j)*q1+R(j)*l1*prop)*qbRb_a(2,1,1)
     .          +(Q(k)*q1+L(k)*l1*prop)*qbRb_b(2,1,1)
         a221=(Q(j)*q1+R(j)*l1*prop)*qbRb_a(2,2,1)
     .          +(Q(k)*q1+R(k)*l1*prop)*qbRb_b(2,2,1)

         a212=(Q(j)*q1+R(j)*r1*prop)*qbRb_a(2,1,2)
     .          +(Q(k)*q1+L(k)*r1*prop)*qbRb_b(2,1,2)
         a222=(Q(j)*q1+R(j)*r1*prop)*qbRb_a(2,2,2)
     .          +(Q(k)*q1+R(k)*r1*prop)*qbRb_b(2,2,2)

         msq_2(j,k)=
     .        faclo*(abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     .        +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)


c---Q-Qb case
         a111=(Q(+j)*q1+L(+j)*l1*prop)*qRb_a(1,1,1)
     .          +(Q(k)*q1+L(k)*l1*prop)*qRb_b(1,1,1)
         a112=(Q(+j)*q1+L(+j)*r1*prop)*qRb_a(1,1,2)
     .          +(Q(k)*q1+L(k)*r1*prop)*qRb_b(1,1,2)
         a221=(Q(+j)*q1+R(+j)*l1*prop)*qRb_a(2,2,1)
     .          +(Q(k)*q1+R(k)*l1*prop)*qRb_b(2,2,1)
         a222=(Q(+j)*q1+R(+j)*r1*prop)*qRb_a(2,2,2)
     .          +(Q(k)*q1+R(k)*r1*prop)*qRb_b(2,2,2)

         a121=(Q(+j)*q1+L(+j)*l1*prop)*qRb_a(1,2,1)
     .          +(Q(k)*q1+R(k)*l1*prop)*qRb_b(1,2,1)
         a122=(Q(+j)*q1+L(+j)*r1*prop)*qRb_a(1,2,2)
     .          +(Q(k)*q1+R(k)*r1*prop)*qRb_b(1,2,2)
         a211=(Q(+j)*q1+R(+j)*l1*prop)*qRb_a(2,1,1)
     .          +(Q(k)*q1+L(k)*l1*prop)*qRb_b(2,1,1)
         a212=(Q(+j)*q1+R(+j)*r1*prop)*qRb_a(2,1,2)
     .          +(Q(k)*q1+L(k)*r1*prop)*qRb_b(2,1,2)

         msq_3(j,k)=
     .        faclo*(abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     .        +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)

C---Qb-q case
         a111=(Q(j)*q1+L(j)*l1*prop)*qbR_a(1,1,1)
     .          +(Q(+k)*q1+L(+k)*l1*prop)*qbR_b(1,1,1)
         a121=(Q(j)*q1+L(j)*l1*prop)*qbR_a(1,2,1)
     .          +(Q(+k)*q1+R(+k)*l1*prop)*qbR_b(1,2,1)
         a112=(Q(j)*q1+L(j)*r1*prop)*qbR_a(1,1,2)
     .          +(Q(+k)*q1+L(+k)*r1*prop)*qbR_b(1,1,2)
         a122=(Q(j)*q1+L(j)*r1*prop)*qbR_a(1,2,2)
     .          +(Q(+k)*q1+R(+k)*r1*prop)*qbR_b(1,2,2)
         a211=(Q(j)*q1+R(j)*l1*prop)*qbR_a(2,1,1)
     .          +(Q(+k)*q1+L(+k)*l1*prop)*qbR_b(2,1,1)
         a221=(Q(j)*q1+R(j)*l1*prop)*qbR_a(2,2,1)
     .          +(Q(+k)*q1+R(+k)*l1*prop)*qbR_b(2,2,1)
         a212=(Q(j)*q1+R(j)*r1*prop)*qbR_a(2,1,2)
     .          +(Q(+k)*q1+L(+k)*r1*prop)*qbR_b(2,1,2)
         a222=(Q(j)*q1+R(j)*r1*prop)*qbR_a(2,2,2)
     .          +(Q(+k)*q1+R(+k)*r1*prop)*qbR_b(2,2,2)

         msq_4(j,k)=
     .        faclo*(abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     .        +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)
      enddo
      enddo

cc    do not add diagonal pieces
      do j=1,nf
         msq_3(j,j)=0d0
         msq_4(j,j)=0d0
      enddo

      do j=1,nf
         do k=1,nf
            msq(j,k)=msq_1(j,k)
            msq(-j,-k)=msq_2(j,k)
            msq(j,-k)=msq(j,-k)+msq_3(j,k)
            msq(-j,k)=msq(-j,k)+msq_4(j,k)
         enddo
      enddo


c     diagonal pieces      
c----QQ case
      do j=1,nf
         k=j
         a111=(Q(j)*q1+L(j)*l1*prop)*(qR_a(1,1,1)+qR_b(1,1,1))
         b111=(Q(j)*q1+L(j)*l1*prop)*(qq_a(1,1,1)+qq_b(1,1,1))
         a112=(Q(j)*q1+L(j)*r1*prop)*(qR_a(1,1,2)+qR_b(1,1,2))
         b112=(Q(j)*q1+L(j)*r1*prop)*(qq_a(1,1,2)+qq_b(1,1,2))
         a221=(Q(j)*q1+R(j)*l1*prop)*(qR_a(2,2,1)+qR_b(2,2,1))
         b221=(Q(j)*q1+R(j)*l1*prop)*(qq_a(2,2,1)+qq_b(2,2,1))
         a222=(Q(j)*q1+R(j)*r1*prop)*(qR_a(2,2,2)+qR_b(2,2,2))
         b222=(Q(j)*q1+R(j)*r1*prop)*(qq_a(2,2,2)+qq_b(2,2,2))

         a121=(Q(j)*q1+L(j)*l1*prop)*qR_a(1,2,1)
     .          +(Q(k)*q1+R(k)*l1*prop)*qR_b(1,2,1)
         b121=(Q(j)*q1+L(j)*l1*prop)*qq_a(1,2,1)
     .          +(Q(k)*q1+R(k)*l1*prop)*qq_b(1,2,1)
         a122=(Q(j)*q1+L(j)*r1*prop)*qR_a(1,2,2)
     .          +(Q(k)*q1+R(k)*r1*prop)*qR_b(1,2,2)
         b122=(Q(j)*q1+L(j)*r1*prop)*qq_a(1,2,2)
     .          +(Q(k)*q1+R(k)*r1*prop)*qq_b(1,2,2)
         a211=(Q(j)*q1+R(j)*l1*prop)*qR_a(2,1,1)
     .          +(Q(k)*q1+L(k)*l1*prop)*qR_b(2,1,1)
         b211=(Q(j)*q1+R(j)*l1*prop)*qq_a(2,1,1)
     .          +(Q(k)*q1+L(k)*l1*prop)*qq_b(2,1,1)
         a212=(Q(j)*q1+R(j)*r1*prop)*qR_a(2,1,2)
     .          +(Q(k)*q1+L(k)*r1*prop)*qR_b(2,1,2)
         b212=(Q(j)*q1+R(j)*r1*prop)*qq_a(2,1,2)
     .          +(Q(k)*q1+L(k)*r1*prop)*qq_b(2,1,2)

         mqq0=half*faclo*(
     .        +Dble(a111*Dconjg(b111))
     .        +Dble(a112*Dconjg(b112))
     .        +Dble(a221*Dconjg(b221))
     .        +Dble(a222*Dconjg(b222)))
     .        *two/xn
         mqq1=half*faclo*
     .        (abs(a111)**2+abs(a112)**2
     .        +abs(a221)**2+abs(a222)**2
     .        +abs(a122)**2+abs(a212)**2
     .        +abs(a121)**2+abs(a211)**2)
         mqq2=half*faclo*(
     .        +abs(b111)**2+abs(b112)**2
     .        +abs(b221)**2+abs(b222)**2
     .        +abs(b122)**2+abs(b212)**2
     .        +abs(b121)**2+abs(b211)**2)
         msqd_1(j)=mqq0+mqq1+mqq2
      enddo

c----QbQb case
      do j=1,nf
         k=j
         a111=(Q(j)*q1+L(j)*l1*prop)*(qbRb_a(1,1,1)+qbRb_b(1,1,1))
         b111=(Q(j)*q1+L(j)*l1*prop)*(qbqb_a(1,1,1)+qbqb_b(1,1,1))
         a112=(Q(j)*q1+L(j)*r1*prop)*(qbRb_a(1,1,2)+qbRb_b(1,1,2))
         b112=(Q(j)*q1+L(j)*r1*prop)*(qbqb_a(1,1,2)+qbqb_b(1,1,2))
         a221=(Q(j)*q1+R(j)*l1*prop)*(qbRb_a(2,2,1)+qbRb_b(2,2,1))
         b221=(Q(j)*q1+R(j)*l1*prop)*(qbqb_a(2,2,1)+qbqb_b(2,2,1))
         a222=(Q(j)*q1+R(j)*r1*prop)*(qbRb_a(2,2,2)+qbRb_b(2,2,2))
         b222=(Q(j)*q1+R(j)*r1*prop)*(qbqb_a(2,2,2)+qbqb_b(2,2,2))


         a121=(Q(j)*q1+L(j)*l1*prop)*qbRb_a(1,2,1)
     .          +(Q(k)*q1+R(k)*l1*prop)*qbRb_b(1,2,1)
         a122=(Q(j)*q1+L(j)*r1*prop)*qbRb_a(1,2,2)
     .          +(Q(k)*q1+R(k)*r1*prop)*qbRb_b(1,2,2)
         a211=(Q(j)*q1+R(j)*l1*prop)*qbRb_a(2,1,1)
     .          +(Q(k)*q1+L(k)*l1*prop)*qbRb_b(2,1,1)
         a212=(Q(j)*q1+R(j)*r1*prop)*qbRb_a(2,1,2)
     .          +(Q(k)*q1+L(k)*r1*prop)*qbRb_b(2,1,2)

         b121=(Q(j)*q1+L(j)*l1*prop)*qbqb_a(1,2,1)
     .          +(Q(k)*q1+R(k)*l1*prop)*qbqb_b(1,2,1)
         b122=(Q(j)*q1+L(j)*r1*prop)*qbqb_a(1,2,2)
     .          +(Q(k)*q1+R(k)*r1*prop)*qbqb_b(1,2,2)
         b211=(Q(j)*q1+R(j)*l1*prop)*qbqb_a(2,1,1)
     .          +(Q(k)*q1+L(k)*l1*prop)*qbqb_b(2,1,1)
         b212=(Q(j)*q1+R(j)*r1*prop)*qbqb_a(2,1,2)
     .          +(Q(k)*q1+L(k)*r1*prop)*qbqb_b(2,1,2)


         mqq0=half*faclo*
     .        (dble(a111*dconjg(b111))
     .        +dble(a112*dconjg(b112))
     .        +dble(a221*dconjg(b221))
     .        +dble(a222*dconjg(b222)))*two/xn
         mqq1=half*faclo*
     .        (abs(a111)**2+abs(a112)**2
     .        +abs(a221)**2+abs(a222)**2
     .        +abs(a122)**2+abs(a212)**2
     .        +abs(a121)**2+abs(a211)**2)
         mqq2=half*faclo*
     .        (abs(b111)**2+abs(b112)**2
     .        +abs(b221)**2+abs(b222)**2
     .        +abs(b122)**2+abs(b212)**2
     .        +abs(b121)**2+abs(b211)**2)
         msqd_2(j)=mqq0+mqq1+mqq2
      enddo

      
      do j=1,nf
         k=-j
c--case where final state from annihilation diagrams is the same quark
         a111=(Q(j)*q1+L(j)*l1*prop)*(qRb_a(1,1,1)+qRb_b(1,1,1))
         b111=(Q(j)*q1+L(j)*l1*prop)*(qqb_a(1,1,1)+qqb_b(1,1,1))

         a112=(Q(j)*q1+L(j)*r1*prop)*(qRb_a(1,1,2)+qRb_b(1,1,2))
         b112=(Q(j)*q1+L(j)*r1*prop)*(qqb_a(1,1,2)+qqb_b(1,1,2))

         a221=(Q(j)*q1+R(j)*l1*prop)*(qRb_a(2,2,1)+qRb_b(2,2,1))
         b221=(Q(j)*q1+R(j)*l1*prop)*(qqb_a(2,2,1)+qqb_b(2,2,1))

         a222=(Q(j)*q1+R(j)*r1*prop)*(qRb_a(2,2,2)+qRb_b(2,2,2))
         b222=(Q(j)*q1+R(j)*r1*prop)*(qqb_a(2,2,2)+qqb_b(2,2,2))

         a121=(Q(+j)*q1+L(+j)*l1*prop)*qRb_a(1,2,1)
     .          +(Q(-k)*q1+R(-k)*l1*prop)*qRb_b(1,2,1)
         a122=(Q(+j)*q1+L(+j)*r1*prop)*qRb_a(1,2,2)
     .          +(Q(-k)*q1+R(-k)*r1*prop)*qRb_b(1,2,2)
         a211=(Q(+j)*q1+R(+j)*l1*prop)*qRb_a(2,1,1)
     .          +(Q(-k)*q1+L(-k)*l1*prop)*qRb_b(2,1,1)
         a212=(Q(+j)*q1+R(+j)*r1*prop)*qRb_a(2,1,2)
     .          +(Q(-k)*q1+L(-k)*r1*prop)*qRb_b(2,1,2)

         b121=(Q(+j)*q1+L(+j)*l1*prop)*qqb_a(1,2,1)
     .          +(Q(-k)*q1+R(-k)*l1*prop)*qqb_b(1,2,1)
         b122=(Q(+j)*q1+L(+j)*r1*prop)*qqb_a(1,2,2)
     .          +(Q(-k)*q1+R(-k)*r1*prop)*qqb_b(1,2,2)
         b211=(Q(+j)*q1+R(+j)*l1*prop)*qqb_a(2,1,1)
     .          +(Q(-k)*q1+L(-k)*l1*prop)*qqb_b(2,1,1)
         b212=(Q(+j)*q1+R(+j)*r1*prop)*qqb_a(2,1,2)
     .          +(Q(-k)*q1+L(-k)*r1*prop)*qqb_b(2,1,2)

         mqq0=faclo*(
     .      +Dble(a111*Dconjg(b111))+Dble(a112*Dconjg(b112))
     .      +Dble(a221*Dconjg(b221))+Dble(a222*Dconjg(b222)))*two/xn
         mqq1=faclo*(
     .      +abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
     .      +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)
         mqq2=faclo*
     .      (abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     .      +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)

         nup=mod(k,2)+1
         ndo=nf-2-mod(k,2)
         if (nf.le. 4) ndo=ndo-1
         if (nf.le. 3) nup=nup-1
         
          b111=(Q(+j)*q1+L(+j)*l1*prop)*qqb_a(1,1,1)
     .          +(Q(+1)*q1+L(+1)*l1*prop)*qqb_b(1,1,1)
          b112=(Q(+j)*q1+L(+j)*r1*prop)*qqb_a(1,1,2)
     .          +(Q(+1)*q1+L(+1)*r1*prop)*qqb_b(1,1,2)
          b221=(Q(+j)*q1+R(+j)*l1*prop)*qqb_a(2,2,1)
     .          +(Q(+1)*q1+R(+1)*l1*prop)*qqb_b(2,2,1)
          b222=(Q(+j)*q1+R(+j)*r1*prop)*qqb_a(2,2,2)
     .          +(Q(+1)*q1+R(+1)*r1*prop)*qqb_b(2,2,2)
          b121=(Q(+j)*q1+L(+j)*l1*prop)*qqb_a(1,2,1)
     .          +(Q(+1)*q1+R(+1)*l1*prop)*qqb_b(1,2,1)
          b122=(Q(+j)*q1+L(+j)*r1*prop)*qqb_a(1,2,2)
     .          +(Q(+1)*q1+R(+1)*r1*prop)*qqb_b(1,2,2)
          b211=(Q(+j)*q1+R(+j)*l1*prop)*qqb_a(2,1,1)
     .          +(Q(+1)*q1+L(+1)*l1*prop)*qqb_b(2,1,1)
          b212=(Q(+j)*q1+R(+j)*r1*prop)*qqb_a(2,1,2)
     .          +(Q(+1)*q1+L(+1)*r1*prop)*qqb_b(2,1,2)
            
      tdo=faclo*(abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
     .             +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)

          b111=(Q(+j)*q1+L(+j)*l1*prop)*qqb_a(1,1,1)
     .          +(Q(+2)*q1+L(+2)*l1*prop)*qqb_b(1,1,1)
          b112=(Q(+j)*q1+L(+j)*r1*prop)*qqb_a(1,1,2)
     .          +(Q(+2)*q1+L(+2)*r1*prop)*qqb_b(1,1,2)
          b221=(Q(+j)*q1+R(+j)*l1*prop)*qqb_a(2,2,1)
     .          +(Q(+2)*q1+R(+2)*l1*prop)*qqb_b(2,2,1)
          b222=(Q(+j)*q1+R(+j)*r1*prop)*qqb_a(2,2,2)
     .          +(Q(+2)*q1+R(+2)*r1*prop)*qqb_b(2,2,2)
          b121=(Q(+j)*q1+L(+j)*l1*prop)*qqb_a(1,2,1)
     .          +(Q(+2)*q1+R(+2)*l1*prop)*qqb_b(1,2,1)
          b122=(Q(+j)*q1+L(+j)*r1*prop)*qqb_a(1,2,2)
     .          +(Q(+2)*q1+R(+2)*r1*prop)*qqb_b(1,2,2)
          b211=(Q(+j)*q1+R(+j)*l1*prop)*qqb_a(2,1,1)
     .          +(Q(+2)*q1+L(+2)*l1*prop)*qqb_b(2,1,1)
          b212=(Q(+j)*q1+R(+j)*r1*prop)*qqb_a(2,1,2)
     .          +(Q(+2)*q1+L(+2)*r1*prop)*qqb_b(2,1,2)
            
      tup=faclo*(abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
     .          +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)

          mqq1=mqq1+dfloat(nup)*tup+dfloat(ndo)*tdo
          msqd_3(j)=mqq0+mqq1+mqq2
       enddo
       
      do k=1,nf
         j=-k
         a111=(Q(-j)*q1+L(-j)*l1*prop)*(qbR_a(1,1,1)+qbR_b(1,1,1))
         b111=(Q(-j)*q1+L(-j)*l1*prop)*(qbq_a(1,1,1)+qbq_b(1,1,1))
         a112=(Q(-j)*q1+L(-j)*r1*prop)*(qbR_a(1,1,2)+qbR_b(1,1,2))
         b112=(Q(-j)*q1+L(-j)*r1*prop)*(qbq_a(1,1,2)+qbq_b(1,1,2))
         a221=(Q(-j)*q1+R(-j)*l1*prop)*(qbR_a(2,2,1)+qbR_b(2,2,1))
         b221=(Q(-j)*q1+R(-j)*l1*prop)*(qbq_a(2,2,1)+qbq_b(2,2,1))
         a222=(Q(-j)*q1+R(-j)*r1*prop)*(qbR_a(2,2,2)+qbR_b(2,2,2))
         b222=(Q(-j)*q1+R(-j)*r1*prop)*(qbq_a(2,2,2)+qbq_b(2,2,2))

         a121=(Q(-j)*q1+L(-j)*l1*prop)*qbR_a(1,2,1)
     .          +(Q(+k)*q1+R(+k)*l1*prop)*qbR_b(1,2,1)
         a122=(Q(-j)*q1+L(-j)*r1*prop)*qbR_a(1,2,2)
     .          +(Q(+k)*q1+R(+k)*r1*prop)*qbR_b(1,2,2)
         a211=(Q(-j)*q1+R(-j)*l1*prop)*qbR_a(2,1,1)
     .          +(Q(+k)*q1+L(+k)*l1*prop)*qbR_b(2,1,1)
         a212=(Q(-j)*q1+R(-j)*r1*prop)*qbR_a(2,1,2)
     .          +(Q(+k)*q1+L(+k)*r1*prop)*qbR_b(2,1,2)

         b121=(Q(-j)*q1+L(-j)*l1*prop)*qbq_a(1,2,1)
     .          +(Q(+k)*q1+R(+k)*l1*prop)*qbq_b(1,2,1)
         b122=(Q(-j)*q1+L(-j)*r1*prop)*qbq_a(1,2,2)
     .          +(Q(+k)*q1+R(+k)*r1*prop)*qbq_b(1,2,2)
         b211=(Q(-j)*q1+R(-j)*l1*prop)*qbq_a(2,1,1)
     .          +(Q(+k)*q1+L(+k)*l1*prop)*qbq_b(2,1,1)
         b212=(Q(-j)*q1+R(-j)*r1*prop)*qbq_a(2,1,2)
     .          +(Q(+k)*q1+L(+k)*r1*prop)*qbq_b(2,1,2)

         mqq0=faclo*(
     .      +Dble(a111*Dconjg(b111))+Dble(a112*Dconjg(b112))
     .      +Dble(a221*Dconjg(b221))+Dble(a222*Dconjg(b222)))*two/xn
         mqq2=faclo*
     .      (abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     .      +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)
         mqq1=faclo*(
     .      +abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
     .      +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)

c--Here we must also add the contribution of other final state quarks
c  unequal to initial annihilating quarks
         nup=mod(k,2)+1
         ndo=nf-2-mod(k,2)
         if (nf.le. 4) ndo=ndo-1
         if (nf.le. 3) nup=nup-1
         
          b111=(Q(-j)*q1+L(-j)*l1*prop)*qbq_a(1,1,1)
     .          +(Q(+3)*q1+L(+3)*l1*prop)*qbq_b(1,1,1)
          b112=(Q(-j)*q1+L(-j)*r1*prop)*qbq_a(1,1,2)
     .          +(Q(+3)*q1+L(+3)*r1*prop)*qbq_b(1,1,2)
          b221=(Q(-j)*q1+R(-j)*l1*prop)*qbq_a(2,2,1)
     .          +(Q(+3)*q1+R(+3)*l1*prop)*qbq_b(2,2,1)
          b222=(Q(-j)*q1+R(-j)*r1*prop)*qbq_a(2,2,2)
     .          +(Q(+3)*q1+R(+3)*r1*prop)*qbq_b(2,2,2)
          b121=(Q(-j)*q1+L(-j)*l1*prop)*qbq_a(1,2,1)
     .          +(Q(+3)*q1+R(+3)*l1*prop)*qbq_b(1,2,1)
          b122=(Q(-j)*q1+L(-j)*r1*prop)*qbq_a(1,2,2)
     .          +(Q(+3)*q1+R(+3)*r1*prop)*qbq_b(1,2,2)
          b211=(Q(-j)*q1+R(-j)*l1*prop)*qbq_a(2,1,1)
     .          +(Q(+3)*q1+L(+3)*l1*prop)*qbq_b(2,1,1)
          b212=(Q(-j)*q1+R(-j)*r1*prop)*qbq_a(2,1,2)
     .          +(Q(+3)*q1+L(+3)*r1*prop)*qbq_b(2,1,2)
      tdo=faclo*(abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
     .          +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)

          b111=(Q(-j)*q1+L(-j)*l1*prop)*qbq_a(1,1,1)
     .          +(Q(+2)*q1+L(+2)*l1*prop)*qbq_b(1,1,1)
          b112=(Q(-j)*q1+L(-j)*r1*prop)*qbq_a(1,1,2)
     .          +(Q(+2)*q1+L(+2)*r1*prop)*qbq_b(1,1,2)
          b221=(Q(-j)*q1+R(-j)*l1*prop)*qbq_a(2,2,1)
     .          +(Q(+2)*q1+R(+2)*l1*prop)*qbq_b(2,2,1)
          b222=(Q(-j)*q1+R(-j)*r1*prop)*qbq_a(2,2,2)
     .          +(Q(+2)*q1+R(+2)*r1*prop)*qbq_b(2,2,2)
          b121=(Q(-j)*q1+L(-j)*l1*prop)*qbq_a(1,2,1)
     .          +(Q(+2)*q1+R(+2)*l1*prop)*qbq_b(1,2,1)
          b122=(Q(-j)*q1+L(-j)*r1*prop)*qbq_a(1,2,2)
     .          +(Q(+2)*q1+R(+2)*r1*prop)*qbq_b(1,2,2)
          b211=(Q(-j)*q1+R(-j)*l1*prop)*qbq_a(2,1,1)
     .          +(Q(+2)*q1+L(+2)*l1*prop)*qbq_b(2,1,1)
          b212=(Q(-j)*q1+R(-j)*r1*prop)*qbq_a(2,1,2)
     .          +(Q(+2)*q1+L(+2)*r1*prop)*qbq_b(2,1,2)
      tup=faclo*(abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
     .          +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)

          mqq1=mqq1+dfloat(nup)*tup+dfloat(ndo)*tdo
          msqd_4(k)=mqq0+mqq1+mqq2
      enddo

       do j=1,nf
          msq(j,j)=msqd_1(j)
          msq(-j,-j)=msqd_2(j)
          msq(j,-j)=msq(j,-j)+msqd_3(j)
          msq(-j,j)=msq(-j,j)+msqd_4(j)
       enddo
      
      return
      end
          
    

      
     
