      subroutine qqb_z1jet_gs(p,msq)
************************************************************************
*     Author: R.K. Ellis                                               *
*     September, 1999.                                                 *
*    Matrix element SUBTRACTION squared averag'd over init'l colors    *
*    and spins                                                         *
c     q(-p1)+qbar(-p2) -->  Z + parton(p5) + parton(p6)                *
c                           |                                          *
c                            -->l(p3)+a(p4)                            *
************************************************************************

      implicit none 
      include 'constants.f'
      include 'ptilde.f'
      include 'qqgg.f'
      integer j,k,nd
c --- remember: nd will count the dipoles
      
      double precision p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      double precision 
     & msq15_2(-nf:nf,-nf:nf),msq25_1(-nf:nf,-nf:nf),
     & msq16_2(-nf:nf,-nf:nf),msq26_1(-nf:nf,-nf:nf),
     & msq15_6(-nf:nf,-nf:nf),msq26_5(-nf:nf,-nf:nf),
     & msq16_5(-nf:nf,-nf:nf),msq25_6(-nf:nf,-nf:nf),
     & msq56_1v(-nf:nf,-nf:nf),msq56_2v(-nf:nf,-nf:nf),
     & msq26_5v(-nf:nf,-nf:nf),msq26_1v(-nf:nf,-nf:nf),
     & msq15_6v(-nf:nf,-nf:nf),msq16_2v(-nf:nf,-nf:nf),
     & msq16_5v(-nf:nf,-nf:nf),msq25_6v(-nf:nf,-nf:nf),
     & msq25_1v(-nf:nf,-nf:nf),
     & msq15_2v(-nf:nf,-nf:nf),
     & dummy(-nf:nf,-nf:nf),
     & sub15_2(4),sub25_1(4),sub16_2(4),sub26_1(4),
     & sub15_6(4),sub16_5(4),sub25_6(4),sub26_5(4),
     & sub56_1(4),sub56_2(4),sub56_1v,sub56_2v,
     & sub26_5v,sub25_1v,sub26_1v,sub16_5v,sub16_2v,sub15_2v,sub15_6v,
     & sub25_6v
      external qqb_z1jet,qqb_z_gvec
      ndmax=6

c--- calculate all the initial-initial dipoles
      call dips(1,p,1,5,2,sub15_2,sub15_2v,msq15_2,msq15_2v,
     . qqb_z1jet,qqb_z_gvec)
      call dips(2,p,2,5,1,sub25_1,sub25_1v,msq25_1,msq25_1v,
     . qqb_z1jet,qqb_z_gvec)
      call dips(3,p,1,6,2,sub16_2,sub16_2v,msq16_2,msq16_2v,
     . qqb_z1jet,qqb_z_gvec)
      call dips(4,p,2,6,1,sub26_1,sub26_1v,msq26_1,msq26_1v,
     . qqb_z1jet,qqb_z_gvec)

c--- now the basic initial final ones
      call dips(5,p,1,5,6,sub15_6,sub15_6v,msq15_6,msq15_6v,
     . qqb_z1jet,qqb_z_gvec)
c--- called for final initial the routine only supplies new values for
c--- sub... and sub...v and msqv
      call dips(5,p,5,6,1,sub56_1,sub56_1v,dummy,msq56_1v,
     . qqb_z1jet,qqb_z_gvec)
      call dips(5,p,1,6,5,sub16_5,sub16_5v,msq16_5,msq16_5v,
     . qqb_z1jet,qqb_z_gvec)

      call dips(6,p,2,6,5,sub26_5,sub26_5v,msq26_5,msq26_5v,
     . qqb_z1jet,qqb_z_gvec)
      call dips(6,p,5,6,2,sub56_2,sub56_2v,dummy,msq56_2v,
     . qqb_z1jet,qqb_z_gvec)
      call dips(6,p,2,5,6,sub25_6,sub25_6v,msq25_6,msq25_6v,
     . qqb_z1jet,qqb_z_gvec)

      msq=0d0
c      do j=-nf,nf
c      do k=-nf,nf      
c      do nd=1,ndmax
c        msq(nd,j,k)=0d0
c      enddo
c      enddo
c      enddo
      
cc--- Gflag subtraction pieces
c      do j=-nf,nf
c      do k=-nf,nf
c      
c      if ((j .ne. 0) .and. (k .ne. 0) .and. (j.ne.-k)) goto 19
c
cc--- do only q-qb and qb-q cases      
c      if (  ((j .gt. 0).and.(k .lt. 0))
c     . .or. ((j .lt. 0).and.(k .gt. 0))) then
cC-----half=statistical factor
c      msq(1,j,k)=-half*msq15_2(j,k)*sub15_2(qq)/xn
c      msq(2,j,k)=-half*msq25_1(j,k)*sub25_1(qq)/xn
c      msq(3,j,k)=-half*msq16_2(j,k)*sub16_2(qq)/xn
c      msq(4,j,k)=-half*msq26_1(j,k)*sub26_1(qq)/xn
c      msq(5,j,k)=half*xn*(
c     .  msq15_6(j,k)*(sub15_6(qq)+0.5d0*sub56_1(gg))
c     . +0.5d0*msq56_1v(j,k)*sub56_1v
c     . +msq16_5(j,k)*(sub16_5(qq)+0.5d0*sub56_1(gg))
c     . +0.5d0*msq56_1v(j,k)*sub56_1v)
c      msq(6,j,k)=half*xn*(
c     .  msq26_5(j,k)*(sub26_5(qq)+0.5d0*sub56_2(gg))
c     . +0.5d0*msq56_2v(j,k)*sub56_2v
c     . +msq25_6(j,k)*(sub25_6(qq)+0.5d0*sub56_2(gg))
c     . +0.5d0*msq56_2v(j,k)*sub56_2v)
c      elseif ((k .eq. 0).and.(j.ne.0)) then
cc--- q-g and qb-g cases
c      msq(2,j,k)=2d0*tr*msq25_1(j,-j)*sub25_1(qg)
c      msq(3,j,k)=xn*msq16_2(j,k)*sub16_2(qq)
c      msq(4,j,k)=xn*(msq26_1(j,k)*sub26_1(gg)+msq26_1v(j,k)*sub26_1v)
c      msq(5,j,k)=-(msq16_5(j,k)*sub16_5(qq)+msq16_5(j,k)*sub56_1(qq))/xn
c      msq(6,j,k)=xn*(msq26_5(j,k)*sub26_5(gg)+msq26_5v(j,k)*sub26_5v
c     .              +msq26_5(j,k)*sub56_2(qq))
c 
c      elseif ((j .eq. 0).and.(k.ne.0)) then
cc--- g-q and g-qb cases
c      msq(1,j,k)=2d0*tr*msq15_2(-k,k)*sub15_2(qg)
c      msq(3,j,k)=xn*(msq16_2(j,k)*sub16_2(gg)+msq16_2v(j,k)*sub16_2v)
c      msq(4,j,k)=xn*msq26_1(j,k)*sub26_1(qq)
c      msq(5,j,k)=xn*(msq16_5(j,k)*sub16_5(gg)+msq16_5v(j,k)*sub16_5v
c     .              +msq15_6(j,k)*sub56_1(qq))
c      msq(6,j,k)=-(msq26_5(j,k)*sub26_5(qq)+msq26_5(j,k)*sub56_2(qq))/xn
c
c      elseif ((j .eq. 0).and.(k .eq. 0)) then
cc--- g-g case (real process is g(p1)+g(p2) --> qb(p5)+q(p6)
cc---Hence 15 split multiplies q(15)+g(p2)-->Z+q(p6)
cc---Hence 25 split multiplies g(p1)+q(p25)-->Z+q(p6)
c      msq(1,j,k)=(msq15_2(+1,k)+msq15_2(+2,k)+msq15_2(+3,k)
c     .           +msq15_2(+4,k)+msq15_2(+5,k))*sub15_2(qg)*2d0*tr
c      msq(2,j,k)=(msq25_1(k,+1)+msq25_1(k,+2)+msq25_1(k,+3)
c     .           +msq25_1(k,+4)+msq25_1(k,+5))*sub25_1(qg)*2d0*tr
c      msq(3,j,k)=(msq16_2(-5,k)+msq16_2(-4,k)+msq16_2(-3,k)
c     .           +msq16_2(-2,k)+msq16_2(-1,k))*sub16_2(qg)*2d0*tr
c      msq(4,j,k)=(msq26_1(k,-5)+msq26_1(k,-4)+msq26_1(k,-3)
c     .           +msq26_1(k,-2)+msq26_1(k,-1))*sub26_1(qg)*2d0*tr
c
c      endif
c
c 19   continue
c      enddo
c      enddo



c--- Gflag subtraction pieces

c---  do only q-qb and qb-q cases
      do j=1,nf
         k=-j
C-----half=statistical factor
         msq(1,j,k)=-half*msq15_2(j,k)*sub15_2(qq)/xn
         msq(2,j,k)=-half*msq25_1(j,k)*sub25_1(qq)/xn
         msq(3,j,k)=-half*msq16_2(j,k)*sub16_2(qq)/xn
         msq(4,j,k)=-half*msq26_1(j,k)*sub26_1(qq)/xn
         msq(5,j,k)=half*xn*(
     .        msq15_6(j,k)*(sub15_6(qq)+0.5d0*sub56_1(gg))
     .        +0.5d0*msq56_1v(j,k)*sub56_1v
     .        +msq16_5(j,k)*(sub16_5(qq)+0.5d0*sub56_1(gg))
     .        +0.5d0*msq56_1v(j,k)*sub56_1v)
         msq(6,j,k)=half*xn*(
     .        msq26_5(j,k)*(sub26_5(qq)+0.5d0*sub56_2(gg))
     .        +0.5d0*msq56_2v(j,k)*sub56_2v
     .        +msq25_6(j,k)*(sub25_6(qq)+0.5d0*sub56_2(gg))
     .        +0.5d0*msq56_2v(j,k)*sub56_2v)
      enddo
      do j=-nf,-1
         k=-j
C-----half=statistical factor
         msq(1,j,k)=-half*msq15_2(j,k)*sub15_2(qq)/xn
         msq(2,j,k)=-half*msq25_1(j,k)*sub25_1(qq)/xn
         msq(3,j,k)=-half*msq16_2(j,k)*sub16_2(qq)/xn
         msq(4,j,k)=-half*msq26_1(j,k)*sub26_1(qq)/xn
         msq(5,j,k)=half*xn*(
     .        msq15_6(j,k)*(sub15_6(qq)+0.5d0*sub56_1(gg))
     .        +0.5d0*msq56_1v(j,k)*sub56_1v
     .        +msq16_5(j,k)*(sub16_5(qq)+0.5d0*sub56_1(gg))
     .        +0.5d0*msq56_1v(j,k)*sub56_1v)
         msq(6,j,k)=half*xn*(
     .        msq26_5(j,k)*(sub26_5(qq)+0.5d0*sub56_2(gg))
     .        +0.5d0*msq56_2v(j,k)*sub56_2v
     .        +msq25_6(j,k)*(sub25_6(qq)+0.5d0*sub56_2(gg))
     .        +0.5d0*msq56_2v(j,k)*sub56_2v)
      enddo
      
      k=0
      do j=1,nf
c--- q-g and qb-g cases
      msq(2,j,k)=2d0*tr*msq25_1(j,-j)*sub25_1(qg)
      msq(3,j,k)=xn*msq16_2(j,k)*sub16_2(qq)
      msq(4,j,k)=xn*(msq26_1(j,k)*sub26_1(gg)+msq26_1v(j,k)*sub26_1v)
      msq(5,j,k)=-(msq16_5(j,k)*sub16_5(qq)+msq16_5(j,k)*sub56_1(qq))/xn
      msq(6,j,k)=xn*(msq26_5(j,k)*sub26_5(gg)+msq26_5v(j,k)*sub26_5v
     .              +msq26_5(j,k)*sub56_2(qq))
      enddo
      do j=-nf,-1
c--- q-g and qb-g cases
      msq(2,j,k)=2d0*tr*msq25_1(j,-j)*sub25_1(qg)
      msq(3,j,k)=xn*msq16_2(j,k)*sub16_2(qq)
      msq(4,j,k)=xn*(msq26_1(j,k)*sub26_1(gg)+msq26_1v(j,k)*sub26_1v)
      msq(5,j,k)=-(msq16_5(j,k)*sub16_5(qq)+msq16_5(j,k)*sub56_1(qq))/xn
      msq(6,j,k)=xn*(msq26_5(j,k)*sub26_5(gg)+msq26_5v(j,k)*sub26_5v
     .              +msq26_5(j,k)*sub56_2(qq))
      enddo
      
      j=0
      do k=1,nf
c--- g-q and g-qb cases
      msq(1,j,k)=2d0*tr*msq15_2(-k,k)*sub15_2(qg)
      msq(3,j,k)=xn*(msq16_2(j,k)*sub16_2(gg)+msq16_2v(j,k)*sub16_2v)
      msq(4,j,k)=xn*msq26_1(j,k)*sub26_1(qq)
      msq(5,j,k)=xn*(msq16_5(j,k)*sub16_5(gg)+msq16_5v(j,k)*sub16_5v
     .              +msq15_6(j,k)*sub56_1(qq))
      msq(6,j,k)=-(msq26_5(j,k)*sub26_5(qq)+msq26_5(j,k)*sub56_2(qq))/xn

      enddo
      do k=-nf,-1
c--- g-q and g-qb cases
      msq(1,j,k)=2d0*tr*msq15_2(-k,k)*sub15_2(qg)
      msq(3,j,k)=xn*(msq16_2(j,k)*sub16_2(gg)+msq16_2v(j,k)*sub16_2v)
      msq(4,j,k)=xn*msq26_1(j,k)*sub26_1(qq)
      msq(5,j,k)=xn*(msq16_5(j,k)*sub16_5(gg)+msq16_5v(j,k)*sub16_5v
     .              +msq15_6(j,k)*sub56_1(qq))
      msq(6,j,k)=-(msq26_5(j,k)*sub26_5(qq)+msq26_5(j,k)*sub56_2(qq))/xn

      enddo

      
      j=0
      k=0
c--- g-g case (real process is g(p1)+g(p2) --> qb(p5)+q(p6)
c---Hence 15 split multiplies q(15)+g(p2)-->Z+q(p6)
c---Hence 25 split multiplies g(p1)+q(p25)-->Z+q(p6)
      msq(1,j,k)=(msq15_2(+1,k)+msq15_2(+2,k)+msq15_2(+3,k)
     .           +msq15_2(+4,k)+msq15_2(+5,k))*sub15_2(qg)*2d0*tr
      msq(2,j,k)=(msq25_1(k,+1)+msq25_1(k,+2)+msq25_1(k,+3)
     .           +msq25_1(k,+4)+msq25_1(k,+5))*sub25_1(qg)*2d0*tr
      msq(3,j,k)=(msq16_2(-5,k)+msq16_2(-4,k)+msq16_2(-3,k)
     .           +msq16_2(-2,k)+msq16_2(-1,k))*sub16_2(qg)*2d0*tr
      msq(4,j,k)=(msq26_1(k,-5)+msq26_1(k,-4)+msq26_1(k,-3)
     .           +msq26_1(k,-2)+msq26_1(k,-1))*sub26_1(qg)*2d0*tr

      
cc--- Qflag subtraction pieces
c      do j=-nf,nf
c      do k=-nf,nf      
c
c      if (((j .gt. 0).and.(k .gt. 0)) .or. 
c     .    ((j .lt. 0).and.(k .lt. 0))) then
cc--q-q or qb-qb
c      if (j.eq.k) then
c      msq(1,j,k)=msq(1,j,k)+0.5d0*(xn-1d0/xn)
c     .  *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
c      msq(2,j,k)=msq(2,j,k)+0.5d0*(xn-1d0/xn)
c     .  *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
c      msq(3,j,k)=msq(3,j,k)+0.5d0*(xn-1d0/xn)
c     .  *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
c      msq(4,j,k)=msq(4,j,k)+0.5d0*(xn-1d0/xn)
c     .  *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
c      else
c      msq(1,j,k)=msq(1,j,k)+(xn-1d0/xn)
c     .  *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
c      msq(4,j,k)=msq(4,j,k)+(xn-1d0/xn)
c     .  *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
c      endif
c      elseif ((j .gt. 0).and.(k .lt. 0)) then
cc q-qbar
c      if (j.eq.-k) then
c      msq(1,j,k)=msq(1,j,k)+(xn-1d0/xn)
c     .  *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
c      msq(4,j,k)=msq(4,j,k)+(xn-1d0/xn)
c     .  *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
c      msq(6,j,k)=msq(6,j,k)+2d0*tr*dfloat(nf)
c     .  *(msq26_5(j,k)*sub56_2(gq)-msq56_2v(j,k)*sub56_2v)
c      else 
c      msq(1,j,k)=msq(1,j,k)+(xn-1d0/xn)
c     .  *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
c      msq(4,j,k)=msq(4,j,k)+(xn-1d0/xn)
c     .  *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
c      endif
cc--qbar-q
c      elseif ((j .lt. 0).and.(k .gt. 0)) then
c      if (j.eq.-k) then
c      msq(2,j,k)=msq(2,j,k)+(xn-1d0/xn)
c     .  *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
c      msq(3,j,k)=msq(3,j,k)+(xn-1d0/xn)
c     .  *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
c      msq(6,j,k)=msq(6,j,k)+2d0*tr*dfloat(nf)
c     . *(msq26_5(j,k)*sub56_2(gq)-msq56_2v(j,k)*sub56_2v)
c      else 
c      msq(2,j,k)=msq(2,j,k)+(xn-1d0/xn)
c     .  *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
c      msq(3,j,k)=msq(3,j,k)+(xn-1d0/xn)
c     .  *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
c
c      endif
c      endif
c
c
c      enddo
c      enddo

c--- Qflag subtraction pieces

c--q-q or qb-qb
      do j=1,nf
         do k=1,nf
            msq(1,j,k)=(xn-1d0/xn)
     .           *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
            msq(4,j,k)=(xn-1d0/xn)
     .           *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
         enddo
      enddo
      do j=-nf,-1
         do k=-nf,-1
            msq(1,j,k)=(xn-1d0/xn)
     .           *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
            msq(4,j,k)=(xn-1d0/xn)
     .           *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
         enddo
      enddo

      
      do j=1,nf
         k=j
         msq(1,j,k)=0.5d0*(xn-1d0/xn)
     .        *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
         msq(2,j,k)=0.5d0*(xn-1d0/xn)
     .        *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
         msq(3,j,k)=0.5d0*(xn-1d0/xn)
     .        *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
         msq(4,j,k)=0.5d0*(xn-1d0/xn)
     .        *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
      enddo
      do j=-nf,-1
         k=j
         msq(1,j,k)=0.5d0*(xn-1d0/xn)
     .        *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
         msq(2,j,k)=0.5d0*(xn-1d0/xn)
     .        *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
         msq(3,j,k)=0.5d0*(xn-1d0/xn)
     .        *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
         msq(4,j,k)=0.5d0*(xn-1d0/xn)
     .        *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
      enddo

c q-qbar
      do j=1,nf
         do k=-nf,-1
            msq(1,j,k)=msq(1,j,k)+(xn-1d0/xn)
     .           *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
            msq(4,j,k)=msq(4,j,k)+(xn-1d0/xn)
     .           *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
         enddo
      enddo
      do j=1,nf                 !subtract diagonal
         k=-j
         msq(1,j,k)=msq(1,j,k)-(xn-1d0/xn)
     .        *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
         msq(4,j,k)=msq(4,j,k)-(xn-1d0/xn)
     .        *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
      enddo
      do j=1,nf
         k=-j
         msq(1,j,k)=msq(1,j,k)+(xn-1d0/xn)
     .        *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
         msq(4,j,k)=msq(4,j,k)+(xn-1d0/xn)
     .        *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
         msq(6,j,k)=msq(6,j,k)+2d0*tr*dfloat(nf)
     .        *(msq26_5(j,k)*sub56_2(gq)-msq56_2v(j,k)*sub56_2v)
      enddo

c--qbar-q
      do j=-nf,-1
         do k=1,nf
            msq(2,j,k)=msq(2,j,k)+(xn-1d0/xn)
     .           *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
            msq(3,j,k)=msq(3,j,k)+(xn-1d0/xn)
     .           *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
         enddo
      enddo
      do j=-nf,-1               !subtract diagonal
         k=-j
         msq(2,j,k)=msq(2,j,k)-(xn-1d0/xn)
     .        *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
         msq(3,j,k)=msq(3,j,k)-(xn-1d0/xn)
     .        *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
      enddo
      do j=-nf,-1
         k=-j
         msq(2,j,k)=msq(2,j,k)+(xn-1d0/xn)
     .        *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
         msq(3,j,k)=msq(3,j,k)+(xn-1d0/xn)
     .        *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
         msq(6,j,k)=msq(6,j,k)+2d0*tr*dfloat(nf)
     .        *(msq26_5(j,k)*sub56_2(gq)-msq56_2v(j,k)*sub56_2v)
      enddo

      
      return
      end
      
