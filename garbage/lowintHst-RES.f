       double precision function lowintHst(r,wgt)
       implicit none
       include 'npart.f'
       include 'constants.f'
       include 'masses.f'
       include 'vegas_common.f'
       include 'limits.f'
       include 'zerowidth.f'
       include 'dynamicscale.f'
!       include 'jetlabel.f'
       double precision r(mxdim),sqrts,resumm,tempp,tempp2,BrnRat,ran2
       integer ii,j,ij,azloopmax,azloop,jstop
       real *8 jac,m,m2,qt,qt2,s,cosh2y,cutoff,val,wgt
       real *8 p(mxpart,4),pV(4),p4cm(4),ptemp(4),kap1(4),kap1b(4)
       real *8 y,ymax,phi,xx1,xx2,xxtemp,cos_th,phi_lep,dexpy,dexpmy
       logical bin
       double precision qtmax,qtmin,lowintHst0
!       common/qtcut/xqtcut
       common/bin/bin
       integer n2,n3
       double precision mass2,width2,mass3,width3,x1,msq,wt,switch,kk
     & ,zeta1,zeta1b,kt1,kt2,mt2,costh_CS,qP1,qP2,ritmx
       common/breit/n2,n3,mass2,width2,mass3,width3
       common/energy/sqrts
       common/BrnRat/BrnRat
      common/ritmx/ritmx
       external resumm
       logical cuts

      integer mysw
      common/mysw/mysw
!      logical bin2
!      common/bin2/bin2

       npart=2

!       qtcutmaxNEW=qtcutmax
       lowintHst=0d0
       lowintHst0=0d0

       jstop=0

       azloopmax=500 !10000
       azloop=0

C      The number of jets is zero for this piece
!       jets=0

       ptemp(1)=0d0
       ptemp(2)=0d0
       ptemp(3)=0d0
       ptemp(4)=0d0
       do ii=1,mxpart
         p(ii,1)=0d0
         p(ii,2)=0d0
         p(ii,3)=0d0
         p(ii,4)=0d0
       enddo

! Generate the vector boson (dilepton) invariant mass

        jac=1d0
        x1=r(1)
        call  breitw(x1,wsqmin,wsqmax,mass3,width3,msq,wt) 
        m2=msq
        m=dsqrt(m2)
        jac=jac*wt


C     Dynamic scale
      if(dynamicscale) call scaleset(m2)


       ymax=0.5d0*dlog(sqrts**2/m2)
       y=-ymax+2d0*ymax*r(2)
       jac=jac*2d0*ymax  

       dexpy=dexp(y)
       dexpmy=dexp(-y)
       cosh2y=((dexpy+dexpmy)*0.5d0)**2

!   qT kinematical limit
       qtmax=dsqrt( (sqrts**2+m2)**2/(4d0*sqrts**2*cosh2y) - m2 )

!       if(qtmax.lt.qtcutmax) qtcutmaxNEW=qtmax
!       qtcutmaxNEW=qtmax

       qtmin=0.1d0

       qt=qtmin+qtmax*r(3)
       jac=jac*(qtmax)
       qt2=qt**2

       xx1=dsqrt(m2/sqrts**2)*dexpy
       xx2=dsqrt(m2/sqrts**2)*dexpmy


! incoming quarks
       p(1,1)=0d0
       p(1,2)=0d0
       p(1,3)=xx1*0.5d0*sqrts
       p(1,4)=xx1*0.5d0*sqrts
       p(2,1)=0d0
       p(2,2)=0d0
       p(2,3)=-xx2*0.5d0*sqrts
       p(2,4)=xx2*0.5d0*sqrts

! First lepton direction: Cos of the polar angle
       cos_th=-1d0+2d0*r(4)
       jac=jac*2d0 

       mt2=m2+qt2

! LOOP over (vector boson and lepton) azimuthal angles 
       do j=1,azloopmax,1

! Vector boson azimuthal angle
        phi=twopi*ran2()
!          phi=twopi*r(5)
!       write(*,*) "phi",phi
! Lepton decay in the center of mass frame 
! First lepton direction:  azimuthal angle
       phi_lep=twopi*ran2()
!       write(*,*) "phi_lep",phi_lep
!       write(*,*) "diff",phi-phi_lep

!  vector boson momentum: pV(4)^2-pV(1)^2-pV(2)^2-pV(3)^2=m2
       pV(1)=qt*dcos(phi)
       pV(2)=qt*dsin(phi)
       pV(3)=0.5d0*dsqrt(mt2)*(dexp(y)-dexp(-y))
       pV(4)=0.5d0*dsqrt(mt2)*(dexp(y)+dexp(-y))
!       pV(3)=0.5d0*dexpmy*(dexp(2d0*y)-1d0)*dsqrt(m2+qt*qt)
!       pV(4)=pV(3)*(dexpy+dexpmy)/(dexpy-dexpmy)


! momentum of the first lepton 

       p4cm(4)=m/2d0
       p4cm(1)=p4cm(4)*dsin(dacos(cos_th))*dsin(phi_lep)
       p4cm(2)=p4cm(4)*dsin(dacos(cos_th))*dcos(phi_lep)
       p4cm(3)=p4cm(4)*cos_th

! Boost to go in the Lab frame
       call boost(m,pV,p4cm,ptemp)
       p(4,1)=ptemp(1)
       p(4,2)=ptemp(2)
       p(4,3)=ptemp(3)
       p(4,4)=ptemp(4)


!  momentum of the second lepton
       p(3,4)=pV(4)-p(4,4)
       p(3,1)=pV(1)-p(4,1)
       p(3,2)=pV(2)-p(4,2)
       p(3,3)=pV(3)-p(4,3)


!!!!! CHECK STEFANO'S FORMULA


!       kt1=(1d0+pV(3)/(dsqrt(m2)+pV(4)))*pV(1)/2d0   !
!       kt2=(1d0+pV(3)/(dsqrt(m2)+pV(4)))*pV(2)/2d0   !  MY PRESCRIPTION

       kt1=pV(1)/2d0   !
       kt2=pV(2)/2d0   !  CS FRAME PRESCRIPTION

!       mt2=m2+qt2
       zeta1=1d0/m2/2d0*(m2+2d0*(pV(1)*kt1+pV(2)*kt2)
     &+dsqrt((m2+2d0*(pV(1)*kt1+pV(2)*kt2))**2-4d0*mt2*(kt1**2+kt2**2)))

       qP1=(pV(4)-pV(3))*sqrts/2d0
       qP2=(pV(4)+pV(3))*sqrts/2d0

       kap1(4)=sqrts/2d0*(zeta1*m2/2d0/qP1
     & +(kt1**2+kt2**2)/zeta1*qP1/m2/sqrts**2*2d0)
       kap1(1)=kt1
       kap1(2)=kt2
       kap1(3)=sqrts/2d0*(zeta1*m2/2d0/qP1
     & -(kt1**2+kt2**2)/zeta1*qP1/m2/sqrts**2*2d0)

!       zeta1b=(2d0-qt2*sqrts/(2d0*(m+pV(4))*qP2))*mt2/2d0/m2
!       kap1b(4)=(1d0-qt2*sqrts/(4d0*(m+pV(4))*qP2))*qP2/sqrts
!     & +qt2/(4d0*(m+pV(4)))
!       kap1b(1)=kt1
!       kap1b(2)=kt2
!       kap1b(3)=(1d0-qt2*sqrts/(4d0*(m+pV(4))*qP2))*qP2/sqrts
!     & -qt2/(4d0*(m+pV(4)))

       costh_CS=1d0-4d0*(kap1(4)*p(4,4)-kap1(3)*p(4,3)-kap1(2)*p(4,2)
     & -kap1(1)*p(4,1))/m2

!       write(*,*) "CHECK STEFANO", cos_th,costh_CS,cos_th/costh_CS
!     &,"k",kap1b(1)/kap1(1),kap1b(2)/kap1(2),kap1b(3)/kap1(3)
!     &,kap1b(4)/kap1(4)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead bail out

      if(cuts(p,0) .eqv. .true.) goto 999


!  SWITCHING FUNCTIONS
      switch=1d0
      if(qt.ge.m*3/4d0)  switch=dexp(-(m*3/4d0-qt)**2/(m/2d0)**2)           ! GAUSS SWITCH
!!!!!!
!      if(qt.ge.m)  switch=dexp(-(m-qt)**2/(m/2d0)**2)           ! GAUSS SWITCH
!      if(qt.ge.m/2d0)  switch=dexp(-(m/2d0-qt)**2/(m/2d0)**2)              ! GAUSS SWITCH
!      if(qt.ge.m/2d0)        switch=dexp(-2d0*(m/2d0-qt)**2/(m/2d0)**2)    ! GAUSS SWITCH FASTER
!      if(qt.ge.m/2d0)        switch=dexp((m2/4d0-qt2)/m2)                  ! EXP SWITCH
!      if(qt.ge.m/2d0)        switch=dexp((m2/4d0-qt2)/m2*4d0)              ! EXP SWITCH FASTER
!      if(qt.ge.m/2d0)        switch=(dcos(pi/50d0*(qt-45.6d0))+1d0)/2d0    ! COS SWITCH
!      if(qt.ge.95.6)         switch=0d0                                    ! COS SWITCH

      if(switch.le.0.01d0)    goto 999    ! 


      if(jstop.ne.1) then
       jstop=1

!--Call to the resummation part
       if((qt.lt.qtmin).or.(switch.lt.0.01)) then
        tempp=0d0
       else
!        tempp=resumm(cos_th,m,qt,y)/(8d0/3d0)   ! MY PRESCRIPTION
        tempp=resumm(costh_CS,m,qt,y)/(8d0/3d0)  ! CS PRESCRIPTION
!        tempp=dexp(-qt**2)
        if(tempp.ne.tempp) tempp=0d0    !  TO AVOID NAN 
       endif

       lowintHst0=jac*tempp/BrnRat
       lowintHst0=lowintHst0*switch    ! SWITCHING

      endif

      azloop=azloop+1      
      val=lowintHst0*wgt/dfloat(azloopmax)
      
       if (bin) then
          val=val/dfloat(itmx)*ritmx   ! This factor needed becouse we have changed the number of itmx in integrate-RES.f
!          write(*,*) "itmx",itmx
          call plotter(p,val,0)
          mysw=mysw+1
       endif

 999  continue


      enddo !end azloop


      lowintHst=lowintHst0*dfloat(azloop)/dfloat(azloopmax)


       RETURN
       END
