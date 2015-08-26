      double precision function realint(vector,wgt)
      implicit none
      include 'constants.f'
      include 'realonly.f'
      include 'virtonly.f'
      include 'noglue.f'
      include 'vegas_common.f'
      include 'ptilde.f'
      include 'npart.f'
      include 'scale.f'
      include 'facscale.f'
      include 'efficiency.f'
      include 'maxwt.f'
      include 'process.f'
      include 'dynamicscale.f'
      include 'dipolescale.f'
      integer ih1,ih2,j,k,nd,nmax,nmin,nvec,doFill
      double precision vector(mxdim),W,val,xint
      double precision sqrts,fx1(-nf:nf),fx2(-nf:nf)
      double precision p(mxpart,4),pjet(mxpart,4),p1ext(4),p2ext(4)
      double precision pswt,rscalestart,fscalestart
      double precision s(mxpart,mxpart),wgt,msq(-nf:nf,-nf:nf)
      double precision msqc(maxd,-nf:nf,-nf:nf),xmsq(0:maxd),xmsqjk
      double precision flux,BrnRat,xreal,xreal2
      double precision xx1,xx2,q(mxpart,4),dot,q2,qt2,xqtcut
      integer n2,n3
      double precision mass2,width2,mass3,width3
      common/breit/n2,n3,mass2,width2,mass3,width3
      common/xreal/xreal,xreal2
      logical bin,first,failed
      logical incldip(0:maxd),includedipole,includereal
      external qqb_z2jet,qqb_z1jet_gs,qqb_w2jet,qqb_w1jet_gs
      common/density/ih1,ih2
      common/energy/sqrts
      common/bin/bin
      common/Pext/p1ext,p2ext
      common/nmax/nmax
      common/BrnRat/BrnRat
      common/nmin/nmin
      common/incldip/incldip
      integer nproc
      common/nproc/nproc
      common/qtcut/xqtcut

      common/doFill/doFill
      external hists_fill

      logical binner
      external binner

      integer ii,jj,kk

      data p/48*0d0/
      data first/.true./
      save first,rscalestart,fscalestart


      if (first) then
         first=.false.
         rscalestart=scale
         fscalestart=facscale
      endif
      ntotshot=ntotshot+1
      pswt=0d0
      realint=0d0      

      W=sqrts**2
      
      if (first) then
         write(6,*)
         write(6,*) 'nmin=',nmin,',nmax=',nmax
         write(6,*)
         first=.false.
      endif
      
      npart=4
      call gen4(vector,p,pswt,*999)

c      print*,'phase space in real'
c      print*,p(3,1),p(3,2),p(3,3),p(3,4)
c      print*,p(4,1),p(4,2),p(4,3),p(4,4)
c      print*,'mass',sqrt((p(4,4)+p(3,4))**2
c     +     - (p(4,1)+p(3,1))**2
c     +     - (p(4,2)+p(3,2))**2
c     +     - (p(4,3)+p(3,3))**2)
c      print*,'y',0.5d0*log((p(4,4)+p(3,4) + (p(4,3)+p(3,3)))/
c     +     (p(4,4)+p(3,4) - (p(4,3)+p(3,3))))
c      print*,'pt',sqrt((p(4,1)+p(3,1))**2 + (p(4,2)+p(3,2))**2)
c      print*
      
      nvec=npart+2

      q2=2*dot(p,3,4)
      qt2=(p(3,1)+p(4,1))**2+(p(3,2)+p(4,2))**2

      call dotem(nvec,p,s)
      
c---impose cuts on final state
      call masscuts(s,*999)


c----reject event if any s(i,j) is too small
      call smalls(s,npart,*999)
      

c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead set to zero
      includereal=includedipole(0,p)
      incldip(0)=includereal 
      
CC   Dynamic scale: set it only if point passes cuts

      if(dynamicscale.and.includereal) then
       call scaleset(q2)
       dipscale(0)=facscale
      endif


      if (includereal .eqv. .false.) then
        do j=-nf,nf
        do k=-nf,nf
          msq(j,k)=0d0
        enddo
        enddo
      endif
      
C---- binner cut
      if (binner(p(3,:),p(4,:)).eqv..false.) goto 999
C---- min qt cut
      if(dsqrt(qt2/q2).lt.xqtcut) goto 999
      
c---- generate collinear points that satisy the jet cuts (for checking)
c      call singgen(p,s,*998)
            
c----calculate the x's for the incoming partons from generated momenta

      xx1=two*(p(1,4)*p2ext(4)-p(1,3)*p2ext(3))/W
      xx2=two*(p(2,4)*p1ext(4)-p(2,3)*p1ext(3))/W

      
      if ((xx1 .gt. 1d0) .or. (xx2 .gt. 1d0)) then
         realint=0d0
         return
      endif



c--- Calculate the required matrix elements    

      if(nproc.eq.3) then
       if (includereal) call qqb_z2jet(p,msq)
       call qqb_z1jet_gs(p,msqc)
      else
       if (includereal) call qqb_w2jet(p,msq)
       call qqb_w1jet_gs(p,msqc)
      endif 

      do nd=0,ndmax
      xmsq(nd)=0d0
      enddo
      
            
      flux=fbGeV2/(two*xx1*xx2*W)


  777 continue    
      do nd=0,ndmax
      xmsq(nd)=0d0
      enddo
           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      do nd=0,ndmax

      call fdist(ih1,xx1,dipscale(nd),fx1)
      call fdist(ih2,xx2,dipscale(nd),fx2)
      
      do j=-nf,nf
      do k=-nf,nf

      if (ggonly) then
      if ((j.ne.0) .or. (k.ne.0)) goto 20
      endif

      if (gqonly) then
      if (((j.eq.0).and.(k.eq.0)) .or. ((j.ne.0).and.(k.ne.0))) goto 20
      endif      
      
      if (noglue) then 
      if ((j.eq.0) .or. (k.eq.0)) goto 20
      endif

      if (realonly) then 
        if(nd.eq.0) then
         xmsq(0)=xmsq(0)+fx1(j)*fx2(k)*msq(j,k)
        else
         xmsq(nd)=0d0
        endif
      elseif (virtonly) then
        if(nd.eq.0) then
         xmsq(0)=0d0
        else
         xmsq(nd)=xmsq(nd)+fx1(j)*fx2(k)*(-msqc(nd,j,k))
        endif
      else

        if(nd.eq.0) then
         xmsqjk=fx1(j)*fx2(k)*msq(j,k)
        else
         xmsqjk=fx1(j)*fx2(k)*(-msqc(nd,j,k))
        endif

        xmsq(nd)=xmsq(nd)+xmsqjk         

      endif

 20   continue
      enddo
      enddo

      enddo

      realint=0d0
      xint=0d0

c---trial with weight of real alone
c---first set up all dipole contributions
c---this is the value of integral including subtractions
      do nd=0,ndmax
        xmsq(nd)=xmsq(nd)*flux*pswt/BrnRat
        failed=.false.
        
c--- if this dipole has no contribution, go to end of loop
c        if (xmsq(nd) .eq. 0d0) goto 997         
         
        if (nd .eq. 0) then
c---if there's no real contribution, record the event as failing to pass cuts
          if (xmsq(nd) .eq. 0d0) then
             failed=.true.
             goto 996
          endif
        else
c--- if this dipole has no contribution, go to end of loop
          if (xmsq(nd) .eq. 0d0) goto 997         
c---check whether each counter-event passes the cuts
          do j=1,mxpart
          do k=1,4
          q(j,k)=ptilde(nd,j,k)
          enddo
          enddo
          incldip(nd)=includedipole(nd,q)
          if (incldip(nd) .eqv. .false.) failed=.true.
        endif

 996    if (failed) then
          if (nd .eq. 0) then
            ncutzero=ncutzero+1
            ntotzero=ntotzero+1
          endif
          call dotem(nvec,p,s)
          xmsq(nd)=0d0
          goto 997         
        endif
c---if it does, add to total
        xint=xint+xmsq(nd)

        val=xmsq(nd)*wgt
                
c--- update the maximum weight so far, if necessary
        if (dabs(val) .gt. wtmax) then
          wtmax=dabs(val)
        endif

c---if we're binning, add to histo too
        if (bin) then
          call getptildejet(nd,pjet)
          call dotem(nvec,pjet,s)
          val=val/dfloat(itmx)
          if (nd .eq. 0) then
            call plotter(pjet,val,3+nd)
          else
              npart=npart-1
           call plotter(pjet,val,3+nd)
              npart=npart+1
         endif
        endif
c---otherwise, skip contribution
 997    continue
      enddo

      call dotem(nvec,p,s)

c 998  continue


C     Fill only if it's last iteration
      if (doFill.ne.0) then
          call hists_fill(p(3,:),p(4,:),xint*wgt)
      endif
      realint=xint

      xreal=xreal+xint*wgt/dfloat(itmx)
      xreal2=xreal2+(xint*wgt)**2/dfloat(itmx)
      

      return

 999  realint=0d0

      ntotzero=ntotzero+1
 
      return
      end
















