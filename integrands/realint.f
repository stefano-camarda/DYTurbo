      double precision function realint(vector,wgt,f)
      implicit none
      include 'constants.f'
      include 'noglue.f'
      include 'vegas_common.f'
      include 'ptilde.f'
      include 'npart.f'
      include 'scale.f'
      include 'facscale.f'
      include 'process.f'
      include 'dynamicscale.f'
      include 'dipolescale.f'
      include 'qcdcouple.f'
      include 'options.f'
      integer ih1,ih2,j,k,nd,nvec
      double precision f(*)
      double precision vector(mxdim),W,val
      double precision sqrts,fx1(-nf:nf),fx2(-nf:nf)
      double precision dipfx1(0:maxd,-nf:nf),dipfx2(0:maxd,-nf:nf)
      double precision p(mxpart,4),pjet(mxpart,4),p1ext(4),p2ext(4)
      double precision ptrans(mxpart,4)
      double precision pswt
      double precision s(mxpart,mxpart),wgt,msq(-nf:nf,-nf:nf)
      double precision msqc(maxd,-nf:nf,-nf:nf),xmsq(0:maxd)
      double precision flux,BrnRat
      double precision xx1,xx2,dot,q2
      logical incldip(0:maxd),includedipole,includereal
      external qqb_z2jet,qqb_z1jet_gs,qqb_w2jet,qqb_w1jet_gs
      common/density/ih1,ih2
      common/energy/sqrts
      common/Pext/p1ext,p2ext
      common/BrnRat/BrnRat
      common/incldip/incldip
      integer nproc
      common/nproc/nproc

C      external hists_fill
      external hists_setpdf
      external hists_real_dipole
      external hists_real_event

      double precision x,omx,sij,sik,sjk
      integer ip,jp,kp

      data p/48*0d0/

      integer npdf,maxpdf
      double precision gsqcentral

      pswt=0d0
      realint=0d0      
      do npdf=0,totpdf-1
         f(npdf+1)=0d0
      enddo

      W=sqrts**2
      
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
c      qt2=(p(3,1)+p(4,1))**2+(p(3,2)+p(4,2))**2

      call dotem(nvec,p,s)
      
c---impose cuts on final state
      call masscuts(s,*999)

c----reject event if any s(i,j) is too small
      call smalls(s,npart,*999)
      
c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead set to zero
      includereal=includedipole(0,p)
      incldip(0)=includereal 

c     check which dipoles are to be included
      incldip(5)=incldip(0)
      incldip(6)=incldip(0) 
      do nd=1,6
         if (nd.eq.1) then
            ip = 1
            jp = 5
            kp = 2
         endif
         if (nd.eq.2) then
            ip = 2
            jp = 5
            kp = 1
         endif
         if (nd.eq.3) then
            ip = 1
            jp = 6
            kp = 2
         endif
         if (nd.eq.4) then
            ip = 2
            jp = 6
            kp = 1
         endif
         if (nd.eq.5) then
            ip = 1
            jp = 5
            kp = 6
         endif
         if (nd.eq.6) then
            ip = 2
            jp = 6
            kp = 5
         endif
         sij=two*dot(p,ip,jp)
         sik=two*dot(p,ip,kp)
         sjk=two*dot(p,jp,kp)

         if (nd.le.4) then
            omx=-(sij+sjk)/sik
         else
            omx=-sjk/(sij+sik)
         endif
         x=one-omx
                     
c     ptrans is the transformation of p into the dipole kinematic
         call transform(p,ptrans,x,ip,jp,kp)
c     ptrans is stored into ptilde(nd), to be used by dips in  dipolesub.f
         call storeptilde(nd,ptrans)
         if (nd.le.4) then
            incldip(nd)=includedipole(nd,ptrans)
         endif
      enddo
      
c     If the real and all the dipoles fail cuts, then exit
      if (.not.incldip(0).and.
     +     .not.incldip(1).and.
     +     .not.incldip(2).and.
     +     .not.incldip(3).and.
     +     .not.incldip(4).and.
     +     .not.incldip(5).and.
     +     .not.incldip(6)) then
         return
      endif
      
CC   Dynamic scale: set it only if point passes cuts 
c     !!!this does not look correct, should be only if(dynamicscale)!!!
      if(dynamicscale.and.includereal) then
       call scaleset(q2)
       dipscale(0)=facscale
      endif

c----calculate the x's for the incoming partons from generated momenta
      xx1=two*(p(1,4)*p2ext(4)-p(1,3)*p2ext(3))/W
      xx2=two*(p(2,4)*p1ext(4)-p(2,3)*p1ext(3))/W

      if ((xx1 .gt. 1d0) .or. (xx2 .gt. 1d0)) then
         return
      endif

c--- Calculate the required matrix elements  (dipscale(nd) are set appropriately in dipolesub.f)
      if (pdferr) then
         call setpdf(0)
      endif
      gsqcentral=gsq
      if(nproc.eq.3) then
       if (includereal) call qqb_z2jet(p,msq)
       call qqb_z1jet_gs(p,msqc)
      else
       if (includereal) call qqb_w2jet(p,msq)
       call qqb_w1jet_gs(p,msqc)
      endif 

      flux=fbGeV2/(two*xx1*xx2*W)

c     skip PDF loop in the preconditioning phase
      maxpdf=0
      if (dofill.ne.0) maxpdf = totpdf-1
      
c     start PDF loop
      do npdf=0,maxpdf
         call setpdf(npdf)
         call hists_setpdf(npdf)
c     intitialise xmsq to 0 for the real and all dipoles
         do nd=0,ndmax
            xmsq(nd)=0d0
         enddo

c     evaluate PDFs
         if (dynamicscale) then
            do nd=0,6
               if (incldip(nd)) then
                  call fdist(ih1,xx1,dipscale(nd),dipfx1(nd,:))
                  call fdist(ih2,xx2,dipscale(nd),dipfx2(nd,:))
               endif
            enddo
         else
            call fdist(ih1,xx1,dipscale(0),fx1)
            call fdist(ih2,xx2,dipscale(0),fx2)
            do nd=0,6
               do j=-nf,nf
                  dipfx1(nd,j)=fx1(j)
                  dipfx2(nd,j)=fx2(j)
               enddo
            enddo
         endif

c     calculate xmsq for the real event
         if (includereal) then
            do j=-nf,nf
               do k=-nf,nf
                  xmsq(0)=xmsq(0)+dipfx1(0,j)*dipfx2(0,k)*msq(j,k)
     .                 *(gsq/gsqcentral)**2
               enddo
            enddo
         endif
      
c     calculate xmsq for the dipole contributions
         do nd=1,6
            if (incldip(nd)) then
               do j=-nf,nf
                  do k=-nf,nf
                     xmsq(nd)=xmsq(nd)
     .                    +dipfx1(nd,j)*dipfx2(nd,k)*(-msqc(nd,j,k))
     .                    *(gsq/gsqcentral)**2
                  enddo
               enddo
            endif
         enddo

c     Sum up the real and all the dipole contributions
         xmsq(0)=xmsq(0)+xmsq(5)+xmsq(6)

c---trial with weight of real alone
c---first set up all dipole contributions
c---this is the value of integral including subtractions
c     
         do nd=0,4              !Start loop on real+dipoles contributions
c---  if this dipole has no contribution, go to end of loop
            if (xmsq(nd).eq.0d0) cycle
            xmsq(nd)=xmsq(nd)*flux*pswt/BrnRat
         
c---  add to total
            f(npdf+1)=f(npdf+1)+xmsq(nd)

C---  Fill only if it's last iteration
            if (doFill.ne.0) then
c     this call to getptildejet can be removed, because pjet and ptildejet are both equal to ptrans,
c     which is the transformation of p into the dipole kinematic
c     ptrans is stored in ptilde(nd), which should be used for filling
               call getptildejet(nd,pjet)
               val=xmsq(nd)*wgt
C---        store information per each dipole
              call hists_real_dipole(pjet(3,:),pjet(4,:),val,nd)
C             call hists_fill(ptilde(nd,3,:),ptilde(nd,4,:),val)
            endif
         enddo                  !End loop on real+dipoles contributions

C---  Fill only if it's last iteration
         if (doFill.ne.0) then
C---     fill the dipole contribution to each bin separatelly
            call hists_real_event()
         endif
      enddo ! end of PDF loop
      realint = f(1)
      return

 999  realint=0d0

      return
      end
