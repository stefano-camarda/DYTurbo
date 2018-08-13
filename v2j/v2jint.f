      double precision function v2jint(vector,wgt,f)
      implicit none
      include 'constants.f'
      include 'noglue.f'
      include 'mxdim.f'
      include 'npart.f'
      include 'scale.f'
      include 'facscale.f'
      include 'process.f'
      include 'dynamicscale.f'
      include 'qcdcouple.f'
      include 'options.f'
      integer ih1,ih2,i,j,k,nd,nvec
      double precision f(*)
      double precision vector(mxdim),W,val
      double precision sqrts,fx1(-nf:nf),fx2(-nf:nf)
      double precision p(mxpart,4),p1ext(4),p2ext(4)
      double precision pswt
      double precision s(mxpart,mxpart),wgt,msq(-nf:nf,-nf:nf)
      double precision xmsq
      double precision flux,BrnRat
      double precision xx1,xx2,dot,q2,qt2
      logical includereal
      external qqb_z2jet,qqb_w2jet
      common/density/ih1,ih2
      common/energy/sqrts
      common/Pext/p1ext,p2ext
      common/BrnRat/BrnRat
      integer nproc
      common/nproc/nproc

      external hists_fill
      external hists_setpdf

      data p/pdim*0d0/

      integer npdf,maxpdf
      double precision gsqcentral

      logical cuts,failedcuts,makecuts
      common/makecuts/makecuts
      logical binner
      external binner
      
      pswt=0d0
      v2jint=0d0      
      do npdf=0,totpdf-1
         f(npdf+1)=0d0
      enddo

      W=sqrts**2
      
      npart=4
      call gen4(vector,p,pswt,*999)

      nvec=npart+2

c     Compute qt2,q2
      q2=2*dot(p,3,4)
      qt2=(p(3,1)+p(4,1))**2+(p(3,2)+p(4,2))**2

      call dotem(nvec,p,s)
      
c---impose cuts on final state
      call masscuts(s,*999)

c----reject event if any s(i,j) is too small
      call smalls(s,npart,*999)
      
c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead set to zero
      includereal=.true.

      if(sqrt(qt2/Q2).lt.xqtcut.or.sqrt(qt2).lt.qtcut) then
         includereal=.false.
         return
      endif

c     apply cuts on jets
      if(sqrt(p(5,1)**2+p(5,2)**2).lt.20) then
         includereal=.false.
         return
      endif

      if(sqrt(p(6,1)**2+p(6,2)**2).lt.20) then
         includereal=.false.
         return
      endif

c--- check the lepton cuts
      failedcuts=(cuts(p,0))
      if (failedcuts) includereal=.false.
      if (.not.binner(p(3,:),p(4,:))) includereal=.false.
      
c     If cuts fail, then exit
      if (.not.includereal) then
         return
      endif
      
CC   Dynamic scale: set it only if point passes cuts 
      if (dynamicscale) then
         call scaleset(q2)
      endif

c----calculate the x's for the incoming partons from generated momenta
      xx1=two*(p(1,4)*p2ext(4)-p(1,3)*p2ext(3))/W
      xx2=two*(p(2,4)*p1ext(4)-p(2,3)*p1ext(3))/W

      if ((xx1 .gt. 1d0) .or. (xx2 .gt. 1d0)) then
         return
      endif

c--- Calculate the required matrix elements
      if (pdferr) then
         call dysetpdf(0)
      endif
      gsqcentral=gsq
      if(nproc.eq.3) then
         call qqb_z2jet(p,msq)
      else
         call qqb_w2jet(p,msq)
      endif 

      flux=fbGeV2/(two*xx1*xx2*W)

c     skip PDF loop in the preconditioning phase
      maxpdf=0
      if (dofill.ne.0) maxpdf = totpdf-1
      
c     start PDF loop
      do npdf=0,maxpdf
         call dysetpdf(npdf)
         call hists_setpdf(npdf)
         
c     intitialise xmsq to 0
         xmsq=0d0

c     evaluate PDFs
         call fdist(ih1,xx1,facscale,fx1)
         call fdist(ih2,xx2,facscale,fx2)

c     calculate xmsq for the real event
         do j=-nf,nf
            do k=-nf,nf
               xmsq=xmsq+fx1(j)*fx2(k)*msq(j,k)
     .              *(gsq/gsqcentral)**2
            enddo
         enddo
      
         xmsq=xmsq*flux*pswt/BrnRat
         
c     save pdf results
         f(npdf+1)=xmsq

c     Fill only in the last iteration
         if (doFill.ne.0) then
            val=xmsq*wgt
            call hists_fill(p(3,:),p(4,:),val)
         endif

      enddo                     ! end of PDF loop
      v2jint = f(1)
      return

 999  v2jint=0d0

      return
      end
