CC    Used to compute Higgs or W(Z) cross section at NLO only

      double precision function lowint(r,wgt,f)
      implicit none
      include 'constants.f'
      include 'npart.f'
      include 'mxdim.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'facscale.f'
      include 'noglue.f'
      include 'dynamicscale.f'
      include 'options.f'
      include 'qcdcouple.f'

      double precision f(*)
      double precision qq2
      double precision dot
      external dot

      integer ih1,ih2,j,k,nvec
      double precision r(mxdim),W,sqrts,xmsq,val,
     . fx1(-nf:nf),fx2(-nf:nf),p(mxpart,4),pswt
      double precision wgt,msq(-nf:nf,-nf:nf)
      double precision xx(2),flux,BrnRat
      logical includedipole
      integer nproc
      common/density/ih1,ih2
      common/energy/sqrts
      common/x1x2/xx
      common/BrnRat/BrnRat
      common/nproc/nproc
      data p/48*0d0/

      external hists_fill

      integer npdf,maxpdf
      double precision gsqcentral
      
      lowint=0d0
      do npdf=0,totpdf-1
         f(npdf+1)=0d0
      enddo

      W=sqrts**2
  
      npart=3     
      call gen3(r,p,pswt,*999)
c      print*,p(1,1),p(1,2),p(1,3),p(1,4)
c      print*,p(2,1),p(2,2),p(2,3),p(2,4)
c      print*,p(3,1),p(3,2),p(3,3),p(3,4)
c      print*,p(4,1),p(4,2),p(4,3),p(4,4)
c      print*,'mass',sqrt((p(4,4)+p(3,4))**2
c     +     - (p(4,1)+p(3,1))**2
c     +     - (p(4,2)+p(3,2))**2
c     +     - (p(4,3)+p(3,3))**2)
c      print*,'y',0.5d0*log((p(4,4)+p(3,4) + (p(4,3)+p(3,3)))/
c     +     (p(4,4)+p(3,4) - (p(4,3)+p(3,3))))
c      print*,'pt',sqrt((p(4,1)+p(3,1))**2 + (p(4,2)+p(3,2))**2)

C     Dynamic scale
      if(dynamicscale) then
         qq2=2*dot(p,3,4)
         call scaleset(qq2)
      endif

      nvec=npart+2
      call dotem(nvec,p,s)
      
      call masscuts(s,*999)
      
c----reject event if any s(i,j) is too small
      call smalls(s,npart,*999)                                                 

c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead bail out
      if (includedipole(0,p) .eqv. .false.) then
         return
      endif
      
      xx(1)=-2d0*p(1,4)/sqrts
      xx(2)=-2d0*p(2,4)/sqrts

c--- Calculate the required matrix elements      
      if (pdferr) then
         call setpdf(0)
      endif
      gsqcentral=gsq

       if(nproc.eq.3) then
          call qqb_z_g(p,msq)
       else
          call qqb_w_g(p,msq)
       endif      

      flux=fbGeV2/(2d0*xx(1)*xx(2)*W)

c     skip PDF loop in the preconditioning phase
      maxpdf=0
      if (doFill.ne.0) maxpdf = totpdf-1
      
c     start PDF loop
      do npdf=0,maxpdf
         call setpdf(npdf)
         call hists_setpdf(npdf)
c     intitialise xmsq to 0
         xmsq=0d0
      
c--- calculate PDF's  
         call fdist(ih1,xx(1),facscale,fx1)
         call fdist(ih2,xx(2),facscale,fx2)

         do j=-nf,nf
         do k=-nf,nf
               
            if (msq(j,k).eq.0d0) cycle
         if (ggonly) then
            if ((j.ne.0) .or. (k.ne.0)) cycle
         endif

         if (gqonly) then
       if (((j.eq.0).and.(k.eq.0)) .or. ((j.ne.0).and.(k.ne.0))) cycle
         endif
         
         if (noglue) then 
            if ((j.eq.0) .or. (k.eq.0)) cycle
         endif
         
c     gsq/gsqcentral correct for a possibly different value of alphas in the PDF (at O(alphas))
         xmsq=xmsq+fx1(j)*fx2(k)*msq(j,k)*(gsq/gsqcentral)
         
         enddo
         enddo

         xmsq=flux*pswt*xmsq/BrnRat
         f(npdf+1)=xmsq

         if (doFill.ne.0) then
            val=xmsq*wgt
            call hists_fill(p(3,:),p(4,:),val)
         endif
      enddo

      lowint = f(1)
      
      return

 999  continue
      
      return
      end
