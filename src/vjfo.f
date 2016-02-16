      function vjfo(m,pt,y)
      implicit none
c     input
      double precision m, pt, y
      double precision q2

      double precision ud,us,ub,cd,cs,cb
      common/cabib/ud,us,ub,
     &             cd,cs,cb

      double precision gGf,gw,xw,gwsq,esq,vevsq
      common/ewcouple/gGf,gw,xw,gwsq,esq,vevsq

      double precision 
     & md,mu,ms,mc,mb,mt,
     & mel,mmu,mtau,
     & hmass,hwidth,
     & wmass,wwidth,
     & zmass,zwidth,
     & twidth,
     & tauwidth,
     & mtausq,mcsq,mbsq
      common/dymasses/
     & md,mu,ms,mc,mb,mt,
     & mel,mmu,mtau,
     & hmass,hwidth,
     & wmass,wwidth,
     & zmass,zwidth,
     & twidth,
     & tauwidth,
     & mtausq,mcsq,mbsq
      
      double precision brz,brw
      
      double precision sroot
      common/energy/sroot
      double precision aemmz
      common/em/aemmz

      double precision gevpb
      common/gevpb/gevpb

      integer iih1,iih2
      common/density/iih1,iih2

      integer order
      common/nnlo/order
      integer nproc
      common/nproc/nproc

      double precision facscale 
      common/facscale/facscale

      double precision scale,musq
      common/scale/scale,musq

      double precision alphasPDF
      external alphasPDF
c     output
      double precision as
      common/asNEW/as
      integer ih1,ih2,nf,ic
      integer iter,locall,nlocall,iord,flagch
      common/flagch/flagch
      common/nf/nf
      double precision xmur,xmuf,xmur2,xmuf2
      common/scales2/xmur,xmuf,xmur2,xmuf2
      common/order/iord
      common/pdf/ih1,ih2
      common/vegint/iter,locall,nlocall
      double precision vud,vus,vub,vcd,vcs,vcb,vtd,vts,vtb
      common/ckm/vud,vus,vub,vcd,vcs,vcb,vtd,vts,vtb
      double precision amv, y1, y2, qtbis, gf, ppi, ssroot, sw2, aem,cw2
      common/cdyqt/amv,y1,y2,qtbis,gf,ppi,ssroot,sw2,aem,ic
      integer prodflag
      common/prodflag/prodflag
      double precision siggamma,sigint,sigz,sigw
      common/sigs/siggamma,sigint,sigz,sigw

      double precision yv,expyp,expym
      common/yv/yv,expyp,expym

c
      double precision vjfo
      double precision res

      gevpb=3.8937966d8
      ppi=dacos(-1d0)
      
      nf=5
      flagch=0

      iter=20
      locall=10000
      nlocall=20000
      
      ssroot=sroot
      ih1 = iih1
      ih2 = iih2

      iord=order-1
      xmur=scale
      xmuf=facscale
      xmur2=scale**2
      xmuf2=facscale**2
      as=alphasPDF(xmur)/ppi

      if (nproc.eq.3) then
         prodflag=5 !prodflag = 3
      else if (nproc.eq.1) then
         prodflag=21
      else if (nproc.eq.2) then
         prodflag=22
      endif
      
      vud=ud
      vus=us
      vub=ub
      vcd=cd
      vcs=cs
      vcb=cb
      vtd=0
      vts=0
      vtb=0

      sw2=xw
      cw2=1d0-sw2

      aemmz=sqrt(2d0)*gGf*(wmass**2)*xw/ppi
      aem=aemmz

      brz=1d0/zwidth*aem/24d0*zmass
     /  *((-1d0+2d0*sw2)**2+(2d0*sw2)**2)/(sw2*cw2) !=0.033638
      brw=1d0/wwidth*aem*wmass/(12d0*sw2)  !=0.10906 
      
      q2=m**2
      sigz=brz*zwidth*q2/(ppi*zmass)
     .     /((q2-zmass**2)**2+zmass**2*zwidth**2) !/(16d0*cw2)
      sigw=brw*wwidth*q2/(ppi*wmass)
     .     /((q2-wmass**2)**2+wmass**2*wwidth**2) !/4d0
      siggamma=aem/(3d0*ppi*q2)
      sigint=-aem/(6d0*ppi)*(q2-zmass**2)
     .     /((q2-zmass**2)**2+zmass**2*zwidth**2)*
     /        (-1d0+4d0*sw2)/(2d0*sqrt(sw2*cw2))               !sqrt(sw2/cw2)/2d0


      
      amv=m
      y1=y
      y2=y
      yv=y
      expyp=exp(y)
      expym=exp(-y)
      qtbis=pt

c      if (fnwa.eq.1) then       ! 
c      sigz=brz                  !/(16d0*cw2)
c      sigw=brw                  !/4d0
c      prodflag=3                !prodflag = 3
c      amv=91.1876
c      endif

c     y1=-5d0
c      y2=5d0
      call fodyqt(res)
c     conversion factors:
c     ds/dqt = ds/dqt2 * dqt2/dqt
c     fb = 1000 * pb
      vjfo = 2*pt*res*1000d0
      return
      end
