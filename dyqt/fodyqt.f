*****************************************************************************
c.....written by giuseppe bozzi
c.....This program computes dsigma/dqT^2 at NLO accuracy
c.....for (gamma,W,Z)+jet production at hadron colliders 
c.....Based on: Gonsalves, Pawlowski, Wai - Phys.Rev.D40 (1989)
c*****************************************************************************

      subroutine fodyqt(ris)
c      implicit none
c      double precision ris
c      double precision aem !amv, y1, y2, qtbis, gf, ppi, ssroot, sw2, aem,cw2
c      integer i 
      implicit real *8 (a-h,o-z)
      integer ih1,ih2,ic,isetproton,nloop,j,prodflag,nf,ord,nc
      integer qtpts,iter,locall,nlocall,iord,flagch
      common/flagch/flagch
      common/nf/nf
      common/vjorder/iord
      common/isetproton/isetproton
      common/pdf/ih1,ih2
      common/nloop/nloop
      common/para/s,ss
c      common/internal/qt,q,q2
      common/couplings/xw,cw,sw,alpha0
      common/const2/pi,cf,ca,tr,xnc
      common/gevpb/gevpb
      common/prodflag/prodflag
      common/vegint/iter,locall,nlocall
      common/ckm/vud,vus,vub,vcd,vcs,vcb,vtd,vts,vtb
      double precision amv, y1, y2, qtbis, gf, ppi, ssroot, sw2, aem,cw2
      common/cdyqt/amv,y1,y2,qtbis,gf,ppi,ssroot,sw2,aem,ic
      include 'fodyqt_inc.f'
      include 'internal_inc.f'

c.....definition of constants and couplings
      ca=3d0
      xnc=3d0
      cf=4/3d0
      tr=nf/2d0
      pi=ppi
      alpha0=aem
      xw=sw2
      sw=sqrt(xw)                       ! sin_w ! 
      cw=sqrt(1-xw)                     ! cos_w !

cc.....read the input file
c      if(ic.eq.1) then
c         ih1=1                      !beam 1 @ LHC!
c         ih2=1                      !beam 2 @ LHC!
c      elseif(ic.eq.-1) then
c         ih1=1                      !beam 1 @ Tevatron!
c         ih2=-1                     !beam 2 @ Tevatron!
c      endif
      ss=ssroot
      s=ss**2
      q=amv
      q2=q**2
      ord=iord
      qt=qtbis
      call qtdy(ris,error,chi2a,y1,y2,ord)
c      print *,q,qt,y1,y2,ris,error
      return
      end

c******************************************************
c.....this subroutine returns the value of dsigma/dqt^2
c******************************************************

      subroutine qtdy(ris,error,chi2a,y1,y2,ord)
c      implicit none
c      double precision ris,error,chi2a,y1,y2
      implicit real *8 (a-h,o-z)
      integer nf,ih1,ih2,ic,isetproton,prodflag
      integer ord,iord,iwseed,ndim,nloop,n3,ncl3
      integer iter,locall,nlocall
      external xdelta,sing
c      common/asNEW/as
c      common/asp/asp
      common/nf/nf
      common/nloop/nloop
c      common/scales2/xmur,xmuf,xmur2,xmuf2
      common/rapidities/yy1,yy2
c      double precision yv,expyp,expym
c      common/yv/yv,expyp,expym
      common/isetproton/isetproton
      common/pdf/ih1,ih2
      common/vjorder/iord
      common/para/s,ss
c      common/tm/tm
c      common/internal/qt,q,q2
      common/const2/pi,cf,ca,tr,xnc
c      common/mand/sh,th,uh,s2
      common/gevpb/gevpb
      common/prodflag/prodflag
      integer flagch
      common/flagch/flagch
      common/vegint/iter,locall,nlocall
      common/couplings/xw,cw,sw,alpha0
      include 'fodyqt_inc.f'
      include 'scales2_inc.f'
      include 'internal_inc.f'

      include 'gauss.f'
      integer i,j
      double precision xmin,xmax
      integer xrule,xintervals
      double precision ax,bx,cx,mx,x,t,jac
      real *8 xx(1:1)

      integer ii,jj,half
      double precision z1min,z1max
      integer z1rule,z1intervals
      double precision az1,bz1,cz1,mz1,z1,jacz1
      double precision z2min,z2max
      integer z2rule,z2intervals
      double precision az2,bz2,cz2,mz2,z2,jacz2
      real *8 zz(1:2)
      
!
      tm=sqrt(q2+qt**2)
c.....y2 not greater than y1
      if(y1.gt.y2) then
         ris=0d0
         return
      endif
      
c.....kinematical limits on qt
      z=q2/s
      xr=(1-z)**2-4*z*(qt/q)**2
      if(xr.lt.0)then
         ris=0d0
         return
      endif

cc.....kinematical limits on y
      tmpx=(q2+s)/ss/tm
      ymax=log((tmpx+sqrt(tmpx**2-4))/2)
c      if(y1.gt.ymax.or.y2.lt.(-ymax)) then
c         ris=0
c         return
c      elseif(y1.lt.-ymax.and.y2.gt.ymax) then
c         yy1=-ymax
c         yy2=ymax     
c      elseif(y1.lt.-ymax.and.y2.gt.-ymax.and.y2.lt.ymax) then
c         yy1=-ymax
c         yy2=y2
c      elseif(y1.gt.-ymax.and.y1.lt.ymax.and.y2.gt.ymax) then
c         yy1=y1
c         yy2=ymax
c      elseif(y1.gt.-ymax.and.y1.lt.ymax.and.y2.lt.ymax) then
c         yy1=y1
c         yy2=y2
c      endif

      if(yv.gt.ymax.or.yv.lt.(-ymax)) then
         ris=0
         return
      endif      
      
c.....compute as at order=nloop
      asp=as*pi

      if(qt/q.le.0.1d0) then
      ncl3=2*locall
      n3=2*iter
      else
      ncl3=locall
      n3=iter
      endif


c.....integration of delta(s2) terms
      ndim=1
!      n3=iter
!      ncl3=locall
      iwseed=1
      if (flagch.eq.3.or.flagch.eq.12) then
      rdelta=0d0
      error=0d0
      else
c      call xinteg(xdelta,ndim,n3,ncl3,av3,d3,chi2a)
c      rdelta=av3
c      error=d3
c      print *,rdelta,d3
      
c     boundaries of integration      
      xmin = 0
      xmax = 1

c     initialise
      rdelta = 0d0

      xintervals=1 !probably overdoing
      xrule=64
      
      do i=1,xintervals
         ax = xmin+(xmax-xmin)*(i-1)/xintervals
         bx = xmin+(xmax-xmin)*i/xintervals
         cx=0.5d0*(ax+bx)
         mx=0.5d0*(bx-ax)
         do j=1,xrule
            x=cx+mx*xxx(xrule,j)
            jac=1d0
            xx(1)=x
            rdelta = rdelta+xdelta(xx)*www(xrule,j)*jac*mx
c            print *,mx,bx,ax
         enddo
      enddo
      endif

c      do i=1,100
c         x=i/100d0
c         xx(1)=x
c         rdelta=xdelta(xx)
c      enddo
      


c.....integration of non-delta(s2) terms (only at NLO)
      if (iord.eq.1) then
c     *************** vegas integration *******************
c         ndim=2
c         if(qt/q.le.0.1d0) then
c            ncl3=2*nlocall
c            n3=2*iter
c         else
c            ncl3=nlocall
c            n3=iter
c         endif
c         call xinteg(sing,ndim,n3,ncl3,av3,d3,chi2a)
c         rsing=av3
c         error=error+d3
c         print *,(asp/2d0/pi)*rsing,(asp/2d0/pi)*d3
c     *****************************************************
         

c     initialise
         rsing = 0d0

c     boundaries of integration      
         z1min = 1D-12
         z1max = 1d0
      
         z1intervals=1
         z1rule=64

c     boundaries of integration      
         z2min = 1D-15
         z2max = 1d0

         z2intervals=1
         z2rule=64

         do i=1,z1intervals
            az1 = z1min+(z1max-z1min)*(i-1)/z1intervals
            bz1 = z1min+(z1max-z1min)*i/z1intervals
            cz1=0.5d0*(az1+bz1)
            mz1=0.5d0*(bz1-az1)
            do j=1,z1rule
               z1=cz1+mz1*xxx(z1rule,j)
c               t=z1
c               jacz1=1d0
c               t=1-z1min*(z1max/z1min)**z1
c               jacz1=(t-1)*log(z1max/z1min)
               t=z1min*(z1max/z1min)**z1
               jacz1=t*log(z1max/z1min)
c               esp=0.1d0
c               t=z1**esp
c               jacz1=esp*z1**(esp-1)
               
               zz(1)=t
               do ii=1,z2intervals
                  az2 = z2min+(z2max-z2min)*(ii-1)/z2intervals
                  bz2 = z2min+(z2max-z2min)*ii/z2intervals
                  cz2=0.5d0*(az2+bz2)
                  mz2=0.5d0*(bz2-az2)
                  do jj=1,z2rule
                     z2=cz2+mz2*xxx(z2rule,jj)
c                     t=z2
c                     jacz2=1d0
c                     t=1-z2min*(z2max/z2min)**z2
c                     jacz2=(t-1)*log(z2max/z2min)
                     t=z2min*(z2max/z2min)**z2
                     jacz2=t*log(z2max/z2min)
c                     esp=0.1d0
c                     t=z2**esp
c                     jacz2=esp*z2**(esp-1)
                     
                     zz(2)=t
                     rsing=rsing+sing(zz)*www(z1rule,j)*www(z2rule,jj)
     .                    *jacz1*jacz2*mz1*mz2
                  enddo
               enddo
            enddo
         enddo

c         print *,yv
c         print *,'{'
c         print *,'TGraph *gy = new TGraph();'
c         do i=0,1000
c            t=z1min*(z1max/z1min)**i/1000d0
c            jacz1=t*log(z1max/z1min)
c            zz(2)=i/1000d0
c            zz(1)=0.5
c            rsing=sing(zz)
c            print *,'gy->SetPoint(gy->GetN(), ',i/1000d0,', ',rsing,');'
cc            print *,zz(1),rsing
c         enddo
c         print *,'gy->Draw();'
c         print *,'}'
         
         rsing=rsing*(asp/2d0/pi)
c         print *,rsing
         
      elseif(iord.eq.0) then
        rsing=0d0
      endif

c.....final result
      ris=rdelta+rsing
c      print *,q,qt,yv,ris,error
c      print *,q,qt,yv,rdelta,rsing
      return
      end

c************************************************************
c.....this function includes all the delta(s2) contributions
c************************************************************

      function xdelta(xx)
      implicit real *8 (a-h,o-z)
      real *8 xx(1:1)
      real *8 fh1(-5:5),fh2(-5:5)
      integer nf,ih1,ih2,ic,isetproton,nloop,ord,iord,prodflag
      integer flagch
c      common/asp/asp
      common/nf/nf
c      common/scales2/xmur,xmuf,xmur2,xmuf2
      common/rapidities/y1,y2
c      double precision yv,expyp,expym
c      common/yv/yv,expyp,expym
      common/isetproton/isetproton
      common/pdf/ih1,ih2
      common/vjorder/iord
      common/para/s,ss
c      common/tm/tm
c      common/fractions/x1,x2
c      common/internal/qt,q,q2
      common/const2/pi,cf,ca,tr,xnc
c      common/mand/sh,th,uh,s2
      common/quarks/eq(5),alq(5),arq(5),ckm(6,6),delta(5,5),tau3(5,5)
      common/couplings/xw,cw,sw,alpha0
      common/gevpb/gevpb
      common/prodflag/prodflag
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
      common/flagch/flagch
      double precision ytot,expytotp,expytotm
      include 'fodyqt_inc.f'
      include 'internal_inc.f'
      include 'luminosities_inc.f'
      include 'scales2_inc.f'
c...........................................................
c.....definition of quantitites for phase space integration
c...........................................................
      xjac=1d0                  !jacobian "initialization"
      t=xx(1)                   !variable introduced to later allow 0<x2<1
      
c.....vector boson rapidity
c      if(y1.ne.y2) then
c         yv=y1+(y2-y1)*xx(2)
c         xjac=xjac*(y2-y1)
c      else
c         yv=y1
c     endif

c.....lower limit of integration for x2
      x2min=(q2-ss*tm*expym)/(ss*tm*expyp-s)
c      print *,'x2min', q,qt,y,x2min
c      print *,q2,ss,tm,expym,expyp,s
      if(x2min.gt.1d0.or.x2min.lt.0d0) then
         write(*,*) 'error in x2min'
         write(*,*) x2min
         write(*,*) 'm',q,'pt',qt,'y',yv
         xdelta=0d0
         return
c         stop
      endif

c.....change of variable to allow integration over 0-1
c....."esp" introduced for better behavior at small qt
      esp=8d0
      x2=dexp((1-t**esp)*log(x2min))
      xjac=xjac*x2*log(x2min)*esp*t**(esp-1d0)
c      x2=x2min*exp(1d0/x2min*t)
c      xjac = xjac*x2(1d0-x2min)
      
c.....definition of x1 in function of x2
      x1=(ss*tm*expyp*x2-q2)/(x2*s-ss*tm*expym)

c.....imposing x1,x2<1
      tiny=0d0 !1d-8
      if (x1.gt.1-tiny.or.x2.gt.1-tiny) then
         xdelta=0d0
         return
      endif
      
c.....definition of partonic mandelstam invariants and jacobian      
      sh=x1*x2*s
c      print *,'x1 x2 s',q,qt,y,x1,x2,s
c      print *,'q2, tm, sh',q,qt,y,sqrt(q2),tm,sh,q2/sh
c     bug fix in DYqT: th <-> uh
c      th=q2-ss*x2*tm*expyp
c      uh=q2-ss*x1*tm*expym
      th=q2-ss*x1*tm*expym
      uh=q2-ss*x2*tm*expyp
      xjj=abs(x2*s-ss*tm*expym)

c     1/xjj is the jacobian of the integration in dx1 of delta(q2-sh-uh-th) dx1
c     following from the change of variable t = q2-sh-uh-th
c     dx1 = dt * dx1(t)/dt  -> 1/xjj = dx1(t)/dt
      
c.....common factor for all the contributions
      factor=gevpb*alpha0*asp*cf/sh
c.....compute parton luminosity
      call flavour
      call utilities(sh,th,uh,q2,0d0)
c.....leading order delta(s2) contributions
      xloqg=(factor/(xnc**2-1d0))*(Aqg0(sh,th,uh,q2)+Agq0(sh,th,uh,q2))
      xloqqb=(factor/xnc)*Aqqb0(sh,th,uh,q2)
            if (flagch.eq.0) then
               xlo=xloqg+xloqqb 
            elseif (flagch.eq.1.or.flagch.eq.11.or.flagch.eq.12) then
               xlo=xloqqb 
            elseif (flagch.eq.2) then
               xlo=xloqg 
            elseif (flagch.eq.3) then
               xlo=0d0
            endif
            
c.....next to leading order delta(s2) contributions
      if (iord.eq.1) then
         xnloqg=(factor/(xnc**2-1d0))*(
     /     Bqg1(sh,th,uh,q2)+Bgq1(sh,th,uh,q2)+
     /     Bqg2(sh,th,uh,q2)+Bgq2(sh,th,uh,q2)+
     /     Cqg1(sh,th,uh,q2)+Cgq1(sh,th,uh,q2)+
     /     Cqg2(sh,th,uh,q2)+Cgq2(sh,th,uh,q2)+
     /     Bqg3(sh,th,uh,q2)+Bgq3(sh,th,uh,q2)         
     /        )
         xnloqqb=(factor/xnc)*(
     /        Bqqb1(sh,th,uh,q2)+Bqqb2(sh,th,uh,q2)+
     /        Cqqb1(sh,th,uh,q2)+D0aa(sh,th,uh,q2)+
     /        Bqqb3(sh,th,uh,q2)
     /        )

            if (flagch.eq.0) then
               xnlo=xnloqg+xnloqqb
            elseif (flagch.eq.1.or.flagch.eq.11.or.flagch.eq.12) then
               xnlo=xnloqqb 
            elseif (flagch.eq.2) then
               xnlo=xnloqg 
            elseif (flagch.eq.3) then
               xnlo=0d0 
            endif
      else
         xnlo=0d0
      endif

c.....total result
      xdelta=xlo+(asp/2/pi)*xnlo
      xdelta=-pi*xdelta*xjac/xjj
      return
      end
      
c**********************************************************
c.....this function includes all the singular contributions
c**********************************************************

      function sing(zz)
      implicit real *8 (a-h,o-z)
      real *8 zz(1:2),lb
      integer nf,ih1,ih2,ic,isetproton,nloop,ord,iord,prodflag
c      common/asp/asp
      common/nf/nf
c      common/scales2/xmur,xmuf,xmur2,xmuf2
      common/rapidities/y1,y2
c      double precision yv,expyp,expym
c      common/yv/yv,expyp,expym
      common/isetproton/isetproton
      common/pdf/ih1,ih2
      common/vjorder/iord
      common/para/s,ss
c      common/tm/tm
c      common/fractions/x1,x2
c      common/internal/qt,q,q2
      common/const2/pi,cf,ca,tr,xnc
c      common/mand/sh,th,uh,s2
      common/couplings/xw,cw,sw,alpha0
      common/gevpb/gevpb
      common/prodflag/prodflag
      integer flagch
      common/flagch/flagch
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
      include 'fodyqt_inc.f'
      include 'internal_inc.f'
      include 'luminosities_inc.f'
      include 'scales2_inc.f'

c.....definition of quantitites for phase space integration
c.....as defined in formula B.1 of Glosser-Schmidt paper (arXiv:hep-ph/0209248)
c     zz(1) <-> z1      
c     zz(2) <-> z2
c     zz(3) <-> yv
c to improve numerical accuracy, the change of variable z1 -> 1-z1 and z2 -> 1-z2 is applied

      xjacy=1d0                 !.....jacobian initialization

c.....vector boson rapidity
c      if(y1.ne.y2) then
c       yv=y1+(y2-y1)*zz(3)
c       xjacy=xjacy*(y2-y1)
c      elseif(y1.eq.y2) then
c       yv=y1
c      endif

c.....x1_0,x2_0 and modification by Massimiliano (dcut)
      x10=tm/ss*expyp
      x20=tm/ss*expym
      dcut=x10*(qt/tm)**2/(1-x10*(1-(qt/tm)**2))
c      dcut=qt/(tm+qt)
c      dcut=0d0

c.....z1,z2,lambda_b and corresponding jacobian 
c.....(to allow integration from 0 to 1)
      xjac=1d0
c******** changed z1 -> 1-z1 ***
c      z1=x20+(1-dcut-x20)*zz(1)
      z1=dcut+(1d0-dcut-x20)*zz(1)
c******** end ***
      xjac=xjac*(1d0-dcut-x20)

      
c******** changed z1 -> 1-z1 ***
c     lb=qt**2/tm**2*z1/(1-z1)
      lb=qt**2/tm**2*(1d0-z1)/z1
c******** end ***
c******** changed z2 -> 1-z2 ***
c      z2=x10*(1+lb)+(1-x10*(1+lb))*zz(2)
      z2=(1-x10*(1+lb))*zz(2)
c******** end ***
      xjac=xjac*(1-x10*(1+lb))     

c.....first term of the integrand for the plus prescription

c.....x1,x2
c******** changed z2 -> 1-z2 ***
c      x1=x10/z2*(1+lb)
      x1=x10/(1d0-z2)*(1+lb)
c******** end ***
c******** changed z1 -> 1-z1 ***
c      x2=x20/z1
      x2=x20/(1d0-z1)
c******** end ***

c.....definition of partonic Mandelstam invariants
c******** changed z2 -> 1-z2  z1 -> 1-z1 ***
c      sh=tm**2/z1/z2*(1+lb)
      sh=tm**2/(1d0-z1)/(1d0-z2)*(1+lb)
c******** end ***
c     bug fix in DYqT: th <-> uh
c      th=-tm**2/z1*(1-z1)*(1+lb)
c      uh=q2-tm**2/z2*(1+lb)
c******** changed z1 -> 1-z1 ***
c      uh=-tm**2/z1*(1-z1)*(1+lb)
      uh=-tm**2/(1d0-z1)*(z1)*(1+lb)
c******** end ***
c******** changed z2 -> 1-z2 ***
c      th=q2-tm**2/z2*(1+lb)
      th=q2-tm**2/(1d0-z2)*(1+lb)
c******** end ***
      
c      print *,'sh',sh,x1*x2*s
c      print *,'th',th,q2-ss*x1*tm*expym
c      print *,'uh',uh,q2-ss*x2*tm*expyp
      
c******** changed z2 -> 1-z2  z1 -> 1-z1 ***
c      s2=tm**2/z1/z2*(1-z2)*(1-z1)*(1+lb)
      s2=tm**2/(1d0-z1)/(1d0-z2)*z2*z1*(1+lb)
c******** end ***
c      print *,'s2',s2,sh+uh+th-q2

      if (s2.eq.0d0) then
         sing=0d0
         return
      endif
      
c.....imposing s2>0
      if(s2.lt.0) then
         write(*,*)'s2 < 0 ! s2 =',s2
         print *,'z1',z1,'z2',z2,'x1',x1,'x2',x2,'sh',sh,'th',th,'uh',uh
         write(*,*) 'm',q,'pt',qt,'y',yv
         xdelta=0d0
         return
c       stop
      endif

c.....phase space prefactor
c******** changed z2 -> 1-z2  z1 -> 1-z1 ***
c      pre1=tm**2*(1+lb)/(z1*z2)**2
      pre1=tm**2*(1+lb)/((1d0-z1)*(1d0-z2))**2
c******** end ***

c.....imposing x1,x2 < 1
      tiny=0d0 !1d-7
      if(x1.gt.1-tiny.or.x2.gt.1-tiny) then
         write(*,*)'x1',x1,'x2',x2
         sing2=0d0
         sing1=0d0
         reg=0d0
         xspur=0d0
         xm10=0
         xdelta=0d0
         return
      else

c.....compute parton luminosity
         call flavour
         call utilities(sh,th,uh,q2,s2)

c.....common factor for all the contributions
         factor=gevpb*alpha0*asp*cf/sh

c.....regular contributions (non-delta and non-singular in s2)
         xrgg=(factor*(xnc/(xnc**2-1d0)**2))*(
     /        Cgg1(sh,th,uh,q2)+Cgg1x(sh,th,uh,q2)
     /        )
         xrqg=(factor/(xnc**2-1d0))*(
     /        Cqg3(sh,th,uh,q2,3)+Cgq3(sh,th,uh,q2,3)+
     /        Cqg3(sh,th,uh,q2,4)+Cgq3(sh,th,uh,q2,4)
     /        )
         xrqqb=(factor/xnc)*(
     /        Cqqb2(sh,th,uh,q2,4)+Cqqb2x(sh,th,uh,q2,4)+
     /        Daa(sh,th,uh,q2,4)+Daax(sh,th,uh,q2,4)+
     /        Dab(sh,th,uh,q2)+Dabx(sh,th,uh,q2)+
     /        Dbb(sh,th,uh,q2)+Dbbx(sh,th,uh,q2)+
     /        Dac(sh,th,uh,q2)+Dad(sh,th,uh,q2)+
     /        Dbc(sh,th,uh,q2)+Dbd(sh,th,uh,q2)+
     /        Dcc(sh,th,uh,q2)+Ddd(sh,th,uh,q2)+
     /        DcdLL(sh,th,uh,q2)+DcdLLx(sh,th,uh,q2)+
     /        DcdLR(sh,th,uh,q2)+DcdLRx(sh,th,uh,q2)
     /        )
         xrqq=(factor/xnc)*(1/2d0)*(
     /        Eac(sh,th,uh,q2)+Ebd(sh,th,uh,q2) +
     /        Ead(sh,th,uh,q2)+Ebc(sh,th,uh,q2)
c     missing pieces in qq
c     /        + Eaa(sh,th,uh,q2)+Ecc(sh,th,uh,q2) ! -> =2*Dcc
c     /        + Ebb(sh,th,uh,q2)+Edd(sh,th,uh,q2) ! -> =2*Ddd
c     /        + EabLL+EabLLx
c     /        + EabLR+EabLRx
c     /        + EcdLL+EcdLLx !-DcdLR
c     /        + EcdLR+EcdLRx !-DcdLL
     /        )
         
            if (flagch.eq.0) then
         xreg=xrgg+xrqg+xrqqb+xrqq
            elseif (flagch.eq.1.or.flagch.eq.11.or.flagch.eq.12) then
         xreg=xrqqb+xrqq
            elseif (flagch.eq.2) then
         xreg=xrqg
            elseif (flagch.eq.3) then
         xreg=xrgg
            endif
      if (xrqg.ne.xrqg) then
         print *,'nan in xrqg, s2 = ',s2
         xrqg = 0d0
      endif

      if (xrqqb.ne.xrqqb) then
         print *,'nan in xrqqb, s2 = ',s2
         xrqqb = 0d0
      endif
c       print *,xreg,factor,sh,th,uh,q2,s2

c.....(log(s2)/s2)_A+ contributions
         xsqg1=(factor/(xnc**2-1d0))*
     /        (Cqg3(sh,th,uh,q2,1)+Cgq3(sh,th,uh,q2,1))
         xsqqb1=(factor/xnc)*
     /        (Cqqb2(sh,th,uh,q2,1)+Cqqb2x(sh,th,uh,q2,1))
            if (flagch.eq.0) then
         xs1=xsqg1+xsqqb1
            elseif (flagch.eq.1.or.flagch.eq.11.or.flagch.eq.12) then
         xs1=xsqqb1
            elseif (flagch.eq.2) then
         xs1=xsqg1
            elseif (flagch.eq.3) then
          xs1=0d0
            endif

c.....(1/s2)_A+ contributions
         xsqg2=(factor/(xnc**2-1d0))*
     /        (Cqg3(sh,th,uh,q2,2)+Cgq3(sh,th,uh,q2,2))
         xsqqb2=(factor/xnc)*(
     /        Cqqb2(sh,th,uh,q2,2)+Cqqb2x(sh,th,uh,q2,2)+
     /        Daa(sh,th,uh,q2,2)+Daax(sh,th,uh,q2,2)
     /        )
            if (flagch.eq.0) then
         xs2=xsqg2+xsqqb2
            elseif (flagch.eq.1.or.flagch.eq.11.or.flagch.eq.12) then
         xs2=xsqqb2
            elseif (flagch.eq.2) then
         xs2=xsqg2
            elseif (flagch.eq.3) then
         xs2=0d0
            endif

      endif 

c     *********************************************************
c     This piece depends only on yv and z1, phase space is the same, with z2=1 (which leads to s2=0)
c     It corresponds to fs(0) in Eq. (2.5) of [Gonsalves, Pawlowsky, Wai]
      
c.....subtraction of z2=1 term according to plus prescription     

c.....x1,x2
      x1=x10*(1+lb)
c******** changed z1 -> 1-z1 ***
c      x2=x20/z1
      x2=x20/(1d0-z1)
c******** end ***

c.....definition of partonic Mandelstam invariants
c******** changed z1 -> 1-z1 ***
c      sh=tm**2/z1*(1+lb)
      sh=tm**2/(1d0-z1)*(1+lb)
c******** end ***
c     bug fix in DYqT: th <-> uh
c      th=-tm**2/z1*(1-z1)*(1+lb)
c      uh=q2-tm**2*(1+lb)
c******** changed z1 -> 1-z1 ***
c      uh=-tm**2/z1*(1-z1)*(1+lb)
      uh=-tm**2/(1d0-z1)*z1*(1+lb)
c******** end ***
      th=q2-tm**2*(1+lb)
      
c      print *,sh,x1*x2*s,th,q2-ss*x1*tm*expym,uh,q2-ss*x2*tm*expyp

c.....phase space prefactor
c******** changed z1 -> 1-z1 ***
c      pre10=tm**2*(1+lb)/(z1)**2
      pre10=tm**2*(1+lb)/(1d0-z1)**2
c******** end ***

c.....imposing x1,x2 < 1
      if(x1.gt.1-tiny.or.x2.gt.1-tiny) then
         sing2=0d0
         sing1=0d0
         reg=0d0
         xspur=0d0
         xm10=0
      else

c.....compute parton luminosity
         call flavour
c     set s2=0d0 explicitly in utilities
         call utilities(sh,th,uh,q2,0d0)
         
c.....common factor for all the contributions
         factor=gevpb*alpha0*asp*cf/sh

c.....(log(s2)/s2)_A+ contributions
         xsqg10=(factor/(xnc**2-1d0))*
     /        (Cqg3(sh,th,uh,q2,1)+Cgq3(sh,th,uh,q2,1))
         xsqqb10=(factor/xnc)*
     /        (Cqqb2(sh,th,uh,q2,1)+Cqqb2x(sh,th,uh,q2,1))
            if (flagch.eq.0) then
         xs10=xsqg10+xsqqb10
            elseif (flagch.eq.1.or.flagch.eq.11.or.flagch.eq.12) then
         xs10=xsqqb10
            elseif (flagch.eq.2) then
         xs10=xsqg10
            elseif (flagch.eq.3) then
         xs10=0d0
            endif

c.....(1/s2)_A+ contributions
         xsqg20=(factor/(xnc**2-1d0))*
     /        (Cqg3(sh,th,uh,q2,2)+Cgq3(sh,th,uh,q2,2))
         xsqqb20=(factor/xnc)*(
     /        Cqqb2(sh,th,uh,q2,2)+Cqqb2x(sh,th,uh,q2,2)+
     /        Daa(sh,th,uh,q2,2)+Daax(sh,th,uh,q2,2)
     /        )
            if (flagch.eq.0) then
         xs20=xsqg20+xsqqb20
            elseif (flagch.eq.1.or.flagch.eq.11.or.flagch.eq.12) then
         xs20=xsqqb20
            elseif (flagch.eq.2) then
         xs20=xsqg20
            elseif (flagch.eq.3) then
         xs20=0d0
            endif
c     *********************************************************


c.....final result.................................................
c.....note: here I am passing from s2 to z2 distributions..........
c.....according to the formulas:...................................
c.....(1/s2)_A+ = (1/-t)*((1/(1-z2))_+ -1 
c.....                   + log(z2min/(1-z2min))*delta(1-z2))
c.....(log(s2/q2)/s2)_A+ = log(-t/q2)*(1/s2)_A+ - 1/t*(
c.....                   (log(1-z2)/(1-z2))_+ + log(z2/(1-z2))
c.....                   - log(z2)/(1-z2) - 1/2*delta(1-z2)*
c.....                   log^2((1-z2min)/z2min))
c.....as a consequence, there appear some "spurious" terms.........
c.....in nlo delta, 1/s2 and regular contributions.................     
c.....I do not include spurious delta terms: they cancel against...
c.....similar contributions in subroutine 'xdelta'.................
c.....provided that A=-th in Gonsalves' formulas...................
c.....I also add mismatch terms (xmism) coming from 0-z2min region.

c     bug fix in DYqT: th <-> uh
c      a=-th
      a=-uh
         z2min=x10*(1+lb)   

         xmism=-0.5d0*xs10*(log(1-z2min))**2-
     /        (xs10*log(a/q2)+xs20)*log(1-z2min)

         xmism=xmism*pre10*(1-dcut-x20)/(a)

c******** changed z2 -> 1-z2 ***
c     sing1=log(1-z2)/(1-z2)*(pre1*xs1-pre10*xs10)*xjac/(a)
         sing1=log(z2)/(z2)*(pre1*xs1-pre10*xs10)*xjac/(a)
c******** end ***

c******** changed z2 -> 1-z2 ***
c         sing2=(pre1*(xs1*log(a/q2)+xs2)-
c     /        pre10*(xs10*log(a/q2)+xs20))/(1-z2)*xjac/(a)
         sing2=(pre1*(xs1*log(a/q2)+xs2)-
     /        pre10*(xs10*log(a/q2)+xs20))/(z2)*xjac/(a)
c******** end ***

c******** changed z2 -> 1-z2 ***
c         reg=pre1*(xreg+
c     /        (xs1*(log(z2/(1-z2))-(log(z2))/(1-z2)-log(a/q2))-
c     /        xs2)/(a))*xjac
         reg=pre1*(xreg+
     /        (xs1*(log((1d0-z2)/(z2))-(log(1d0-z2))/(z2)-log(a/q2))-
     /        xs2)/(a))*xjac
c******** end ***
         
      endif
      
      sing=(sing1+sing2+reg-xmism)/s
!      write(*,*) reg
      sing=pi*xjacy*sing
      if (sing.ne.sing) then
         print *,'nan in xsing', sing1,sing2,reg,xmism
         sing=0d0
      endif
      return
      end

