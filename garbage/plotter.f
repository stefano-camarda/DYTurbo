      subroutine bookplot(n,tag,titlex,var,wt,xmin,xmax,dx,llplot) 
      implicit none
      include 'nplot.f'
      integer n
      character*(*) titlex
      character*3 llplot
      character*4 tag,mypart
      double precision var,wt,xmin,xmax,dx
      common/mypart/mypart

      if (tag.eq.'book') then
          call mbook(n,titlex,dx,xmin,xmax)
          call mbook(20+n,titlex,dx,xmin,xmax)
          call mbook(40+n,titlex,dx,xmin,xmax)
          call mbook(60+n,titlex,dx,xmin,xmax)
          call mbook(80+n,titlex,dx,xmin,xmax)
          call mbook(100+n,titlex,dx,xmin,xmax)
          call mbook(120+n,titlex,dx,xmin,xmax)
      elseif (tag .eq. 'plot') then
          call mfill(n,var,wt)
        linlog(n)=llplot
        titlearray(n)=titlex
      endif

      return
      end

    
      subroutine plotter(p,wt,switch)
      implicit none
      include 'clustering.f'
      include 'constants.f'
      include 'cutoff.f'
      include 'jetlabel.f'
      include 'npart.f'
      include 'mxdim.f'
      include 'process.f'
      include 'removebr.f'
      include 'masses.f'


      character*2 plabel(mxpart)
      common/plabel/plabel

      integer n,switch,nplotmax,i,j
      character tag*4
  
      double precision wt
      double precision m34,p(mxpart,4),fphi,HT
      double precision pt3,pt4,pt5,pt6,pt7,pt8
      double precision eta3,eta4,eta5,eta6,eta7,eta8
      double precision pt34,y34,pjm,pt3dpt4,tmass,ptmin,ptmax,pt34om34
      double precision pt34y01,pt34y12,pt34y224
      double precision m34y01,m34y12,m34y224
      double precision pt34aa,pt34a,pt34b,pt34c,pt34d,pt34e,pt34f
      double precision pt340,pt341,pt342,pt343,pt344
      double precision pt4a,pt4b,pt4c,pt4d
      double precision eta4a,eta4b,eta4c,eta4d,eta4e,eta4f
      double precision pt,etarap,yraptwo,yrapfour,pttwo,R
      double precision cosphi34,deltaphi,pto(1:4),tmp(1:2),costh_CS
      integer eventpart,nqcdjets,nqcdstart
      double precision phi3,phi4,cosdelphi,delphi,phi_acop,theta_star
     &,phi_st


      logical first,jetmerge
      character*30 runstring
      common/runstring/runstring
      common/nplotmax/nplotmax
      common/nqcdjets/nqcdjets,nqcdstart
      common/jetmerge/jetmerge
      double precision realeventp(mxpart,4)
      common/realeventp/realeventp
      integer order,nproc,ndec
      common/nproc/nproc
      common/nnlo/order
      data first/.true./
      save first
       if (first) then
        tag='book'
c--- ensure we initialize all possible histograms
        eventpart=npart+3
        eta3=0d0
        pt3=0d0
        eta4=0d0
        pt4=0d0
        eta5=0d0
        pt5=0d0
        eta6=0d0
        pt6=0d0
        eta7=0d0
        pt7=0d0
        eta8=0d0
        pt8=0d0

        y34=0d0
        pt34=0d0
        pt34om34=0d0
        pt340=0d0
        pt341=0d0
        pt342=0d0
        pt343=0d0
        pt344=0d0
        pt34aa=0d0
        pt34a=0d0
        pt34b=0d0
        pt34c=0d0
        pt34d=0d0
        pt34e=0d0
        pt34f=0d0
        pt4a=0d0
        pt4b=0d0
        pt4c=0d0
        pt4d=0d0
        eta4a=0d0
        eta4b=0d0
        eta4c=0d0
        eta4d=0d0
        eta4e=0d0
        eta4f=0d0
        m34=0d0
        cosphi34=0d0
        costh_CS=0d0

        pt34y01=0d0
        pt34y12=0d0
        pt34y224=0d0
        m34y01=0d0
        m34y12=0d0
        m34y224=0d0

        HT=0d0

        deltaphi=0d0
        phi3=0d0
        phi4=0d0
        cosdelphi=0d0
        delphi=0d0
        phi_acop=0d0
        theta_star=0d0
        phi_st=0d0

        pto(1)=0d0
        pto(2)=0d0
        pto(3)=0d0
        pto(4)=0d0


        jetmerge=.true.
CC      Here set jets to the maximum number of jets 
CC      to book the necessary histograms: 0 at LO, 1 at NLO and 2 at NNLO
        jets=order
CC
        goto 99
      else
        tag='plot'
      endif

C     ndec is the number of decay products of the vector boson

      ndec=2

C     enventpart is the total number of four momenta in the event
C     2 for initial state + ndec=2 for the V decay + number of jets 

      eventpart=2+ndec+jets
    

      eta3=etarap(3,p)
      pt3=pt(3,p)
      eta4=etarap(4,p)
      pt4=pt(4,p)        
      y34=yraptwo(3,4,p)
      pt34=pttwo(3,4,p)
      m34=dsqrt((p(3,4)+p(4,4))**2-(p(3,1)+p(4,1))**2
     .          -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2)

      pt34om34=pt34/m34

      pt34y01=-10000d0
      m34y01=-10000d0
      pt34y12=-10000d0
      m34y12=-10000d0
      pt34y224=-10000d0
      m34y224=-10000d0

      if (abs(y34).lt.1d0) then
         pt34y01=pt34
         m34y01=m34
      endif
      if (abs(y34).ge.1d0.and.abs(y34).lt.2d0) then
         pt34y12=pt34
         m34y12=m34
      endif
      if (abs(y34).ge.2d0.and.abs(y34).lt.2.4d0) then
         pt34y224=pt34
         m34y224=m34
      endif


      HT=pt3+pt4

C     Transverse mass

      pt3dpt4=p(3,1)*p(4,1)+p(3,2)*p(4,2)

      if(pt3*pt4.ge.pt3dpt4) then
      tmass=dsqrt(2d0*(pt3*pt4-pt3dpt4))
      else
      tmass=-999d0  ! do not bin
      write(*,*) "WARNING! tmass^2=",2d0*(pt3*pt4-pt3dpt4)
      endif


      ptmin=min(pt3,pt4)
      ptmax=max(pt3,pt4)


!  WRITE EVENTS ON FILE
      call evtwriter(p,wt,switch)


 99   continue

      n=1                  
      
      if(nproc.eq.3) then
         call bookplot(n,tag,'pt34',pt34,wt,0d0,100d0,1d0,'lin')
         n=n+1
         call bookplot(n,tag,'pt34y01',pt34y01,wt,0d0,100d0,1d0,'lin')
         n=n+1
         call bookplot(n,tag,'pt34y12',pt34y12,wt,0d0,100d0,1d0,'lin')
         n=n+1
         call bookplot(n,tag,'pt34y224',pt34y224,wt,0d0,100d0,1d0,'lin')
         n=n+1
         call bookplot(n,tag,'m34',m34,wt,60d0,120d0,1d0,'lin')
         n=n+1
         call bookplot(n,tag,'m34y01',m34y01,wt,60d0,120d0,1d0,'lin')
         n=n+1
         call bookplot(n,tag,'m34y12',m34y12,wt,60d0,120d0,1d0,'lin')
         n=n+1
         call bookplot(n,tag,'m34y224',m34y224,wt,60d0,120d0,1d0,'lin')
         n=n+1
      endif

      if(nproc.eq.1) then
         call bookplot(n,tag,'pte',pt4,wt,30d0,50d0,0.5d0,'lin')
         n=n+1
         call bookplot(n,tag,'ptmiss',pt3,wt,0d0,100d0,1d0,'lin')
         n=n+1
         call bookplot(n,tag,'pt34',pt34,wt,0d0,100d0,1d0,'lin')
         n=n+1
         call bookplot(n,tag,'tmass',tmass,wt,60d0,110d0,1d0,'lin')
         n=n+1
         call bookplot(n,tag,'etae',eta4,wt,-5d0,5d0,0.5d0,'lin')
         n=n+1
         call bookplot(n,tag,'etanu',eta3,wt,-5d0,5d0,0.5d0,'lin')
         n=n+1
      elseif(nproc.eq.2) then
         call bookplot(n,tag,'pte',pt3,wt,30d0,50d0,0.5d0,'lin')
         n=n+1
         call bookplot(n,tag,'ptmiss',pt4,wt,0d0,100d0,1d0,'lin')
         n=n+1
         call bookplot(n,tag,'pt34',pt34,wt,0d0,100d0,1d0,'lin')
         n=n+1
         call bookplot(n,tag,'tmass',tmass,wt,60d0,110d0,1d0,'lin')
         n=n+1
         call bookplot(n,tag,'etae',eta3,wt,-5d0,5d0,0.5d0,'lin')
         n=n+1
         call bookplot(n,tag,'etanu',eta4,wt,-5d0,5d0,0.5d0,'lin')
         n=n+1
      endif



    
      n=n-1

      if (n .gt. 20) then
        write(6,*) 'WARNING - TOO MANY HISTOGRAMS!'
        write(6,*) n,' > 20, which is the built-in maximum'
        stop
      endif

c--- set the maximum number of plots, on the first call
      if (first) then
        first=.false.
        nplotmax=n
      endif
      

      return 
      end
      

      
