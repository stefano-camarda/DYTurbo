c**********************************************
c.....LO contributions for QG->V+X and GQ->V+X
c**********************************************

c.....Real contribution to QG->V+X
      function Aqg0(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
      include 'luminosities_inc.f'
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR

      temp=(-sh/th-th/sh-(2d0*uh*q2)/(sh*th))
      Aqg0=temp*xlumqg
      
      return
      end

c.....Real contribution to GQ->V+X
      function Agq0(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
      include 'luminosities_inc.f'
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR

      temp=(-sh/uh-uh/sh-(2d0*th*q2)/(sh*uh))
      Agq0=temp*xlumgq

      return
      end

c**********************************************
c.....LO contribution for QQB->V+X (=QBQ->V+X)
c**********************************************

c.....Real contribution to QQB->V+X
      function Aqqb0(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
      include 'luminosities_inc.f'
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR

      temp=(uh/th+th/uh+(2d0*sh*q2)/(uh*th))
      Aqqb0=temp*xlumqqb

      return
      end

c******************************************************
c.....NLO virtual contributions to QG->V+X and GQ->V+X
c******************************************************

c.....Virtual (one-loop) contribution to QG->V+X
      function Bqg1(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
      
      temp=(sh/th+th/sh+(2d0*uh*q2)/(sh*th))*
     /     (-8*cf-cf*fu**2-pi**2/3d0*(cf-ca)+
     /     1/2d0*ca*(fu**2-fs**2-ft**2)+ca*(11/6d0*fmu2+f1t-f2t))+
     /     cf*(uh/(th+uh)+uh/(sh+uh)+(uh+th)/sh+(uh+sh)/th)+
     /     ft*(cf*(4*uh**2+2*uh*th+4*sh*uh+sh*th)/(sh+uh)**2+
     /     (ca*th)/(sh+uh))+
     /     fs*(cf*(4*uh**2+2*sh*uh+4*uh*th+th*sh)/(th+uh)**2+
     /     (ca*sh)/(th+uh))+
     /     (2*cf-ca)*(2*fu*(uh**2/(sh+th)**2+2*uh/(sh+th))-q2/(sh*th)*
     /     (sh**2+th**2)/(sh+th)+(uh**2+(uh+sh)**2)/(sh*th)*(f1t+f1u-
     /     ft*fu+pi**2/2d0)+(uh**2+(th+uh)**2)/(sh*th)*(f1u-f2u))
      Bqg1=-temp*xlumqg
      
      return
      end

c.....Virtual (one-loop) contribution to GQ->V+X
      function Bgq1(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
      
      temp=(sh/uh+uh/sh+(2d0*th*q2)/(sh*uh))*
     /     (-8*cf-cf*ft**2-pi**2/3d0*(cf-ca)+
     /     1/2d0*ca*(ft**2-fs**2-fu**2)+ca*(11/6d0*fmu2+f1u-f2u))+
     /     cf*(th/(uh+th)+th/(sh+th)+(th+uh)/sh+(th+sh)/uh)+
     /     fu*(cf*(4*th**2+2*th*uh+4*sh*th+sh*uh)/(sh+th)**2+
     /     (ca*uh)/(sh+th))+
     /     fs*(cf*(4*th**2+2*sh*th+4*th*uh+uh*sh)/(uh+th)**2+
     /     (ca*sh)/(uh+th))+
     /     (2*cf-ca)*(2*ft*(th**2/(sh+uh)**2+2*th/(sh+uh))-q2/(sh*uh)*
     /     (sh**2+uh**2)/(sh+uh)+(th**2+(th+sh)**2)/(sh*uh)*(f1u+f1t-
     /     fu*ft+pi**2/2d0)+(th**2+(uh+th)**2)/(sh*uh)*(f1t-f2t))
      Bgq1=-temp*xlumgq

      return
      end

c.....Virtual (renormalization counterterm) contribution to QG->V+X
      function Bqg2(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
      
      temp=(sh/th+th/sh+(2d0*uh*q2)/(sh*th))*fmu2/3d0*nf
      Bqg2=temp*xlumqg

      return
      end

c.....Virtual (renormalization counterterm) contribution to GQ->V+X
      function Bgq2(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
      
      temp=(sh/uh+uh/sh+(2d0*th*q2)/(sh*uh))*fmu2/3d0*nf
      Bgq2=temp*xlumgq

      return
      end

c.....Virtual (triangular quark loops) contribution to QG->V+X
      function Bqg3(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
      
      temp=(uh+q2)/(uh-q2)*(1-(q2*fu)/(uh-q2))
      Bqg3=temp*xlumqgtr

      return
      end

c.....Virtual (triangular quark loops) contribution to GQ->V+X
      function Bgq3(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
      
      temp=(th+q2)/(th-q2)*(1-(q2*ft)/(th-q2))
      Bgq3=temp*xlumgqtr

      return
      end

c********************************************************
c.....NLO virtual contributions for QQB->V+X (=QBQ->V+X)
c********************************************************

c.....Virtual (one-loop) contribution to QQB->V+X
      function Bqqb1(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
      
      temp=(uh/th+th/uh+(2d0*sh*q2)/(uh*th))*
     /     (-8*cf-cf*fs**2+pi**2/6d0*(4*cf-ca)+
     /     1/2d0*ca*(fs**2-(ft+fu)**2)+ca*(11/6d0*fmu2+f1t+f1u))+
     /     cf*(sh/(th+sh)+sh/(sh+uh)+(sh+th)/uh+(uh+sh)/th)+
     /     ft*(cf*(4*sh**2+2*sh*th+4*sh*uh+uh*th)/(sh+uh)**2+
     /     (ca*th)/(sh+uh))+
     /     fu*(cf*(4*sh**2+2*sh*uh+4*sh*th+th*uh)/(th+sh)**2+
     /     (ca*uh)/(th+sh))+
     /     (2*cf-ca)*(2*fs*(sh**2/(uh+th)**2+2*sh/(uh+th))-q2/(uh*th)*
     /     (uh**2+th**2)/(uh+th)+(sh**2+(uh+sh)**2)/(uh*th)*(f1t-f2t)+
     /     (sh**2+(th+sh)**2)/(uh*th)*(f1u-f2u))
      Bqqb1=temp*xlumqqb

      return
      end

c.....Virtual (renormalization counterterm) contribution to QQB->V+X
      function Bqqb2(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
      
      temp=-(uh/th+th/uh+(2d0*sh*q2)/(uh*th))*fmu2/3d0*nf
      Bqqb2=temp*xlumqqb

      return
      end

c.....Virtual (triangular quark loops) contribution to QQB->V+X
      function Bqqb3(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
      real *8 fh1(-5:5),fh2(-5:5)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
      
      temp=-(sh+q2)/(sh-q2)*(1-(q2*fs)/(sh-q2))
      Bqqb3=temp*xlumqqbtr

      return
      end

c***************************************************
c.....NLO real contributions to QG->VQG and GQ->VQG
c***************************************************

c.....First Real contribution to QG->VQG (proportional to delta(s2))
      function Cqg1(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
      
      temp=(-sh/th-th/sh-(2d0*uh*q2)/(sh*th))*
     /     (cf*(7/2d0+2*fm2*(fu-fa)+fa**2-3/2d0*(fm2+fa))+
     /     ca*(pi**2/6d0+1/2d0*(fs-ft-fu)**2+2*ft*(fm2-fa)+
     /     2*fa*(fs-fu-fm2+fa)-11/6d0*fm2))
      Cqg1=temp*xlumqg

      return
      end

c.....First Real contribution to GQ->VQG (proportional to delta(s2))
      function Cgq1(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
      
      temp=(-sh/uh-uh/sh-(2d0*th*q2)/(sh*uh))*
     /     (cf*(7/2d0+2*fm2*(ft-fa)+fa**2-3/2d0*(fm2+fa))+
     /     ca*(pi**2/6d0+1/2d0*(fs-fu-ft)**2+2*fu*(fm2-fa)+
     /     2*fa*(fs-ft-fm2+fa)-11/6d0*fm2))
      Cgq1=temp*xlumgq

      return
      end

c.....Second Real contribution to QG->VQG (proportional to delta(s2))
      function Cqg2(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
      
      temp=(-sh/th-th/sh-(2d0*uh*q2)/(sh*th))*fm2/3d0*nf
      Cqg2=temp*xlumqg

      return
      end

c.....Second Real contribution to GQ->VQG (proportional to delta(s2))
      function Cgq2(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
      
      temp=(-sh/uh-uh/sh-(2d0*th*q2)/(sh*uh))*fm2/3d0*nf
      Cgq2=temp*xlumgq

      return
      end

c.....Third Real contribution to QG->VQG (for nonzero values of s2)
      function Cqg3(sh,th,uh,q2,what)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton,what
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
c      s2=sh+th+uh-q2

      if (what.eq.1) then
c     term proportional to (fs2/s2)_A+
      temp1=(-sh/th-th/sh-(2d0*uh*q2)/(sh*th))*(2d0*cf+4d0*ca)
      temp=temp1
      elseif (what.eq.2) then
c     term proportional to (1/s2)_A+
      temp1=(-sh/th-th/sh-(2d0*uh*q2)/(sh*th))*(
     /     cf*(-3/2d0+fsu+2d0*ftu-2d0*fm2+(th+uh)*fla/la)+
     /     ca*(2d0*fstu+(fst-fsu)/2d0-ftu-2d0*fm2))
      temp=temp1
      elseif (what.eq.3) then
c     term proportional to (1/s2), FINITE in the limit s2->0
      temp1=(cf-ca/2d0)*((sh+2d0*uh)/th*((th+uh)*fla/la+fsu)+
     /     (2d0*th+4d0*uh)*ftu/sh)
      temp=temp1/s2
      elseif (what.eq.4) then
c     regular term
      temp1=fla/la**3*(cf*((3*sh)*(sh+q2-s2)*(uh**2-th**2)*
     /     (1-uh/th)/(4*la**2)+(uh**2-th**2)*((q2-s2**2)/(4*sh)+
     /     11/4d0+s2/th-uh/(2*th))-2*sh*(th+uh**2/th)-
     /     sh/2d0*(sh-s2)*(uh/th-3d0)+4*s2*(th-uh)+5*sh*uh)+
     /     ca*(1-uh/th)*((th+uh)/4d0*(sh+q2-s2)+sh*(sh-s2)))+
     /     fla/la*(cf*((th-uh-2*s2)/(4*sh)+(11*sh-2*uh-4*s2)/(4*th)+
     /     2d0/sh*(1+uh/th)*(5*q2+s2)-(16d0*q2*s2)/(sh*th))-
     /     ca/2d0*((th+uh)/(2d0*th)+5d0*(uh/sh-(uh+s2)/th)-
     /     7d0/th*(sh+uh*s2/sh)+3d0/sh*(th-3*s2)+
     /     2d0/(sh*th)*(uh**2+3*s2**2)))+
     /     1/la**2*(cf*(3*sh/(2*la**2)*(th-uh)**2*(1+uh/th)+
     /     sh/(2*th)*(uh-3*th)+(1-uh/th)*(s2-7d0*uh/4d0)+
     /     (th-uh)*(11/4d0+(3*th+3*uh)/(2*sh)-2*s2/sh))+
     /     ca*(uh/th-1)*(sh+s2))+
     /     2*dt**3*sh*th*(4*cf-3*ca-(cf-ca)*(fs2-fm2))+
     /     dt**2*(cf*(3*sh-8*th-4*uh)-ca*(9*sh-4*th+4*uh)+
     /     2d0*(fs2-fm2)*(cf*(th+uh)+ca*(3*sh+2*uh)))+
     /     dt*((cf*(fs2-fm2)-ca*(fs2-fm2+(fstu+fs2-ftu)/2d0+
     /     (fst+flat)/4d0))*(4d0*(1+uh/sh)*(1-uh/th)-
     /     2d0*(sh/th+th/sh))-cf*(fs2-fm2)*(sh+2*uh)/th+
     /     ca*((2*flat+2*fst)*(uh+sh)/th+2d0*(fs2-fm2)*
     /     ((5*sh+4*uh)/th-2))-cf*(4+(sh-2*uh)/th+(2*uh-3*th)/sh-
     /     uh**2/(sh*th))+ca*(7-(3*sh+2*uh)/th))+
     /     du*(2*cf*(fm2-fs2-sh/th)+ca*(sh/th-1))+
     /     dst*cf*(fsu-2*ftu-flat)*(1-s2/th)*(s2**2-2*uh*(s2-uh))/sh**2+
     /     dtu*cf*(2d0*(fm2-fs2-ftu)*(2*sh/th*(q2+2*uh)+th-s2+
     /     4*uh/th*(2*uh-s2)+2d0/sh*(uh-s2)*(th-2*uh+2*uh**2/th))+
     /     s2-th+4*uh*(q2/th+(1-uh/th)*(th-q2)/sh))+
     /     cf*(1/(sh*th)*((2*uh*(s2-uh)-s2**2)/sh-s2)*(fsu-2*ftu-flat)-
     /     2d0*(1/th-2/sh+4*uh/(sh*th))*ftu+(uh-3*s2+5*q2)/(sh*th)*flat+
     /     (5/sh+1/uh*(1+2*th/sh-3*s2/sh)-s2*q2/(sh*uh**2))*fsu+
     /     (8/sh+3/th+1/uh-4*uh/(sh*th)+(2*th-3*s2)/(sh*uh)-
     /     q2*(1/th**2+s2/(sh*uh**2)))*(fs2-fm2)-2/sh+3/(4*th)-1/uh+
     /     3*uh/(sh*th)-(th+s2)/(sh*uh)+q2*(s2/(sh*uh**2)-1/(2*th**2)))+
     /     ca*(2d0*(fs2-fm2)*(6/sh-2/th-2*q2/(sh*th)-3*s2/(sh*th))+
     /     (15/(2*sh)-9*s2/(2*sh*th)-2*q2/th**2)*fst+
     /     2d0*(1/th**2*(sh+s2**2/sh)+2*s2*q2/th**3-4*q2/(sh*th)*
     /     (1-q2/th)+2*s2*q2/(sh*th**2)*(1-s2/th)*(1-q2/th))*
     /     (fs2-fm2+fst)+1/sh*(4-2*uh/th-s2/th)*(fs2+fstu-(fst+fsu)/2)+
     /     1/sh*(1-uh/th)*fsu+1/sh*(2*uh/th-1)*ftu-2d0*(3/th+1/sh+
     /     2d0*(uh-s2)/(sh*th))*flat-1/sh-1/th+2*s2/(sh*th)+
     /     1/th**2*(3*q2+4*s2*(1-3*q2/th))-2/th**2*(sh+s2**2/sh)+
     /     2*q2/(sh*th)*(1-q2/th)-(12*s2*q2)/(sh*th**2)*(1-s2/th)*
     /     (1-q2/th))
      temp=temp1
      endif
      Cqg3=temp*xlumqg

      return
      end

c.....Third Real contribution to GQ->VQG (for nonzero values of s2)
      function Cgq3(sh,th,uh,q2,what)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton,what
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
c      s2=sh+th+uh-q2

      if (what.eq.1) then
c     term proportional to (fs2/s2)_A+
      temp1x=(-sh/uh-uh/sh-(2d0*th*q2)/(sh*uh))*(2d0*cf+4d0*ca)
      temp=temp1x
      elseif (what.eq.2) then
c     term proportional to (1/s2)_A+
      temp1x=(-sh/uh-uh/sh-(2d0*th*q2)/(sh*uh))*(
     /     cf*(-3/2d0+fst+2d0*ftu-2d0*fm2+(uh+th)*fla/la)+
     /     ca*(2d0*fstu+(fsu-fst)/2d0-ftu-2d0*fm2))
      temp=temp1x
      elseif (what.eq.3) then
c     term proportional to (1/s2), FINITE in the limit s2->0
      temp1x=(cf-ca/2d0)*((sh+2d0*th)/uh*((uh+th)*fla/la+fst)+
     /     (2d0*uh+4d0*th)*ftu/sh)
      temp=temp1x/s2
      elseif (what.eq.4) then
c     regular term
      temp1x=fla/la**3*(cf*((3*sh)*(sh+q2-s2)*(th**2-uh**2)*
     /     (1-th/uh)/(4*la**2)+(th**2-uh**2)*((q2-s2**2)/(4*sh)+
     /     11/4d0+s2/uh-th/(2*uh))-2*sh*(uh+th**2/uh)-
     /     sh/2d0*(sh-s2)*(th/uh-3d0)+4*s2*(uh-th)+5*sh*th)+
     /     ca*(1-th/uh)*((uh+th)/4d0*(sh+q2-s2)+sh*(sh-s2)))+
     /     fla/la*(cf*((uh-th-2*s2)/(4*sh)+(11*sh-2*th-4*s2)/(4*uh)+
     /     2d0/sh*(1+th/uh)*(5*q2+s2)-(16d0*q2*s2)/(sh*uh))-
     /     ca/2d0*((uh+th)/(2d0*uh)+5d0*(th/sh-(th+s2)/uh)-
     /     7d0/uh*(sh+th*s2/sh)+3d0/sh*(uh-3*s2)+
     /     2d0/(sh*uh)*(th**2+3*s2**2)))+
     /     1/la**2*(cf*(3*sh/(2*la**2)*(uh-th)**2*(1+th/uh)+
     /     sh/(2*uh)*(th-3*uh)+(1-th/uh)*(s2-7d0*th/4d0)+
     /     (uh-th)*(11/4d0+(3*uh+3*th)/(2*sh)-2*s2/sh))+
     /     ca*(th/uh-1)*(sh+s2))+
     /     2*du**3*sh*uh*(4*cf-3*ca-(cf-ca)*(fs2-fm2))+
     /     du**2*(cf*(3*sh-8*uh-4*th)-ca*(9*sh-4*uh+4*th)+
     /     2d0*(fs2-fm2)*(cf*(uh+th)+ca*(3*sh+2*th)))+
     /     du*((cf*(fs2-fm2)-ca*(fs2-fm2+(fstu+fs2-ftu)/2d0+
     /     (fsu+flau)/4d0))*(4d0*(1+th/sh)*(1-th/uh)-
     /     2d0*(sh/uh+uh/sh))-cf*(fs2-fm2)*(sh+2*th)/uh+
     /     ca*((2*flau+2*fsu)*(th+sh)/uh+2d0*(fs2-fm2)*
     /     ((5*sh+4*th)/uh-2))-cf*(4+(sh-2*th)/uh+(2*th-3*uh)/sh-
     /     th**2/(sh*uh))+ca*(7-(3*sh+2*th)/uh))+
     /     dt*(2*cf*(fm2-fs2-sh/uh)+ca*(sh/uh-1))+
     /     dsu*cf*(fst-2*ftu-flau)*(1-s2/uh)*(s2**2-2*th*(s2-th))/sh**2+
     /     dtu*cf*(2d0*(fm2-fs2-ftu)*(2*sh/uh*(q2+2*th)+uh-s2+
     /     4*th/uh*(2*th-s2)+2d0/sh*(th-s2)*(uh-2*th+2*th**2/uh))+
     /     s2-uh+4*th*(q2/uh+(1-th/uh)*(uh-q2)/sh))+
     /     cf*(1/(sh*uh)*((2*th*(s2-th)-s2**2)/sh-s2)*(fst-2*ftu-flau)-
     /     2d0*(1/uh-2/sh+4*th/(sh*uh))*ftu+(th-3*s2+5*q2)/(sh*uh)*flau+
     /     (5/sh+1/th*(1+2*uh/sh-3*s2/sh)-s2*q2/(sh*th**2))*fst+
     /     (8/sh+3/uh+1/th-4*th/(sh*uh)+(2*uh-3*s2)/(sh*th)-
     /     q2*(1/uh**2+s2/(sh*th**2)))*(fs2-fm2)-2/sh+3/(4*uh)-1/th+
     /     3*th/(sh*uh)-(uh+s2)/(sh*th)+q2*(s2/(sh*th**2)-1/(2*uh**2)))+
     /     ca*(2d0*(fs2-fm2)*(6/sh-2/uh-2*q2/(sh*uh)-3*s2/(sh*uh))+
     /     (15/(2*sh)-9*s2/(2*sh*uh)-2*q2/uh**2)*fsu+
     /     2d0*(1/uh**2*(sh+s2**2/sh)+2*s2*q2/uh**3-4*q2/(sh*uh)*
     /     (1-q2/uh)+2*s2*q2/(sh*uh**2)*(1-s2/uh)*(1-q2/uh))*
     /     (fs2-fm2+fsu)+1/sh*(4-2*th/uh-s2/uh)*(fs2+fstu-(fsu+fst)/2)+
     /     1/sh*(1-th/uh)*fst+1/sh*(2*th/uh-1)*ftu-2d0*(3/uh+1/sh+
     /     2d0*(th-s2)/(sh*uh))*flau-1/sh-1/uh+2*s2/(sh*uh)+
     /     1/uh**2*(3*q2+4*s2*(1-3*q2/uh))-2/uh**2*(sh+s2**2/sh)+
     /     2*q2/(sh*uh)*(1-q2/uh)-(12*s2*q2)/(sh*uh**2)*(1-s2/uh)*
     /     (1-q2/uh))
      temp=temp1x
      endif
      Cgq3=temp*xlumgq

      return
      end

c****************************************
c.....NLO real contributions to GG->VQQB
c****************************************

c.....First Real contribution to GG->VQQB 
      function Cgg1(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
c      s2=sh+th+uh-q2

      temp=ca*ds*fla/la*(th**2/la**2*(3/(2*la**2)*(2*sh-th-uh)*
     /     (uh**2-th**2)+2*sh-3*th+uh+(th**2-uh**2)/sh)-3*sh/2-
     /     th*(9/2d0+5*th/sh+3*uh/sh))+
     /     ca*fla/la*(th**2/la**2*(3d0/(2*la**2)*(uh**2-th**2)+1)*
     /     (q2-s2)/sh-11d0/(4*sh)*(q2-s2)-4*(1+th/sh))+
     /     cf*fla/la*(2*ds*(th+uh-2*sh)+6)+ca*th/la**2*
     /     (3*(th-uh)**2/la**2*(1+(th+uh)/(2*sh)-2*sh*ds)+
     /     ds*(2*sh-3*th+3*uh)-1-uh/sh+s2*(th-uh)/sh**2+
     /     (ds-1/(2*sh))*(th-uh)**2/sh)+
     /     4*cf*dtu*((fs2-fm2+ftu)*(sh/2+2*th*(1+th/sh))-
     /     th*(1+(th+uh)/sh))+
     /     2*ds*((ca-2*cf)*dt*th**2/sh*(fst+flat)+cf*dst*
     /     (fsu-flat-2*ftu)*(sh+2*uh*(1+uh/sh)))+
     /     ds*(2*cf*(fst-(3+4*(th+uh)/sh)*flat-(4+8*th/sh)*ftu)+
     /     ca/2d0*((fst+flat)*(2+(uh+3*th)/sh)+1+th*(th-uh)/sh**2))+
     /     2*cf*dst*(2*ftu+flat-fsu)*(2+(th+uh)/sh)+8*th*dt**2*
     /     (cf*(1-fs2+fm2)+ca)+
     /     (4*cf-2*ca)*dt*du*th*uh/sh*(fs2-ftu+fstu)+
     /     dt*((ca-2*cf)*((2+(th+uh)/sh)*(fst+flat)+2/sh*(uh-th)*
     /     (fs2+fstu-ftu))-8*cf*(fs2-fm2-1)+4*ca*(2+(uh-th)/sh))+
     /     cf*(2*(fs2-fm2)*(2/sh-2/th-1/(sh*th)*(uh+s2*q2/th))+
     /     4/sh*(fs2+fstu)-2/th*(2+1/sh*(uh+s2*q2/th))*fst+
     /     1/(sh*th)*(4*q2*(1+s2/th)-2*uh))-
     /     ca/sh*(2*(fs2+fstu)+fst/2+5*flat/2d0+15/4d0)
      Cgg1=temp*xlumgg

      return
      end

c.....Crossed Real contribution to GG->VQQB
      function Cgg1x(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
c      s2=sh+th+uh-q2

      temp=ca*ds*fla/la*(uh**2/la**2*(3/(2*la**2)*(2*sh-uh-th)*
     /     (th**2-uh**2)+2*sh-3*uh+th+(uh**2-th**2)/sh)-3*sh/2-
     /     uh*(9/2d0+5*uh/sh+3*th/sh))+
     /     ca*fla/la*(uh**2/la**2*(3d0/(2*la**2)*(th**2-uh**2)+1)*
     /     (q2-s2)/sh-11d0/(4*sh)*(q2-s2)-4*(1+uh/sh))+
     /     cf*fla/la*(2*ds*(uh+th-2*sh)+6)+ca*uh/la**2*
     /     (3*(uh-th)**2/la**2*(1+(uh+th)/(2*sh)-2*sh*ds)+
     /     ds*(2*sh-3*uh+3*th)-1-th/sh+s2*(uh-th)/sh**2+
     /     (ds-1/(2*sh))*(uh-th)**2/sh)+
     /     4*cf*dtu*((fs2-fm2+ftu)*(sh/2+2*uh*(1+uh/sh))-
     /     uh*(1+(uh+th)/sh))+
     /     2*ds*((ca-2*cf)*du*uh**2/sh*(fsu+flau)+cf*dsu*
     /     (fst-flau-2*ftu)*(sh+2*th*(1+th/sh)))+
     /     ds*(2*cf*(fsu-(3+4*(uh+th)/sh)*flau-(4+8*uh/sh)*ftu)+
     /     ca/2d0*((fsu+flau)*(2+(th+3*uh)/sh)+1+uh*(uh-th)/sh**2))+
     /     2*cf*dsu*(2*ftu+flau-fst)*(2+(uh+th)/sh)+8*uh*du**2*
     /     (cf*(1-fs2+fm2)+ca)+
     /     (4*cf-2*ca)*du*dt*uh*th/sh*(fs2-ftu+fstu)+
     /     du*((ca-2*cf)*((2+(uh+th)/sh)*(fsu+flau)+2/sh*(th-uh)*
     /     (fs2+fstu-ftu))-8*cf*(fs2-fm2-1)+4*ca*(2+(th-uh)/sh))+
     /     cf*(2*(fs2-fm2)*(2/sh-2/uh-1/(sh*uh)*(th+s2*q2/uh))+
     /     4/sh*(fs2+fstu)-2/uh*(2+1/sh*(th+s2*q2/uh))*fsu+
     /     1/(sh*uh)*(4*q2*(1+s2/uh)-2*th))-
     /     ca/sh*(2*(fs2+fstu)+fsu/2+5*flau/2d0+15/4d0)
      Cgg1x=temp*xlumgg

      return
      end

c****************************************
c.....NLO real contributions to QQB->VGG
c****************************************

c.....First Real contribution to QQB->VGG (proportional to delta(s2))
      function Cqqb1(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
      
      temp=(uh/th+th/uh+(2d0*sh*q2)/(uh*th))*
     /     (ca*(67/18d0-11*fa/6d0+fa**2)+cf*(2*ft+2*fu-4*fa-3)*fm2+
     /     (cf-ca/2d0)*(pi**2/3d0+(2*fa+fs-ft-fu)**2))
      Cqqb1=temp*xlumqqb

      return
      end

c.....Second Real contribution to QQB->VGG (for nonzero values of s2)
      function Cqqb2(sh,th,uh,q2,what)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton,what
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
c      s2=sh+th+uh-q2
      
      if (what.eq.1) then
c     term proportional to (fs2/s2)_A+
      temp1=1/2d0*(uh/th+th/uh+2*q2*sh/(th*uh))*(8*cf-2*ca)
      temp=temp1
      elseif (what.eq.2) then
c     term proportional to (1/s2)_A+
      temp1=1/2d0*(uh/th+th/uh+2*q2*sh/(th*uh))*
     /     (-11*ca/6d0+2*cf*(ftu-2*fm2)+(2*cf-ca)*(2*fstu-ftu))
      temp=temp1
      elseif (what.eq.4) then
c     regular term
      temp1=cf*(sh*dt**2+2*dt-sh/(th*uh)+dtu*(s2+2*q2*(uh-s2)/th)-
     /     2*dt*q2*(1/uh-1/th))+
     /     ca*(-11*sh/(6*th*uh)+dt**2*sh**2/uh*(3*s2/(2*th)-2)+
     /     2*dt*sh/uh+q2/(3*th**2))+
     /     (2*cf-ca)*(dt/uh*(q2-uh)**2*(du-1/th)*(fstu+fs2-ftu)+
     /     (sh-q2)/(th*uh)*(fstu+fs2))+
     /     ftu*(cf*dtu*(4*q2/(th*uh)*(q2-th)**2+2*(q2+sh)-s2)+
     /     cf*(s2-2*sh)/(th*uh)-ca*q2/(th*uh))+
     /     cf*(fm2-fs2)*(dtu*(4*(uh-q2)-s2-4*q2/(th*uh)*(uh-q2)**2)+
     /     dt**2*(q2-uh)-dt/th*(2*q2-uh)-2/th+q2/th**2+4*q2/(th*uh))
      temp=temp1
      endif
      Cqqb2=temp*xlumqqb
      if (Cqqb2.ne.Cqqb2) print *,'Cqqb2',what,s2,fs2
      return
      end

c.....Crossed Second Real contribution to QQB->VGG (for nonzero values of s2)
      function Cqqb2x(sh,th,uh,q2,what)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton,what
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
c      s2=sh+th+uh-q2

      if (what.eq.1) then
c     term proportional to (fs2/s2)_A+
      temp1x=1/2d0*(th/uh+uh/th+2*q2*sh/(uh*th))*(8*cf-2*ca)
      temp=temp1x
      elseif (what.eq.2) then
c     term proportional to (1/s2)_A+
      temp1x=1/2d0*(th/uh+uh/th+2*q2*sh/(uh*th))*
     /     (-11*ca/6+2*cf*(ftu-2*fm2)+(2*cf-ca)*(2*fstu-ftu))
      temp=temp1x
      elseif (what.eq.4) then
c     regular term
      temp1x=cf*(sh*du**2+2*du-sh/(uh*th)+dtu*(s2+2*q2*(th-s2)/uh)-
     /     2*du*q2*(1/th-1/uh))+
     /     ca*(-11*sh/(6*uh*th)+du**2*sh**2/th*(3*s2/(2*uh)-2)+
     /     2*du*sh/th+q2/(3*uh**2))+
     /     (2*cf-ca)*(du/th*(q2-th)**2*(dt-1/uh)*(fstu+fs2-ftu)+
     /     (sh-q2)/(uh*th)*(fstu+fs2))+
     /     ftu*(cf*dtu*(4*q2/(uh*th)*(q2-uh)**2+2*(q2+sh)-s2)+
     /     cf*(s2-2*sh)/(uh*th)-ca*q2/(uh*th))+
     /     cf*(fm2-fs2)*(dtu*(4*(th-q2)-s2-4*q2/(uh*th)*(th-q2)**2)+
     /     du**2*(q2-th)-du/uh*(2*q2-th)-2/uh+q2/uh**2+4*q2/(uh*th))
      temp=temp1x
      endif
      Cqqb2x=temp*xlumqqb

      return
      end

c********************************************************
c.....NLO real contributions for QQB->VQQB (=QBQ->VQQB)
c********************************************************

c.....First Real contribution to QQB->VQQB
      function D0aa(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
      
      temp=(uh/th+th/uh+(2d0*sh*q2)/(uh*th))*(fa/3d0-5/9d0)*nf
      D0aa=temp*xlumqqb

      return
      end

c.....Second Real contribution to QQB->VQQB
      function Daa(sh,th,uh,q2,what)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton,what
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
      
      if (what.eq.2) then
c     term proportional to (1/s2)_A+
      temp1=1/6d0*(uh/th+th/uh+(2*sh*q2)/(th*uh))
      temp=temp1*nf
      elseif (what.eq.4) then
c     regular term
      temp1=1/3d0*(sh/(th*uh)-q2/th**2)
      temp=temp1*nf
      endif
      Daa=temp*xlumqqb

      return
      end

c.....Second Real contribution to QQB->VQQB (crossed)
      function Daax(sh,th,uh,q2,what)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton,what
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
      
      if (what.eq.2) then
c     term proportional to (1/s2)_A+
      temp1=1/6d0*(th/uh+uh/th+(2*sh*q2)/(uh*th))
      temp=temp1*nf
      elseif (what.eq.4) then
c     regular term
      temp1=1/3d0*(sh/(uh*th)-q2/uh**2)
      temp=temp1*nf
      endif
      Daax=temp*xlumqqb

      return
      end

c.....Third real (only for Z production) contribution to QQB->VQQB
      function Dab(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
c      s2=sh+th+uh-q2

      temp=1/2d0*(fla*q2/(la*th)*((2*sh+th+uh)/la**2*(3*sh*(th-uh)**2/
     /     la**2+uh-th-sh)-1)+(2*sh+th+uh)/(2*la**2*th)*(-3/la**2*
     /     (sh+q2-s2)*(th-uh)**2+3*th-uh+2*sh-2*s2))
      Dab=temp*xlumqqbtr

      return
      end
      
c.....Third real (only for Z production) contribution to QQB->VQQB (crossed)
      function Dabx(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
c      s2=sh+th+uh-q2

      temp=1/2d0*(fla*q2/(la*uh)*((2*sh+uh+th)/la**2*(3*sh*(uh-th)**2/
     /     la**2+th-uh-sh)-1)+(2*sh+uh+th)/(2*la**2*uh)*(-3/la**2*
     /     (sh+q2-s2)*(uh-th)**2+3*uh-th+2*sh-2*s2))
      Dabx=temp*xlumqqbtr

      return
      end
      
c.....Fourth real contribution to QQB->VQQB
      function Dbb(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
c      s2=sh+th+uh-q2

      temp=1/2d0*(ds/2d0*(uh/sh-1d0)-5d0/(4*sh)+1/la**2*(uh*ds*(-2*sh+
     /     3*uh/(2*sh)*(th-uh)+4*uh-2*th)+uh/(2*sh)*(2*sh+2*s2+th-uh))+
     /     3*uh**2/la**4*(uh-th)*(ds*(2*sh-uh-th)-(sh+s2)/sh)+fla/la*
     /     (ds*(2*uh**2/sh+5*uh/2d0+3*sh/2d0)+3/4d0+uh/sh-s2/(2*sh))+
     /     fla/la**3*(uh**2*ds*(3*uh-th-uh**2/sh+th**2/sh-2*sh)+uh/sh*
     /     (2*s2*th-uh*th-2*s2**2+4*uh*s2-3*uh**2+2*sh*s2-uh*sh))+
     /     fla/la**5*(3*ds*uh**2*(uh**2-th**2)*(sh-q2)+3*uh**2*q2/sh*
     /     (uh-th)*(uh+th-2*s2)))
      Dbb=temp*xlumqqbdbb

      return
      end

c.....Fourth real contribution to QQB->VQQB (crossed)
      function Dbbx(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
c      s2=sh+th+uh-q2

      temp=1/2d0*(ds/2d0*(th/sh-1d0)-5d0/(4*sh)+1/la**2*(th*ds*(-2*sh+
     /     3*th/(2*sh)*(uh-th)+4*th-2*uh)+th/(2*sh)*(2*sh+2*s2+uh-th))+
     /     3*th**2/la**4*(th-uh)*(ds*(2*sh-th-uh)-(sh+s2)/sh)+fla/la*
     /     (ds*(2*th**2/sh+5*th/2d0+3*sh/2d0)+3/4d0+th/sh-s2/(2*sh))+
     /     fla/la**3*(th**2*ds*(3*th-uh-th**2/sh+uh**2/sh-2*sh)+th/sh*
     /     (2*s2*uh-th*uh-2*s2**2+4*th*s2-3*th**2+2*sh*s2-th*sh))+
     /     fla/la**5*(3*ds*th**2*(th**2-uh**2)*(sh-q2)+3*th**2*q2/sh*
     /     (th-uh)*(th+uh-2*s2)))
      Dbbx=temp*xlumqqbdbb

      return
      end

c.....Fifth real contribution to QQB->VQQB
      function Dac(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
c      s2=sh+th+uh-q2

      temp=(cf-ca/2d0)*(fla/la*(3*sh**2*q2/la**4*(th-uh)**2*(1/th+1/uh)+
     /     sh*q2/la**2*(5*sh/th-7*sh/uh+4*uh/th-4)+
     /     1/s2*(2*sh**2/uh+(q2*(2*sh+uh)-uh*sh)/th)+(2*s2-uh)/th-2)+
     /     (3*sh*q2)/(la**4*th*uh)*(2*s2-uh-th)*(th-uh)**2+
     /     1/la**2*(sh*(sh-s2)*(7/uh-5/th)+1/2d0*(th/uh-uh/th)*
     /     (10*sh+th+3*uh-2*s2)+2*sh-2*s2*(1-uh/th))+
     /     dt*sh**2/uh*(dt+3/th)+fst/th*(1/s2*(2*sh**2/uh+2*sh+uh)-
     /     2*q2/th)-1/(2*uh)+3/(2*th)-3*q2/th**2)
      Dac=temp*xlumqqb

      return
      end

c.....Sixth real contribution to QQB->VQQB
      function Dbc(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
c      s2=sh+th+uh-q2

      temp=(cf-ca/2d0)*(fla/la*(6*sh*s2*q2/(la**4*th)*(th-uh)**2+
     /     1/la**2*(2*sh*uh+(s2*(1+th/sh)-th/(2*sh)*(th+uh))*(th+uh)*
     /     (1-uh/th)+3*sh/2d0*(th-uh)*(1-uh/th)-4*s2*q2*(1+(sh-uh)/th))+
     /     ds*(6*sh+2*uh+(th**2-uh**2)/(2*sh))-1+(sh+4*uh-2*s2)/(2*th)+
     /     (3*th+2*uh-6*s2)/sh+(uh**2-2*uh*s2+4*s2**2)/(sh*th))+
     /     3*s2/la**4*(th-uh)**2*(1+uh/th-2*q2/th)+1/la**2*(-th/2+3*uh+
     /     s2+(th**2-uh**2)/sh-uh/(2*th)*(uh+6*s2))+2*dt+ds*(fst+flat)*
     /     (2*(sh+uh)/th+1/(2*sh)*(th+uh**2/th))+fst/2*(2/th+3/sh+
     /     (uh-2*s2)/(sh*th))+flat/2*(-2/th+1/sh-(uh+2*s2)/(sh*th))+
     /     1/(2*th)-1/sh+1/(sh*th)*(2*uh-4*s2*q2/th))
      Dbc=temp*xlumqqbdbc

      return
      end

c.....Seventh real contribution to QQB->VQQB (Dac crossed)
      function Dad(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
c      s2=sh+th+uh-q2

      temp=(cf-ca/2d0)*(fla/la*(3*sh**2*q2/la**4*(uh-th)**2*(1/uh+1/th)+
     /     sh*q2/la**2*(5*sh/uh-7*sh/th+4*th/uh-4)+
     /     1/s2*(2*sh**2/th+(q2*(2*sh+th)-th*sh)/uh)+(2*s2-th)/uh-2)+
     /     (3*sh*q2)/(la**4*uh*th)*(2*s2-th-uh)*(uh-th)**2+
     /     1/la**2*(sh*(sh-s2)*(7/th-5/uh)+1/2d0*(uh/th-th/uh)*
     /     (10*sh+uh+3*th-2*s2)+2*sh-2*s2*(1-th/uh))+
     /     du*sh**2/th*(du+3/uh)+fsu/uh*(1/s2*(2*sh**2/th+2*sh+th)-
     /     2*q2/uh)-1/(2*th)+3/(2*uh)-3*q2/uh**2)
      Dad=temp*xlumqqb

      return
      end

c.....Eigth real contribution to QQB->VQQB (Dbc crossed)
      function Dbd(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
      real *8 fh1(-5:5),fh2(-5:5)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
c      s2=sh+th+uh-q2

      temp=(cf-ca/2d0)*(fla/la*(6*sh*s2*q2/(la**4*uh)*(uh-th)**2+
     /     1/la**2*(2*sh*th+(s2*(1+uh/sh)-uh/(2*sh)*(uh+th))*(uh+th)*
     /     (1-th/uh)+3*sh/2d0*(uh-th)*(1-th/uh)-4*s2*q2*(1+(sh-th)/uh))+
     /     ds*(6*sh+2*th+(uh**2-th**2)/(2*sh))-1+(sh+4*th-2*s2)/(2*uh)+
     /     (3*uh+2*th-6*s2)/sh+(th**2-2*th*s2+4*s2**2)/(sh*uh))+
     /     3*s2/la**4*(uh-th)**2*(1+th/uh-2*q2/uh)+1/la**2*(-uh/2+3*th+
     /     s2+(uh**2-th**2)/sh-th/(2*uh)*(th+6*s2))+2*du+ds*(fsu+flau)*
     /     (2*(sh+th)/uh+1/(2*sh)*(uh+th**2/uh))+fsu/2*(2/uh+3/sh+
     /     (th-2*s2)/(sh*uh))+flau/2*(-2/uh+1/sh-(th+2*s2)/(sh*uh))+
     /     1/(2*uh)-1/sh+1/(sh*uh)*(2*th-4*s2*q2/uh))
      Dbd=temp*xlumqqbdbc

      return
      end

c.....Nineth real contribution to QQB->VQQB
      function Dcc(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
c      s2=sh+th+uh-q2

      temp=1/2d0*(fla/la*(sh/(la**2*th)*(uh**2-th**2)+(3*sh+2*uh)/th)+
     /     (th+uh-2*s2)/la**2*(1-uh/th)+dt*(fs2-fm2)*(sh*dt+4*sh/th+
     /     2*uh/th)+(fst-fm2+fs2)*((sh+s2-q2)**2/(sh*th**2)+1/sh*
     /     (uh/th-2*s2*q2/th**2)**2+4/th*(q2/th-1))+
     /     2*fst/th*(1-q2/th)-dt*(sh*dt-1-uh/th)+2/th*(uh/th+s2*(1/th+
     /     1/sh-s2/(sh*th)))+1/th*(q2/th-1)+12*s2*q2/(sh*th**3)*
     /     (uh-s2*q2/th))
      Dcc=temp*xlumqqbdcc

      return
      end

c.....Tenth real contribution to QQB->VQQB (Dcc crossed)
      function Ddd(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
       
c      s2=sh+th+uh-q2

      temp=1/2d0*(fla/la*(sh/(la**2*uh)*(th**2-uh**2)+(3*sh+2*th)/uh)+
     /     (uh+th-2*s2)/la**2*(1-th/uh)+du*(fs2-fm2)*(sh*du+4*sh/uh+
     /     2*th/uh)+(fsu-fm2+fs2)*((sh+s2-q2)**2/(sh*uh**2)+1/sh*
     /     (th/uh-2*s2*q2/uh**2)**2+4/uh*(q2/uh-1))+
     /     2*fsu/uh*(1-q2/uh)-du*(sh*du-1-th/uh)+2/uh*(th/uh+s2*(1/uh+
     /     1/sh-s2/(sh*uh)))+1/uh*(q2/uh-1)+12*s2*q2/(sh*uh**3)*
     /     (th-s2*q2/uh))
c      Ddd=temp*xlumqqbdcc
      Ddd=temp*xlumqqbddd !!! -> Bug Fix, (see Formula 2.18 of http://journals.aps.org/prd/pdf/10.1103/PhysRevD.40.2245)

      return
      end

c.....Eleventh real contribution to QQB->VQQB (cd interference/same helicity)
      function DcdLL(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
c      s2=sh+th+uh-q2

      temp=1/2d0*(fla/la*(-sh/(la**2*th)*(th-uh)**2+2*ds*(sh-th)+
     /     (3*sh+2*s2)/th-1)+(uh-2*s2)/la**2*(uh/th-1)+ds*dst/th*
     /     (fsu-flat-2*ftu)*(2*sh*(sh+uh)+uh**2)+
     /     ds*((fst-ftu)*(1+uh/th+2*th/uh+4*sh*(q2+s2)/(th*uh))-
     /     (flat+ftu)*(1+uh/th))+
     /     dst*(2*ftu-fsu+flat)*(1+(2*sh+uh)/th)+2*dt*sh/uh+ftu/(th*uh)*
     /     (sh+s2)-1/th*(fst+flat+1-sh/uh))
      DcdLL=temp*xlumqqbLL

      return
      end
c.....Eleventh real contribution to QQB->VQQB (crossed)
      function DcdLLx(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
c      s2=sh+th+uh-q2

      temp=1/2d0*(fla/la*(-sh/(la**2*uh)*(uh-th)**2+2*ds*(sh-uh)+
     /     (3*sh+2*s2)/uh-1)+(th-2*s2)/la**2*(th/uh-1)+ds*dsu/uh*
     /     (fst-flau-2*ftu)*(2*sh*(sh+th)+th**2)+
     /     ds*((fsu-ftu)*(1+th/uh+2*uh/th+4*sh*(q2+s2)/(uh*th))-
     /     (flau+ftu)*(1+th/uh))+
     /     dsu*(2*ftu-fst+flau)*(1+(2*sh+th)/uh)+2*du*sh/th+ftu/(uh*th)*
     /     (sh+s2)-1/uh*(fsu+flau+1-sh/th))
      DcdLLx=temp*xlumqqbLL

      return
      end

c.....Twelth real contribution to QQB->VQQB (cd interference/different helicity)
      function DcdLR(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
c      s2=sh+th+uh-q2

      temp=1/2d0*(2*ds*sh**2/th*(dst*(fsu-flat-2*ftu)+2/uh*(fst-ftu))+
     /     dst*(2*ftu-fsu+flat)*(1+2*sh/th-uh/th)+2d0/(th*uh)*
     /     (ftu-fsu)*(sh+s2-uh)-2d0/th*(ftu+flat))
      DcdLR=temp*xlumqqbLR

      return
      end

c.....Twelth real contribution to QQB->VQQB (crossed)
      function DcdLRx(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
c      s2=sh+th+uh-q2

      temp=1/2d0*(2*ds*sh**2/uh*(dsu*(fst-flau-2*ftu)+2/th*(fsu-ftu))+
     /     dsu*(2*ftu-fst+flau)*(1+2*sh/uh-th/uh)+2d0/(uh*th)*
     /     (ftu-fst)*(sh+s2-th)-2d0/uh*(ftu+flau))
      DcdLRx=temp*xlumqqbLR

      return
      end

c********************************************************
c.....NLO real contributions for QQ->VQQ (=QBQB->VQQ)
c********************************************************

c.....First contribution which is different from D functions
      function Eac(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
      temp=(cf-ca/2d0)*((fst+flat)*(4*dt*ds*sh**2/th-2*dt*(1-uh/th)-
     /     4/th)+4*q2/th**2*fst)
      Eac=temp*xlumqq

      return
      end

c.....Second contribution which is different from D functions
      function Ead(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
c      s2=sh+th+uh-q2

      temp=(cf-ca/2d0)*(1/sh*(2*fs2-fst-fsu+2*fstu)*(-2-uh/th-th/uh+
     /     2*s2/(th*uh)*(q2-sh))-2/sh*(2+uh/th+th/uh-2*s2*q2*
     /     (1/th**2+1/uh**2)))
      Ead=temp*xlumqqead

      return
      end

c.....Third contribution which is different from D functions (Eac crossed)
      function Ebd(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
      temp=(cf-ca/2d0)*((fsu+flau)*(4*du*ds*sh**2/uh-2*du*(1-th/uh)-
     /     4/uh)+4*q2/uh**2*fsu)
      Ebd=temp*xlumqq

      return
      end

c.....Fourth contribution which is different from D functions (Ead crossed)
      function Ebc(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
c      s2=sh+th+uh-q2

      temp=(cf-ca/2d0)*(1/sh*(2*fs2-fsu-fst+2*fstu)*(-2-th/uh-uh/th+
     /     2*s2/(uh*th)*(q2-sh))-2/sh*(2+th/uh+uh/th-2*s2*q2*
     /     (1/uh**2+1/th**2)))
      Ebc=temp*xlumqqead

      return
      end


!!! MIO
c.....Nineth real contribution to QQB->VQQB
      function Eaa(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
c      s2=sh+th+uh-q2

      temp=1/2d0*(fla/la*(sh/(la**2*th)*(uh**2-th**2)+(3*sh+2*uh)/th)+
     /     (th+uh-2*s2)/la**2*(1-uh/th)+dt*(fs2-fm2)*(sh*dt+4*sh/th+
     /     2*uh/th)+(fst-fm2+fs2)*((sh+s2-q2)**2/(sh*th**2)+1/sh*
     /     (uh/th-2*s2*q2/th**2)**2+4/th*(q2/th-1))+
     /     2*fst/th*(1-q2/th)-dt*(sh*dt-1-uh/th)+2/th*(uh/th+s2*(1/th+
     /     1/sh-s2/(sh*th)))+1/th*(q2/th-1)+12*s2*q2/(sh*th**3)*
     /     (uh-s2*q2/th))
      Eaa=temp*xlumqqeaa

      return
      end

c.....Tenth real contribution to QQB->VQQB (Dcc crossed)
      function Ebb(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
       
c      s2=sh+th+uh-q2

      temp=1/2d0*(fla/la*(sh/(la**2*uh)*(th**2-uh**2)+(3*sh+2*th)/uh)+
     /     (uh+th-2*s2)/la**2*(1-th/uh)+du*(fs2-fm2)*(sh*du+4*sh/uh+
     /     2*th/uh)+(fsu-fm2+fs2)*((sh+s2-q2)**2/(sh*uh**2)+1/sh*
     /     (th/uh-2*s2*q2/uh**2)**2+4/uh*(q2/uh-1))+
     /     2*fsu/uh*(1-q2/uh)-du*(sh*du-1-th/uh)+2/uh*(th/uh+s2*(1/uh+
     /     1/sh-s2/(sh*uh)))+1/uh*(q2/uh-1)+12*s2*q2/(sh*uh**3)*
     /     (th-s2*q2/uh))
      Ebb=temp*xlumqqeaa
c      Ebb=temp*xlumqqebb !!! -> Bug Fix, but this piece is not used (see Formula 2.21 of http://journals.aps.org/prd/pdf/10.1103/PhysRevD.40.2245)

      return
      end
c.....Eleventh real contribution to QQB->VQQB (cd interference/same helicity)
      function EabLR(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
c      s2=sh+th+uh-q2

      temp=-(1/2d0*(fla/la*(-sh/(la**2*th)*(th-uh)**2+2*ds*(sh-th)+
     /     (3*sh+2*s2)/th-1)+(uh-2*s2)/la**2*(uh/th-1)+ds*dst/th*
     /     (fsu-flat-2*ftu)*(2*sh*(sh+uh)+uh**2)+
     /     ds*((fst-ftu)*(1+uh/th+2*th/uh+4*sh*(q2+s2)/(th*uh))-
     /     (flat+ftu)*(1+uh/th))+
     /     dst*(2*ftu-fsu+flat)*(1+(2*sh+uh)/th)+2*dt*sh/uh+ftu/(th*uh)*
     /     (sh+s2)-1/th*(fst+flat+1-sh/uh)))
      EabLR=temp*xlumqqLR

      return
      end
c.....Eleventh real contribution to QQB->VQQB (crossed)
      function EabLRx(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
c      s2=sh+th+uh-q2

      temp=-(1/2d0*(fla/la*(-sh/(la**2*uh)*(uh-th)**2+2*ds*(sh-uh)+
     /     (3*sh+2*s2)/uh-1)+(th-2*s2)/la**2*(th/uh-1)+ds*dsu/uh*
     /     (fst-flau-2*ftu)*(2*sh*(sh+th)+th**2)+
     /     ds*((fsu-ftu)*(1+th/uh+2*uh/th+4*sh*(q2+s2)/(uh*th))-
     /     (flau+ftu)*(1+th/uh))+
     /     dsu*(2*ftu-fst+flau)*(1+(2*sh+th)/uh)+2*du*sh/th+ftu/(uh*th)*
     /     (sh+s2)-1/uh*(fsu+flau+1-sh/th)))
      EabLRx=temp*xlumqqLR

      return
      end

c.....Twelth real contribution to QQB->VQQB (cd interference/different helicity)
      function EabLL(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
c      s2=sh+th+uh-q2

      temp=-(1/2d0*(2*ds*sh**2/th*(dst*(fsu-flat-2*ftu)+2/uh*(fst-ftu))+
     /     dst*(2*ftu-fsu+flat)*(1+2*sh/th-uh/th)+2d0/(th*uh)*
     /     (ftu-fsu)*(sh+s2-uh)-2d0/th*(ftu+flat)))
      EabLL=temp*xlumqqLL

      return
      end

c.....Twelth real contribution to QQB->VQQB (crossed)
      function EabLLx(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
c      s2=sh+th+uh-q2

      temp=-(1/2d0*(2*ds*sh**2/uh*(dsu*(fst-flau-2*ftu)+2/th*(fsu-ftu))+
     /     dsu*(2*ftu-fst+flau)*(1+2*sh/uh-th/uh)+2d0/(uh*th)*
     /     (ftu-fst)*(sh+s2-th)-2d0/uh*(ftu+flau)))
      EabLLx=temp*xlumqqLL

      return
      end

!!! FINE MIO

c     Need to add missing contributions of Eq. (2.21) of [Gonsalves, Pawlowsky, Wai]
      function Ecc(sh,th,uh,q2)
      implicit real *8 (a-h,o-z)
c      real *8 la
      integer prodflag,nf,ih1,ih2,isetproton
      common/nf/nf
      common/const2/pi,cf,ca,tr,xnc
      include 'functions_inc.f'
      include 'luminosities_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
c     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
c     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
c     /     xlumqqLL,xlumqqLR
c      common/recmass/s2
      
c      s2=sh+th+uh-q2

      temp=1/2d0*(fla/la*(sh/(la**2*th)*(uh**2-th**2)+(3*sh+2*uh)/th)+
     /     (th+uh-2*s2)/la**2*(1-uh/th)+dt*(fs2-fm2)*(sh*dt+4*sh/th+
     /     2*uh/th)+(fst-fm2+fs2)*((sh+s2-q2)**2/(sh*th**2)+1/sh*
     /     (uh/th-2*s2*q2/th**2)**2+4/th*(q2/th-1))+
     /     2*fst/th*(1-q2/th)-dt*(sh*dt-1-uh/th)+2/th*(uh/th+s2*(1/th+
     /     1/sh-s2/(sh*th)))+1/th*(q2/th-1)+12*s2*q2/(sh*th**3)*
     /     (uh-s2*q2/th))
      Ecc=temp*xlumqqecc

      return
      end
c end missing      
      
c********************************
c.....Dilogarithmic function Li2 
c********************************
     
c.....Function Li2(x) for arguments x < = 1.0                             
      function dyLI2(x)                                                     
      implicit none                                                           
      real*8 X,Y,T,S,A,PI3,PI6,ZERO,ONE,HALF,MALF,MONE,MTWO  
      real*8 C(0:18),H,ALFA,B0,B1,B2,LI2OLD                 
      real*8 dyLi2                                             
      integer i                                                       
      
      data ZERO /0.0d0/, ONE /1.0d0/                               
      data HALF /0.5d0/, MALF /-0.5d0/                             
      data MONE /-1.0d0/, MTWO /-2.0d0/                            
      data PI3 /3.289868133696453d0/, PI6 /1.644934066848226d0/               
      data C( 0) / 0.4299669356081370d0/
      data C( 1) / 0.4097598753307711d0/                              
      data C( 2) /-0.0185884366501460d0/                              
      data C( 3) / 0.0014575108406227d0/                              
      data C( 4) /-0.0001430418444234d0/                              
      data C( 5) / 0.0000158841554188d0/                              
      data C( 6) /-0.0000019078495939d0/                              
      data C( 7) / 0.0000002419518085d0/                              
      data C( 8) /-0.0000000319334127d0/                              
      data C( 9) / 0.0000000043454506d0/                              
      data C(10) /-0.0000000006057848d0/                              
      data C(11) / 0.0000000000861210d0/                              
      data C(12) /-0.0000000000124433d0/                              
      data C(13) / 0.0000000000018226d0/                              
      data C(14) /-0.0000000000002701d0/                              
      data C(15) / 0.0000000000000404d0/                              
      data C(16) /-0.0000000000000061d0/                              
      data C(17) / 0.0000000000000009d0/                              
      data C(18) /-0.0000000000000001d0/                              
      
      if(x .gt. 1.00000000001d0) then                                    
         write(6,*)'problems in dyLI2'
         write(6,*)'x=',x 
         stop                                               
      elseif(x .gt. 1.0d0) then                                          
         x = 1.d0                                                      
      endif                                                              
      if(X .eq. ONE) then                                                
         LI2OLD=PI6
         dyLI2=LI2OLD                                                       
         return                                                            
      else if(X .eq. MONE) then                                          
         LI2OLD=MALF*PI6
         dyLI2=LI2OLD                                                  
         return                                                            
      end if                                                             
      T=-X                                                               
      if(T .le. MTWO) then                                               
         Y=MONE/(ONE+T)                                                    
         S=ONE                                                             
         A=-PI3+HALF*(LOG(-T)**2-LOG(ONE+ONE/T)**2)                        
      else if(T .lt. MONE) then                                          
         Y=MONE-T                                                          
         S=MONE                                                            
         A=LOG(-T)                                                         
         A=-PI6+A*(A+LOG(ONE+ONE/T))                                       
      else if(T .le. MALF) then                                          
         Y=(MONE-T)/T                                                      
         S=ONE                                                             
         A=LOG(-T)                                                         
         A=-PI6+A*(MALF*A+LOG(ONE+T))                                      
      else if(T .lt. ZERO) then                                          
         Y=-T/(ONE+T)                                                      
         S=MONE                                                            
         A=HALF*LOG(ONE+T)**2                                              
      else if(T .le. ONE) then                                           
         Y=T                                                               
         S=ONE                                                             
         A=ZERO                                                            
      else                                                               
         Y=ONE/T                                                           
         S=MONE                                                            
         A=PI6+HALF*LOG(T)**2                                              
      end if                                                             
      
      H=Y+Y-ONE                                                          
      ALFA=H+H                                                           
      B1=ZERO                                                            
      B2=ZERO                                                            
      do I = 18,0,-1                                                    
         B0=C(I)+ALFA*B1-B2                                               
         B2=B1                                                            
         B1=B0                                                            
      enddo                                                              
      LI2OLD=-(S*(B0-H*B2)+A) 
      dyLI2=LI2OLD
      end     
      
c*************************************************************
c.....Frequently occurring functions used in defining B,C,D,E
c*************************************************************
      subroutine utilities3(q2)
      implicit none
      real *8 q2,a
      include 'scales2_inc.f'
      include 'functions_inc.f'
      fm2=log(xmuf2/q2)
      fmu2=log(xmur2/q2)
      end

      subroutine utilities2(uh,q2)
      implicit none
      real *8 uh,q2,a
      include 'scales2_inc.f'
      include 'functions_inc.f'
      a=-uh
      fu=log(-uh/q2)
      fa=fu !log(a/q2)
      end
      
      subroutine utilities(sh,th,uh,q2,ss2)
      implicit real *8 (a-h,o-z)
      real *8 dyLI2
c      common/fractions/x1,x2
      include 'scales2_inc.f'
c      common/scales2/xmur,xmuf,xmur2,xmuf2
      include 'functions_inc.f'
c      common/denominators/dt,du,ds,dst,dsu,dtu,la
c      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
c     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
c      common/recmass/s2
      INTEGER :: nthreads, myid
      INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
c     bug fix in DYqT: th <-> uh
c      a=-th
      a=-uh

c     pass s2 as a parameter to set it to 0 when appropriate, and improve numerical precision
      s2=ss2
c      s2=sh+th+uh-q2
c      if (s2.lt.1d-7) then
c         s2=0d0
c      else
c         print *,s2,(s2.lt.1d-7)
c      endif
c.....denominator factors
      dt=1d0/(s2-th)
      du=1d0/(s2-uh)
      ds=1d0/(sh+q2-s2)
      dst=1d0/(sh+th-s2)
      dsu=1d0/(sh+uh-s2)
      dtu=1d0/(th*uh-s2*q2)
      la=sqrt((uh+th)**2-4d0*s2*q2)
c.....transcendental functions
      fs=log(sh/q2)
      ft=log(-th/q2)
!      fu=log(-uh/q2)
c      fs2=log(s2/q2)
      if (s2.gt.0d0) fs2=log(s2/q2)
      else fs2 = 0d0       !when s2 is 0, fs2 is infinity!
!      fm2=log(xmuf2/q2)
!      fmu2=log(xmur2/q2)
!      fa=log(a/q2)
      fst=log(sh*th**2/(q2*(s2-th)**2))
      fsu=log(sh*uh**2/(q2*(s2-uh)**2))
      fla=log((sh+q2-s2+la)/(sh+q2-s2-la))
      fstu=log(sh*q2/((s2-th)*(s2-uh)))
      ftu=log((th*uh-s2*q2)/((s2-th)*(s2-uh)))
      flat=log(sh*q2*(s2-th)**2/(s2*(2*q2-uh)-q2*th)**2)
      flau=log(sh*q2*(s2-uh)**2/(s2*(2*q2-th)-q2*uh)**2)
!      f1t=dyLI2(q2/(q2-th))+1d0/2d0*(log(q2/(q2-th)))**2
!      f2t=dyLI2(q2/sh)+1d0/2d0*fs**2+fs*log(-th/(sh-q2))
!      f1u=dyLI2(q2/(q2-uh))+1d0/2d0*(log(q2/(q2-uh)))**2
!      f2u=dyLI2(q2/sh)+1d0/2d0*fs**2+fs*log(-uh/(sh-q2))
!      print *,OMP_GET_THREAD_NUM(),'mand',sh,uh,th,q2,s2
!      print *,OMP_GET_THREAD_NUM(),'dens',dt,du,ds,dst,dsu,dtu
!      print *,OMP_GET_THREAD_NUM(),'trans1',fs,ft,fu,fs2,fm2,fmu2,fa
!      print *,OMP_GET_THREAD_NUM(),'trans2',fst,fsu,fla,fstu,ftu,flat,flau
!      print *,OMP_GET_THREAD_NUM(),'s',s2,q2,s2,q2,fs2,sh,fs,f1t
      return
      end

      subroutine utilities_dilog(sh,th,uh,q2)
      implicit none
      real *8 sh,th,uh,q2
      real *8 dyLI2
      include 'functions_inc.f'
      f1t=dyLI2(q2/(q2-th))+1d0/2d0*(log(q2/(q2-th)))**2
      f2t=dyLI2(q2/sh)+1d0/2d0*fs**2+fs*log(-th/(sh-q2))
      f1u=dyLI2(q2/(q2-uh))+1d0/2d0*(log(q2/(q2-uh)))**2
      f2u=dyLI2(q2/sh)+1d0/2d0*fs**2+fs*log(-uh/(sh-q2))
!      print *,OMP_GET_THREAD_NUM(),'dilog',f1t,f2t,f1u,f2u
      return
      end
      
