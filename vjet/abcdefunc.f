c**********************************************
c.....LO contributions for QG->V+X and GQ->V+X
c**********************************************

c.....Real contribution to QG->V+X
      function Aqg0(sh,th,uh,q2)
      implicit none
      double precision aqg0
      double precision sh,th,uh,q2
      double precision temp
      include 'luminosities_inc.f'

      temp=(-sh/th-th/sh-(2d0*uh*q2)/(sh*th))
      Aqg0=temp*xlumqg
      
      return
      end

c.....Real contribution to GQ->V+X
      function Agq0(sh,th,uh,q2)
      implicit none
      double precision agq0
      double precision sh,th,uh,q2
      double precision temp
      include 'luminosities_inc.f'

      temp=(-sh/uh-uh/sh-(2d0*th*q2)/(sh*uh))
      Agq0=temp*xlumgq

      return
      end

c**********************************************
c.....LO contribution for QQB->V+X (=QBQ->V+X)
c**********************************************

c.....Real contribution to QQB->V+X
      function Aqqb0(sh,th,uh,q2)
      implicit none
      double precision aqqb0
      double precision sh,th,uh,q2
      double precision temp
      include 'luminosities_inc.f'

      temp=(uh/th+th/uh+(2d0*sh*q2)/(uh*th))
      Aqqb0=temp*xlumqqb

      return
      end

c******************************************************
c.....NLO virtual contributions to QG->V+X and GQ->V+X
c******************************************************

c.....Virtual (one-loop) contribution to QG->V+X
      function Bqg1(sh,th,uh,q2)
      implicit none
      double precision bqg1
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
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
      implicit none
      double precision bgq1
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
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
      implicit none
      double precision bqg2
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
      temp=(sh/th+th/sh+(2d0*uh*q2)/(sh*th))*fmu2/3d0*nf
      Bqg2=temp*xlumqg

      return
      end

c.....Virtual (renormalization counterterm) contribution to GQ->V+X
      function Bgq2(sh,th,uh,q2)
      implicit none
      double precision bgq2
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
      temp=(sh/uh+uh/sh+(2d0*th*q2)/(sh*uh))*fmu2/3d0*nf
      Bgq2=temp*xlumgq

      return
      end

c.....Virtual (triangular quark loops) contribution to QG->V+X
      function Bqg3(sh,th,uh,q2)
      implicit none
      double precision bqg3
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
      temp=(uh+q2)/(uh-q2)*(1-(q2*fu)/(uh-q2))
      Bqg3=temp*xlumqgtr

      return
      end

c.....Virtual (triangular quark loops) contribution to GQ->V+X
      function Bgq3(sh,th,uh,q2)
      implicit none
      double precision bgq3
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
      temp=(th+q2)/(th-q2)*(1-(q2*ft)/(th-q2))
      Bgq3=temp*xlumgqtr

      return
      end

c********************************************************
c.....NLO virtual contributions for QQB->V+X (=QBQ->V+X)
c********************************************************

c.....Virtual (one-loop) contribution to QQB->V+X
      function Bqqb1(sh,th,uh,q2)
      implicit none
      double precision bqqb1
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
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
      implicit none
      double precision bqqb2
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
      temp=-(uh/th+th/uh+(2d0*sh*q2)/(uh*th))*fmu2/3d0*nf
      Bqqb2=temp*xlumqqb

      return
      end

c.....Virtual (triangular quark loops) contribution to QQB->V+X
      function Bqqb3(sh,th,uh,q2)
      implicit none
      double precision bqqb3
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
      temp=-(sh+q2)/(sh-q2)*(1-(q2*fs)/(sh-q2))
      Bqqb3=temp*xlumqqbtr

      return
      end

c***************************************************
c.....NLO real contributions to QG->VQG and GQ->VQG
c***************************************************

c.....First Real contribution to QG->VQG (proportional to delta(s2))
      function Cqg1(sh,th,uh,q2)
      implicit none
      double precision cqg1
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
      temp=(-sh/th-th/sh-(2d0*uh*q2)/(sh*th))*
     /     (cf*(7/2d0+2*fm2*(fu-fa)+fa**2-3/2d0*(fm2+fa))+
     /     ca*(pi**2/6d0+1/2d0*(fs-ft-fu)**2+2*ft*(fm2-fa)+
     /     2*fa*(fs-fu-fm2+fa)-11/6d0*fm2))
      Cqg1=temp*xlumqg

      return
      end

c.....First Real contribution to GQ->VQG (proportional to delta(s2))
      function Cgq1(sh,th,uh,q2)
      implicit none
      double precision cgq1
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
      temp=(-sh/uh-uh/sh-(2d0*th*q2)/(sh*uh))*
     /     (cf*(7/2d0+2*fm2*(ft-fa)+fa**2-3/2d0*(fm2+fa))+
     /     ca*(pi**2/6d0+1/2d0*(fs-fu-ft)**2+2*fu*(fm2-fa)+
     /     2*fa*(fs-ft-fm2+fa)-11/6d0*fm2))
      Cgq1=temp*xlumgq

      return
      end

c.....Second Real contribution to QG->VQG (proportional to delta(s2))
      function Cqg2(sh,th,uh,q2)
      implicit none
      double precision cqg2
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
      temp=(-sh/th-th/sh-(2d0*uh*q2)/(sh*th))*fm2/3d0*nf
      Cqg2=temp*xlumqg

      return
      end

c.....Second Real contribution to GQ->VQG (proportional to delta(s2))
      function Cgq2(sh,th,uh,q2)
      implicit none
      double precision cgq2
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
      temp=(-sh/uh-uh/sh-(2d0*th*q2)/(sh*uh))*fm2/3d0*nf
      Cgq2=temp*xlumgq

      return
      end

c.....Third Real contribution to QG->VQG (for nonzero values of s2)
      function Cqg3(sh,th,uh,q2,what)
      implicit none
      double precision cqg3
      double precision sh,th,uh,q2
      double precision temp
      integer what
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
      if (what.eq.1) then
c     term proportional to (fs2/s2)_A+
      temp=(-sh/th-th/sh-(2d0*uh*q2)/(sh*th))*(2d0*cf+4d0*ca)
      elseif (what.eq.2) then
c     term proportional to (1/s2)_A+
      temp=(-sh/th-th/sh-(2d0*uh*q2)/(sh*th))*(
     /     cf*(-3/2d0+fsu+2d0*ftu-2d0*fm2+(th+uh)*fla/la)+
     /     ca*(2d0*fstu+(fst-fsu)/2d0-ftu-2d0*fm2))
      elseif (what.eq.3) then
c     term proportional to (1/s2), FINITE in the limit s2->0
      temp=(cf-ca/2d0)*((sh+2d0*uh)/th*((th+uh)*fla/la+fsu)+
     /     (2d0*th+4d0*uh)*ftu/sh)
      temp=temp/s2
      elseif (what.eq.4) then
c     regular term
      temp=fla/la**3*(cf*((3*sh)*(sh+q2-s2)*(uh**2-th**2)*
     /     (1-uh/th)/(4*la**2)+(uh**2-th**2)*((q2-s2)/(4*sh)+ ! bug fix !!! !     /     (1-uh/th)/(4*la**2)+(uh**2-th**2)*((q2-s2**2)/(4*sh)+
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
      endif
      Cqg3=temp*xlumqg

      return
      end

c.....Third Real contribution to GQ->VQG (for nonzero values of s2)
      function Cgq3(sh,th,uh,q2,what)
      implicit none
      double precision cgq3
      double precision sh,th,uh,q2
      double precision temp
      integer what
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'

      if (what.eq.1) then
c     term proportional to (fs2/s2)_A+
      temp=(-sh/uh-uh/sh-(2d0*th*q2)/(sh*uh))*(2d0*cf+4d0*ca)
      elseif (what.eq.2) then
c     term proportional to (1/s2)_A+
      temp=(-sh/uh-uh/sh-(2d0*th*q2)/(sh*uh))*(
     /     cf*(-3/2d0+fst+2d0*ftu-2d0*fm2+(uh+th)*fla/la)+
     /     ca*(2d0*fstu+(fsu-fst)/2d0-ftu-2d0*fm2))
      elseif (what.eq.3) then
c     term proportional to (1/s2), FINITE in the limit s2->0
      temp=(cf-ca/2d0)*((sh+2d0*th)/uh*((uh+th)*fla/la+fst)+
     /     (2d0*uh+4d0*th)*ftu/sh)
      temp=temp/s2
      elseif (what.eq.4) then
c     regular term
      temp=fla/la**3*(cf*((3*sh)*(sh+q2-s2)*(th**2-uh**2)*
     /     (1-th/uh)/(4*la**2)+(th**2-uh**2)*((q2-s2)/(4*sh)+ ! bug fix !!!!     /     (1-th/uh)/(4*la**2)+(th**2-uh**2)*((q2-s2**2)/(4*sh)+
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
      endif
      Cgq3=temp*xlumgq

      return
      end

c****************************************
c.....NLO real contributions to GG->VQQB
c****************************************

c.....First Real contribution to GG->VQQB 
      function Cgg1(sh,th,uh,q2)
      implicit none
      double precision cgg1
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'

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
      implicit none
      double precision cgg1x
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'

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
      implicit none
      double precision cqqb1
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
      temp=(uh/th+th/uh+(2d0*sh*q2)/(uh*th))*
     /     (ca*(67/18d0-11*fa/6d0+fa**2)+cf*(2*ft+2*fu-4*fa-3)*fm2+
     /     (cf-ca/2d0)*(pi**2/3d0+(2*fa+fs-ft-fu)**2))
      Cqqb1=temp*xlumqqb

      return
      end

c.....Second Real contribution to QQB->VGG (for nonzero values of s2)
      function Cqqb2(sh,th,uh,q2,what)
      implicit none
      double precision cqqb2
      double precision sh,th,uh,q2
      double precision temp
      integer what
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
      if (what.eq.1) then
c     term proportional to (fs2/s2)_A+
      temp=1/2d0*(uh/th+th/uh+2*q2*sh/(th*uh))*(8*cf-2*ca)
      elseif (what.eq.2) then
c     term proportional to (1/s2)_A+
      temp=1/2d0*(uh/th+th/uh+2*q2*sh/(th*uh))*
     /     (-11*ca/6d0+2*cf*(ftu-2*fm2)+(2*cf-ca)*(2*fstu-ftu))
      elseif (what.eq.4) then
c     regular term
      temp=cf*(sh*dt**2+2*dt-sh/(th*uh)+dtu*(s2+2*q2*(uh-s2)/th)-
     /     2*dt*q2*(1/uh-1/th))+
     /     ca*(-11*sh/(6*th*uh)+dt**2*sh**2/uh*(3*s2/(2*th)-2)+
     /     2*dt*sh/uh+q2/(3*th**2))+
     /     (2*cf-ca)*(dt/uh*(q2-uh)**2*(du-1/th)*(fstu+fs2-ftu)+
     /     (sh-q2)/(th*uh)*(fstu+fs2))+
     /     ftu*(cf*dtu*(4*q2/(th*uh)*(q2-th)**2+2*(q2+sh)-s2)+
     /     cf*(s2-2*sh)/(th*uh)-ca*q2/(th*uh))+
     /     cf*(fm2-fs2)*(dtu*(4*(uh-q2)-s2-4*q2/(th*uh)*(uh-q2)**2)+
     /     dt**2*(q2-uh)-dt/th*(2*q2-uh)-2/th+q2/th**2+4*q2/(th*uh))
      endif
      Cqqb2=temp*xlumqqb

      return
      end

c.....Crossed Second Real contribution to QQB->VGG (for nonzero values of s2)
      function Cqqb2x(sh,th,uh,q2,what)
      implicit none
      double precision cqqb2x
      double precision sh,th,uh,q2
      double precision temp
      integer what
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'

      if (what.eq.1) then
c     term proportional to (fs2/s2)_A+
      temp=1/2d0*(th/uh+uh/th+2*q2*sh/(uh*th))*(8*cf-2*ca)
      elseif (what.eq.2) then
c     term proportional to (1/s2)_A+
      temp=1/2d0*(th/uh+uh/th+2*q2*sh/(uh*th))*
     /     (-11*ca/6+2*cf*(ftu-2*fm2)+(2*cf-ca)*(2*fstu-ftu))
      elseif (what.eq.4) then
c     regular term
      temp=cf*(sh*du**2+2*du-sh/(uh*th)+dtu*(s2+2*q2*(th-s2)/uh)-
     /     2*du*q2*(1/th-1/uh))+
     /     ca*(-11*sh/(6*uh*th)+du**2*sh**2/th*(3*s2/(2*uh)-2)+
     /     2*du*sh/th+q2/(3*uh**2))+
     /     (2*cf-ca)*(du/th*(q2-th)**2*(dt-1/uh)*(fstu+fs2-ftu)+
     /     (sh-q2)/(uh*th)*(fstu+fs2))+
     /     ftu*(cf*dtu*(4*q2/(uh*th)*(q2-uh)**2+2*(q2+sh)-s2)+
     /     cf*(s2-2*sh)/(uh*th)-ca*q2/(uh*th))+
     /     cf*(fm2-fs2)*(dtu*(4*(th-q2)-s2-4*q2/(uh*th)*(th-q2)**2)+
     /     du**2*(q2-th)-du/uh*(2*q2-th)-2/uh+q2/uh**2+4*q2/(uh*th))
      endif
      Cqqb2x=temp*xlumqqb

      return
      end

c********************************************************
c.....NLO real contributions for QQB->VQQB (=QBQ->VQQB)
c********************************************************

c.....First Real contribution to QQB->VQQB
      function D0aa(sh,th,uh,q2)
      implicit none
      double precision d0aa
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
      temp=(uh/th+th/uh+(2d0*sh*q2)/(uh*th))*(fa/3d0-5/9d0)*nf
      D0aa=temp*xlumqqb

      return
      end

c.....Second Real contribution to QQB->VQQB
      function Daa(sh,th,uh,q2,what)
      implicit none
      double precision daa
      double precision sh,th,uh,q2
      double precision temp
      integer what
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
      if (what.eq.2) then
c     term proportional to (1/s2)_A+
      temp=1/6d0*(uh/th+th/uh+(2*sh*q2)/(th*uh))
      temp=temp*nf
      elseif (what.eq.4) then
c     regular term
      temp=1/3d0*(sh/(th*uh)-q2/th**2)
      temp=temp*nf
      endif
      Daa=temp*xlumqqb

      return
      end

c.....Second Real contribution to QQB->VQQB (crossed)
      function Daax(sh,th,uh,q2,what)
      implicit none
      double precision daax
      double precision sh,th,uh,q2
      double precision temp
      integer what
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
      if (what.eq.2) then
c     term proportional to (1/s2)_A+
      temp=1/6d0*(th/uh+uh/th+(2*sh*q2)/(uh*th))
      temp=temp*nf
      elseif (what.eq.4) then
c     regular term
      temp=1/3d0*(sh/(uh*th)-q2/uh**2)
      temp=temp*nf
      endif
      Daax=temp*xlumqqb

      return
      end

c.....Third real (only for Z production) contribution to QQB->VQQB
      function Dab(sh,th,uh,q2)
      implicit none
      double precision dab
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
      temp=1/2d0*(fla*q2/(la*th)*((2*sh+th+uh)/la**2*(3*sh*(th-uh)**2/
     /     la**2+uh-th-sh)-1)+(2*sh+th+uh)/(2*la**2*th)*(-3/la**2*
     /     (sh+q2-s2)*(th-uh)**2+3*th-uh+2*sh-2*s2))
      Dab=temp*xlumqqbtr

      return
      end
      
c.....Third real (only for Z production) contribution to QQB->VQQB (crossed)
      function Dabx(sh,th,uh,q2)
      implicit none
      double precision dabx
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
      temp=1/2d0*(fla*q2/(la*uh)*((2*sh+uh+th)/la**2*(3*sh*(uh-th)**2/
     /     la**2+th-uh-sh)-1)+(2*sh+uh+th)/(2*la**2*uh)*(-3/la**2*
     /     (sh+q2-s2)*(uh-th)**2+3*uh-th+2*sh-2*s2))
      Dabx=temp*xlumqqbtr

      return
      end
      
c.....Fourth real contribution to QQB->VQQB
      function Dbb(sh,th,uh,q2)
      implicit none
      double precision dbb
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
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
      implicit none
      double precision dbbx
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
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
      implicit none
      double precision dac
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
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
      implicit none
      double precision dbc
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
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
      implicit none
      double precision dad
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
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
      implicit none
      double precision dbd
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
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
      implicit none
      double precision dcc
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
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
      implicit none
      double precision ddd
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
       
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
      implicit none
      double precision dcdll
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
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
      implicit none
      double precision dcdllx
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
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
      implicit none
      double precision dcdlr
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
      temp=1/2d0*(2*ds*sh**2/th*(dst*(fsu-flat-2*ftu)+2/uh*(fst-ftu))+
     /     dst*(2*ftu-fsu+flat)*(1+2*sh/th-uh/th)+2d0/(th*uh)*
     /     (ftu-fsu)*(sh+s2-uh)-2d0/th*(ftu+flat))
      DcdLR=temp*xlumqqbLR

      return
      end

c.....Twelth real contribution to QQB->VQQB (crossed)
      function DcdLRx(sh,th,uh,q2)
      implicit none
      double precision dcdlrx
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
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
      implicit none
      double precision eac
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
      temp=(cf-ca/2d0)*((fst+flat)*(4*dt*ds*sh**2/th-2*dt*(1-uh/th)-
     /     4/th)+4*q2/th**2*fst)
      Eac=temp*xlumqq

      return
      end

c.....Second contribution which is different from D functions
      function Ead(sh,th,uh,q2)
      implicit none
      double precision ead
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
      temp=(cf-ca/2d0)*(1/sh*(2*fs2-fst-fsu+2*fstu)*(-2-uh/th-th/uh+
     /     2*s2/(th*uh)*(q2-sh))-2/sh*(2+uh/th+th/uh-2*s2*q2*
     /     (1/th**2+1/uh**2)))
      Ead=temp*xlumqqead

      return
      end

c.....Third contribution which is different from D functions (Eac crossed)
      function Ebd(sh,th,uh,q2)
      implicit none
      double precision ebd
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
      temp=(cf-ca/2d0)*((fsu+flau)*(4*du*ds*sh**2/uh-2*du*(1-th/uh)-
     /     4/uh)+4*q2/uh**2*fsu)
      Ebd=temp*xlumqq

      return
      end

c.....Fourth contribution which is different from D functions (Ead crossed)
      function Ebc(sh,th,uh,q2)
      implicit none
      double precision ebc
      double precision sh,th,uh,q2
      double precision temp
      include 'dyqcd.f'
      include 'functions_inc.f'
      include 'luminosities_inc.f'
      
      temp=(cf-ca/2d0)*(1/sh*(2*fs2-fsu-fst+2*fstu)*(-2-th/uh-uh/th+
     /     2*s2/(uh*th)*(q2-sh))-2/sh*(2+th/uh+uh/th-2*s2*q2*
     /     (1/uh**2+1/th**2)))
      Ebc=temp*xlumqqead

      return
      end

c     Eaa=Ecc=Dcc (A28)
c     Ebb=Edd=Ddd (A29)
c     EabLL=EcdLL=-DcdLR (A30)
c     EabLR=EcdLR=-DcdLL (A31)
