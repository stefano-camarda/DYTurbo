c.....subroutine calling vegas
      subroutine xinteg(sig,indim,initn,incalls,avg,stddev,chi2a)
      implicit real * 8 (a-h,o-z)
      common/bveg1/xl(10),xu(10),acc,ndim,ncall,itmx,nprn
      common/bveg2/xi(50,10),si,si2,swgt,schi,ndo,it
      common/cwww/w1max,w1vgs,w1evt
      common/ciweight/iweight
      external sig

      do j=1,10
        xl(j)=0.d0
        xu(j)=1.d0
      enddo
      nprn=0
      acc=-1d0
      ndim=indim
      ncall=incalls
      itmx=3
      do j=1,initn
        if(j.eq.1)then
          call vegas(sig,avgi,sd,chi2a)
        else
          call vegas3(sig,avgi,sd,chi2a)
        endif
      enddo
      avg=avgi
      stddev=sd
      return
      end

c.....VEGAS subroutine
      subroutine vegas(fxn,avgi,sd,chi2a)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      common/bveg1/xl(10),xu(10),acc,ndim,ncall,itmx,nprn
      common/bveg2/xi(50,10),si,si2,swgt,schi,ndo,it
      COMMON/SEED/NUM,NUM2
      dimension d(50,10),di(50,10),xin(50),r(50),dx(10),dt(10),
     1     x(10),kg(10),ia(10)
      dimension RAND(10)
      data ndmx/50/,alph/1.5d0/,one/1.d0/,mds/1/
     
      NUM=1
c     NUM2 e' irrilevante
      ndo=1
      do 1 j=1,ndim
 1       xi(1,j)=one
         entry vegas1(fxn,avgi,sd,chi2a)
c........initialises  cumulative  variables but not grid
         it=0
         si=0.
         si2=si
         swgt=si
         schi=si
         entry vegas2(fxn,avgi,sd,chi2a)
c........no initialisation
         nd=ndmx
         ng=1
         if(mds.eq.0)go to 2
         ng=(ncall/2.)**(1./ndim)
         mds=1
         if((2*ng-ndmx).lt.0)go to 2
         mds=-1
         npg=ng/ndmx+1
         nd=ng/npg
         ng=npg*nd
 2       k=ng**ndim
         npg=ncall/k
         if(npg.lt.2)npg=2
         calls=npg*k
         dxg=one/ng
         dv2g=(calls*dxg**ndim)**2/npg/npg/(npg-one)
         xnd=nd
         ndm=nd-1
         dxg=dxg*xnd
         xjac=one/calls
         do 3 j=1,ndim
         dx(j)=xu(j)-xl(j)
 3       xjac=xjac*dx(j)
c........rebin preserving bin density
         if(nd.eq.ndo)go to 8
         rc=ndo/xnd
         do 7 J=1,ndim
         k=0
         xn=0.
         dr=xn
         i=k
 4       k=k+1
         dr=dr+one
         xo=xn
         xn=xi(k,j)
 5       if(rc.gt.dr)go to 4
         i=i+1
         dr=dr-rc
         xin(i)=xn-(xn-xo)*dr
         if(i.lt.ndm)go to 5
         do 6 i=1,ndm
 6       xi(i,j)=xin(i)
 7       xi(nd,j)=one
         ndo=nd
 8       if(nprn.ne.0)write(6,200)ndim,calls,it,itmx,acc
     1   ,mds,nd,(xl(j),xu(j),j=1,ndim)
         entry vegas3(fxn,avgi,sd,chi2a)
c........main integration loop
 9       it=it+1
          ti=0.
         tsi=ti
         do 10 j=1,ndim
         kg(j)=1
         do 10 i=1,nd
         d(i,j)=ti
 10      di(i,j)=ti
 11      fb=0.
         f2b=fb
         k=0
 12      k=k+1
       call randa(ndim,rand)
         wgt=xjac
         do 15 j=1,ndim
         xn=(kg(j)-rand(j))*dxg+one
         ia(j)=xn
         if(ia(j).gt.1)go to 13
         xo=xi(ia(j),j)
         rc=(xn-ia(j))*xo
         go to 14
13       xO=xi(ia(j),j)-xi(ia(j)-1,j)
         rc=xi(ia(j)-1,j)+(xn-ia(j))*xo
 14      x(j)=xl(j)+rc*dx(j)
 15      wgt=wgt*xo*xnd
c
         f=wgt
         f=f*fxn(x,wgt)
         f2=f*f
         fb=fb+f
         f2b=f2b+f2
         do 16 j=1,ndim
         di(ia(j),j)=di(ia(j),j)+f
 16      if(mds.ge.0)d(ia(j),J)=d(ia(j),J)+f2
         if(k.lt.npg) go to 12
888    FORMAT(1X,'F',G14.6,'F2',G14.6,'FB',G14.6,'F2B',G14.6)
         f2b= sqrt(f2b*      NPG)
         f2b=(f2b-fb)*(f2b+fb)
1661   FORMAT(1X,'F2B',G14.6,'NPG',  I10)
         ti=ti+fb
         tsi=tsi+f2b
33     FORMAT(1X,'TSI',G14.6,'F2B',G14.6)
         if(mds.ge.0)go to 18
         do 17 j=1,ndim
 17      d(ia(j),j)=d(ia(j),j)+f2b
 18      k=ndim
 19      kg(k)=mod(kg(k),ng)+1
         if(kg(k).ne.1)go to 11
         k=k-1
         if(k.gt.0)go to 19
c........final results for this iteration
        tsi=tsi*dv2g
        ti2=ti*ti
88     format(1x,'tsi',g14.6)
       if(tsi.eq.0.d0)then 
         rpp=1.d15
       else
         rpp=abs(ti2/tsi)
       endif
       if(rpp.lt.1.d14) then
        wgt=ti2/tsi
        si=si+ti*wgt
        si2=si2+ti2
        swgt=swgt+wgt
        schi=schi+ti2*wgt
995    FORMAT(1X,'SWGT',G14.6,'SI2',G14.6)
        avgi=si/swgt
        sd=swgt*it/si2
        chi2a=sd*(schi/swgt-avgi*avgi)/(it-.999)
        sd=dsqrt(one/sd)
       else
        write(*,*) '***VEGAS WARNING***: zero error!'
        write(*,*) ' we guess that integral is exact '
        avgi=ti
        si=ti
        si2=0.d0
        sd=0
        chi2a=-1
        return
       endif
c
        if(nprn.eq.0)go to 21
        tsi=dsqrt(tsi)
        write(6,201)it,ti,tsi,avgi,sd,chi2a
        if(nprn.ge.0)go to 21
        do 20 j=1,ndim
 20     write(6,202) j,(xi(i,j),di(i,j),d(i,j),i=1,nd)
c.......refine grid
 21     do 23 j=1,ndim
        xo=d(1,j)
        xn=d(2,j)
        d(1,j)=(xo+xn)/2.
        dt(j)=d(1,j)
        do 22 i=2,ndm
        d(i,j)=xo+xn
        xo=xn
        xn=d(i+1,j)
        d(i,j)=(d(i,j)+xn)/3.
 22     dt(j)=dt(j)+d(i,j)
        d(nd,j)=(xn+xo)/2.
 23     dt(j)=dt(j)+d(nd,j)
        do 28 j=1,ndim
        rc=0.
        do 24 i=1,nd
        r(i)=0.
        if(d(i,j).le.0.)go to 24
        xo=dt(j)/d(i,j)
        r(i)=((xo-one)/xo/dlog(xo))**alph
 24     rc=rc+r(i)
        rc=rc/xnd
        k=0
        xn=0.
        dr=xn
        i=k
 25     k=k+1
        dr=dr+r(k)
        xo=xn
        xn=xi(k,j)
 26     if(rc.gt.dr)go to 25
        i=i+1
        dr=dr-rc
        xin(i)=xn-(xn-xo)*dr/r(k)
        if(i.lt.ndm)go to 26
        do 27 i=1,ndm
 27     xi(i,j)=xin(i)
 28     xi(nd,j)=one
        if(it.lt.itmx.and.acc*dabs(avgi).lt.sd)go to 9
 200    format(1X,'0input parameters for vegas:  ndim=',i3,
     1  '   ncall=',f8.0/28x,'  it=',i5,'    itmx=',i5/28x,
     2  '  acc=',g9.3/28x,'  mds=',i3,'     nd=',i4/28x,
     3  '  (xl,xu)=',(t40,'( ',g12.6,' , ',g12.6,' )'))
 201    format(///' integration by vegas' / '0iteration no.',i5,
     1  ':  integral=',g14.8/21x,'std dev =',g10.4 /
     2  ' accumulated results:   integral=',g14.8/
     3  24x,'std dev =',g10.4 / 24x,'chi**2 per it''n =',g10.4)
 202    format(1X,'0data for axis',i2,/,' ',6x,'x',7x,'  delt i ',
     1  2x,'conv','ce   ',11x,'x',7x,'  delt i ',2x,'conv','ce  '
     2  ,11x,'x',7x,'   delt i ',2x,'conv','CE  ',/,
     3  (1X,' ',3g12.4,5x,3g12.4,5x,3g12.4))
        return
        entry vegas4(fxn,avgi,sd,chi2a)
        if(si2.eq.0.d0)then
          avgi=si
          sd=0
          chi2a=-1
        else
          avgi=si/swgt
          sd=swgt*it/si2
          chi2a=sd*(schi/swgt-avgi*avgi)/(it-.999)
          sd=dsqrt(one/sd)
        endif
        if(nprn.ne.0) write(6,201)it,0.d0,0.d0,avgi,sd,chi2a
        return
        end

c        subroutine save(ndim)
c        implicit real*8 (a-h,o-z)
c       implicit integer*4 (i-n)
c        common/bveg2/xi(50,10),si,si2,swgt,schi,ndo,it
c
c      stores vegas data   (unit 7) for later initialisation
c
c        write(7,200) ndo,it,si,si2,swgt,schi,
c     1       ((xi(i,j),i=1,ndo),j=1,ndim)
c        return
c        entry restr(ndim)
c
c         enters initialisation data for vegas
c
c        read(7,200) ndo,it,si,si2,swgt,schi,
c     1    ((xi(i,j),i= 1,ndo),j=1,ndim)
c 200    format(2i8,4z16/(5z16))
c        return
c        end

      
      FUNCTION RANDOM(SEED)
*     -----------------
* Ref.: K. Park and K.W. Miller, Comm. of the ACM 31 (1988) p.1192
* Use seed = 1 as first value.
*
      IMPLICIT INTEGER*4 (A-Z)
      DOUBLE PRECISION MINV,RANDOM
      SAVE
      PARAMETER(M=2147483647,A=16807,Q=127773,R=2836)
      PARAMETER(MINV=0.46566128752458d-09)
      HI = SEED/Q
      LO = MOD(SEED,Q)
      SEED = A*LO - R*HI
      IF(SEED.LE.0) SEED = SEED + M
      RANDOM = SEED*MINV
      END

      subroutine randa(n,rand)
      implicit double precision (a-h,o-z)
      COMMON/SEED/NUM,NUM2
      common/caso/caso(5)
      dimension rand(10)
      do 1 i=1,n
      rand(i)=random(NUM)
1     continue
      do 2 i=1,5
      caso(i)=random(NUM)
2     continue
      return
      end


      subroutine sysdep
c use newver='NEW' for vaxes, 'UNKNOWN' in other machines,
      character * 7 newver
      common/newver/newver
      newver = 'UNKNOWN'
      end

      subroutine delete(fname)
      character * 80 str
      character * (*) fname
      l = len(fname)
      k = 1
      dowhile(fname(k:k).eq.' '.and.k.lt.l)
         k = k+1
      enddo
      dowhile(fname(l:l).eq.' '.and.l.gt.k)
         l = l-1
      enddo
      if(l-k.gt.70) then
         write(*,*) 'delete: filename > 70 chars not allowed'
         stop
      endif
      if(l.eq.k) then
         write(*,*) 'delete: void filename'
         stop
      endif
      str(1:) = '\\rm '
      str(7:) = fname(k:l)
      call system(str)
      end
