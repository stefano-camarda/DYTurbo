C-----------------------------------------------------------------------
C                              Vegas
C
      BLOCK DATA vegdat
C
C     makes default parameter assignments for vegas
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      COMMON/bveg1/ncall,itmx,nprn,ndev,xl(10),xu(10),acc
      COMMON/bveg2/it,ndo,si,swgt,schi,xi(50,10)
      COMMON/bveg3/alph,ndmx,mds
      COMMON/rnsd/iseed
      DATA ncall/5000/,itmx/5/,nprn/5/,acc/-1./,
     1     xl/0.,0.,0.,0.,0.,0.,0.,0.,0.,0./,
     2     xu/1.,1.,1.,1.,1.,1.,1.,1.,1.,1./,
     3     alph/1.5/,ndmx/50/,mds/1/,ndev/6/,
     4     ndo/1/,xi/500*1./,it/0/,si,swgt,schi/3*0./

      END
C-----------------------------------------------------------------------
      SUBROUTINE VEGAS(ndim,fxn,avgi,sd,chi2a)
C
C     Subroutine performs ndim-dimensional Monte Carlo integ'n
C     - by G.P. Lepage    Sept 1976/(rev)Aug 1979
C     - algorithm described in J Comp Phys 27,192(1978)
C
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      COMMON/bveg1/ncall,itmx,nprn,ndev,xl(10),xu(10),acc
      COMMON/bveg2/it,ndo,si,swgt,schi,xi(50,10)
      COMMON/bveg3/alph,ndmx,mds
      COMMON/bveg4/calls,ti,tsi
      DIMENSION d(50,10),di(50,10),xin(50),r(50),dx(10),ia(10),
     1          kg(10),dt(10),x(10)
C      DIMENSION rand(10)
      double precision rand(10)
      DATA one/1d0/
      sqrt(a)=dsqrt(a)
      alog(a)=dlog(a)
      abs(a)=dabs(a)
C
      ndo=1
      DO 1 j=1,ndim
1     xi(1,j)=one
C
      ENTRY vegas1(ndim,fxn,avgi,sd,chi2a)
C        - initializes cumulative variables, but not grid
      it=0
      si=0d0
      swgt=si
      schi=si
C
      ENTRY vegas2(ndim,fxn,avgi,sd,chi2a)
C        - no initialization
      nd=ndmx
      ng=1
      IF(mds.EQ.0) GO TO 2
      ng=(ncall/2d0)**(1d0/ndim)
      mds=1
      IF((2*ng-ndmx).LT.0) GO TO 2
      mds=-1
      npg=ng/ndmx+1
      nd=ng/npg
      ng=npg*nd
2     k=ng**ndim
      npg=ncall/k
      IF(npg.LT.2) npg=2
      calls=npg*k
      dxg=one/ng
      dv2g=(calls*dxg**ndim)**2/npg/npg/(npg-one)
      xnd=nd
      ndm=nd-1
      dxg=dxg*xnd
      xjac=one/calls
      DO 3 j=1,ndim
      dx(j)=xu(j)-xl(j)
3     xjac=xjac*dx(j)
C
C   rebin, preserving bin density
      IF(nd.EQ.ndo) GO TO 8
      rc=ndo/xnd
      DO 7 j=1, ndim
      k=0
      xn=0.
      dr=xn
      i=k
4     k=k+1
      dr=dr+one
      xo=xn
      xn=xi(k,j)
5     IF(rc.GT.dr) GO TO 4
      i=i+1
      dr=dr-rC
      xin(i)=xn-(xn-xo)*dr
      IF(i.LT.ndm) GO TO 5
      DO 6 i=1,ndm
6     xi(i,j)=xin(i)
7     xi(nd,j)=one
      ndo=nd
C
8     IF(nprn.ge.0) WRITE(ndev,200) ndim,calls,it,itmx,acc,nprn,
     1                    alph,mds,nd,(xl(j),xu(j),j=1,ndim)
      IF(nprn.ge.0) WRITE(ndev,222)
C
      ENTRY vegas3(ndim,fxn,avgi,sd,chi2a)
C        - main integration loop
9     it=it+1
      ti=0d0
      tsi=ti
      DO 10 j=1,ndim
      kg(j)=1
      DO 10 i=1,nd
      d(i,j)=ti
10    di(i,j)=ti
C
11    fb=0d0
      f2b=fb
      k=0
12    k=k+1
      CALL randa(ndim,rand)
      wgt=xjaC
      DO 15 j=1,ndim
      xn=(kg(j)-rand(j))*dxg+one
      ia(j)=xn
      IF(ia(j).GT.1) GO TO 13
      xo=xi(ia(j),j)
      rc=(xn-ia(j))*xo
      GO TO 14
13    xo=xi(ia(j),j)-xi(ia(j)-1,j)
      rc=xi(ia(j)-1,j)+(xn-ia(j))*xo
14    x(j)=xl(j)+rc*dx(j)
15    wgt=wgt*xo*xnd
C
      f=wgt
      f=f*fxn(x,wgt)
      f2=f*f
      fb=fb+f
      f2b=f2b+f2
      DO 16 j=1,ndim
      di(ia(j),j)=di(ia(j),j)+f
16    IF(mds.ge.0) d(ia(j),j)=d(ia(j),j)+f2
      IF(k.LT.npg) GO TO 12
C
      f2b=sqrt(f2b*npg)
      f2b=(f2b-fb)*(f2b+fb)
      ti=ti+fb
      tsi=tsi+f2b
      IF(mds.ge.0) GO TO 18
      DO 17 j=1,ndim
17    d(ia(j),j)=d(ia(j),j)+f2b
18    k=ndim
19    kg(k)=mod(kg(k),ng)+1
      IF(kg(k).ne.1) GO TO 11
      k=k-1
      IF(k.GT.0) GO TO 19
C
C   compute final results for this iteration
      tsi=tsi*dv2g
      ti2=ti*ti
      wgt=one/tsi
      si=si+ti*wgt
      swgt=swgt+wgt
      schi=schi+ti2*wgt
      avgi=si/swgt
      chi2a=(schi-si*avgi)/(it-.9999d0)
      sd=sqrt(one/swgt)
C
      IF(nprn.LT.0) GO TO 21
      tsi=sqrt(tsi)
      WRITE(ndev,201) it,ti,tsi,avgi,sd,chi2a
      IF(nprn.EQ.0) GO TO 21
      DO 20 j=1,ndim
20    WRITE(ndev,202) j,(xi(i,j),di(i,j),i=1+nprn/2,nd,nprn)
C
c  refine grid
21    DO 23 j=1,ndim
      xo=d(1,j)
      xn=d(2,j)
      d(1,j)=(xo+xn)/2.
      dt(j)=d(1,j)
      DO 22 i=2,ndm
      d(i,j)=xo+xn
      xo=xn
      xn=d(i+1,j)
      d(i,j)=(d(i,j)+xn)/3.
22    dt(j)=dt(j)+d(i,j)
      d(nd,j)=(xo+xn)/2.
23    dt(j)=dt(j)+d(nd,j)
C
      DO 28 j=1,ndim
      rc=0.
      DO 24 i=1,nd
      r(i)=0d0
      IF(d(i,j).le.0.) GO TO 24
      xo=dt(j)/d(i,j)
      r(i)=((xo-one)/xo/alog(xo))**alph
24    rc=rc+r(i)
      rc=rc/xnd
      k=0
      xn=0d0
      dr=xn
      i=k
25    k=k+1
      dr=dr+r(k)
      xo=xn
      xn=xi(k,j)
26    IF(rc.GT.dr) GO TO 25
      i=i+1
      dr=dr-rC
      xin(i)=xn-(xn-xo)*dr/r(k)
      IF(i.LT.ndm) GO TO 26
      DO 27 i=1,ndm
27    xi(i,j)=xin(i)
28    xi(nd,j)=one
C
      IF(it.LT.itmx.AND.acc*abs(avgi).LT.sd) GO TO 9
200   FORMAT(/35h Input parameters for Vegas:  ndim=,i3,8h  ncall=,f8.0
     1  /28x,5h  it+,i5,7h  itmx=,i5/28x,6h  acc=,g9.3
     2  /28x,7h  nprn=,i3,7h  alph=,f5.2/28x,6h  mds=,i3,6h   nd=,i4
     3  /28x,10h  (xl,xu)=,(t40,2h( ,g12.6,3h , ,g12.6,2h )))
222   FORMAT(
     4  //1x,'Iter.#. Integral    Std.-Dev.   Accum.-Int. Std.-Dev.   ',
     5   'Chi2a',/)
201   FORMAT(1x,i3,5x,5(g9.3,3x))
*201   FORMAT(///21h integration by vegas//14h iteration no.,i3,
*     1  14h:   integral =,g14.8/24x,10hstd dev  =,g10.4/
*     2  34h accumulated results:   integral =,g14.8/
*     3  24x,10hstd dev  =,g10.4/24x,17hchi**2 per it'n =,g10.4)
202   FORMAT(/15h DATA for axis ,i2/25h    x       delta i       ,
     1  24h   x       delta i      ,18h   x       delta i,
     2  /(1h ,f7.6,1x,g11.4,5x,f7.6,1x,g11.4,5x,f7.6,1x,g11.4))
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE randa(n,rand)
C     subroutine generates uniformly distributed random no's x(i),i=1,n
      IMPLICIT NONE
      INTEGER n
      DOUBLE PRECISION rand(n),ran2,ranf
      EXTERNAL ran2
      INTEGER i,idum
      COMMON /ranseed/ idum
      DO i=1,n
C        rand(i)=ran2(idum)
         rand(i)=ranf(idum)
      ENDDO
      RETURN
      END
C-----------------------------------------------------------------------
      BLOCK DATA seeddata
      IMPLICIT NONE
      INTEGER idum
      COMMON /ranseed/ idum
      DATA idum /-123456789/
      END
C-----------------------------------------------------------------------
      SUBROUTINE getseed(seed)
      IMPLICIT NONE
      INTEGER seed,idum
      COMMON /ranseed/ idum
      seed=idum
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE setseed(seed)
      IMPLICIT NONE
      INTEGER seed,idum
      COMMON /ranseed/ idum
      IF (seed.GT.0) THEN
          idum=-seed
      ELSE
          idum=seed
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      double precision FUNCTION ran2(idum)
*     from Numerical recipes
*     converted to double precision
      implicit none
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      double precision AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=3d-16,RNMX=1d0-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION RANF(ISEED)
C     Park-Miller generator using Schrage's algorithm
      IMPLICIT NONE
      INTEGER ISEED
      INTEGER IA,IC,IQ,IR
      PARAMETER (IA=16807,IC=2147483647,IQ=127773,IR=2836)
      INTEGER IH,IL,IT
      IH=ISEED/IQ
      IL=MOD(ISEED,IQ)
      IT=IA*IL-IR*IH
      IF(IT.GT.0) THEN
          ISEED=IT
      ELSE
          ISEED=IC+IT
      ENDIF
      RANF=ISEED/DBLE(IC)
      END
