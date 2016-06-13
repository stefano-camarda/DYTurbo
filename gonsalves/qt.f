************************************************************************
*     Numerical Program for Inclusive gamma*, W & Z production         *
*     at large transverse momentum in hadron-hadron collisions         *
*     Includes next-to-leading order QCD corrections given in          *
*     R.J. Gonsalves, J. Pawlowski, C.-F. Wai, Phys. Rev. D40,         *
*     2245 (1989).                                                     *
*----------------------------------------------------------------------*
*     UBVM FILE  :  QT FORTRAN                                         *
*     Created    :  1989                                               *
*     Revised    :  June 30, 1990                                      *
*     Updated    :  July-August 2004:                                  *
*                   Renamed qt.f                                       *
*                   Updated for latest MRST, CTEQ parton distributions *
*                   Changed alpha 1/137 -> alpha(M_Z)                  *
*                   Updated top quark mass                             *
*                   September 2011:                                    *
*                   Corrected xlu,xld in setcon                        *
*                   Corrected triangle loop code in plumin             *
*----------------------------------------------------------------------*
*     Example input file "qt.inp"                                      *
*       computes dsigma/dQ_T**2 in p-pbar collisions at 1.8 TeV        *
*       for W+ + W- production using MRST2002 NLO parton densities     *
*       factorization/renormalization scales are set to Q_T            *
*       Q_T values: 10,30,50,70,90,110,130,150,170,190 GeV             *
*---- Begin: Cut here and remove leading and trailing * ---------------*
*     Input parameters for qt.f: Inclusive gamma*, W, Z production     *
*     Hadrons:    lppb    lns    e       pdset                         *
*                 t       f      1800    mrs02nlo                      *
*     Bosons:     npwz    lwpm   q     qt    y                         *
*                 2       t      80    50    0                         *
*     Partons:    lallpr  1  2  3  4  5  6  7  8  9  10 11 12          *
*                 t       f  t  t  f  f  f  f  f  f  f  f  f           *
*     QCD:        nscale  nalpha  lord1 lord2                          *
*                 1       1       t     t                              *
*     Scales:     lscaleq  lscalm  lscalmu  scalm  scalmu              *
*                 t        f       f        0.0    0.0                 *
*     Thresholds: linctri                                              *
*                 f                                                    *
*     Integrate:  ndim    ntot   begin   step   npoints                *
*                 3       0      10.0    20.0   10                     *
*     Vegas:      ncall   itmx0   itmx   ndev   nprn   alph0  alph     *
*                 40000   3       7      2      -1     1.5    1.5      *
*     Ran Seed:   iseed                                                *
*                 13579                                                *
*     Output:     lterm  lauto  ofname  comment                        *
*                 t      f      qt.out  #                              *
*---- End: Example input file -----------------------------------------*
*     Input file parameters:                                           *
* ------ line 3 ------                                                 *
*     lppb: t = proton proton, f = proton antiproton                   *
*     lns:  t = non-singlet,   f = non-singlet and singlet             *
*     e = hadron-hadron center of mass energy in GeV                   *
*     pdset: parton density set (8 character string)                   *
*         See BLOCK DATA pddata for sets implemented                   *
* ----- line 5 ------                                                  *
*     npwz: vector boson V produced                                    *
*         1 = gamma*, 2 = W-, 3 = W+, 4 = Z                            *
*     lwpm: t = sum W+ and W-  f = single W- (npwz=2) or W+ (npwz=3)   *
*     q  = virtual photon mass in GeV [program sets q=M_(W/Z) for W/Z] *
*     qt = vector boson transverse momentum in GeV                     *
*     y  = vector boson rapidity                                       *
* ------ line 7 ------                                                 *
*     lallpr: t = include all subprocesses                             *
*             f = include subprocesses marked t                        *
*     nproc:  1  2  3  4  5  6  7  8  9  10 11 12 in function xcfin    *
* ------ line 9 ------                                                 *
*     linctri:  t = include virtual triangle diagram, f = exclude      *
* ------ line 11 ------                                                *
*     nscale: energy scale used for renormalization/factorization      *
*         1 = q_transverse,  2 = Q=M_V,  3 = sqrt(q_t**2 + Q**2)       *
*     nalpha: use of NLO and LO QCD coupling alpha_s                   *
*         1 = NLO + NLO**2, 2 = NLO + LO*NLO, 3 = NLO + LO**2          *
*     lord1:  t = include,  f = exclude leading order                  *
*     lord2:  t = include,  f = exclude next-leading order             *
* ------ line 13 ------                                                *
*     lscaleq: t = set renormalization factorization scales equal      *
*     lscalm: t = cross section as function of factorization scale     *
*     lscalmu: t = cross section as function of renormalization scale  *
*     scalm: factorization scale is multiplied by sqrt(10**scalm)      *
*     scalmu: renormalization scale is multiplied by sqrt(10**scalmu)  *
* ------ line 15 ------                                                *
*     ndim: dimension of phase space integration                       *
*         2 = dsigma/dq_t**2/dy                                        *
*         3 = dsigma/dq_t**2                                           *
*         4 = sigma_total with (q_t > q_t_min)                         *
*     ntot: normalize by total cross section                           *
*         1 = compute and divide by total cross section                *
*     begin: starting value for chosen variable                        *
*     step:  step in chosen variable                                   *
*     npoints: number of points in chosen variable                     *
* ------ line 17 ------                                                *
*     ncall = number of integrand evaluations                          *
*     itmx0 = maximum number of warmup (discarded) iterations          *
*     itmx  = maximum number of iterations                             *
*     ndev  = Fortran device for output from Vegas                     *
*     nprn  : Vegas prints following on ndev for each iteration        *
*           > 0 : integral, std dev, chi^2, grid information           *
*           = 0 : integral, std dev, chi^2                             *
*           < 0 : nothing                                              *
*     alph0 = rate at which grid is modified during warmup             *
*     alph  = rate at which grid is modified during data taking        *
* ------ line 19 ------                                                *
*     iseed  = random number seed for ran2                             *
* ------ line 21 ------                                                *
*     lterm:  t = write numerical values on term, f = don't write      *
*     lauto:  t = derive output file name from input file name         *
*                 by replacing .inp with .out                          *
*             f = use ofname for output file                           *
*     ofname: name of output file                                      *
*     comment: comment character for output file                       *
************************************************************************
C-----------------------------------------------------------------------
      BLOCK DATA pddata
      IMPLICIT NONE
      INTEGER maxset,i
      PARAMETER (maxset=13)
      CHARACTER pdset(maxset)*8,pdname(maxset)*45
      COMMON /pdsets/ pdset,pdname
      DOUBLE PRECISION lambda(maxset)
      COMMON /lambdas/ lambda
      DATA (pdset(i),i=1,maxset)/
     & 'do84set1',
     & 'do84set2',
     & 'ehlqset1',
     & 'ehlqset2',
     & 'mrs123m1',
     & 'mrs123m2',
     & 'mrs123m3',
     & 'mrs89em1',
     & 'mrs89bm2',
     & 'mrs01lo1',     
     & 'mrs02nlo',     
     & 'mrs02nnl',     
     & 'ctq61nlo'
     & /
      DATA (pdname(i),i=1,maxset)/
     & 'Duke-Owens 1984 Set 1 Lambda=0.2 GeV         ',
     & 'Duke-Owens 1984 Set 2 Lambda=0.4 GeV         ',
     & 'EHLQ Set 1 Lambda=0.2 GeV                    ',
     & 'EHLQ Set 2 Lambda=0.29 GeV                   ',
     & 'MRS123 Mode 1 (soft glue) Lambda=0.107 GeV   ',
     & 'MRS123 Mode 2 (hard glue) Lambda=0.250 GeV   ',
     & 'MRS123 Mode 3 (1/RTX glue) Lambda=0.178 GeV  ',
     & 'MRSEB Mode 1 (EMC FIT) Lambda=0.091 GeV      ',
     & 'MRSEB Mode 2 (BCDMS FIT) Lambda=0.228 GeV    ',
     & 'MRST2001 Mode 1 LO Lambda_4=0.220 GeV        ',     
     & 'MRST2002 Mode 1 NLO Lambda_4=0.323 GeV       ',     
     & 'MRST2002 Mode 2 NNLO Lambda_4=0.234 GeV      ',     
     & 'Cteq61 Mode 1 NLO Lambda_4=0.326 GeV         '
     & /
      DATA (lambda(i),i=1,maxset)/
     & 0.200d0, 0.400d0, 0.200d0, 0.290d0,
     & 0.107d0, 0.250d0, 0.178d0, 0.091d0, 0.228d0,
     & 0.220d0, 0.323d0, 0.234d0, 0.326d0
     & /
      END
C-----------------------------------------------------------------------
      PROGRAM qtmain
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      CHARACTER today*30,proces*50,numb(12)*3,string*50,ifname*50
     &       ,ofname*50
      LOGICAL yes,no,lppb,lns,lwpm,lord1,lord2,lallpr,lproc(12)
     &       ,lscaleq,lscalm,lscalmu
      PARAMETER (yes=.true.,no=.false.)
      DIMENSION result(5,1000)
      COMMON /switch/ npwz,nscale,nalpha,nset,ndim,ntot
     &         ,lppb,lns,lwpm,lord1,lord2,lproc
      COMMON /const/ pi,tonbs,alpha,sin2tw,cf,ca,xnc,qw
     & ,xmw,xmz,xmc,xmb,xmt,qu,qd,xlu,xld,xru,xrd,v2(3,3)
      COMMON /variab/ e,q,qt,y,scalm,scalmu,qq,qt2,eqt,xlambd,escal2
     & ,sh,th,uh,aslo,asnlo,xnf,sqf2,sv2,sxrml,szf2,sv2d(3),sv2u(3)
     & ,x1,s2,s,t,u,d1,d2,d3,d4,d5,d6,d9,d10
     & ,fa,fm,fmu,fs,ft,fu,fs2,fst,fsu,ftu,fstu,fla,flt,flu
      LOGICAL lasnnlo
      COMMON /nnlo/ asnnlo,lasnnlo
      COMMON /bveg1/ ncall,itmx,nprn,ndev,xl(10),xu(10),acc
      COMMON /bveg3/ alph,ndmx,mds
      COMMON /counter/ iread
      COMMON /iquad/ alph0,itmx0
      DATA (numb(i),i=1,12) /' 1 ',' 2 ',' 3 ',' 4 ',' 5 ',' 6 ',
     & ' 7 ',' 8 ',' 9 ','10 ','11 ','12 '/
      INTEGER maxset,i
      PARAMETER (maxset=13)
      CHARACTER pdset(maxset)*8,pdname(maxset)*45
      COMMON /pdsets/ pdset,pdname
      DOUBLE PRECISION lambda(maxset)
      COMMON /lambdas/ lambda
      CHARACTER pdsetname*8,comment*2
      COMMON /pdsetn/ pdsetname
      LOGICAL lterm,lauto
      LOGICAL linctri
      COMMON /triangle/ linctri
 
      CALL ctime(time(),today)
*     use iargc and getarg to get input file name from command line
      IF (iargc().GT.0) THEN
          CALL getarg(1,ifname)
      ELSE
          ifname='qt.inp'
      ENDIF

      ninput=1
      OPEN (unit=ninput,file=ifname,status='old')
      PRINT *,'Reading input from file ',ifname
      READ (ninput,*)
      READ (ninput,*)
      READ (ninput,*) lppb,lns,e,string
      READ (ninput,*)
      READ (ninput,*) npwz,lwpm,q,qt,y
      READ (ninput,*)
      READ (ninput,*) lallpr,(lproc(i),i=1,12)
      READ (ninput,*)
      READ (ninput,*) nscale,nalpha,lord1,lord2
      READ (ninput,*)
      READ (ninput,*) lscaleq,lscalm,lscalmu,scalm,scalmu
      READ (ninput,*)
      READ (ninput,*) linctri
      READ (ninput,*)
      READ (ninput,*) ndim,ntot,begin,step,npoints
      READ (ninput,*)
      READ (ninput,*) ncall,itmx0,itmx,ndev,nprn,alph0,alph
      READ (ninput,*)
      READ (ninput,*) iseed
      READ (ninput,*)
      READ (ninput,*) lterm,lauto,ofname,comment

      nset=0
      DO i=1,maxset
          IF (string(1:8).EQ.pdset(i)) THEN
              nset=i
              xlambd=lambda(i)
              pdsetname=pdset(i)
          ENDIF
      ENDDO
      IF (pdsetname.EQ.'mrs02nnl') THEN
          lasnnlo=.TRUE.
      ELSE
          lasnnlo=.FALSE.
      ENDIF
      IF (nset.EQ.0) THEN
          PRINT *,'Bad parton density set ',string
          STOP
      ENDIF
      IF (lallpr) THEN
          DO i=1,12
              lproc(i)=yes
          ENDDO
      ENDIF
      scalm=10d0**scalm
      scalmu=10d0**scalmu
      idum=iseed
      CALL setseed(idum)
      IF (lauto) THEN
          i=index(ifname,'.inp')
          IF (i.EQ.0) i=index(ifname,' ')
          ofname=ifname(1:i-1)//'.out'
      ENDIF
      OPEN (unit=2,file=ofname,status='unknown')
      IF (lterm) PRINT *,'Output will go to file ',ofname
 
      IF (npwz.EQ.1) lproc(12)=no
      IF (npwz.EQ.2.OR.npwz.EQ.3) THEN
          lproc(10)=no
          lproc(11)=no
          lproc(12)=no
      ENDIF
      IF (lwpm.AND.npwz.EQ.3) npwz=2
      IF (lns) THEN
          lproc(2)=no
          lproc(3)=no
          lproc(5)=no
      ENDIF
      IF (lord1.AND..NOT.lord2) THEN
          DO i=3,12
              lproc(i)=no
          ENDDO
      ENDIF
      proces=' '
      DO i=1,12
          IF (lproc(13-i)) proces=(numb(13-i)//proces)
      ENDDO
      IF (lord2) proces=('nlo '//proces)
      IF (lord1) proces=('lo '//proces)
 
      CALL setcon
      IF (ntot.EQ.1) THEN
          CALL settot (totxc,totsd,totchi)
      ENDIF
      nres=0
      IF (lscalm.OR.lscalmu) THEN
          IF (lscalm.AND.lscalmu.AND..NOT.lscaleq) THEN
              DO i=1,npoints
                  scalm=10d0**(begin+step*(i-1))
                  DO j=1,npoints
                      scalmu=10d0**(begin+step*(j-1))
                      nres=nres+1
                      result(1,nres)=scalm
                      result(2,nres)=scalmu
                      CALL setup(result(3,nres),result(4,nres)
     &                          ,result(5,nres))
                      IF (ntot.EQ.1) THEN
                          result(3,nres)=result(3,nres)/totxc*2.
                          result(4,nres)=result(4,nres)/totxc*2.
                      ENDIF
                      IF (lterm) WRITE (6,105) (result(k,nres),k=1,5)
                  ENDDO
              ENDDO
          ELSE
              DO i=1,npoints
                  scale=10d0**(begin+step*(i-1))
                  nres=nres+1
                  result(1,nres)=scale
                  IF (lscalm.OR.lscaleq) scalm=scale
                  IF (lscalmu.OR.lscaleq) scalmu=scale
                  CALL setup (result(2,nres),result(3,nres)
     &                       ,result(4,nres))
                  IF (ntot.EQ.1) THEN
                      result(2,nres)=result(2,nres)/totxc*2.
                      result(3,nres)=result(3,nres)/totxc*2.
                  ENDIF
                  IF (lterm) WRITE (6,100) (result(j,nres),j=1,4)
              ENDDO
          ENDIF
      ELSE
          DO qtp=begin,begin+step*(npoints-1),step
              qt=qtp
              nres=nres+1
              result(1,nres)=qt
              CALL setup (result(2,nres),result(3,nres),result(4,nres))
              IF (ntot.EQ.1) THEN
                  result(2,nres)=result(2,nres)/totxc*2.
                  result(3,nres)=result(3,nres)/totxc*2.
              ENDIF
              IF (lterm) WRITE (6,100) (result(j,nres),j=1,4)
          ENDDO
      ENDIF

      WRITE (ndev,110) comment,today,ifname
      WRITE (ndev,120) comment
      IF (lppb)       string='p pbar collision'
      IF (.NOT. lppb) string='p p collision'
      IF (lns)        string='Non-singlet: ppb - pp'
      WRITE (ndev,130) comment,string,e
      WRITE (ndev,140) comment,pdname(nset)
      IF (npwz.EQ.1) string='Virtual photon'
      IF (npwz.EQ.2) THEN
          string='W^- Production'
          IF (lwpm)  string=' (W^-) + (W^+)'
      ENDIF
      IF (npwz.EQ.3) string='W^+ Production'
      IF (npwz.EQ.4) string='Z_0 Production'
      WRITE (ndev,150) comment,string,q
      WRITE (ndev,160) comment,proces
      WRITE (ndev,170) comment,nscale,nalpha
      IF (lscalm.OR.lscalmu) THEN
          IF (lscalm) WRITE (ndev,180) comment,qt,scalmu
          IF (lscalmu) WRITE (ndev,190) comment,qt,scalm
      ELSE
          WRITE (ndev,200) comment,scalm,scalmu
      ENDIF
      IF (npwz.EQ.4) THEN
          IF (linctri) THEN
              WRITE (ndev,205) comment,' included'
          ELSE
              WRITE (ndev,205) comment,' not included'
          ENDIF
      ENDIF
      IF (ndim.EQ.2) string='d sigma/(dqt**2 dy)'
      IF (ndim.EQ.3) string='d sigma/dqt**2'
      IF (ndim.EQ.4) string='sigma(qt > qtmin)'
      WRITE (ndev,210) comment,string
      IF (ntot.EQ.1) WRITE (ndev,220) comment,totxc,totsd,totchi
      WRITE (ndev,230) comment,ncall,itmx0,itmx,alph0,alph
      WRITE (ndev,240) comment,iseed
      WRITE (ndev,250) comment
      IF (lscalm.AND.lscalmu.AND..NOT.lscaleq) THEN
          DO i=1,nres
              IF (i.GT.1.AND.MOD(i-1,npoints).EQ.0) WRITE(ndev,106)
              WRITE (ndev,105) (result(j,i),j=1,5)
          ENDDO
      ELSE
          DO i=1,nres
              WRITE (ndev,100) (result(j,i),j=1,4)
          ENDDO
      ENDIF

      CLOSE (UNIT=ninput)
      CLOSE (UNIT=2)

 100  FORMAT (1x,4(g9.3,3x))
 105  FORMAT (1x,5(g9.3,3x))
 106  FORMAT (1x)
 110  FORMAT (a2,a30,'Input file: ',a30)
 120  FORMAT (a2,'Electroweak Boson Production at Large Q_T')
 130  FORMAT (a2,a30,'sqrt(S) = ',f7.1,' GeV')
 140  FORMAT (a2,'Parton Densities: ',a45)
 150  FORMAT (a2,'Vector Boson: ',a20,'Q = ',f7.2,' GeV')
 160  FORMAT (a2,'Parton Processes: ',a50)
 170  FORMAT (a2,'QCD: nscale = ',i2,2x,'nalpha = ',i2)
 180  FORMAT (a2,'Scalm dependence at Q_T = ',f7.2,' scalmu = ',f7.2)
 190  FORMAT (a2,'Scalmu dependence at Q_T = ',f7.2,' scalm = ',f7.2)
 200  FORMAT (a2,'scalm = ',f7.2,' scalmu = ',f7.2)
 205  FORMAT (a2,'Z_0 virtual triangle diagrams',a14)
 210  FORMAT (a2,'Cross section computed: ',a30)
 220  FORMAT (a2,'Divided by totxc/2 = ',3(g9.3,3x))
 230  FORMAT (a2,'Vegas: ncall = ',i8,' itmx0 = ',i4,' itmx = ',i4,
     &        '  alph0 = ',f4.2,'  alph = ',f4.2) 
 240  FORMAT (a2,'Seed for random number generator = ',i12)
 250  FORMAT (a2,'Variable   Integral    Std. Dev.   ChiSqd/dof')
      END
C-----------------------------------------------------------------------
      SUBROUTINE setcon
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      COMMON /const/ pi,tonbs,alpha,sin2tw,cf,ca,xnc,qw
     & ,xmw,xmz,xmc,xmb,xmt,qu,qd,xlu,xld,xru,xrd,v2(3,3)
 
      pi=3.141592654
      tonbs=.389386d6
*      alpha=1./137.0360
      alpha=1./127.918
      sin2tw=0.2312
      sw=dsqrt(sin2tw)
      cw=dsqrt(1.d0-sin2tw)
      cf=4./3.
      ca=3.
      xnc=3.
      xmz=91.1876
*      xmw=xmz*cw
      xmw=80.425
      xmc=1.25
      xmb=4.25
*      xmt=178.1
      xmt=1e6
      qu=2./3.
      qd=-1./3.
      qw=1./dsqrt(2.d0)/sw
      xlu=1./2./sw/cw-qu*sw/cw
      xld=-1./2./sw/cw-qd*sw/cw
      xru=-qu*sw/cw
      xrd=-qd*sw/cw
      v2(1,1)=(0.9738)**2
      v2(1,2)=(0.2200)**2
      v2(1,3)=1.-v2(1,1)-v2(1,2)
      v2(2,1)=(0.224)**2
      v2(2,2)=(0.996)**2
      v2(2,3)=1.-v2(2,1)-v2(2,2)
      v2(3,1)=1.-v2(1,1)-v2(2,1)
      v2(3,2)=1.-v2(1,2)-v2(2,2)
      v2(3,3)=1.-v2(1,3)-v2(2,3)
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE setup (ans,sd,chi2a)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      LOGICAL lppb,lns,lwpm,lord1,lord2,lallpr,lproc(12)
      COMMON /switch/ npwz,nscale,nalpha,nset,ndim,ntot
     &         ,lppb,lns,lwpm,lord1,lord2,lproc
      COMMON /const/ pi,tonbs,alpha,sin2tw,cf,ca,xnc,qw
     & ,xmw,xmz,xmc,xmb,xmt,qu,qd,xlu,xld,xru,xrd,v2(3,3)
      COMMON /variab/ e,q,qt,y,scalm,scalmu,qq,qt2,eqt,xlambd,escal2
     & ,sh,th,uh,aslo,asnlo,xnf,sqf2,sv2,sxrml,szf2,sv2d(3),sv2u(3)
     & ,x1,s2,s,t,u,d1,d2,d3,d4,d5,d6,d9,d10
     & ,fa,fm,fmu,fs,ft,fu,fs2,fst,fsu,ftu,fstu,fla,flt,flu
      LOGICAL lasnnlo
      COMMON /nnlo/ asnnlo,lasnnlo
      COMMON /bveg1/ ncall,itmx,nprn,ndev,xl(10),xu(10),acc
      COMMON /bveg3/ alph,ndmx,mds
      COMMON /iquad/ alph0,itmx0
      EXTERNAL f2veg,f2dim,f3veg,f3dim,s2lims,s2llim,s2ulim
      EXTERNAL f4dim,f4veg
 
      sh=e**2
      IF (npwz.EQ.2.OR.npwz.EQ.3) q=xmw
      IF (npwz.EQ.4) q=xmz
      qq=q**2
      qt2=qt**2
      eqt=dsqrt(qt2+qq)
      IF (nscale.EQ.1) escal2=qt2
      IF (nscale.EQ.2) escal2=qq
      IF (nscale.EQ.3) escal2=qt2+qq
      fm=dlog(scalm*escal2/qq)
      fmu=dlog(scalmu*escal2/qq)
      CALL flavor
      IF (nset.LE.9) THEN
          aslo=alphas_old(1)/2./pi
      ELSE
          scale=dsqrt(escal2*scalmu)
          aslo=alphas(scale,xlambd,0)/2d0/pi
      ENDIF
      asnlo=aslo
      asnnlo=aslo
      IF (lord2) THEN
          IF (nset.LE.9) THEN
              asnlo=alphas_old(2)/2./pi
          ELSE
              scale=dsqrt(escal2*scalmu)
              asnlo=alphas(scale,xlambd,1)/2d0/pi
              asnnlo=alphas(scale,xlambd,2)/2d0/pi
          ENDIF
      ENDIF
      surd=dsqrt((sh-qq)**2-4.*sh*qt2)
      ymax=1./2.*dlog((sh+qq+surd)/(sh+qq-surd))
      qtmax=(sh-qq)/2d0/e
      qtmin=qt
      xl(1)=0d0
      xu(1)=sh-qq
      xl(2)=dlog(qq/sh)
      xu(2)=0d0
      xl(3)=-ymax
      xu(3)=ymax
      xl(4)=dlog(qtmin**2/sh)
      xu(4)=dlog(qtmax**2/sh)
      npr=nprn
      itm=itmx
      alp=alph
      nprn=-1
      itmx=itmx0
      alph=alph0
      IF (ndim.EQ.2) CALL vegas (2,f2veg,ans,sd,chi2a)
      IF (ndim.EQ.3) CALL vegas (3,f3veg,ans,sd,chi2a)
      IF (ndim.EQ.4) CALL vegas (4,f4veg,ans,sd,chi2a)
      nprn=npr
      itmx=itm
      alph=alp
      IF (ndim.EQ.2) CALL vegas1 (2,f2veg,ans,sd,chi2a)
      IF (ndim.EQ.3) CALL vegas (3,f3veg,ans,sd,chi2a)
      IF (ndim.EQ.4) CALL vegas1 (4,f4veg,ans,sd,chi2a)
      omult=4.*pi*alpha*cf*xnc*tonbs
      ans=ans*omult
      sd=sd*omult
      IF (npwz.EQ.1) THEN
          emult=alpha/3./pi/qq
          ans=ans*emult
          sd=sd*emult
      ENDIF
      ans=ans*pi
      sd=sd*pi
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE settot (ans,sd,chi2a)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      LOGICAL lppb,lns,lwpm,lord1,lord2,lallpr,lproc(12)
      COMMON /switch/ npwz,nscale,nalpha,nset,ndim,ntot
     &         ,lppb,lns,lwpm,lord1,lord2,lproc
      COMMON /const/ pi,tonbs,alpha,sin2tw,cf,ca,xnc,qw
     & ,xmw,xmz,xmc,xmb,xmt,qu,qd,xlu,xld,xru,xrd,v2(3,3)
      COMMON /variab/ e,q,qt,y,scalm,scalmu,qq,qt2,eqt,xlambd,escal2
     & ,sh,th,uh,aslo,asnlo,xnf,sqf2,sv2,sxrml,szf2,sv2d(3),sv2u(3)
     & ,x1,s2,s,t,u,d1,d2,d3,d4,d5,d6,d9,d10
     & ,fa,fm,fmu,fs,ft,fu,fs2,fst,fsu,ftu,fstu,fla,flt,flu
      COMMON /bveg1/ ncall,itmx,nprn,ndev,xl(10),xu(10),acc
      COMMON /bveg3/ alph,ndmx,mds
      COMMON /iquad/ alph0,itmx0
      EXTERNAL totveg
 
      sh=e**2
      IF (npwz.EQ.2.OR.npwz.EQ.3) q=xmw
      IF (npwz.EQ.4) q=xmz
      qq=q**2
      escal2=qq
      CALL flavor
      saveit=scalmu
      scalmu=1.
      IF (nset.LE.9) THEN
          aslo=alphas_old(1)/2./pi
      ELSE
          aslo=alphas(q,xlambd,1)/2d0/pi
      ENDIF
      scalmu=saveit
      ymax=1./2.*dlog(sh/qq)
      xl(1)=0d0
      xu(1)=1.
      xl(2)=-ymax
      xu(2)=ymax
      npr=nprn
      itm=itmx
      alp=alph
      nprn=-1
      itmx=itmx0
      alph=alph0
      CALL vegas (2,totveg,ans,sd,chi2a)
      nprn=npr
      itmx=itm
      alph=alp
      CALL vegas1 (2,totveg,ans,sd,chi2a)
      omult=4.*pi**2*alpha*xnc*tonbs
      ans=ans*omult
      sd=sd*omult
      IF (npwz.EQ.1) THEN
          emult=alpha/3./pi/qq
          ans=ans*emult
          sd=sd*emult
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE flavor
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      LOGICAL lppb,lns,lwpm,lord1,lord2,lallpr,lproc(12)
      COMMON /switch/ npwz,nscale,nalpha,nset,ndim,ntot
     &         ,lppb,lns,lwpm,lord1,lord2,lproc
      COMMON /const/ pi,tonbs,alpha,sin2tw,cf,ca,xnc,qw
     & ,xmw,xmz,xmc,xmb,xmt,qu,qd,xlu,xld,xru,xrd,v2(3,3)
      COMMON /variab/ e,q,qt,y,scalm,scalmu,qq,qt2,eqt,xlambd,escal2
     & ,sh,th,uh,aslo,asnlo,xnf,sqf2,sv2,sxrml,szf2,sv2d(3),sv2u(3)
     & ,x1,s2,s,t,u,d1,d2,d3,d4,d5,d6,d9,d10
     & ,fa,fm,fmu,fs,ft,fu,fs2,fst,fsu,ftu,fstu,fla,flt,flu
      CHARACTER pdsetname*8
      COMMON /pdsetn/ pdsetname

      qsq=escal2
      xnf=3.
      sqf2=qu**2+2.*qd**2
      sxrml=xrd-xld
      DO 1 i=1,3
          sv2u(i)=v2(i,1)+v2(i,2)
          sv2d(i)=v2(1,i)
1     CONTINUE
      sv2=sv2u(1)
      szf2=(xlu**2+xru**2)/2.+xld**2+xrd**2
      IF (qsq.GT.4*xmc**2) THEN
          xnf=4.
          sqf2=sqf2+qu**2
          DO 2 i=1,3
              sv2d(i)=sv2d(i)+v2(2,i)
2         CONTINUE
          sv2=sv2u(1)+sv2u(2)
          szf2=szf2+(xlu**2+xru**2)/2.
          sxrml=0d0
      ENDIF
      IF (pdsetname(1:4).EQ.'do84') RETURN
      IF (qsq.GT.4*xmb**2) THEN
          xnf=5.
          sqf2=sqf2+qd**2
          DO 3 i=1,3
              sv2u(i)=sv2u(i)+v2(i,3)
3         CONTINUE
          sv2=sv2u(1)+sv2u(2)
          szf2=szf2+(xld**2+xrd**2)/2.
          sxrml=xrd-xld
      ENDIF
      IF (pdsetname(1:4).NE.'ehlq') RETURN
      IF (qsq.GT.4*xmt**2) THEN
          xnf=6.
          sqf2=sqf2+qu**2
          DO 4 i=1,3
              sv2d(i)=sv2d(i)+v2(3,i)
4         CONTINUE
          sv2=sv2u(1)+sv2u(2)+sv2u(3)
          szf2=szf2+(xlu**2+xru**2)/2.
          sxrml=0d0
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
C     QCD Effective coupling - Bardeen et al. definition
C
      DOUBLE PRECISION FUNCTION alphas_old (norder)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      LOGICAL lppb,lns,lwpm,lord1,lord2,lproc(12)
      COMMON /const/ pi,tonbs,alpha,sin2tw,cf,ca,xnc,qw
     & ,xmw,xmz,xmc,xmb,xmt,qu,qd,xlu,xld,xru,xrd,v2(3,3)
      COMMON /variab/ e,q,qt,y,scalm,scalmu,qq,qt2,eqt,xlambd,escal2
     & ,sh,th,uh,aslo,asnlo,xnf,sqf2,sv2,sxrml,szf2,sv2d(3),sv2u(3)
     & ,x1,s2,s,t,u,d1,d2,d3,d4,d5,d6,d9,d10
     & ,fa,fm,fmu,fs,ft,fu,fs2,fst,fsu,ftu,fstu,fla,flt,flu
      qsq=escal2*scalmu
      t=dlog(qsq/xlambd**2)
      b1=12.*pi/(33.-2.*xnf)
      alphas=b1/t
      IF (norder.EQ.2) THEN
          b2=24.*pi**2/(153.-19.*xnf)
          alphas=alphas*(1.-b1**2*dlog(t)/b2/t)
      ENDIF
      alphas_old=alphas
      RETURN
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION f2dim (xlogx1,s2p)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      COMMON /variab/ e,q,qt,y,scalm,scalmu,qq,qt2,eqt,xlambd,escal2
     & ,sh,th,uh,aslo,asnlo,xnf,sqf2,sv2,sxrml,szf2,sv2d(3),sv2u(3)
     & ,x1,s2,s,t,u,d1,d2,d3,d4,d5,d6,d9,d10
     & ,fa,fm,fmu,fs,ft,fu,fs2,fst,fsu,ftu,fstu,fla,flt,flu
      x1=dexp(xlogx1)
      s2=s2p
      f2dim=preigd ()
      RETURN
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION f3dim (yp,xlogx1,s2p)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      COMMON /variab/ e,q,qt,y,scalm,scalmu,qq,qt2,eqt,xlambd,escal2
     & ,sh,th,uh,aslo,asnlo,xnf,sqf2,sv2,sxrml,szf2,sv2d(3),sv2u(3)
     & ,x1,s2,s,t,u,d1,d2,d3,d4,d5,d6,d9,d10
     & ,fa,fm,fmu,fs,ft,fu,fs2,fst,fsu,ftu,fstu,fla,flt,flu
      y=yp
      x1=dexp(xlogx1)
      s2=s2p
      f3dim=preigd ()
      RETURN
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION f4dim (xlqt2,yp,xlogx1,s2p)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      COMMON /variab/ e,q,qt,y,scalm,scalmu,qq,qt2,eqt,xlambd,escal2
     & ,sh,th,uh,aslo,asnlo,xnf,sqf2,sv2,sxrml,szf2,sv2d(3),sv2u(3)
     & ,x1,s2,s,t,u,d1,d2,d3,d4,d5,d6,d9,d10
     & ,fa,fm,fmu,fs,ft,fu,fs2,fst,fsu,ftu,fstu,fla,flt,flu
      qt=e*dexp(xlqt2/2.)
      qt2=qt*qt
      eqt=dsqrt(qt2+qq)
      y=yp
      x1=dexp(xlogx1)
      s2=s2p
      f4dim=qt2*preigd ()
      RETURN
      END
C----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION f2tot (y,z)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      f2tot=pretot(y,z)
      RETURN
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION f2veg (x,wgt)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION x(2)
      f2veg=f2dim(x(2),x(1))
      RETURN
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION f3veg (x,wgt)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION x(3)
      f3veg=f3dim(x(3),x(2),x(1))
      RETURN
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION f4veg (x,wgt)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION x(4)
      f4veg=f4dim(x(4),x(3),x(2),x(1))
      RETURN
      END
C----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION totveg (x,wgt)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION x(2)
      totveg=f2tot(x(2),x(1))
      RETURN
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION s2lims (xlogx1)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      COMMON /variab/ e,q,qt,y,scalm,scalmu,qq,qt2,eqt,xlambd,escal2
     & ,sh,th,uh,aslo,asnlo,xnf,sqf2,sv2,sxrml,szf2,sv2d(3),sv2u(3)
     & ,x1,s2,s,t,u,d1,d2,d3,d4,d5,d6,d9,d10
     & ,fa,fm,fmu,fs,ft,fu,fs2,fst,fsu,ftu,fstu,fla,flt,flu
 
      ENTRY s2ulim (xlogx1)
      s2ulim = uh+dexp(xlogx1)*(sh+th-qq)
      RETURN
      ENTRY s2llim (xlogx1)
      s2llim = 0d0
      RETURN
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION preigd ()
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      LOGICAL lppb,lns,lwpm,lord1,lord2,lallpr,lproc(12)
      DIMENSION pdd(13,2),pdd0(13,2)
      COMMON /switch/ npwz,nscale,nalpha,nset,ndim,ntot
     &         ,lppb,lns,lwpm,lord1,lord2,lproc
      COMMON /const/ pi,tonbs,alpha,sin2tw,cf,ca,xnc,qw
     & ,xmw,xmz,xmc,xmb,xmt,qu,qd,xlu,xld,xru,xrd,v2(3,3)
      COMMON /variab/ e,q,qt,y,scalm,scalmu,qq,qt2,eqt,xlambd,escal2
     & ,sh,th,uh,aslo,asnlo,xnf,sqf2,sv2,sxrml,szf2,sv2d(3),sv2u(3)
     & ,x1,s2,s,t,u,d1,d2,d3,d4,d5,d6,d9,d10
     & ,fa,fm,fmu,fs,ft,fu,fs2,fst,fsu,ftu,fstu,fla,flt,flu
      LOGICAL lasnnlo
      COMMON /nnlo/ asnnlo,lasnnlo
      LOGICAL linctri
      COMMON /triangle/ linctri
 
      th=dexp(y)
      uh=qq-e*eqt*th
      th=qq-e*eqt/th
 
C     Phase space boundaries for Vegas
      IF (ndim.EQ.4) THEN
          surd=dsqrt((sh-qq)**2-4.*sh*qt2)
          ymax=1./2.*dlog((sh+qq+surd)/(sh+qq-surd))
          IF(dabs(y).GT.ymax) GOTO 100
      ENDIF
      b=-uh/(sh+th-qq)
      IF (x1.LT.b) GOTO 100
      a=uh+x1*(sh+th-qq)
      IF (s2.GT.a) GOTO 100
 
      xjac=x1*sh+uh-qq
      x20=(-x1*th-qq*(1.-x1))/xjac
      x2=s2/xjac+x20
      fs2=dlog(s2/qq)
      fa=dlog(a/qq)
      fscale=escal2*scalm
      CALL compdd (x1,x2,x20,fscale,pdd,pdd0)
      alow=0d0
      ads2=0d0
      as2a=0d0
      afin=0d0
      CALL comvar (0,x20)
      check=(x1+x20)/2.*e-eqt*dcosh(y)
      IF (check.LT.0.D0) STOP 'PREIGD error: Insufficient Energy'
      DO 2 i=1,2
          IF (lord1) THEN
              IF (lproc(1).AND.i.EQ.1) alow=alow+xclow(1)*pdd0(1,i)
              IF (lproc(2)) alow=alow+xclow(2)*pdd0(2,i)
          ENDIF
          IF (lord2) THEN
              IF (lproc(1).AND.i.EQ.1) THEN
                  ads2=ads2+xcds2(1)*pdd0(1,i)
                  IF (linctri) THEN
                      ads2=ads2+triang(1)*pdd0(12,i)
                  ENDIF
                  as2a=as2a-xcs2a(1)*pdd0(1,i)
              ENDIF
              IF (lproc(2)) THEN
                  ads2=ads2+xcds2(2)*pdd0(2,i)
                  IF (linctri) THEN
                      ads2=ads2+triang(2)*pdd0(13,i)
                  ENDIF
                  as2a=as2a-xcs2a(2)*pdd0(2,i)
              ENDIF
          ENDIF
      CALL exchtu
2     CONTINUE
      alow=alow/s
      ads2=ads2/s
      as2a=as2a/s
      IF (.NOT.lord2) GOTO 5
      CALL comvar (1,x2)
      DO 3 i=1,2
          IF (lproc(1).AND.i.EQ.1) as2a=as2a+xcs2a(1)*pdd(1,i)/s
          IF (lproc(2)) as2a=as2a+xcs2a(2)*pdd(2,i)/s
          DO 4 j=1,12
              IF (lproc(j)) afin=afin+xcfin(j)*pdd(j,i)
4         CONTINUE
          CALL exchtu
3     CONTINUE
      afin=afin/s+as2a
5     CONTINUE
      IF (lasnnlo) THEN
          as23=asnnlo
      ELSE
          as23=asnlo
      ENDIF
      IF (nalpha.EQ.1) preigd=as23*alow/a+as23**2*(ads2/a+afin)
      IF (nalpha.EQ.2) preigd=as23*((alow+aslo*ads2)/a+aslo*afin)
      IF (nalpha.EQ.3) preigd=as23*alow/a+aslo**2*(ads2/a+afin)
      preigd=preigd*x1/xjac
      RETURN
100   preigd=0d0
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION pretot (yq,z)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      LOGICAL lppb,lns,lwpm,lord1,lord2,lallpr,lproc(12)
      DIMENSION pdd(13,2),pdd0(13,2),pdd1(13,2)
      COMMON /switch/ npwz,nscale,nalpha,nset,ndim,ntot
     &         ,lppb,lns,lwpm,lord1,lord2,lproc
      COMMON /const/ pi,tonbs,alpha,sin2tw,cf,ca,xnc,qw
     & ,xmw,xmz,xmc,xmb,xmt,qu,qd,xlu,xld,xru,xrd,v2(3,3)
      COMMON /variab/ e,q,qt,y,scalm,scalmu,qq,qt2,eqt,xlambd,escal2
     & ,sh,th,uh,aslo,asnlo,xnf,sqf2,sv2,sxrml,szf2,sv2d(3),sv2u(3)
     & ,x1,s2,s,t,u,d1,d2,d3,d4,d5,d6,d9,d10
     & ,fa,fm,fmu,fs,ft,fu,fs2,fst,fsu,ftu,fstu,fla,flt,flu
 
      tau=qq/sh
 
C     Phase space boundary for Vegas
      zmin=tau*dexp(2.*dabs(yq))
      IF (z.LT.zmin) GOTO 100
 
      x11=dsqrt(tau)
      x21=x11
      dexpyq=dexp(yq)
      x11=x11*dexpyq
      x21=x21/dexpyq
      x2=dsqrt(z)
      x1=x11/x2
      x2=x21/x2
      CALL compdd (x11,x21,x21,qq,pdd1,pdd0)
      CALL compdd (x1,x2,x2,qq,pdd,pdd0)
 
      pretot=0d0
      IF (lproc(1)) THEN
          qqb=1./(1.-zmin)*(1.+aslo*cf*(2.*pi**2/3.-8.))*pdd1(1,1)
     &        +aslo*cf*(4.*dlog(1.-z)/(1.-z)*((1.+z**2)/z
     &                 *pdd(1,1)-2.*pdd1(1,1))
     &                 -2.*(1.+z**2)/z*dlog(z)*pdd(1,1))
          pretot=pretot+qqb
      ENDIF
      IF (lproc(2)) THEN
          qg=aslo*cf*((z**2+(1.-z)**2)*(dlog((1.-z)**2/z)-3./2.)
     &                +2.-z**2/2.)*(pdd(2,1)+pdd(2,2))
          pretot=pretot+qg
      ENDIF
      pretot=pretot/sh
      RETURN
100   pretot=0d0
      END
C-----------------------------------------------------------------------
      SUBROUTINE comvar (ns2,x2)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      COMMON /variab/ e,q,qt,y,scalm,scalmu,qq,qt2,eqt,xlambd,escal2
     & ,sh,th,uh,aslo,asnlo,xnf,sqf2,sv2,sxrml,szf2,sv2d(3),sv2u(3)
     & ,x1,s2,s,t,u,d1,d2,d3,d4,d5,d6,d9,d10
     & ,fa,fm,fmu,fs,ft,fu,fs2,fst,fsu,ftu,fstu,fla,flt,flu
 
      s=x1*x2*sh
      t=x1*th+(1.-x1)*qq
      u=x2*uh+(1.-x2)*qq
      s2p=s2*ns2
      d1=1./(s2p-t)
      d2=1./(s2p-u)
      d3=1./(s+t-s2p)
      d4=1./(s+u-s2p)
      d5=1./(s+qq-s2p)
      d6=1./(u*t-s2p*qq)
      d10=1./((u+t)**2-4.*s2p*qq)
      d9=dsqrt(d10)
      fs=dlog(s/qq)
      ft=dlog(-t/qq)
      fu=dlog(-u/qq)
      fst=fs-2.*dlog(-1./d1/t)
      fsu=fs-2.*dlog(-1./d2/u)
      ftu=dlog(d1*d2/d6)
      fstu=dlog(s*qq*d1*d2)
      fla=dlog((d9+d5)/(d9-d5))
      flt=dlog(s*qq*(1./(s2p*(2.*qq-u)-qq*t)/d1)**2)
      flu=dlog(s*qq*(1./(s2p*(2.*qq-t)-qq*u)/d2)**2)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE exchtu
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      COMMON /variab/ e,q,qt,y,scalm,scalmu,qq,qt2,eqt,xlambd,escal2
     & ,sh,th,uh,aslo,asnlo,xnf,sqf2,sv2,sxrml,szf2,sv2d(3),sv2u(3)
     & ,x1,s2,s,t,u,d1,d2,d3,d4,d5,d6,d9,d10
     & ,fa,fm,fmu,fs,ft,fu,fs2,fst,fsu,ftu,fstu,fla,flt,flu
 
      CALL exch (t,u)
      CALL exch (d1,d2)
      CALL exch (d3,d4)
      CALL exch (ft,fu)
      CALL exch (fst,fsu)
      CALL exch (flt,flu)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE exch (a,b)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      temp=a
      a=b
      b=temp
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE compdd (x1,x2,x20,qsq,pdd,pdd0)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      LOGICAL lppb,lns,lwpm,lord1,lord2,lallpr,lproc(12)
      DIMENSION pdd(13,2), pdd0(13,2), a(-6:6), b(-6:6), b0(-6:6)
      DIMENSION qdd(13,2), qdd0(13,2)
      COMMON /const/ pi,tonbs,alpha,sin2tw,cf,ca,xnc,qw
     & ,xmw,xmz,xmc,xmb,xmt,qu,qd,xlu,xld,xru,xrd,v2(3,3)
      COMMON /switch/ npwz,nscale,nalpha,nset,ndim,ntot
     &         ,lppb,lns,lwpm,lord1,lord2,lproc
      COMMON /counter/ iread
      CHARACTER pdsetname*8
      COMMON /pdsetn/ pdsetname

      IF (pdsetname(1:4).EQ.'do84') THEN
          IF (pdsetname(8:8).EQ.'1') mode=1
          IF (pdsetname(8:8).EQ.'2') mode=2
          CALL dopden (mode,x1,qsq,a)
          CALL dopden (mode,x2,qsq,b)
          CALL dopden (mode,x20,qsq,b0)
      ELSEIF (pdsetname(1:4).EQ.'ehlq') THEN
          IF (pdsetname(8:8).EQ.'1') mode=1
          IF (pdsetname(8:8).EQ.'2') mode=2
          CALL ehlq (mode,x1,qsq,a)
          CALL ehlq (mode,x2,qsq,b)
          CALL ehlq (mode,x20,qsq,b0)
      ELSEIF (pdsetname(1:6).EQ.'mrs123') THEN
          IF (pdsetname(8:8).EQ.'1') mode=1
          IF (pdsetname(8:8).EQ.'2') mode=2
          IF (pdsetname(8:8).EQ.'3') mode=3
          CALL mrs (mode,x1,qsq,a)
          CALL mrs (mode,x2,qsq,b)
          CALL mrs (mode,x20,qsq,b0)
      ELSEIF (pdsetname(1:5).EQ.'mrs89') THEN
          IF (pdsetname(6:6).EQ.'e') mode=4
          IF (pdsetname(6:6).EQ.'b') mode=5
          CALL mrs (mode,x1,qsq,a)
          CALL mrs (mode,x2,qsq,b)
          CALL mrs (mode,x20,qsq,b0)
      ELSEIF (pdsetname(1:5).EQ.'mrs01') THEN
          IF (pdsetname(6:8).EQ.'lo1') mode=0
          CALL mrst (mode,x1,qsq,a)
          CALL mrst (mode,x2,qsq,b)
          CALL mrst (mode,x20,qsq,b0)
      ELSEIF (pdsetname(1:5).EQ.'mrs02') THEN
          IF (pdsetname(6:8).EQ.'nlo') mode=1
          IF (pdsetname(6:8).EQ.'nnl') mode=2
          CALL mrst (mode,x1,qsq,a)
          CALL mrst (mode,x2,qsq,b)
          CALL mrst (mode,x20,qsq,b0)
      ELSEIF (pdsetname(1:5).EQ.'ctq61') THEN
          IF (pdsetname(6:8).EQ.'nlo') mode=1
          CALL cteq (mode,x1,qsq,a)
          CALL cteq (mode,x2,qsq,b)
          CALL cteq (mode,x20,qsq,b0)
      ENDIF
C     Initial color averages
      DO 1 i=-6,6
          a(i) = a(i)/xnc
          b(i) = b(i)/xnc
          b0(i) = b0(i)/xnc
1     CONTINUE
      glue = xnc/(xnc**2-1d0)
      a(0) = a(0)*glue
      b(0) = b(0)*glue
      b0(0) = b0(0)*glue
C
      CALL plumin (lppb,lns,npwz,a,b,pdd)
      CALL plumin (lppb,lns,npwz,a,b0,pdd0)
      IF (lwpm.AND.npwz.EQ.2.AND..NOT.lns) THEN
          CALL plumin (lppb,lns,3,a,b,qdd)
          CALL plumin (lppb,lns,3,a,b0,qdd0)
          DO 2 i=1,13
          DO 2 j=1,2
              pdd(i,j)=pdd(i,j)+qdd(i,j)
2             pdd0(i,j)=pdd0(i,j)+qdd0(i,j)
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE plumin (lppb,lns,npwz,a,b,pdd)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      LOGICAL lppb,lns
      DIMENSION a(-6:6),b(-6:6),pdd(13,2),h1u(3),h1ub(3),h2u(3),h2ub(3)
     &         ,h1d(3),h1db(3),h2d(3),h2db(3)
      COMMON /const/ pi,tonbs,alpha,sin2tw,cf,ca,xnc,qw
     & ,xmw,xmz,xmc,xmb,xmt,qu,qd,xlu,xld,xru,xrd,v2(3,3)
      COMMON /variab/ e,q,qt,y,scalm,scalmu,qq,qt2,eqt,xlambd,escal2
     & ,sh,th,uh,aslo,asnlo,xnf,sqf2,sv2,sxrml,szf2,sv2d(3),sv2u(3)
     & ,x1,s2,s,t,u,d1,d2,d3,d4,d5,d6,d9,d10
     & ,fa,fm,fmu,fs,ft,fu,fs2,fst,fsu,ftu,fstu,fla,flt,flu
      EXTERNAL sum
 
      qu2=qu**2
      qd2=qd**2
      w2=qw**2/2.
      zu2=(xlu**2+xru**2)/2.
      zd2=(xld**2+xrd**2)/2.
 
      CALL vector (a(1),a(4),a(6),h1u)
      CALL vector (a(-1),a(-4),a(-6),h1ub)
      CALL vector (a(2),a(3),a(5),h1d)
      CALL vector (a(-2),a(-3),a(-5),h1db)
      IF (lppb) THEN
          CALL vector (b(1),b(4),b(6),h2ub)
          CALL vector (b(-1),b(-4),b(-6),h2u)
          CALL vector (b(2),b(3),b(5),h2db)
          CALL vector (b(-2),b(-3),b(-5),h2d)
      ELSE
          CALL vector (b(1),b(4),b(6),h2u)
          CALL vector (b(-1),b(-4),b(-6),h2ub)
          CALL vector (b(2),b(3),b(5),h2d)
          CALL vector (b(-2),b(-3),b(-5),h2db)
      ENDIF
      sh1u=sum(h1u)+sum(h1ub)
      sh2u=sum(h2u)+sum(h2ub)
      sh1d=sum(h1d)+sum(h1db)
      sh2d=sum(h2d)+sum(h2db)
      sh1=sh1u+sh1d
      sh2=sh2u+sh2d
 
      IF (lns) THEN
          pdd(2,1)=0d0
          pdd(3,1)=0d0
          pdd(5,1)=0d0
          pdd(13,1)=0d0
          uv1=a(1)-a(-1)
          dv1=a(2)-a(-2)
          uv2=b(1)-b(-1)
          dv2=b(2)-b(-2)
          IF (npwz.EQ.1) THEN
              pdd(1,1)=qu2*uv1*uv2+qd2*dv1*dv2
              pdd(4,1)=sqf2*(uv1*uv2+dv1*dv2)
              pdd(7,1)=pdd(1,1)
              pdd(10,1)=2.*(qu*uv1+qd*dv1)*(qu*uv2+qd*dv2)
              pdd(12,1)=0d0
          ELSEIF (npwz.EQ.2.OR.npwz.EQ.3) THEN
              pdd(1,1)=w2*v2(1,1)*(uv1*dv2+dv1*uv2)
              pdd(4,1)=w2*sv2*(uv1*uv2+dv1*dv2)
              pdd(7,1)=w2*(sv2u(1)*dv1*dv2+sv2d(1)*uv1*uv2)
              pdd(10,1)=0d0
              pdd(12,1)=0d0
          ELSEIF (npwz.EQ.4) THEN
              pdd(1,1)=zu2*uv1*uv2+zd2*dv1*dv2
              pdd(4,1)=szf2*(uv1*uv2+dv1*dv2)
              pdd(7,1)=pdd(1,1)
              pdd(10,1)=((xlu+xru)*uv1+(xld+xrd)*dv1)/2.
              pdd(10,1)=pdd(10,1)*((xlu+xru)*uv2+(xld+xrd)*dv2)
              pdd(12,1)=(uv1*uv2+dv1*dv2)*(xru-xlu)*sxrml/2.
          ENDIF
          pdd(6,1)=pdd(1,1)
          pdd(8,1)=-pdd(1,1)/2.
          pdd(9,1)=-pdd(7,1)/2.
          pdd(11,1)=pdd(10,1)
          DO 55 i=1,13
55        pdd(i,2)=pdd(i,1)
          RETURN
      ENDIF
 
C     1. ann + qqb12
      IF (npwz.EQ.1) THEN
          pdd(1,1)=qu2*(dot(h1u,h2ub)+dot(h1ub,h2u))
     &            +qd2*(dot(h1d,h2db)+dot(h1db,h2d))
      ELSEIF (npwz.EQ.2) THEN
          pdd(1,1)=w2*(prod(h2ub,v2,h1d)+prod(h1ub,v2,h2d))
      ELSEIF (npwz.EQ.3) THEN
          pdd(1,1)=w2*(prod(h1u,v2,h2db)+prod(h2u,v2,h1db))
      ELSEIF (npwz.EQ.4) THEN
          pdd(1,1)=zu2*(dot(h1u,h2ub)+dot(h1ub,h2u))
     &            +zd2*(dot(h1d,h2db)+dot(h1db,h2d))
      ENDIF
      pdd(1,2)=pdd(1,1)
 
C     2. com
      IF (npwz.EQ.1) THEN
          pdd(2,1)=qu2*sh1u+qd2*sh1d
          pdd(2,2)=qu2*sh2u+qd2*sh2d
      ELSEIF (npwz.EQ.2) THEN
          pdd(2,1)=w2*(dot(sv2u,h1d)+dot(h1ub,sv2d))
          pdd(2,2)=w2*(dot(sv2u,h2d)+dot(h2ub,sv2d))
      ELSEIF (npwz.EQ.3) THEN
          pdd(2,1)=w2*(dot(h1u,sv2d)+dot(sv2u,h1db))
          pdd(2,2)=w2*(dot(h2u,sv2d)+dot(sv2u,h2db))
      ELSEIF (npwz.EQ.4) THEN
          pdd(2,1)=zu2*sh1u+zd2*sh1d
          pdd(2,2)=zu2*sh2u+zd2*sh2d
      ENDIF
      pdd(2,1)=pdd(2,1)*b(0)
      pdd(2,2)=a(0)*pdd(2,2)
 
C     3. fus
      IF (npwz.EQ.1) pdd(3,1)=sqf2
      IF (npwz.EQ.2.OR.npwz.EQ.3) pdd(3,1)=w2*sv2
      IF (npwz.EQ.4) pdd(3,1)=szf2
      pdd(3,1)=pdd(3,1)*a(0)*b(0)
      pdd(3,2)=pdd(3,1)
 
C     4. qqb34 squared
      dd=dot(h1u,h2ub)+dot(h1d,h2db)+dot(h1ub,h2u)+dot(h1db,h2d)
      IF (npwz.EQ.1) pdd(4,1)=sqf2*dd
      IF (npwz.EQ.2.OR.npwz.EQ.3) pdd(4,1)=w2*sv2*dd
      IF (npwz.EQ.4) pdd(4,1)=szf2*dd
      pdd(4,2)=pdd(4,1)
 
C     5. qqb56, qqb78, qq12, qq34, qq56, qq78, squared
      IF (npwz.EQ.1) THEN
          pdd(5,1)=qu2*sh1u+qd2*sh1d
          pdd(5,2)=qu2*sh2u+qd2*sh2d
      ELSEIF (npwz.EQ.2) THEN
          pdd(5,1)=w2*(dot(sv2u,h1d)+dot(h1ub,sv2d))
          pdd(5,2)=w2*(dot(sv2u,h2d)+dot(h2ub,sv2d))
      ELSEIF (npwz.EQ.3) THEN
          pdd(5,1)=w2*(dot(h1u,sv2d)+dot(sv2u,h1db))
          pdd(5,2)=w2*(dot(h2u,sv2d)+dot(sv2u,h2db))
      ELSEIF (npwz.EQ.4) THEN
          pdd(5,1)=zu2*sh1u+zd2*sh1d
          pdd(5,2)=zu2*sh2u+zd2*sh2d
      ENDIF
      pdd(5,1)=pdd(5,1)*sh2
      pdd(5,2)=sh1*pdd(5,2)
 
C     6. qqb1256 and qqb1278
      pdd(6,1)=pdd(1,1)
      pdd(6,2)=pdd(1,2)
 
C     7. qqb3456 and qqb3478
      IF (npwz.EQ.1.OR.npwz.EQ.4) THEN
          pdd(7,1)=pdd(1,1)
          pdd(7,2)=pdd(1,2)
      ELSEIF (npwz.EQ.2) THEN
          pdd(7,1)=w2*(dot3(sv2u,h1d,h2db)+dot3(h1ub,h2u,sv2d))
          pdd(7,2)=w2*(dot3(sv2u,h1db,h2d)+dot3(h1u,h2ub,sv2d))
      ELSEIF (npwz.EQ.3) THEN
          pdd(7,1)=w2*(dot3(sv2u,h1db,h2d)+dot3(h1u,h2ub,sv2d))
          pdd(7,2)=w2*(dot3(sv2u,h1d,h2db)+dot3(h1ub,h2u,sv2d))
      ENDIF
 
C     8. qq1256 and qq3478
      IF (npwz.EQ.1) THEN
          pdd(8,1)=(qu2*(dot(h1u,h2u)+dot(h1ub,h2ub))
     &             +qd2*(dot(h1d,h2d)+dot(h1db,h2db)))/2.
          pdd(8,2)=pdd(8,1)
      ELSEIF (npwz.EQ.2) THEN
          pdd(8,1)=w2*(prod(h2u,v2,h1d)+prod(h1ub,v2,h2db))/2.
          pdd(8,2)=w2*(prod(h1u,v2,h2d)+prod(h2ub,v2,h1db))/2.
      ELSEIF (npwz.EQ.3) THEN
          pdd(8,1)=w2*(prod(h1u,v2,h2d)+prod(h2ub,v2,h1db))/2.
          pdd(8,2)=w2*(prod(h2u,v2,h1d)+prod(h1ub,v2,h2db))/2.
      ELSEIF (npwz.EQ.4) THEN
          pdd(8,1)=(zu2*(dot(h1u,h2u)+dot(h1ub,h2ub))
     &             +zd2*(dot(h1d,h2d)+dot(h1db,h2db)))/2.
          pdd(8,2)=pdd(8,1)
      ENDIF
 
C     9. qq1278 and qq3456
      IF (npwz.EQ.1.OR.npwz.EQ.4) THEN
          pdd(9,1)=pdd(8,1)
      ELSEIF (npwz.EQ.2) THEN
          pdd(9,1)=w2*(dot3(sv2u,h1d,h2d)+dot3(h1ub,h2ub,sv2d))/2.
      ELSEIF (npwz.EQ.3) THEN
          pdd(9,1)=w2*(dot3(h1u,h2u,sv2d)+dot3(sv2u,h1db,h2db))/2.
      ENDIF
      pdd(9,2)=pdd(9,1)
 
C     10. qqb5678LL and qq1234LR
      IF (npwz.EQ.1) THEN
          tmp1=qu*(sum(h1u)-sum(h1ub))+qd*(sum(h1d)-sum(h1db))
          tmp2=qu*(sum(h2ub)-sum(h2u))+qd*(sum(h2db)-sum(h2d))
          pdd(10,1)=tmp1*tmp2
      ELSEIF (npwz.EQ.2.OR.npwz.EQ.3) THEN
          pdd(10,1)=0d0
      ELSEIF (npwz.EQ.4) THEN
          tmp1=xlu*sum(h1u)+xld*sum(h1d)-xru*sum(h1ub)-xrd*sum(h1db)
          tmp2=xlu*sum(h2ub)+xld*sum(h2db)-xru*sum(h2u)-xrd*sum(h2d)
          pdd(10,1)=tmp1*tmp2/2.
          tmp1=xlu*sum(h1ub)+xld*sum(h1db)-xru*sum(h1u)-xrd*sum(h1d)
          tmp2=xlu*sum(h2u)+xld*sum(h2d)-xru*sum(h2ub)-xrd*sum(h2db)
          pdd(10,1)=pdd(10,1)+tmp1*tmp2/2.
      ENDIF
      pdd(10,2)=pdd(10,1)
 
C     11. qqb5678LR and qq1234LL
      IF (npwz.LT.4) THEN
          pdd(11,1)=pdd(10,1)
      ELSEIF (npwz.EQ.4) THEN
          tmp1=xlu*sum(h1u)+xld*sum(h1d)-xru*sum(h1ub)-xrd*sum(h1db)
          tmp2=xru*sum(h2ub)+xrd*sum(h2db)-xlu*sum(h2u)-xld*sum(h2d)
          pdd(11,1)=tmp1*tmp2/2.
          tmp1=xlu*sum(h1ub)+xld*sum(h1db)-xru*sum(h1u)-xrd*sum(h1d)
          tmp2=xru*sum(h2u)+xrd*sum(h2d)-xlu*sum(h2ub)-xld*sum(h2db)
          pdd(11,1)=pdd(11,1)+tmp1*tmp2/2.
      ENDIF
      pdd(11,2)=pdd(11,1)
 
C     12. qqb1234
      IF (npwz.LT.4) THEN
          pdd(12,1)=0d0
      ELSEIF (npwz.EQ.4) THEN
          pdd(12,1)=dot(h1u,h2ub)+dot(h1ub,h2u)
          pdd(12,1)=pdd(12,1)-dot(h1d,h2db)-dot(h1db,h2d)
          pdd(12,1)=pdd(12,1)*(xru-xlu)*sxrml/2.
      ENDIF
      pdd(12,2)=pdd(12,1)
 
C     13. COM -- triangle loop contributions -- corrected 2011
      IF (npwz.LT.4) THEN
          pdd(13,1)=0d0
          pdd(13,2)=0d0
      ELSEIF (npwz.EQ.4) THEN
          pdd(13,1)=sum(h1u)+sum(h1ub)-sum(h1d)-sum(h1db)
          pdd(13,1)=pdd(13,1)*(xru-xlu)*sxrml*b(0)/2.
          pdd(13,2)=sum(h2u)+sum(h2ub)-sum(h2d)-sum(h2db)
          pdd(13,2)=a(0)*pdd(13,2)*(xru-xlu)*sxrml/2.
      ENDIF
 
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE vector (a,b,c,vec)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION vec(3)
      vec(1)=a
      vec(2)=b
      vec(3)=c
      RETURN
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION sum (a)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION a(3)
      sum=a(1)+a(2)+a(3)
      RETURN
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION dot (a,b)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION a(3),b(3)
      dot=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      RETURN
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION dot3 (a,b,c)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION a(3),b(3),c(3)
      dot3=a(1)*b(1)*c(1)+a(2)*b(2)*c(2)+a(3)*b(3)*c(3)
      RETURN
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION prod (a,b,c)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION a(3),b(3,3),c(3)
      prod=0d0
      DO 1 i=1,3
      DO 1 j=1,3
1     prod=prod+a(i)*b(i,j)*c(j)
      RETURN
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION xclow (nproc)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      COMMON /const/ pi,tonbs,alpha,sin2tw,cf,ca,xnc,qw
     & ,xmw,xmz,xmc,xmb,xmt,qu,qd,xlu,xld,xru,xrd,v2(3,3)
      COMMON /variab/ e,q,qt,y,scalm,scalmu,qq,qt2,eqt,xlambd,escal2
     & ,sh,th,uh,aslo,asnlo,xnf,sqf2,sv2,sxrml,szf2,sv2d(3),sv2u(3)
     & ,x1,s2,s,t,u,d1,d2,d3,d4,d5,d6,d9,d10
     & ,fa,fm,fmu,fs,ft,fu,fs2,fst,fsu,ftu,fstu,fla,flt,flu
 
      GOTO (1,2) nproc
C     q qbar -> v G
1     xclow=u/t+t/u+2.*s*qq/t/u
      RETURN
C     q G -> V q
2     xclow=-s/t-t/s-2.*u*qq/s/t
      RETURN
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION xcds2 (nproc)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      COMMON /const/ pi,tonbs,alpha,sin2tw,cf,ca,xnc,qw
     & ,xmw,xmz,xmc,xmb,xmt,qu,qd,xlu,xld,xru,xrd,v2(3,3)
      COMMON /variab/ e,q,qt,y,scalm,scalmu,qq,qt2,eqt,xlambd,escal2
     & ,sh,th,uh,aslo,asnlo,xnf,sqf2,sv2,sxrml,szf2,sv2d(3),sv2u(3)
     & ,x1,s2,s,t,u,d1,d2,d3,d4,d5,d6,d9,d10
     & ,fa,fm,fmu,fs,ft,fu,fs2,fst,fsu,ftu,fstu,fla,flt,flu
 
      sti=d3
      sui=d4
      tui=-d9
      GOTO (1,2) nproc
1     CONTINUE
C     q qbar -> v G
      t0ut = u/t+t/u+2.*qq*s/u/t
      realp =
     &  t0ut*(-(11./6.*ca-xnf/3.)*fa+67./18.*ca-5./9.*xnf
     & +ca*fa**2+2.*cf*(ft+fu-2.*fa)*fm
     & +(cf-ca/2.)*(pi**2./3.+(2.*fa+fs-fu-ft)**2)-3.*cf*fm)
      virtp =
     &  t0ut*(-8.*cf-cf*fs**2+(2.*cf/3.-ca/6.)*pi**2
     &    +ca/2.*(fs**2-(ft+fu)**2)+(11./6.*ca-xnf/3.)*fmu)
     & +cf*(s*(sti+sui)+(s+t)/u+(s+u)/t)
     & +ft*(cf*(4.*s**2+2.*s*t+4.*s*u+u*t)*sui**2+ca*t*sui)
     & +fu*(cf*(4.*s**2+2.*s*u+4.*s*t+t*u)*sti**2+ca*u*sti)
     & +(2.*cf-ca)*(2.*fs*(s**2*tui**2+2.*s*tui)
     &    -qq*(u**2+t**2)/u/t*tui)
     & -(2.*cf-ca)*((s**2+(s+u)**2.)/u/t*fr1(qq,s,t)
     &    +(s**2+(s+t)**2)/t/u*fr1(qq,s,u))
     & +ca*t0ut*fr2(qq,t,u)  + 0.0
      xcds2 = realp + virtp
      RETURN
2     CONTINUE
C     q G -> v G
      t0st = s/t+t/s+2.*qq*u/s/t
      realp =
     &  -t0st*(cf*(7./2.+2*fm*(fu-fa)+fa**2-3./2.*(fm+fa))
     &      +ca*(pi**2/6.+(fs-ft-fu)**2/2.+2.*ft*(fm-fa)
     &      +2.*fa*(fs-fu+fa-fm))-fm*(11./6.*ca-xnf/3.))
      virtp =
     &  t0st*(-8.*cf-cf*fu**2-1./3.*(cf-ca)*pi**2
     &    +ca/2.*(fu**2-fs**2-ft**2)+(11./6.*ca-xnf/3.)*fmu)
     & +cf*(u*(tui+sui)+(u+t)/s+(s+u)/t)
     & +ft*(cf*(4.*u**2+2.*u*t+4.*s*u+s*t)*sui**2+ca*t*sui)
     & +fs*(cf*(4.*u**2+2.*s*u+4.*u*t+t*s)*tui**2+ca*s*tui)
     & +(2.*cf-ca)*(2.*fu*(u**2*sti**2+2.*u*sti)
     &    -qq*(s**2+t**2)/s/t*sti)
     & -(2.*cf-ca)*((u**2+(u+s)**2)/s/t*(ft*fu-pi**2
     &    +pi**2/2.-fr2(qq,t,u))+(u**2+(t+u)**2)/s/t*fr1(qq,s,u))
     & -ca*t0st*fr1(qq,s,t)
      virtp = - virtp
      xcds2 = realp + virtp
      RETURN
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION triang (nproc)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      COMMON /const/ pi,tonbs,alpha,sin2tw,cf,ca,xnc,qw
     & ,xmw,xmz,xmc,xmb,xmt,qu,qd,xlu,xld,xru,xrd,v2(3,3)
      COMMON /variab/ e,q,qt,y,scalm,scalmu,qq,qt2,eqt,xlambd,escal2
     & ,sh,th,uh,aslo,asnlo,xnf,sqf2,sv2,sxrml,szf2,sv2d(3),sv2u(3)
     & ,x1,s2,s,t,u,d1,d2,d3,d4,d5,d6,d9,d10
     & ,fa,fm,fmu,fs,ft,fu,fs2,fst,fsu,ftu,fstu,fla,flt,flu
 
      GOTO (1,2) nproc
1     CONTINUE
C     q qbar -> v G
      triang = -(s+qq)/(s-qq)*(1.-qq/(s-qq)*fs)
      RETURN
 
2     CONTINUE
C     q G -> v G
      triang = (u+qq)/(u-qq)*(1.-qq/(u-qq)*fu)
      RETURN
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION fr1 (qq,s,t)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      fr1 = dlog(s/qq)*dlog(t/(qq-s))+0.5*dlog(qq/s)**2
     &   -0.5*dlog((qq-t)/qq)**2+dilog(qq/s)-dilog(qq/(qq-t))
      RETURN
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION fr2 (qq,t,u)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      fr2 = 0.5*dlog((qq-t)/qq)**2+0.5*dlog((qq-u)/qq)**2
     &       +dilog(qq/(qq-t))+dilog(qq/(qq-u))
      RETURN
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION xcs2a (nproc)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      COMMON /const/ pi,tonbs,alpha,sin2tw,cf,ca,xnc,qw
     & ,xmw,xmz,xmc,xmb,xmt,qu,qd,xlu,xld,xru,xrd,v2(3,3)
      COMMON /variab/ e,q,qt,y,scalm,scalmu,qq,qt2,eqt,xlambd,escal2
     & ,sh,th,uh,aslo,asnlo,xnf,sqf2,sv2,sxrml,szf2,sv2d(3),sv2u(3)
     & ,x1,s2,s,t,u,d1,d2,d3,d4,d5,d6,d9,d10
     & ,fa,fm,fmu,fs,ft,fu,fs2,fst,fsu,ftu,fstu,fla,flt,flu
 
      GOTO (1,2) nproc
1     CONTINUE
      t0ut = u/t+t/u+2.*qq*s/u/t
      anns2a =
     &  t0ut*(xnf/3.-11./6.*ca+2.*cf*(ftu-2.*fm)
     &    +(cf-ca/2.)*(4.*fstu-2.*ftu)
     &    +fs2*(8.*cf-2.*ca))
      xcs2a = anns2a/s2
      RETURN
2     CONTINUE
      t0st = s/t+t/s+2.*qq*u/s/t
      coms2a =
     &     -t0st*(cf*(-3./2.+fsu+2.*ftu-2.*fm+(t+u)*fla*d9)
     &            +ca*(2.*fstu+(fst-fsu)/2.-ftu-2.*fm))
     &     -fs2*t0st*(2.*cf+4.*ca)
     &     +(cf-ca/2)*((fla*(t+u)*d9+fsu)*(s+2.*u)/t
     &                +2.*ftu*(t+2.*u)/s)
      xcs2a = coms2a/s2
      RETURN
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION xcfin (nproc)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      COMMON /const/ pi,tonbs,alpha,sin2tw,cf,ca,xnc,qw
     & ,xmw,xmz,xmc,xmb,xmt,qu,qd,xlu,xld,xru,xrd,v2(3,3)
      COMMON /variab/ e,q,qt,y,scalm,scalmu,qq,qt2,eqt,xlambd,escal2
     & ,sh,th,uh,aslo,asnlo,xnf,sqf2,sv2,sxrml,szf2,sv2d(3),sv2u(3)
     & ,x1,s2,s,t,u,d1,d2,d3,d4,d5,d6,d9,d10
     & ,fa,fm,fmu,fs,ft,fu,fs2,fst,fsu,ftu,fstu,fla,flt,flu
 
      GOTO (1,2,3,4,5,6,7,8,9,10,11,12) nproc
 
1     CONTINUE
C     1. ann and qqb12
      xcfin =
     &  cf*(s*d1**2+2.*d1-s/u/t+(2.*qq*u+t*s2-2.*qq*s2)*d6/t
     &    -2.*qq*(t-u)*d1/t/u)
     & +ca*(-11./6.*s/u/t+s**2*(3.*s2-4.*t)*d1**2/2./t/u
     &    +2.*s*d1/u+qq/3./t**2)
     & +xnf/3.*(s/u/t-qq/t**2)
     & +2.*(cf-ca/2.)*(fstu+fs2-ftu)
     &    *(qq-u)**2*d1/u*(d2-1/t)
     & +2.*(cf-ca/2)*(fstu+fs2)*(s-qq)/t/u
     & +cf*ftu*(4.*qq*(qq-t)**2*d6/u/t+(2.*(qq+s)-s2)*d6
     &    +(s2-2.*s)/u/t)
     & -ca*ftu*qq/u/t
     & +cf*(fm-fs2)*(d6/u/t*(4.*u*t*(u-qq)-4.*qq*(u-qq)**2
     &    -u*t*s2)+(qq-u)*d1**2-(2.*qq-u)*d1/t-2./t+qq/t**2
     &    +4.*qq/u/t)  + 0.0
      RETURN
 
2     CONTINUE
C     2. com
      xcfin =
     &  fla*d9*d10*(d10*3./4.*cf*s*(s+qq-s2)*(u**2-t**2)*(1-u/t)
     &     +cf*((u**2-t**2)*((qq-s2)/4./s+11./4.+s2/t-u/2./t)
     &       -2.*s*(t+u**2/t)-s*(s-s2)*(u/t-3.)/2.+4.*s2*(t-u)+5.*s*u)
     &     +ca*(1.-u/t)*((t+u)*(s+qq-s2)/4.+s*(s-s2)))
     & +fla*d9*(cf*((t-u-2.*s2)/4./s+(11.*s-2.*u-4.*s2)/4./t
     &             +(2.*(1.+u/t)*(5.*qq+s2)-16.*qq*s2/t)/s)
     &        -ca/2.*((1.+u/t)/2.+5.*(u/s-(u+s2)/t)-7.*(s+u*s2/s)/t
     &             +3.*(t-3.*s2)/s+2.*(u**2+3.*s2**2)/s/t))
     & +d10*(cf*(3./2.*d10*s*(1.+u/t)*(t-u)**2+s*(u-3.*t)/2./t
     &    +(1.-u/t)*(s2-7.*u/4.)+(t-u)*(11./4.+(3./2.*(t+u)-2.*s2)/s))
     &      +ca*(u/t-1.)*(s+s2))
     & +2.*d1**3*s*t*(4.*cf-3.*ca-(fs2-fm)*(cf-ca))
     & +d1**2*(cf*(3.*s-8.*t-4.*u)-ca*(9.*s-4.*(t-u))
     &      +2.*(fs2-fm)*(cf*(t+u)+ca*(3.*s+2.*u)))
      xcfin=xcfin
     & +d1*((cf*(fs2-fm)-ca*(fs2-fm+1./2.*(fstu+fs2
     &    -ftu+(fst+flt)/2.)))*(4.*(1.+u/s)*(1.-u/t)-2.*(s/t+t/s))
     &    -cf*(fs2-fm)*(s+2.*u)/t
     &    +ca*((2.*flt+2.*fst)*(u+s)/t+2.*(fs2-fm)*((5.*s+4.*u)/t-2.))
     &    -cf*(4.+(s-2.*u)/t+(2.*u-3.*t)/s-u**2/s/t)
     &    +ca*(7.-(3.*s+2.*u)/t))
     & +d2*(2.*cf*(fm-fs2-s/t)+ca*(s/t-1.))
     & +d3*cf*(fsu-2.*ftu-flt)*(1.-s2/t)*(s2**2-2.*u*(s2-u))/s**2
     & +d6*cf*(2.*(fm-fs2-ftu)*(2.*s/t*(qq+2.*u)+t-s2+4.*u/t*(2.*u-s2)
     &      +2.*(u-s2)*(t-2.*u+2.*u**2/t)/s)
     &      +s2-t+4.*u*(qq/t+(1.-u/t)*(t-qq)/s))
      xcfin=xcfin
     & +cf*((fsu-2.*ftu-flt)*((2.*u*(s2-u)-s2**2)/s-s2)/s/t
     &      -2.*ftu*(1./t-2./s+4.*u/s/t)+flt*(u-3.*s2+5.*qq)/s/t
     &      +fsu*(5./s+(1.+2.*t/s-3.*s2/s)/u-s2*qq/s/u**2)
     &      +(fs2-fm)*(8./s+3./t+1./u-4.*u/s/t+(2.*t-3.*s2)/s/u
     &        -qq*(1./t**2+s2/s/u**2))-2./s+3./4./t-1./u+3.*u/s/t
     &      -(t+s2)/s/u+qq*(s2/s/u**2-1./2./t**2))
     & +ca*(2.*(fs2-fm)*(6./s-2./t-2.*qq/s/t-3.*s2/s/t)
     &      -2.*fst*(-15./4./s+9.*s2/4./s/t+qq/t**2)
     &      +2.*(fs2-fm+fst)*((s+s2**2/s)/t**2
     &      +2.*s2*qq/t**3-4.*qq/t*(1.-qq/t)/s
     &      +2.*s2/s/t*qq/t*(1.-s2/t)*(1.-qq/t))
     &      +(fs2+fstu-(fst+fsu)/2.)*(4.-2.*u/t-s2/t)/s+fsu*(1.-u/t)/s
     &      +ftu*(2.*u/t-1.)/s-2.*flt*(3./t+1./s+2.*(u-s2)/s/t)
     &      -1./s-1./t+2.*s2/s/t+(3.*qq+4.*s2*(1.-3.*qq/t))/t**2
     &      -2.*(s+s2**2./s)/t**2-12.*s2/s/t*qq/t*(1.-s2/t)*(1.-qq/t)
     &      +2./s*qq/t*(1.-qq/t))   +  0.0
      RETURN
 
3     CONTINUE
C     3. fus
      xcfin =
     & ca*d5*d9*fla*(d10*t**2*(3./2.*d10*(2.*s-t-u)*(u**2-t**2)
     &  +(2.*s-3.*t+u+(t**2-u**2)/s))-3.*s/2.-t*(9./2.+5.*t/s+3.*u/s))
     &+ca*d9*fla*(d10*t**2*(3./2.*d10*(u**2-t**2)+1.)*(qq-s2)/s
     &     -11.*(qq-s2)/s/4.-4.*(1.+t/s))
     &+cf*d9*fla*(2.*d5*(t+u-2.*s)+6.)
     &+ca*d10*t*(3.*d10*(t-u)**2*(1.+(t+u)/2./s-2.*s*d5)
     &     +d5*(2.*s-3.*t+3.*u)-1.-u/s+s2*(t-u)/s**2
     &     +(d5-1./2./s)*(t-u)**2/s)
     &+cf*4.*d6*((fs2-fm+ftu)*(s/2.+2.*t*(1.+t/s))-t*(1.+(t+u)/s))
     &+d5*2*(d1*(ca-2.*cf)*t**2*(fst+flt)/s
     &      +d3*cf*(fsu-flt-2.*ftu)*(s+2.*u*(1.+u/s)))
     &+d5*(2.*cf*(fst-(3.+4.*(t+u)/s)*flt-(4.+8.*t/s)*ftu)
     &    +ca/2.*((fst+flt)*(2.+(u+3.*t)/s)+1.+t*(t-u)/s**2))
      xcfin = xcfin
     &+2.*cf*d3*(2.*ftu+flt-fsu)*(2.+(t+u)/s)
     &+8.*t*d1**2*(cf*(1.-fs2+fm)+ca)
     &+4.*(cf-ca/2.)*d1*d2*(fs2-ftu+fstu)*t*u/s
     &+d1*((ca-2.*cf)*((2.+(t+u)/s)*(fst+flt)+2.*(u-t)/s
     &    *(fs2+fstu-ftu))-8.*cf*(fs2-fm-1.)+4.*ca*(2.+(u-t)/s))
     &+cf*(2.*(fs2-fm)*(2./s-2./t-(u+s2*qq/t)/s/t)
     &    +4.*(fs2+fstu)/s-2.*fst/t*(2.+(u+s2*qq/t)/s)
     &    +(4.*qq*(1.+s2/t)-2.*u)/s/t)
     &-ca/s*(2.*(fs2+fstu)+fst/2.+5./2.*flt+15./4.)  + 0.0
      RETURN
 
4     CONTINUE
C     4. qqb34
      xcfin =
     &  d5/2.*(u/s-1.)-5./4./s
     & +d10*(u*d5*(-2.*s+3./2.*u/s*(t-u)+4.*u-2.*t)
     &    +u/2./s*(2.*s+2.*s2+t-u))
     & +d10**2*3.*u**2*(u-t)*(d5*(2.*s-u-t)-(s+s2)/s)
     & +d9*fla*(d5*(2.*u**2/s+5./2.*u+3./2.*s)+(3./4.*s+u-s2/2.)/s)
     & +d9*d10*fla*(u**2*d5*(3.*u-t-u**2/s+t**2/s-2.*s)
     &    +u/s*(2.*s2*t-u*t-2.*s2**2+4.*u*s2-3.*u**2+2.*s*s2-u*s))
     & +d9*d10**2*fla*(3.*d5*(u**2-t**2)*u**2*(s-qq)
     &    +3.*u**2*qq/s*(u-t)*(u+t-2.*s2))  + 0.0
      xcfin=xcfin/2.
      RETURN
 
5     CONTINUE
C     5. qqb56
      xcfin =
     &  d9*d10*s*fla*(u**2-t**2)/t
     & +d9/t*fla*(3.*s+2.*u)
     & +d10*(t+u-2.*s2)*(1.-u/t)
     & +(fs2-fm)*d1*(s*d1+4.*s/t+2.*u/t)
     &  +(fst-fm+fs2)*((s+s2-qq)**2/s/t**2
     &    +(u/t-2.*s2*qq/t**2)**2/s+4.*(qq/t-1.)/t)
     &  +fst*2*(1.-qq/t)/t
     &  -d1*(s*d1-1.-u/t)+2.*(u/t+s2*(1./t+1./s-s2/t/s))/t
     &  +(qq/t-1.)/t+12.*s2*qq*(u-s2*qq/t)/s/t**3  + 0.0
      xcfin=xcfin/2.
      RETURN
 
6     CONTINUE
C     6. qqb1256
      xcfin =
     &  d9*fla*(d10**2*3.*s**2*qq*(t-u)**2*(1./t+1./u)
     &         +d10*s*qq*(5.*s/t-7.*s/u+4.*u/t-4.)
     &         +(2.*s**2/u+(qq*(2.*s+u)-u*s)/t)/s2+(2.*s2-u)/t-2.)
     & +d10**2*3.*s*qq*(2.*s2-u-t)*(t-u)**2/u/t
     & +d10*(s*(s-s2)*(7./u-5./t)+(t/u-u/t)*(7.*s+t+3.*u-2.*s2)/2.
     &         +3.*s*(t/u-u/t)/2.+2.*s-2.*s2*(1.-u/t))
     & +d1*s**2*(d1+3./t)/u+fst*((2.*s*(1.+s/u)+u)/s2-2.*qq/t)/t
     &         -1./2./u+3./2./t-3.*qq/t**2
      xcfin=xcfin*(cf-ca/2.)
      RETURN
 
7     CONTINUE
C     7. qqb3456
      xcfin =
     &  d9*fla*(6.*s*s2*qq*d10**2*(t-u)**2/t
     &       +d10*(2.*s*u+(s2*(1.+t/s)-t*(t+u)/2./s)
     &           *(t+u)*(1.-u/t)+3.*s/2.*(t-u)*(1.-u/t)
     &           -4.*s2*qq*(1.+(s-u)/t))
     &       +d5*(6.*s+2.*u+(t**2-u**2)/2./s)
     &       -1.+(s/2.+2.*u-s2)/t+(3.*t+2.*u-6.*s2)/s
     &           +(u**2-2*u*s2+4.*s2**2)/s/t)
     & +d10**2*3.*s2*(t-u)**2*(1.+u/t-2.*qq/t)
     & +d10*(-t/2.+3.*u+s2+(t**2-u**2)/s-u*(u+6.*s2)/2/t)
     & +d5*(fst+flt)*(2.*(s+u)/t+(t+u**2/t)/2./s)  +2.*d1
     & +fst/2.*(2./t+3./s+(u-2.*s2)/s/t)
     & +flt/2.*(-2./t+1./s-(u+2.*s2)/s/t)
     & +1./2./t-1./s+(2.*u-4.*s2*qq/t)/s/t  + 0.0
      xcfin=xcfin*(cf-ca/2.)
      RETURN
 
8     CONTINUE
C     8. qq1256
      xcfin =
     &  (fst+flt)*(4.*d1*d5*s**2/t-2.*d1*(1.-u/t)-4./t)
     &  +4.*qq*fst/t**2   + 0.0
      xcfin=xcfin*(cf-ca/2.)
      RETURN
 
9     CONTINUE
C     9. qq1278
      xcfin =
     &  (2.*fs2-fst-fsu+2.*fstu)*(-2.-u/t-t/u+2.*s2*(1./u+1./t
     &      -s2/u/t))/s
     & -4.*(1.+(u/t+t/u)/2.-s2*qq*(1./t**2+1./u**2))/s
      xcfin=xcfin*(cf-ca/2.)
      RETURN
 
10    CONTINUE
C     10. qqb5678LL
      xcfin =
     &  d9*fla*(-d10*s/t*(t-u)**2+2.*d5*(s-t)+3.*s/t+2.*s2/t-1.)
     & +d10*(u/t-1.)*(u-2.*s2)
     & +d3*d5*(fsu-flt-2.*ftu)*(2.*s*(s+u)+u**2)/t
     & +d5*((fst-ftu)*(1.+u/t+2.*t/u+4.*s*(1./t+1./u+s/t/u))
     &     -(flt+ftu)*(1+u/t))
     & +d3*(2.*ftu-fsu+flt)*(1.+2.*s/t+u/t) + 2.*d1*s/u
     & +ftu*(s+s2)/t/u-(fst+flt+1.-s/u)/t
      xcfin = xcfin/2.
      RETURN
 
11    CONTINUE
C     11. qqb5678LR
      xcfin =
     &  2.*d5*(d3*(fsu-flt-2.*ftu)+2.*(fst-ftu)/u)*s**2/t
     & +d3*(1.+2.*s/t-u/t)*(2.*ftu-fsu+flt)
     & +2.*(ftu-fsu)*(s+s2-u)/t/u - 2.*(ftu+flt)/t
      xcfin = xcfin/2.
      RETURN
 
12    CONTINUE
C     12. qqb1234
      xcfin =
     &  d9*fla*(d10*(d10*3.*s*qq*(2.*s+t+u)*(t-u)**2/t
     &       +qq*(2.*s+t+u)*(u-t-s)/t) - qq/t)
     & +d10*(-d10*3./2./t*(s+qq-s2)*(2.*s+t+u)*(t-u)**2
     &       +(2.*s+t+u)*(3./2.+(s-s2)/t-u/2./t))
      xcfin = xcfin/2.
      RETURN
 
      END
C-----------------------------------------------------------------------
C     Dilogarithm Li2(x) for -1 <= x <= 1
C
      DOUBLE PRECISION FUNCTION dilog (x)
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      EXTERNAL dilitg
      pi=3.141592654
      IF (DABS(x).GT.1.d0) PAUSE
     &  'Dilog: Argument outside range [-1,1]'
      IF (x.EQ.1.d0) THEN
          dilog = pi**2/6.
          RETURN
      ENDIF
      IF (x.EQ.0.d0) THEN
          dilog=0.
          RETURN
      ENDIF
      xx=x
      IF (x.GT..5d0) xx=1.d0-x
      IF (x.LT.0.d0) xx=-x/(1.-x)
      CALL dilgau (dilitg,0.d0,xx,dilog)
      IF (x.GT..5d0) THEN
          dilog=-dilog+pi**2/6.-dlog(x)*dlog(xx)
      ELSEIF (x.LT.0.d0) THEN
          dilog=-dilog-.5d0*(dlog(1.-x))**2
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION dilitg (x)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      dilitg = -dlog(1.d0-x)/x
      RETURN
      END
C-----------------------------------------------------------------------
C      Gaussian Integrator - from Numerical Recipes - page 122
C
      SUBROUTINE dilgau (func,a,b,ss)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION x(5), w(5)
      DATA x /.1488743389d0,.4333953941d0,.6794095682d0,
     &        .8650633666d0,.9739065285d0/
      DATA w /.2955242247d0,.2692667193d0,.2190863625d0,
     &        .1494513491d0,.0666713443d0/
      xm=0.5d0*(b+a)
      xr=0.5d0*(b-a)
      ss=0.d0
      DO 11 j=1,5
          dx=xr*x(j)
          ss=ss+w(j)*(func(xm+dx)+func(xm-dx))
11    CONTINUE
      ss=xr*ss
      RETURN
      END
