c Written by Daniel de Florian. Last update 18/05/2011 

C COMPUTES THE QT RESUMMATION INCLUDING RAPIDITY
C  USES DOUBLE MELLIN INVERSION and Bessel Quadratures
C Gives dsigma/dqt/dy and the integral over y
C Attention, where it says eta it actually means RAPIDITY=y

!     PROGRAM MAIN
      double precision function resumm(costh,mm,qtt,yy,mode)
      implicit none
c     IMPLICIT DOUBLE PRECISION (A-I,L-Z)
      double precision costh,mm,qtt,yy
      integer mode
      integer mod              !mode 0: differential mode 1: integrated in costh mode 2: integrated in costh and y
      common /mod/mod
      double precision resu
      double precision cthmom0,cthmom1,cthmom2
      INTEGER N, NVP, ISET, NAORD, NNF, ih,flag1,IER,flag
      double precision mv2,g_param !  real *16 mm,qtt,yy,theta
      integer icoll,ih1,ih2,inorm
      integer cc,mord,kk,sii,sjj
C     Parameters for Bessel quadratures
      integer lenaw
      parameter (lenaw = 8000)
      double precision tiny
      parameter (tiny = 1.0d-307)
      double precision aw(0 : lenaw - 1)
      COMMON/QUADRATURES/aw
c     alphas(Mz) from lhapdf
      double precision amz
      common/couple/amz
      COMMON/collider/ih1,ih2
      double precision Q2i, Q2MUR, Q2MUF
      COMMON / SCALESbis / Q2i, Q2MUR, Q2MUF
      double precision alpqf,alpqr
      COMMON / COUPL  / ALPQF, ALPQR
      double precision c,co,si,ax
      COMMON / CONT   / C, CO, SI, AX
      double precision eta
      COMMON / ETA / ETA
      double precision etamax,etamin
      COMMON/ETAlim/ETAmax,etamin
      double precision etalim
      COMMON/etalimite/etalim
      double precision etam
      COMMON / NAORD / NAORD
c      double precision v
c      COMMON/v/v
      COMMON/NFLAVORS/nnF
      COMMON/morder/mord
      COMMON/flag1/flag1
      COMMON/flag/flag
      COMMON/g_param/g_param
c     COMMON/binteg/phi2,min,max
      double precision Vud,Vus,Vub,Vcd,Vcs,Vcb
      common/cabib/Vud,Vus,Vub,Vcd,Vcs,Vcb
      include 'sigmaij_inc.f'
      double precision g
      common/NP/g
      common/dyphoton/phot
      double complex loga,logmuf2q2,logq2muf2,logq2mur2
      common/clogs/loga,logmuf2q2,logq2muf2,logq2mur2
      
c     cached for invres and cachecoeff
      DOUBLE PRECISION ax1,ax2,xx1,xx2
      integer nmax1,nmax2
      common /cxx/ax1,ax2,xx1,xx2,nmax1,nmax2

      COMPLEX*16 cCEX1(136),cCEX2p(136),cCEX2m(136)
      COMMON/ccex/cCEX1,cCEX2p,cCEX2m

      INTEGER I
      COMPLEX*16 CCp,CCm, Np(136),Nm(136)
      COMMON / MOMS2    / Np,Nm,CCP,CCm

      logical changedmass
      double precision facZ,facW,chi1,chi2
      common/vcoup/facZ,facW,chi1,chi2

      double precision q2s,x
      common/q2s/q2s,x

      double precision xborn,xborn2,xnorm,xnormal

      double precision dyalphas_mcfm,dyalphas_lhapdf
      external dyalphas_mcfm,dyalphas_lhapdf

      double precision xsection
      external xsection
      double precision xsection2
      external xsection2

      integer flag5,brflag,fnwa
      common/flags2/flag5,brflag,fnwa

      double precision cqt,cq2,cy
      common/kinematic/cqt,cq2,cy
C     
c      include 'ewinput.f' 
      include 'masses.f' 
      include 'const.h' 
c     include 'constants.f' 
      include 'sudakov_inc.f'
      include 'scales_inc.f' 
      double precision gevpb
      data gevpb/3.8937966d8/
      double precision gevfb
      data gevfb/3.8937966d11/
      include 'zerowidth.f'
      integer order,nproc,iih1,iih2,phot
      common/nnlo/order
      common/nproc/nproc
      double precision sqs
      common/energy/SQS 
      common/density/iih1,iih2
      include 'scale.f'
      include 'facscale.f'
      include 'quadrules.f'

      cqt=qtt
      cq2=mm*mm
      cy=yy
      
      mod = mode
      if(flag.eq.0)  then       !! ONE TIME INITIALIZATION
         !print *, 'Warming up: first call to resum ... '
         write(*,'(A)',advance='no') 'First call to resum ... '
         !write(*,'()',advance='no') 'First call to resum ... '
         call resummconst
         call reno2const !this seems to be unused
     
C     Initialization for Gauss inversion
         if (approxpdf.eq.0) then
            call initmellingauss ! higher order quadrature rule
         else
            call INITO           ! dyres original quadrature rule
         endif
         print *,'Done'
c     Initialization of redundant variables       
         flag1=order            ! set flag1 to order of calculation (carbon copy of order, 1=NLO+NLL, 2=NNLO+NNLL)
c     Choose pp or ppbar collider
         ih1=iih1               !1
         ih2=iih2               !1

c     Set factorization and renormalization scales (to work with dynamic scale need to move this outside init stage)
         mur=scale
         muf=facscale
C     Scales           
         mur2=mur**2
         muf2=muf**2

         q2mur=mur2
         q2muf=muf2

C non-perturbative parameter
         g=g_param

c     precompute log of scales
         loga=log(a_param)
         rloga=log(a_param)

c     narrow width, no branching ratio, always 0
         fnwa=0
         brflag=0

c     choose real axis (complex plane) integration of bstar (b) (always 0)
c         flagrealcomplex=0

c     set flag5 = 3 for Z, 21 for W+, 22 for W- (if phot is on (-1) flag5 is 5)
         if(nproc.eq.3) then
            flag5=3-phot*2
         elseif(nproc.eq.1.or.nproc.eq.2) then
            flag5=nproc+20
         endif
         
c     narrow width approximation (never used)      
         if (fnwa.eq.1) then
            if (flag5.eq.3)  q=Mz
            if (flag5.eq.2.or.flag5.eq.21.or.flag5.eq.22) q=Mw
         endif

c     normal (imod=0) or modified (imod=1) sudakov      
c         imod=1

         nnf=nf                 !number of flavours (why duplicated?)
         b0p=b0*a_param
       
c     check on remove branching ratio, not used
         if (brflag.eq.1.and.fnwa.eq.0) then
            write(*,*)"ERROR: Remove BR flag can be true, only if Narrow Width
     /Approximation flag is true!"
            stop
         endif

c     flag1 is the order of calculation (carbon copy of order, 1=NLO+NLL, 2=NNLO+NNLL)
         if(flag1.eq.1) then
            mord=1              ! Matching at LO
            iord=0              ! LO evolution
            naord=1
         elseif(flag1.eq.2) then
            mord=2              ! Matching at NLO
            iord=1              ! NLO evolution
            naord=2
         endif

      endif                     ! end initialization


C   ALPQR = ALPHA AT RENORMALIZATION SCALE
      if (approxpdf.eq.1) then
         ALPQR=dyalphas_mcfm(dsqrt(q2mur),amz,3)/4d0/pi
      else
         ALPQR=dyalphas_lhapdf(dsqrt(q2mur))/4d0/pi
      endif
C as = ALPHAS/PI
      aass = ALPQR*4d0

      g=g_param
      loga=log(a_param)
      rloga=log(a_param)
      
c******************************************
c     mass dependent part (cache mass value)
      if (q.eq.mm) then
         changedmass=.false.
      else
         changedmass=.true.
         q=mm
         q2=q**2
      endif
c******************************************


c******************************************
c     Initialise the born-level matrix elements sigmaij
c     ---> sigmaij depends on mass and costh
c     when costh is integrated and costh moments are provided, sigmaij
c     are convoluted with the rapidity depent expression
c      and they depend also on rapidity
c.....Born cross sections
      if (mod.eq.0) then
         call initsigma(mm,costh)
      elseif (mod.eq.1) then
         call cthmoments(cthmom0,cthmom1,cthmom2)
         call initsigmacth(mm,cthmom0,cthmom1,cthmom2)
      endif
c******************************************

c*****************************************
c     mass dependent part
      if(changedmass) then      !.or.changedeta) then

C     squared of resummation scale (used for the evolution)       
         q2s=q2/a_param**2
C     
c     careful, changed name of common block energy/sroot/ to energy/SQS/
         X=Q2/SQS**2
         sroot=SQS 
         shad=sroot**2
C...  PARAMETERS OF THE INTEGRATION CONTOUR IN THE COMPLEX N PLANE :
         AX=dlog(X)

C     Kinematical limit
         etalim=-0.5*ax
C     Limit eta_max to avoid reaching the end of phase space       
         etam=IDINT(etalim*10)/10d0

C...  COUPLING CONSTANTS AT INPUT SCALE = ALPHAS/4/PI
C     ALPQF = ALPHA AT RESUMMATION SCALE (used to start evolution,

      if (approxpdf.eq.1) then
         ALPQF=dyalphas_mcfm(dsqrt(q2s),amz,3)/4d0/pi
      else
         ALPQF=dyalphas_lhapdf(dsqrt(q2s))/4d0/pi
      endif
         
c     precompute scales
         logmuf2q2=log(muf2/q2)
         logq2muf2=log(q2/muf2)
         logq2mur2=log(q2/mur2)
         rlogq2mur2=log(q2/mur2)
         rblim=b0p*(1/q)*exp(1/(2*aass*beta0)) ! avoid Landau pole     
         if(q2.gt.2.2d0*mur2) then ! AVOID LARGE VALE OF b IN f2(y) WHEN q2>mur2 
            rblim=b0p*(1/q)*exp(1/(2*aass*beta0))*dsqrt(mur2/q2)
         endif
         cblim=b0p*(1/q)*exp(1/(2*aass*beta0))   
         if(q2.gt.2.2d0*mur2) then
            cblim=b0p*(1/q)*exp(1/(2*aass*beta0))*dsqrt(mur2/q2)
         endif
      endif
c************************

c*****************************************
c     initialization part
       if(flag.eq.0) then
             
          if (mod.eq.0.and.approxpdf.eq.1) then
C     Call fitting routine
             CALL fiteador(x,muf2,etam)
C     CALL initialization subroutine for PDFS and Anomalous dimensions
             CALL INITOFIT
          else
C     initialization of exact PDF moments and anomalous dimensions
             call initmoments
          endif

C     define points for quadratures integration
          call intdeoini(lenaw, tiny, 1.0d-2, aw)

          flag=1
      endif

c*****************************************
!      write(*,*) flag
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC       
C
C introduce qt and eta from here   
C START Phase Space generation from here!!!!
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Notice that "Higgs Mass" will also change when producing the Phase Space
C
 
       qt=qtt
C Generate rapidity of gauge-boson
c*****************************************
c     qt and rapidity dependent part, qt and eta are stored in a common block
       eta=yy !etta       
C Evaluate kinematics to choose best fit and logs       
       if (flag5.eq.3)  mv2=Mz**2 !missing flag5=5, but this has not impact since mv2 in evaluatekin is not used
         if (flag5.eq.2.or.flag5.eq.21.or.flag5.eq.22) mv2=Mw**2
c     evaluatekin select IFIT, which rapidity dependent,
c     be careful, evaluatekin also set qt in the common block
c     mv2 is not used         
       call evaluatekin(eta,qt,mv2)
c*****************************************
c       print *,IFIT
c*****************************************
c AX=Log[z1*z2]
C eta=1/2 (log[z1]-Log[z2])
c cache ax1, ax2, xx1, xx2 for invres(b)
        AX1 = (AX+2*ETA)/2d0 
        AX2 = (AX-2*ETA)/2d0 
        xx1=dExp(ax1)
        xx2=dExp(ax2)

c to avoid issues at x>0.87d0 (very large rapidities) suppress the result       
        if (xx1.gt.0.87d0) xnormal=(1-xx1)**3/(1-0.87d0)**3
        if (xx2.gt.0.87d0) xnormal=(1-xx2)**3/(1-0.87d0)**3

c in rapidity integrated mode nmax1 and nmax2 are evaluate for the maximum and minimum rapidity        
        if (mod.eq.0.and.approxpdf.eq.1) then
        if ( xX1 .le. 0.001 ) then                          
           NMAX1 = 40           ! zmax = 2 
        else if ( xX1 .le. 0.05 ) then                   
           NMAX1 = 40           ! zmax = 4  
        else if ( Xx1 .le. 0.2 ) then                   
           NMAX1 = 56           ! zmax = 8  
        else if ( xX1 .le. 0.4 ) then                   
           NMAX1 = 72           ! zmax = 12 
        else if ( xX1 .le. 0.7 ) then                   
           NMAX1 = 88           ! zmax = 18 
        else                                                          
c     NMAX1 = 136                ! zmax = 36
           NMAX1 = 88           ! zmax = 18 
        end if                  
        
        if ( xX2 .le. 0.001 ) then                          
           NMAX2 = 40           ! zmax = 2 
        else if ( xX2 .le. 0.05 ) then                   
           NMAX2 = 40           ! zmax = 4  
        else if ( xX2 .le. 0.2 ) then                   
           NMAX2 = 56           ! zmax = 8  
        else if ( xX2 .le. 0.4 ) then                   
           NMAX2 = 72           ! zmax = 12 
        else if ( xX2 .le. 0.7 ) then                   
           NMAX2 = 88           ! zmax = 18 
        else                                                          
c     NMAX2 = 136                ! zmax = 36
           NMAX2 = 88           ! zmax = 18 
        end if                  
      else
c     In costh and rapidity integrated modes remove rapidity dependences
c     Also the integration path in the complex plane is in this case along the imaginary axis,
c     need to integrate up to higher z max 18 (22)
         NMAX1 = mdim
         NMAX2 = mdim
      end if
c***************************************
c rapidity (and mass) dependence in ax1 ax2        
      do I = 1, max(NMAX1,NMAX2) !136
      cCEX1(I) = EXP (-Np(I) * AX1) / pi * CCp
      cCEX2p(I) = EXP (-Np(I) * AX2) / pi * CCp
      cCEX2m(I) = EXP (-Nm(I) * AX2) / pi * CCm
c     here can cache the products cCEX1(I)*cCEX2p(I) and cCEX1(I)*cCEX2p(I)
c and perform I1*I2 integrations in rapidity      
      enddo

c***************************************

c***************************************
c     mass dependence in logq2muf2 logq2mur2, and (avoidable) rapidity dependence in nnmax
c      if(changedmass) then !.or.changedeta) then
         call cachecoeff
c      endif
c***************************************

c*****************************************

c*****************************************
c     mass, qt dependent part
C Normalization
      xnorm=x*qt/2/Q2
c*****************************************

c     Apply rescaling from approximate to original PDF only in differential mode mod=0
c     to reproduce the original dyres result
      if (mod.eq.0.and.approxpdf.eq.1) then
c*****************************************
c     mass, rapidity costh dependent part (costh is in sigmaij for the the xsections)
C     xsection provides the born xsection computed using exact pdfs
         xborn=xsection(ax,eta,muf,ih1,ih2)
C     xsection2 provides the born xsection computed using approx. pdfs
         xborn2=xsection2(ax,eta,muf,ih1,ih2)
C     Normalize PDF's to born
         xnormal=xborn/xborn2
!     write(*,*) "y,qT,Q,xnormal",yy,qT,dsqrt(Q2),xnormal
C     Avoid problems at extreme kinematics (and eventual divisions by 0!)     
         if ((xnormal.gt.1.25).or.(xnormal.lt.0.75)) then
            xnormal=1d0
         endif
c     xnormal=1d0
c     print *,yy,xnormal
c*****************************************
      else
         xnormal = 1d0
      endif

      if (mod.eq.2.or.mod.eq.1) then
         xnormal = 1d0
      endif
      if (abs(eta) .gt. etalim) then
         resumm=0
         return
      endif
C     COMPUTE RESUMMED CONTRIBUTION (b integral 0->infinity)
c     print*,'phase space point in resumm', qtt, yy, mm, costh

cc*****************************************
cc     C++ rewritten check
c      call setmesq_expy(mod, mm, costh, yy)
c      call hcoeff_calc(aass,logmuf2q2,logq2muf2,logq2mur2,loga)
cc*****************************************
c     dependence on qt, m, also y, and costh unless integrated
      call resummation(resu)
c*****************************************
      resumm=resu*xnorm * xnormal
      END
C      
      Subroutine resummation(resu)
      IMPLICIT DOUBLE PRECISION (A - Z)
      integer lenaw
      parameter (lenaw = 8000)
      dimension aw(0 : lenaw - 1)
      external INVRES
      COMMON/transverse/qt
      COMMON/QUADRATURES/aw
C Compute integral using quadratures for Oscillating functions
C as described in intde2.f
      call intdeo(invres, 0.0d0, qt, aw, resu, errt)
c     print *,'dequad result of inverse bessel transform',resu, errt
      return
      end
 
 
      FUNCTION INVres (b)    
      IMPLICIT NONE
      DOUBLE PRECISION INVRES,b,dbesj0,xj0,ax,eta,
     .     FZ,C,CO,SI,qt,FUN
      double precision funp1m1,funp2m2,funm1p1,funm2p2
c     cached from resumm
      DOUBLE PRECISION ax1,ax2,xx1,xx2
      integer nmax1,nmax2
      common /cxx/ax1,ax2,xx1,xx2,nmax1,nmax2

      INTEGER NMAX, I1, I2, I
      DOUBLE PRECISION  WN(136),WZ(8),Zs(8),DOWN(17),UP(17)
      COMPLEX*16 CCp,CCm, Np(136),Nm(136),scale2,sudak,S,factorfin
      COMPLEX*16 XNM1,XNM2,XN1,XN2,CEX1,CEX2,HCRN,INT1,INT2,bb
      COMPLEX*16 int1p1m1,int1p2m2,int1m1p1,int1m2p2
      COMPLEX*16 int2p1m1,int2p2m2,int2m1p1,int2m2p2
      COMPLEX*16 cCEX1(136),cCEX2p(136),cCEX2m(136)
      COMMON/ccex/cCEX1,cCEX2p,cCEX2m
      COMMON / CONT   / C, CO, SI, AX
      COMMON / ETA / ETA
      COMMON / WEIGHTS2 / WN
      COMMON / MOMS2    / Np,Nm,CCp,CCm
      COMMON/transverse/qt
      external S       

      DOUBLE COMPLEX alphasl
      DOUBLE PRECISION XL, XL1, SALP
      DOUBLE COMPLEX alpq,ALPr
      COMMON/alphasldata/XL,XL1,SALP,alpq,ALPr
      DOUBLE PRECISION  ALPQF, ALPQR
      COMMON / COUPL  / ALPQF, ALPQR
      INTEGER NAORD
      COMMON / NAORD / NAORD

      integer iq
      DOUBLE COMPLEX fn1(9),fn2(9),fn3(9)
      COMPLEX *16 fx1(-5:5), fx2(-5:5), fx3(-5:5)
      include 'pdfn_inc.f'

      integer ih1,ih2
      COMMON/collider/ih1,ih2

      DOUBLE COMPLEX utemp,dtemp
      COMPLEX *16 sHCRN(-2:2,-2:2)
      common/chcrn/sHCRN
      integer nproc
      common/nproc/nproc

      include 'sigmaij_inc.f'

      complex *16 sigmaintijp(-5:5,-5:5,136,136)
      complex *16 sigmaintijm(-5:5,-5:5,136,136)
      common/sigmaint/sigmaintijp,sigmaintijm

      integer mod              !mode 0: differential mode 1: integrated in costh mode 2: integrated in costh and y
      common /mod/mod

      include 'sudakov_inc.f'

      COMPLEX *16 mellinint_integrand

      double precision cqt,cq2,cy
      common/kinematic/cqt,cq2,cy

      double precision besselint
      external besselint
      
      include 'constants.f' 
      scale2=complex(b0p**2/b**2,0d0)
      bb=complex(b,0d0)
C     USES BESSEL FUNCTION SINCE INTEGRATION IS DONE ALONG THE REAL AXIS
c      print *, b
c     ********************
c     qt and b dependence (bessel function)
      xj0= 2*dbesj0((qt*b))
c     ********************
c     ********************
c     Sudakov is only mass and b dependent
      sudak=S(bb)
      if (sudak.eq.complex(0, 0)) then
         invres = 0
         return
      endif     
c     ********************
c     ********************
c     qt and mass dependence
      factorfin=bb*xj0*sudak
c     ********************

c     ********************
c     ax1, ax2, xx1, xx2 are cached in resumm, they are eta and mass dependence (AX is mass dependent)
c     ********************

c**************************************
c     b-dependence
*...  alphasl gives the LL/NLL evolution of alpha from Qres=Q/a_param to
*     q2=bo^2/b^2
      alpq= alpqf * alphasl(scale2)
      XL = ALPQF / ALPQ 
      SALP   = LOG (XL)
      XL1 = 1.- XL
C     SELECT ORDER FOR EVOLUTION LO/NLO
      ALPr= ALPQ * complex(1d0,0d0)*(NAORD-1)
c**************************************

c     ax1, ax2, xx1, xx2, ccex1, ccex2p, ccex2m are cached in resumm, they are eta and mass dependence (AX is mass dependent)


c     ***************************
      FUN = 0.D0
      funp1m1 = 0d0
      funp2m2 = 0d0
      funm1p1 = 0d0
      funm2p2 = 0d0
c       print *,'Start gauss loop',NMAX1*NMAX2
c This is a n=8 gaussian quadrature in the complex plane, at nodes I1, I2 of INTERESnew (I1, I2 are points in z space)
c I1 and I2 are 8*17=136 (17 intervals) for z = [0,36]
c dependence on b in SCALE2
c dependence on mass and eta in I1, I2 n points (can be avoided by setting NMAX1 and NMAX2, probably not useful)
c dependence on mass and eta in ax1, ax2
c dependence on mass and costh in sigmaij inside INTERESnew
c apparently no dependence on qt -> cache it and perform qt integration first? yes, but you have a b dependence
c The set of b points are determined by qt! can fully cache b points, but qt integration must performed at last
c linearise this loop to what is only I dependent, that is do a preliminary I loop to cache I-dependent values, and than the I1 I2 loop

c ******************************
c cache I dependence of fp1, fp2+, fp2-
      do I = 1, max(nmax1,nmax2) ! 136
c     output of reno2: gammaiab_j, depends on Ij z point, rapidity (in IFIT), depend on b (in scale2)
c     output of reno2: FP1 (FP2), depends on I1, (I2), IFIT, b
*  U   UB   D   DB   S   C    B    G 
*  1    2   3   4    5   6    7    8 
c     the subroutine reno2 performs the evolution of the mellin moments of the PDF, fx,
c     from muf to the scale corresponding to the impact parameter b
c     reno2 (fx, I, ISIGN, IBEAM, ALPS, SCALE2)
c
         call reno2 (cfx1(:,I), I, 1, 1, ALPQF, SCALE2)
         call reno2 (cfx2p(:,I), I, 1, 2, ALPQF, SCALE2)
         call reno2 (cfx2m(:,I), I, -1, 2, ALPQF, SCALE2)
c     convert proton into antiproton (need to switch only u and d, since the sea is assumed symmetric)
         if (ih1.eq.-1) then
            UTEMP=cfx1(1,I)
            cfx1(1,I)=cfx1(-1,I)
            cfx1(-1,I)=UTEMP
            DTEMP=cfx1(2,I)
            cfx1(2,I)=cfx1(-2,I)
            cfx1(-2,I)=DTEMP
         endif
         if (ih2.eq.-1) then
            UTEMP=cfx2p(1,I)
            cfx2p(1,I)=cfx2p(-1,I)
            cfx2p(-1,I)=UTEMP
            DTEMP=cfx2p(2,I)
            cfx2p(2,I)=cfx2p(-2,I)
            cfx2p(-2,I)=DTEMP

            UTEMP=cfx2m(1,I)
            cfx2m(1,I)=cfx2m(-1,I)
            cfx2m(-1,I)=UTEMP
            DTEMP=cfx2m(2,I)
            cfx2m(2,I)=cfx2m(-2,I)
            cfx2m(-2,I)=DTEMP
         endif                
c     print *,'fortran',b,I,cfx2m(0,I)
c     Cache the positive and negative branch of coefficients which depend only on one I index
         call cachehcoeff(I,1)
         call cachehcoeff(I,-1)
      enddo
cc     ******************** rapidity integrated cross section
c      if (mod.eq.0.or.mod.eq.1) then
c         do I = 1, NMAX1
c            CALL INTERESnew (HCRN, I, I, 1)
c            INT1= HCRN*exp(-Np(I)*AX)/pi*CCp
c
c            INT2 = 0
c            FZ=-DBLE(1./2*(INT1-INT2)*WN(I))
c            FUN= FUN + FZ
c         enddo
cc     ********************************************************      
      
c ******************************
c     C++ rewritten check
c      call hcoeff_calcb(aass,logmuf2q2,loga,alpq,aexp,aexpb)
c ******************************
      DO 100 I1 = 1, NMAX1
         DO 10 I2 = 1, NMAX2
c here scale2 is fixed (b-dependent), and the function is called many times at I1 I2 points       
c part of the coefficients calculation is hoisted in the previous I loop

c     merge positive and negative branch
c *************************************
c orig fortran code            
             if (mod.eq.0.or.mod.eq.1) then
C     COMPUTE POSITIVE BRANCH
           CALL INTERESnew (HCRN, I1, I2, 1)
           INT1= (HCRN * cCEX2p(I2))
           
C     COMPUTE NEGATIVE BRANCH
           CALL INTERESnew (HCRN, I1, I2, -1)
           INT2= (HCRN * cCEX2m(I2))
           
           FZ=-DBLE( 1./2*(INT1-INT2)*cCEX1(I1)*WN(I1)*WN(I2))

           FUN= FUN + FZ
        elseif (mod.eq.2) then
c     Rapidity integrated mode:
c     sigmaij are fatorised from HCRN and numerical integration in y is performed in rapintegrals
c     the full expression is HCRN(I1,I2)_ij * ccex(I1,I2) * sigma_ij
c     HCRN(I1,I2)_ij is only b dependent
c     ccex(I1,I2) is rapidity and mass dependent
c     sigma_ij is costh and mass dependent, but becomes rapidity dependent after integration of the costh moments
c     The integrals are solved analitically when no cuts on the leptons are applied
            
            CALL INTERESnew (HCRN, I1, I2, 1)
            int1p1m1=shcrn(1,-1)*sigmaintijp(1,-1,I1,I2)
            int1p2m2=shcrn(2,-2)*sigmaintijp(2,-2,I1,I2)
            int1m1p1=shcrn(-1,1)*sigmaintijp(-1,1,I1,I2)
            int1m2p2=shcrn(-2,2)*sigmaintijp(-2,2,I1,I2)
            
            CALL INTERESnew (HCRN, I1, I2, -1)
            int2p1m1=shcrn(1,-1)*sigmaintijm(1,-1,I1,I2)
            int2p2m2=shcrn(2,-2)*sigmaintijm(2,-2,I1,I2)
            int2m1p1=shcrn(-1,1)*sigmaintijm(-1,1,I1,I2)
            int2m2p2=shcrn(-2,2)*sigmaintijm(-2,2,I1,I2)

            funp1m1=funp1m1-DBLE((int1p1m1-int2p1m1)) !*WN(I1)*WN(I2))
            funp2m2=funp2m2-DBLE((int1p2m2-int2p2m2)) !*WN(I1)*WN(I2))
            funm1p1=funm1p1-DBLE((int1m1p1-int2m1p1)) !*WN(I1)*WN(I2))
            funm2p2=funm2p2-DBLE((int1m2p2-int2m2p2)) !*WN(I1)*WN(I2))
        endif    
c end of orig fortran code
c *************************************
ccc     C++ rewritten check
c            call pdfevol(I1,I2,1)
c            call mellinint_pdf_mesq_expy(I1,I2,1)
c            INT1=mellinint_integrand(I1,I2,1)
c            call pdfevol(I1,I2,2)
c            call mellinint_pdf_mesq_expy(I1,I2,2)
c            INT2=mellinint_integrand(I1,I2,2)
c            FZ=-DBLE(0.5d0*(INT1-INT2))
cc     FZ=-0.5d0*(DBLE(INT1)-DBLE(INT2))
c            FUN= FUN + FZ
cc            print *,'fortran',I1,I2,INT1,INT2
ccc *************************************


 10   CONTINUE
 100  CONTINUE
c     ***************************
      if (mod.eq.0.or.mod.eq.1) then
         INVRES = fun*factorfin
      elseif (mod.eq.2) then
c     ***************************
c orig fortran code            
         INVRES= 0.5d0*(funp1m1+funp2m2
     1        +funm1p1+funm2p2)*factorfin
c     ***************************
c     C++ rewritten check
c         INVRES = fun*factorfin
c     ***************************
      endif
      if (invres.ne.invres) then
         print *, 'Warning, invres =', invres, 'b =', b, 'qt=', qt
         invres = 0
      endif
c      print *,'fortran',b, INVRES, fun, factorfin
c      print *,b, INVRES, besselint(b,cqt,cq2)
      RETURN                           
      END


      subroutine cachehcoeff(I,ISIGN)
      implicit none
      integer I, SIG, ISIGN
      double complex H1q      
      
      DOUBLE PRECISION aassh,aasshsq
      DOUBLE PRECISION XL, XL1, SALP
      DOUBLE COMPLEX alpq,ALPr
      COMMON/alphasldata/XL,XL1,SALP,alpq,ALPr
      double complex loga,logmuf2q2,logq2muf2,logq2mur2
      common/clogs/loga,logmuf2q2,logq2muf2,logq2mur2
      
      include 'sudakov_inc.f'
      
      integer flag1
      COMMON/flag1/flag1

c     Cached values
      DOUBLE COMPLEX cgamma1qq(136,2),cgamma1qg(136,2),cgamma1gq(136,2),
     1     cgamma1gg(136,2),cgamma2qq(136,2),cgamma2qqV(136,2),
     2     cgamma2qqbV(136,2),cgamma2qqS(136,2),cgamma2qqbS(136,2),
     3     cgamma2qg(136,2),cgamma2gq(136,2),cgamma2gg(136,2)
      DOUBLE COMPLEX cC1QQ(136,2),cC1QG(136,2),cC1GQ(136,2),cC1GG
      DOUBLE COMPLEX cans(136,2),cam(136,2), cap(136,2), cal(136,2), 
     1     cbe(136,2),cab(136,2)
      DOUBLE COMPLEX crmin(136,2),crplus(136,2),
     1     crqq(136,2), crqg(136,2), crgq(136,2), crgg(136,2)
      DOUBLE COMPLEX cAC(136,2),cNMP(136,2),cNPM(136,2),cDMQQ(136,2),
     1     cDMQG(136,2),cDMGQ(136,2),cDMGG(136,2),cDPQQ(136,2),
     2     cDPQG(136,2),cDPGQ(136,2),cDPGG(136,2),cRMMQQ(136,2),
     3     cRMMQG(136,2),cRMMGQ(136,2),
     3     cRMMGG(136,2),cRMPQQ(136,2),cRMPQG(136,2),cRMPGQ(136,2),
     4     cRMPGG(136,2),cRPMQQ(136,2),cRPMQG(136,2),cRPMGQ(136,2),
     5     cRPMGG(136,2),cRPPQQ(136,2),cRPPQG(136,2),cRPPGQ(136,2),
     6     cRPPGG(136,2)
      COMPLEX*16 cC2qgM(136,2),cC2NSqqM(136,2),cC2SqqbM(136,2),
     1     cC2NSqqbM(136,2)
      COMMON /ANOMCACHE/cans,cam,cap,cal,cbe,cab,crmin,crplus,crqq,
     1     crqg,crgq,crgg,
     2     cAC,cNMP,cNPM,cDMQQ,
     3     cDMQG,cDMGQ,cDMGG,cDPQQ,cDPQG,
     4     cDPGQ,cDPGG,cRMMQQ,cRMMQG,cRMMGQ,
     5     cRMMGG,cRMPQQ,cRMPQG,cRMPGQ,
     6     cRMPGG,cRPMQQ,cRPMQG,cRPMGQ,
     7     cRPMGG,cRPPQQ,cRPPQG,cRPPGQ,
     8     cRPPGG,cC1QQ,cC1QG,cC1GQ,cC1GG,
     9     cgamma1qq,cgamma1qg,cgamma1gq,
     1     cgamma1gg,cgamma2qq,cgamma2qqV,
     2     cgamma2qqbV,cgamma2qqS,cgamma2qqbS,
     3     cgamma2qg,cgamma2gq,cgamma2gg,
     4     cC2qgM,cC2NSqqM,cC2SqqbM,cC2NSqqbM

      DOUBLE COMPLEX cHqqb(136,136,2)
      DOUBLE COMPLEX cH1stqqb(136,136,2),cH1stqg(136,2),cH1stgg,
     1     cH2stqqp(136,2),cH2stqq(136,2),cH2stqqb(136,136,2),
     2     cH2stqg_1(136,136,2),cH2stqg_2(136,136,2),cH2stgg(136,136,2)
      COMMON /ccoefficients/cHqqb,
     1     cH1stqqb,cH1stqg,cH1stgg,
     2     cH2stqqp,cH2stqq,cH2stqqb,cH2stqg_1,cH2stqg_2,cH2stgg

c     output
      double complex cHqgnll(136,2),cHqqnnll(136,2),cHqqpnnll(136,2),
     1     caexpqq(136,2),caexpqg(136,2)
      COMMON /chcoeff/cHqgnll,cHqqnnll,cHqqpnnll,caexpqq,caexpqg
      include 'const.h'


c     **************************************
c     retrieve values cached in cacheanom for
C     NORMALIZED ANOMALOUS DIMENSIONS AND COEFFICIENTS
      SIG = (-ISIGN+1)/2+1
      H1q=dcmplx(0d0,0d0)
      aassh=aass/2
      aasshsq=(aass/2)**2

c     begin NLL
      if(flag1.eq.1) then 
         
         cHqgnll(I,SIG) = alpq*2d0*cC1qg(I,SIG) 
     .        +(aass/2)*(-cgamma1qg(I,SIG))*(logmuf2q2+2*loga)
      endif                     ! end NLL
c     begin NNLL
c     these NNLL coefficients depend on I, SIG
c     and also on b through aexpb, aexp
      if(flag1.eq.2) then

         caexpqq(I,SIG)=aexpB**(aassh*(cC1qq(I,SIG)-pisq329))
         caexpqg(I,SIG)=aexpB**(aassh*(cC2qgM(I,SIG)/cC1qg(I,SIG)
     .        -pisq329))
         cHqqnnll(I,SIG)=aasshsq*cH2stqq(I,SIG)
     .        *caexpqq(I,SIG)*(aexp)**2
         cHqqpnnll(I,SIG)=aasshsq*cH2stqqp(I,SIG)
     .        *caexpqq(I,SIG)*(aexp)**2
         
      endif                     ! end NNLL
c     *********************************************************
      return
      end
        
 
c**********************************
c input I1 I2, scale2(b)
c output HCRN (function of what? b, rapidity, I1, I2)
c do positive and negative in the same loop
c break coefficient calculation, which is not dependent on sigmaij
       SUBROUTINE INTERESnew (HCRN, I1,I2, isign)
c       IMPLICIT DOUBLE COMPLEX (A - Z)
       implicit none
       double complex HCRN
       DOUBLE PRECISION aassh,aasshsq
       DOUBLE COMPLEX aexpqq1,aexpqq2,aexpqg1,aexpqg2

       complex *16 FX1(-5:5), FX2(-5:5)
       complex *16 GGN,QGN_1,QGN_2,QQBN_1,QQBN_2,QQBN_3,QQBN_4
       complex *16 sGGN(-2:2,-2:2),sQGN_1(-2:2,-2:2),sQGN_2(-2:2,-2:2)
       complex *16 sQQBN_1(-2:2,-2:2),sQQBN_2(-2:2,-2:2)
       complex *16 sQQBN_3(-2:2,-2:2),sQQBN_4(-2:2,-2:2)
       complex *16 Hqqb,Hqg_1,Hqg_2,Hgg,Hqq_1,Hqq_2,Hqq,Hqqp_1,Hqqp_2

       complex *16 sHCRN(-2:2,-2:2)
       common/chcrn/sHCRN

       INTEGER I,J,flag1,I1,I2,isign,ih1,ih2,sig
       INTEGER SI,SJ,SK
       include 'const.h'
       include 'scales_inc.f'

       COMMON/flag1/flag1

       COMMON/collider/ih1,ih2

       include 'sudakov_inc.f'
       include 'sigmaij_inc.f'

      integer nproc
      common/nproc/nproc

c     Cached values
      DOUBLE COMPLEX cgamma1qq(136,2),cgamma1qg(136,2),cgamma1gq(136,2),
     1     cgamma1gg(136,2),cgamma2qq(136,2),cgamma2qqV(136,2),
     2     cgamma2qqbV(136,2),cgamma2qqS(136,2),cgamma2qqbS(136,2),
     3     cgamma2qg(136,2),cgamma2gq(136,2),cgamma2gg(136,2)
      DOUBLE COMPLEX cC1QQ(136,2),cC1QG(136,2),cC1GQ(136,2),cC1GG
      DOUBLE COMPLEX cans(136,2),cam(136,2), cap(136,2), cal(136,2), 
     1     cbe(136,2),cab(136,2)
      DOUBLE COMPLEX crmin(136,2),crplus(136,2),
     1     crqq(136,2), crqg(136,2), crgq(136,2), crgg(136,2)
      DOUBLE COMPLEX cAC(136,2),cNMP(136,2),cNPM(136,2),cDMQQ(136,2),
     1     cDMQG(136,2),cDMGQ(136,2),cDMGG(136,2),cDPQQ(136,2),
     2     cDPQG(136,2),cDPGQ(136,2),cDPGG(136,2),cRMMQQ(136,2),
     3     cRMMQG(136,2),cRMMGQ(136,2),
     3     cRMMGG(136,2),cRMPQQ(136,2),cRMPQG(136,2),cRMPGQ(136,2),
     4     cRMPGG(136,2),cRPMQQ(136,2),cRPMQG(136,2),cRPMGQ(136,2),
     5     cRPMGG(136,2),cRPPQQ(136,2),cRPPQG(136,2),cRPPGQ(136,2),
     6     cRPPGG(136,2)
      COMPLEX*16 cC2qgM(136,2),cC2NSqqM(136,2),cC2SqqbM(136,2),
     1     cC2NSqqbM(136,2)
      COMMON /ANOMCACHE/cans,cam,cap,cal,cbe,cab,crmin,crplus,crqq,
     1     crqg,crgq,crgg,
     2     cAC,cNMP,cNPM,cDMQQ,
     3     cDMQG,cDMGQ,cDMGG,cDPQQ,cDPQG,
     4     cDPGQ,cDPGG,cRMMQQ,cRMMQG,cRMMGQ,
     5     cRMMGG,cRMPQQ,cRMPQG,cRMPGQ,
     6     cRMPGG,cRPMQQ,cRPMQG,cRPMGQ,
     7     cRPMGG,cRPPQQ,cRPPQG,cRPPGQ,
     8     cRPPGG,cC1QQ,cC1QG,cC1GQ,cC1GG,
     9     cgamma1qq,cgamma1qg,cgamma1gq,
     1     cgamma1gg,cgamma2qq,cgamma2qqV,
     2     cgamma2qqbV,cgamma2qqS,cgamma2qqbS,
     3     cgamma2qg,cgamma2gq,cgamma2gg,
     4     cC2qgM,cC2NSqqM,cC2SqqbM,cC2NSqqbM

      DOUBLE COMPLEX cHqqb(136,136,2)
      DOUBLE COMPLEX cH1stqqb(136,136,2),cH1stqg(136,2),cH1stgg,
     1     cH2stqqp(136,2),cH2stqq(136,2),cH2stqqb(136,136,2),
     2     cH2stqg_1(136,136,2),cH2stqg_2(136,136,2),cH2stgg(136,136,2)
      COMMON /ccoefficients/cHqqb,
     1     cH1stqqb,cH1stqg,cH1stgg,
     2     cH2stqqp,cH2stqq,cH2stqqb,cH2stqg_1,cH2stqg_2,cH2stgg

c      DIMENSION FP1(9), fp2(9), FN(15)

c     cached fps
      integer iq
      include 'pdfn_inc.f'

c      cached hcoefficients
      double complex cHqgnll(136,2),cHqqnnll(136,2),cHqqpnnll(136,2),
     1     caexpqq(136,2),caexpqg(136,2)
      COMMON /chcoeff/cHqgnll,cHqqnnll,cHqqpnnll,caexpqq,caexpqg

c     integration mode
      integer mod               !mode 0: differential mode 1: integrated in costh mode 2: integrated in costh and y
      common /mod/mod
      
cc**************************************
c retrieve values cached in cacheanom for
C NORMALIZED ANOMALOUS DIMENSIONS AND COEFFICIENTS
      SIG = (-ISIGN+1)/2+1

c     begin LL
      if(flag1.eq.0) then
         Hqqb=0d0 !cHqqb(I1,I2,SIG)
C     all other channels are =0 at LL
         Hqg_1= 0d0
         Hqg_2= 0d0
         Hgg = 0d0                !  Q Q -> Qb Q  = Qb Qb -> Q Qb  
         Hqq_1 = 0d0            !  Qb Qb -> Qb Q =  Q Q -> Q Qb      
         Hqq_2 = 0d0            ! Average QQ->QQb  and QbQb->QQb
         Hqq = 0d0               !  qqp_1  means   Q' Q -> Qb Q  flavor in sigmaQQb determined by "second parton"
         Hqqp_1 = 0d0           !  qqp_2  means   Q Q' -> Q Qb  flavor in sigmaQQb determined by "first parton"
         Hqqp_2 = 0d0
      endif
      
c begin NLL
c the Hqqb NLL coefficients depend only on I1 and I2 (and isign)
c Hqg_1(2) also on b through alpq in Cqg_1(2)
      if(flag1.eq.1) then 
         
         Hqqb=cHqqb(I1,I2,SIG)
c     *******************************
c     b dependence
         Hqg_1=cHqgnll(I1,1)
         Hqg_2=cHqgnll(I2,SIG)
c     *******************************
C     all other channels are =0 at NLL
C     GG
         Hgg = 0d0                !  Q Q -> Qb Q  = Qb Qb -> Q Qb  
         Hqq_1 = 0d0            !  Qb Qb -> Qb Q =  Q Q -> Q Qb      
         Hqq_2 = 0d0            ! Average QQ->QQb  and QbQb->QQb
         Hqq = 0d0               !  qqp_1  means   Q' Q -> Qb Q  flavor in sigmaQQb determined by "second parton"
         Hqqp_1 = 0d0           !  qqp_2  means   Q Q' -> Q Qb  flavor in sigmaQQb determined by "first parton"
         Hqqp_2 = 0d0
      endif
c     END NLL

C
C            qg_1 means GQ initial state
C            qg_2 means QG initial state
C
c begin NNLL
c **************************************************
c the NNLL coefficients depend on I1 I2
c and some of them also on b through aexpb
       if(flag1.eq.2) then
c      H1q=dcmplx(0d0,0d0)
          aassh=aass/2
          aasshsq=(aass/2)**2

c     cached I1 I2 and ISIGN dependence
c     also b dependence

C  QQb
       Hqqb = (1+aassh*cH1stqqb(I1,I2,SIG)
     .     +aasshsq*cH2stqqb(I1,I2,SIG))
     .  *caexpqq(I1,1)              !leg q
     .  *caexpqq(I2,SIG)             !leg qb

C  qg_1 means GQ initial state
       Hqg_1 = (aassh*cH1stqg(I1,1)+aasshsq*cH2stqg_1(I1,I2,SIG))
     .   *aexp
     .   *caexpqg(I1,1)  
     .   *caexpqq(I1,1)

C  qg_2 means QG initial state          
       Hqg_2 = (aassh*cH1stqg(I2,SIG)+aasshsq*cH2stqg_2(I1,I2,SIG))
     .   *aexp
     .   *caexpqg(I2,SIG)  
     .   *caexpqq(I2,SIG) 

C  GG
       Hgg = aasshsq*cH2stgg(I1,I2,SIG)
     .   *aexp*caexpqg(I1,1)  
     .   *aexp*caexpqg(I2,SIG)              

       Hqq_1 = cHqqnnll(I1,1)   !  Q Q -> Qb Q  = Qb Qb -> Q Qb
       Hqq_2 = cHqqnnll(I2,SIG)   !  Qb Qb -> Qb Q =  Q Q -> Q Qb
       Hqq= (Hqq_1+Hqq_2)/2d0     ! Average QQ->QQb  and QbQb->QQb
       Hqqp_1 = cHqqpnnll(I1,1)   ! qqp_1  means   Q' Q -> Qb Q  flavor in sigmaQQb determined by "second parton"
       Hqqp_2 = cHqqpnnll(I2,SIG) ! qqp_2  means   Q Q' -> Q Qb  flavor in sigmaQQb determined by "first parton"
      endif ! end NNLL


      FX1( 1) = cfx1( 1,I1)
      FX1( 2) = cfx1( 2,I1)
      FX1( 3) = cfx1( 3,I1)
      FX1( 4) = cfx1( 4,I1)
      FX1( 5) = cfx1( 5,I1)
      FX1( 0) = cfx1( 0,I1)
      FX1(-1) = cfx1(-1,I1)
      FX1(-2) = cfx1(-2,I1)
      FX1(-3) = cfx1(-3,I1)
      FX1(-4) = cfx1(-4,I1)
      FX1(-5) = cfx1(-5,I1)
      if (ISIGN.eq.1) then
         FX2( 1) = cfx2p( 1,I2)
         FX2( 2) = cfx2p( 2,I2)
         FX2( 3) = cfx2p( 3,I2)
         FX2( 4) = cfx2p( 4,I2)
         FX2( 5) = cfx2p( 5,I2)
         FX2( 0) = cfx2p( 0,I2)
         FX2(-1) = cfx2p(-1,I2)
         FX2(-2) = cfx2p(-2,I2)
         FX2(-3) = cfx2p(-3,I2)
         FX2(-4) = cfx2p(-4,I2)
         FX2(-5) = cfx2p(-5,I2)
      else
         FX2( 1) = cfx2m( 1,I2)
         FX2( 2) = cfx2m( 2,I2)
         FX2( 3) = cfx2m( 3,I2)
         FX2( 4) = cfx2m( 4,I2)
         FX2( 5) = cfx2m( 5,I2)
         FX2( 0) = cfx2m( 0,I2)
         FX2(-1) = cfx2m(-1,I2)
         FX2(-2) = cfx2m(-2,I2)
         FX2(-3) = cfx2m(-3,I2)
         FX2(-4) = cfx2m(-4,I2)
         FX2(-5) = cfx2m(-5,I2)
      endif

c The following loops are the main core of the resummed calculation
c FX1 (FX2) depend on I1 (I2), b, mass, rapidity
c sigmaij depends on mass, costh


       GGN=0
c       QGN_1=0
c       QGN_2=0
c       QQBN_1=0
       QQBN_2=0
       QQBN_3=0
       QQBN_4=0
       if (mod.eq.2) then
          sGGN   (1,-1)=0
          sQQBN_2(1,-1)=0
          sQQBN_3(1,-1)=0
          sQQBN_4(1,-1)=0
          sGGN   (2,-2)=0
          sQQBN_2(2,-2)=0
          sQQBN_3(2,-2)=0
          sQQBN_4(2,-2)=0
          sGGN   (-1,1)=0
          sQQBN_2(-1,1)=0
          sQQBN_3(-1,1)=0
          sQQBN_4(-1,1)=0
          sGGN   (-2,2)=0
          sQQBN_2(-2,2)=0
          sQQBN_3(-2,2)=0
          sQQBN_4(-2,2)=0
       endif
c Speed up integration for Z/gamma, calculate only antidiagonal sj=-si terms, unrolled loops, factorise u-type and d-type quarks
       if (nproc.eq.3) then     !
          if (mod.eq.0.or.mod.eq.1) then
cc**************************************
cc NLL part
      QGN_1 = (FX2(-1)+FX2(-4))*sigmaij(1,-1)
     2      + (FX2(-2)+FX2(-3)+FX2(-5))*sigmaij(2,-2)
     1      + (FX2(1)+FX2(4))*sigmaij(-1,1)
     2      + (FX2(2)+FX2(3)+FX2(5))*sigmaij(-2,2)

      QGN_2 = (FX1(1)+FX1(4))*sigmaij(1,-1)
     2      + (FX1(2)+FX1(3)+FX1(5))*sigmaij(2,-2)
     1      + (FX1(-1)+FX1(-4))*sigmaij(-1,1)
     2      + (FX1(-2)+FX1(-3)+FX1(-5))*sigmaij(-2,2)

      QQBN_1 = (FX1(1)*FX2(-1)+FX1(4)*FX2(-4))*sigmaij(1,-1)
     2  + (FX1(2)*FX2(-2)+FX1(3)*FX2(-3)+FX1(5)*FX2(-5))*sigmaij(2,-2)
     1  + (FX1(-1)*FX2(1)+FX1(-4)*FX2(4))*sigmaij(-1,1)
     2  + (FX1(-2)*FX2(2)+FX1(-3)*FX2(3)+FX1(-5)*FX2(5))*sigmaij(-2,2)

c **************************************

c **************************************
c     NNLL
          if(flag1.eq.2) then 
             GGN = FX1(0)*FX2(0)*
     1            (2*sigmaij(1,-1)
     2           + 3*sigmaij(2,-2)
     1           + 2*sigmaij(-1,1)
     2            + 3*sigmaij(-2,2))
             
             QQBN_2 = (FX1(1)*FX2(1)+FX1(4)*FX2(4))*sigmaij(1,-1)
     2     + (FX1(2)*FX2(2)+FX1(3)*FX2(3)+FX1(5)*FX2(5))*sigmaij(2,-2)
     1     + (FX1(-1)*FX2(-1)+FX1(-4)*FX2(-4))*sigmaij(-1,1)
     2     + (FX1(-2)*FX2(-2)+FX1(-3)*FX2(-3)
     3       +FX1(-5)*FX2(-5))*sigmaij(-2,2)

         QQBN_3 = FX1(1)*( FX2(-2)*sigmaij(2,-2)
     3                   + FX2(-3)*sigmaij(2,-2)
     4                   + FX2(-4)*sigmaij(1,-1)
     5                   + FX2(-5)*sigmaij(2,-2)
     2                   + FX2(2)*sigmaij(-2,2)
     3                   + FX2(3)*sigmaij(-2,2)
     4                   + FX2(4)*sigmaij(-1,1)
     5                   + FX2(5)*sigmaij(-2,2))
     1          + FX1(2)*(FX2(-1)*sigmaij(1,-1)
     3                   + FX2(-3)*sigmaij(2,-2)
     4                   + FX2(-4)*sigmaij(1,-1)
     5                   + FX2(-5)*sigmaij(2,-2)
     1                   + FX2(1)*sigmaij(-1,1)
     3                   + FX2(3)*sigmaij(-2,2)
     4                   + FX2(4)*sigmaij(-1,1)
     5                   + FX2(5)*sigmaij(-2,2))
     1          + FX1(3)*(FX2(-1)*sigmaij(1,-1)
     2                   + FX2(-2)*sigmaij(2,-2)
     4                   + FX2(-4)*sigmaij(1,-1)
     5                   + FX2(-5)*sigmaij(2,-2)
     1                   + FX2(1)*sigmaij(-1,1)
     2                   + FX2(2)*sigmaij(-2,2)
     4                   + FX2(4)*sigmaij(-1,1)
     5                   + FX2(5)*sigmaij(-2,2))
     1          + FX1(4)*(FX2(-1)*sigmaij(1,-1)
     2                   + FX2(-2)*sigmaij(2,-2)
     3                   + FX2(-3)*sigmaij(2,-2)
     5                   + FX2(-5)*sigmaij(2,-2)
     1                   + FX2(1)*sigmaij(-1,1)
     2                   + FX2(2)*sigmaij(-2,2)
     3                   + FX2(3)*sigmaij(-2,2)
     5                   + FX2(5)*sigmaij(-2,2))
     1          + FX1(5)*(FX2(-1)*sigmaij(1,-1)
     2                   + FX2(-2)*sigmaij(2,-2)
     3                   + FX2(-3)*sigmaij(2,-2)
     4                   + FX2(-4)*sigmaij(1,-1)
     1                   + FX2(1)*sigmaij(-1,1)
     2                   + FX2(2)*sigmaij(-2,2)
     3                   + FX2(3)*sigmaij(-2,2)
     4                   + FX2(4)*sigmaij(-1,1))
     1          + FX1(-1)*(FX2(-2)*sigmaij(2,-2)
     3                   + FX2(-3)*sigmaij(2,-2)
     4                   + FX2(-4)*sigmaij(1,-1)
     5                   + FX2(-5)*sigmaij(2,-2)
     2                   + FX2(2)*sigmaij(-2,2)
     3                   + FX2(3)*sigmaij(-2,2)
     4                   + FX2(4)*sigmaij(-1,1)
     5                   + FX2(5)*sigmaij(-2,2))
     1          + FX1(-2)*(FX2(-1)*sigmaij(1,-1)
     3                   + FX2(-3)*sigmaij(2,-2)
     4                   + FX2(-4)*sigmaij(1,-1)
     5                   + FX2(-5)*sigmaij(2,-2)
     1                   + FX2(1)*sigmaij(-1,1)
     3                   + FX2(3)*sigmaij(-2,2)
     4                   + FX2(4)*sigmaij(-1,1)
     5                   + FX2(5)*sigmaij(-2,2))
     1          + FX1(-3)*(FX2(-1)*sigmaij(1,-1)
     2                   + FX2(-2)*sigmaij(2,-2)
     4                   + FX2(-4)*sigmaij(1,-1)
     5                   + FX2(-5)*sigmaij(2,-2)
     1                   + FX2(1)*sigmaij(-1,1)
     2                   + FX2(2)*sigmaij(-2,2)
     4                   + FX2(4)*sigmaij(-1,1)
     5                   + FX2(5)*sigmaij(-2,2))
     1          + FX1(-4)*(FX2(-1)*sigmaij(1,-1)
     2                   + FX2(-2)*sigmaij(2,-2)
     3                   + FX2(-3)*sigmaij(2,-2)
     5                   + FX2(-5)*sigmaij(2,-2)
     1                   + FX2(1)*sigmaij(-1,1)
     2                   + FX2(2)*sigmaij(-2,2)
     3                   + FX2(3)*sigmaij(-2,2)
     5                   + FX2(5)*sigmaij(-2,2))
     1          + FX1(-5)*(FX2(-1)*sigmaij(1,-1)
     2                   + FX2(-2)*sigmaij(2,-2)
     3                   + FX2(-3)*sigmaij(2,-2)
     4                   + FX2(-4)*sigmaij(1,-1)
     1                   + FX2(1)*sigmaij(-1,1)
     2                   + FX2(2)*sigmaij(-2,2)
     3                   + FX2(3)*sigmaij(-2,2)
     4                   + FX2(4)*sigmaij(-1,1))

         QQBN_4 = FX2(1)*( FX1(2)*sigmaij(2,-2)
     3                   + FX1(3)*sigmaij(2,-2)
     4                   + FX1(4)*sigmaij(1,-1)
     5                   + FX1(5)*sigmaij(2,-2)
     2                   + FX1(-2)*sigmaij(-2,2)
     3                   + FX1(-3)*sigmaij(-2,2)
     4                   + FX1(-4)*sigmaij(-1,1)
     5                   + FX1(-5)*sigmaij(-2,2))
     1          + FX2(2)* (FX1(1)*sigmaij(1,-1)
     3                   + FX1(3)*sigmaij(2,-2)
     4                   + FX1(4)*sigmaij(1,-1)
     5                   + FX1(5)*sigmaij(2,-2)
     1                   + FX1(-1)*sigmaij(-1,1)
     3                   + FX1(-3)*sigmaij(-2,2)
     4                   + FX1(-4)*sigmaij(-1,1)
     5                   + FX1(-5)*sigmaij(-2,2))
     1          + FX2(3)* (FX1(1)*sigmaij(1,-1)
     2                   + FX1(2)*sigmaij(2,-2)
     4                   + FX1(4)*sigmaij(1,-1)
     5                   + FX1(5)*sigmaij(2,-2)
     1                   + FX1(-1)*sigmaij(-1,1)
     2                   + FX1(-2)*sigmaij(-2,2)
     4                   + FX1(-4)*sigmaij(-1,1)
     5                   + FX1(-5)*sigmaij(-2,2))
     1          + FX2(4)* (FX1(1)*sigmaij(1,-1)
     2                   + FX1(2)*sigmaij(2,-2)
     3                   + FX1(3)*sigmaij(2,-2)
     5                   + FX1(5)*sigmaij(2,-2)
     1                   + FX1(-1)*sigmaij(-1,1)
     2                   + FX1(-2)*sigmaij(-2,2)
     3                   + FX1(-3)*sigmaij(-2,2)
     5                   + FX1(-5)*sigmaij(-2,2))
     1          + FX2(5)* (FX1(1)*sigmaij(1,-1)
     2                   + FX1(2)*sigmaij(2,-2)
     3                   + FX1(3)*sigmaij(2,-2)
     4                   + FX1(4)*sigmaij(1,-1)
     1                   + FX1(-1)*sigmaij(-1,1)
     2                   + FX1(-2)*sigmaij(-2,2)
     3                   + FX1(-3)*sigmaij(-2,2)
     4                   + FX1(-4)*sigmaij(-1,1))
     1          + FX2(-1)*(FX1(2)*sigmaij(2,-2)
     3                   + FX1(3)*sigmaij(2,-2)
     4                   + FX1(4)*sigmaij(1,-1)
     5                   + FX1(5)*sigmaij(2,-2)
     2                   + FX1(-2)*sigmaij(-2,2)
     3                   + FX1(-3)*sigmaij(-2,2)
     4                   + FX1(-4)*sigmaij(-1,1)
     5                   + FX1(-5)*sigmaij(-2,2))
     1          + FX2(-2)*(FX1(1)*sigmaij(1,-1)
     3                   + FX1(3)*sigmaij(2,-2)
     4                   + FX1(4)*sigmaij(1,-1)
     5                   + FX1(5)*sigmaij(2,-2)
     1                   + FX1(-1)*sigmaij(-1,1)
     3                   + FX1(-3)*sigmaij(-2,2)
     4                   + FX1(-4)*sigmaij(-1,1)
     5                   + FX1(-5)*sigmaij(-2,2))
     1          + FX2(-3)*(FX1(1)*sigmaij(1,-1)
     2                   + FX1(2)*sigmaij(2,-2)
     4                   + FX1(4)*sigmaij(1,-1)
     5                   + FX1(5)*sigmaij(2,-2)
     1                   + FX1(-1)*sigmaij(-1,1)
     2                   + FX1(-2)*sigmaij(-2,2)
     4                   + FX1(-4)*sigmaij(-1,1)
     5                   + FX1(-5)*sigmaij(-2,2))
     1          + FX2(-4)*(FX1(1)*sigmaij(1,-1)
     2                   + FX1(2)*sigmaij(2,-2)
     3                   + FX1(3)*sigmaij(2,-2)
     5                   + FX1(5)*sigmaij(2,-2)
     1                   + FX1(-1)*sigmaij(-1,1)
     2                   + FX1(-2)*sigmaij(-2,2)
     3                   + FX1(-3)*sigmaij(-2,2)
     5                   + FX1(-5)*sigmaij(-2,2))
     1          + FX2(-5)*(FX1(1)*sigmaij(1,-1)
     2                   + FX1(2)*sigmaij(2,-2)
     3                   + FX1(3)*sigmaij(2,-2)
     4                   + FX1(4)*sigmaij(1,-1)
     1                   + FX1(-1)*sigmaij(-1,1)
     2                   + FX1(-2)*sigmaij(-2,2)
     3                   + FX1(-3)*sigmaij(-2,2)
     4        + FX1(-4)*sigmaij(-1,1))
      endif
      elseif(mod.eq.2) then
      sQGN_1(1,-1)=FX1(0)*(FX2(-1)+FX2(-4))        
      sQGN_1(2,-2)=FX1(0)*(FX2(-2)+FX2(-3)+FX2(-5))
      sQGN_1(-1,1)=FX1(0)*(FX2(1)+FX2(4))          
      sQGN_1(-2,2)=FX1(0)*(FX2(2)+FX2(3)+FX2(5))   
      
      sQGN_2(1,-1)=FX2(0)*(FX1(1)+FX1(4))          
      sQGN_2(2,-2)=FX2(0)*(FX1(2)+FX1(3)+FX1(5))   
      sQGN_2(-1,1)=FX2(0)*(FX1(-1)+FX1(-4))        
      sQGN_2(-2,2)=FX2(0)*(FX1(-2)+FX1(-3)+FX1(-5))
      
      sQQBN_1(1,-1)=FX1(1)*FX2(-1)+FX1(4)*FX2(-4)
      sQQBN_1(2,-2)=FX1(2)*FX2(-2)+FX1(3)*FX2(-3)+FX1(5)*FX2(-5)
      sQQBN_1(-1,1)=FX1(-1)*FX2(1)+FX1(-4)*FX2(4)
      sQQBN_1(-2,2)=FX1(-2)*FX2(2)+FX1(-3)*FX2(3)+FX1(-5)*FX2(5)

c     NNLL
      if(flag1.eq.2) then
      sGGN(1,-1)=FX1(0)*FX2(0)*2			
      sGGN(2,-2)=FX1(0)*FX2(0)*3			
      sGGN(-1,1)=FX1(0)*FX2(0)*2			
      sGGN(-2,2)=FX1(0)*FX2(0)*3			

      sQQBN_2(1,-1)=FX1(1)*FX2(1)+FX1(4)*FX2(4)
      sQQBN_2(2,-2)=FX1(2)*FX2(2)+FX1(3)*FX2(3)+FX1(5)*FX2(5)
      sQQBN_2(-1,1)=FX1(-1)*FX2(-1)+FX1(-4)*FX2(-4)
      sQQBN_2(-2,2)=FX1(-2)*FX2(-2)+FX1(-3)*FX2(-3)+FX1(-5)*FX2(-5)

      sQQBN_3(1,-1)=FX1(1)*FX2(-4)			               
     .             +(FX1(2)*FX2(-1)+FX1(2)*FX2(-4))
     .             +(FX1(3)*FX2(-1)+FX1(3)*FX2(-4))
     .             +FX1(4)*FX2(-1)                                
     .             +(FX1(5)*FX2(-1)+FX1(5)*FX2(-4))
     .             +FX1(-1)*FX2(-4)                                  
     .             +(FX1(-2)*FX2(-1)+FX1(-2)*FX2(-4))
     .             +(FX1(-3)*FX2(-1)+FX1(-3)*FX2(-4))
     .             +FX1(-4)*FX2(-1)                                  
     .             +(FX1(-5)*FX2(-1)+FX1(-5)*FX2(-4))

      sQQBN_3(2,-2)=FX1(1)*(FX2(-2)+FX2(-3)+FX2(-5))
     .		   +(FX1(2)*FX2(-3)+FX1(2)*FX2(-5))
     .		   +(FX1(3)*FX2(-2)+FX1(3)*FX2(-5))
     .		   +FX1(4)*(FX2(-2)+FX2(-3)+FX2(-5))
     .		   +(FX1(5)*FX2(-2)+FX1(5)*FX2(-3))
     .		   +FX1(-1)*(FX2(-2)+FX2(-3)+FX2(-5))
     .		   +(FX1(-2)*FX2(-3)+FX1(-2)*FX2(-5))
     .		   +(FX1(-3)*FX2(-2)+FX1(-3)*FX2(-5))
     .		   +FX1(-4)*(FX2(-2)+FX2(-3)+FX2(-5))
     .		   +(FX1(-5)*FX2(-2)+FX1(-5)*FX2(-3))

      sQQBN_3(-1,1)=FX1(1)*FX2(4)                                    
     .		   +(FX1(2)*FX2(1)+FX1(2)*FX2(4))
     .		   +(FX1(3)*FX2(1)+FX1(3)*FX2(4))
     .		   +FX1(4)*FX2(1)                                 
     .		   +(FX1(5)*FX2(1)+FX1(5)*FX2(4))
     .		   +FX1(-1)*FX2(4)                                   
     .		   +(FX1(-2)*FX2(1)+FX1(-2)*FX2(4))
     .		   +(FX1(-3)*FX2(1)+FX1(-3)*FX2(4))
     .		   +FX1(-4)*FX2(1)                                   
     .		   +(FX1(-5)*FX2(1)+FX1(-5)*FX2(4))

      sQQBN_3(-2,2)=FX1(1)*(FX2(2)+FX2(3)+FX2(5))
     .		   +(FX1(2)*FX2(3)+FX1(2)*FX2(5))
     .		   +(FX1(3)*FX2(2)+FX1(3)*FX2(5))
     .		   +FX1(4)*(FX2(2)+FX2(3)+FX2(5))
     .		   +(FX1(5)*FX2(2)+FX1(5)*FX2(3))
     .		   +FX1(-1)*(FX2(2)+FX2(3)+FX2(5))
     .		   +(FX1(-2)*FX2(3)+FX1(-2)*FX2(5))
     .		   +(FX1(-3)*FX2(2)+FX1(-3)*FX2(5))
     .		   +FX1(-4)*(FX2(2)+FX2(3)+FX2(5))
     .		   +(FX1(-5)*FX2(2)+FX1(-5)*FX2(3))

      sQQBN_4(1,-1)=FX2(1)*FX1(4)
     .            +(FX2(2)*FX1(1)+FX2(2)*FX1(4))
     .            +(FX2(3)*FX1(1)+FX2(3)*FX1(4))
     .            +FX2(4)*FX1(1)
     .            +(FX2(5)*FX1(1)+FX2(5)*FX1(4))
     .            +FX2(-1)*FX1(4)
     .            +(FX2(-2)*FX1(1)+FX2(-2)*FX1(4))
     .            +(FX2(-3)*FX1(1)+FX2(-3)*FX1(4))
     .            +FX2(-4)*FX1(1)
     .            +(FX2(-5)*FX1(1)+FX2(-5)*FX1(4))

      sQQBN_4(2,-2)=FX2(1)*(FX1(2)+FX1(3)+FX1(5))
     .            +(FX2(2)*FX1(3)+FX2(2)*FX1(5))
     .            +(FX2(3)*FX1(2)+FX2(3)*FX1(5))
     .            +FX2(4)*(FX1(2)+FX1(3)+FX1(5))
     .            +(FX2(5)*FX1(2)+FX2(5)*FX1(3))
     .            +FX2(-1)*(FX1(2)+FX1(3)+FX1(5))
     .            +(FX2(-2)*FX1(3)+FX2(-2)*FX1(5))
     .            +(FX2(-3)*FX1(2)+FX2(-3)*FX1(5))
     .            +FX2(-4)*(FX1(2)+FX1(3)+FX1(5))
     .            +(FX2(-5)*FX1(2)+FX2(-5)*FX1(3))

      sQQBN_4(-1,1)=FX2(1)*FX1(-4)
     .            +(FX2(2)*FX1(-1)+FX2(2)*FX1(-4))
     .            +(FX2(3)*FX1(-1)+FX2(3)*FX1(-4))
     .            +FX2(4)*FX1(-1)
     .            +(FX2(5)*FX1(-1)+FX2(5)*FX1(-4))
     .            +FX2(-1)*FX1(-4)
     .            +(FX2(-2)*FX1(-1)+FX2(-2)*FX1(-4))
     .            +(FX2(-3)*FX1(-1)+FX2(-3)*FX1(-4))
     .            +FX2(-4)*FX1(-1)
     .            +(FX2(-5)*FX1(-1)+FX2(-5)*FX1(-4))

      sQQBN_4(-2,2)=FX2(1)*(FX1(-2)+FX1(-3)+FX1(-5))
     .            +(FX2(2)*FX1(-3)+FX2(2)*FX1(-5))
     .            +(FX2(3)*FX1(-2)+FX2(3)*FX1(-5))
     .            +FX2(4)*(FX1(-2)+FX1(-3)+FX1(-5))
     .            +(FX2(5)*FX1(-2)+FX2(5)*FX1(-3))
     .            +FX2(-1)*(FX1(-2)+FX1(-3)+FX1(-5))
     .            +(FX2(-2)*FX1(-3)+FX2(-2)*FX1(-5))
     .            +(FX2(-3)*FX1(-2)+FX2(-3)*FX1(-5))
     .            +FX2(-4)*(FX1(-2)+FX1(-3)+FX1(-5))
     .            +(FX2(-5)*FX1(-2)+FX2(-5)*FX1(-3))
      endif
      endif
c **************************************
c     Loop unrolling for W+ to speed up integration
      else if (nproc.eq.1) then     ! 
cc **************************************
cc NLL part
      QGN_1 = FX2(-2)*sigmaij(1,-2)
     2      + FX2(-3)*sigmaij(1,-3)
     3      + FX2(-5)*sigmaij(1,-5)
     1      + FX2(-2)*sigmaij(4,-2)
     2      + FX2(-3)*sigmaij(4,-3)
     3      + FX2(-5)*sigmaij(4,-5)
     2      + FX2(1)*sigmaij(-2,1)
     3      + FX2(4)*sigmaij(-2,4)
     2      + FX2(1)*sigmaij(-3,1)
     3      + FX2(4)*sigmaij(-3,4)
     2      + FX2(1)*sigmaij(-5,1)
     3      + FX2(4)*sigmaij(-5,4)
      
      QGN_2 = FX1(1)*sigmaij(1,-2)
     2      + FX1(1)*sigmaij(1,-3)
     3      + FX1(1)*sigmaij(1,-5)
     1      + FX1(4)*sigmaij(4,-2)
     2      + FX1(4)*sigmaij(4,-3)
     3      + FX1(4)*sigmaij(4,-5)
     2      + FX1(-2)*sigmaij(-2,1)
     3      + FX1(-2)*sigmaij(-2,4)
     2      + FX1(-3)*sigmaij(-3,1)
     3      + FX1(-3)*sigmaij(-3,4)
     2      + FX1(-5)*sigmaij(-5,1)
     3      + FX1(-5)*sigmaij(-5,4)

      QQBN_1 = FX1(1)*FX2(-2)*sigmaij(1,-2)
     2       + FX1(1)*FX2(-3)*sigmaij(1,-3)
     3       + FX1(1)*FX2(-5)*sigmaij(1,-5)
     1       + FX1(4)*FX2(-2)*sigmaij(4,-2)
     2       + FX1(4)*FX2(-3)*sigmaij(4,-3)
     3       + FX1(4)*FX2(-5)*sigmaij(4,-5)
     2       + FX1(-2)*FX2(1)*sigmaij(-2,1)
     3       + FX1(-2)*FX2(4)*sigmaij(-2,4)
     2       + FX1(-3)*FX2(1)*sigmaij(-3,1)
     3       + FX1(-3)*FX2(4)*sigmaij(-3,4)
     2       + FX1(-5)*FX2(1)*sigmaij(-5,1)
     3       + FX1(-5)*FX2(4)*sigmaij(-5,4)
cc **************************************
c **************************************
c     NNLL
      if(flag1.eq.2) then 
         GGN = FX1(0)*FX2(0)*
     1        (sigmaij(1,-2)
     2        +   sigmaij(1,-3)
     3        +   sigmaij(1,-5)
     4        +   sigmaij(4,-2)
     5        +   sigmaij(4,-3)
     6        +   sigmaij(4,-5)
     7        +   sigmaij(-2,1)
     8        +   sigmaij(-2,4)
     9        +   sigmaij(-3,1)
     1        +   sigmaij(-3,4)
     2        +   sigmaij(-5,1)
     3        +   sigmaij(-5,4))

         QQBN_2 = FX1(1)*FX2(2)*sigmaij(1,-2)
     2          + FX1(1)*FX2(3)*sigmaij(1,-3)
     3          + FX1(1)*FX2(5)*sigmaij(1,-5)
     1          + FX1(4)*FX2(2)*sigmaij(4,-2)
     2          + FX1(4)*FX2(3)*sigmaij(4,-3)
     3          + FX1(4)*FX2(5)*sigmaij(4,-5)
     2          + FX1(-2)*FX2(-1)*sigmaij(-2,1)
     3          + FX1(-2)*FX2(-4)*sigmaij(-2,4)
     2          + FX1(-3)*FX2(-1)*sigmaij(-3,1)
     3          + FX1(-3)*FX2(-4)*sigmaij(-3,4)
     2          + FX1(-5)*FX2(-1)*sigmaij(-5,1)
     3          + FX1(-5)*FX2(-4)*sigmaij(-5,4)

         QQBN_3 =
     1    (FX1(1)+FX1(-1))*
     1        ( FX2(-2)*sigmaij(4,-2)
     2        + FX2(-3)*sigmaij(4,-3)
     3        + FX2(-5)*sigmaij(4,-5)
     2        + FX2(1)*sigmaij(-2,1)
     3        + FX2(4)*sigmaij(-2,4)
     2        + FX2(1)*sigmaij(-3,1)
     3        + FX2(4)*sigmaij(-3,4)
     2        + FX2(1)*sigmaij(-5,1)
     3        + FX2(4)*sigmaij(-5,4))
     1  + (FX1(2)+FX1(-2))*
     1        ( FX2(-2)*sigmaij(1,-2)
     2        + FX2(-3)*sigmaij(1,-3)
     3        + FX2(-5)*sigmaij(1,-5)
     1        + FX2(-2)*sigmaij(4,-2)
     2        + FX2(-3)*sigmaij(4,-3)
     3        + FX2(-5)*sigmaij(4,-5)
     2        + FX2(1)*sigmaij(-3,1)
     3        + FX2(4)*sigmaij(-3,4)
     2        + FX2(1)*sigmaij(-5,1)
     3        + FX2(4)*sigmaij(-5,4))
     1        + (FX1(3)+FX1(-3))*
     1        ( FX2(-2)*sigmaij(1,-2)
     2        + FX2(-3)*sigmaij(1,-3)
     3        + FX2(-5)*sigmaij(1,-5)
     1        + FX2(-2)*sigmaij(4,-2)
     2        + FX2(-3)*sigmaij(4,-3)
     3        + FX2(-5)*sigmaij(4,-5)
     2        + FX2(1)*sigmaij(-2,1)
     3        + FX2(4)*sigmaij(-2,4)
     2        + FX2(1)*sigmaij(-5,1)
     3        + FX2(4)*sigmaij(-5,4))
     1  + (FX1(4)+FX1(-4))*
     1        ( FX2(-2)*sigmaij(1,-2)
     2        + FX2(-3)*sigmaij(1,-3)
     3        + FX2(-5)*sigmaij(1,-5)
     2        + FX2(1)*sigmaij(-2,1)
     3        + FX2(4)*sigmaij(-2,4)
     2        + FX2(1)*sigmaij(-3,1)
     3        + FX2(4)*sigmaij(-3,4)
     2        + FX2(1)*sigmaij(-5,1)
     3        + FX2(4)*sigmaij(-5,4))
     1  + (FX1(5)+FX1(-5))*
     1        ( FX2(-2)*sigmaij(1,-2)
     2        + FX2(-3)*sigmaij(1,-3)
     3        + FX2(-5)*sigmaij(1,-5)
     1        + FX2(-2)*sigmaij(4,-2)
     2        + FX2(-3)*sigmaij(4,-3)
     3        + FX2(-5)*sigmaij(4,-5)
     2        + FX2(1)*sigmaij(-2,1)
     3        + FX2(4)*sigmaij(-2,4)
     2        + FX2(1)*sigmaij(-3,1)
     3        + FX2(4)*sigmaij(-3,4))

         QQBN_4 = 
     1     (FX2(1)+FX2(-1))*
     1       (FX1(1)*sigmaij(1,-2)
     2      + FX1(1)*sigmaij(1,-3)
     3      + FX1(1)*sigmaij(1,-5)
     1      + FX1(4)*sigmaij(4,-2)
     2      + FX1(4)*sigmaij(4,-3)
     3      + FX1(4)*sigmaij(4,-5)
     3      + FX1(-2)*sigmaij(-2,4)
     3      + FX1(-3)*sigmaij(-3,4)
     3      + FX1(-5)*sigmaij(-5,4))
     1  + (FX2(2)+FX2(-2))*
     2       (FX1(1)*sigmaij(1,-3)
     3      + FX1(1)*sigmaij(1,-5)
     2      + FX1(4)*sigmaij(4,-3)
     3      + FX1(4)*sigmaij(4,-5)
     2      + FX1(-2)*sigmaij(-2,1)
     3      + FX1(-2)*sigmaij(-2,4)
     2      + FX1(-3)*sigmaij(-3,1)
     3      + FX1(-3)*sigmaij(-3,4)
     2      + FX1(-5)*sigmaij(-5,1)
     3      + FX1(-5)*sigmaij(-5,4))
     1  + (FX2(3)+FX2(-3))*
     1       (FX1(1)*sigmaij(1,-2)
     3      + FX1(1)*sigmaij(1,-5)
     1      + FX1(4)*sigmaij(4,-2)
     3      + FX1(4)*sigmaij(4,-5)
     2      + FX1(-2)*sigmaij(-2,1)
     3      + FX1(-2)*sigmaij(-2,4)
     2      + FX1(-3)*sigmaij(-3,1)
     3      + FX1(-3)*sigmaij(-3,4)
     2      + FX1(-5)*sigmaij(-5,1)
     3      + FX1(-5)*sigmaij(-5,4))
     1  + (FX2(4)+FX2(-4))*
     1       (FX1(1)*sigmaij(1,-2)
     2      + FX1(1)*sigmaij(1,-3)
     3      + FX1(1)*sigmaij(1,-5)
     1      + FX1(4)*sigmaij(4,-2)
     2      + FX1(4)*sigmaij(4,-3)
     3      + FX1(4)*sigmaij(4,-5)
     2      + FX1(-2)*sigmaij(-2,1)
     2      + FX1(-3)*sigmaij(-3,1)
     2      + FX1(-5)*sigmaij(-5,1))
     1  + (FX2(5)+FX2(-5))*
     1       (FX1(1)*sigmaij(1,-2)
     2      + FX1(1)*sigmaij(1,-3)
     1      + FX1(4)*sigmaij(4,-2)
     2      + FX1(4)*sigmaij(4,-3)
     2      + FX1(-2)*sigmaij(-2,1)
     3      + FX1(-2)*sigmaij(-2,4)
     2      + FX1(-3)*sigmaij(-3,1)
     3      + FX1(-3)*sigmaij(-3,4)
     2      + FX1(-5)*sigmaij(-5,1)
     3      + FX1(-5)*sigmaij(-5,4))
         
      endif
c     **************************************


c     Loop unrolling for W- to speed up integration
      else if (nproc.eq.2) then     ! 
cc **************************************
cc NLL part
      QGN_1 = FX2(2)*sigmaij(-1,2)
     2      + FX2(3)*sigmaij(-1,3)
     3      + FX2(5)*sigmaij(-1,5)
     1      + FX2(2)*sigmaij(-4,2)
     2      + FX2(3)*sigmaij(-4,3)
     3      + FX2(5)*sigmaij(-4,5)
     2      + FX2(-1)*sigmaij(2,-1)
     3      + FX2(-4)*sigmaij(2,-4)
     2      + FX2(-1)*sigmaij(3,-1)
     3      + FX2(-4)*sigmaij(3,-4)
     2      + FX2(-1)*sigmaij(5,-1)
     3      + FX2(-4)*sigmaij(5,-4)
      
      QGN_2 = FX1(-1)*sigmaij(-1,2)
     2      + FX1(-1)*sigmaij(-1,3)
     3      + FX1(-1)*sigmaij(-1,5)
     1      + FX1(-4)*sigmaij(-4,2)
     2      + FX1(-4)*sigmaij(-4,3)
     3      + FX1(-4)*sigmaij(-4,5)
     2      + FX1(2)*sigmaij(2,-1)
     3      + FX1(2)*sigmaij(2,-4)
     2      + FX1(3)*sigmaij(3,-1)
     3      + FX1(3)*sigmaij(3,-4)
     2      + FX1(5)*sigmaij(5,-1)
     3      + FX1(5)*sigmaij(5,-4)

      QQBN_1 = FX1(-1)*FX2(2)*sigmaij(-1,2)
     2       + FX1(-1)*FX2(3)*sigmaij(-1,3)
     3       + FX1(-1)*FX2(5)*sigmaij(-1,5)
     1       + FX1(-4)*FX2(2)*sigmaij(-4,2)
     2       + FX1(-4)*FX2(3)*sigmaij(-4,3)
     3       + FX1(-4)*FX2(5)*sigmaij(-4,5)
     2       + FX1(2)*FX2(-1)*sigmaij(2,-1)
     3       + FX1(2)*FX2(-4)*sigmaij(2,-4)
     2       + FX1(3)*FX2(-1)*sigmaij(3,-1)
     3       + FX1(3)*FX2(-4)*sigmaij(3,-4)
     2       + FX1(5)*FX2(-1)*sigmaij(5,-1)
     3       + FX1(5)*FX2(-4)*sigmaij(5,-4)
cc **************************************
c **************************************
c     NNLL
      if(flag1.eq.2) then 
         GGN = FX1(0)*FX2(0)*
     1           (sigmaij(-1,2)
     2        +   sigmaij(-1,3)
     3        +   sigmaij(-1,5)
     4        +   sigmaij(-4,2)
     5        +   sigmaij(-4,3)
     6        +   sigmaij(-4,5)
     7        +   sigmaij(2,-1)
     8        +   sigmaij(2,-4)
     9        +   sigmaij(3,-1)
     1        +   sigmaij(3,-4)
     2        +   sigmaij(5,-1)
     3        +   sigmaij(5,-4))

         QQBN_2 = FX1(-1)*FX2(-2)*sigmaij(-1,2)
     2          + FX1(-1)*FX2(-3)*sigmaij(-1,3)
     3          + FX1(-1)*FX2(-5)*sigmaij(-1,5)
     1          + FX1(-4)*FX2(-2)*sigmaij(-4,2)
     2          + FX1(-4)*FX2(-3)*sigmaij(-4,3)
     3          + FX1(-4)*FX2(-5)*sigmaij(-4,5)
     2          + FX1(2)*FX2(1)*sigmaij(2,-1)
     3          + FX1(2)*FX2(4)*sigmaij(2,-4)
     2          + FX1(3)*FX2(1)*sigmaij(3,-1)
     3          + FX1(3)*FX2(4)*sigmaij(3,-4)
     2          + FX1(5)*FX2(1)*sigmaij(5,-1)
     3          + FX1(5)*FX2(4)*sigmaij(5,-4)

         QQBN_3 =
     1    (FX1(-1)+FX1(1))*
     1        ( FX2(2)*sigmaij(-4,2)
     2        + FX2(3)*sigmaij(-4,3)
     3        + FX2(5)*sigmaij(-4,5)
     2        + FX2(-1)*sigmaij(2,-1)
     3        + FX2(-4)*sigmaij(2,-4)
     2        + FX2(-1)*sigmaij(3,-1)
     3        + FX2(-4)*sigmaij(3,-4)
     2        + FX2(-1)*sigmaij(5,-1)
     3        + FX2(-4)*sigmaij(5,-4))
     1  + (FX1(-2)+FX1(2))*
     1        ( FX2(2)*sigmaij(-1,2)
     2        + FX2(3)*sigmaij(-1,3)
     3        + FX2(5)*sigmaij(-1,5)
     1        + FX2(2)*sigmaij(-4,2)
     2        + FX2(3)*sigmaij(-4,3)
     3        + FX2(5)*sigmaij(-4,5)
     2        + FX2(-1)*sigmaij(3,-1)
     3        + FX2(-4)*sigmaij(3,-4)
     2        + FX2(-1)*sigmaij(5,-1)
     3        + FX2(-4)*sigmaij(5,-4))
     1        + (FX1(-3)+FX1(3))*
     1        ( FX2(2)*sigmaij(-1,2)
     2        + FX2(3)*sigmaij(-1,3)
     3        + FX2(5)*sigmaij(-1,5)
     1        + FX2(2)*sigmaij(-4,2)
     2        + FX2(3)*sigmaij(-4,3)
     3        + FX2(5)*sigmaij(-4,5)
     2        + FX2(-1)*sigmaij(2,-1)
     3        + FX2(-4)*sigmaij(2,-4)
     2        + FX2(-1)*sigmaij(5,-1)
     3        + FX2(-4)*sigmaij(5,-4))
     1  + (FX1(-4)+FX1(4))*
     1        ( FX2(2)*sigmaij(-1,2)
     2        + FX2(3)*sigmaij(-1,3)
     3        + FX2(5)*sigmaij(-1,5)
     2        + FX2(-1)*sigmaij(2,-1)
     3        + FX2(-4)*sigmaij(2,-4)
     2        + FX2(-1)*sigmaij(3,-1)
     3        + FX2(-4)*sigmaij(3,-4)
     2        + FX2(-1)*sigmaij(5,-1)
     3        + FX2(-4)*sigmaij(5,-4))
     1  + (FX1(-5)+FX1(5))*
     1        ( FX2(2)*sigmaij(-1,2)
     2        + FX2(3)*sigmaij(-1,3)
     3        + FX2(5)*sigmaij(-1,5)
     1        + FX2(2)*sigmaij(-4,2)
     2        + FX2(3)*sigmaij(-4,3)
     3        + FX2(5)*sigmaij(-4,5)
     2        + FX2(-1)*sigmaij(2,-1)
     3        + FX2(-4)*sigmaij(2,-4)
     2        + FX2(-1)*sigmaij(3,-1)
     3        + FX2(-4)*sigmaij(3,-4))

         QQBN_4 = 
     1     (FX2(-1)+FX2(1))*
     1       (FX1(-1)*sigmaij(-1,2)
     2      + FX1(-1)*sigmaij(-1,3)
     3      + FX1(-1)*sigmaij(-1,5)
     1      + FX1(-4)*sigmaij(-4,2)
     2      + FX1(-4)*sigmaij(-4,3)
     3      + FX1(-4)*sigmaij(-4,5)
     3      + FX1(2)*sigmaij(2,-4)
     3      + FX1(3)*sigmaij(3,-4)
     3      + FX1(5)*sigmaij(5,-4))
     1  + (FX2(-2)+FX2(2))*
     2       (FX1(-1)*sigmaij(-1,3)
     3      + FX1(-1)*sigmaij(-1,5)
     2      + FX1(-4)*sigmaij(-4,3)
     3      + FX1(-4)*sigmaij(-4,5)
     2      + FX1(2)*sigmaij(2,-1)
     3      + FX1(2)*sigmaij(2,-4)
     2      + FX1(3)*sigmaij(3,-1)
     3      + FX1(3)*sigmaij(3,-4)
     2      + FX1(5)*sigmaij(5,-1)
     3      + FX1(5)*sigmaij(5,-4))
     1  + (FX2(-3)+FX2(3))*
     1       (FX1(-1)*sigmaij(-1,2)
     3      + FX1(-1)*sigmaij(-1,5)
     1      + FX1(-4)*sigmaij(-4,2)
     3      + FX1(-4)*sigmaij(-4,5)
     2      + FX1(2)*sigmaij(2,-1)
     3      + FX1(2)*sigmaij(2,-4)
     2      + FX1(3)*sigmaij(3,-1)
     3      + FX1(3)*sigmaij(3,-4)
     2      + FX1(5)*sigmaij(5,-1)
     3      + FX1(5)*sigmaij(5,-4))
     1  + (FX2(-4)+FX2(4))*
     1       (FX1(-1)*sigmaij(-1,2)
     2      + FX1(-1)*sigmaij(-1,3)
     3      + FX1(-1)*sigmaij(-1,5)
     1      + FX1(-4)*sigmaij(-4,2)
     2      + FX1(-4)*sigmaij(-4,3)
     3      + FX1(-4)*sigmaij(-4,5)
     2      + FX1(2)*sigmaij(2,-1)
     2      + FX1(3)*sigmaij(3,-1)
     2      + FX1(5)*sigmaij(5,-1))
     1  + (FX2(-5)+FX2(5))*
     1       (FX1(-1)*sigmaij(-1,2)
     2      + FX1(-1)*sigmaij(-1,3)
     1      + FX1(-4)*sigmaij(-4,2)
     2      + FX1(-4)*sigmaij(-4,3)
     2      + FX1(2)*sigmaij(2,-1)
     3      + FX1(2)*sigmaij(2,-4)
     2      + FX1(3)*sigmaij(3,-1)
     3      + FX1(3)*sigmaij(3,-4)
     2      + FX1(5)*sigmaij(5,-1)
     3      + FX1(5)*sigmaij(5,-4))
         

      endif
c     **************************************
      endif
      if (mod.eq.2) then
      sHCRN(1,-1) = sGGN(1,-1)*Hgg+sQGN_1(1,-1)*Hqg_1+sQGN_2(1,-1)*Hqg_2
     1 	   + sQQBN_1(1,-1)*Hqqb+sQQBN_2(1,-1)*Hqq
     1 	   + sQQBN_3(1,-1)*Hqqp_1+sQQBN_4(1,-1)*Hqqp_2
      sHCRN(2,-2) = sGGN(2,-2)*Hgg+sQGN_1(2,-2)*Hqg_1+sQGN_2(2,-2)*Hqg_2
     1      + sQQBN_1(2,-2)*Hqqb+sQQBN_2(2,-2)*Hqq
     1 	   + sQQBN_3(2,-2)*Hqqp_1+sQQBN_4(2,-2)*Hqqp_2
      sHCRN(-1,1) = sGGN(-1,1)*Hgg+sQGN_1(-1,1)*Hqg_1+sQGN_2(-1,1)*Hqg_2
     1      + sQQBN_1(-1,1)*Hqqb+sQQBN_2(-1,1)*Hqq
     1 	   + sQQBN_3(-1,1)*Hqqp_1+sQQBN_4(-1,1)*Hqqp_2
      sHCRN(-2,2) = sGGN(-2,2)*Hgg+sQGN_1(-2,2)*Hqg_1+sQGN_2(-2,2)*Hqg_2
     1      + sQQBN_1(-2,2)*Hqqb+sQQBN_2(-2,2)*Hqq
     1 	   + sQQBN_3(-2,2)*Hqqp_1+sQQBN_4(-2,2)*Hqqp_2

c         HCRN = sHCRN(1,-1)*sigmaij(1,-1)
c     1    + sHCRN(2,-2)*sigmaij(2,-2)
c     1    + sHCRN(-1,1)*sigmaij(-1,1)
c     1    + sHCRN(-2,2)*sigmaij(-2,2)
      else
         HCRN = GGN*Hgg + FX1(0)*QGN_1*Hqg_1 + FX2(0)*QGN_2*Hqg_2
     1        + QQBN_1*Hqqb + QQBN_2*Hqq + QQBN_3*Hqqp_1 + QQBN_4*Hqqp_2
      endif
c      print *,'fortran',I1,I2,GGN,Hgg
      RETURN
      END


C
C PROVIDES MOMENTS OF DENSITIES AT A GIVEN SCALE (INCLUDES EVOLUTION)
C AND ANOMALOUS DIMENSIONS
C EVERYTHING IN MELLIN SPACE
C
c     input: I point in z grid, ISIGN (+ or -), IBEAM (1 or 2), SCALE2 (b), alps = alpqf (mass dependent)
c     output: FN, alpq = alps * alphasl(scale2)
c     dependence: FN(b, mass, I)

c ***************************************
c Can completely cache the output FN in FN(I,IFIT,ISIG), and calculate in the init? No because of b
c reno2 is cached and looped only on I, before entering the I1 I2 doube loop into cfx1 cfx2p cfx2m
       subroutine reno2 (fx, I, ISIGN, IBEAM, ALPS, SCALE2)
     
       IMPLICIT DOUBLE COMPLEX (A - Z)
       DIMENSION fx(-5:5)

       DOUBLE PRECISION ALPS
       DOUBLE PRECISION S

       DOUBLE PRECISION XL, XL1, SALP
       DOUBLE COMPLEX alpq,ALPr
       COMMON/alphasldata/XL,XL1,SALP,alpq,ALPr

       DOUBLE PRECISION  ALPQF, ALPQR
       COMMON / COUPL  / ALPQF, ALPQR

      INTEGER NMX
      PARAMETER (NMX = 512)

       INTEGER F, nnF, NAORD, I, ISIGN, IBEAM, IFIT, SIG
         COMPLEX*16 UVP(NMX,30),DVP(NMX,30),USP(NMX,30),DSP(NMX,30),
     .         SSP(NMX,30),GLP(NMX,30),CHP(NMX,30),BOP(NMX,30)
       COMPLEX*16 UVM(NMX,30),DVM(NMX,30),USM(NMX,30),DSM(NMX,30),
     .           SSM(NMX,30),GLM(NMX,30),CHM(NMX,30),BOM(NMX,30)
       COMPLEX*16 UVP2(NMX,30),DVP2(NMX,30),USP2(NMX,30),DSP2(NMX,30),
     .           SSP2(NMX,30),GLP2(NMX,30),CHP2(NMX,30),BOP2(NMX,30)
       COMPLEX*16 UVM2(NMX,30),DVM2(NMX,30),USM2(NMX,30),DSM2(NMX,30),
     .           SSM2(NMX,30),GLM2(NMX,30),CHM2(NMX,30),BOM2(NMX,30)
          
       COMPLEX*16 QQIP(NMX),QGFP(NMX), GQIP(NMX), GGIP(NMX), GGFP(NMX),
     1            NS1MIP(NMX), NS1PIP(NMX), NS1FP(NMX),QQ1FP(NMX), 
     2       QG1FP(NMX), GQ1IP(NMX), GQ1FP(NMX), GG1IP(NMX), GG1FP(NMX) 
       COMPLEX*16 QQIM(NMX),QGFM(NMX), GQIM(NMX), GGIM(NMX), GGFM(NMX),
     1           NS1MIM(NMX), NS1PIM(NMX), NS1FM(NMX),QQ1FM(NMX), 
     2       QG1FM(NMX), GQ1IM(NMX), GQ1FM(NMX), GG1IM(NMX), GG1FM(NMX)
       COMPLEX*16 C2qgMp(NMX),C2NSqqMp(NMX),C2SqqbMp(NMX),
     .            C2NSqqbMp(NMX)
       COMPLEX*16 C2qgMm(NMX),C2NSqqMm(NMX),C2SqqbMm(NMX),
     .            C2NSqqbMm(NMX)
      
       COMMON / ANOMP/QQIp, QGFp, GQIp, GGIp, GGFp, NS1MIp, NS1PIp, 
     1          NS1Fp, QQ1Fp, QG1Fp, GQ1Ip, GQ1Fp, GG1Ip, GG1Fp
       COMMON / ANOMM/QQIm, QGFm, GQIm, GGIm, GGFm, NS1MIm, NS1PIm, 
     1          NS1Fm, QQ1Fm, QG1Fm, GQ1Im, GQ1Fm, GG1Im, GG1Fm
       COMMON / DISTP1/ UVP,DVP,USP,DSP,SSP,GLP,CHP,BOP
       COMMON / DISTM1/ UVM,DVM,USM,DSM,SSM,GLM,CHM,BOM
       COMMON / DISTP2/ UVP2,DVP2,USP2,DSP2,SSP2,GLP2,CHP2,BOP2
       COMMON / DISTM2/ UVM2,DVM2,USM2,DSM2,SSM2,GLM2,CHM2,BOM2
       COMMON / NAORD / NAORD
       COMMON/NFLAVORS/nnF
       COMMON/COEFF/BETA0N,BETA1N,BETA2N,A1GN,A2GN,A3GN,B1GN,B2GN,H1G
       COMMON/ IFIT/ IFIT    
       COMMON / H2COEF /C2qgMp,C2NSqqMp,C2SqqbMp,C2NSqqbMp,
     .                  C2qgMm,C2NSqqMm,C2SqqbMm,C2NSqqbMm
      COMMON/exponent/aexp,aexpB

c     Cached values
      DOUBLE COMPLEX cgamma1qq(136,2),cgamma1qg(136,2),cgamma1gq(136,2),
     1     cgamma1gg(136,2),cgamma2qq(136,2),cgamma2qqV(136,2),
     2     cgamma2qqbV(136,2),cgamma2qqS(136,2),cgamma2qqbS(136,2),
     3     cgamma2qg(136,2),cgamma2gq(136,2),cgamma2gg(136,2)
      DOUBLE COMPLEX cC1QQ(136,2),cC1QG(136,2),cC1GQ(136,2),cC1GG
      DOUBLE COMPLEX cans(136,2),cam(136,2), cap(136,2), cal(136,2), 
     1     cbe(136,2),cab(136,2)
      DOUBLE COMPLEX crmin(136,2),crplus(136,2),
     1     crqq(136,2), crqg(136,2), crgq(136,2), crgg(136,2)
      DOUBLE COMPLEX cAC(136,2),cNMP(136,2),cNPM(136,2),cDMQQ(136,2),
     1     cDMQG(136,2),cDMGQ(136,2),cDMGG(136,2),cDPQQ(136,2),
     2     cDPQG(136,2),cDPGQ(136,2),cDPGG(136,2),cRMMQQ(136,2),
     3     cRMMQG(136,2),cRMMGQ(136,2),
     3     cRMMGG(136,2),cRMPQQ(136,2),cRMPQG(136,2),cRMPGQ(136,2),
     4     cRMPGG(136,2),cRPMQQ(136,2),cRPMQG(136,2),cRPMGQ(136,2),
     5     cRPMGG(136,2),cRPPQQ(136,2),cRPPQG(136,2),cRPPGQ(136,2),
     6     cRPPGG(136,2)
      COMPLEX*16 cC2qgM(136,2),cC2NSqqM(136,2),cC2SqqbM(136,2),
     1     cC2NSqqbM(136,2)
      COMMON /ANOMCACHE/cans,cam,cap,cal,cbe,cab,crmin,crplus,crqq,
     1     crqg,crgq,crgg,
     2     cAC,cNMP,cNPM,cDMQQ,
     3     cDMQG,cDMGQ,cDMGG,cDPQQ,cDPQG,
     4     cDPGQ,cDPGG,cRMMQQ,cRMMQG,cRMMGQ,
     5     cRMMGG,cRMPQQ,cRMPQG,cRMPGQ,
     6     cRMPGG,cRPMQQ,cRPMQG,cRPMGQ,
     7     cRPMGG,cRPPQQ,cRPPQG,cRPPGQ,
     8     cRPPGG,cC1QQ,cC1QG,cC1GQ,cC1GG,
     9     cgamma1qq,cgamma1qg,cgamma1gq,
     1     cgamma1gg,cgamma2qq,cgamma2qqV,
     2     cgamma2qqbV,cgamma2qqS,cgamma2qqbS,
     3     cgamma2qg,cgamma2gq,cgamma2gg,
     4     cC2qgM,cC2NSqqM,cC2SqqbM,cC2NSqqbM
     
       include 'const.h'

c**************************************
c     IFIT is only rapidity dependent, I is the point in the z grid of the gaussian quadrature loop
       IF (IBEAM.EQ.1) THEN

       IF (ISIGN.EQ.1) THEN
        UVI = UVp(I,IFIT)
        DVI = DVp(I,IFIT)
        USI = USp(I,IFIT)
        DSI = DSp(I,IFIT)
        SSI = SSp(I,IFIT)
        GLI = GLp(I,IFIT)
        CHI = CHp(I,IFIT)
        BOI = BOp(I,IFIT)

        ELSE 
        UVI = UVm(I,IFIT)
        DVI = DVm(I,IFIT)
        USI = USm(I,IFIT)
        DSI = DSm(I,IFIT)
        SSI = SSm(I,IFIT)
        GLI = GLm(I,IFIT)
        CHI = CHm(I,IFIT)
        BOI = BOm(I,IFIT)
        
         ENDIF      

       ELSEIF (IBEAM.EQ.2) THEN
       IF (ISIGN.EQ.1) THEN
        UVI = UVp2(I,IFIT)
        DVI = DVp2(I,IFIT)
        USI = USp2(I,IFIT)
        DSI = DSp2(I,IFIT)
        SSI = SSp2(I,IFIT)
        GLI = GLp2(I,IFIT)
        CHI = CHp2(I,IFIT)
        BOI = BOp2(I,IFIT)

        ELSE 
        UVI = UVm2(I,IFIT)
        DVI = DVm2(I,IFIT)
        USI = USm2(I,IFIT)
        DSI = DSm2(I,IFIT)
        SSI = SSm2(I,IFIT)
        GLI = GLm2(I,IFIT)
        CHI = CHm2(I,IFIT)
        BOI = BOm2(I,IFIT)
        
         ENDIF      
       ELSE
        write(6,*)"wrong IBEAM in RENOSUS= ",IBEAM
        STOP
       ENDIF      
             
c ***********************
c     this part can be precomputed
c     also add IFIT dependence
      UVN = UVI
      DVN = DVI
      NS3N = UVI + 2*USI - DVI - 2.*DSI
      NS8N = UVI + 2*USI + DVI + 2.*DSI - 4.*SSI
      GLN = GLI
      SIN = UVI + DVI + 2*USI + 2*DSI
     .     + 2*SSI + 2*CHI + 2*BOI 
      NS15N = UVI + DVI + 2*USI + 2*DSI + 2*SSI - 6*CHI
      NS24N = UVI + DVI + 2*USI + 2*DSI +
     .     2*SSI + 2*CHI - 8*BOI
      NS35N = SIN
 
      IF (nnf.eq.3) then
        SIN = UVI + DVI + 2*USI + 2*DSI + 2*SSI 
        NS15N = SIN
        NS24N = SIN
        NS35N = SIN
      ENDIF
      F = 5
      SG = SIN
      GL = GLN
c ***********************
c**************************************

c**************************************
c retrieve cached values
      S=SALP

c retrieved values cached in cacheanom
      SIG = (-ISIGN+1)/2+1
      ANS=cANS(I,SIG)
      AM=cAM(I,SIG)
      AP=cAP(I,SIG)
      AL=cAL(I,SIG)
      BE=cBE(I,SIG)
      AB=cAB(I,SIG)
      RMIN=cRMIN(I,SIG)
      RPLUS=cRPLUS(I,SIG)
      AC  = cAC(I,SIG)
      RMMQQ = cRMMQQ(I,SIG)
      RMMQG = cRMMQG(I,SIG)
      RMMGQ = cRMMGQ(I,SIG)
      RMMGG = cRMMGG(I,SIG)
      RMPQQ = cRMPQQ(I,SIG)
      RMPQG = cRMPQG(I,SIG)
      RMPGQ = cRMPGQ(I,SIG)
      RMPGG = cRMPGG(I,SIG)
      RPMQQ = cRPMQQ(I,SIG)
      RPMQG = cRPMQG(I,SIG)
      RPMGQ = cRPMGQ(I,SIG)
      RPMGG = cRPMGG(I,SIG)
      RPPQQ = cRPPQQ(I,SIG)
      RPPQG = cRPPQG(I,SIG)
      RPPGQ = cRPPGQ(I,SIG)
      RPPGG = cRPPGG(I,SIG)
c**************************************

c**************************************
c     b-dependence
       ENS = EXP (-ANS*S)
C
       EM  = EXP (-AM*S)
       EP  = EXP (-AP*S)
       EMP = EM/EP
       EPM = EP/EM
c**************************************

c**************************************
c     b-dependence
C...EVOLUTION OF LIGHT PARTON DENSITIES
         UVN  = UVN  * ENS * (1.+  ALPr * XL1 * RMIN)
         DVN  = DVN  * ENS * (1.+  ALPr * XL1 * RMIN)
         NS3N = NS3N * ENS * (1.+  ALPr * XL1 * RPLUS)
         NS8N = NS8N * ENS * (1.+  ALPr * XL1 * RPLUS)
c**************************************
C
c**************************************
c     b-dependence
       SIN = EM * ((AL + ALPr * (RMMQQ * XL1 + RMPQQ * (EPM-XL)))* SG
     1           + (BE + ALPr * (RMMQG * XL1 + RMPQG * (EPM-XL))) * GL)
     2     + EP * ((AC + ALPr * (RPPQQ * XL1 + RPMQQ * (EMP-XL))) * SG
     3           +(-BE + ALPr * (RPPQG * XL1 + RPMQG * (EMP-XL))) * GL)
       GLN = EM * ((AB + ALPr * (RMMGQ * XL1 + RMPGQ * (EPM-XL))) * SG
     1           + (AC + ALPr * (RMMGG * XL1 + RMPGG * (EPM-XL))) * GL)
     2     + EP *((-AB + ALPr * (RPPGQ * XL1 + RPMGQ * (EMP-XL))) * SG
     3           + (AL + ALPr * (RPPGG * XL1 + RPMGG * (EMP-XL))) * GL)
C
      NS15N = NS15N * ENS * (1.+  ALPr * XL1 * RPLUS)
      NS24N = NS24N * ENS * (1.+  ALPr * XL1 * RPLUS)
      NS35N = SIN
c**************************************

c**************************************
c     b-dependence
C...  FLAVOUR DECOMPOSITION OF THE QUARK SEA :
      SSN = (10.* SIN + 2.* NS35N + 3.* NS24N + 5.* NS15N - 20.* NS8N)
     1     / 120.
      DSN = (10.* SIN + 2.* NS35N + 3.* NS24N + 5.* NS15N + 10.* NS8N
     1     - 30.* NS3N - 60.* DVN) / 120.
      USN = (10.* SIN + 2.* NS35N + 3.* NS24N + 5.* NS15N + 10.* NS8N
     1     + 30.* NS3N - 60.* UVN) / 120.
      CHN = (UVN + DVN + 2*USN + 2*DSN + 2*SSN - NS15N)/6.
      BON = (UVN + DVN + 2*USN + 2*DSN + 2*SSN + 2*CHN - NS24N)/8.
      
      if (nnf.eq.3) then        !GRV
         SSN= (20.* SIN - 20.* NS8N)/120.
         DSN = (20.* SIN + 10.* NS8N - 30.* NS3N - 60.* DVN) / 120.
         USN = (20.* SIN + 10.* NS8N + 30.* NS3N - 60.* UVN) / 120.
         CHN=CMPLX(0D0,0D0)
         BON=CMPLX(0D0,0D0)
      endif
c**************************************

*...  OUTPUT  
*  U   UB   D   DB   S   C    B    G 
*  1    2   3   4    5   6    7    8 
c************************************
c     b dependence and I dependence and ISIGN
      fx(0) = GLN
      fx(1) = UVN + USN
      fx(-1) = USN
      fx(2) = DVN + DSN
      fx(-2) = DSN
      fx(3) = SSN
      fx(-3) = SSN
      if (nnf.ge.4) then 
         fx(4) = CHN
         fx(-4) = CHN
      else 
         fx(4) = 0D0
         fx(-4) = 0D0
      endif
      if(nf.ge.5) then 
         fx(5) = BON
         fx(-5) = BON
      else 
         fx(5) = 0D0
         fx(-5) = 0D0
      endif
      return
      end

C
C
C...BETA FUNCTION FOR COMPLEX ARGUMENT :
       DOUBLE COMPLEX FUNCTION CBETA (Z1, Z2)
       IMPLICIT DOUBLE COMPLEX (A - Z)
       LNGAM (X) = (X - 0.5) * LOG (X) - X + 0.91893853 + 1./(12.* X)
     1              * (1.- 1./(30.* X*X) * (1.- 1./(3.5 * X*X)
     2              * (1.- 4./(3.* X*X))))
       SUB = CMPLX (0., 0.)
       ZZ1 = Z1
  1    CONTINUE
       IF ( DREAL (ZZ1) .LT. 15.) THEN
          SUB = SUB + LOG ((ZZ1+Z2) / ZZ1)
          ZZ1 = ZZ1 + 1.
          GOTO 1
       END IF
       ZZ2 = Z2
  2    CONTINUE
       IF ( DREAL (ZZ2) .LT. 15.) THEN
          SUB = SUB + LOG ((ZZ1+ZZ2) / ZZ2)
          ZZ2 = ZZ2 + 1.
          GOTO 2
       END IF
       LG1 = LNGAM (ZZ1)
       LG2 = LNGAM (ZZ2)
       LG12 = LNGAM (ZZ1 + ZZ2)
       CBETA = EXP (LG1 + LG2 - LG12 + SUB)
       RETURN
       END
*



C     Computes the moments of pdfs for beam 1
      SUBROUTINE F0MOMENTS1(N,UV,DV,US,DS,SS,GL,CH,BO)
      implicit DOUBLE COMPLEX (A-Z)
      integer flagi,IK,NFITMAX
      double precision A1UV,A2UV,A3UV,A4UV,A5UV,A6UV,A7UV,A8UV
      double precision A1DV,A2DV,A3DV,A4DV,A5DV,A6DV,A7DV,A8DV
      double precision A1US,A2US,A3US,A4US,A5US,A6US,A7US,A8US 
      double precision A1DS,A2DS,A3DS,A4DS,A5DS,A6DS,A7DS,A8DS
      double precision A1SS,A2SS,A3SS,A4SS,A5SS,A6SS,A7SS,A8SS
      double precision A1GL,A2GL,A3GL,A4GL,A5GL,A6GL,A7GL,A8GL
      double precision A1CH,A2CH,A3CH,A4CH,A5CH,A6CH,A7CH,A8CH
      double precision A1BO,A2BO,A3BO,A4BO,A5BO,A6BO,A7BO,A8BO
      dimension UV(30),DV(30),US(30),DS(30),SS(30),GL(30),CH(30),BO(30)
c
       dimension A1UV(30),A2UV(30),A3UV(30),A4UV(30),
     .       A5UV(30),A6UV(30),A7UV(30),A8UV(30)
       dimension A1DV(30),A2DV(30),A3DV(30),A4DV(30),
     .       A5DV(30),A6DV(30),A7DV(30),A8DV(30)
       dimension A1US(30),A2US(30),A3US(30),A4US(30),
     .       A5US(30),A6US(30),A7US(30),A8US(30)
       dimension A1DS(30),A2DS(30),A3DS(30),A4DS(30),
     .       A5DS(30),A6DS(30),A7DS(30),A8DS(30)
       dimension A1SS(30),A2SS(30),A3SS(30),A4SS(30),
     .       A5SS(30),A6SS(30),A7SS(30),A8SS(30)
       dimension A1GL(30),A2GL(30),A3GL(30),A4GL(30),
     .       A5GL(30),A6GL(30),A7GL(30),A8GL(30)
       dimension A1CH(30),A2CH(30),A3CH(30),A4CH(30),
     .       A5CH(30),A6CH(30),A7CH(30),A8CH(30)
       dimension A1BO(30),A2BO(30),A3BO(30),A4BO(30),
     .       A5BO(30),A6BO(30),A7BO(30),A8BO(30)
c
      double precision GE,ZETA2,ZETA3,PI,aa 
      COMMON/ CUV1/ A1UV,A2UV,A3UV,A4UV,A5UV,A6UV,A7UV,A8UV
      COMMON/ CDV1/ A1DV,A2DV,A3DV,A4DV,A5DV,A6DV,A7DV,A8DV
      COMMON/ CUS1/ A1US,A2US,A3US,A4US,A5US,A6US,A7US,A8US
      COMMON/ CDS1/ A1DS,A2DS,A3DS,A4DS,A5DS,A6DS,A7DS,A8DS
      COMMON/ CSS1/ A1SS,A2SS,A3SS,A4SS,A5SS,A6SS,A7SS,A8SS
      COMMON/ CGL1/ A1GL,A2GL,A3GL,A4GL,A5GL,A6GL,A7GL,A8GL
      COMMON/ CCH1/ A1CH,A2CH,A3CH,A4CH,A5CH,A6CH,A7CH,A8CH
      COMMON/ CBO1/ A1BO,A2BO,A3BO,A4BO,A5BO,A6BO,A7BO,A8BO
      COMMON/expp/aa
      COMMON/NFITMAX/NFITMAX
 
      data GE/ 0.577216d0 /
      data ZETA2/1.64493d0 / 
      data ZETA3/1.20206d0 /
      data PI/3.1415926536d0 /
      AAA=CMPLX(aa,0.d0)

      DO IK=1,NFITMAX
C     U-VALENCE
      CA1=CMPLX(A1UV(ik),0.d0)
      CA2=CMPLX(A2UV(ik),0.d0)
      CA3=CMPLX(A3UV(ik),0.d0)
      CA4=CMPLX(A4UV(ik),0.d0)
      CA5=CMPLX(A5UV(ik),0.d0)
      CA6=CMPLX(A6UV(ik),0.d0)
      CA7=CMPLX(A7UV(ik),0.d0)
      CA8=CMPLX(A8UV(ik),0.d0)
      UV(ik)=CA1*CBETA(CA2-1+N,CA3+1) + CA1*CA4*CBETA(CA2+N,CA3+1)+
     .     CA1*CA5*CBETA(CA2-1+0.5+N,CA3+1)
     .     + CA1*CA6*CBETA(CA2-1+1.5+N,CA3+1) 
     .     + CA1*CA7*CBETA(CA2-1+2+N,CA3+1)
     .     + CA1*CA8*CBETA(CA2-1+AAA+N,CA3+1)
C     D-VALENCE
      CA1=CMPLX(A1DV(ik),0.d0)
      CA2=CMPLX(A2DV(ik),0.d0)
      CA3=CMPLX(A3DV(ik),0.d0)
      CA4=CMPLX(A4DV(ik),0.d0)
      CA5=CMPLX(A5DV(ik),0.d0)
      CA6=CMPLX(A6DV(ik),0.d0)
      CA7=CMPLX(A7DV(ik),0.d0)
      CA8=CMPLX(A8DV(ik),0.d0)
      DV(ik)=CA1*CBETA(CA2-1+N,CA3+1) + CA1*CA4*CBETA(CA2+N,CA3+1)+
     .     CA1*CA5*CBETA(CA2-1+0.5+N,CA3+1) 
     .     + CA1*CA6*CBETA(CA2-1+1.5+N,CA3+1) 
     .     + CA1*CA7*CBETA(CA2-1+2+N,CA3+1)
     .     + CA1*CA8*CBETA(CA2-1+AAA+N,CA3+1)
C     U-BAR
      CA1=CMPLX(A1US(ik),0.d0)
      CA2=CMPLX(A2US(ik),0.d0)
      CA3=CMPLX(A3US(ik),0.d0)
      CA4=CMPLX(A4US(ik),0.d0)
      CA5=CMPLX(A5US(ik),0.d0)
      CA6=CMPLX(A6US(ik),0.d0)
      CA7=CMPLX(A7US(ik),0.d0)
      CA8=CMPLX(A8US(ik),0.d0)
      US(ik)=CA1*CBETA(CA2-1+N,CA3+1) + CA1*CA4*CBETA(CA2+N,CA3+1)+
     .     CA1*CA5*CBETA(CA2-1+0.5+N,CA3+1) 
     .     + CA1*CA6*CBETA(CA2-1+1.5+N,CA3+1) 
     .     + CA1*CA7*CBETA(CA2-1+2+N,CA3+1)
     .     + CA1*CA8*CBETA(CA2-1+AAA+N,CA3+1)
C     D-BAR
      CA1=CMPLX(A1DS(ik),0.d0)
      CA2=CMPLX(A2DS(ik),0.d0)
      CA3=CMPLX(A3DS(ik),0.d0)
      CA4=CMPLX(A4DS(ik),0.d0)
      CA5=CMPLX(A5DS(ik),0.d0)
      CA6=CMPLX(A6DS(ik),0.d0)
      CA7=CMPLX(A7DS(ik),0.d0)
      CA8=CMPLX(A8DS(ik),0.d0)
      DS(ik)=CA1*CBETA(CA2-1+N,CA3+1) + CA1*CA4*CBETA(CA2+N,CA3+1)+
     .     CA1*CA5*CBETA(CA2-1+0.5+N,CA3+1) 
     .     + CA1*CA6*CBETA(CA2-1+1.5+N,CA3+1) 
     .     + CA1*CA7*CBETA(CA2-1+2+N,CA3+1)
     .     + CA1*CA8*CBETA(CA2-1+AAA+N,CA3+1)
C     S-BAR
      CA1=CMPLX(A1SS(ik),0.d0)
      CA2=CMPLX(A2SS(ik),0.d0)
      CA3=CMPLX(A3SS(ik),0.d0)
      CA4=CMPLX(A4SS(ik),0.d0)
      CA5=CMPLX(A5SS(ik),0.d0)
      CA6=CMPLX(A6SS(ik),0.d0)
      CA7=CMPLX(A7SS(ik),0.d0)
      CA8=CMPLX(A8SS(ik),0.d0)
      SS(ik)=CA1*CBETA(CA2-1+N,CA3+1) + CA1*CA4*CBETA(CA2+N,CA3+1)+
     .     CA1*CA5*CBETA(CA2-1+0.5+N,CA3+1)
     .     + CA1*CA6*CBETA(CA2-1+1.5+N,CA3+1) 
     .     + CA1*CA7*CBETA(CA2-1+2+N,CA3+1)
     .     + CA1*CA8*CBETA(CA2-1+AAA+N,CA3+1)
C     GLUON
      CA1=CMPLX(A1GL(ik),0.d0)
      CA2=CMPLX(A2GL(ik),0.d0)
      CA3=CMPLX(A3GL(ik),0.d0)
      CA4=CMPLX(A4GL(ik),0.d0)
      CA5=CMPLX(A5GL(ik),0.d0)
      CA6=CMPLX(A6GL(ik),0.d0)
      CA7=CMPLX(A7GL(ik),0.d0)
      CA8=CMPLX(A8GL(ik),0.d0)
      GL(ik)=CA1*CBETA(CA2-1+N,CA3+1) + CA1*CA4*CBETA(CA2+N,CA3+1)+
     .     CA1*CA5*CBETA(CA2-1+0.5+N,CA3+1)
     .     + CA1*CA6*CBETA(CA2-1+1.5+N,CA3+1) 
     .     + CA1*CA7*CBETA(CA2-1+2+N,CA3+1)
     .     + CA1*CA8*CBETA(CA2-1+AAA+N,CA3+1)
C     CHARM
      CA1=CMPLX(A1CH(ik),0.d0)
      CA2=CMPLX(A2CH(ik),0.d0)
      CA3=CMPLX(A3CH(ik),0.d0)
      CA4=CMPLX(A4CH(ik),0.d0)
      CA5=CMPLX(A5CH(ik),0.d0)
      CA6=CMPLX(A6CH(ik),0.d0)
      CA7=CMPLX(A7CH(ik),0.d0)
      CA8=CMPLX(A8CH(ik),0.d0)
      CH(ik)=CA1*CBETA(CA2-1+N,CA3+1) + CA1*CA4*CBETA(CA2+N,CA3+1)+
     .     CA1*CA5*CBETA(CA2-1+0.5+N,CA3+1)
     .     + CA1*CA6*CBETA(CA2-1+1.5+N,CA3+1) 
     .     + CA1*CA7*CBETA(CA2-1+2+N,CA3+1)
     .     + CA1*CA8*CBETA(CA2-1+AAA+N,CA3+1)
C     BOTTOM
      CA1=CMPLX(A1BO(ik),0.d0)
      CA2=CMPLX(A2BO(ik),0.d0)
      CA3=CMPLX(A3BO(ik),0.d0)
      CA4=CMPLX(A4BO(ik),0.d0)
      CA5=CMPLX(A5BO(ik),0.d0)
      CA6=CMPLX(A6BO(ik),0.d0)
      CA7=CMPLX(A7BO(ik),0.d0)
      CA8=CMPLX(A8BO(ik),0.d0)
      BO(ik)=CA1*CBETA(CA2-1+N,CA3+1) + CA1*CA4*CBETA(CA2+N,CA3+1)+
     .     CA1*CA5*CBETA(CA2-1+0.5+N,CA3+1)
     .     + CA1*CA6*CBETA(CA2-1+1.5+N,CA3+1) 
     .     + CA1*CA7*CBETA(CA2-1+2+N,CA3+1)
     .     + CA1*CA8*CBETA(CA2-1+AAA+N,CA3+1)
      ENDDO
      RETURN
      END      
      
      
      
C     Computes the moments of pdfs for beam 2
      SUBROUTINE F0MOMENTS2(N,UV,DV,US,DS,SS,GL,CH,BO)
      implicit DOUBLE COMPLEX (A-Z)
      integer flagi,IK,NFITMAX
      double precision A1UV,A2UV,A3UV,A4UV,A5UV,A6UV,A7UV,A8UV
      double precision A1DV,A2DV,A3DV,A4DV,A5DV,A6DV,A7DV,A8DV
      double precision A1US,A2US,A3US,A4US,A5US,A6US,A7US,A8US 
      double precision A1DS,A2DS,A3DS,A4DS,A5DS,A6DS,A7DS,A8DS
      double precision A1SS,A2SS,A3SS,A4SS,A5SS,A6SS,A7SS,A8SS
      double precision A1GL,A2GL,A3GL,A4GL,A5GL,A6GL,A7GL,A8GL
      double precision A1CH,A2CH,A3CH,A4CH,A5CH,A6CH,A7CH,A8CH
      double precision A1BO,A2BO,A3BO,A4BO,A5BO,A6BO,A7BO,A8BO
      dimension UV(30),DV(30),US(30),DS(30),SS(30),GL(30),CH(30),BO(30)
c
       dimension A1UV(30),A2UV(30),A3UV(30),A4UV(30),
     .       A5UV(30),A6UV(30),A7UV(30),A8UV(30)
       dimension A1DV(30),A2DV(30),A3DV(30),A4DV(30),
     .       A5DV(30),A6DV(30),A7DV(30),A8DV(30)
       dimension A1US(30),A2US(30),A3US(30),A4US(30),
     .       A5US(30),A6US(30),A7US(30),A8US(30)
       dimension A1DS(30),A2DS(30),A3DS(30),A4DS(30),
     .       A5DS(30),A6DS(30),A7DS(30),A8DS(30)
       dimension A1SS(30),A2SS(30),A3SS(30),A4SS(30),
     .       A5SS(30),A6SS(30),A7SS(30),A8SS(30)
       dimension A1GL(30),A2GL(30),A3GL(30),A4GL(30),
     .       A5GL(30),A6GL(30),A7GL(30),A8GL(30)
       dimension A1CH(30),A2CH(30),A3CH(30),A4CH(30),
     .       A5CH(30),A6CH(30),A7CH(30),A8CH(30)
       dimension A1BO(30),A2BO(30),A3BO(30),A4BO(30),
     .       A5BO(30),A6BO(30),A7BO(30),A8BO(30)
c
      double precision GE,ZETA2,ZETA3,PI,aa 
      COMMON/ CUV2/ A1UV,A2UV,A3UV,A4UV,A5UV,A6UV,A7UV,A8UV
      COMMON/ CDV2/ A1DV,A2DV,A3DV,A4DV,A5DV,A6DV,A7DV,A8DV
      COMMON/ CUS2/ A1US,A2US,A3US,A4US,A5US,A6US,A7US,A8US
      COMMON/ CDS2/ A1DS,A2DS,A3DS,A4DS,A5DS,A6DS,A7DS,A8DS
      COMMON/ CSS2/ A1SS,A2SS,A3SS,A4SS,A5SS,A6SS,A7SS,A8SS
      COMMON/ CGL2/ A1GL,A2GL,A3GL,A4GL,A5GL,A6GL,A7GL,A8GL
      COMMON/ CCH2/ A1CH,A2CH,A3CH,A4CH,A5CH,A6CH,A7CH,A8CH
      COMMON/ CBO2/ A1BO,A2BO,A3BO,A4BO,A5BO,A6BO,A7BO,A8BO
      COMMON/expp/aa
      COMMON/NFITMAX/NFITMAX
 
      data GE/ 0.577216d0 /
      data ZETA2/1.64493d0 / 
      data ZETA3/1.20206d0 /
      data PI/3.1415926536d0 /
      AAA=CMPLX(aa,0.d0)

      DO IK=1,NFITMAX
C     U-VALENCE
      CA1=CMPLX(A1UV(ik),0.d0)
      CA2=CMPLX(A2UV(ik),0.d0)
      CA3=CMPLX(A3UV(ik),0.d0)
      CA4=CMPLX(A4UV(ik),0.d0)
      CA5=CMPLX(A5UV(ik),0.d0)
      CA6=CMPLX(A6UV(ik),0.d0)
      CA7=CMPLX(A7UV(ik),0.d0)
      CA8=CMPLX(A8UV(ik),0.d0)
      UV(ik)=CA1*CBETA(CA2-1+N,CA3+1) + CA1*CA4*CBETA(CA2+N,CA3+1)+
     .     CA1*CA5*CBETA(CA2-1+0.5+N,CA3+1)
     .     + CA1*CA6*CBETA(CA2-1+1.5+N,CA3+1) 
     .     + CA1*CA7*CBETA(CA2-1+2+N,CA3+1)
     .     + CA1*CA8*CBETA(CA2-1+AAA+N,CA3+1)
C     D-VALENCE
      CA1=CMPLX(A1DV(ik),0.d0)
      CA2=CMPLX(A2DV(ik),0.d0)
      CA3=CMPLX(A3DV(ik),0.d0)
      CA4=CMPLX(A4DV(ik),0.d0)
      CA5=CMPLX(A5DV(ik),0.d0)
      CA6=CMPLX(A6DV(ik),0.d0)
      CA7=CMPLX(A7DV(ik),0.d0)
      CA8=CMPLX(A8DV(ik),0.d0)
      DV(ik)=CA1*CBETA(CA2-1+N,CA3+1) + CA1*CA4*CBETA(CA2+N,CA3+1)+
     .     CA1*CA5*CBETA(CA2-1+0.5+N,CA3+1) 
     .     + CA1*CA6*CBETA(CA2-1+1.5+N,CA3+1) 
     .     + CA1*CA7*CBETA(CA2-1+2+N,CA3+1)
     .     + CA1*CA8*CBETA(CA2-1+AAA+N,CA3+1)
C     U-BAR
      CA1=CMPLX(A1US(ik),0.d0)
      CA2=CMPLX(A2US(ik),0.d0)
      CA3=CMPLX(A3US(ik),0.d0)
      CA4=CMPLX(A4US(ik),0.d0)
      CA5=CMPLX(A5US(ik),0.d0)
      CA6=CMPLX(A6US(ik),0.d0)
      CA7=CMPLX(A7US(ik),0.d0)
      CA8=CMPLX(A8US(ik),0.d0)
      US(ik)=CA1*CBETA(CA2-1+N,CA3+1) + CA1*CA4*CBETA(CA2+N,CA3+1)+
     .     CA1*CA5*CBETA(CA2-1+0.5+N,CA3+1) 
     .     + CA1*CA6*CBETA(CA2-1+1.5+N,CA3+1) 
     .     + CA1*CA7*CBETA(CA2-1+2+N,CA3+1)
     .     + CA1*CA8*CBETA(CA2-1+AAA+N,CA3+1)
C     D-BAR
      CA1=CMPLX(A1DS(ik),0.d0)
      CA2=CMPLX(A2DS(ik),0.d0)
      CA3=CMPLX(A3DS(ik),0.d0)
      CA4=CMPLX(A4DS(ik),0.d0)
      CA5=CMPLX(A5DS(ik),0.d0)
      CA6=CMPLX(A6DS(ik),0.d0)
      CA7=CMPLX(A7DS(ik),0.d0)
      CA8=CMPLX(A8DS(ik),0.d0)
      DS(ik)=CA1*CBETA(CA2-1+N,CA3+1) + CA1*CA4*CBETA(CA2+N,CA3+1)+
     .     CA1*CA5*CBETA(CA2-1+0.5+N,CA3+1) 
     .     + CA1*CA6*CBETA(CA2-1+1.5+N,CA3+1) 
     .     + CA1*CA7*CBETA(CA2-1+2+N,CA3+1)
     .     + CA1*CA8*CBETA(CA2-1+AAA+N,CA3+1)
C     S-BAR
      CA1=CMPLX(A1SS(ik),0.d0)
      CA2=CMPLX(A2SS(ik),0.d0)
      CA3=CMPLX(A3SS(ik),0.d0)
      CA4=CMPLX(A4SS(ik),0.d0)
      CA5=CMPLX(A5SS(ik),0.d0)
      CA6=CMPLX(A6SS(ik),0.d0)
      CA7=CMPLX(A7SS(ik),0.d0)
      CA8=CMPLX(A8SS(ik),0.d0)
      SS(ik)=CA1*CBETA(CA2-1+N,CA3+1) + CA1*CA4*CBETA(CA2+N,CA3+1)+
     .     CA1*CA5*CBETA(CA2-1+0.5+N,CA3+1)
     .     + CA1*CA6*CBETA(CA2-1+1.5+N,CA3+1) 
     .     + CA1*CA7*CBETA(CA2-1+2+N,CA3+1)
     .     + CA1*CA8*CBETA(CA2-1+AAA+N,CA3+1)
C     GLUON
      CA1=CMPLX(A1GL(ik),0.d0)
      CA2=CMPLX(A2GL(ik),0.d0)
      CA3=CMPLX(A3GL(ik),0.d0)
      CA4=CMPLX(A4GL(ik),0.d0)
      CA5=CMPLX(A5GL(ik),0.d0)
      CA6=CMPLX(A6GL(ik),0.d0)
      CA7=CMPLX(A7GL(ik),0.d0)
      CA8=CMPLX(A8GL(ik),0.d0)
      GL(ik)=CA1*CBETA(CA2-1+N,CA3+1) + CA1*CA4*CBETA(CA2+N,CA3+1)+
     .     CA1*CA5*CBETA(CA2-1+0.5+N,CA3+1)
     .     + CA1*CA6*CBETA(CA2-1+1.5+N,CA3+1) 
     .     + CA1*CA7*CBETA(CA2-1+2+N,CA3+1)
     .     + CA1*CA8*CBETA(CA2-1+AAA+N,CA3+1)
C     CHARM
      CA1=CMPLX(A1CH(ik),0.d0)
      CA2=CMPLX(A2CH(ik),0.d0)
      CA3=CMPLX(A3CH(ik),0.d0)
      CA4=CMPLX(A4CH(ik),0.d0)
      CA5=CMPLX(A5CH(ik),0.d0)
      CA6=CMPLX(A6CH(ik),0.d0)
      CA7=CMPLX(A7CH(ik),0.d0)
      CA8=CMPLX(A8CH(ik),0.d0)
      CH(ik)=CA1*CBETA(CA2-1+N,CA3+1) + CA1*CA4*CBETA(CA2+N,CA3+1)+
     .     CA1*CA5*CBETA(CA2-1+0.5+N,CA3+1)
     .     + CA1*CA6*CBETA(CA2-1+1.5+N,CA3+1) 
     .     + CA1*CA7*CBETA(CA2-1+2+N,CA3+1)
     .     + CA1*CA8*CBETA(CA2-1+AAA+N,CA3+1)
C     BOTTOM
      CA1=CMPLX(A1BO(ik),0.d0)
      CA2=CMPLX(A2BO(ik),0.d0)
      CA3=CMPLX(A3BO(ik),0.d0)
      CA4=CMPLX(A4BO(ik),0.d0)
      CA5=CMPLX(A5BO(ik),0.d0)
      CA6=CMPLX(A6BO(ik),0.d0)
      CA7=CMPLX(A7BO(ik),0.d0)
      CA8=CMPLX(A8BO(ik),0.d0)
      BO(ik)=CA1*CBETA(CA2-1+N,CA3+1) + CA1*CA4*CBETA(CA2+N,CA3+1)+
     .     CA1*CA5*CBETA(CA2-1+0.5+N,CA3+1)
     .     + CA1*CA6*CBETA(CA2-1+1.5+N,CA3+1) 
     .     + CA1*CA7*CBETA(CA2-1+2+N,CA3+1)
     .     + CA1*CA8*CBETA(CA2-1+AAA+N,CA3+1)
      ENDDO
      RETURN
      END      
      
    
      subroutine rapintegrals(ymin,ymax,m,nolepcuts)
      implicit none
      double precision ymin,ymax,m
      logical nolepcuts

      include 'quadrules.f'
      include 'gauss.f' 
      
c     store integration results in the common block, they are passed to initsigmacthy
      integer I1,I2
      complex *16 Ith0p(136,136)
      complex *16 Ith1p(136,136)
      complex *16 Ith2p(136,136)
      complex *16 Ith0m(136,136)
      complex *16 Ith1m(136,136)
      complex *16 Ith2m(136,136)
      common/ITHMOM/Ith0p,Ith1p,Ith2p,Ith0m,Ith1m,Ith2m
      
c nodes for the gaussian quadrature for the double Mellin inversion
      complex*16 CCp,CCm, Np(136),Nm(136)
      common /MOMS2/ Np,Nm,CCP,CCm
      
      double precision ax1,ax2,xx1,xx2
      integer nmax1,nmax2
      common /cxx/ax1,ax2,xx1,xx2,nmax1,nmax2

      include 'const.h'
      double precision C,CO,SI,ax
      COMMON /CONT/ C,CO,SI,AX

      double precision xsection
      external xsection
      double precision xsection2
      external xsection2
      double precision sqs
      common/energy/SQS 
      integer ih1,ih2
      common/collider/ih1,ih2
      include 'scales_inc.f'
      
      complex *16 yintp,yintm
      complex *16 fpy,fmy
      double precision cthmom0,cthmom1,cthmom2
      double precision y,xc,xm
      integer i,j
      double precision ya,yb

      complex *16 cfpm(136,136)
      complex *16 cfmm(136,136)

c     weights of the gaussian quadrature
      double precision  WN(136)
      COMMON /WEIGHTS2/ WN      

c     cached values from cacheyrapint
      complex *16 cfpy(136,136,yrule*yintervals)
      complex *16 cfmy(136,136,yrule*yintervals)
      common /cachedrapint/ cfpy,cfmy
      
c     print *,'in rapintegrals'
      ax=dlog(m**2/sqs**2)

      NMAX1 = mdim
      NMAX2 = mdim
      

!!!   Important, notice that in all the expressions, the dependence on I1, I2 is of the type I2-I1
!!!   This means that the double loops on I1 I2 can be reduced to single loops on I2-I1
!!!   This is true for both analytical and numerical integration
      
c     If there are no cuts on the leptons, calculate the integrals analitically
      if (nolepcuts.eqv..true.) then
         do I1 = 1, NMAX1          !136
            do I2 = 1, NMAX2       !136
               if (I1.eq.I2) then
                  yintp=(CCp/pi)**2*exp(-(Np(I1)+Np(I2))*ax/2d0)
     .                 *(ymax-ymin)
               else
                  yintp=1d0/(-Np(I1)+Np(I2))
     +        *(CCp/pi)**2*exp(-(Np(I1)+Np(I2))*ax/2d0)
     +        *(exp((-Np(I1)+Np(I2))*ymax)-exp((-Np(I1)+Np(I2))*ymin))
               endif
               if (Np(I1).eq.Nm(I2)) then
                  yintm=(CCp/pi)*(CCm/pi)*exp(-(Np(I1)+Nm(I2))*ax/2d0)
     +                 *(ymax-ymin)
               else
                  yintm=1d0/(-Np(I1)+Nm(I2))
     +         *(CCp/pi)*(CCm/pi) *exp(-(Np(I1)+Nm(I2))*ax/2d0)
     +         *(exp((-Np(I1)+Nm(I2))*ymax)-exp((-Np(I1)+Nm(I2))*ymin))
               endif
               Ith0p(I1,I2) = 2d0*yintp*WN(I1)*WN(I2)
               Ith1p(I1,I2) = 0d0
               Ith2p(I1,I2) = 2d0/3d0*yintp*WN(I1)*WN(I2)
               Ith0m(I1,I2) = 2d0*yintm*WN(I1)*WN(I2)
               Ith1m(I1,I2) = 0d0
               Ith2m(I1,I2) = 2d0/3d0*yintm*WN(I1)*WN(I2)
            enddo
         enddo
         call initsigmacthy(m)
         return
      endif
     
c     Initialize integrals      
      do I1 = 1, NMAX1             !136
         do I2 = 1, NMAX2          !136
            Ith0p(I1,I2) = 0d0
            Ith1p(I1,I2) = 0d0
            Ith2p(I1,I2) = 0d0
            Ith0m(I1,I2) = 0d0
            Ith1m(I1,I2) = 0d0
            Ith2m(I1,I2) = 0d0
         enddo
      enddo

c     cache the mass dependent part (ax) of the exponential
      do I1 = 1, NMAX1          !136
         do I2 = 1, NMAX2       !136
      cfpm(I1,I2)=(CCp/pi)**2*exp(-(Np(I1)+Np(I2))*ax/2)
      cfmm(I1,I2)=CCp*CCm/pi**2*exp(-(Np(I1)+Nm(I2))*ax/2)
            enddo
      enddo
      
c     start integration
      xc=0.5d0*(ymin+ymax)
      xm=0.5d0*(ymax-ymin)
      do i=1,yintervals
         ya = ymin+(ymax-ymin)*(i-1)/yintervals
         yb = ymin+(ymax-ymin)*i/yintervals
         xc=0.5d0*(ya+yb)
         xm=0.5d0*(yb-ya)
         do j=1,yrule
            y=xc+xm*xxx(yrule,j)

c     calculate costheta moments as a function of y
         call sety(y)
         call genV4p()
         call cthmoments(cthmom0,cthmom1,cthmom2)
c      print *,cthmom0,cthmom1,cthmom2
         do I1 = 1, NMAX1  !136
            do I2 = 1, NMAX2 !136
c     functions to integrate (rapidity part is cached at initialisation)
               fpy=cfpm(I1,I2)*cfpy(I1,I2,j+(i-1)*yrule)
               fmy=cfmm(I1,I2)*cfmy(I1,I2,j+(i-1)*yrule)
c     integrals
               Ith0p(I1,I2)=Ith0p(I1,I2)+fpy*cthmom0
               Ith1p(I1,I2)=Ith1p(I1,I2)+fpy*cthmom1
               Ith2p(I1,I2)=Ith2p(I1,I2)+fpy*cthmom2

               Ith0m(I1,I2)=Ith0m(I1,I2)+fmy*cthmom0
               Ith1m(I1,I2)=Ith1m(I1,I2)+fmy*cthmom1
               Ith2m(I1,I2)=Ith2m(I1,I2)+fmy*cthmom2
            enddo
         enddo
      enddo
      enddo

c      print *,
c      do I1 = 1, 3              !136
c         do I2 = 1, 3           !136
c            print*,I1,I2,Ith1p(I1,I2),cfpm(I1,I2)
c         enddo
c      enddo
      call initsigmacthy(m)
      return
      end

c     setup the sigmaij integrated in costh and rapidity
      subroutine initsigmacthy(m)
      implicit none
      double precision m

      integer flag5,brflag,fnwa
      common/flags2/flag5,brflag,fnwa

      include 'const.h'
      double precision q2

      double precision facZ,facW,chi1,chi2
      common/vcoup/facZ,facW,chi1,chi2

      double precision Vud,Vus,Vub,Vcd,Vcs,Vcb
      common/cabib/Vud,Vus,Vub,Vcd,Vcs,Vcb

      double precision gevfb
      data gevfb/3.8937966d11/

      double precision ax1,ax2,xx1,xx2
      integer nmax1,nmax2
      common /cxx/ax1,ax2,xx1,xx2,nmax1,nmax2
      
c     Input: integrated costh moments
      integer I1,I2
      complex *16 Ith0p(136,136)
      complex *16 Ith1p(136,136)
      complex *16 Ith2p(136,136)
      complex *16 Ith0m(136,136)
      complex *16 Ith1m(136,136)
      complex *16 Ith2m(136,136)
      common/ITHMOM/Ith0p,Ith1p,Ith2p,Ith0m,Ith1m,Ith2m

c     Output: rapidity and costh integrated born-level matrix elements
      complex *16 sigmaintijp(-5:5,-5:5,136,136)
      complex *16 sigmaintijm(-5:5,-5:5,136,136)
      common/sigmaint/sigmaintijp,sigmaintijm

      q2=m**2
      facZ=1/9d0/pi*q2/((q2-Mz**2)**2+Mz**2*zw**2)*gevfb
      facW=1/9d0/pi*q2/((q2-Mw**2)**2+Mw**2*ww**2)*gevfb
      chi1=(q2-Mz**2)/q2
      chi2=((q2-Mz**2)**2+Mz**2*zw**2)/q2/q2

      do I1 = 1, NMAX1    !136
         do I2 = 1, NMAX2  !136
            if (flag5.eq.3) then
               sigmaintijp(1,-1,I1,I2)=facZ*( 
     \              (gLZu**2+gRZu**2)*(fLZ**2+fRZ**2)
     \              *   (Ith0p(I1,I2)+Ith2p(I1,I2))
     \              -   (gLZu**2-gRZu**2)*(fLZ**2-fRZ**2)
     \              *     (2d0*Ith1p(I1,I2)))
               sigmaintijp(4,-4,I1,I2)=sigmaintijp(1,-1,I1,I2)

               sigmaintijp(-1,1,I1,I2)=facZ*( 
     \              (gLZu**2+gRZu**2)*(fLZ**2+fRZ**2)
     \              *     (Ith0p(I1,I2)+Ith2p(I1,I2))
     \              +   (gLZu**2-gRZu**2)*(fLZ**2-fRZ**2)
     \              *     (2d0*Ith1p(I1,I2)))
               sigmaintijp(-4,4,I1,I2)=sigmaintijp(-1,1,I1,I2)
c     
               sigmaintijp(2,-2,I1,I2)=facZ*( 
     \              (gLZd**2+gRZd**2)*(fLZ**2+fRZ**2)
     \              *     (Ith0p(I1,I2)+Ith2p(I1,I2))
     \              -   (gLZd**2-gRZd**2)*(fLZ**2-fRZ**2)
     \              *     (2d0*Ith1p(I1,I2)))
               sigmaintijp(3,-3,I1,I2)=sigmaintijp(2,-2,I1,I2)
               sigmaintijp(5,-5,I1,I2)=sigmaintijp(2,-2,I1,I2)
c     
               sigmaintijp(-2,2,I1,I2)=facZ*( 
     \              (gLZd**2+gRZd**2)*(fLZ**2+fRZ**2)
     \              *     (Ith0p(I1,I2)+Ith2p(I1,I2))
     \              +   (gLZd**2-gRZd**2)*(fLZ**2-fRZ**2)
     \              *     (2d0*Ith1p(I1,I2)))
               sigmaintijp(-3,3,I1,I2)=sigmaintijp(-2,2,I1,I2)
               sigmaintijp(-5,5,I1,I2)=sigmaintijp(-2,2,I1,I2)
c     negative
               sigmaintijm(1,-1,I1,I2)=facZ*( 
     \              (gLZu**2+gRZu**2)*(fLZ**2+fRZ**2)
     \              *   (Ith0m(I1,I2)+Ith2m(I1,I2))
     \              -   (gLZu**2-gRZu**2)*(fLZ**2-fRZ**2)
     \              *     (2d0*Ith1m(I1,I2)))
               sigmaintijm(4,-4,I1,I2)=sigmaintijm(1,-1,I1,I2)

               sigmaintijm(-1,1,I1,I2)=facZ*( 
     \              (gLZu**2+gRZu**2)*(fLZ**2+fRZ**2)
     \              *     (Ith0m(I1,I2)+Ith2m(I1,I2))
     \              +   (gLZu**2-gRZu**2)*(fLZ**2-fRZ**2)
     \              *     (2d0*Ith1m(I1,I2)))
               sigmaintijm(-4,4,I1,I2)=sigmaintijm(-1,1,I1,I2)
c     
               sigmaintijm(2,-2,I1,I2)=facZ*( 
     \              (gLZd**2+gRZd**2)*(fLZ**2+fRZ**2)
     \              *     (Ith0m(I1,I2)+Ith2m(I1,I2))
     \              -   (gLZd**2-gRZd**2)*(fLZ**2-fRZ**2)
     \              *     (2d0*Ith1m(I1,I2)))
               sigmaintijm(3,-3,I1,I2)=sigmaintijm(2,-2,I1,I2)
               sigmaintijm(5,-5,I1,I2)=sigmaintijm(2,-2,I1,I2)
c     
               sigmaintijm(-2,2,I1,I2)=facZ*( 
     \              (gLZd**2+gRZd**2)*(fLZ**2+fRZ**2)
     \              *     (Ith0m(I1,I2)+Ith2m(I1,I2))
     \              +   (gLZd**2-gRZd**2)*(fLZ**2-fRZ**2)
     \              *     (2d0*Ith1m(I1,I2)))
               sigmaintijm(-3,3,I1,I2)=sigmaintijm(-2,2,I1,I2)
               sigmaintijm(-5,5,I1,I2)=sigmaintijm(-2,2,I1,I2)
            elseif(flag5.eq.21) then
               sigmaintijp(1,-2,I1,I2)=facW/16d0*vud**2*gLW**2*fLW**2
     .              *(Ith0p(I1,I2)-2*Ith1p(I1,I2)+Ith2p(I1,I2))
               sigmaintijp(-2,1,I1,I2)=facW/16d0*vud**2*gLW**2*fLW**2
     .              *(Ith0p(I1,I2)+2*Ith1p(I1,I2)+Ith2p(I1,I2))
               sigmaintijp(1,-3,I1,I2)=facW/16d0*vus**2*gLW**2*fLW**2
     .              *(Ith0p(I1,I2)-2*Ith1p(I1,I2)+Ith2p(I1,I2))
               sigmaintijp(-3,1,I1,I2)=facW/16d0*vus**2*gLW**2*fLW**2
     .              *(Ith0p(I1,I2)+2*Ith1p(I1,I2)+Ith2p(I1,I2))
               sigmaintijp(1,-5,I1,I2)=facW/16d0*vub**2*gLW**2*fLW**2
     .              *(Ith0p(I1,I2)-2*Ith1p(I1,I2)+Ith2p(I1,I2))
               sigmaintijp(-5,1,I1,I2)=facW/16d0*vub**2*gLW**2*fLW**2
     .              *(Ith0p(I1,I2)+2*Ith1p(I1,I2)+Ith2p(I1,I2))
               sigmaintijp(4,-3,I1,I2)=facW/16d0*vcs**2*gLW**2*fLW**2
     .              *(Ith0p(I1,I2)-2*Ith1p(I1,I2)+Ith2p(I1,I2))
               sigmaintijp(-3,4,I1,I2)=facW/16d0*vcs**2*gLW**2*fLW**2
     .              *(Ith0p(I1,I2)+2*Ith1p(I1,I2)+Ith2p(I1,I2))
               sigmaintijp(4,-2,I1,I2)=facW/16d0*vcd**2*gLW**2*fLW**2
     .              *(Ith0p(I1,I2)-2*Ith1p(I1,I2)+Ith2p(I1,I2))
               sigmaintijp(-2,4,I1,I2)=facW/16d0*vcd**2*gLW**2*fLW**2
     .              *(Ith0p(I1,I2)+2*Ith1p(I1,I2)+Ith2p(I1,I2))
               sigmaintijp(4,-5,I1,I2)=facW/16d0*vcb**2*gLW**2*fLW**2
     .              *(Ith0p(I1,I2)-2*Ith1p(I1,I2)+Ith2p(I1,I2))
               sigmaintijp(-5,4,I1,I2)=facW/16d0*vcb**2*gLW**2*fLW**2
     .              *(Ith0p(I1,I2)+2*Ith1p(I1,I2)+Ith2p(I1,I2))
            elseif(flag5.eq.22) then ! CHECK IT    !!!!!
               sigmaintijp(2,-1,I1,I2)=facW/16d0*vud**2*gLW**2*fLW**2
     .              *(Ith0p(I1,I2)-2*Ith1p(I1,I2)+Ith2p(I1,I2))
               sigmaintijp(-1,2,I1,I2)=facW/16d0*vud**2*gLW**2*fLW**2
     .              *(Ith0p(I1,I2)+2*Ith1p(I1,I2)+Ith2p(I1,I2))
               sigmaintijp(3,-1,I1,I2)=facW/16d0*vus**2*gLW**2*fLW**2
     .              *(Ith0p(I1,I2)-2*Ith1p(I1,I2)+Ith2p(I1,I2))
               sigmaintijp(-1,3,I1,I2)=facW/16d0*vus**2*gLW**2*fLW**2
     .              *(Ith0p(I1,I2)+2*Ith1p(I1,I2)+Ith2p(I1,I2))
               sigmaintijp(5,-1,I1,I2)=facW/16d0*vub**2*gLW**2*fLW**2
     .              *(Ith0p(I1,I2)-2*Ith1p(I1,I2)+Ith2p(I1,I2))
               sigmaintijp(-1,5,I1,I2)=facW/16d0*vub**2*gLW**2*fLW**2
     .              *(Ith0p(I1,I2)+2*Ith1p(I1,I2)+Ith2p(I1,I2))
               sigmaintijp(3,-4,I1,I2)=facW/16d0*vcs**2*gLW**2*fLW**2
     .              *(Ith0p(I1,I2)-2*Ith1p(I1,I2)+Ith2p(I1,I2))
               sigmaintijp(-4,3,I1,I2)=facW/16d0*vcs**2*gLW**2*fLW**2
     .              *(Ith0p(I1,I2)+2*Ith1p(I1,I2)+Ith2p(I1,I2))
               sigmaintijp(2,-4,I1,I2)=facW/16d0*vcd**2*gLW**2*fLW**2
     .              *(Ith0p(I1,I2)-2*Ith1p(I1,I2)+Ith2p(I1,I2))
               sigmaintijp(-4,2,I1,I2)=facW/16d0*vcd**2*gLW**2*fLW**2
     .              *(Ith0p(I1,I2)+2*Ith1p(I1,I2)+Ith2p(I1,I2))
               sigmaintijp(5,-4,I1,I2)=facW/16d0*vcb**2*gLW**2*fLW**2
     .              *(Ith0p(I1,I2)-2*Ith1p(I1,I2)+Ith2p(I1,I2))
               sigmaintijp(-4,5,I1,I2)=facW/16d0*vcb**2*gLW**2*fLW**2
     .              *(Ith0p(I1,I2)+2*Ith1p(I1,I2)+Ith2p(I1,I2))
            elseif (flag5.eq.5) then !  T=TYPO in hep-ph/9704239 !!!!!
               sigmaintijp(1,-1,I1,I2)=facZ*( 
     \              ( (gLZu**2+gRZu**2)*(fLZ**2+fRZ**2) 
     \              +1/4d0*(4d0*pi*aem)**2*(eequ)**2*chi2 
     \              -1/2d0*(4d0*pi*aem)*eequ*(gLZu+gRZu)*(fLZ+fRZ)*chi1)
!     
     \              *(Ith0p(I1,I2)+Ith2p(I1,I2))
     \              -  ( (gLZu**2-gRZu**2)*(fLZ**2-fRZ**2)
     \              -1/2d0*(4d0*pi*aem)*eequ*(gLZu-gRZu)*(fLZ-fRZ)*chi1)
!     
     \              *(2d0*Ith1p(I1,I2)) 
     \              )
               sigmaintijp(4,-4,I1,I2)=sigmaintijp(1,-1,I1,I2)
               sigmaintijp(-1,1,I1,I2)=facZ*( 
     \              ( (gLZu**2+gRZu**2)*(fLZ**2+fRZ**2) 
     \              +1/4d0*(4d0*pi*aem)**2*(eequ)**2*chi2 
     \              -1/2d0*(4d0*pi*aem)*eequ*(gLZu+gRZu)*(fLZ+fRZ)*chi1)
!     
     \              *(Ith0p(I1,I2)+Ith2p(I1,I2))
     \              +  ( (gLZu**2-gRZu**2)*(fLZ**2-fRZ**2)
     \              -1/2d0*(4d0*pi*aem)*eequ*(gLZu-gRZu)*(fLZ-fRZ)*chi1)
!     
     \              *(2d0*Ith1p(I1,I2)) 
     \              )
               sigmaintijp(-4,4,I1,I2)=sigmaintijp(-1,1,I1,I2)
c     
               sigmaintijp(2,-2,I1,I2)=facZ*( 
     \              ( (gLZd**2+gRZd**2)*(fLZ**2+fRZ**2) 
     \              +1/4d0*(4d0*pi*aem)**2*(eeqd)**2*chi2 
     \              -1/2d0*(4d0*pi*aem)*eeqd*(gLZd+gRZd)*(fLZ+fRZ)*chi1)
!     
     \              *(Ith0p(I1,I2)+Ith2p(I1,I2))
     \              -  ( (gLZd**2-gRZd**2)*(fLZ**2-fRZ**2)
     \              -1/2d0*(4d0*pi*aem)*eeqd*(gLZd-gRZd)*(fLZ-fRZ)*chi1)
!     
     \              *(2d0*Ith1p(I1,I2)) 
     \              )
               sigmaintijp(3,-3,I1,I2)=sigmaintijp(2,-2,I1,I2)
               sigmaintijp(5,-5,I1,I2)=sigmaintijp(2,-2,I1,I2)
               sigmaintijp(-2,2,I1,I2)=facZ*( 
     \              ( (gLZd**2+gRZd**2)*(fLZ**2+fRZ**2) 
     \              +1/4d0*(4d0*pi*aem)**2*(eeqd)**2*chi2 
     \              -1/2d0*(4d0*pi*aem)*eeqd*(gLZd+gRZd)*(fLZ+fRZ)*chi1)
!     
     \              *(Ith0p(I1,I2)+Ith2p(I1,I2))
     \              +  ( (gLZd**2-gRZd**2)*(fLZ**2-fRZ**2)
     \              -1/2d0*(4d0*pi*aem)*eeqd*(gLZd-gRZd)*(fLZ-fRZ)*chi1)
!     
     \              *(2d0*Ith1p(I1,I2)) 
     \              )
               sigmaintijp(-3,3,I1,I2)=sigmaintijp(-2,2,I1,I2)
               sigmaintijp(-5,5,I1,I2)=sigmaintijp(-2,2,I1,I2)
c     negative
               sigmaintijm(1,-1,I1,I2)=facZ*( 
     \              ( (gLZu**2+gRZu**2)*(fLZ**2+fRZ**2) 
     \              +1/4d0*(4d0*pi*aem)**2*(eequ)**2*chi2 
     \              -1/2d0*(4d0*pi*aem)*eequ*(gLZu+gRZu)*(fLZ+fRZ)*chi1)
!     
     \              *(Ith0m(I1,I2)+Ith2m(I1,I2))
     \              -  ( (gLZu**2-gRZu**2)*(fLZ**2-fRZ**2)
     \              -1/2d0*(4d0*pi*aem)*eequ*(gLZu-gRZu)*(fLZ-fRZ)*chi1)
!     
     \              *(2d0*Ith1m(I1,I2)) 
     \              )
               sigmaintijm(4,-4,I1,I2)=sigmaintijm(1,-1,I1,I2)
               sigmaintijm(-1,1,I1,I2)=facZ*( 
     \              ( (gLZu**2+gRZu**2)*(fLZ**2+fRZ**2) 
     \              +1/4d0*(4d0*pi*aem)**2*(eequ)**2*chi2 
     \              -1/2d0*(4d0*pi*aem)*eequ*(gLZu+gRZu)*(fLZ+fRZ)*chi1)
!     
     \              *(Ith0m(I1,I2)+Ith2m(I1,I2))
     \              +  ( (gLZu**2-gRZu**2)*(fLZ**2-fRZ**2)
     \              -1/2d0*(4d0*pi*aem)*eequ*(gLZu-gRZu)*(fLZ-fRZ)*chi1)
!     
     \              *(2d0*Ith1m(I1,I2)) 
     \              )
               sigmaintijm(-4,4,I1,I2)=sigmaintijm(-1,1,I1,I2)
c     
               sigmaintijm(2,-2,I1,I2)=facZ*( 
     \              ( (gLZd**2+gRZd**2)*(fLZ**2+fRZ**2) 
     \              +1/4d0*(4d0*pi*aem)**2*(eeqd)**2*chi2 
     \              -1/2d0*(4d0*pi*aem)*eeqd*(gLZd+gRZd)*(fLZ+fRZ)*chi1)
!     
     \              *(Ith0m(I1,I2)+Ith2m(I1,I2))
     \              -  ( (gLZd**2-gRZd**2)*(fLZ**2-fRZ**2)
     \              -1/2d0*(4d0*pi*aem)*eeqd*(gLZd-gRZd)*(fLZ-fRZ)*chi1)
!     
     \              *(2d0*Ith1m(I1,I2)) 
     \              )
               sigmaintijm(3,-3,I1,I2)=sigmaintijm(2,-2,I1,I2)
               sigmaintijm(5,-5,I1,I2)=sigmaintijm(2,-2,I1,I2)
               sigmaintijm(-2,2,I1,I2)=facZ*( 
     \              ( (gLZd**2+gRZd**2)*(fLZ**2+fRZ**2) 
     \              +1/4d0*(4d0*pi*aem)**2*(eeqd)**2*chi2 
     \              -1/2d0*(4d0*pi*aem)*eeqd*(gLZd+gRZd)*(fLZ+fRZ)*chi1)
!     
     \              *(Ith0m(I1,I2)+Ith2m(I1,I2))
     \              +  ( (gLZd**2-gRZd**2)*(fLZ**2-fRZ**2)
     \              -1/2d0*(4d0*pi*aem)*eeqd*(gLZd-gRZd)*(fLZ-fRZ)*chi1)
!     
     \              *(2d0*Ith1m(I1,I2)) 
     \              )
               sigmaintijm(-3,3,I1,I2)=sigmaintijm(-2,2,I1,I2)
               sigmaintijm(-5,5,I1,I2)=sigmaintijm(-2,2,I1,I2)
            endif      
         enddo
      enddo
      return
      end

c     setup the sigmaij integrated in costh
      subroutine initsigmacth(m,cthmom0,cthmom1,cthmom2)
      implicit none
      double precision m
      double precision cthmom0,cthmom1,cthmom2
      
      integer flag5,brflag,fnwa
      common/flags2/flag5,brflag,fnwa

      include 'const.h'
      double precision q2

      double precision facZ,facW,chi1,chi2
      common/vcoup/facZ,facW,chi1,chi2

      double precision Vud,Vus,Vub,Vcd,Vcs,Vcb
      common/cabib/Vud,Vus,Vub,Vcd,Vcs,Vcb

      double precision gevfb
      data gevfb/3.8937966d11/

      include 'sigmaij_inc.f'

      q2=m**2
      facZ=1/9d0/pi*q2/((q2-Mz**2)**2+Mz**2*zw**2)*gevfb
      facW=1/9d0/pi*q2/((q2-Mw**2)**2+Mw**2*ww**2)*gevfb
      chi1=(q2-Mz**2)/q2
      chi2=((q2-Mz**2)**2+Mz**2*zw**2)/q2/q2

      if (flag5.eq.3) then
         sigmaij(1,-1)=facZ*( 
     \        (gLZu**2+gRZu**2)*(fLZ**2+fRZ**2)*(cthmom0+cthmom2)
     \        -   (gLZu**2-gRZu**2)*(fLZ**2-fRZ**2)*(2d0*cthmom1))
         sigmaij(4,-4)=sigmaij(1,-1)
c     
         sigmaij(-1,1)=facZ*( 
     \        (gLZu**2+gRZu**2)*(fLZ**2+fRZ**2)*(cthmom0+cthmom2)
     \        +   (gLZu**2-gRZu**2)*(fLZ**2-fRZ**2)*(2d0*cthmom1))
         sigmaij(-4,4)=sigmaij(-1,1)
c     
         sigmaij(2,-2)=facZ*( 
     \        (gLZd**2+gRZd**2)*(fLZ**2+fRZ**2)*(cthmom0+cthmom2)
     \        -   (gLZd**2-gRZd**2)*(fLZ**2-fRZ**2)*(2d0*cthmom1))
         sigmaij(3,-3)=sigmaij(2,-2)
         sigmaij(5,-5)=sigmaij(2,-2)
c     
         sigmaij(-2,2)=facZ*( 
     \        (gLZd**2+gRZd**2)*(fLZ**2+fRZ**2)*(cthmom0+cthmom2)
     \        +   (gLZd**2-gRZd**2)*(fLZ**2-fRZ**2)*(2d0*cthmom1))
         sigmaij(-3,3)=sigmaij(-2,2)
         sigmaij(-5,5)=sigmaij(-2,2)
      elseif(flag5.eq.21) then
         sigmaij(1,-2)=facW/16d0*vud**2*gLW**2*fLW**2
     .        *(cthmom0-2*cthmom1+cthmom2)
         sigmaij(-2,1)=facW/16d0*vud**2*gLW**2*fLW**2
     .        *(cthmom0+2*cthmom1+cthmom2)
         sigmaij(1,-3)=facW/16d0*vus**2*gLW**2*fLW**2
     .        *(cthmom0-2*cthmom1+cthmom2)
         sigmaij(-3,1)=facW/16d0*vus**2*gLW**2*fLW**2
     .        *(cthmom0+2*cthmom1+cthmom2)
         sigmaij(1,-5)=facW/16d0*vub**2*gLW**2*fLW**2
     .        *(cthmom0-2*cthmom1+cthmom2)
         sigmaij(-5,1)=facW/16d0*vub**2*gLW**2*fLW**2
     .        *(cthmom0+2*cthmom1+cthmom2)
         sigmaij(4,-3)=facW/16d0*vcs**2*gLW**2*fLW**2
     .        *(cthmom0-2*cthmom1+cthmom2)
         sigmaij(-3,4)=facW/16d0*vcs**2*gLW**2*fLW**2
     .        *(cthmom0+2*cthmom1+cthmom2)
         sigmaij(4,-2)=facW/16d0*vcd**2*gLW**2*fLW**2
     .        *(cthmom0-2*cthmom1+cthmom2)
         sigmaij(-2,4)=facW/16d0*vcd**2*gLW**2*fLW**2
     .        *(cthmom0+2*cthmom1+cthmom2)
         sigmaij(4,-5)=facW/16d0*vcb**2*gLW**2*fLW**2
     .        *(cthmom0-2*cthmom1+cthmom2)
         sigmaij(-5,4)=facW/16d0*vcb**2*gLW**2*fLW**2
     .        *(cthmom0+2*cthmom1+cthmom2)
      elseif(flag5.eq.22) then  ! CHECK IT    !!!!!
         sigmaij(2,-1)=facW/16d0*vud**2*gLW**2*fLW**2
     .        *(cthmom0-2*cthmom1+cthmom2)
         sigmaij(-1,2)=facW/16d0*vud**2*gLW**2*fLW**2
     .        *(cthmom0+2*cthmom1+cthmom2)
         sigmaij(3,-1)=facW/16d0*vus**2*gLW**2*fLW**2
     .        *(cthmom0-2*cthmom1+cthmom2)
         sigmaij(-1,3)=facW/16d0*vus**2*gLW**2*fLW**2
     .        *(cthmom0+2*cthmom1+cthmom2)
         sigmaij(5,-1)=facW/16d0*vub**2*gLW**2*fLW**2
     .        *(cthmom0-2*cthmom1+cthmom2)
         sigmaij(-1,5)=facW/16d0*vub**2*gLW**2*fLW**2
     .        *(cthmom0+2*cthmom1+cthmom2)
         sigmaij(3,-4)=facW/16d0*vcs**2*gLW**2*fLW**2
     .        *(cthmom0-2*cthmom1+cthmom2)
         sigmaij(-4,3)=facW/16d0*vcs**2*gLW**2*fLW**2
     .        *(cthmom0+2*cthmom1+cthmom2)
         sigmaij(2,-4)=facW/16d0*vcd**2*gLW**2*fLW**2
     .        *(cthmom0-2*cthmom1+cthmom2)
         sigmaij(-4,2)=facW/16d0*vcd**2*gLW**2*fLW**2
     .        *(cthmom0+2*cthmom1+cthmom2)
        sigmaij(5,-4)=facW/16d0*vcb**2*gLW**2*fLW**2
     .        *(cthmom0-2*cthmom1+cthmom2)
        sigmaij(-4,5)=facW/16d0*vcb**2*gLW**2*fLW**2
     .        *(cthmom0+2*cthmom1+cthmom2)
      elseif (flag5.eq.5) then
         sigmaij(1,-1)=facZ*( 
     \        ( (gLZu**2+gRZu**2)*(fLZ**2+fRZ**2) 
     \        +1/4d0*(4d0*pi*aem)**2*(eequ)**2*chi2 
     \        -1/2d0*(4d0*pi*aem)*eequ*(gLZu+gRZu)*(fLZ+fRZ)*chi1 )
!     
     \        *(cthmom0+cthmom2)
     \        -  ( (gLZu**2-gRZu**2)*(fLZ**2-fRZ**2)
     \        -1/2d0*(4d0*pi*aem)*eequ*(gLZu-gRZu)*(fLZ-fRZ)*chi1 )
!     
     \        *(2d0*cthmom1) 
     \        )
         sigmaij(4,-4)=sigmaij(1,-1)
         sigmaij(-1,1)=facZ*( 
     \        ( (gLZu**2+gRZu**2)*(fLZ**2+fRZ**2) 
     \        +1/4d0*(4d0*pi*aem)**2*(eequ)**2*chi2 
     \        -1/2d0*(4d0*pi*aem)*eequ*(gLZu+gRZu)*(fLZ+fRZ)*chi1 )
!     
     \        *(cthmom0+cthmom2)
     \        +  ( (gLZu**2-gRZu**2)*(fLZ**2-fRZ**2)
     \        -1/2d0*(4d0*pi*aem)*eequ*(gLZu-gRZu)*(fLZ-fRZ)*chi1 )
!     
     \        *(2d0*cthmom1) 
     \        )
         sigmaij(-4,4)=sigmaij(-1,1)
c     
         sigmaij(2,-2)=facZ*( 
     \        ( (gLZd**2+gRZd**2)*(fLZ**2+fRZ**2) 
     \        +1/4d0*(4d0*pi*aem)**2*(eeqd)**2*chi2 
     \        -1/2d0*(4d0*pi*aem)*eeqd*(gLZd+gRZd)*(fLZ+fRZ)*chi1 )
!     
     \        *(cthmom0+cthmom2)
     \        -  ( (gLZd**2-gRZd**2)*(fLZ**2-fRZ**2)
     \        -1/2d0*(4d0*pi*aem)*eeqd*(gLZd-gRZd)*(fLZ-fRZ)*chi1 )
!     
     \        *(2d0*cthmom1) 
     \        )
         sigmaij(3,-3)=sigmaij(2,-2)
         sigmaij(5,-5)=sigmaij(2,-2)
         sigmaij(-2,2)=facZ*( 
     \        ( (gLZd**2+gRZd**2)*(fLZ**2+fRZ**2) 
     \        +1/4d0*(4d0*pi*aem)**2*(eeqd)**2*chi2 
     \        -1/2d0*(4d0*pi*aem)*eeqd*(gLZd+gRZd)*(fLZ+fRZ)*chi1 )
!     
     \        *(cthmom0+cthmom2)
     \        +  ( (gLZd**2-gRZd**2)*(fLZ**2-fRZ**2)
     \        -1/2d0*(4d0*pi*aem)*eeqd*(gLZd-gRZd)*(fLZ-fRZ)*chi1 )
!     
     \        *(2d0*cthmom1) 
     \        )
         sigmaij(-3,3)=sigmaij(-2,2)
         sigmaij(-5,5)=sigmaij(-2,2)
      endif
      return
      end


c     setup the sigmaij integrated in costh
      subroutine initsigma(m,costh)
      implicit none
      double precision m,costh
      
      integer flag5,brflag,fnwa
      common/flags2/flag5,brflag,fnwa

      include 'const.h'
      double precision q2

      double precision facZ,facW,chi1,chi2
      common/vcoup/facZ,facW,chi1,chi2

      double precision Vud,Vus,Vub,Vcd,Vcs,Vcb
      common/cabib/Vud,Vus,Vub,Vcd,Vcs,Vcb

      double precision gevfb
      data gevfb/3.8937966d11/

      include 'sigmaij_inc.f'

      q2=m**2
      facZ=1/9d0/pi*q2/((q2-Mz**2)**2+Mz**2*zw**2)*gevfb
      facW=1/9d0/pi*q2/((q2-Mw**2)**2+Mw**2*ww**2)*gevfb
      chi1=(q2-Mz**2)/q2
      chi2=((q2-Mz**2)**2+Mz**2*zw**2)/q2/q2
      if (flag5.eq.3) then
         sigmaij(1,-1)=facZ*( 
     \        (gLZu**2+gRZu**2)*(fLZ**2+fRZ**2)*(1d0+costh**2)
     \        -   (gLZu**2-gRZu**2)*(fLZ**2-fRZ**2)*(2d0*costh)   )
         sigmaij(4,-4)=sigmaij(1,-1)
c     
         sigmaij(-1,1)=facZ*( 
     \        (gLZu**2+gRZu**2)*(fLZ**2+fRZ**2)*(1d0+costh**2)
     \        +   (gLZu**2-gRZu**2)*(fLZ**2-fRZ**2)*(2d0*costh)   )
         sigmaij(-4,4)=sigmaij(-1,1)
c     
         sigmaij(2,-2)=facZ*( 
     \        (gLZd**2+gRZd**2)*(fLZ**2+fRZ**2)*(1d0+costh**2)
     \        -   (gLZd**2-gRZd**2)*(fLZ**2-fRZ**2)*(2d0*costh)   )
         sigmaij(3,-3)=sigmaij(2,-2)
         sigmaij(5,-5)=sigmaij(2,-2)
c
         sigmaij(-2,2)=facZ*( 
     \        (gLZd**2+gRZd**2)*(fLZ**2+fRZ**2)*(1d0+costh**2)
     \        +   (gLZd**2-gRZd**2)*(fLZ**2-fRZ**2)*(2d0*costh)   )
         sigmaij(-3,3)=sigmaij(-2,2)
         sigmaij(-5,5)=sigmaij(-2,2)

      elseif(flag5.eq.21) then
         sigmaij(1,-2)=facW/16d0*vud**2*gLW**2*fLW**2*(1d0-costh)**2
         sigmaij(-2,1)=facW/16d0*vud**2*gLW**2*fLW**2*(1d0+costh)**2
         sigmaij(1,-3)=facW/16d0*vus**2*gLW**2*fLW**2*(1d0-costh)**2
         sigmaij(-3,1)=facW/16d0*vus**2*gLW**2*fLW**2*(1d0+costh)**2
         sigmaij(1,-5)=facW/16d0*vub**2*gLW**2*fLW**2*(1d0-costh)**2
         sigmaij(-5,1)=facW/16d0*vub**2*gLW**2*fLW**2*(1d0+costh)**2
         sigmaij(4,-3)=facW/16d0*vcs**2*gLW**2*fLW**2*(1d0-costh)**2
         sigmaij(-3,4)=facW/16d0*vcs**2*gLW**2*fLW**2*(1d0+costh)**2
         sigmaij(4,-2)=facW/16d0*vcd**2*gLW**2*fLW**2*(1d0-costh)**2
         sigmaij(-2,4)=facW/16d0*vcd**2*gLW**2*fLW**2*(1d0+costh)**2
         sigmaij(4,-5)=facW/16d0*vcb**2*gLW**2*fLW**2*(1d0-costh)**2
         sigmaij(-5,4)=facW/16d0*vcb**2*gLW**2*fLW**2*(1d0+costh)**2
      elseif(flag5.eq.22) then  ! CHECK IT    !!!!!
         sigmaij(2,-1)=facW/16d0*vud**2*gLW**2*fLW**2*(1d0-costh)**2
         sigmaij(-1,2)=facW/16d0*vud**2*gLW**2*fLW**2*(1d0+costh)**2
         sigmaij(3,-1)=facW/16d0*vus**2*gLW**2*fLW**2*(1d0-costh)**2
         sigmaij(-1,3)=facW/16d0*vus**2*gLW**2*fLW**2*(1d0+costh)**2
         sigmaij(5,-1)=facW/16d0*vub**2*gLW**2*fLW**2*(1d0-costh)**2
         sigmaij(-1,5)=facW/16d0*vub**2*gLW**2*fLW**2*(1d0+costh)**2
         sigmaij(3,-4)=facW/16d0*vcs**2*gLW**2*fLW**2*(1d0-costh)**2
         sigmaij(-4,3)=facW/16d0*vcs**2*gLW**2*fLW**2*(1d0+costh)**2
         sigmaij(2,-4)=facW/16d0*vcd**2*gLW**2*fLW**2*(1d0-costh)**2
         sigmaij(-4,2)=facW/16d0*vcd**2*gLW**2*fLW**2*(1d0+costh)**2
         sigmaij(5,-4)=facW/16d0*vcb**2*gLW**2*fLW**2*(1d0-costh)**2
         sigmaij(-4,5)=facW/16d0*vcb**2*gLW**2*fLW**2*(1d0+costh)**2
      elseif(flag5.eq.2) then   ! CHECK IT     !!!!!
         sigmaij(1,-2)=facW/16d0*vud**2*gLW**2*fLW**2*(1d0-costh)**2
         sigmaij(-2,1)=facW/16d0*vud**2*gLW**2*fLW**2*(1d0+costh)**2
         sigmaij(1,-3)=facW/16d0*vus**2*gLW**2*fLW**2*(1d0-costh)**2
         sigmaij(-3,1)=facW/16d0*vus**2*gLW**2*fLW**2*(1d0+costh)**2
         sigmaij(1,-5)=facW/16d0*vub**2*gLW**2*fLW**2*(1d0-costh)**2
         sigmaij(-5,1)=facW/16d0*vub**2*gLW**2*fLW**2*(1d0+costh)**2
         sigmaij(4,-3)=facW/16d0*vcs**2*gLW**2*fLW**2*(1d0-costh)**2
         sigmaij(-3,4)=facW/16d0*vcs**2*gLW**2*fLW**2*(1d0+costh)**2
         sigmaij(4,-2)=facW/16d0*vcd**2*gLW**2*fLW**2*(1d0-costh)**2
         sigmaij(-2,4)=facW/16d0*vcd**2*gLW**2*fLW**2*(1d0+costh)**2
         sigmaij(4,-5)=facW/16d0*vcb**2*gLW**2*fLW**2*(1d0-costh)**2
         sigmaij(-5,4)=facW/16d0*vcb**2*gLW**2*fLW**2*(1d0+costh)**2
         sigmaij(2,-1)=facW/16d0*vud**2*gLW**2*fLW**2*(1d0-costh)**2
         sigmaij(-1,2)=facW/16d0*vud**2*gLW**2*fLW**2*(1d0+costh)**2
         sigmaij(3,-1)=facW/16d0*vus**2*gLW**2*fLW**2*(1d0-costh)**2
         sigmaij(-1,3)=facW/16d0*vus**2*gLW**2*fLW**2*(1d0+costh)**2
         sigmaij(5,-1)=facW/16d0*vub**2*gLW**2*fLW**2*(1d0-costh)**2
         sigmaij(-1,5)=facW/16d0*vub**2*gLW**2*fLW**2*(1d0+costh)**2
         sigmaij(3,-4)=facW/16d0*vcs**2*gLW**2*fLW**2*(1d0-costh)**2
         sigmaij(-4,3)=facW/16d0*vcs**2*gLW**2*fLW**2*(1d0+costh)**2
         sigmaij(2,-4)=facW/16d0*vcd**2*gLW**2*fLW**2*(1d0-costh)**2
         sigmaij(-4,2)=facW/16d0*vcd**2*gLW**2*fLW**2*(1d0+costh)**2
         sigmaij(5,-4)=facW/16d0*vcb**2*gLW**2*fLW**2*(1d0-costh)**2
         sigmaij(-4,5)=facW/16d0*vcb**2*gLW**2*fLW**2*(1d0+costh)**2
      elseif (flag5.eq.5) then  !  T=TYPO in hep-ph/9704239 !!!!!
            sigmaij(1,-1)=facZ*( 
     \         ( (gLZu**2+gRZu**2)*(fLZ**2+fRZ**2) 
!T     \          +4d0*(4d0*pi*aem)**2*(eequ)**2*chi2 
!T     \          -2d0*(4d0*pi*aem)*eequ*(gLZu+gRZu)*(fLZ+fRZ)*chi1 )
     \          +1/4d0*(4d0*pi*aem)**2*(eequ)**2*chi2 
     \          -1/2d0*(4d0*pi*aem)*eequ*(gLZu+gRZu)*(fLZ+fRZ)*chi1 )
!
     \          *(1d0+costh**2)
     \      -  ( (gLZu**2-gRZu**2)*(fLZ**2-fRZ**2)
!T     \          -2d0*(4d0*pi*aem)*eequ*(gLZu-gRZu)*(fLZ-fRZ)*chi1 )
     \          -1/2d0*(4d0*pi*aem)*eequ*(gLZu-gRZu)*(fLZ-fRZ)*chi1 )
!
     \          *(2d0*costh) 
     \           )
            sigmaij(4,-4)=sigmaij(1,-1)
            sigmaij(-1,1)=facZ*( 
     \         ( (gLZu**2+gRZu**2)*(fLZ**2+fRZ**2) 
!     \          +4d0*(4d0*pi*aem)**2*(eequ)**2*chi2 
!     \          -2d0*(4d0*pi*aem)*eequ*(gLZu+gRZu)*(fLZ+fRZ)*chi1 )
     \          +1/4d0*(4d0*pi*aem)**2*(eequ)**2*chi2 
     \          -1/2d0*(4d0*pi*aem)*eequ*(gLZu+gRZu)*(fLZ+fRZ)*chi1 )
!
     \          *(1d0+costh**2)
     \      +  ( (gLZu**2-gRZu**2)*(fLZ**2-fRZ**2)
!T     \          -2d0*(4d0*pi*aem)*eequ*(gLZu-gRZu)*(fLZ-fRZ)*chi1 )
     \          -1/2d0*(4d0*pi*aem)*eequ*(gLZu-gRZu)*(fLZ-fRZ)*chi1 )
!
     \          *(2d0*costh) 
     \           )
            sigmaij(-4,4)=sigmaij(-1,1)
c
            sigmaij(2,-2)=facZ*( 
     \         ( (gLZd**2+gRZd**2)*(fLZ**2+fRZ**2) 
!     \          +4d0*(4d0*pi*aem)**2*(eeqd)**2*chi2 
!     \          -2d0*(4d0*pi*aem)*eeqd*(gLZd+gRZd)*(fLZ+fRZ)*chi1 )
     \          +1/4d0*(4d0*pi*aem)**2*(eeqd)**2*chi2 
     \          -1/2d0*(4d0*pi*aem)*eeqd*(gLZd+gRZd)*(fLZ+fRZ)*chi1 )
!
     \          *(1d0+costh**2)
     \      -  ( (gLZd**2-gRZd**2)*(fLZ**2-fRZ**2)
!     \          -2d0*(4d0*pi*aem)*eeqd*(gLZd-gRZd)*(fLZ-fRZ)*chi1 )
     \          -1/2d0*(4d0*pi*aem)*eeqd*(gLZd-gRZd)*(fLZ-fRZ)*chi1 )
!
     \          *(2d0*costh) 
     \           )
            sigmaij(3,-3)=sigmaij(2,-2)
            sigmaij(5,-5)=sigmaij(2,-2)
            sigmaij(-2,2)=facZ*( 
     \         ( (gLZd**2+gRZd**2)*(fLZ**2+fRZ**2) 
!     \          +4d0*(4d0*pi*aem)**2*(eeqd)**2*chi2 
!     \          -2d0*(4d0*pi*aem)*eeqd*(gLZd+gRZd)*(fLZ+fRZ)*chi1 )
     \          +1/4d0*(4d0*pi*aem)**2*(eeqd)**2*chi2 
     \          -1/2d0*(4d0*pi*aem)*eeqd*(gLZd+gRZd)*(fLZ+fRZ)*chi1 )
!
     \          *(1d0+costh**2)
     \      +  ( (gLZd**2-gRZd**2)*(fLZ**2-fRZ**2)
!     \          -2d0*(4d0*pi*aem)*eeqd*(gLZd-gRZd)*(fLZ-fRZ)*chi1 )
     \          -1/2d0*(4d0*pi*aem)*eeqd*(gLZd-gRZd)*(fLZ-fRZ)*chi1 )
!
     \          *(2d0*costh) 
     \           )
            sigmaij(-3,3)=sigmaij(-2,2)
            sigmaij(-5,5)=sigmaij(-2,2)
      endif
      return
      end
