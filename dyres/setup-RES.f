      subroutine setup
      implicit none
      include 'constants.f'
      include 'virtonly.f'
      include 'realonly.f'
      include 'noglue.f'
      include 'lc.f'
      include 'cutoff.f'
      include 'maxwt.f'
      include 'masses.f'
      include 'process.f'
      include 'scale.f'
      include 'facscale.f'
      include 'zerowidth.f'
      include 'removebr.f'
      include 'clustering.f'
      include 'gridinfo.f'
      include 'limits.f'
      include 'alfacut.f'
      include 'pdfiset.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'nlooprun.f'
      include 'rescoeff.f'
      include 'jetcuts.f'
      include 'flags.f'
      include 'dipolescale.f'
CC
      include 'vegas_common.f'
      include 'nwz.f'
      include 'ewcharge.f'
      include 'dynamicscale.f'
      include 'lhapdf.f'
      character *2 plabel(mxpart)
      common/plabel/plabel
      integer order,notag,nqcdjets,nqcdstart,isub,nproc,ndec,flag,phot
     .,nd
      common/nnlo/order
      common/notag/notag
      common/nqcdjets/nqcdjets,nqcdstart
      common/isub/isub
      double precision BrnRat,gamgambr,wwbr,zzbr,br0,vmass
      common/BrnRat/BrnRat
      common/nproc/nproc
!      common/mlimits/wsqmin,wsqmax
      common/flag/flag
      common/dyphoton/phot
!
      double precision xqtcut
      common/qtcut/xqtcut

!resummation scale
      COMMON/a_param/a_param,b0p
      COMMON/g_param/g_param
      double precision a_param,b0p,g_param

      double precision beta1,H2qqdelta,H2qqD0
      common/Hstcoeff/beta1,H2qqdelta,H2qqD0

      logical int
      common/int/int

      logical dorebin
      common/dorebin/dorebin

      character *12 format_string
      character *30  runstring2
      character *30  runstring3
      character *30  runstring4
      character *30  runstring5
      character *30  runstring6

      character *50 prefix
      character *36 pdfstring
      integer nset
      common/prefix/nset,prefix
      common/pdfstring/pdfstring

      character*4 part
      character*30 runstring
      integer j
      logical makecuts
      integer nmin,nmax,n2,n3
      integer ih1,ih2,itmx1,itmx2,ncall1,ncall2,idum,rseed
      double precision rtsmin,sroot,LT
      double precision Mwmin,Mwmax
      common/mwminmax/Mwmin,Mwmax
      double precision Rcut
      double precision ran2,randummy
      double precision cmass,bmass
      double precision mass2,width2,mass3,width3
      double precision amz,dyalphas
      double precision brwen,brzee,brtau,brtop
      
      character *3 str1
      character *10 str2
      character *38 str3
      character *50 string
      character *30 mycard
      character *30 mytopfile

      integer pr,mysw
      common/pr/pr
      common/mysw/mysw


      common/couple/amz
      
      common/breit/n2,n3,mass2,width2,mass3,width3
      
      common/nmin/nmin
      common/nmax/nmax
      common/rtsmin/rtsmin
 

      common/part/part
      common/runstring/runstring
      common/energy/sroot
      common/density/ih1,ih2
      common/iterat/itmx1,ncall1,itmx2,ncall2
      common/ranno/idum

      
      common/Rcut/Rcut
      common/makecuts/makecuts

      common/qmass/cmass,bmass

      common/rseed/rseed
      save /ranno/

      logical lhapdfs
      common/lhapdfs/lhapdfs

!      a_param=2d0   ! resummation scale
      

      flag=0

      lhapdfs=.false.

      virtonly=.false.
      realonly=.false.

      noglue=.false.
      ggonly=.false.
      gqonly=.false.

      nmin=1
      nmax=2

      clustering=.true.
      colourchoice=0
      rtsmin=40d0
      cutoff=0.001d0

CC    

      Qflag=.true.
      Gflag=.true.

CC
!      aii=1d-2
!      aif=1d-2
!      afi=1d-2
!      aff=1d-2

      aii=1d0
      aif=1d0
      afi=1d0
      aff=1d0


      inclusive=.true.
      algorithm='none'
 
CC   Parameters used to define jets

      ptjetmin=0d0
      etajetmin=0d0
      etajetmax=20d0

      Rcut=0.4d0

CC    Dynamic scale (if true muf=mur=q)

      dynamicscale=.false.

      removebr=.false.
      makecuts=.true.     


CC    Adjust the grid at each iteration

      dorebin=.true.


CC    Read a previously saved grid

!      readin=.true.
      readin=.true.



      writeout=.false.
      ingridfile='vegas'
      outgridfile='vegas'

CC    Read inputfile

C      read(*,*) sroot
C      read(*,*) ih1,ih2
C      read(*,*) nproc !decay mode
C      read(*,*) scale,facscale  ! mur,muf!      read(*,*) a_param
C      read(*,*) a_param
C      read(*,*) g_param
C      read(*,*) order        
C      read(*,*) part
C      read(*,*) zerowidth
C      read(*,*) Mwmin,Mwmax
C      read(*,*) itmx1,ncall1
C      read(*,*) itmx2,ncall2
C      read(*,*) rseed          ! seed
C      read(*,*) iset,nset     
C      read(*,*) PDFname,PDFmember
C      read(*,*) runstring
C      read(*,*) pr

CC    Set all factorization scales to facscale
CC    to avoid problems when dynamicscale=.false.

      do nd=0,40
       dipscale(nd)=facscale
      enddo


      if(pr.ne.0) then! only if we want to write "ntuples"
cc      if(order.eq.0) runstring3=trim(runstring)//'_born.evt'
cc      if(order.ge.1) then
        runstring3=trim(runstring)//'_res.evt'
        runstring4=trim(runstring)//'_fin.evt'
      open (unit=121, file=runstring4, action='write', status='unknown')
cc      endif
      open (unit=120, file=runstring3, action='write', status='unknown')

      mycard=trim(runstring)//'_evt'
      mytopfile=trim(runstring)//'_top'
      open(unit=200,file=mycard,action='write',status='unknown')
      write(200,"(I1,6x,A24)") order,"[order]"
      write(200,"(I2,6x,A24)") pr,"[iterations recorded]"
cc      write(200,"(I2,6x,A24)") itmx2,"[iterations acc. run]"
      write(200,"(A4,6x,A24)") part,"[part]"
      write(200,*) mytopfile, "  [generated top file]"
      write(200,*) runstring3, "  [event file for resummed]"
      if(order.ge.1) write(200,*) runstring4, " [event file for finite]"
      close(200)
      endif


CC    Cut on qt/Q
      xqtcut=0.008d0  !qtcut=1d0  

C     Limits on invariant mass of vector boson decay products
C     (irrelevant if zerowidth=.true.)


      wsqmin=Mwmin**2
      wsqmax=Mwmax**2

C     Check if the limits are compatible with sroot
      if(wsqmax.gt.(sroot**2)) wsqmax=sroot**2

      do j=1,mxpart
      plabel(j)=''
      enddo


      plabel(1)='pp'
      plabel(2)='pp'

c--- the default behaviour is to remove no branching ratio

      BrnRat=1d0

      call dycoupling

      call cstring(pdfstring)

C      if(lhapdfs.eqv.(.false.)) then
C       write(6,*)'CCCCCCCCCC      Parton Distributions    CCCCCCCCCCCC'
C       write(6,*)'C                                                  C'
C       write(6,*)'C       ', pdfstring,'       C'
C       write(6,*)'C                                                  C'
C      endif
C      write(6,*)'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'

      notag=0
      nqcdjets=0
      isub=0


      
      mb=0
    
    
     
      if((ih1.eq.1).and.(ih2.eq.-1)) then
       str1='ppb'
      elseif((ih1.eq.1).and.(ih2.eq.1)) then
       str1='pp'
      else
       write(*,*)'Initial state not allowed'
      endif


CC    Here ndim is the number of dimensions for Born

      if(nproc.eq.1) then
C
C     W+->l+nubar
C
       str2=' -> W+ -> ' 
       ndec=2
       ndim=4    
       plabel(3)='nl'
       plabel(4)='ea'
       plabel(5)='pp'
       plabel(6)='pp'
       n2=0
       n3=1
       nwz=1
       mass3=wmass
       width3=wwidth

       if (removebr) then
        call branch(brwen,brzee,brtau,brtop)
        BrnRat=brwen
       endif

      str3=' nu(p3)+e+(p4)'

      elseif(nproc.eq.2) then
C
C      W-=>l-nu
C
       str2=' -> W- -> ' 
       ndec=2
       ndim=4
       plabel(3)='el'
       plabel(4)='na'
       plabel(5)='pp'
       plabel(6)='pp'
       n2=0
       n3=1
       nwz=-1
       mass3=wmass
       width3=wwidth
     
       
     

       if (removebr) then
        call branch(brwen,brzee,brtau,brtop)
        BrnRat=brwen
       endif

      str3=' e^-(p3)+nu~(p4)'

      elseif(nproc.eq.3) then
C
C      Z->e+e-
C
       str2=' -> Z -> ' 
       ndec=2
       ndim=4
       plabel(3)='el'
       plabel(4)='ea'
       plabel(5)='pp'
       plabel(6)='pp'
       n2=0
       n3=1
       nwz=0
       mass3=zmass
       width3=zwidth

       l1=le     
       r1=re


C       q1=-1
!       q1=0 !switch off the photon

       phot=q1

       int=.false.

       str3=' e-(p3)+e+(p4)'
       
     

       if (removebr) then
        call branch(brwen,brzee,brtau,brtop)
        BrnRat=brzee
       endif




      else
       write(*,*)'Wrong decay channel'
       stop
      endif


!     qtcut max
!      qtcutmax=5d0*mass3

C     New: decide if using Breit-Wigner or not

      n3=0
      vmass=wmass
      if(nproc.eq.3) vmass=zmass

      if(mwmin.lt.vmass.and.mwmax.gt.vmass) n3=1

CCCCCCCCCCCCC

      call strcat(str1,str2,string)
      call strcat(string,str3,string)

     
      call cstring(string)      
     

      write(6,*)'C                                                  C'

      if(order.eq.0) then
!      write(6,*)'C         Computing LO cross section for           C'
      write(6,*)'WRONG ORDER (0): CHOOSE NLL+LO (1) OR NNLL+NLO (2)'
      stop
      elseif(order.eq.1) then
      write(6,*)'C        Computing NLO cross section for           C'
      elseif(order.eq.2) then
      write(6,*)'C        Computing NNLO cross section for          C'
      else
      write(*,*)'Order can be 0,1 or 2 !'
       stop
      endif
  
      write(6,*)'C                                                  C' 
      write(*,*)'C',string,'C'
      write(6,*)'C                                                  C' 
      write(6,96)sroot
      write(6,*)'C                                                  C'
      write(6,*)'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC' 


 96   format(' C           at  Sqrt(s)=',f8.2,' GeV               C')

      nqcdjets=0
    

      call ckmfill(nwz)

    

CCCCCCCCCCCCCCCCCC


c--- set-up the random number generator with a negative seed
      idum=-abs(rseed)
      randummy=ran2()

c--- initialize masses for alpha_s routine
      cmass=dsqrt(mcsq)
      bmass=dsqrt(mbsq)


c--- check that we have a valid value of 'part'
      if ( (part .ne. 'lord') .and. (part .ne. 'real') .and.
     .     (part .ne. 'virt') .and. (part .ne. 'tota') ) then
          write(6,*) 'part=',part,' is not a valid option'
          stop     
      endif      



        as=dyalphas(scale,amz,nlooprun)
        ason2pi=as/twopi
        ason4pi=as/fourpi
        gsq=fourpi*as
        musq=scale**2

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CC    Resummation coefficients

      beta0=(33-2*nf)/12d0

      beta1=(153d0-19*nf)/24d0

      Kappa=67/6d0-(pi**2)/2d0-5d0/9d0*nf

      A1q=4d0/3
      A2q=0.5d0*A1q*Kappa
      B1q=-2d0


      B2q=4d0/9*(pi**2-3d0/4-12*Z3)+(11d0/9*pi**2-193d0/12+6*Z3)
     & +nf/6d0*(17d0/3-4d0/9*pi**2)

      

C     Delta term in c1qq coefficient

      C1qqdelta=(pi**2-8)/3d0

C     Delta term in P2qq splitting function (as/pi normalization)

      Delta2qq=16d0/9*(3d0/8-pi**2/2+6*Z3)
     &   +4*(17d0/24+11d0*pi**2/18-3*Z3)-2d0/3*nf*(1d0/6+2*pi**2/9d0)

      Delta2qq=Delta2qq/4d0

CC    Coefficients of D0 and D1 in P*P (as/pi normalization)

      D0qqqq=8d0/3
      D1qqqq=32d0/9


CC    Coefficients of delta(1-z) in P*P

      Deltaqqqq=4d0/9*(9d0/4-2*pi**2/3d0)

C     H2qq contribution: coefficient of delta(1-z)

      H2qqdelta=-2561d0/144+127d0*nf/72+3*pi**2/2-19d0*nf*Pi**2/81+
     &          49d0*Pi**4/324 +58d0*Z3/9 + 8d0*nf*Z3/27

C     H2qq contribution: coefficient of D0(z)

      H2qqD0=-404d0/27+(56d0*nf)/81+14*Z3


      return

      end
      
