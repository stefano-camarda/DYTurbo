C     Computes the complex N Mellin moments of pdfs
      subroutine pdfmoments(beam,N,UV,DV,US,DS,SS,GL,CH,BO)
      implicit none
      complex(8) N
      integer beam,hadron

      integer ih1,ih2
      common/collider/ih1,ih2

      complex(8) UV,DV,US,DS,SS,GL,CH,BO

      real(8) fx(-5:5)
      real(8) fxtemp(-2:2)
      real (8) mu

      real(8) x

      real(8) xmin,xmax,xc,xm,xa,xb

      include 'gauss.inc'

      integer intervals,i,j
      double precision muf

      integer approxpdf,pdfintervals
      common/opts/approxpdf,pdfintervals
      
c     select beam
      if (beam.eq.1) then
         hadron = ih1
      elseif (beam.eq.2) then
         hadron = ih2
      endif
c     factorization scale
      muf=91.1876d0
c     muf=2D0

c     boundaries of integration      
      xa = 1D-8
      xb = 1

c     initialise
      uv = 0d0
      dv = 0d0
      us = 0d0
      ds = 0d0
      ss = 0d0
      gl = 0d0
      ch = 0d0
      bo = 0d0
      

      intervals=pdfintervals
      do i=1,intervals
         xmin = xa*((xb/xa)**(real(i-1, 8)/real(intervals,8)))
         xmax = xa*((xb/xa)**(real(i,8)/real(intervals,8)))
c         xmin = xa+(xb-xa)*(i-1)/intervals
c         xmax = xa+(xb-xa)*i/intervals
         xc=0.5d0*(xmin+xmax)
         xm=0.5d0*(xmax-xmin)

         do j=1,24
            x=xc+xm*xxx24(j)
            call fdist(hadron,x,muf,fx)
            
c     call distributions1(x,fx(1),fx(2),fx(-1),fx(-2),
c     .           fx(3),fx(0),fx(4),fx(5))

c     integral_0^1{ x^(N-1) fx dx}
            uv = uv+x**(N-1)*(fx(2)-fx(-2))*www24(j)*xm
            dv = dv+x**(N-1)*(fx(1)-fx(-1))*www24(j)*xm
            us = us+x**(N-1)*(fx(-1))*www24(j)*xm
            ds = ds+x**(N-1)*(fx(-2))*www24(j)*xm
            ss = ss+x**(N-1)*(fx(-3))*www24(j)*xm
            gl = gl+x**(N-1)*(fx(0))*www24(j)*xm
            ch = ch+x**(N-1)*(fx(-4))*www24(j)*xm
            bo = bo+x**(N-1)*(fx(-4))*www24(j)*xm
         enddo
      enddo
c      print *,'doub'
c      print *,uv
c      print *,dv
c      print *,us
c      print *,ds
c      print *,ss
c      print *,gl
c      print *,ch
c      print *,bo
      
cc **************************************
cc     t = log(x) change of variable
c      uv = 0
c      dv = 0
c      us = 0
c      xa = 0d0
c      xb = -log(1D-8)
c      intervals=10000
c      do i=1,intervals
c         xmin = xa+(xb-xa)*(real(i-1, 8)/real(intervals,8))
c         xmax = xa+(xb-xa)*(real(i,8)/real(intervals,8))
c         xc=0.5d0*(xmin+xmax)
c         xm=0.5d0*(xmax-xmin)
c         do j=1,24
c            x=xc+xm*xxx24(j)
c            call fdist(hadron,exp(-x),xmu,fx)
c            fxtemp(1)=fx(2)
c            fxtemp(-1)=fx(-2)
c            fxtemp(2)=fx(1)
c            fxtemp(-2)=fx(-1)
c            fx(1)=fxtemp(1)
c            fx(-1)=fxtemp(-1)
c            fx(2)=fxtemp(2)
c            fx(-2)=fxtemp(-2)
c            fx(1) =exp(-x)*fx(1) 
c            fx(-1)=exp(-x)*fx(-1)
c            fx(2) =exp(-x)*fx(2) 
c            fx(-2)=exp(-x)*fx(-2)
c            fx(3) =exp(-x)*fx(3) 
c            fx(-3)=exp(-x)*fx(-3)
c            fx(3) =exp(-x)*fx(3) 
c            fx(-3)=exp(-x)*fx(-3)
c            fx(4) =exp(-x)*fx(4) 
c            fx(-4)=exp(-x)*fx(-4)
c            fx(5) =exp(-x)*fx(5) 
c            fx(-5)=exp(-x)*fx(-5)
c
cc     call distributions1(exp(-x),fx(1),fx(2),fx(-1),fx(-2),
cc     .           fx(3),fx(0),fx(4),fx(5))
c
c            uv = uv+exp(-(N-2+1)*x)*(fx(1)-fx(-1))*www24(j)*xm
c            dv = dv+exp(-(N-2+1)*x)*(fx(2)-fx(-2))*www24(j)*xm
c            us = us+exp(-(N-2+1)*x)*(fx(-1))*www24(j)*xm
c         enddo
c      enddo
c      print *,'logi',dv
c **************************************
      return
      end


      subroutine initmoments
c     IMPLICIT DOUBLE PRECISION (A - Z)
      implicit none
      INTEGER k,ik,NFITMAX
      COMPLEX*16 uval,dval,usea,dsea,ssea,glu,charm,bot
      COMPLEX*16 MellinH2qq,MellinH2gg,MellinH2gq
      
      COMPLEX*16 QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F,
     1     QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F
      
      COMPLEX*16 QQIP(136),QGFP(136), GQIP(136), GGIP(136), GGFP(136),
     1     NS1MIP(136), NS1PIP(136), NS1FP(136),QQ1FP(136), 
     2     QG1FP(136), GQ1IP(136), GQ1FP(136), GG1IP(136), GG1FP(136)
      
      COMPLEX*16 QQIM(136),QGFM(136), GQIM(136), GGIM(136), GGFM(136),
     1     NS1MIM(136), NS1PIM(136), NS1FM(136),QQ1FM(136), 
     2     QG1FM(136), GQ1IM(136), GQ1FM(136), GG1IM(136), GG1FM(136)
      
      COMPLEX*16 C2qgMp(136),C2NSqqMp(136),C2SqqbMp(136),
     .     C2NSqqbMp(136)
      COMPLEX*16 C2qgMm(136),C2NSqqMm(136),C2SqqbMm(136),
     .     C2NSqqbMm(136)
      COMPLEX*16 C2qg,C2NSqq,C2Sqqb,C2NSqqb
     
      COMMON / ANOMP/QQIp, QGFp, GQIp, GGIp, GGFp, NS1MIp, NS1PIp, 
     1          NS1Fp, QQ1Fp, QG1Fp, GQ1Ip, GQ1Fp, GG1Ip, GG1Fp
      COMMON / ANOMM/QQIm, QGFm, GQIm, GGIm, GGFm, NS1MIm, NS1PIm, 
     1          NS1Fm, QQ1Fm, QG1Fm, GQ1Im, GQ1Fm, GG1Im, GG1Fm
       
      COMMON / H2COEF /C2qgMp,C2NSqqMp,C2SqqbMp,C2NSqqbMp,
     .     C2qgMm,C2NSqqMm,C2SqqbMm,C2NSqqbMm

c     Common blocks of PDFs mellin moments, where the moments are stored
      complex*16 UVP(136,30),DVP(136,30),USP(136,30),DSP(136,30),
     .     SSP(136,30),GLP(136,30),CHP(136,30),BOP(136,30)
      complex*16 UVM(136,30),DVM(136,30),USM(136,30),DSM(136,30),
     .     SSM(136,30),GLM(136,30),CHM(136,30),BOM(136,30)
      complex*16 UVP2(136,30),DVP2(136,30),USP2(136,30),DSP2(136,30),
     .     SSP2(136,30),GLP2(136,30),CHP2(136,30),BOP2(136,30)
      complex*16 UVM2(136,30),DVM2(136,30),USM2(136,30),DSM2(136,30),
     .     SSM2(136,30),GLM2(136,30),CHM2(136,30),BOM2(136,30)

      common / DISTP1/ UVP,DVP,USP,DSP,SSP,GLP,CHP,BOP
      common / DISTM1/ UVM,DVM,USM,DSM,SSM,GLM,CHM,BOM
      common / DISTP2/ UVP2,DVP2,USP2,DSP2,SSP2,GLP2,CHP2,BOP2
      common / DISTM2/ UVM2,DVM2,USM2,DSM2,SSM2,GLM2,CHM2,BOM2

c     Gaussian nodes of the integration contour in the complex plane
      COMPLEX*16 CCp,CCm, Np(136),Nm(136),XN
      COMMON / MOMS2    / Np,Nm,CCP,CCm

      common/NFITMAX/NFITMAX

      Write(6,*)'Initialise PDF moments with numerical integration'
      NFITMAX = 14
c     calculate Mellin moments of PDFs
c     Beam 1        
      do k=1,88
C     positive branch
         XN=Np(k)  
         call pdfmoments(1,XN,uval,dval,usea,dsea,ssea,glu,charm,bot)
c         print *,'beam 1 positive'
c         print *,'moment',k,XN
c         print *,'uval  ',uval
c         print *,'dval  ',dval
c         print *,'usea  ',usea
c         print *,'dsea  ',dsea
c         print *,'gluon ',glu
c         print *,'charm ',charm
c         print *,'bottom',bot
         do ik=1,NFITMAX
            UVp(k,ik)= uval
            DVp(k,ik)= dval
            USp(k,ik)= usea
            DSp(k,ik)= dsea
            SSp(k,ik)= ssea
            GLp(k,ik)= glu
            CHp(k,ik)= charm
            BOp(k,ik)= bot
         enddo 
c     negative branch
         XN=Nm(k)
         call pdfmoments(1,XN,uval,dval,usea,dsea,ssea,glu,charm,bot)
         do ik=1,NFITMAX
            UVm(k,ik)= uval
            DVm(k,ik)= dval
            USm(k,ik)= usea
            DSm(k,ik)= dsea
            SSm(k,ik)= ssea
            GLm(k,ik)= glu
            CHm(k,ik)= charm
            BOm(k,ik)= bot
         enddo
      enddo

c     Beam 2
      do k=1,88
c     positive branch
         XN=Np(k)  
         call pdfmoments(2,XN,uval,dval,usea,dsea,ssea,glu,charm,bot)
         do ik=1,NFITMAX
            UVp2(k,ik)= uval 
            DVp2(k,ik)= dval 
            USp2(k,ik)= usea 
            DSp2(k,ik)= dsea 
            SSp2(k,ik)= ssea 
            GLp2(k,ik)= glu  
            CHp2(k,ik)= charm
            BOp2(k,ik)= bot  
         enddo
c     negative branch
         XN=Nm(k)
         call pdfmoments(2,XN,uval,dval,usea,dsea,ssea,glu,charm,bot)
         do ik=1,NFITMAX
            UVm2(k,ik)= uval 
            DVm2(k,ik)= dval 
            USm2(k,ik)= usea 
            DSm2(k,ik)= dsea 
            SSm2(k,ik)= ssea 
            GLm2(k,ik)= glu  
            CHm2(k,ik)= charm
            BOm2(k,ik)= bot  
         enddo
      enddo

c     calculate Mellin moments of anomalous dimensions
      do k=1,88
c     positive branch
         XN=Np(k)  
         CALL ANCALC (QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F,
     1        QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, XN)
         QQIp(k) = QQI 
         QGFp(k) = QGF 
         GQIp(k) = GQI 
         GGIp(k) = GGI 
         GGFp(k) = GGF 
         NS1MIp(k) = NS1MI  
         NS1PIp(k) = NS1PI 
         NS1Fp(k) = NS1F            
         QQ1Fp(k) = QQ1F 
         QG1Fp(k) = QG1F 
         GQ1Ip(k) = GQ1I 
         GQ1Fp(k) = GQ1F 
         GG1Ip(k) = GG1I 
         GG1Fp(k) = GG1F 
         
         call H2calc(C2qg,C2NSqqb,C2NSqq,C2Sqqb,XN)

         C2qgMp(k)= C2qg
         C2NSqqMp(k)= C2NSqq
         C2SqqbMp(k)= C2Sqqb
         C2NSqqbMp(k)=  C2NSqqb      

c     negative branch
         XN=Nm(k)
         CALL ANCALC (QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F,
     1        QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, XN)
         QQIm(k) = QQI 
         QGFm(k) = QGF 
         GQIm(k) = GQI 
         GGIm(k) = GGI 
         GGFm(k) = GGF
         NS1MIm(k) = NS1MI  
         NS1PIm(k) = NS1PI 
         NS1Fm(k) = NS1F            
         QQ1Fm(k) = QQ1F 
         QG1Fm(k) = QG1F 
         GQ1Im(k) = GQ1I 
         GQ1Fm(k) = GQ1F 
         GG1Im(k) = GG1I 
         GG1Fm(k) = GG1F 

         call H2calc(C2qg,C2NSqqb,C2NSqq,C2Sqqb,XN)

         C2qgMm(k)= C2qg
         C2NSqqMm(k)= C2NSqq
         C2SqqbMm(k)= C2Sqqb
         C2NSqqbMm(k)=  C2NSqqb              
      enddo

      call cacheanom
      Write(6,*)'End initialization'
      return
      end
