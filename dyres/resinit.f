C...  ANOMALOUS DIMENSIONS FOR LEADING AND NEXT TO LEADING ORDER
C...  EVOLUTION OF PARTON DENSITIES AND WILSON COEFFICIENTS FOR
C...  NLO STRUCTURE FUNCTIONS :

c************************************************
c     all the output values are function of I (z point) and are cached
c************************************************

c     settings: f = 5 (flavours)
c     inputs: QQI, QGF, GQI, GGI, GGF (only function of I)
c     outputs ans, am, ap, al, be, ab are also only function of I, they are cached together with the inputs
c     outputs rmin, rplus, rqq, rgg same as above
      SUBROUTINE CACHEANOM
c       IMPLICIT DOUBLE PRECISION (A-I,L-Z)
      IMPLICIT none
c     Input from INITO and INITOFIT
      COMPLEX*16 CCp,CCm, Np(136),Nm(136),XNN
      COMMON / MOMS2    / Np,Nm,CCP,CCm
      COMPLEX*16 QQIP(136),QGFP(136), GQIP(136), GGIP(136), GGFP(136),
     1     NS1MIP(136), NS1PIP(136), NS1FP(136),QQ1FP(136), 
     2     QG1FP(136), GQ1IP(136), GQ1FP(136), GG1IP(136), GG1FP(136) 
      COMPLEX*16 QQIM(136),QGFM(136), GQIM(136), GGIM(136), GGFM(136),
     1     NS1MIM(136), NS1PIM(136), NS1FM(136),QQ1FM(136), 
     2     QG1FM(136), GQ1IM(136), GQ1FM(136), GG1IM(136), GG1FM(136)
      COMMON / ANOMP/QQIp, QGFp, GQIp, GGIp, GGFp, NS1MIp, NS1PIp, 
     1     NS1Fp, QQ1Fp, QG1Fp, GQ1Ip, GQ1Fp, GG1Ip, GG1Fp
      COMMON / ANOMM/QQIm, QGFm, GQIm, GGIm, GGFm, NS1MIm, NS1PIm, 
     1     NS1Fm, QQ1Fm, QG1Fm, GQ1Im, GQ1Fm, GG1Im, GG1Fm
      COMPLEX*16 C2qgMp(136),C2NSqqMp(136),C2SqqbMp(136),
     1     C2NSqqbMp(136)
      COMPLEX*16 C2qgMm(136),C2NSqqMm(136),C2SqqbMm(136),
     1     C2NSqqbMm(136)
      COMMON / H2COEF /C2qgMp,C2NSqqMp,C2SqqbMp,C2NSqqbMp,
     1     C2qgMm,C2NSqqMm,C2SqqbMm,C2NSqqbMm

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

c     temporary variables      
      DOUBLE PRECISION B04, B1, B02, B10
      INTEGER N, F, S, I
      DOUBLE COMPLEX QQ, QG, GQ, GG, SQ, GP, GM
      DOUBLE COMPLEX QQI,QGF,GQI,GGI,GGF,
     1     NS1MI,NS1PI,NS1F,QQ1F,QG1F,GQ1I,GQ1F,GG1I,GG1F
      DOUBLE COMPLEX NS1M,NS1P,QQ1,QG1,GQ1,GG1

      double precision pi
      parameter(pi=3.14159265358979d0)
      integer nf
c      include 'const.h' 
c      include 'constants.f' 

C...  ANOMALOUS DIMENSIONS AND RELATED QUANTITIES IN LEADING ORDER :
      F = 5
      NF = F
      N = 1
      B04 = 11.- 2./3.* F
      B02 = 2.* B04
      do s=1,2
         do I=1,136
            if (s.eq.1) then
               QQI = QQIp(I)                    
               QGF = QGFp(I) 
               GQI = GQIp(I) 
               GGI = GGIp(I) 
               GGF = GGFp(I) 
               NS1MI = NS1MIp(I)  
               NS1PI = NS1PIp(I) 
               NS1F = NS1Fp(I)            
               QQ1F = QQ1Fp(I) 
               QG1F = QG1Fp(I) 
               GQ1I = GQ1Ip(I) 
               GQ1F = GQ1Fp(I) 
               GG1I = GG1Ip(I) 
               GG1F = GG1Fp(I) 
            else
               QQI = QQIm(I) 
               QGF = QGFm(I)
               GQI = GQIm(I)
               GGI = GGIm(I)
               GGF = GGFm(I)
               NS1MI = NS1MIm(I)
               NS1PI = NS1PIm(I)
               NS1F = NS1Fm(I)            
               QQ1F = QQ1Fm(I) 
               QG1F = QG1Fm(I) 
               GQ1I = GQ1Im(I) 
               GQ1F = GQ1Fm(I) 
               GG1I = GG1Im(I) 
               GG1F = GG1Fm(I) 
            endif      
            QQ = QQI
            QG = F * QGF
            GQ = GQI
            GG = GGI + F * GGF
            SQ = SQRT ((GG - QQ) * (GG - QQ) + 4.* QG * GQ)
            GP = 0.5 * (QQ + GG + SQ)
            GM = 0.5 * (QQ + GG - SQ)
            cANS(I,s) = QQ / B02
            cAM(I,s) = GM / B02
            cAP(I,s) = GP / B02
            cAL(I,s) = (QQ - GP) / (GM - GP)
            cBE(I,s) = QG / (GM - GP)
            cAB(I,s) = GQ / (GM - GP)
C...  NEXT TO LEADING ORDER : ANOMALOUS DIMENSIONS AND WILSON COEFFICIENTS
C...  IN THE MS-BAR FACTORIZATION SCHEME OF BARDEEN ET AL. (1981) :
            NS1M = NS1MI + F * NS1F
            NS1P = NS1PI + F * NS1F
            QQ1 = NS1P + F * QQ1F
            QG1 = F * QG1F
            GQ1 = GQ1I + F * GQ1F
            GG1 = GG1I + F * GG1F
C...  COMBINATIONS OF ANOMALOUS DIMENSIONS FOR NLO - SINGLET EVOLUTION :
            B1 = 102 - 38./3.* F
            B10 = B1 / B04
            cRMIN(I,s) = (NS1M - QQ * B10) / B02
            cRPLUS(I,s) = (NS1P - QQ * B10) / B02
            cRQQ(I,s) = (QQ1 - QQ * B10) / B02
            cRQG(I,s) = (QG1 - QG * B10) / B02
            cRGQ(I,s) = (GQ1 - GQ * B10) / B02
            cRGG(I,s) = (GG1 - GG * B10) / B02

c            additional caching for reno2
      cAC(I,s)  = 1.- cAL(I,s)
      cNMP(I,s) = 1.- cAM(I,s) + cAP(I,s)
      cNPM(I,s) = 1.- cAP(I,s) + cAM(I,s)
      cDMQQ(I,s) =  cAL(I,s) * cRQQ(I,s) + cBE(I,s) * cRGQ(I,s)
      cDMQG(I,s) =  cAL(I,s) * cRQG(I,s) + cBE(I,s) * cRGG(I,s)
      cDMGQ(I,s) =  cAB(I,s) * cRQQ(I,s) + cAC(I,s) * cRGQ(I,s)
      cDMGG(I,s) =  cAB(I,s) * cRQG(I,s) + cAC(I,s) * cRGG(I,s)
      cDPQQ(I,s) =  cAC(I,s) * cRQQ(I,s) - cBE(I,s) * cRGQ(I,s)
      cDPQG(I,s) =  cAC(I,s) * cRQG(I,s) - cBE(I,s) * cRGG(I,s)
      cDPGQ(I,s) = -cAB(I,s) * cRQQ(I,s) + cAL(I,s) * cRGQ(I,s)
      cDPGG(I,s) = -cAB(I,s) * cRQG(I,s) + cAL(I,s) * cRGG(I,s)
      cRMMQQ(I,s)=  cAL(I,s)*cDMQQ(I,s)+cAB(I,s)*cDMQG(I,s)
      cRMMQG(I,s)=  cBE(I,s)*cDMQQ(I,s)+cAC(I,s)*cDMQG(I,s)
      cRMMGQ(I,s)=  cAL(I,s)*cDMGQ(I,s)+cAB(I,s)*cDMGG(I,s)
      cRMMGG(I,s)=  cBE(I,s)*cDMGQ(I,s)+cAC(I,s)*cDMGG(I,s)
      cRMPQQ(I,s)= (cAC(I,s)*cDMQQ(I,s)-cAB(I,s)*cDMQG(I,s))/cNMP(I,s)
      cRMPQG(I,s)=(-cBE(I,s)*cDMQQ(I,s)+cAL(I,s)*cDMQG(I,s))/cNMP(I,s)
      cRMPGQ(I,s)= (cAC(I,s)*cDMGQ(I,s)-cAB(I,s)*cDMGG(I,s))/cNMP(I,s)
      cRMPGG(I,s)=(-cBE(I,s)*cDMGQ(I,s)+cAL(I,s)*cDMGG(I,s))/cNMP(I,s)
      cRPMQQ(I,s)= (cAL(I,s)*cDPQQ(I,s)+cAB(I,s)*cDPQG(I,s))/cNPM(I,s)
      cRPMQG(I,s)= (cBE(I,s)*cDPQQ(I,s)+cAC(I,s)*cDPQG(I,s))/cNPM(I,s)
      cRPMGQ(I,s)= (cAL(I,s)*cDPGQ(I,s)+cAB(I,s)*cDPGG(I,s))/cNPM(I,s)
      cRPMGG(I,s)= (cBE(I,s)*cDPGQ(I,s)+cAC(I,s)*cDPGG(I,s))/cNPM(I,s)
      cRPPQQ(I,s)=  cAC(I,s)*cDPQQ(I,s)-cAB(I,s)*cDPQG(I,s)
      cRPPQG(I,s)= -cBE(I,s)*cDPQQ(I,s)+cAL(I,s)*cDPQG(I,s)
      cRPPGQ(I,s)=  cAC(I,s)*cDPGQ(I,s)-cAB(I,s)*cDPGG(I,s)
      cRPPGG(I,s)= -cBE(I,s)*cDPGQ(I,s)+cAL(I,s)*cDPGG(I,s)


      
c cache C1 coefficients      
      if (s.eq.1) then
         XNN = Np(I)
      else
         XNN = Nm(I)
      endif
      cC1QG(I,s)=1d0/((XNN+1)*(XNN+2))
      cC1GQ(I,s)=4/3d0/(XNN+1)
      cC1qq(I,s)=2*pi**2/3d0-16/3d0+4/3d0/(XNN*(XNN+1))
      cC1gg=pi**2/2d0+11/2d0+pi**2

c     cache gamma1 gamma2: NORMALIZED ANOMALOUS DIMENSIONS AND COEFFICIENTS
      cgamma1qq(I,s)=-1d0*(QQI/4d0)
      cgamma1qg(I,s)=-1d0*(QGF/8d0)
      cgamma1gq(I,s)=-1d0*(GQI/4d0)
      cgamma1gg(I,s)=-1d0*((GGI+nf*GGF)/4d0)
      cgamma2qq(I,s)=-1d0*(((NS1PI+nf*NS1F)+nf*QQ1F)/8d0)
      cgamma2qqV(I,s)=-1d0*(NS1PI+2*nf*NS1F+NS1MI)/16d0
      cgamma2qqbV(I,s)=-1d0*(NS1PI-NS1MI)/16d0
      cgamma2qqS(I,s)=-1d0*(QQ1F/16d0)
      cgamma2qqbS(I,s)=cgamma2qqS(I,s)
      cgamma2qg(I,s)=-1d0*(QG1F/16d0)
      cgamma2gq(I,s)=-1d0*((GQ1I+nf*GQ1F)/8d0)
      cgamma2gg(I,s)=-1d0*((GG1I+nf*GG1F)/8d0)

      if (s.eq.1) then
         cC2qgM(I,s)    = C2qgMp(I)
         cC2NSqqM(I,s)  = C2NSqqMp(I)
         cC2SqqbM(I,s)  = C2SqqbMp(I)
         cC2NSqqbM(I,s) = C2NSqqbMp(I)
      else
         cC2qgM(I,s)    = C2qgMm(I)
         cC2NSqqM(I,s)  = C2NSqqMm(I)
         cC2SqqbM(I,s)  = C2SqqbMm(I)
         cC2NSqqbM(I,s) = C2NSqqbMm(I)
      endif

         enddo
      enddo

      return
      end

      SUBROUTINE cachecoeff
      implicit none
c     Input from resumm
       double complex loga,logmuf2q2,logq2muf2,logq2mur2
       common/clogs/loga,logmuf2q2,logq2muf2,logq2mur2
       double precision aass
       COMMON/aass/aass

c     input from cacheanom
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

      DOUBLE COMPLEX H1q
c     cached coefficients
      DOUBLE COMPLEX cHqqb(136,136,2)
      DOUBLE COMPLEX cH1stqqb(136,136,2),cH1stqg(136,2),cH1stgg,
     1     cH2stqqp(136,2),cH2stqq(136,2),cH2stqqb(136,136,2),
     2     cH2stqg_1(136,136,2),cH2stqg_2(136,136,2),cH2stgg(136,136,2)
      COMMON /ccoefficients/cHqqb,
     1     cH1stqqb,cH1stqg,cH1stgg,
     2     cH2stqqp,cH2stqq,cH2stqqb,cH2stqg_1,cH2stqg_2,cH2stgg

      INTEGER S, I1, I2, I
      integer flag1
      COMMON/flag1/flag1

c     cached from resumm
       DOUBLE PRECISION ax1,ax2,xx1,xx2
       integer nmax1,nmax2
       common /cxx/ax1,ax2,xx1,xx2,nmax1,nmax2

       integer nmax

      include 'const.h' 

      H1q=dcmplx(0d0,0d0)

      nmax = max(nmax1,nmax2)
c      nmax = 88 !136

      if(flag1.eq.1) then
      do s=1,2
         do I1=1,nmax
            do I2=1,nmax
      cHqqb(I1,I2,s)=1d0+aass/2d0*(H1q+cC1qq(I1,1)+cC1qq(I2,s))
     .   -aass/2d0*(cgamma1qq(I1,1)+cgamma1qq(I2,s))*(logmuf2q2+2*loga)
     .   +aass/2d0*(-4*loga)*(B1q+A1q*loga)
      enddo
      enddo
      enddo
      endif


      if(flag1.eq.2) then
         cH1stgg=0d0
         do s=1,2
            do I=1,nmax
               cH1stqg(I,s) = cC1qg(I,s) 
     .              +(-cgamma1qg(I,s))*(logmuf2q2+2*loga)
               cH2stqq(I,s)= cC2NSqqbM(I,s)+ cC2SqqbM(I,s)
               cH2stqqp(I,s)= cC2SqqbM(I,s)


C  qq(I,s)  means   Qb Qb -> Qb Q =  Q Q -> Q Qb
      cH2stqq(I,s)=cH2stqq(I,s) + 4d0*(
     / -2d0*loga*(
     / +(1/4d0*(cgamma2qqbV(I,s)+cgamma2qqbS(I,s))))
     / +(1/4d0*(cgamma2qqbV(I,s)+cgamma2qqbS(I,s)))*logq2muf2
     /+1/2d0*(cH1stqg(I,s)+cC1qg(I,s))/2d0
     /     *(1/2d0*(cgamma1gq(I,s))*(logq2muf2
     / -2d0*loga)
     /))     


C  qqp_2  means   Q Q' -> Q Qb  flavor in sigmaQQb determined by "first parton"
      cH2stqqp(I,s)=cH2stqqp(I,s)+4d0*(
     / -2d0*loga*(
     / +((1/4d0*cgamma2qqbS(I,s))))  
     / +((1/4d0*cgamma2qqbS(I,s)))*logq2muf2
     /+1/2d0*(cH1stqg(I,s)+cC1qg(I,s))/2d0
     /     *(1/2d0*(cgamma1gq(I,s))*(logq2muf2
     / -2d0*loga)
     /))
      enddo
      enddo

         do s=1,2
            do I1=1,nmax
               do I2=1,nmax

       cH1stqqb(I1,I2,s)=(H1q+cC1qq(I1,1)+cC1qq(I2,s))
     .   -(cgamma1qq(I1,1)+cgamma1qq(I2,s))*(logmuf2q2+2*loga)
     .   +(-4*loga)*(B1q+A1q*loga)

C  All *4 because of normalization (as/2pi)^2 en lugar de as/pi
C normalization of as/2pi  implies  gamma^1->gamma^1/2
C                                   C^1, H^1 -> C^1/2, H^1/2 
C                                   gamma^2-> gamma^2/4
C                              beta,A,B -> beta,A,B are in as/pi already


      cH2stqqb(I1,I2,s) = cC2NSqqM(I1,1) 
     .     + cC2NSqqM(I2,s) + cC2SqqbM(I1,1)+ cC2SqqbM(I2,s)+
     .            cC1qq(I1,1) * cC1qq(I2,s)

       cH2stqqb(I1,I2,s) = cH2stqqb(I1,I2,s) + 4d0*(
     .  + 1d0/6d0*A1q*beta0*8*loga**3  
     .  + 1d0/2d0*4*loga**2*(A2q-beta0*(B1q+2*A1q*loga 
     .                  +cgamma1qq(I1,1)/2d0 +cgamma1qq(I2,s)/2d0 ))     
     .  - 2*loga*(B2q+2*A2q*loga-beta0*(cC1qq(I1,1)+cC1qq(I2,s))/2d0 
     .    + cgamma2qqV(I1,1)/4d0 + cgamma2qqS(I1,1)/4d0 
     .    + cgamma2qqV(I2,s)/4d0 + cgamma2qqS(I2,s)/4d0 )      
     .  + beta0/2d0*(cgamma1qq(I1,1)+cgamma1qq(I2,s))/2d0*logq2muf2**2
     .    + (cgamma2qqV(I1,1)/4d0 + cgamma2qqS(I1,1)/4d0 
     .    + cgamma2qqV(I2,s)/4d0 + cgamma2qqS(I2,s)/4d0 )           
     .   *logq2muf2
     .  - cH1stqqb(I1,I2,s)/2d0*beta0*logq2mur2
     .  + 1d0/2d0*(cH1stqqb(I1,I2,s)+H1q+cC1qq(I1,1)+cC1qq(I2,s))/2d0*(
     .       (cgamma1qq(I1,1)+cgamma1qq(I2,s))/2d0*(logq2muf2-2d0*loga)
     .        -((B1q+A1q*loga)*2*loga ))
     . +1d0/4d0*(cH1stqg(I1,1)+cC1qg(I1,1))
     .     *cgamma1gq(I1,1)/2*(logq2muf2-2d0*loga)
     . +1d0/4d0*(cH1stqg(I2,s)+cC1qg(I2,s))
     .     *cgamma1gq(I2,s)/2*(logq2muf2-2d0*loga)
     .          ) 

       cH2stqg_1(I1,I2,s) = cC2qgM(I1,1) + cC1qg(I1,1) * cC1qq(I2,s) 
     .      + 4*(  
     .      + 1d0/2d0*beta0*4*loga**2*(-cgamma1qg(I1,1)/2d0)
     .      - 2*loga*(-beta0 * cC1qg(I1,1)/2d0 + cgamma2qg(I1,1)/4d0)
     .      + 1d0/2d0*beta0*logq2muf2**2*(cgamma1qg(I1,1)/2d0)
     .      + cgamma2qg(I1,1)/4d0*logq2muf2
     .      - beta0*logq2mur2*cH1stqg(I1,1)/2d0
     .      + 1d0/2d0*(cH1stqqb(I1,I2,s) 
     .      + H1q + cC1qq(I1,1) + cC1qq(I2,s))/2d0*
     .             (logq2muf2-2*loga) * cgamma1qg(I1,1)/2d0
     .      + 1d0/2d0 *  ( cH1stqg(I1,1) + cC1qg(I1,1))/2d0 *
     .       ( (logq2muf2-2*loga) *(cgamma1qq(I2,s) 
     .      + cgamma1gg(I1,1) )/2d0-((B1q+A1q*loga)*2*loga ) ) )

       cH2stqg_2(I1,I2,s) = cC2qgM(I2,s) + cC1qg(I2,s) * cC1qq(I1,1) 
     .      + 4*(  
     .      + 1d0/2d0*beta0*4*loga**2*(-cgamma1qg(I2,s)/2d0)
     .      - 2*loga*(-beta0 * cC1qg(I2,s)/2d0 + cgamma2qg(I2,s)/4d0)
     .      + 1d0/2d0*beta0*logq2muf2**2*(cgamma1qg(I2,s)/2d0)
     .      + cgamma2qg(I2,s)/4d0*logq2muf2
     .      - beta0*logq2mur2*cH1stqg(I2,s)/2d0
     .      + 1d0/2d0*(cH1stqqb(I1,I2,s) 
     .      + H1q + cC1qq(I2,s) + cC1qq(I1,1))/2d0*
     .             (logq2muf2-2*loga) * cgamma1qg(I2,s)/2d0
     .      + 1d0/2d0 *  ( cH1stqg(I2,s) + cC1qg(I2,s))/2d0 *
     .      ( (logq2muf2-2*loga) *(cgamma1qq(I1,1) 
     .      + cgamma1gg(I2,s) )/2d0-((B1q+A1q*loga)*2*loga ))  )

C  GG done
       cH2stgg(I1,I2,s) = cC1qg(I2,s)*cC1qg(I1,1)
     .      -4*( 1d0/2d0 * (logmuf2q2+2*loga)*
     .   ( (cH1stqg(I1,1) + cC1qg(I1,1) )/2d0*cgamma1qg(I2,s)/2d0 + 
     .     (cH1stqg(I2,s) + cC1qg(I2,s) )/2d0*cgamma1qg(I1,1)/2d0 ) )


      enddo
      enddo
      enddo
      endif


      return
      end

c**************************************
c All these constants must be initialised only once
      subroutine reno2const
      implicit none
      DOUBLE COMPLEX BETA0N,BETA1N,BETA2N,
     1     A1GN,A2GN,A3GN,
     2     B1GN,B2GN,H1G
      COMMON/COEFF/BETA0N,BETA1N,BETA2N,A1GN,A2GN,A3GN,B1GN,B2GN,H1G
      include 'const.h' 
c.....this is because SIGMA IS NORMALIZED TO AS/2*PI
c.....while our coefficients are normalized to AS/PI
c**************************************
      beta0N=beta0*2d0
      beta1N=beta1*4d0
      beta2N=beta2*8d0
      A1gN=A1g*2d0
      A2gN=A2g*4d0
      A3gN=A3g*8d0
      B1gN=B1g*2d0
      B2gN=B2g*4d0
      end
c**************************************

c**************************************
c All these constants must be initialised only once
      subroutine resummconst
      implicit none
      include 'const.h' 
      include 'masses.f' 
      include 'ewinput.f' 
cc******************************************
      pi=dacos(-1d0)!Expensive way of calculating pi?
      nf=5
c      nnf=nf !number of flavours (why duplicated?)
      CA=3d0
      Cf=4/3d0
      Euler=0.57721566d0
      Z2=1.644934d0
      Z3=1.202057d0
      b0=2*dexp(-Euler)
c      b0p=b0*a_param
      beta0=(33-2*nf)/12d0
      beta1=(153-19*nf)/24d0
      beta2=2857/128d0-5033*nf/1152d0+325*nf**2/3456d0
c......gluon coefficients
      A1g=CA
      A2g=CA/2d0*(67/6d0-(pi**2)/2d0-5d0/9*nf)
      A3g=CA*(13.81-2.15*nf-nf**2/108d0)
      B1g=-(11*CA-2*nf)/6d0
      B2g=CA**2*(23/24d0+(11*pi**2)/18d0-3*Z3/2d0)+
     /        Cf*nf/2d0-CA*nf*(1/12d0+pi**2/9d0)-11/8d0*Cf*CA 
      C1ggn=(pi**2/2d0+11/2d0+pi**2)/2d0
C
c.....quark coefficients
      A1q=Cf
      A2q=Cf/2d0*(67/6d0-(pi**2)/2d0-5/9d0*nf)
      A3q=Cf*(13.81-2.15*nf-nf**2/108d0)            ! 
     /+Cf*(CA*(29.9259d0-28d0*Z3)-8.2963d0*nf/2d0)  ! 
     /*2d0*(beta0*4d0)/64d0                         ! A3 from Becher & Neubert
      B1q=-(3d0*Cf)/2d0

      B2q=Cf**2*(pi**2/4d0-3/16d0-3*Z3)+
     /        CA*Cf*(11*pi**2/36d0-193/48d0+3*Z3/2d0)+
     /        Cf*nf*(17/24d0-pi**2/18d0)
       C1qqn=Cf/2d0*(pi**2/2d0-4d0)! Only delta(1-z) part, i.e. N independent part
c**************************************


c******************************************
c     Copy or set physics constants, inputs from mdata.f
      gf=Gf_inp !1.16639d-5
      Mz=zmass_inp !91.1876d0                        ! z mass ! 
      Mw=wmass_inp !80.399d0                         ! w mass !
      zw=zwidth !2.4952d0                         ! z width ! 
      ww=wwidth !2.085d0                          ! w width !
      gZ=dsqrt(dsqrt(2d0)*gf*Mz**2)
      gW=dsqrt(4d0*dsqrt(2d0)*gf*Mw**2)
c******************************************

c******************************************
c     Copy or set physics constants, inputs from mdata.f (can group with above, do only once in init)
C  ELECTROWEAK PARAMETER SCHEME: SAME SCHEME OF DYNNLO (G_mu scheme)
c      aem=1d0/128.89d0  
c      sw2=pi*aem/(dsqrt(2d0)*gf*Mw**2)
      sw2=1d0-(Mw/Mz)**2
      cw2=1d0-sw2
      aem=dsqrt(2d0)*gf*Mw**2*sw2/pi
      aweak=aem/sw2     
c******************************************

c******************************************
c     More settings of set physics constants, (should do only once in init)
      gLZu=gZ*(1/2d0-eequ*sw2)
      gLZd=gZ*(-1/2d0-eeqd*sw2)
      gRZu=-gZ*eequ*sw2
      gRZd=-gZ*eeqd*sw2
      fLZ=gZ*(-1/2d0+sw2)
      fRZ=gZ*sw2

      gLW=gW/dsqrt(2d0)
      fLW=gW/dsqrt(2d0)
!      gRW=0d0
!      fRW=0d0
c******************************************

      end
