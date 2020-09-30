C ADDED: Computes the various C2 coefficients at a given value of XN

      subroutine H2calc(C2qg,C2NSqqb,C2NSqq,C2Sqqb,XN)
      implicit complex*16 (A - Z) 
      REAL*8 zeta2,zeta3,zeta4,plog41f2,log2,pi,GE,CF,CA
      integer nf
      common/nf/nf

c       zeta2=   1.64493406685d0
c       zeta3=   1.20205690316d0
c     zeta4=   1.08232323371d0
      zeta2 = 1.64493406684822643647d0
      zeta3 = 1.20205690315959428540d0
      zeta4 = 1.08232323371113819152d0
       
       plog41f2=0.5174790616738993d0
       log2=0.693147180559945d0
       pi=3.141592653589793d0

       GE   = 0.57721566490153D0
       ONE=CMPLX(1.0D0,0.0D0)
       CF=4/3d0
       CA=3d0


       XNP1= XN+ONE
       XNP2= XN+ONE+ONE
       XNP3= XN+ONE+ONE+ONE
       XNM1= XN-ONE
       XNf2P1= XN/2d0+ONE
       XNf2P2= XN/2d0+ONE+ONE
       XNf2P3f2= XN/2d0+3/2d0*ONE
       XNf2P1f2= XN/2d0+1/2d0*ONE
       XNf2= XN/2d0


!     Function needed to write the H2 coefficient in the Mellin space

      CALL PSI0(XNP1,PS0NP1)            
      CALL PSI0(XNP2,PS0NP2)
      CALL PSI0(XNP3,PS0NP3)
      CALL PSI0(XN,PS0N)
      CALL PSI0(XNf2P1,PS0Nf2P1)
      CALL PSI0(XNf2P3f2,PS0Nf2P3f2)
      CALL PSI0(XNf2P1f2,PS0Nf2P1f2)
      CALL PSI0(XNf2,PS0Nf2)
      CALL PSI0(XNf2P2,PS0Nf2P2)
!
      CALL PSI1(XNP1,PS1NP1)
      CALL PSI1(XNP2,PS1NP2)
      CALL PSI1(XNP3,PS1NP3)
      CALL PSI1(XN,PS1N)
      CALL PSI1(XNf2,PS1Nf2)
      CALL PSI1(XNf2P1,PS1Nf2P1)
      CALL PSI1(XNf2P3f2,PS1Nf2P3f2)
      CALL PSI1(XNf2P1f2,PS1Nf2P1f2)
      CALL PSI1(XNf2P2,PS1Nf2P2)

!
      CALL PSI2(XNP1,PS2NP1)
      CALL PSI2(XNP2,PS2NP2)
      CALL PSI2(XNP3,PS2NP3)
      CALL PSI2(XN,PS2N)
      CALL PSI2(XNf2,PS2Nf2)
      CALL PSI2(XNf2P1,PS2Nf2P1)
      CALL PSI2(XNf2P3f2,PS2Nf2P3f2)
      CALL PSI2(XNf2P1f2,PS2Nf2P1f2)
      CALL PSI2(XNf2P2,PS2Nf2P2)
!
      CALL PSI3(XN,PS3N)
      CALL PSI3(XNf2,PS3Nf2)
      CALL PSI3(XNf2P1f2,PS3Nf2P1f2)
      CALL PSI3(XNf2,PS3Nf2)
      
      CALL BET(XN,BETA0N)
      CALL BET(XNP1,BETA0NP1)
      CALL BET(XNP2,BETA0NP2)
      CALL BET(XNP3,BETA0NP3)
      CALL BET1(XN,BETA1N)
      CALL BET2(XN,BETA2N)
      CALL BET3(XN,BETA3N)

      S1N= PS0NP1+GE        
      S1NP1= PS0NP2+GE      
      S1NP2= PS0NP3+GE
      S1NM1= PS0N+GE
      S1Nf2= PS0Nf2P1+GE
      S1NM1f2= PS0Nf2P1f2+GE
      S1NP1f2= PS0Nf2P3f2+GE
      S1NM2f2= PS0Nf2+GE
      S1NP2f2= PS0Nf2P2+GE

      S2N= -PS1NP1+zeta2
      S2NP1= -PS1NP2+zeta2  
      S2NP2= -PS1NP3+zeta2  
      S2NM1= -PS1N+zeta2    
      S2Nf2= -PS1Nf2P1+zeta2
      S2NP1f2= -PS1Nf2P3f2+zeta2
      S2NM1f2= -PS1Nf2P1f2+zeta2
      S2NP2f2= -PS1Nf2P2+zeta2

      S3N= PS2NP1/2d0+zeta3
      S3NP1= PS2NP2/2d0+zeta3
      S3NP2=  PS2NP3/2d0+zeta3
!      S2NP2= PS1NP3+zeta2
      S3NM1= PS2N/2d0+zeta3
      S3Nf2= PS2Nf2P1/2d0+zeta3
      S3NP1f2= PS2Nf2P3f2/2d0+zeta3
      S3NM1f2= PS2Nf2P1f2/2d0+zeta3
      S3NP2f2= PS2Nf2P2/2d0+zeta3


      S4NM1= -PS3N/6d0+zeta4
      S4NM1f2= -PS3Nf2P1f2/6d0+zeta4
      S4NM2f2= -PS3Nf2/6d0+zeta4
      S211NM1= -ACG21(XNM1)+6d0/5*zeta2**2
!      S112NM1=0
!      S121NM1=0

      S31NM1= ACG20(XNM1)+zeta2*S2NM1-1d0/2*zeta2**2


      Sm1NM1= (1)*BETA0N-log2
      Sm1N=   (-1)*BETA0NP1-log2
      Sm1NP1= (1)*BETA0NP2-log2
      Sm1NP2= (-1)*BETA0NP3-log2
      Sm2NM1= (-1)*BETA1N-zeta2/2d0
      Sm3NM1= (1)*BETA2N/2d0-zeta3*3/4d0
      Sm4NM1= (-1)*BETA3N/6d0-zeta4*7/8d0

      Sm21N= ACG3(XN)*(1)+zeta2*Sm1N
     /       -5/8d0*zeta3+zeta2*log2
      Sm21NP1= ACG3(XNP1)*(-1)+zeta2*Sm1NP1
     /       -5/8d0*zeta3+zeta2*log2
      Sm21NP2= ACG3(XNP2)*(1)+zeta2*Sm1NP2
     /       -5/8d0*zeta3+zeta2*log2

      Sm31NM1= ACG6(XNM1)*(1)+zeta2*Sm2NM1-zeta3*Sm1NM1
     /       -3/5d0*zeta2**2+2*plog41f2+3/4d0*zeta3*log2
     /       -zeta2/2d0*log2**2+log2**4/12d0
      Sm22NM1= ACG5(XNM1)*(1)-2d0*Sm31NM1
     /         +2d0*zeta2*Sm2NM1+3/40d0*zeta2**2
      Sm211NM1= ACG8(XNM1)*(-1)+zeta3*Sm1NM1-plog41f2+zeta2**2/8d0
     /     +zeta3*log2/8d0+zeta2*log2**2/4d0-log2**4/24d0


      ACG13XNM1 = ACG13(XNM1)
      ACG1XNM1  = ACG1(XNM1)
      ACG1pXN   = ACG1p(XN)
      ACG1pXNP1 = ACG1p(XNP1)
      ACG1ppXNM1= ACG1pp(XNM1)
      ACG2pXNM1 = ACG2p(XNM1)
      ACG4XN    = ACG4(XN)
      ACG4XNP1  = ACG4(XNP1)
      ACG4pXNM1 = ACG4p(XNM1)
      ACG5XNM1  = ACG5(XNM1)
      ACG6XNM1  = ACG6(XNM1)
      ACG7XNM1  = ACG7(XNM1)
      ACG9XNM1  = ACG9(XNM1)


!      C2 coefficient in the Mellin space

       C2Sqqb=
     / (CF*(172/XNM1-54/XN**4-27/XN**3
     / -126/XN**2-315/XN-54/XNP1**4- 27/XNP1**3+180/XNP1**2
     / +279/XNP1 -72/XNP2**3-192/XNP2**2
     / -136/XNP2-(72*zeta2)/XNM1 +(108*zeta2)/XN-
     / (108*zeta2)/XNP1+(72*zeta2)/XNP2
     / +(72*S2NM1)/XNM1 -(108*S2N)/XN
     / +(108*S2NP1)/XNP1 -(72*S2NP2)/XNP2))/432


      C2NSqq =
!
     / CF*nf*(
     / 127/192d0 -1/(24*XN**3)+5/(72*XN**2)-37/(216*XN) -1/(24*XNP1**3)
     /  +5/(72*XNP1**2)-19/(216*XNP1)-(2*zeta2)/3+(7*zeta3)/36
     /  -(7*S1NM1)/27+(5*S2NM1)/36-S3NM1/12)+
     / CA*CF*(
     / -1535/384d0-1/(8*XN**4)+11/(48*XN**3)-47/(72*XN**2)
     / +50/(27*XN)
     / - 1/(8*XNP1**4)-1/(48*XNP1**3)-83/(72*XNP1**2)+1/(54*XNP1)+
     / (97*zeta2)/24-(3*zeta2)/(8*XN)+(3*zeta2)/(8*XNP1)+(85*zeta3)/72-
     /(5*zeta3)/(8*XN)-(5*zeta3)/(8*XNP1)
     /-(29*zeta4)/16+(101*S1NM1)/54-
     /   (5*zeta3*S1NM1)/4+(zeta2*S1NM1**2)/4+(zeta2*S1N)/(4*XN)-
     /   S1NP1/(8*XNP1)+(zeta2*S1NP1)/(4*XNP1)-(19*S2NM1)/18
     /   +(zeta2*S2NM1)/4-(S1NM1**2*S2NM1)/4+S2N/(4*XN**2)+
     /   S2N/(4*XN)-(S1N*S2N)/(4*XN)+S2NP1/(4*XNP1**2)
     /   -S2NP1/(4*XNP1)-(S1NP1*S2NP1)/(4*XNP1)
     /  + S211NM1/2 +(11*S3NM1)/24-(S1NM1*S3NM1)/2
     /  - S3N/(4*XN)-S3NP1/(4*XNP1)-S4NM1/2)+
!
     / CF**2*(
     / 255/128d0-1/(8*XN**4)+3/(8*XN**2)-19/(8*XN)-1/(8*XNP1**4)+
     /   1/(4*XNP1**3)+2/XNP1**2+19/(8*XNP1)-(35*zeta2)/16
     /+(3*zeta2)/(2*XN)-(3*zeta2)/(2*XNP1)-(3*zeta3)/2-(3*zeta3)/(4*XN)-
     /  (3*zeta3)/(4*XNP1)+(61*zeta4)/16-(3*zeta3*S1NM1)/2-
     /  (zeta2*S1NM1**2)/2-S1N/(4*XN**2)
     /   -(zeta2*S1N)/(2*XN)+S1N**2/(8*XN**2)+
     /   S1NP1/(4*XNP1**2)+S1NP1/(8*XNP1)-
     /   (zeta2*S1NP1)/(2*XNP1)+S1NP1**2/(8*XNP1**2)+S2NM1
     /   -(zeta2*S2NM1)/2+(S1NM1**2*S2NM1)/2-(3*S2N)/(8*XN**2)-
     /   (3*S2N)/(4*XN)+(S1N*S2N)/(2*XN)-(3*S2NP1)/(8*XNP1**2)+
     /   (3*S2NP1)/(4*XNP1)+(S1NP1*S2NP1)/(2*XNP1)
     /  - S211NM1/2 -(3*S3NM1)/8+(3*S1NM1*S3NM1)/2+
     /   (3*S3N)/(4*XN)+(3*S3NP1)/(4*XNP1)-S31NM1/2+S4NM1)


      C2NSqqb = (-CA/2d0+CF)*CF/24d0*(
     /168d0+4d0*XN*(180d0+3*XN*(97d0+XN*(79d0+XN*(32d0+XN*(7d0+XN))))-
     /2*XN*XNP1**3*(1+2*XN)*pi**2)-XN**2*XNP1**2*
     /(192*XN**2*XNP1**2*ACG13XNM1+16*XN**2*XNP1**2*pi**2*ACG1XNM1+
     /48*(PS1Nf2-PS1N)+XN*(96*XNP1**2*ACG1pXN
     /-96*XN*XNP1*ACG1pXNP1-48*XN*ACG1ppXNM1-96*XN**2*ACG1ppXNM1
     /-48*XN**3*ACG1ppXNM1+96*XN*ACG2pXNM1+192*XN**2*ACG2pXNM1+
     /96*XN**3*ACG2pXNM1+96*ACG4XN+192*XN*ACG4XN
     /+96*XN**2*ACG4XN-96*XN*ACG4XNP1-96*XN**2*ACG4XNP1+
     /96*XN*ACG4pXNM1+192*XN**2*ACG4pXNM1+96*XN**3*ACG4pXNM1+
     /96*XN*ACG5XNM1+192*XN**2*ACG5XNM1+96*XN**3*ACG5XNM1-
     /192*XN*ACG6XNM1-384*XN**2*ACG6XNM1-192*XN**3*ACG6XNM1
     /-288*XN*ACG7XNM1-576*XN**2*ACG7XNM1-288*XN**3*ACG7XNM1
     /+192*XN*ACG9XNM1+384*XN**2*ACG9XNM1+
     /192*XN**3*ACG9XNM1-4*XN*XNP1*pi**2*PS0Nf2P3f2-
     /12*(-6+(-3+XN)*XN)*PS1Nf2+48*XN*XNP1*PS1N-
     /12*XN*(3+XN)*PS1Nf2P3f2-3*XNP1*(2+3*XN)*
     /PS2Nf2+24*XNP1**2*PS2N+
     /3*XN*XNP1*PS2Nf2P3f2-XN*XNP1**2*PS3Nf2+
     /8*XN*XNP1**2*PS3N+8*XNP1**2*log2*
     /(pi**2+6*XN*zeta3)-8*XNP1**2*PS0N*
     /(pi**2+6*XN*zeta3)
     /+4*XNP1*PS0Nf2
     /*((2+3*XN)*pi**2+12*XN*XNP1*zeta3)
     /))
     /)/(4d0*XN**4*XNP1**4)

       C2qg=
!
     /CF*(
     /1/(16*XN**4)+1/(32*XN**3)-1/(4*XN**2)-13/(32*XN)
     /-1/(8*XNP1**4)+
     /3/(8*XNP1**3)-15/(32*XNP1**2)+43/(32*XNP1)+1/(4*XNP2**4)-
     /1/(4*XNP2**3)+1/(4*XNP2**2)-5/(4*XNP2)+(3*zeta2)/(4*XNP1)-
     /(3*zeta2)/(4*XNP2)+zeta3/XN-(2*zeta3)/XNP1+(2*zeta3)/XNP2-
     /S1N**2/(16*XN**2)+S1N**3/(48*XN)-S1NP1/(4*XNP1**2)+
     /(3*S1NP1)/(16*XNP1)+S1NP1**2/(8*XNP1**2)
     /+S1NP1**2/(8*XNP1)-S1NP1**3/(24*XNP1)+
     /S1NP2/(4*XNP2**2)-S1NP2/(4*XNP2)-
     /S1NP2**2/(8*XNP2**2)-S1NP2**2/(8*XNP2)+
     /S1NP2**3/(24*XNP2)-S2N/(16*XN**2)+(S1N*S2N)/(16*XN)
     /+S2NP1/(8*XNP1**2)-S2NP1/(8*XNP1)-
     /(S1NP1*S2NP1)/(8*XNP1)-S2NP2/(8*XNP2**2)
     /+S2NP2/(8*XNP2)+(S1NP2*S2NP2)/(8*XNP2)-S3N/(12*XN)
     /+S3NP1/(6*XNP1)-S3NP2/(6*XNP2))+
!
     /CA*(
     /43/(108*XNM1)-1/(8*XN**4)-1/(16*XN**3)
     /-7/(24*XN**2)-35/(48*XN)-
     /1/(4*XNP1**4)+1/(4*XNP1**3)+5/(12*XNP1**2)+43/(48*XNP1)-
     /11/(12*XNP2**3)-17/(18*XNP2**2)-149/(216*XNP2)-
     /zeta2/(6*XNM1)+zeta2/(8*XN**2)
     /-((-1)*zeta2)/(8*XN**2)+zeta2/(4*XN)+
     /zeta2/(4*XNP1**2)+((-1)*zeta2)/(4*XNP1**2)-(7*zeta2)/(8*XNP1)-
     /((-1)*zeta2)/(8*XNP1)+zeta2/(4*XNP2**2)
     /-((-1)*zeta2)/(4*XNP2**2)+
     /(19*zeta2)/(24*XNP2)+((-1)*zeta2)/(8*XNP2)-(5*zeta3)/(16*XN)-
     /((-1)*zeta3)/(16*XN)+(5*zeta3)/(8*XNP1)
     /+((-1)*zeta3)/(8*XNP1)-
     /(5*zeta3)/(8*XNP2)-((-1)*zeta3)/(8*XNP2)
     /-(zeta2*log2)/(8*XN)+
     /((1)*zeta2*log2)/(8*XN)-(zeta2*log2)/(4*XNP1)+
     /((1)*zeta2*log2)/(4*XNP1)
     /-(zeta2*log2)/(4*XNP2)+
     /((1)*zeta2*log2)/(4*XNP2)+log2**2/(4*XN**2)-
     /((1)*log2**2)/(4*XN**2)+log2**2/(2*XNP1**2)-
     /((1)*log2**2)/(2*XNP1**2)+log2**2/(2*XNP2**2)-
     /((1)*log2**2)/(2*XNP2**2)
     /+(log2*S1NM1f2)/(4*XN**2)-
     /((1)*log2*S1NM1f2)/(4*XN**2)
     /-(log2*S1Nf2)/(4*XN**2)+
     /((1)*log2*S1Nf2)/(4*XN**2)
     /+(log2*S1Nf2)/(2*XNP1**2)-
     /((1)*log2*S1Nf2)/(2*XNP1**2)+(zeta2*S1N)/(8*XN)+
     /((-1)*zeta2*S1N)/(8*XN)-(S1NM1f2*S1N)/(4*XN**2)+
     /((1)*S1NM1f2*S1N)/(4*XN**2)+(S1Nf2*S1N)/(4*XN**2)-
     /((1)*S1Nf2*S1N)/(4*XN**2)-S1N**3/(48*XN)-
     /(log2*S1NP1f2)/(2*XNP1**2)+((1)*log2*S1NP1f2)/
     /(2*XNP1**2)+(log2*S1NP1f2)/(2*XNP2**2)-
     /((1)*log2*S1NP1f2)/(2*XNP2**2)-
     /(3*S1NP1)/(16*XNP1)-(zeta2*S1NP1)/(4*XNP1)-
     /((-1)*zeta2*S1NP1)/(4*XNP1)-(S1Nf2*S1NP1)/(2*XNP1**2)+
     /((1)*S1Nf2*S1NP1)/(2*XNP1**2)+
     /(S1NP1f2*S1NP1)/(2*XNP1**2)-
     /((1)*S1NP1f2*S1NP1)/(2*XNP1**2)-
     /S1NP1**2/(8*XNP1)+S1NP1**3/(24*XNP1)-
     /(log2*S1NP2f2)/(2*XNP2**2)+((1)*log2*S1NP2f2)/
     /(2*XNP2**2)+S1NP2/(4*XNP2)+(zeta2*S1NP2)/(4*XNP2)+
     /((-1)*zeta2*S1NP2)/(4*XNP2)-(S1NP1f2*S1NP2)/
     /(2*XNP2**2)+((1)*S1NP1f2*S1NP2)/(2*XNP2**2)+
     /(S1NP2f2*S1NP2)/(2*XNP2**2)-
     /((1)*S1NP2f2*S1NP2)/(2*XNP2**2)+
     /S1NP2**2/(8*XNP2)-S1NP2**3/(24*XNP2)-
     /S2NM1f2/(16*XN**2)-((-1)*S2NM1f2)/(16*XN**2)+
     /((1)*S2NM1f2)/(8*XN**2)+(log2*S2NM1f2)/(16*XN)-
     /((1)*log2*S2NM1f2)/(16*XN)-(S1N*S2NM1f2)/
     /(16*XN)+((-1)*S1N*S2NM1f2)/(16*XN)+S2NM1/(6*XNM1)+
     /S2Nf2/(16*XN**2)-((-1)*S2Nf2)/(16*XN**2)-
     /((1)*S2Nf2)/(8*XN**2)-S2Nf2/(8*XNP1**2)+
     /((-1)*S2Nf2)/(8*XNP1**2)+((1)*S2Nf2)/(4*XNP1**2)-
     /((-1)*S2Nf2)/(16*XNP1)-((1)*S2Nf2)/(16*XNP1)-
     /(log2*S2Nf2)/(16*XN)+((1)*log2*S2Nf2)/(16*XN)+
     /(log2*S2Nf2)/(8*XNP1)
     /-((1)*log2*S2Nf2)/(8*XNP1)+
     /(S1N*S2Nf2)/(16*XN)+((-1)*S1N*S2Nf2)/(16*XN)-
     /(S1NP1*S2Nf2)/(8*XNP1)-((-1)*S1NP1*S2Nf2)/(8*XNP1)+
     /((-1)*S2N)/(4*XN**2)-S2N/(4*XN)-(S1NM1f2*S2N)/(8*XN)+
     /((1)*S1NM1f2*S2N)/(8*XN)+(S1Nf2*S2N)/(8*XN)-
     /((1)*S1Nf2*S2N)/(8*XN)-(3*S1N*S2N)/(16*XN)-
     /((-1)*S1N*S2N)/(4*XN)+S2NP1f2/(8*XNP1**2)+
     /((-1)*S2NP1f2)/(8*XNP1**2)-((1)*S2NP1f2)/
     /(4*XNP1**2)-((-1)*S2NP1f2)/(16*XNP1)+
     /((1)*S2NP1f2)/(16*XNP1)-S2NP1f2/(8*XNP2**2)-
     /((-1)*S2NP1f2)/(8*XNP2**2)+((1)*S2NP1f2)/
     /(4*XNP2**2)+((-1)*S2NP1f2)/(16*XNP2)-
     /((1)*S2NP1f2)/(16*XNP2)-(log2*S2NP1f2)/
     /(8*XNP1)+((1)*log2*S2NP1f2)/(8*XNP1)+
     /(log2*S2NP1f2)/(8*XNP2)-((1)*log2*S2NP1f2)/
     /(8*XNP2)+(S1NP1*S2NP1f2)/(8*XNP1)
     /-((-1)*S1NP1*S2NP1f2)/(8*XNP1)-
     /(S1NP2*S2NP1f2)/(8*XNP2)+((-1)*S1NP2*S2NP1f2)/
     /(8*XNP2)-S2NP1/(2*XNP1**2)-((-1)*S2NP1)/
     /(2*XNP1**2)+(7*S2NP1)/(8*XNP1)+
     /((-1)*S2NP1)/(4*XNP1)-(S1Nf2*S2NP1)/(4*XNP1)+
     /((1)*S1Nf2*S2NP1)/(4*XNP1)+(S1NP1f2*S2NP1)/
     /(4*XNP1)-((1)*S1NP1f2*S2NP1)/(4*XNP1)+
     /(3*S1NP1*S2NP1)/(8*XNP1)+((-1)*S1NP1*S2NP1)/
     /(2*XNP1)+S2NP2f2/(8*XNP2**2)-((-1)*S2NP2f2)/
     /(8*XNP2**2)-((1)*S2NP2f2)/(4*XNP2**2)
     /+((-1)*S2NP2f2)/(16*XNP2)+((1)*S2NP2f2)/
     /(16*XNP2)-(log2*S2NP2f2)/(8*XNP2)+
     /((1)*log2*S2NP2f2)/(8*XNP2)+
     /(S1NP2*S2NP2f2)/(8*XNP2)+((-1)*S1NP2*S2NP2f2)/
     /(8*XNP2)+((-1)*S2NP2)/(2*XNP2**2)
     /-(19*S2NP2)/(24*XNP2)-((-1)*S2NP2)/(4*XNP2)-
     /(S1NP1f2*S2NP2)/(4*XNP2)+
     /((1)*S1NP1f2*S2NP2)/(4*XNP2)+
     /(S1NP2f2*S2NP2)/(4*XNP2)
     /-((1)*S1NP2f2*S2NP2)/(4*XNP2)-
     /(3*S1NP2*S2NP2)/(8*XNP2)-((-1)*S1NP2*S2NP2)/
     /(2*XNP2)
     /-(3*S3NM1f2)/(64*XN)+((-1)*S3NM1f2)/
     /(64*XN)
     /+((1)*S3NM1f2)/(32*XN)+(3*S3Nf2)/(64*XN)
     /+((-1)*S3Nf2)/(64*XN)-((1)*S3Nf2)/(32*XN)
     /-(3*S3Nf2)/(32*XNP1)-((-1)*S3Nf2)/(32*XNP1)
     /+((1)*S3Nf2)/(16*XNP1)-S3N/(6*XN)-((-1)*S3N)/(8*XN)
     /+(3*S3NP1f2)/(32*XNP1)-((-1)*S3NP1f2)/(32*XNP1)-
     /((1)*S3NP1f2)/(16*XNP1)-(3*S3NP1f2)/(32*XNP2)+
     /((-1)*S3NP1f2)/(32*XNP2)+((1)*S3NP1f2)/
     /(16*XNP2)+S3NP1/(3*XNP1)+((-1)*S3NP1)/(4*XNP1)
     /+(3*S3NP2f2)/(32*XNP2)+((-1)*S3NP2f2)/(32*XNP2)-
     /((1)*S3NP2f2)/(16*XNP2)-S3NP2/(3*XNP2)-
     /((-1)*S3NP2)/(4*XNP2)-((-1)*Sm21N)/(4*XN)+
     /((-1)*Sm21NP1)/(2*XNP1)-((-1)*Sm21NP2)/(2*XNP2))



! NOTE: C2 coefficients have aS/Pi normalizations. Change to aS/(2*Pi)

       C2qg=4d0*C2qg
       C2NSqqb=4d0*C2NSqqb
       C2Sqqb=4d0*C2Sqqb
       C2NSqq=4d0*C2NSqq

      return
      end
      

C
C    BLUMLEIN FUNCTIONS
C
      SUBROUTINE PSI0(ZZ,RES)
C     -----------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  PSI(Z) FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 Z,T,ONE,Y,T0,RES,ZZ
      COMPLEX*16 Y2 !!!bug fix
C
      ONE=DCMPLX(1.0D0,0.0D0)
C
      Z=ZZ
      T=DCMPLX(0.0D0,0.0D0)
2     R=DSQRT(DREAL(Z)**2+DIMAG(Z)**2)
      IF(R.GT.10.0D0) GOTO 1
      T=T-ONE/Z
      Z=Z+ONE
      GOTO 2
1     Y=ONE/Z
      Y2=Y*Y
      T0 = (-1.D0/2.D0+(-1.D0/12.D0+(1.D0/120.D0+(-1.D0/252.D0+(1.D0/240
     #.D0+(-1.D0/132.D0+(691.D0/32760.D0+(-ONE/12.0D0+ONE*3617.0D0
     #             /8160.D0*Y2
     #  )*Y2  )*Y2  )*Y2  )*Y2  )*Y2  )*Y2  )*Y)*Y-LOG(Y)
C
      RES=T+T0
C
      RETURN
      END
      SUBROUTINE PSI1(ZZ,RES)
C     -----------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  PSI'(Z) FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 Z,T,ONE,Y,T0,RES,ZZ
      COMPLEX*16 Y2 !!!bug fix
C
      ONE=DCMPLX(1.0D0,0.0D0)
C
      Z=ZZ
      T=DCMPLX(0.0D0,0.0D0)
2     R=DSQRT(DREAL(Z)**2+DIMAG(Z)**2)
      IF(R.GT.10.0D0) GOTO 1
      T=T+ONE/Z**2
      Z=Z+ONE
      GOTO 2
1     Y=ONE/Z
      Y2=Y*Y
      T0 = (1.D0+(1.D0/2.D0+(1.D0/6.D0+(-1.D0/30.D0+
     &(1.D0/42.D0+(-1.D0/30.D0+(5.D0/66.D0-691.D0/2730.D0*Y2)
     &*Y2)*Y2)*Y2)*Y2)*Y)*Y)*Y
C
      RES=T+T0
C
      RETURN
      END
      SUBROUTINE PSI2(ZZ,RES)
C     -----------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  PSI''(Z) FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 Z,T,ONE,Y,T0,RES,ZZ
      COMPLEX*16 Y2 !!!bug fix
C
      ONE=DCMPLX(1.0D0,0.0D0)
      TWO=ONE*2.0D0
C
      Z=ZZ
      T=DCMPLX(0.0D0,0.0D0)
2     R=DSQRT(DREAL(Z)**2+DIMAG(Z)**2)
      IF(R.GT.10.0D0) GOTO 1
      T=T-TWO/Z**3
      Z=Z+ONE
      GOTO 2
1     Y=ONE/Z
      Y2=Y*Y
      T0 =(-1.D0+(-1.D0+(-1.D0/2.D0+(1.D0/6.D0+(-1.D0/6.D0+(3.D0/
     &10.D0+(-5.D0/6.D0+691.D0/210.D0*Y2)*Y2)*Y2)*Y2)*Y2)*Y)*Y)*Y2
C
      RES=T+T0
C
      RETURN
      END
      SUBROUTINE PSI3(ZZ,RES)
C     -----------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  PSI'''(Z) FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 Z,T,ONE,Y,T0,RES,ZZ
      COMPLEX*16 Y2 !!!bug fix
C
      ONE=DCMPLX(1.0D0,0.0D0)
      SIX=ONE*6.0D0
C
      Z=ZZ
      T=DCMPLX(0.0D0,0.0D0)
2     R=DSQRT(DREAL(Z)**2+DIMAG(Z)**2)
      IF(R.GT.10.0D0) GOTO 1
      T=T+SIX/Z**4
      Z=Z+ONE
      GOTO 2
1     Y=ONE/Z
      Y2=Y*Y
C
      T0 = (2.D0+(3.D0+(2.D0+(-1.D0+(4.D0/3.D0+(-3.D0+(10.D0+(-691.D0/15
     #.D0+(280.D0-10851.D0/5.D0*Y2  )*Y2  )*Y2  )*Y2  )*Y2  )*Y2  )*Y2
     #)*Y)*Y)*Y2*Y
C
      RES=T+T0
C
      RETURN
      END

      SUBROUTINE BET(ZZ,RES)
C     ----------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  BET(Z) FOR COMPLEX ARGUMENT  $\beta(z)$
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 Z,ONE,RES,V1,V2,Z1,Z2,ZZ
C
      ONE=DCMPLX(1.0D0,0.0D0)
C
      Z=ZZ
      Z1=(Z+ONE)/2.0D0
      Z2=Z/2.0D0
      CALL PSI0(Z1,V1)
      CALL PSI0(Z2,V2)
C
      RES=(V1-V2)/2.0D0
C
      RETURN
      END
      SUBROUTINE BET1(ZZ,RES)
C     -----------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  BET'(Z) FOR COMPLEX ARGUMENT  $\beta'(z)$
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 Z,ONE,RES,V1,V2,Z1,Z2,ZZ
C
      ONE=DCMPLX(1.0D0,0.0D0)
C
      Z=ZZ
      Z1=(Z+ONE)/2.0D0
      Z2=Z/2.0D0
      CALL PSI1(Z1,V1)
      CALL PSI1(Z2,V2)
C
      RES=(V1-V2)/4.0D0
C
      RETURN
      END
      SUBROUTINE BET2(ZZ,RES)
C     -----------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  BET''(Z) FOR COMPLEX ARGUMENT  $\beta''(z)$
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 Z,ONE,RES,V1,V2,Z1,Z2,ZZ
C
      ONE=DCMPLX(1.0D0,0.0D0)
C
      Z=ZZ
      Z1=(Z+ONE)/2.0D0
      Z2=Z/2.0D0
      CALL PSI2(Z1,V1)
      CALL PSI2(Z2,V2)
C
      RES=(V1-V2)/8.0D0
C
      RETURN
      END
      SUBROUTINE BET3(ZZ,RES)
C     -----------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  BET'''(Z) FOR COMPLEX ARGUMENT  $\beta'''(z)$
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 Z,ONE,RES,V1,V2,Z1,Z2,ZZ
C
      ONE=DCMPLX(1.0D0,0.0D0)
C
      Z=ZZ
      Z1=(Z+ONE)/2.0D0
      Z2=Z/2.0D0
      CALL PSI3(Z1,V1)
      CALL PSI3(Z2,V2)
C
!!!!!!!!!  TYPO IN BLUMLEIN!!!!!!!!
!!!!!!!!!      RES=(V1-V2)/8.0D0
      RES=(V1-V2)/16.0D0
C
      RETURN
      END
      COMPLEX*16 FUNCTION ACG1(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(1+X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 AK2,DL,ZERO,ONE
C
!      COMMON/ACLOG2/ AK2(10)
!      COMMON/ACCON1/ DL
!      COMMON/VAL   / ZERO,ONE
!      DIMENSION AK2(10)
!      DATA AK2/0.999999980543793D+0,
!     &        -0.999995797779624D+0,
!     &         0.916516447393493D+0,
!     &        -0.831229921350708D+0,
!     &         0.745873737923571D+0,
!     &        -0.634523908078600D+0,
!     &         0.467104011423750D+0,
!     &        -0.261348046799178D+0,
!     &         0.936814286867420D-1,
!     &        -0.156249375012462D-1/

      DIMENSION AK2(23)
      DATA AK2/1.0000000000000000D-0,
     & -0.9999999999999985D-0,
     &  0.9166666666663948D-0,
     & -0.8333333333136118D-0,
     &  0.7611111103508889D-0,
     & -0.6999999819735105D-0,
     &  0.6482139985629993D-0,
     & -0.6039649964806160D-0,
     &  0.5657662306410356D-0,
     & -0.5323631718571445D-0,
     &  0.5024238774786239D-0,
     & -0.4738508288315496D-0,
     &  0.4427472719775835D-0,
     & -0.4029142806330511D-0,
     &  0.3476841543351489D-0,
     & -0.2748590021353420D-0,
     &  0.1915627642585285D-0,
     & -0.1130763066428224D-0,
     &  0.5415661067306229D-1,
     & -0.1999877298940919D-1,
     &  0.5303624439388411D-2,
     & -0.8944156375768203D-3,
     &  0.7179502917974332D-4/
      
      DL = DLOG(2.0D0)
      ZERO=0.0D0
      ONE =1.0D0
C
      T=DCMPLX(ZERO,ZERO)
      DO 1 L=2,24!11
      K=L-1
      T=T+AK2(K)/(ZN+DBLE(K+1))
1     CONTINUE
C
      ACG1=(DL*DL- ZN*T)/2.0D0
C
      RETURN
      END

      COMPLEX*16 FUNCTION ACG1p(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(1+X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 AK2,DL,ZERO,ONE
C
!      COMMON/ACLOG2/ AK2(10)
!      COMMON/ACCON1/ DL
!      COMMON/VAL   / ZERO,ONE
!      DIMENSION AK2(10)
!      DATA AK2/0.999999980543793D+0,
!     &        -0.999995797779624D+0,
!     &         0.916516447393493D+0,
!     &        -0.831229921350708D+0,
!     &         0.745873737923571D+0,
!     &        -0.634523908078600D+0,
!     &         0.467104011423750D+0,
!     &        -0.261348046799178D+0,
!     &         0.936814286867420D-1,
!     &        -0.156249375012462D-1/
      DIMENSION AK2(23)
      DATA AK2/1.0000000000000000D-0,
     & -0.9999999999999985D-0,
     &  0.9166666666663948D-0,
     & -0.8333333333136118D-0,
     &  0.7611111103508889D-0,
     & -0.6999999819735105D-0,
     &  0.6482139985629993D-0,
     & -0.6039649964806160D-0,
     &  0.5657662306410356D-0,
     & -0.5323631718571445D-0,
     &  0.5024238774786239D-0,
     & -0.4738508288315496D-0,
     &  0.4427472719775835D-0,
     & -0.4029142806330511D-0,
     &  0.3476841543351489D-0,
     & -0.2748590021353420D-0,
     &  0.1915627642585285D-0,
     & -0.1130763066428224D-0,
     &  0.5415661067306229D-1,
     & -0.1999877298940919D-1,
     &  0.5303624439388411D-2,
     & -0.8944156375768203D-3,
     &  0.7179502917974332D-4/

      DL = DLOG(2.0D0)
      ZERO=0.0D0
      ONE =1.0D0
C
      T=DCMPLX(ZERO,ZERO)
      T1=DCMPLX(ZERO,ZERO)
      DO 1 L=2,24!11
      K=L-1
      T=T+AK2(K)/(ZN+DBLE(K+1))
      T1=T1+AK2(K)/(ZN+DBLE(K+1))/(ZN+DBLE(K+1))
1     CONTINUE
C
      ACG1p=(- T + ZN*T1)/2.0D0
C
      RETURN
      END
      COMPLEX*16 FUNCTION ACG1pp(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(1+X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 AK2,DL,ZERO,ONE
C
!      COMMON/ACLOG2/ AK2(10)
!      COMMON/ACCON1/ DL
!      COMMON/VAL   / ZERO,ONE
!      DIMENSION AK2(10)
!      DATA AK2/0.999999980543793D+0,
!     &        -0.999995797779624D+0,
!     &         0.916516447393493D+0,
!     &        -0.831229921350708D+0,
!     &         0.745873737923571D+0,
!     &        -0.634523908078600D+0,
!     &         0.467104011423750D+0,
!     &        -0.261348046799178D+0,
!     &         0.936814286867420D-1,
!     &        -0.156249375012462D-1/
      DIMENSION AK2(23)
      DATA AK2/1.0000000000000000D-0,
     & -0.9999999999999985D-0,
     &  0.9166666666663948D-0,
     & -0.8333333333136118D-0,
     &  0.7611111103508889D-0,
     & -0.6999999819735105D-0,
     &  0.6482139985629993D-0,
     & -0.6039649964806160D-0,
     &  0.5657662306410356D-0,
     & -0.5323631718571445D-0,
     &  0.5024238774786239D-0,
     & -0.4738508288315496D-0,
     &  0.4427472719775835D-0,
     & -0.4029142806330511D-0,
     &  0.3476841543351489D-0,
     & -0.2748590021353420D-0,
     &  0.1915627642585285D-0,
     & -0.1130763066428224D-0,
     &  0.5415661067306229D-1,
     & -0.1999877298940919D-1,
     &  0.5303624439388411D-2,
     & -0.8944156375768203D-3,
     &  0.7179502917974332D-4/

      DL = DLOG(2.0D0)
      ZERO=0.0D0
      ONE =1.0D0
C
      T=DCMPLX(ZERO,ZERO)
      T1=DCMPLX(ZERO,ZERO)
      DO 1 L=2,24!11
      K=L-1
      T=T+AK2(K)/(ZN+DBLE(K+1))/(ZN+DBLE(K+1))
      T1=T1+AK2(K)/(ZN+DBLE(K+1))/(ZN+DBLE(K+1))/(ZN+DBLE(K+1))
1     CONTINUE
C
      ACG1pp=( T - ZN*T1)
C
      RETURN
      END

      COMPLEX*16 FUNCTION ACG2(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(1+X)**2/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 AK3,DL,ZERO,ONE
C
!      COMMON/ACLOG3/ AK3(11)
!      COMMON/ACCON1/ DL
!      COMMON/VAL   / ZERO,ONE
      DIMENSION AK3(11)
      DATA AK3/9.99999989322696D-1,
     &        -1.49999722020708D+0,
     &         1.74988008499745D+0,
     &        -1.87296689068405D+0,
     &         1.91539974617231D+0,
     &        -1.85963744001295D+0,
     &         1.62987195424434D+0,
     &        -1.17982353224299D+0,
     &         6.28710122994999D-1,
     &        -2.11307487211713D-1,
     &         3.28953352932140D-2/
      DL = DLOG(2.0D0)
      ZERO=0.0D0
      ONE =1.0D0
C
      T=DCMPLX(ZERO,ZERO)
      DO 1 L=3,13
      K=L-2
      T=T+AK3(K)/(ZN+DBLE(K+2))
1     CONTINUE
C
      ACG2=(DL**3- ZN*T)/3.0D0
C
      RETURN
      END
      COMPLEX*16 FUNCTION ACG2p(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(1+X)**2/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 AK3,DL,ZERO,ONE
C
!      COMMON/ACLOG3/ AK3(11)
!      COMMON/ACCON1/ DL
!      COMMON/VAL   / ZERO,ONE
      DIMENSION AK3(11)
      DATA AK3/9.99999989322696D-1,
     &        -1.49999722020708D+0,
     &         1.74988008499745D+0,
     &        -1.87296689068405D+0,
     &         1.91539974617231D+0,
     &        -1.85963744001295D+0,
     &         1.62987195424434D+0,
     &        -1.17982353224299D+0,
     &         6.28710122994999D-1,
     &        -2.11307487211713D-1,
     &         3.28953352932140D-2/
      DL = DLOG(2.0D0)
      ZERO=0.0D0
      ONE =1.0D0
C
      T=DCMPLX(ZERO,ZERO)
      T1=DCMPLX(ZERO,ZERO)
      DO 1 L=3,13
      K=L-2
      T=T+AK3(K)/(ZN+DBLE(K+2))
      T1=T1+AK3(K)/(ZN+DBLE(K+2))/(ZN+DBLE(K+2))
1     CONTINUE
C
      ACG2p=(- T + ZN*T1)/3.0D0
C
      RETURN
      END

      COMPLEX*16 FUNCTION ACG3(ZN1)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LI2(X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 AK1,DL,ZERO,ONE,ZET2,GE
C
!      COMMON/ACLOG1/ AK1(9)
!      COMMON/IAPP  / IAPP
!      COMMON/VAL   / ZERO,ONE
!      COMMON/ACCON1/ DL
!      COMMON/ACCON2/ ZET2
!      COMMON/ACCON4/ GE
C
c     DIMENSION AK1(9)
      DIMENSION AK1(20)
C
c      ZET2=1.64493406685d0
      ZET2 = 1.64493406684822643647d0
      ZERO=0.0D0
      ONE =1.0D0
      GE   = 0.57721566490153D+0
      DL = DLOG(2.0D0)
      IAPP=1 !3
C
c      DATA AK1/0.999999974532240D+0,
c     &        -0.499995525890027D+0,
c     &         0.333203435554182D+0,
c     &        -0.248529457735332D+0,
c     &         0.191451164493502D+0,
c     &        -0.137466222203386D+0,
c     &         0.792107405737825D-1,
c     &        -0.301109652783781D-1,
c     &     0.538406198111749D-2/



      DATA AK1/0.9999999999999925D-0,
     &     -0.4999999999988568D-0,
     &      0.3333333332641123D-0,
     &     -0.2499999977763199D-0,
     &      0.1999999561535526D-0,
     &     -0.1666660875051348D-0,
     &      0.1428517138099479D-0,
     &     -0.1249623936313475D-0,
     &      0.1109128496887138D-0,
     &     -0.9918652787800788D-1,
     &      0.8826572954250856D-1,
     &     -0.7643209265133132D-1,
     &      0.6225829212455825D-1,
     &     -0.4572477090315515D-1,
     &      0.2890194939889559D-1,
     &     -0.1496621145891488D-1,
     &      0.6003156359511387D-2,
     &     -0.1731328252868496D-2,
     &      0.3172112728405899D-3,
     &     -0.2760099875146713D-4/
      
      
C
!
!
      ZN=ZN1+ONE
C
      IF(IAPP.EQ.1) GOTO 10
      IF(IAPP.EQ.2) GOTO 20
      IF(IAPP.EQ.3) GOTO 30
      WRITE(6,*) '*** ERROR IN  ACG3, IAPP=',IAPP,' WRONG, STOP ***'
      STOP
C
10    T=DCMPLX(DL*ZET2,ZERO)
      Z=ZN1
      Z1=Z+ONE
      DO 1 L=1,20 !9
      ZL =DCMPLX(DBLE(L),ZERO)
      ZL1=Z+ZL+ONE
      CALL PSI0(ZL1,PS)
      S1=PS+GE
      T=T-AK1(L)*(ZET2*Z/(Z+ZL)+ZL/(Z+ZL)**2*S1)
1     CONTINUE
      GOTO 100
C
20    T=1.01/(ZN+ONE)-0.846/(ZN+ONE*2)+1.155/(ZN+ONE*3)
     &   -1.074/(ZN+ONE*4)+0.55/(ZN+ONE*5)
      GOTO 100
30    T=1.004/(ZN+ONE)-0.846/(ZN+ONE*2)+1.342/(ZN+ONE*3)
     &   -1.532/(ZN+ONE*4)+0.839/(ZN+ONE*5)
100   CONTINUE
      ACG3=T
C
      RETURN
      END
      COMPLEX*16 FUNCTION ACG4(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LI2(-X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 AK1,DL,ZERO,ONE,ZET2
C
!      COMMON/ACLOG1/ AK1(9)
!      COMMON/ACCON1/ DL
!      COMMON/ACCON2/ ZET2
!      COMMON/VAL   / ZERO,ONE
c     DIMENSION AK1(9)
      DIMENSION AK1(20)
C
c      ZET2=1.64493406685d0
      ZET2 = 1.64493406684822643647d0
      ZERO=0.0D0
      ONE =1.0D0
      DL = DLOG(2.0D0)
C
!      DATA AK1/0.999999974532240D+0,
!     &        -0.499995525890027D+0,
!     &         0.333203435554182D+0,
!     &        -0.248529457735332D+0,
!     &         0.191451164493502D+0,
!     &        -0.137466222203386D+0,
!     &         0.792107405737825D-1,
!     &        -0.301109652783781D-1,
!     &         0.538406198111749D-2/
C

      DATA AK1/0.9999999999999925D-0,
     &     -0.4999999999988568D-0,
     &      0.3333333332641123D-0,
     &     -0.2499999977763199D-0,
     &      0.1999999561535526D-0,
     &     -0.1666660875051348D-0,
     &      0.1428517138099479D-0,
     &     -0.1249623936313475D-0,
     &      0.1109128496887138D-0,
     &     -0.9918652787800788D-1,
     &      0.8826572954250856D-1,
     &     -0.7643209265133132D-1,
     &      0.6225829212455825D-1,
     &     -0.4572477090315515D-1,
     &      0.2890194939889559D-1,
     &     -0.1496621145891488D-1,
     &      0.6003156359511387D-2,
     &     -0.1731328252868496D-2,
     &      0.3172112728405899D-3,
     &     -0.2760099875146713D-4/
      
      T=DCMPLX(-ZET2/2.0D0*DL,ZERO)
C
      DO 1 L=1,20 !9
      ZL =DCMPLX(DBLE(L),ZERO)
      ZNL =ZN+ZL
      ZNL1=ZNL+DCMPLX(ONE,ZERO)
      CALL BET(ZNL1,V1)
C
      T=T+AK1(L)*(ZN/ZNL*ZET2/2.0D0+ZL/ZNL**2*(DL-V1))
1     CONTINUE
C
      ACG4=T
C
      RETURN
      END
      COMPLEX*16 FUNCTION ACG4p(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LI2(-X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 AK1,DL,ZERO,ONE,ZET2
C
!      COMMON/ACLOG1/ AK1(9)
!      COMMON/ACCON1/ DL
!      COMMON/ACCON2/ ZET2
!      COMMON/VAL   / ZERO,ONE
c     DIMENSION AK1(9)
      DIMENSION AK1(20)
C
c      ZET2=1.64493406685d0
      ZET2 = 1.64493406684822643647d0
      ZERO=0.0D0
      ONE =1.0D0
      DL = DLOG(2.0D0)
C
!      DATA AK1/0.999999974532240D+0,
!     &        -0.499995525890027D+0,
!     &         0.333203435554182D+0,
!     &        -0.248529457735332D+0,
!     &         0.191451164493502D+0,
!     &        -0.137466222203386D+0,
!     &         0.792107405737825D-1,
!     &        -0.301109652783781D-1,
!     &         0.538406198111749D-2/
C
      DATA AK1/0.9999999999999925D-0,
     &     -0.4999999999988568D-0,
     &      0.3333333332641123D-0,
     &     -0.2499999977763199D-0,
     &      0.1999999561535526D-0,
     &     -0.1666660875051348D-0,
     &      0.1428517138099479D-0,
     &     -0.1249623936313475D-0,
     &      0.1109128496887138D-0,
     &     -0.9918652787800788D-1,
     &      0.8826572954250856D-1,
     &     -0.7643209265133132D-1,
     &      0.6225829212455825D-1,
     &     -0.4572477090315515D-1,
     &      0.2890194939889559D-1,
     &     -0.1496621145891488D-1,
     &      0.6003156359511387D-2,
     &     -0.1731328252868496D-2,
     &      0.3172112728405899D-3,
     &     -0.2760099875146713D-4/
      
      
      T=DCMPLX(ZERO,ZERO)
C
      DO 1 L=1,20!9
      ZL =DCMPLX(DBLE(L),ZERO)
      ZNL =ZN+ZL
      ZNL1=ZNL+DCMPLX(ONE,ZERO)
      CALL BET(ZNL1,V1)
      CALL BET1(ZNL1,V2)
C
      T=T+AK1(L)*((ZNL-ZN)/ZNL**2*ZET2/2.0D0-2.0D0*ZL/ZNL**3*(DL-V1)
     /-ZL/ZNL**2*V2)
1     CONTINUE
C
      ACG4p=T
C
      RETURN
      END

      COMPLEX*16 FUNCTION ACG5(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(X) LI2(X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 AK1,ZERO,ONE,DL,ZET2,GE
C
!      COMMON/ACLOG1/ AK1(9)
!      COMMON/VAL   / ZERO,ONE
!      COMMON/ACCON1/ DL
!      COMMON/ACCON2/ ZET2
!      COMMON/ACCON4/ GE
C
!     DIMENSION AK1(9)
      DIMENSION AK1(20)
C
c      ZET2=1.64493406685d0
      ZET2 = 1.64493406684822643647d0
      ZET3=1.20205690315959428540D+0
      ZERO=0.0D0
      ONE =1.0D0
      GE   = 0.57721566490153D+0
      DL = DLOG(2.0D0)
C
!      DATA AK1/0.999999974532240D+0,
!     &        -0.499995525890027D+0,
!     &         0.333203435554182D+0,
!     &        -0.248529457735332D+0,
!     &         0.191451164493502D+0,
!     &        -0.137466222203386D+0,
!     &         0.792107405737825D-1,
!     &        -0.301109652783781D-1,
!     &         0.538406198111749D-2/
C
      DATA AK1/0.9999999999999925D-0,
     &     -0.4999999999988568D-0,
     &      0.3333333332641123D-0,
     &     -0.2499999977763199D-0,
     &      0.1999999561535526D-0,
     &     -0.1666660875051348D-0,
     &      0.1428517138099479D-0,
     &     -0.1249623936313475D-0,
     &      0.1109128496887138D-0,
     &     -0.9918652787800788D-1,
     &      0.8826572954250856D-1,
     &     -0.7643209265133132D-1,
     &      0.6225829212455825D-1,
     &     -0.4572477090315515D-1,
     &      0.2890194939889559D-1,
     &     -0.1496621145891488D-1,
     &      0.6003156359511387D-2,
     &     -0.1731328252868496D-2,
     &      0.3172112728405899D-3,
     &     -0.2760099875146713D-4/
      
!
!
      T=DCMPLX(ZERO,ZERO)
      Z=ZN
      DO 1 L=1,20 !9
      ZL =DCMPLX(DBLE(L),ZERO)
      ZNL=Z+ZL
      ZL1=ZNL+ONE
      CALL PSI0(ZL1,PS)
      CALL PSI1(ZL1,S1P)
      S1=PS+GE
      T=T-AK1(L)*ZL/ZNL**2*(ZET2+S1P-2.0D0*S1/ZNL)
1     CONTINUE
C
      ACG5=T
C
      RETURN
      END
      COMPLEX*16 FUNCTION ACG6(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LI3(X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 AK1,DL,ZERO,ONE,ZET3,GE,ZET2
C
!      COMMON/ACLOG1/ AK1(9)
!      COMMON/ACCON1/ DL
!      COMMON/ACCON2/ ZET2
!      COMMON/ACCON3/ ZET3
!      COMMON/ACCON4/ GE
!      COMMON/VAL   / ZERO,ONE
C
c     DIMENSION AK1(9)
      DIMENSION AK1(20)
C
c      ZET2=1.64493406685d0
      ZET2 = 1.64493406684822643647d0
      ZET3=1.20205690315959428540D+0
      ZERO=0.0D0
      ONE =1.0D0
      GE   = 0.57721566490153D+0
      DL = DLOG(2.0D0)
C
!      DATA AK1/0.999999974532240D+0,
!     &        -0.499995525890027D+0,
!     &         0.333203435554182D+0,
!     &        -0.248529457735332D+0,
!     &         0.191451164493502D+0,
!     &        -0.137466222203386D+0,
!     &         0.792107405737825D-1,
!     &        -0.301109652783781D-1,
!     &         0.538406198111749D-2/
C
      DATA AK1/0.9999999999999925D-0,
     &     -0.4999999999988568D-0,
     &      0.3333333332641123D-0,
     &     -0.2499999977763199D-0,
     &      0.1999999561535526D-0,
     &     -0.1666660875051348D-0,
     &      0.1428517138099479D-0,
     &     -0.1249623936313475D-0,
     &      0.1109128496887138D-0,
     &     -0.9918652787800788D-1,
     &      0.8826572954250856D-1,
     &     -0.7643209265133132D-1,
     &      0.6225829212455825D-1,
     &     -0.4572477090315515D-1,
     &      0.2890194939889559D-1,
     &     -0.1496621145891488D-1,
     &      0.6003156359511387D-2,
     &     -0.1731328252868496D-2,
     &      0.3172112728405899D-3,
     &     -0.2760099875146713D-4/
      
!
!
      T=DCMPLX(DL*ZET3,ZERO)
      DO 1 L=1,20 !9
      ZL=DCMPLX(DBLE(L),ZERO)
      ZNL=ZN+ZL
      ZNL1=ZNL+ONE
      CALL PSI0(ZNL1,V1)
      S1=V1+GE
C
      T=T-AK1(L)*(ZN/ZNL*ZET3+ZL/ZNL**2*(ZET2-S1/ZNL))
1     CONTINUE
C
      ACG6=T
C
      RETURN
      END
      COMPLEX*16 FUNCTION ACG7(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LI3(-X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 AK1,DL,ZERO,ONE,ZET2,ZET3
C
!      COMMON/ACLOG1/ AK1(9)
!      COMMON/ACCON1/ DL
!      COMMON/ACCON2/ ZET2
!      COMMON/ACCON3/ ZET3
!      COMMON/VAL   / ZERO,ONE
c     DIMENSION AK1(9)
      DIMENSION AK1(20)
C
c      ZET2=1.64493406685d0
      ZET2 = 1.64493406684822643647d0
      ZET3=1.20205690315959428540D+0
      ZERO=0.0D0
      ONE =1.0D0
      DL = DLOG(2.0D0)
C
!      DATA AK1/0.999999974532240D+0,
!     &        -0.499995525890027D+0,
!     &         0.333203435554182D+0,
!     &        -0.248529457735332D+0,
!     &         0.191451164493502D+0,
!     &        -0.137466222203386D+0,
!     &         0.792107405737825D-1,
!     &        -0.301109652783781D-1,
!     &         0.538406198111749D-2/
C
      DATA AK1/0.9999999999999925D-0,
     &     -0.4999999999988568D-0,
     &      0.3333333332641123D-0,
     &     -0.2499999977763199D-0,
     &      0.1999999561535526D-0,
     &     -0.1666660875051348D-0,
     &      0.1428517138099479D-0,
     &     -0.1249623936313475D-0,
     &      0.1109128496887138D-0,
     &     -0.9918652787800788D-1,
     &      0.8826572954250856D-1,
     &     -0.7643209265133132D-1,
     &      0.6225829212455825D-1,
     &     -0.4572477090315515D-1,
     &      0.2890194939889559D-1,
     &     -0.1496621145891488D-1,
     &      0.6003156359511387D-2,
     &     -0.1731328252868496D-2,
     &      0.3172112728405899D-3,
     &     -0.2760099875146713D-4/
      
      T=DCMPLX(-3.0D0*ZET3/4.0D0*DL,ZERO)
C
      DO 1 L=1,20 !9
      ZL =DCMPLX(DBLE(L),ZERO)
      ZNL =ZN+ZL
      ZNL1=ZNL+DCMPLX(ONE,ZERO)
      CALL BET(ZNL1,V1)
C
      T=T+AK1(L)*(ZN/ZNL*3.0D0*ZET3/4.0D0+ZL/ZNL**2/2.0D0*ZET2
     & -ZL/ZNL**3*(DL-V1))
1     CONTINUE
C
      ACG7=T
C
      RETURN
      END
      COMPLEX*16 FUNCTION ACG8(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: S12(X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 AK1,DL,ZERO,ONE,ZET2,GE,ZET3
C
!      COMMON/ACLOG1/ AK1(9)
!      COMMON/ACCON1/ DL
!      COMMON/ACCON2/ ZET2
!      COMMON/ACCON3/ ZET3
!      COMMON/ACCON4/ GE
!      COMMON/VAL   / ZERO,ONE
c      DIMENSION AK1(9)
      DIMENSION AK1(20)
C
c      ZET2=1.64493406685d0
      ZET2 = 1.64493406684822643647d0
      ZET3=1.20205690315959428540D+0
      ZERO=0.0D0
      ONE =1.0D0
      GE   = 0.57721566490153D+0
      DL = DLOG(2.0D0)
C
!      DATA AK1/0.999999974532240D+0,
!     &        -0.499995525890027D+0,
!     &         0.333203435554182D+0,
!     &        -0.248529457735332D+0,
!     &         0.191451164493502D+0,
!     &        -0.137466222203386D+0,
!     &         0.792107405737825D-1,
!     &        -0.301109652783781D-1,
!     &         0.538406198111749D-2/
C
      DATA AK1/0.9999999999999925D-0,
     &     -0.4999999999988568D-0,
     &      0.3333333332641123D-0,
     &     -0.2499999977763199D-0,
     &      0.1999999561535526D-0,
     &     -0.1666660875051348D-0,
     &      0.1428517138099479D-0,
     &     -0.1249623936313475D-0,
     &      0.1109128496887138D-0,
     &     -0.9918652787800788D-1,
     &      0.8826572954250856D-1,
     &     -0.7643209265133132D-1,
     &      0.6225829212455825D-1,
     &     -0.4572477090315515D-1,
     &      0.2890194939889559D-1,
     &     -0.1496621145891488D-1,
     &      0.6003156359511387D-2,
     &     -0.1731328252868496D-2,
     &      0.3172112728405899D-3,
     &     -0.2760099875146713D-4/
      
!
!
      T=DCMPLX(ZERO,ZERO)
      DO 1 L=1,20 !9
      ZL=DCMPLX(DBLE(L),ZERO)
      ZNL=ZN+ZL
      ZNL1=ZNL+ONE
      CALL PSI0(ZNL1,PS0)
      CALL PSI1(ZNL1,PS1)
      S1=PS0+GE
      S2=ZET2-PS1
      T=T-AK1(L)*(ZN*ZET3/ZNL+ZL/ZNL**2/2.0D0*(S1**2+S2))
1     CONTINUE
C
      ACG8=T+DL*ZET3
C
      RETURN
      END
      COMPLEX*16 FUNCTION ACG9(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: S12(-X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 AK1,AK2,AK3,DL,ZERO,ONE,ZET3
C
!      COMMON/ACLOG1/ AK1(9)
!      COMMON/ACLOG2/ AK2(10)
!      COMMON/ACLOG3/ AK3(11)
!      COMMON/ACCON1/ DL
!      COMMON/ACCON3/ ZET3
!      COMMON/VAL   / ZERO,ONE
 !     DIMENSION AK1(9)
 !     DIMENSION AK2(10)
      DIMENSION AK1(20)
      DIMENSION AK2(23)
       DIMENSION AK3(11)
C
      ZET3=1.20205690315959428540D+0
      ZERO=0.0D0
      ONE =1.0D0
      DL = DLOG(2.0D0)
C
!      DATA AK1/0.999999974532240D+0,
!     &        -0.499995525890027D+0,
!     &         0.333203435554182D+0,
!     &        -0.248529457735332D+0,
!     &         0.191451164493502D+0,
!     &        -0.137466222203386D+0,
!     &         0.792107405737825D-1,
!     &        -0.301109652783781D-1,
!     &         0.538406198111749D-2/
!      DATA AK2/0.999999980543793D+0,
!     &        -0.999995797779624D+0,
!     &         0.916516447393493D+0,
!     &        -0.831229921350708D+0,
!     &         0.745873737923571D+0,
!     &        -0.634523908078600D+0,
!     &         0.467104011423750D+0,
!     &        -0.261348046799178D+0,
!     &         0.936814286867420D-1,
!     &        -0.156249375012462D-1/
      DATA AK3/9.99999989322696D-1,
     &        -1.49999722020708D+0,
     &         1.74988008499745D+0,
     &        -1.87296689068405D+0,
     &         1.91539974617231D+0,
     &        -1.85963744001295D+0,
     &         1.62987195424434D+0,
     &        -1.17982353224299D+0,
     &         6.28710122994999D-1,
     &        -2.11307487211713D-1,
     &         3.28953352932140D-2/
C
      DATA AK1/0.9999999999999925D-0,
     &     -0.4999999999988568D-0,
     &      0.3333333332641123D-0,
     &     -0.2499999977763199D-0,
     &      0.1999999561535526D-0,
     &     -0.1666660875051348D-0,
     &      0.1428517138099479D-0,
     &     -0.1249623936313475D-0,
     &      0.1109128496887138D-0,
     &     -0.9918652787800788D-1,
     &      0.8826572954250856D-1,
     &     -0.7643209265133132D-1,
     &      0.6225829212455825D-1,
     &     -0.4572477090315515D-1,
     &      0.2890194939889559D-1,
     &     -0.1496621145891488D-1,
     &      0.6003156359511387D-2,
     &     -0.1731328252868496D-2,
     &      0.3172112728405899D-3,
     &     -0.2760099875146713D-4/
      DATA AK2/1.0000000000000000D-0,
     & -0.9999999999999985D-0,
     &  0.9166666666663948D-0,
     & -0.8333333333136118D-0,
     &  0.7611111103508889D-0,
     & -0.6999999819735105D-0,
     &  0.6482139985629993D-0,
     & -0.6039649964806160D-0,
     &  0.5657662306410356D-0,
     & -0.5323631718571445D-0,
     &  0.5024238774786239D-0,
     & -0.4738508288315496D-0,
     &  0.4427472719775835D-0,
     & -0.4029142806330511D-0,
     &  0.3476841543351489D-0,
     & -0.2748590021353420D-0,
     &  0.1915627642585285D-0,
     & -0.1130763066428224D-0,
     &  0.5415661067306229D-1,
     & -0.1999877298940919D-1,
     &  0.5303624439388411D-2,
     & -0.8944156375768203D-3,
     &  0.7179502917974332D-4/

      
      T=DCMPLX(ZET3*DL/8.0D0,ZERO)
      DO 1 K=1,20 !9
      T1=DCMPLX(ZERO,ZERO)
      DO 2 L=2,24 !11
      L1=L-1
      ZNKL=ZN+DCMPLX(DBLE(K+L),ZERO)
      T1=T1+AK2(L1)/ZNKL
2     CONTINUE
      ZK=DCMPLX(DBLE(K),ZERO)
      ZNK=ZN+ZK
      T=T-AK1(K)*ZN/ZNK*(ZET3/8.0D0-T1/2.0D0)
1     CONTINUE
      T2=DCMPLX(ZERO,ZERO)
      DO 3 K=3,13
      L=K-2
      ZK=DCMPLX(DBLE(K),ZERO)
      ZNK=ZN+ZK
      T2=T2+AK3(L)/ZNK
3     CONTINUE
      T=T-T2/2.0D0
C
      ACG9=T
C
      RETURN
      END
      COMPLEX*16 FUNCTION ACG13(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(1+X)/(1+X)*LI2(-X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 AK2,AK3,DL,ZERO,ONE,ZET2
C
!      COMMON/ACLOG2/ AK2(10)
!      COMMON/ACLOG3/ AK3(11)
!      COMMON/ACCON1/ DL
!      COMMON/ACCON2/ ZET2
!      COMMON/VAL   / ZERO,ONE
c     DIMENSION AK2(10)
      DIMENSION AK2(23)
      DIMENSION AK3(11)
C
c      ZET2=1.64493406685d0
      ZET2 = 1.64493406684822643647d0
      ZERO=0.0D0
      ONE =1.0D0
      DL = DLOG(2.0D0)
C
!      DATA AK2/0.999999980543793D+0,
!     &        -0.999995797779624D+0,
!     &         0.916516447393493D+0,
!     &        -0.831229921350708D+0,
!     &         0.745873737923571D+0,
!     &        -0.634523908078600D+0,
!     &         0.467104011423750D+0,
!     &        -0.261348046799178D+0,
!     &         0.936814286867420D-1,
!     &        -0.156249375012462D-1/
      DATA AK2/1.0000000000000000D-0,
     & -0.9999999999999985D-0,
     &  0.9166666666663948D-0,
     & -0.8333333333136118D-0,
     &  0.7611111103508889D-0,
     & -0.6999999819735105D-0,
     &  0.6482139985629993D-0,
     & -0.6039649964806160D-0,
     &  0.5657662306410356D-0,
     & -0.5323631718571445D-0,
     &  0.5024238774786239D-0,
     & -0.4738508288315496D-0,
     &  0.4427472719775835D-0,
     & -0.4029142806330511D-0,
     &  0.3476841543351489D-0,
     & -0.2748590021353420D-0,
     &  0.1915627642585285D-0,
     & -0.1130763066428224D-0,
     &  0.5415661067306229D-1,
     & -0.1999877298940919D-1,
     &  0.5303624439388411D-2,
     & -0.8944156375768203D-3,
     &  0.7179502917974332D-4/

      DATA AK3/9.99999989322696D-1,
     &        -1.49999722020708D+0,
     &         1.74988008499745D+0,
     &        -1.87296689068405D+0,
     &         1.91539974617231D+0,
     &        -1.85963744001295D+0,
     &         1.62987195424434D+0,
     &        -1.17982353224299D+0,
     &         6.28710122994999D-1,
     &        -2.11307487211713D-1,
     &         3.28953352932140D-2/
C
      T0=DCMPLX(-1.0D0/4.0D0*ZET2*DL**2,ZERO)
C
      T1=DCMPLX(ZERO,ZERO)
      DO 1 L=3,13
      K=L-2
      T1=T1+AK3(K)/(ZN+DBLE(L))
1     CONTINUE
C
      T2=DCMPLX(ZERO,ZERO)
      DO 2 L=2,24 !11
      K=L-1
      ZNK1=ZN+DBLE(L+1)
      CALL BET(ZNK1,V1)
      T2=T2+AK2(K)*ZN/(ZN+DBLE(L))*(ZET2/2.0D0-(DL-V1)/(ZN+DBLE(L)))
2     CONTINUE
C
      ACG13=T0+(T1+T2)/2.0D0
C
      RETURN
      END
      COMPLEX*16 FUNCTION ACG20(ZZ)
C     -----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: (LI3(X)-ZETA3)/(X-1)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 DL,ZERO,ONE,ZET2,ZET3,GE,P24,P34,P22,CK2,CK4,P12,P14
C
      DIMENSION P24(5),P34(5),P22(4),CK2(13),CK4(13),P12(4),P14(5)
!
      DATA CK2/ .0,.0,2.1253046159349207, -1.0523034804772378,
!     & 0.480322239287449D+0,
!     &         -0.168480825099580D+1,
!     &          0.209270571620726D+1,
!     &         -0.101728150275998D+1,
!
     &          0.160179976133047D+0,
     &         -0.351982236677917D+0,
     &          0.141033316846244D+1,
     &         -0.353343997843391D+1,
     &          0.593934696819832D+1,
     &         -0.660019784804042D+1,
     &          0.466330349413063D+1,
     &         -0.189825467489058D+1,
     &          0.339772909487512D+0/
      DATA CK4/ .0,2.215008697869307,-0.9133677154535804,
     /3.4783104357500143,-2.823955592989266,
!     &          0.192962504274437D+0,
!     &          0.000005641557253D+0,
!     &         -0.196891075399448D+1,
!     &          0.392919138747074D+1,
!     &         -0.290306105685546D+1,
!
     &          0.992890266001707D+0,
     &         -0.130026190226546D+1,
     &          0.341870577921103D+1,
     &         -0.576763902370864D+1,
     &          0.645554138192407D+1,
     &         -0.459405622046138D+1,
     &          0.188510809558304D+1,
     &         -0.340476080290674D+0/
c      ZET2=1.64493406685d0
      ZET2 = 1.64493406684822643647d0
      ZET3=1.20205690315959428540D+0
      ZERO=0.0D0
      ONE =1.0D0
      GE   = 0.57721566490153D+0
      DL = DLOG(2.0D0)
!
      P12(1)=ZET3-11.0D0/6.0*ZET2+4.0D0/3.0D0
      P12(2)=3.0D0*ZET2-13.0D0/4.0D0
      P12(3)=-3.0D0/2.0D0*ZET2+5.0D0/2.0D0
      P12(4)=1.0D0/3.0D0*ZET2-7.0D0/12.0D0
!
      P14(1)=257.D0/144.0D0-205.0D0/72.0D0*ZET2+ZET2**2
      P14(2)=-167.0D0/36.0D0+25.0D0/6.0D0*ZET2
      P14(3)=101.0D0/24.0D0-23.0D0/12.0D0*ZET2
      P14(4)=-59.0D0/36.0D0+13.0D0/18.0D0*ZET2
      P14(5)=41.0D0/144.0D0-ZET2/8.0D0

!      CK2(1)=CK2(1)+P12(1)
!      CK2(2)=CK2(2)+P12(2)
!      CK2(3)=CK2(3)+P12(3)
!      CK2(4)=CK2(4)+P12(4)
C
!      CK4(1)=CK4(1)+P14(1)
!      CK4(2)=CK4(2)+P14(2)
!      CK4(3)=CK4(3)+P14(3)
!      CK4(4)=CK4(4)+P14(4)
!      CK4(5)=CK4(5)+P14(5)

C
!      COMMON/ACLOG7/ CK2(13)
!      COMMON/ACLOG9/ CK4(13)
!      COMMON/POLY2 / P22(4)
!      COMMON/POLY4 / P24(5),P34(5)
!      COMMON/ACCON1/ DL
!      COMMON/ACCON2/ ZET2
!      COMMON/ACCON3/ ZET3
!      COMMON/ACCON4/ GE
!      COMMON/VAL   / ZERO,ONE
!
!
      P22(1)=-1.0D0
      P22(2)=5.0D0/2.0D0
      P22(3)=-2.0D0
      P22(4)=1.0D0/2.0D0
!
      P24(1)=-167.0D0/36.0D0+25.0D0/6.0D0*ZET2
      P24(2)=235.0D0/18.0D0-8.0D0*ZET2
      P24(3)=-40.0D0/3.0D0+6.0D0*ZET2
      P24(4)=109.0D0/18.0D0-8.0D0/3.0D0*ZET2
      P24(5)=-41.0D0/36.0D0+ZET2/2.0D0
C-------------------------------                       
      P34(1)=35.0D0/12.0D0
      P34(2)=-26.0D0/3.0D0
      P34(3)=19.0D0/2.0D0
      P34(4)=-14.0D0/3.0D0
      P34(5)=11.0D0/12.0D0
!

C
      T=DCMPLX(ZET2**2/2.0D0,ZERO)
C
      ZN = ZZ
      ZN1=ZN+ONE
      CALL PSI0(ZN1,PS)
      S1=PS+GE
      T=T-ZET3*S1
C
      DO 1 K=1,13
      L=K-1
      ZNK=ZN+DBLE(L)
      ZNK1=ZNK+ONE
      CALL PSI0(ZNK1,PS)
      S1=PS+GE
      R=ZN/ZNK
C
      T=T+CK2(K)*R*S1
1     CONTINUE
C
      DO 3 K=1,4
      L=K-1
      ZNK=ZN+DBLE(L)
      ZNK1=ZNK+ONE
      CALL PSI0(ZNK1,PS)
      CALL PSI1(ZNK1,PS1)
      S1=PS+GE
      S2=-PS1+ZET2
      R=ZN/ZNK
      T1=S1**2+S2
C
      T=T-P22(K)*R*T1
3     CONTINUE
C
      DO 4 K=1,13
      L=K-1
      ZNK=ZN+DBLE(L)
      ZNK1=ZNK+ONE
      CALL PSI0(ZNK1,PS)
      S1=PS+GE
C
      T=T-CK4(K)/ZNK*ZN/2.0D0
4     CONTINUE
C
      DO 5 K=1,5
      L=K-1
      ZNK=ZN+DBLE(L)
      ZNK1=ZNK+ONE
      CALL PSI0(ZNK1,PS)
      CALL PSI1(ZNK1,PS1)
      S1=PS+GE
      S2=-PS1+ZET2
      T1=-S1/ZNK
      T2=(S1**2+S2)/ZNK
C
      T=T-(P24(K)*T1+P34(K)*T2)*ZN/2.0D0
5     CONTINUE
C
C
      ACG20=T

C
      RETURN
      END
!
      COMPLEX*16 FUNCTION ACG21(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: (S12(X)-ZETA3)/(X-1)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 ZERO,ONE,ZET2,ZET3,GE,CK3,P13,P23,P33,ZETA2,ZETA3
C
      DIMENSION P13(5),P23(5),P33(5),CK3(10)
!      COMMON/ACLOG8/ CK3(10)
!      COMMON/POLY3 / P23(5),P33(5)
!      COMMON/ACCON2/ ZET2
!      COMMON/ACCON3/ ZET3
!      COMMON/ACCON4/ GE
!      COMMON/VAL   / ZERO,ONE
C
!
c      ZET2=1.64493406685d0
      ZET2 = 1.64493406684822643647d0
      ZET3=1.20205690315959428540D+0
      ZERO=0.0D0
      ONE =1.0D0
      GE   = 0.57721566490153D+0


      DATA CK3/ .0,1.4236162474052558,-0.08001203559240111,!!!!!!! ERRR ,0.249849075518710,
     /-0.39875367195395994,0.339241791547134,
!      DATA CK3/
!     &         -0.243948949064443D-1,
!     &          0.000005136294145D+0,
!     &          0.249849075518710D+0,
!     &         -0.498290708990997D+0,
!     &          0.354866791547134D+0,
!
     &         -0.522116678353452D-1,
     &         -0.648354706049337D-1,
     &          0.644165053822532D-1,
     &         -0.394927322542075D-1,
     &          0.100879370657869D-1/


      P13(1)=ZET3-2035.0D0/1728.0D0
      P13(2)=205.0D0/144.0D0
      P13(3)=-95.0D0/288.0D0
      P13(4)=43.0D0/432.0D0
      P13(5)=-1.0D0/64.0D0

      P23(1)=205.0D0/144.0D0
      P23(2)=-25.0D0/12.0D0
      P23(3)=23.0D0/24.0D0
      P23(4)=-13.0D0/36.0D0
      P23(5)=1.0D0/16.0D0

      P33(1)=-25.0D0/24.0D0
      P33(2)=2.0D0
      P33(3)=-3.0D0/2.0D0
      P33(4)=2.0D0/3.0D0
      P33(5)=-1.0D0/8.0D0

!      CK3(1)=CK3(1)+P13(1)
!      CK3(2)=CK3(2)+P13(2)
!      CK3(3)=CK3(3)+P13(3)
!      CK3(4)=CK3(4)+P13(4)
!      CK3(5)=CK3(5)+P13(5)

!
      T=DCMPLX(ZERO,ZERO)
      ZN1=ZN+ONE
      CALL PSI0(ZN1,PS)
      CALL PSI1(ZN1,PS1)
      CALL PSI2(ZN1,PS2)
      S1=PS+GE
      S2=-PS1+ZET2
      S3= PS2/2.0D0+ZET3
C
      T=T-ZET3*S1
      T=T+(S1**3+3.0D0*S1*S2+2.0D0*S3)/2.0D0/ZN
C
      DO 1 K=1,10
      L=K-1
      ZNK=ZN+DBLE(L)
      ZNK1=ZNK+ONE
      CALL PSI0(ZNK1,PS)
      S1=PS+GE
      T=T+CK3(K)*S1*ZN/ZNK
1     CONTINUE
C
      DO 2 K=1,5
      L=K-1
      ZNK=ZN+DBLE(L)
      ZNK1=ZNK+ONE
      CALL PSI0(ZNK1,PS)
      CALL PSI1(ZNK1,PS1)
      CALL PSI2(ZNK1,PS2)
      S1=PS+GE
      S2=-PS1+ZET2
      S3= PS2/2.0D0+ZET3
      T2=S1**2+S2
      T3=S1**3+3.0D0*S1*S2+2.0D0*S3
      T=T+ZN/ZNK*(P33(K)*T3-P23(K)*T2)
2     CONTINUE
C
      ACG21=T
C
      RETURN
      END
C     END BLUMLEIN FUNCTIONS

