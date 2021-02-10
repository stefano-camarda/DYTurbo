C Computes the various C2 coefficients at a given value of XN

      subroutine dyH2calc(C2qg,C2NSqqb,C2NSqq,C2Sqqb,XN)
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
       
      plog41f2=0.51747906167389938633d0
      log2=    0.69314718055994530942d0
       pi=3.141592653589793d0

       GE   = 0.57721566490153286061d0
       ONE=CMPLX(1.0D0,0.0D0)
       CF=4d0/3d0
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

      PS0NP1     = DYPSI0(XNP1    )       
      PS0NP2     = DYPSI0(XNP2    )
      PS0NP3     = DYPSI0(XNP3    )
      PS0N       = DYPSI0(XN      )
      PS0Nf2P1   = DYPSI0(XNf2P1  )
      PS0Nf2P3f2 = DYPSI0(XNf2P3f2)
      PS0Nf2P1f2 = DYPSI0(XNf2P1f2)
      PS0Nf2     = DYPSI0(XNf2    )
      PS0Nf2P2   = DYPSI0(XNf2P2  )

      PS1NP1     = DYPSI1(XNP1    )
      PS1NP2     = DYPSI1(XNP2    )
      PS1NP3     = DYPSI1(XNP3    )
      PS1N       = DYPSI1(XN      )
      PS1Nf2     = DYPSI1(XNf2    )
      PS1Nf2P1   = DYPSI1(XNf2P1  )
      PS1Nf2P3f2 = DYPSI1(XNf2P3f2)
      PS1Nf2P1f2 = DYPSI1(XNf2P1f2)
      PS1Nf2P2   = DYPSI1(XNf2P2  )

      PS2NP1     = DYPSI2(XNP1    )
      PS2NP2     = DYPSI2(XNP2    )
      PS2NP3     = DYPSI2(XNP3    )
      PS2N       = DYPSI2(XN      )
      PS2Nf2     = DYPSI2(XNf2    )
      PS2Nf2P1   = DYPSI2(XNf2P1  )
      PS2Nf2P3f2 = DYPSI2(XNf2P3f2)
      PS2Nf2P1f2 = DYPSI2(XNf2P1f2)
      PS2Nf2P2   = DYPSI2(XNf2P2  )

      PS3N       = DYPSI3(XN      )
      PS3Nf2     = DYPSI3(XNf2    )
      PS3Nf2P1f2 = DYPSI3(XNf2P1f2)
      PS3Nf2     = DYPSI3(XNf2    )

      BETA0N     = DYBE0(XN       )
      BETA0NP1   = DYBE0(XNP1     )
      BETA0NP2   = DYBE0(XNP2     )
      BETA0NP3   = DYBE0(XNP3     )
      BETA1N     = DYBE1(XN       )
      BETA2N     = DYBE2(XN       )
      BETA3N     = DYBE3(XN       )

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
      S211NM1= -fun16(XNM1)+6d0/5*zeta2**2
!      S112NM1=0
!      S121NM1=0

      S31NM1= fun15(XNM1)+zeta2*S2NM1-1d0/2*zeta2**2


      Sm1NM1= (1)*BETA0N-log2
      Sm1N=   (-1)*BETA0NP1-log2
      Sm1NP1= (1)*BETA0NP2-log2
      Sm1NP2= (-1)*BETA0NP3-log2
      Sm2NM1= (-1)*BETA1N-zeta2/2d0
      Sm3NM1= (1)*BETA2N/2d0-zeta3*3/4d0
      Sm4NM1= (-1)*BETA3N/6d0-zeta4*7/8d0

      Sm21N= fun6(XN)*(1)+zeta2*Sm1N
     /       -5/8d0*zeta3+zeta2*log2
      Sm21NP1= fun6(XNP1)*(-1)+zeta2*Sm1NP1
     /       -5/8d0*zeta3+zeta2*log2
      Sm21NP2= fun6(XNP2)*(1)+zeta2*Sm1NP2
     /       -5/8d0*zeta3+zeta2*log2

      Sm31NM1= fun10(XNM1)*(1)+zeta2*Sm2NM1-zeta3*Sm1NM1
     /       -3/5d0*zeta2**2+2*plog41f2+3/4d0*zeta3*log2
     /       -zeta2/2d0*log2**2+log2**4/12d0
      Sm22NM1= fun9(XNM1)*(1)-2d0*Sm31NM1
     /         +2d0*zeta2*Sm2NM1+3/40d0*zeta2**2
      Sm211NM1= fun12(XNM1)*(-1)+zeta3*Sm1NM1-plog41f2+zeta2**2/8d0
     /     +zeta3*log2/8d0+zeta2*log2**2/4d0-log2**4/24d0


      ACG13XNM1 = fun14(XNM1)
      ACG1XNM1  = fun1(XNM1)
      ACG1pXN   = fun2(XN)
      ACG1pXNP1 = fun2(XNP1)
      ACG1ppXNM1= fun3(XNM1)
      ACG2pXNM1 = fun5(XNM1)
      ACG4XN    = fun7(XN)
      ACG4XNP1  = fun7(XNP1)
      ACG4pXNM1 = fun8(XNM1)
      ACG5XNM1  = fun9(XNM1)
      ACG6XNM1  = fun10(XNM1)
      ACG7XNM1  = fun11(XNM1)
      ACG9XNM1  = fun13(XNM1)


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
