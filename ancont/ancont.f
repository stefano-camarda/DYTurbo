      PROGRAM ANCONT
C     --------------
C
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C     CODE ANCONT    1.00
C
C     TO EVALUATE THE MELLIN TRANSFORMS OF NIELSEN INTEGRALS AND RELATED
C     FUNCTIONS AND THEIR ANALYTIC CONTINUATIONS INTO THE COMPLEX PLANE
C     FOR ALL HARMONIC SUMS OCCURING UP TO THE 2--LOOP LEVEL 
C     (TRANSCENDENTALITY 4)
C
C************************************************************************
C
C     REFS.: [1] J. BLUEMLEIN AND S. KURTH, DESY 97-160, HEP-PH/9708388
C            [2] J. BLUEMLEIN AND S. KURTH, PHYS. REV. D60 (1999) 014018
C            [3] J. BLUEMLEIN Comput. Phys. Commun. 133 (2000) 76.
C
C   >>>      Conditions of use: Citation of Refs.: [2] and [3]
C   >>>      Compilation g77 ancont.f
C
C************************************************************************
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C--------------     START
C
      CALL ACSTAR
C
C
C--------------     INITIALIZATION
C
      CALL ACINI
C
C--------------     RUNNING
C
      CALL ACRUN
C
C--------------     END
C
      CALL ACEND
C
      STOP
      END
      SUBROUTINE ACSTAR
C     -----------------
C
      WRITE(6,*) '                                   '
      WRITE(6,*) '***** ANCONT: VS 1.00         *****'
      WRITE(6,*) '***** J. BLUEMLEIN 01.10.1999 *****'
      WRITE(6,*) '                                   '
      WRITE(6,*) '                                   '
C
      RETURN
      END
      SUBROUTINE ACINI
C     ----------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  INITIALIZE CONSTANTS
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /KM/KEY,MAX
C
      COMMON/ACCON1/ DL
      COMMON/ACCON2/ ZET2
      COMMON/ACCON3/ ZET3
      COMMON/ACCON4/ GE
      COMMON/ACCON5/ XLI4
      COMMON/IAPP  / IAPP
      COMMON/VAL   / ZERO,ONE
      COMMON/ACLOG1/ AK1(9)
      COMMON/ACLOG2/ AK2(10)
      COMMON/ACLOG3/ AK3(11)
      COMMON/ACLOG4/ BK1(9)
      COMMON/ACLOG5/ BK2(10)
      COMMON/ACLOG6/ CK1(12)
      COMMON/ACLOG7/ CK2(13)
      COMMON/ACLOG8/ CK3(10)
      COMMON/ACLOG9/ CK4(13)
      COMMON/ACLG10/ CK5(9)
      COMMON/ACLG11/ DK5(10)
      COMMON/POLY1 / P21(4)
      COMMON/POLY2 / P22(4)
      COMMON/POLY3 / P23(5),P33(5)
      COMMON/POLY4 / P24(5),P34(5)
C
      DIMENSION P11(4),P12(4),P13(5),P14(5)
C
      DATA AK1/0.999999974532240D+0,
     &        -0.499995525890027D+0,
     &         0.333203435554182D+0,
     &        -0.248529457735332D+0,
     &         0.191451164493502D+0,
     &        -0.137466222203386D+0,
     &         0.792107405737825D-1,
     &        -0.301109652783781D-1,
     &         0.538406198111749D-2/
C
      DATA AK2/0.999999980543793D+0,
     &        -0.999995797779624D+0,
     &         0.916516447393493D+0,
     &        -0.831229921350708D+0,
     &         0.745873737923571D+0,
     &        -0.634523908078600D+0,
     &         0.467104011423750D+0,
     &        -0.261348046799178D+0,
     &         0.936814286867420D-1,
     &        -0.156249375012462D-1/
C
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
      DATA BK1/0.693147166991375D+0,
     &        -0.306850436868254D+0,
     &         0.193078041088284D+0,
     &        -0.139403892894644D+0,
     &         0.105269615988049D+0,
     &        -0.746801353858524D-1,
     &         0.427339135378207D-1,
     &        -0.161809049989783D-1,
     &         0.288664611077007D-2/
C
      DATA BK2/0.480453024731510D+0,
     &         0.480450679641120D+0,
     &        -0.519463586324817D+0,
     &         0.479285947990175D+0,
     &        -0.427765744446172D+0,
     &         0.360855321373065D+0,
     &        -0.263827078164263D+0,
     &         0.146927719341510D+0,
     &        -0.525105367350968D-1,
     &         0.874144396622167D-2/
C
      DATA CK1/-0.283822933724932D+0,
     &          0.999994319023731D+0,
     &         -0.124975762907682D+1,
     &          0.607076808008983D+0,
     &         -0.280403220046588D-1,
     &         -0.181869786537805D+0,
     &          0.532318519269331D+0,
     &         -0.107281686995035D+1,
     &          0.138194913357518D+1,
     &         -0.111100841298484D+1,
     &          0.506649587198046D+0,
     &         -0.100672390783659D+0/
C
      DATA CK2/ 0.480322239287449D+0,
     &         -0.168480825099580D+1,
     &          0.209270571620726D+1,
     &         -0.101728150275998D+1,
     &          0.160179976133047D+0,
     &         -0.351982236677917D+0,
     &          0.141033316846244D+1,
     &         -0.353343997843391D+1,
     &          0.593934696819832D+1,
     &         -0.660019784804042D+1,
     &          0.466330349413063D+1,
     &         -0.189825467489058D+1,
     &          0.339772909487512D+0/
C
      DATA CK3/-0.243948949064443D-1,
     &          0.000005136294145D+0,
     &          0.249849075518710D+0,
     &         -0.498290708990997D+0,
     &          0.354866791547134D+0,
     &         -0.522116678353452D-1,
     &         -0.648354706049337D-1,
     &          0.644165053822532D-1,
     &         -0.394927322542075D-1,
     &          0.100879370657869D-1/
C
      DATA CK4/ 0.192962504274437D+0,
     &          0.000005641557253D+0,
     &         -0.196891075399448D+1,
     &          0.392919138747074D+1,
     &         -0.290306105685546D+1,
     &          0.992890266001707D+0,
     &         -0.130026190226546D+1,
     &          0.341870577921103D+1,
     &         -0.576763902370864D+1,
     &          0.645554138192407D+1,
     &         -0.459405622046138D+1,
     &          0.188510809558304D+1,
     &         -0.340476080290674D+0/
C
      DATA CK5/-0.822467033400775D+0,
     &          0.887664705657325D-1,
     &         -0.241549406045162D-1,
     &          0.965074750946139D-2,
     &         -0.470587487919749D-2,
     &          0.246014308378549D-2,
     &         -0.116431121874067D-2,
     &          0.395705193848026D-3,
     &         -0.664699010014505D-4/
C
      DATA DK5/-0.822467033400775D+0,
     &          0.999999974532240D+0,
     &         -0.249997762945014D+0,
     &          0.111067811851394D+0,
     &         -0.621323644338330D-1,
     &          0.382902328987004D-1,
     &         -0.229110370338977D-1,
     &          0.113158200819689D-1,
     &         -0.376387065979726D-2,
     &          0.598229109013054D-3/
C
      PI = 3.141592653589793238462643D0
      ZETA2 = PI**2/6.0D0       
      ZETA3 = 1.20205690315959428540D0
      ZET2=ZETA2
      ZET3=ZETA3
      ZLI4 = FLI4(0.5D0)
      XLI4 = FLI4(0.5D0)
      D2 = DLOG(2.0D0)
      DL = DLOG(2.0D0)
      GE   = 0.57721566490153D+0
C
      ZERO=0.0D0
      ONE =1.0D0
C
C----------------------------------------------------------------------
      P11(1)=-49.0D0/36.0+ZETA2
      P11(2)=11.0D0/6.0D0
      P11(3)=-7.0D0/12.0D0
      P11(4)=1.0D0/9.0D0
C-------------------------------                       
      P21(1)=11.0D0/6.0
      P21(2)=-3.0D0
      P21(3)=3.0D0/2.0D0
      P21(4)=-1.0D0/3.0D0
C-------------------------------                       
      P12(1)=ZETA3-11.0D0/6.0*ZETA2+4.0D0/3.0D0
      P12(2)=3.0D0*ZETA2-13.0D0/4.0D0
      P12(3)=-3.0D0/2.0D0*ZETA2+5.0D0/2.0D0
      P12(4)=1.0D0/3.0D0*ZETA2-7.0D0/12.0D0
C-------------------------------                       
      P22(1)=-1.0D0
      P22(2)=5.0D0/2.0D0
      P22(3)=-2.0D0
      P22(4)=1.0D0/2.0D0
C-------------------------------                       
      P13(1)=ZETA3-2035.0D0/1728.0D0
      P13(2)=205.0D0/144.0D0
      P13(3)=-95.0D0/288.0D0
      P13(4)=43.0D0/432.0D0
      P13(5)=-1.0D0/64.0D0
C-------------------------------                       
      P23(1)=205.0D0/144.0D0
      P23(2)=-25.0D0/12.0D0
      P23(3)=23.0D0/24.0D0
      P23(4)=-13.0D0/36.0D0
      P23(5)=1.0D0/16.0D0
C-------------------------------                       
      P33(1)=-25.0D0/24.0D0
      P33(2)=2.0D0
      P33(3)=-3.0D0/2.0D0
      P33(4)=2.0D0/3.0D0
      P33(5)=-1.0D0/8.0D0
C-------------------------------                       
      P14(1)=257.D0/144.0D0-205.0D0/72.0D0*ZET2+ZET2**2
      P14(2)=-167.0D0/36.0D0+25.0D0/6.0D0*ZET2
      P14(3)=101.0D0/24.0D0-23.0D0/12.0D0*ZET2
      P14(4)=-59.0D0/36.0D0+13.0D0/18.0D0*ZET2
      P14(5)=41.0D0/144.0D0-ZET2/8.0D0
C-------------------------------                       
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
C----------------------------------------------------------------------
C
C >>> ACCOUNT FOR POLYNOM PARTS
C
      CK1(1)=CK1(1)+P11(1)
      CK1(2)=CK1(2)+P11(2)
      CK1(3)=CK1(3)+P11(3)
      CK1(4)=CK1(4)+P11(4)
C
      CK2(1)=CK2(1)+P12(1)
      CK2(2)=CK2(2)+P12(2)
      CK2(3)=CK2(3)+P12(3)
      CK2(4)=CK2(4)+P12(4)
C
      CK3(1)=CK3(1)+P13(1)
      CK3(2)=CK3(2)+P13(2)
      CK3(3)=CK3(3)+P13(3)
      CK3(4)=CK3(4)+P13(4)
      CK3(5)=CK3(5)+P13(5)
C
      CK4(1)=CK4(1)+P14(1)
      CK4(2)=CK4(2)+P14(2)
      CK4(3)=CK4(3)+P14(3)
      CK4(4)=CK4(4)+P14(4)
      CK4(5)=CK4(5)+P14(5)
C
C----------------------------------------------------------------------
C
      KEY = 2
      MAX = 10000
C
      IAPP=1
C
      CALL INVINI
C
      CALL DEFAUL
C
      CALL UINIT
C
      CALL WROUT
C
      RETURN
      END
      SUBROUTINE INVINI
C     -----------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  INITIALIZE CONSTANTS: MELLIN INVERSION
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
*
      COMMON /OFFS  / CP(26)
      COMMON /PI    / PI
      COMMON /ZET2/ ZETA2
      COMMON /ZET3/ ZETA3
      COMMON /WGAUSS/ XL8(8),WL8(8),XL32(32),WL32(32)
*
      DATA  PI/ 3.14159265358979323846D+0/
*                 1      2      3      4      5       6       7
      DATA CP/-0.7D0,-1.7D0,-0.8D0,-0.7D0,-0.7D0,-0.7D0,-0.7D0,
     &        -1.7D0,-1.7D0,-1.7D0,-1.7D0,-1.7D0,-1.7D0, 0.3D0,
     &        -0.7D0,-0.7D0,-1.7D0, 0.3D0, 0.3D0, 0.3D0, 0.2D0,
     &        -0.5D0, 0.5D0, 0.5D0,0.5D0,-2.8D0/
*
      DATA XL8/-0.960289856497536231684D+0,
     &         -0.796666477413626739592D+0,
     &         -0.525532409916328985818D+0,
     &         -0.183434642495649804939D+0,
     &          0.183434642495649804939D+0,
     &          0.525532409916328985818D+0,
     &          0.796666477413626739592D+0,
     &          0.960289856497536231684D+0/
*
      DATA WL8/ 0.101228536290376259163D+0,
     &          0.222381034453374470518D+0,
     &          0.313706645877887287348D+0,
     &          0.362683783378361982964D+0,
     &          0.362683783378361982964D+0,
     &          0.313706645877887287348D+0,
     &          0.222381034453374470518D+0,
     &          0.101228536290376259163D+0/
*
      DATA XL32/-0.997263861849481563545D+0,
     &          -0.985611511545268335400D+0,
     &          -0.964762255587506430774D+0,
     &          -0.934906075937739689171D+0,
     &          -0.896321155766052123965D+0,
     &          -0.849367613732569970134D+0,
     &          -0.794483795967942406963D+0,
     &          -0.732182118740289680387D+0,
     &          -0.663044266930215200975D+0,
     &          -0.587715757240762329041D+0,
     &          -0.506899908932229390024D+0,
     &          -0.421351276130635345364D+0,
     &          -0.331868602282127649780D+0,
     &          -0.239287362252137074545D+0,
     &          -0.144471961582796493485D+0,
     &          -0.483076656877383162348D-1,
     &           0.483076656877383162348D-1,
     &           0.144471961582796493485D+0,
     &           0.239287362252137074545D+0,
     &           0.331868602282127649780D+0,
     &           0.421351276130635345364D+0,
     &           0.506899908932229390024D+0,
     &           0.587715757240762329041D+0,
     &           0.663044266930215200975D+0,
     &           0.732182118740289680387D+0,
     &           0.794483795967942406963D+0,
     &           0.849367613732569970134D+0,
     &           0.896321155766052123965D+0,
     &           0.934906075937739689171D+0,
     &           0.964762255587506430774D+0,
     &           0.985611511545268335400D+0,
     &           0.997263861849481563545D+0/
*
      DATA WL32/0.701861000943074980416D-2,
     &          0.162743947308230037851D-1,
     &          0.253920653084367920322D-1,
     &          0.342738629139351946052D-1,
     &          0.428358980239763007216D-1,
     &          0.509980592613025930634D-1,
     &          0.586840934788972562042D-1,
     &          0.658222227761892440852D-1,
     &          0.723457941088502071174D-1,
     &          0.781938957870682705896D-1,
     &          0.833119242269464135758D-1,
     &          0.876520930044038163114D-1,
     &          0.911738786957638874594D-1,
     &          0.938443990808045654352D-1,
     &          0.956387200792748594152D-1,
     &          0.965400885147278005676D-1,
     &          0.965400885147278005676D-1,
     &          0.956387200792748594152D-1,
     &          0.938443990808045654352D-1,
     &          0.911738786957638874594D-1,
     &          0.876520930044038163114D-1,
     &          0.833119242269464135758D-1,
     &          0.781938957870682705896D-1,
     &          0.723457941088502071174D-1,
     &          0.658222227761892440852D-1,
     &          0.586840934788972562042D-1,
     &          0.509980592613025930634D-1,
     &          0.428358980239763007216D-1,
     &          0.342738629139351946052D-1,
     &          0.253920653084367920322D-1,
     &          0.162743947308230037851D-1,
     &          0.701861000943074980416D-2/
*
      ZETA2=PI**2/6.0D0
      ZETA3=1.20205690315959428540D+0
*     PHI=3.0D0/4.0D0*PI
      PHI=9.0D0/10.0D0*PI
*
      RETURN
      END
      SUBROUTINE DEFAUL
C     -----------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  DEFAULT PARAMETERS
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /EP/EPS
      COMMON/RUN / IRUN
      COMMON/TEST/ ITEST1,ITEST2,ITEST3
      COMMON/MOMPA/ NMIN,NMAX
      COMMON/FUNPA/ IMIN,IMAX
      COMMON/IAPP  / IAPP
C
      EPS = 1.0D-9
      IRUN = 1
      ITEST1=1
      ITEST2=1
      ITEST3=1
      IMIN=1
      IMAX=26
      NMIN=1
      NMAX=20
      IAPP=1
C
      RETURN
      END
      SUBROUTINE WROUT
C     ----------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  WRITE OUT OF THE DEFAULT PARAMETERS
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /EP/EPS
      COMMON/TEST/ ITEST1,ITEST2,ITEST3
      COMMON/MOMPA/ NMIN,NMAX
      COMMON/FUNPA/ IMIN,IMAX
      COMMON/IAPP  / IAPP
      COMMON/RUN / IRUN
C
      WRITE(6,*) '  '
      WRITE(6,*) '*** INPUT PARAMETERS :'
      WRITE(6,*) '  '
      WRITE(6,*) '*** IRUN   = ', IRUN
      WRITE(6,*) '*** ITEST1 = ', ITEST1
      WRITE(6,*) '*** ITEST2 = ', ITEST2
      WRITE(6,*) '*** ITEST3 = ', ITEST3
      WRITE(6,*) '*** IMIN   = ', IMIN
      WRITE(6,*) '*** IMAX   = ', IMAX
      WRITE(6,*) '*** NMIN   = ', NMIN
      WRITE(6,*) '*** NMAX   = ', NMAX
      WRITE(6,*) '*** IAPP   = ', IAPP
      WRITE(6,*) '*** EPS    = ', EPS
      WRITE(6,*) '        '
      WRITE(6,*) '        '
C
      RETURN
      END
      SUBROUTINE ACRUN
C     ----------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  RUNNING THE CODE
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMPLEX*16 ACG1,ACG2,ACG3,ACG4,ACG5,ACG6,ACG7,ACG8,ACG9,ACG10,
     &           ACG11,ACG12,ACG13,ACG14,ACG15,ACG16,ACG17,ACG18,ACG19,
     &           ACG20,ACG21,ACG22,ACG23,ACG24,ACG25,ACG26,ZZ,R
C
      COMMON/TEST/ ITEST1,ITEST2,ITEST3
      COMMON/MOMPA/ NMIN,NMAX
      COMMON/FUNPA/ IMIN,IMAX
      COMMON/IAPP  / IAPP
      COMMON/RUN / IRUN
C
      DIMENSION XP(26)
C
      DATA XP/1.0D-7,1.0D-6,1.0D-5,1.0D-4,1.0D-3,1.0D-2,5.0D-2,
     &        1.0D-1,1.5D-1,2.0D-1,2.5D-1,3.0D-1,3.5D-1,4.0D-1,
     &        4.5D-1,5.0D-1,5.5D-1,6.0D-1,6.5D-1,7.0D-1,7.5D-1,
     &        8.0D-1,8.5D-1,9.0D-1,9.5D-1,9.9D-1/
C
C---  RUNNING THE CODE
C
      I1=0
      I2=0
      I3=0
      IF(IRUN.LT.0.OR.IRUN.GT.5) GOTO  3000
      IF(IRUN.EQ.0.AND.ITEST1.EQ.1) I1=1
      IF(IRUN.EQ.0.AND.ITEST2.EQ.1) I2=1
      IF(IRUN.EQ.0.AND.ITEST3.EQ.1) I3=1
      IF(I1.NE.1.AND.I2.NE.1.AND.I3.NE.1) GOTO  3001
      GOTO 3002
3000  WRITE(6,*) '***** IRUN=',IRUN,' ** OUT OF RANGE, STOP ***'
      STOP
3002  CONTINUE
      IF(ITEST1.EQ.1) GOTO 101
1001  IF(ITEST2.EQ.1) GOTO 102
1002  IF(ITEST3.EQ.1) GOTO 103
      GOTO 1000
C
C---  TEST 1: POSITIVE INTEGER MOMENTS: NUM. REPRESENTATION VS HARMONIC
C                                                              SUMS
C
101   CONTINUE
C
      DO 1 I=IMIN,IMAX
      WRITE(6,*) '****************************************'
      DO 2 N=NMIN,NMAX
C
      CALL MOMTES(I,N,RES,F1)
C
      WRITE(6,2003) I,N,RES,F1
2003  FORMAT(1X,'TEST1: I,K,R=',I4,1X,I4,1X,E12.3,1X,E12.5)
C
2     CONTINUE
      WRITE(6,*) '****************************************'
1     CONTINUE
C
      GOTO 1001
C
102   CONTINUE
C
C---  TEST 2: POSITIVE INTEGER MOMENTS: HARMONIC SUMS VS COMPLEX MOMENT
C                                                        REPRESENTATION
C
      DO 3 I=IMIN,IMAX
      WRITE(6,*) '****************************************'
      DO 4 N=NMIN,NMAX
C
      CALL ACREL(I,N,RES)
C
      IF(I.NE.3) THEN
      WRITE(6,2001) I,N,RES
2001  FORMAT(1X,'TEST2: I,N,R=',I4,1X,I4,1X,E12.3)
      ELSE
      WRITE(6,2002) IAPP,I,N,RES
2002  FORMAT(1X,'TEST2: IAPP,I,K,R=',I4,1X,I4,1X,I4,1X,D12.3)
      ENDIF
C
4     CONTINUE
      WRITE(6,*) '****************************************'
3     CONTINUE
C
      GOTO 1002
C
103   CONTINUE
C
C---  TEST 3: X DEPENDENCE OF BASIC FUNCTIONS: NUM. RESULTS AGAINST
C                                                   MELLIN INVERSION
C
*
      DO 5 I=IMIN,IMAX
      WRITE(6,*)
     &'**********************************************************'
      WRITE(6,*) '       '
      WRITE(6,*) '*** I = ',I,'  ***'
      WRITE(6,*) '       '
      DO 6 J=1,26
C
      X=XP(J)
      CALL INVERS(X,I,RES,F1)
C
      WRITE(6,9000) X,RES,F1
9000  FORMAT('TEST3:X,RAT,VAL=',D12.5,3X,D12.5,3X,D12.5)
C
6     CONTINUE
      WRITE(6,*)
     &'**********************************************************'
5     CONTINUE
      GOTO 1000
3001  CONTINUE
C
C >>> IRUN=1,2,3,4,5
C
      GOTO(4001,4002,4003,4004,4005), IRUN
C
4001  CONTINUE
C
      DO 5001 II=IMIN,IMAX
      DO 5002 NN=NMIN,NMAX
      IF(II.EQ.1) T=FCT1(NN)
      IF(II.EQ.2) T=FCT2(NN)
      IF(II.EQ.3) T=FCT3(NN)
      IF(II.EQ.4) T=FCT4(NN)
      IF(II.EQ.5) T=FCT5(NN)
      IF(II.EQ.6) T=FCT6(NN)
      IF(II.EQ.7) T=FCT7(NN)
      IF(II.EQ.8) T=FCT8(NN)
      IF(II.EQ.9) T=FCT9(NN)
      IF(II.EQ.10) T=FCT10(NN)
      IF(II.EQ.11) T=FCT11(NN)
      IF(II.EQ.12) T=FCT12(NN)
      IF(II.EQ.13) T=FCT13(NN)
      IF(II.EQ.14) T=FCT14(NN)
      IF(II.EQ.15) T=FCT15(NN)
      IF(II.EQ.16) T=FCT16(NN)
      IF(II.EQ.17) T=FCT17(NN)
      IF(II.EQ.18) T=FCT18(NN)
      IF(II.EQ.19) T=FCT19(NN)
      IF(II.EQ.20) T=FCT20(NN)
      IF(II.EQ.21) T=FCT21(NN)
      IF(II.EQ.22) T=FCT22(NN)
      IF(II.EQ.23) T=FCT23(NN)
      IF(II.EQ.24) T=FCT24(NN)
      IF(II.EQ.25) T=FCT25(NN)
      IF(II.EQ.26) T=FCT26(NN)
C
      WRITE(6,7001) II,NN,T
7001  FORMAT(1X,'IRUN=1 *** I,N, MOM=',I4,1X,I4,1X,E12.5)
5002  CONTINUE
5001  CONTINUE
      GOTO 1000
4002  CONTINUE
C
      DO 5003 II=IMIN,IMAX
      DO 5004 NN=NMIN,NMAX
      IF(II.EQ.1) T=XCG1(NN)
      IF(II.EQ.2) T=XCG2(NN)
      IF(II.EQ.3) T=XCG3(NN)
      IF(II.EQ.4) T=XCG4(NN)
      IF(II.EQ.5) T=XCG5(NN)
      IF(II.EQ.6) T=XCG6(NN)
      IF(II.EQ.7) T=XCG7(NN)
      IF(II.EQ.8) T=XCG8(NN)
      IF(II.EQ.9) T=XCG9(NN)
      IF(II.EQ.10) T=XCG10(NN)
      IF(II.EQ.11) T=XCG11(NN)
      IF(II.EQ.12) T=XCG12(NN)
      IF(II.EQ.13) T=XCG13(NN)
      IF(II.EQ.14) T=XCG14(NN)
      IF(II.EQ.15) T=XCG15(NN)
      IF(II.EQ.16) T=XCG16(NN)
      IF(II.EQ.17) T=XCG17(NN)
      IF(II.EQ.18) T=XCG18(NN)
      IF(II.EQ.19) T=XCG19(NN)
      IF(II.EQ.20) T=XCG20(NN)
      IF(II.EQ.21) T=XCG21(NN)
      IF(II.EQ.22) T=XCG22(NN)
      IF(II.EQ.23) T=XCG23(NN)
      IF(II.EQ.24) T=XCG24(NN)
      IF(II.EQ.25) T=XCG25(NN)
      IF(II.EQ.26) T=XCG26(NN)
C
      WRITE(6,7002) II,NN,T
7002  FORMAT(1X,'IRUN=2 *** I,N, MOM=',I4,1X,I4,1X,E12.5)
5004  CONTINUE
5003  CONTINUE
      GOTO 1000
4003  CONTINUE
C
      DO 5005 II=IMIN,IMAX
      DO 5006 NN=NMIN,NMAX
      R=DCMPLX(DBLE(NN),0.0D0)
      IF(II.EQ.1) ZZ=ACG1(R)
      IF(II.EQ.2) ZZ=ACG2(R)
      IF(II.EQ.3) ZZ=ACG3(R)
      IF(II.EQ.4) ZZ=ACG4(R)
      IF(II.EQ.5) ZZ=ACG5(R)
      IF(II.EQ.6) ZZ=ACG6(R)
      IF(II.EQ.7) ZZ=ACG7(R)
      IF(II.EQ.8) ZZ=ACG8(R)
      IF(II.EQ.9) ZZ=ACG9(R)
      IF(II.EQ.10) ZZ=ACG10(R)
      IF(II.EQ.11) ZZ=ACG11(R)
      IF(II.EQ.12) ZZ=ACG12(R)
      IF(II.EQ.13) ZZ=ACG13(R)
      IF(II.EQ.14) ZZ=ACG14(R)
      IF(II.EQ.15) ZZ=ACG15(R)
      IF(II.EQ.16) ZZ=ACG16(R)
      IF(II.EQ.17) ZZ=ACG17(R)
      IF(II.EQ.18) ZZ=ACG18(R)
      IF(II.EQ.19) ZZ=ACG19(R)
      IF(II.EQ.20) ZZ=ACG20(R)
      IF(II.EQ.21) ZZ=ACG21(R)
      IF(II.EQ.22) ZZ=ACG22(R)
      IF(II.EQ.23) ZZ=ACG23(R)
      IF(II.EQ.24) ZZ=ACG24(R)
      IF(II.EQ.25) ZZ=ACG25(R)
      IF(II.EQ.26) ZZ=ACG26(R)
C
      WRITE(6,*) 'IRUN=3 ***, I,N, MOM=',II,R,ZZ
5006  CONTINUE
5005  CONTINUE
      GOTO 1000
C
4004  CONTINUE
      DO 5007 II=IMIN,IMAX
      DO 5008 IX=1,26
C
      X=XP(IX)
      IF(II.EQ.1) T=FKN1(X)
      IF(II.EQ.2) T=FKN2(X)
      IF(II.EQ.3) T=FKN3(X)
      IF(II.EQ.4) T=FKN4(X)
      IF(II.EQ.5) T=FKN5(X)
      IF(II.EQ.6) T=FKN6(X)
      IF(II.EQ.7) T=FKN7(X)
      IF(II.EQ.8) T=FKN8(X)
      IF(II.EQ.9) T=FKN9(X)
      IF(II.EQ.10) T=FKN10(X)
      IF(II.EQ.11) T=FKN11(X)
      IF(II.EQ.12) T=FKN12(X)
      IF(II.EQ.13) T=FKN13(X)
      IF(II.EQ.14) T=FKN14(X)
      IF(II.EQ.15) T=FKN15(X)
      IF(II.EQ.16) T=FKN16(X)
      IF(II.EQ.17) T=FKN17(X)
      IF(II.EQ.18) T=FKN18(X)
      IF(II.EQ.19) T=FKN19(X)
      IF(II.EQ.20) T=FKN20(X)
      IF(II.EQ.21) T=FKN21(X)
      IF(II.EQ.22) T=FKN22(X)
      IF(II.EQ.23) T=FKN23(X)
      IF(II.EQ.24) T=FKN24(X)
      IF(II.EQ.25) T=FKN25(X)
      IF(II.EQ.26) T=FKN26(X)
C
      WRITE(6,*) II,X,T
7004  FORMAT(1X,'IRUN=4 *** I,X, FUN=',I4,1X,E12.3,1X,E12.5)
5008  CONTINUE
5007  CONTINUE
      GOTO 1000
C
4005  CONTINUE
C
      DO 5009 II=IMIN,IMAX
      DO 5010 IX=1,26
C
      X=XP(IX)
      CALL INV1(X,II,F2)
C
      WRITE(6,*) II,X,F2
7005  FORMAT(1X,'IRUN=5 *** I,X, FUN=',I4,1X,E12.3,1X,E12.5)
5010  CONTINUE
5009  CONTINUE
C
1000  CONTINUE
C
C
      CALL URUN
C
      RETURN
      END
      SUBROUTINE MOMTES(I,N,RES,F1)
C     -----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  CALCULATE THE RELATIVE ACCURACY OF THE Nth MOMENT COMPARING
C---  THE REPRESENTATION THROUGH A NUMERICAL INTEGRAL FCTi WITH THAT
C---  BY (ALTERNATING) NESTED HARMONIC SUMS  XCGi
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON/IAPP  / IAPP
      COMMON/VAL   / ZERO,ONE

C
      EXTERNAL FCT1,FCT2,FCT3,FCT4,FCT5,FCT6,FCT7,FCT8,FCT9,FCT10,
     &         FCT11,FCT12,FCT13,FCT14,FCT15,FCT16,FCT17,FCT18,FCT19,
     &         FCT20,FCT21,FCT22,FCT23,FCT24,FCT25,FCT26
      EXTERNAL XCG1,XCG2,XCG3,XCG4,XCG5,XCG6,XCG7,XCG8,XCG9,XCG10,
     &         XCG11,XCG12,XCG13,XCG14,XCG15,XCG16,XCG17,XCG18,XCG19,
     &         XCG20,XCG21,XCG22,XCG23,XCG24,XCG25,XCG26
C
      ZN=DCMPLX(DBLE(N),ZERO)
C
      GOTO(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     &     23,24,25,26), I
C
      IF(I.LT.1.OR.I.GT.26) 
     &WRITE(6,*) '*** FLAG I WRONG IN SR MOMTES,I=',I,'  *** STOP ***'
      IF(I.LT.1.OR.I.GT.26) STOP
C
1     F1=XCG1(N)
      RES=F1/FCT1(N)-ONE
      RETURN
2     F1=XCG2(N)
      RES=F1/FCT2(N)-ONE
      RETURN
3     F1=XCG3(N)
      RES=F1/FCT3(N)-ONE
      RETURN
4     F1=XCG4(N)
      RES=F1/FCT4(N)-ONE
      RETURN
5     F1=XCG5(N)
      RES=F1/FCT5(N)-ONE
      RETURN
6     F1=XCG6(N)
      RES=F1/FCT6(N)-ONE
      RETURN
7     F1=XCG7(N)
      RES=F1/FCT7(N)-ONE
      RETURN
8     F1=XCG8(N)
      RES=F1/FCT8(N)-ONE
      RETURN
9     F1=XCG9(N)
      RES=F1/FCT9(N)-ONE
      RETURN
10    F1=XCG10(N)
      RES=F1/FCT10(N)-ONE
      RETURN
11    F1=XCG11(N)
      RES=F1/FCT11(N)-ONE
      RETURN
12    F1=XCG12(N)
      RES=F1/FCT12(N)-ONE
      RETURN
13    F1=XCG13(N)
      RES=F1/FCT13(N)-ONE
      RETURN
14    F1=XCG14(N)
      RES=F1/FCT14(N)-ONE
      RETURN
15    F1=XCG15(N)
      RES=F1/FCT15(N)-ONE
      RETURN
16    F1=XCG16(N)
      RES=F1/FCT16(N)-ONE
      RETURN
17    F1=XCG17(N)
      RES=F1/FCT17(N)-ONE
      RETURN
18    F1=XCG18(N)
      RES=F1/FCT18(N)-ONE
      RETURN
19    F1=XCG19(N)
      RES=F1/FCT19(N)-ONE
      RETURN
20    F1=XCG20(N)
      RES=F1/FCT20(N)-ONE
      RETURN
21    F1=XCG21(N)
      RES=F1/FCT21(N)-ONE
      RETURN
22    F1=XCG22(N)
      RES=F1/FCT22(N)-ONE
      RETURN
23    F1=XCG23(N)
      RES=F1/FCT23(N)-ONE
      RETURN
24    F1=XCG24(N)
      RES=F1/FCT24(N)-ONE
      RETURN
25    F1=XCG25(N)
      RES=F1/FCT25(N)-ONE
      RETURN
26    F1=XCG26(N)
      RES=F1/FCT26(N)-ONE
      RETURN
      END
      SUBROUTINE ACREL(I,N,RES)
C     -------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  CALCULATE THE RELATIVE ACCURACY OF THE Nth MOMENT COMPARING
C---  THE REPRESENTATION FOR THE COMPLEX MELLIN MOMENTS ACGi
C---  WITH THOSE DUE TO FINITE (ALTERNATING) HARMONIC SUMS XCGi
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 ZERO,ONE,XCG1,XCG2,XCG3,XCG4,XCG5,XCG6
      REAL*8 XCG7,XCG8,XCG9,XCG10,XCG11,XCG12
      REAL*8 XCG13,XCG14,XCG15,XCG16,XCG17,XCG18
      REAL*8 XCG19,XCG20,XCG21,XCG22,XCG23,XCG24
      REAL*8 XCG25,XCG26
C
      REAL*8      RES
      INTEGER           I,IAPP
      COMMON/IAPP  / IAPP
      COMMON/VAL   / ZERO,ONE

C
      EXTERNAL ACG1,ACG2,ACG3,ACG4,ACG5,ACG6,ACG7,ACG8,ACG9,ACG10,
     &         ACG11,ACG12,ACG13,ACG14,ACG15,ACG16,ACG17,ACG18,ACG19,
     &         ACG20,ACG21,ACG22,ACG23,ACG24,ACG25,ACG26
      EXTERNAL XCG1,XCG2,XCG3,XCG4,XCG5,XCG6,XCG7,XCG8,XCG9,XCG10,
     &         XCG11,XCG12,XCG13,XCG14,XCG15,XCG16,XCG17,XCG18,XCG19,
     &         XCG20,XCG21,XCG22,XCG23,XCG24,XCG25,XCG26
C
      ZN=DCMPLX(DBLE(N),ZERO)
C
      GOTO(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     &     23,24,25,26), I
C
      IF(I.LT.1.OR.I.GT.26) 
     &WRITE(6,*) '*** FLAG I WRONG IN SR ACREL, I=',I,'  *** STOP ***'
      IF(I.LT.1.OR.I.GT.26) STOP
C
1     RES=DREAL(ACG1(ZN)/XCG1(N))-ONE
      RETURN
2     RES=DREAL(ACG2(ZN)/XCG2(N))-ONE
      RETURN
3     RES=DREAL(ACG3(ZN)/XCG3(N))-ONE
      RETURN
4     RES=DREAL(ACG4(ZN)/XCG4(N))-ONE
      RETURN
5     RES=DREAL(ACG5(ZN)/XCG5(N))-ONE
      RETURN
6     RES=DREAL(ACG6(ZN)/XCG6(N))-ONE
      RETURN
7     RES=DREAL(ACG7(ZN)/XCG7(N))-ONE
      RETURN
8     RES=DREAL(ACG8(ZN)/XCG8(N))-ONE
      RETURN
9     RES=DREAL(ACG9(ZN)/XCG9(N))-ONE
      RETURN
10    RES=DREAL(ACG10(ZN)/XCG10(N))-ONE
      RETURN
11    RES=DREAL(ACG11(ZN)/XCG11(N))-ONE
      RETURN
12    RES=DREAL(ACG12(ZN)/XCG12(N))-ONE
      RETURN
13    RES=DREAL(ACG13(ZN)/XCG13(N))-ONE
      RETURN
14    RES=DREAL(ACG14(ZN)/XCG14(N))-ONE
      RETURN
15    RES=DREAL(ACG15(ZN)/XCG15(N))-ONE
      RETURN
16    RES=DREAL(ACG16(ZN)/XCG16(N))-ONE
      RETURN
17    RES=DREAL(ACG17(ZN)/XCG17(N))-ONE
      RETURN
18    RES=DREAL(ACG18(ZN)/XCG18(N))-ONE
      RETURN
19    RES=DREAL(ACG19(ZN)/XCG19(N))-ONE
      RETURN
20    RES=DREAL(ACG20(ZN)/XCG20(N))-ONE
      RETURN
21    RES=DREAL(ACG21(ZN)/XCG21(N))-ONE
      RETURN
22    RES=DREAL(ACG22(ZN)/XCG22(N))-ONE
      RETURN
23    RES=DREAL(ACG23(ZN)/XCG23(N))-ONE
      RETURN
24    RES=DREAL(ACG24(ZN)/XCG24(N))-ONE
      RETURN
25    RES=DREAL(ACG25(ZN)/XCG25(N))-ONE
      RETURN
26    CONTINUE
      RES=DREAL(ACG26(ZN)/XCG26(N))-ONE
      RETURN
      END
      SUBROUTINE INVERS(X,I,RES,F1)
*     -----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  CALCULATE THE RELATIVE ACCURACY OF THE MELLIN INVERSION USING
C---  THE ANALYTIC CONTINUATION ACGi VS THE NUMERICAL REPRESENTATION
C---  OF THE BASUC FUNCTIONS FKNi
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 DIR,STA,SS,ZN,FF,FUNC
*
      COMMON /OFFS  / CP(26)
      COMMON /PI    / PI
      COMMON /WGAUSS/ XL8(8),WL8(8),XL32(32),WL32(32)
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON /FKNLL/  LL
      COMMON /ZET2/ ZETA2
      COMMON /ZET3/ ZETA3
*
      DIMENSION ZINT(51)
*
      LL = I
      IGAUSS = 8
      ZERO=0.0D0
*
      CPAR=CP(I)
      IF(X.LE.0.7D0) RANGE=71.0D0
      IF(X.LE.0.7D0) RANGE=71.0D0
      IF(X.LE.0.7D0) INTE =20
      IF(X.GT.0.7D0) INTE =20
      IF(X.LE.0.7D0) PHI=3.0D0/4.0D0*PI
      IF(X.GT.0.7D0) PHI=7.0D0/8.0D0*PI
      IF(X.GT.0.98D0) PHI=19.0D0/20.0D0*PI
      DIR=DCMPLX(DCOS(PHI),DSIN(PHI))
      STA=DCMPLX(CPAR,ZERO)
      ST   =EXP(LOG(RANGE)/DBLE(INTE))
      ST0  =1.0D0
      ZINT(1)=ZERO
      DO 10 L=1,20
      ST0 =ST0*ST
      L1=L+1
      ZINT(L1)=ST0-1.0D0
10    CONTINUE
C
      SS=DCMPLX(ZERO,ZERO)
      XBL=LOG(X)
*
      DO 20 L=1,20
      A=ZINT(L)
      B=ZINT(L+1)
      DLI=(B-A)/2.0D0
      SLI=(B+A)/2.0D0
*
      IF(IGAUSS.EQ. 8) NM=8
      IF(IGAUSS.EQ.32) NM=32
      IF(IGAUSS.NE.8.AND.IGAUSS.NE.32) 
     &WRITE(6,*) '*** IGAUSS=',IGAUSS,'*** WRONG IN SR INVERSE, STOP'
      IF(IGAUSS.NE.8.AND.IGAUSS.NE.32) STOP
*
      DO 21 M=1,NM
      IF(IGAUSS.EQ. 8) W=WL8(M)
      IF(IGAUSS.EQ.32) W=WL32(M)
      IF(IGAUSS.EQ. 8) ZZ=XL8(M)
      IF(IGAUSS.EQ.32) ZZ=XL32(M)
      Z= DLI*ZZ+SLI
      ZN=STA+Z*DIR
      FF=FUNC(ZN)
      SS=SS+FF*DIR*EXP(-XBL*ZN)*W*DLI
21    CONTINUE
*
20    CONTINUE
*
      F=DIMAG(SS)/PI
      F1=FUNO(X)
      RES=F/F1-1.0D0
1     CONTINUE
*
100   CONTINUE
*
      RETURN
      END
      SUBROUTINE INV1(X,I,F2)
*     -----------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  CALCULATE THE MELLIN INVERSION USING
C---  THE ANALYTIC CONTINUATION ACGi
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 DIR,STA,SS,ZN,FF,FUNC
*
      COMMON /OFFS  / CP(26)
      COMMON /PI    / PI
      COMMON /WGAUSS/ XL8(8),WL8(8),XL32(32),WL32(32)
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON /FKNLL/  LL
      COMMON /ZET2/ ZETA2
      COMMON /ZET3/ ZETA3
*
      DIMENSION ZINT(51)
*
      LL = I
      IGAUSS = 8
      ZERO=0.0D0
*
      CPAR=CP(I)
      IF(X.LE.0.7D0) RANGE=71.0D0
      IF(X.LE.0.7D0) RANGE=71.0D0
      IF(X.LE.0.7D0) INTE =20
      IF(X.GT.0.7D0) INTE =20
      IF(X.LE.0.7D0) PHI=3.0D0/4.0D0*PI
      IF(X.GT.0.7D0) PHI=7.0D0/8.0D0*PI
      IF(X.GT.0.98D0) PHI=19.0D0/20.0D0*PI
      DIR=DCMPLX(DCOS(PHI),DSIN(PHI))
      STA=DCMPLX(CPAR,ZERO)
      ST   =EXP(LOG(RANGE)/DBLE(INTE))
      ST0  =1.0D0
      ZINT(1)=ZERO
      DO 10 L=1,20
      ST0 =ST0*ST
      L1=L+1
      ZINT(L1)=ST0-1.0D0
10    CONTINUE
C
      SS=DCMPLX(ZERO,ZERO)
      XBL=LOG(X)
*
      DO 20 L=1,20
      A=ZINT(L)
      B=ZINT(L+1)
      DLI=(B-A)/2.0D0
      SLI=(B+A)/2.0D0
*
      IF(IGAUSS.EQ. 8) NM=8
      IF(IGAUSS.EQ.32) NM=32
      IF(IGAUSS.NE.8.AND.IGAUSS.NE.32) 
     &WRITE(6,*) '*** IGAUSS=',IGAUSS,'*** WRONG IN SR INVERSE, STOP'
      IF(IGAUSS.NE.8.AND.IGAUSS.NE.32) STOP
*
      DO 21 M=1,NM
      IF(IGAUSS.EQ. 8) W=WL8(M)
      IF(IGAUSS.EQ.32) W=WL32(M)
      IF(IGAUSS.EQ. 8) ZZ=XL8(M)
      IF(IGAUSS.EQ.32) ZZ=XL32(M)
      Z= DLI*ZZ+SLI
      ZN=STA+Z*DIR
      FF=FUNC(ZN)
      SS=SS+FF*DIR*EXP(-XBL*ZN)*W*DLI
21    CONTINUE
*
20    CONTINUE
*
      F=DIMAG(SS)/PI
      F2=F
1     CONTINUE
*
100   CONTINUE
*
      RETURN
      END
      REAL*8 FUNCTION FUNO(X)
*     -----------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  SELECTION ROUTINE FOR: FKNi
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /FKNLL/  LL
*
      GOTO(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     &     23,24,25,26), LL
1     F=FKN1(X)
      GOTO 1000
2     F=FKN2(X)
      GOTO 1000
3     F=FKN3(X)
      GOTO 1000
4     F=FKN4(X)
      GOTO 1000
5     F=FKN5(X)
      GOTO 1000
6     F=FKN6(X)
      GOTO 1000
7     F=FKN7(X)
      GOTO 1000
8     F=FKN8(X)
      GOTO 1000
9     F=FKN9(X)
      GOTO 1000
10    F=FKN10(X)
      GOTO 1000
11    F=FKN11(X)
      GOTO 1000
12    F=FKN12(X)
      GOTO 1000
13    F=FKN13(X)
      GOTO 1000
14    F=FKN14(X)
      GOTO 1000
15    F=FKN15(X)
      GOTO 1000
16    F=FKN16(X)
      GOTO 1000
17    F=FKN17(X)
      GOTO 1000
18    F=FKN18(X)
      GOTO 1000
19    F=FKN19(X)
      GOTO 1000
20    F=FKN20(X)
      GOTO 1000
21    F=FKN21(X)
      GOTO 1000
22    F=FKN22(X)
      GOTO 1000
23    F=FKN23(X)
      GOTO 1000
24    F=FKN24(X)
      GOTO 1000
25    F=FKN25(X)
      GOTO 1000
26    F=FKN26(X)
*
1000  FUNO=F
*
      RETURN
      END
      COMPLEX*16 FUNCTION FUNC(Z)
*     ---------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  SELECTION ROUTINE FOR: ACGi
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      COMMON /FKNLL/  LL
*
      Z1=Z-1.0D0
*
      GOTO(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     &     23,24,25,26), LL
1     F=ACG1(Z1)
      GOTO 1000
2     F=ACG2(Z1)
      GOTO 1000
3     F=ACG3(Z1)
      GOTO 1000
4     F=ACG4(Z1)
      GOTO 1000
5     F=ACG5(Z1)
      GOTO 1000
6     F=ACG6(Z1)
      GOTO 1000
7     F=ACG7(Z1)
      GOTO 1000
8     F=ACG8(Z1)
      GOTO 1000
9     F=ACG9(Z1)
      GOTO 1000
10    F=ACG10(Z1)
      GOTO 1000
11    F=ACG11(Z1)
      GOTO 1000
12    F=ACG12(Z1)
      GOTO 1000
13    F=ACG13(Z1)
      GOTO 1000
14    F=ACG14(Z1)
      GOTO 1000
15    F=ACG15(Z1)
      GOTO 1000
16    F=ACG16(Z1)
      GOTO 1000
17    F=ACG17(Z1)
      GOTO 1000
18    F=ACG18(Z1)
      GOTO 1000
19    F=ACG19(Z1)
      GOTO 1000
20    F=ACG20(Z1)
      GOTO 1000
21    F=ACG21(Z1)
      GOTO 1000
22    F=ACG22(Z1)
      GOTO 1000
23    F=ACG23(Z1)
      GOTO 1000
24    F=ACG24(Z1)
      GOTO 1000
25    F=ACG25(Z1)
      GOTO 1000
26    F=ACG26(Z1)
*
1000  FUNC=F
*
      RETURN
      END
      SUBROUTINE ACEND
C     ----------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C     END OF THE CODE: ANCONT
C
C************************************************************************
C
C
      CALL UOUT
C
      WRITE(6,*) '***** ANCONT COMPLETED *****'
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
      COMMON/ACLOG2/ AK2(10)
      COMMON/ACCON1/ DL
      COMMON/VAL   / ZERO,ONE
C
      T=DCMPLX(ZERO,ZERO)
      DO 1 L=2,11
      K=L-1
      T=T+AK2(K)/(ZN+DBLE(K+1))
1     CONTINUE
C
      ACG1=(DL*DL- ZN*T)/2.0D0
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
      COMMON/ACLOG3/ AK3(11)
      COMMON/ACCON1/ DL
      COMMON/VAL   / ZERO,ONE
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
      COMMON/ACLOG1/ AK1(9)
      COMMON/IAPP  / IAPP
      COMMON/VAL   / ZERO,ONE
      COMMON/ACCON1/ DL
      COMMON/ACCON2/ ZET2
      COMMON/ACCON4/ GE
C
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
      DO 1 L=1,9
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
      COMMON/ACLOG1/ AK1(9)
      COMMON/ACCON1/ DL
      COMMON/ACCON2/ ZET2
      COMMON/VAL   / ZERO,ONE
C
      T=DCMPLX(-ZET2/2.0D0*DL,ZERO)
C
      DO 1 L=1,9
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
      COMMON/ACLOG1/ AK1(9)
      COMMON/VAL   / ZERO,ONE
      COMMON/ACCON1/ DL
      COMMON/ACCON2/ ZET2
      COMMON/ACCON4/ GE
C
      T=DCMPLX(ZERO,ZERO)
      Z=ZN
      DO 1 L=1,9
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
      COMMON/ACLOG1/ AK1(9)
      COMMON/ACCON1/ DL
      COMMON/ACCON2/ ZET2
      COMMON/ACCON3/ ZET3
      COMMON/ACCON4/ GE
      COMMON/VAL   / ZERO,ONE
C
      T=DCMPLX(DL*ZET3,ZERO)
      DO 1 L=1,9
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
      COMMON/ACLOG1/ AK1(9)
      COMMON/ACCON1/ DL
      COMMON/ACCON2/ ZET2
      COMMON/ACCON3/ ZET3
      COMMON/VAL   / ZERO,ONE
C
      T=DCMPLX(-3.0D0*ZET3/4.0D0*DL,ZERO)
C
      DO 1 L=1,9
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
      COMMON/ACLOG1/ AK1(9)
      COMMON/ACCON1/ DL
      COMMON/ACCON2/ ZET2
      COMMON/ACCON3/ ZET3
      COMMON/ACCON4/ GE
      COMMON/VAL   / ZERO,ONE
C
      T=DCMPLX(ZERO,ZERO)
      DO 1 L=1,9
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
      COMMON/ACLOG1/ AK1(9)
      COMMON/ACLOG2/ AK2(10)
      COMMON/ACLOG3/ AK3(11)
      COMMON/ACCON1/ DL
      COMMON/ACCON3/ ZET3
      COMMON/VAL   / ZERO,ONE
C
      T=DCMPLX(ZET3*DL/8.0D0,ZERO)
      DO 1 K=1,9
      T1=DCMPLX(ZERO,ZERO)
      DO 2 L=2,11
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
      COMPLEX*16 FUNCTION ACG10(ZN)
C     -----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: I1(X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 AK1,AK2,DL,ZERO,ONE,ZET3,GE
C
      COMMON/ACLOG1/ AK1(9)
      COMMON/ACLOG2/ AK2(10)
      COMMON/ACCON1/ DL
      COMMON/ACCON4/ GE
      COMMON/ACCON3/ ZET3
      COMMON/VAL   / ZERO,ONE
C
      T=DCMPLX(-5.0D0/8.0D0*ZET3*DL,ZERO)
C
      DO 1 K=2,11
      L=K-1
      ZNK=ZN+DBLE(K)
      ZNK1=ZNK+ONE
      CALL PSI0(ZNK1,PS)
      S1=PS+GE
      T=T+AK2(L)*S1/ZNK
1     CONTINUE
      DO 2 K=1,9
      ZNK=ZN+DBLE(K)
      T2=DCMPLX(ZERO,ZERO)
      DO 3 L=1,9
      ZNKL=ZNK+DBLE(L)
      ZNKL1=ZNKL+ONE
      CALL PSI0(ZNKL1,PS1)
      S2=PS1+GE
      T2=T2-AK1(L)*S2/ZNKL
3     CONTINUE
      T=T+AK1(K)*ZN/ZNK*(5.0D0/8.0D0*ZET3+T2)
2     CONTINUE
C
      ACG10=T
C
      RETURN
      END
      COMPLEX*16 FUNCTION ACG11(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(1-X)/(1+X) LI2(X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 AK1,CK1,P21,DL,ZERO,ONE,ZET2,GE,ZET3
C
      COMMON/ACLOG1/ AK1(9)
      COMMON/ACLOG6/ CK1(12)
      COMMON/POLY1 / P21(4)
      COMMON/ACCON1/ DL
      COMMON/ACCON2/ ZET2
      COMMON/ACCON3/ ZET3
      COMMON/ACCON4/ GE
      COMMON/VAL   / ZERO,ONE
C
      T=DCMPLX(ZERO,ZERO)
      T1=DCMPLX((-ZET2+DL**2)/2.0D0,ZERO)
      T2=DCMPLX(7.0D0/4.0D0*ZET3-ZET2*DL+DL**3/3.0D0,ZERO)
C
      DO 1 K=1,12
      L=K-1
      U2=T1
      ZK=DCMPLX(DBLE(L))
C
      DO 2 M=1,9
      ZM=DCMPLX(DBLE(M))
      ZM1=ZM+ONE
      ZNM=ZN+ZM+ZK
      ZNM1=ZNM+ONE
      CALL PSI0(ZNM1,PS)
      CALL PSI0(ZM1,PS1)
      S1=PS+GE
      S11=PS1+GE
      U2=U2-AK1(M)*(ZM/ZNM*S1-S11)
2     CONTINUE
C
      T=T+CK1(K)*U2
1     CONTINUE
C
      DO 3 K=1,4       
      L=K-1
C
      U3=T2
      ZK=DCMPLX(DBLE(L),ZERO)
      DO 4 M=1,9
      ZM=DCMPLX(DBLE(M),ZERO)
      ZM1=ZM+ONE
      ZMN=ZN+ZK+ZM
      ZMN1=ZMN+ONE
      CALL PSI0(ZMN1,PS5)
      CALL PSI1(ZMN1,PS2)
      CALL PSI0(ZM1,PS3)
      CALL PSI1(ZM1,PS4)
      S1=PS5+GE
      S2=-PS2+ZET2
      S11=PS3+GE
      S21=-PS4+ZET2
      U3=U3+AK1(M)*(ZM/ZMN*(S1**2+S2)-(S11**2+S21))
4     CONTINUE
      T=T+P21(K)*U3
3     CONTINUE
C
      ACG11=T
C
      RETURN
      END
      COMPLEX*16 FUNCTION ACG12(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(1-X)/(1+X)*LI2(-X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 AK1,AK2,DL,ZERO,ONE,GE
C
      COMMON/ACLOG1/ AK1(9)
      COMMON/ACLOG2/ AK2(10)
      COMMON/ACCON1/ DL
      COMMON/ACCON4/ GE
      COMMON/VAL   / ZERO,ONE
C
      T=DCMPLX(ZERO,ZERO)
      DO 1 K=1,9
      ZK=DCMPLX(DBLE(K),ZERO)
      T2=DCMPLX(ZERO,ZERO)
      ZNK=ZN+ZK
      DO 2 L=2,11
      L1=L-1
      ZNKL=ZNK+DBLE(L)
      T2=T2+AK2(L1)*ZNK/ZNKL
2     CONTINUE
      ZNK1=ZNK+ONE
      CALL BET(ZNK1,U1)
      CALL PSI0(ZNK1,U2)
      V1=U1*(U2+GE-DL)
      CALL BET1(ZNK1,V2)
C
      RES=-ONE/2.0D0*(DL**2-T2)+V2-V1
      T=T-AK1(K)/ZK*RES
1     CONTINUE
C
      ACG12=T
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
      COMMON/ACLOG2/ AK2(10)
      COMMON/ACLOG3/ AK3(11)
      COMMON/ACCON1/ DL
      COMMON/ACCON2/ ZET2
      COMMON/VAL   / ZERO,ONE
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
      DO 2 L=2,11
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
      COMPLEX*16 FUNCTION ACG14(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: (LOG^2(1+X)-LOG^2(2))/(X-1)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 BK2,ZERO,ONE
C
      COMMON/ACLOG5/ BK2(10)
      COMMON/VAL   / ZERO,ONE
C
      T=DCMPLX(ZERO,ZERO)
      DO 1 L=1,10
      ZNL=ZN+DBLE(L)
      T=T+BK2(L)/ZNL
1     CONTINUE
C
      ACG14=T
C
      RETURN
      END
      COMPLEX*16 FUNCTION ACG15(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: (LOG(1+x)-LOG(2))/(X-1)*LI2(X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 BK1,ZERO,ONE,ZET2,GE
C
      COMMON/ACLOG4/ BK1(9)
      COMMON/ACCON2/ ZET2
      COMMON/ACCON4/ GE
      COMMON/VAL   / ZERO,ONE
C
      T=ZERO
      DO 1 L=1,9
      ZNK=ZN+DCMPLX(DBLE(L),ZERO)
      ZNK2=ZNK+ONE
      CALL PSI0(ZNK2,PS)
      S1=PS+GE
      T=T+BK1(L)/ZNK*(ZET2-S1/ZNK)
1     CONTINUE
C
      ACG15= T
C
      RETURN
      END
      COMPLEX*16 FUNCTION ACG16(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: (LOG(1+x)-LOG(2))/(X-1)*LI2(-X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 BK1,ZERO,ONE,ZET2,DL
C
      COMMON/ACLOG4/ BK1(9)
      COMMON/ACCON2/ ZET2
      COMMON/ACCON1/ DL
      COMMON/VAL   / ZERO,ONE
C
      T=ZERO
      DO 1 L=1,9
      ZNK=ZN+DCMPLX(DBLE(L),ZERO)
      ZNK2=ZNK+ONE
      CALL BET(ZNK2,V1)
      T1=DL-V1
      T=T+BK1(L)/ZNK*(-ZET2/2.0D0+T1/ZNK)
1     CONTINUE
C
      ACG16= T
C
      RETURN
      END
      COMPLEX*16 FUNCTION ACG17(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(X) LOG^2(1+X)/(X-1)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 AK2,ZERO,ONE,ZET2
C
      COMMON/ACLOG2/ AK2(10)
      COMMON/VAL   / ZERO,ONE
      COMMON/ACCON2/ ZET2
C
      T=DCMPLX(ZERO,ZERO)
      DO 1 K=2,11
      L=K-1
      ZNK1=ZN+DBLE(K+1)
      CALL PSI1(ZNK1,V1)
      T=T+AK2(L)*V1
1     CONTINUE
C
      ACG17=T
C
      RETURN
      END
      COMPLEX*16 FUNCTION ACG18(ZN)
C     -----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: (LI2(X)-ZETA2)/(X-1)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 CK1,P21,ZERO,ONE,ZET2,ZET3,GE
C
      COMMON/ACLOG6/ CK1(12)
      COMMON/POLY1 / P21(4)
      COMMON/ACCON2/ ZET2
      COMMON/ACCON3/ ZET3
      COMMON/ACCON4/ GE
      COMMON/VAL   / ZERO,ONE
C
      ZN1=ZN+ONE
      T=DCMPLX(ZERO,ZERO)
C
      CALL PSI0(ZN1,PS1)
      CALL PSI1(ZN1,PS2)
      SPS1=PS1+GE
      SPS2=-PS2+ZET2
      T=T+(SPS1**2+SPS2)/ZN-ZET2*SPS1
C
      DO 1 L=1,12
      ZNK1=ZN+DBLE(L-1)
      ZNK2=ZNK1+ONE
      CALL PSI0(ZNK2,PS)
      S1=PS+GE
      T=T+CK1(L)*S1/ZNK1*ZN
1     CONTINUE
C
      DO 2 L=1,4
      ZNK1=ZN+DBLE(L-1)
      ZNK2=ZNK1+ONE
      CALL PSI0(ZNK2,PS)
      CALL PSI1(ZNK2,PS1)
      S1=PS+GE
      S2=-PS1+ZET2
      T=T-P21(L)*(S1**2+S2)/ZNK1*ZN
2     CONTINUE
C
      ACG18=T
C
      RETURN
      END
      COMPLEX*16 FUNCTION ACG19(ZN)
C     -----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: (LI2(-X)+ZETA2/2)/(X-1)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 AK1,ZERO,ONE,ZET2,GE
C
      COMMON/ACLOG1/ AK1(9)
      COMMON/ACCON2/ ZET2
      COMMON/ACCON4/ GE
      COMMON/VAL   / ZERO,ONE
C
      T=DCMPLX(ZERO,ZERO)
      ZN1=ZN+ONE
      CALL PSI0(ZN1,PS)
      S1=PS+GE
      T=T+S1*ZET2/2.0D0
      DO 1 L=1,9
      ZK=DCMPLX(DBLE(L))
      ZNK1=ZN+ZK+ONE
      CALL PSI0(ZNK1,PS)
      S1=PS+GE
      T=T-AK1(L)*S1/ZK
1     CONTINUE
C
      ACG19=T
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
      REAL*8 DL,ZERO,ONE,ZET2,ZET3,GE,P24,P34,P22,CK2,CK4
C
      COMMON/ACLOG7/ CK2(13)
      COMMON/ACLOG9/ CK4(13)
      COMMON/POLY2 / P22(4)
      COMMON/POLY4 / P24(5),P34(5)
      COMMON/ACCON1/ DL
      COMMON/ACCON2/ ZET2
      COMMON/ACCON3/ ZET3
      COMMON/ACCON4/ GE
      COMMON/VAL   / ZERO,ONE
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
      REAL*8 ZERO,ONE,ZET2,ZET3,GE,CK3,P23,P33
C
      COMMON/ACLOG8/ CK3(10)
      COMMON/POLY3 / P23(5),P33(5)
      COMMON/ACCON2/ ZET2
      COMMON/ACCON3/ ZET3
      COMMON/ACCON4/ GE
      COMMON/VAL   / ZERO,ONE
C
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
      COMPLEX*16 FUNCTION ACG22(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(X) LI2(X)/(X-1)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 CK1,P21,GE,ZERO,ONE
C
      COMMON/ACLOG6/ CK1(12)
      COMMON/POLY1 / P21(4)
      COMMON/ACCON4/ GE
      COMMON/VAL   / ZERO,ONE
C
      T=DCMPLX(ZERO,ZERO)
C
      DO 1 L=1,12
      ZNK1=ZN+DBLE(L)
      CALL PSI1(ZNK1,PS1)
      T=T+CK1(L)*PS1
1     CONTINUE
C
      DO 2 L=1,4
      ZNK1=ZN+DBLE(L)
      CALL PSI0(ZNK1,PS)
      CALL PSI1(ZNK1,PS1)
      CALL PSI2(ZNK1,PS2)
      TA=(PS+GE)*PS1-PS2/2.0D0
      T=T-P21(L)*TA
2     CONTINUE
C
      ACG22=T
C
      RETURN
      END
      COMPLEX*16 FUNCTION ACG23(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: (LI3(-X)+3/4*ZETA3)/(X-1)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 AK1,GE,ZERO,ONE,ZET2,ZET3,XLI4,DL
C
      COMMON/ACLOG1/ AK1(9)
      COMMON/ACCON1/ DL
      COMMON/ACCON2/ ZET2
      COMMON/ACCON3/ ZET3
      COMMON/ACCON4/ GE
      COMMON/ACCON5/ XLI4
      COMMON/VAL   / ZERO,ONE
C
C >>> M[Li3(x)/(1+x)](N)
C
      U=ACG6(ZN)
C
      T=DCMPLX(ZERO,ZERO)
      T=T-ZET2**2/2.0D0+ZET3*DL
      ZN1=ZN+ONE
      CALL PSI0(ZN1,PS)
      S1=PS+GE
      T=T+3.0D0/4.0D0*ZET3*S1-U
      DO 1 L=1,9
      ZK=DCMPLX(DBLE(L),ZERO)
      ZNK=ZN+ZK
      ZNK1=ZNK+ONE
      CALL PSI0(ZNK1,PS)
      S1=PS+GE
      R=ZN/ZNK
      H=ZET3-ZET2/ZNK-ZET2/ZK+
     &   S1*(1.0D0/ZNK**2+1.0D0/ZK**2+1.0D0/ZNK/ZK)
      T=T-AK1(L)*R*H
1     CONTINUE
C
      ACG23=T
C
      RETURN
      END
      COMPLEX*16 FUNCTION ACG24(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: (I1(X)+5/8*ZETA3)/(X-1)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 AK1,CK5,DK5,DL,ZERO,ONE,ZET3,GE,ZET2
C
      COMMON/ACLOG1/ AK1(9)
      COMMON/ACLG10/ CK5(9)
      COMMON/ACLG11/ DK5(10)
      COMMON/ACCON1/ DL
      COMMON/ACCON2/ ZET2
      COMMON/ACCON3/ ZET3
      COMMON/ACCON4/ GE
      COMMON/VAL   / ZERO,ONE
C
      T=DCMPLX(-2.0D0*ZET3*DL,ZERO)
C
C >>> M[S12(x)/(1+x)](N)
C
      U=ACG8(ZN)
C
      T=T+2.0D0*U
      ZN1=ZN+ONE
      CALL PSI0(ZN1,PS)
      S1=PS+GE
      T=T+5.0D0/8.0D0*ZET3*S1
C
      DO 1 L=1,9
      ZNK=ZN+DBLE(L)
      ZNK1=ZNK+ONE
      CALL PSI0(ZNK1,PS)
      CALL PSI1(ZNK1,PS1)
      S1=PS+GE
      S2=-PS1+ZET2
      T1=ZN/ZNK*(2.0D0*ZET3-(S1**2+S2)/ZNK)
      T=T+AK1(L)*T1
1     CONTINUE
C
      DO 2 L=1,9
      ZNK=ZN+DBLE(L)
      ZNK1=ZNK+ONE
      CALL PSI0(ZNK1,PS)
      S1=PS+GE
      T=T+CK5(L)*S1/ZNK*ZN
2     CONTINUE
C
      DO 3 K=1,10
      L=K-1
      ZNK=ZN+DBLE(L)
      ZNK1=ZNK+ONE
      CALL PSI0(ZNK1,PS)
      CALL PSI1(ZNK1,PS1)
      S1=PS+GE
      S2=-PS1+ZET2
      T1=ZN/ZNK*(S1**2+S2)
      T=T-DK5(K)*T1
3     CONTINUE
C
      ACG24=T
C
      RETURN
      END
      COMPLEX*16 FUNCTION ACG25(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: (S12(-X) -ZETA3/8)/(X-1)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
C
      REAL*8 AK1,AK2,GE,DL,ZET3,ZERO,ONE
C
      COMMON/ACLOG1/ AK1(9)
      COMMON/ACLOG2/ AK2(10)
      COMMON/ACCON1/ DL
      COMMON/ACCON3/ ZET3
      COMMON/ACCON4/ GE
      COMMON/VAL   / ZERO,ONE
C
C >>> M[I1(x)/(1+x)](N)
C
      U=ACG10(ZN)
C
      T=DCMPLX(5.0D0/16.0D0*ZET3*DL,ZERO)+U/2.0D0
      ZN1=ZN+ONE
      CALL PSI0(ZN1,PS)
      S1=PS+GE
      T=T-ZET3/8.0D0*S1
C
      DO 1 L=2,11
      K=L-1
      ZK=DCMPLX(DBLE(L),ZERO)
      ZNK=ZN+ZK
      ZNK1=ZNK+ONE
      CALL PSI0(ZNK1,PS)
      S1=PS+GE
      T=T+AK2(K)*S1/ZK/ZNK*ZN/2.0D0
1     CONTINUE
C
      CO=-5.0D0/8.0D0*ZET3
C
      DO 2 K=1,9
      ZK=DCMPLX(DBLE(K),ZERO)
      ZNK=ZN+ZK
      T1=CO
      DO 3 L=1,9
      ZNKL=ZN+DBLE(K+L)
      ZNKL1=ZNKL+ONE
      CALL PSI0(ZNKL1,PS)
      S1=PS+GE
      T1=T1+AK1(L)*S1/ZNKL
3     CONTINUE
      T=T+AK1(K)/ZNK*ZN/2.0D0*T1
2     CONTINUE
C
      ACG25=T
C
      RETURN
      END
      COMPLEX*16 FUNCTION ACG26(ZN)
C     -----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG^3(1-X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      REAL*8 AK1,ZERO,ONE,XLI4,ZET2,ZET3,GE
C
      COMMON/ACLOG1/ AK1(9)
      COMMON/ACCON2/ ZET2
      COMMON/ACCON3/ ZET3
      COMMON/ACCON4/ GE
      COMMON/ACCON5/ XLI4
      COMMON/VAL   / ZERO,ONE
C
      T=DCMPLX(-6.0D0*XLI4,ZERO)
      DO 1 K=1,9
      ZK=DCMPLX(DBLE(K),ZERO)
      ZK1=ZK+ONE
      ZKN=ZK+ZN
      ZKN1=ZKN+ONE
      CALL PSI0(ZKN1,PS)
      CALL PSI1(ZKN1,PS1)
      CALL PSI2(ZKN1,PS2)
      S1=PS+GE
      S2=-PS1+ZET2
      S3= PS2/2.0D0+ZET3
      CALL PSI0(ZK1,PS)
      CALL PSI1(ZK1,PS1)
      CALL PSI2(ZK1,PS2)
      H1=PS+GE
      H2=-PS1+ZET2
      H3= PS2/2.0D0+ZET3
      T=T-AK1(K)*(ZK/ZKN*(S1**3+3.0D0*S1*S2+2.0D0*S3)
     &                  -(H1**3+3.0D0*H1*H2+2.0D0*H3))
1     CONTINUE
C
      ACG26=T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FLI2(X)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  LI2(X)
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      EXTERNAL FR1
C
      PI = 3.141592653589793238462643D0
      ZETA2 = PI**2/6.0D0       
C
      IF(X.LT.-1.0D0.OR.X.GT.1.0D0) GOTO 100
      IF(X.EQ.1.0D0)  GOTO 5
      IF(X.EQ.-1.0D0) GOTO 6
      IF(X.LT.0.0D0)  GOTO 1
      IF(X.EQ.0.0D0)  GOTO 2
      IF(X.GT.0.5D0)  GOTO 3
      T=FR1(X)
      GOTO 200
100   WRITE(6,*) 'FLI2 -> NOT ALLOWED,X=',X,'STOP ***'
      STOP
1     Y=-X
      IF(Y.GT.0.5D0) GOTO 4
      Y2=Y**2
      T=FR1(Y2)/2.0D0-FR1(Y)
      GOTO 200
2     T=0.0D0
      GOTO 200
3     XM=1.0D0-X
      T=-FR1(XM)-LOG(X)*LOG(XM)+ZETA2
      GOTO 200
4     YM=1.0D0-Y
      T1=-FR1(YM)-LOG(Y)*LOG(YM)+ZETA2
      Y2=Y**2
      IF(Y2.LT.0.5) THEN
      T2=FR1(Y2)
      ELSE
      T2=-FR1(1.0D0-Y2)-LOG(Y2)*LOG(1.0D0-Y2)+ZETA2
      ENDIF
      T=T2/2.0D0-T1
      GOTO 200
5     T=ZETA2
      GOTO 200
6     T=-ZETA2/2.0D0
      GOTO 200
C 
200   FLI2=T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FR1(X)
C     --------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  SUBSIDIARTY ROUTINE FOR  LI2(X)
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      Y=LOG(1.0D0-X)
      Y2=Y*Y
C
      T = (-1.D0+(-1.D0/4.D0+(-1.D0/36.D0+(1.D0/3600.D0+(
     &-1.D0/211680.D0+(1.D0/10886400.D0+(-1.D0/526901760.D0
     &+(691.D0/16999766784000.D0+(-1.D0/1120863744000.D0
     &+3617.D0/0.18140058832896D18*Y2)*Y2)*Y2)*Y2)*Y2)*Y2)
     &*Y2)*Y)*Y)*Y
C
      FR1=T
C
      RETURN
      END
      SUBROUTINE GAMMAL(ZZ,RES)
C     -------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  LOG(GAMMA(Z)) FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 T1,T2,ZZ,Z,X,X2,RES
C
      Z=ZZ
      PI = 3.141592653589793238462643D0
C
      ONE=DCMPLX(1.0D0,0.0D0)
      T=DCMPLX(0.0D0,0.0D0)
2     R=SQRT(DREAL(Z)**2+DIMAG(Z)**2)
      IF(R.GT.10.0D0) GOTO 1
      T=T-LOG(Z)
      Z=Z+ONE
      GOTO 2
1     CONTINUE
C
      T1=Z*(LOG(Z)-1.0D0)+LOG(2.0D0*PI/Z)/2.0D0
C
      X=ONE/Z
      X2=X*X
      T2 = (1.D0/12.D0+(-1.D0/360.D0+(1.D0/1260.D0+(-1.D0/1680.D0+(1.D0/
     #1188.D0+(-691.D0/360360.D0+(1.D0/156.D0-3617.D0/122400.D0*X2)*X2
     #)*X2)*X2)*X2)*X2)*X2)*X
C
      RES=T1+T2+T
C
      RETURN
      END
      SUBROUTINE GAMMA(ZZ,RES)
C     ------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  GAMMA(Z) FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      COMPLEX*16 ZZ,T,RES
C
      CALL GAMMAL(ZZ,T)
C
      RES=EXP(T)
C
      RETURN
      END
      SUBROUTINE BETA(AA,BB,RES)
C     --------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C---  BETA(A,B) FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      COMPLEX*16 T1,T2,T3,T,RES,AA,BB
C
      CALL GAMMAL(AA,T1)
      CALL GAMMAL(BB,T2)
      CALL GAMMAL(AA+BB,T3)
      T=T1+T2-T3
C
      RES=EXP(T)
C
      RETURN
      END
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
C
      ONE=DCMPLX(1.0D0,0.0D0)
C
      Z=ZZ
      T=DCMPLX(0.0D0,0.0D0)
2     R=SQRT(DREAL(Z)**2+DIMAG(Z)**2)
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
C
      ONE=DCMPLX(1.0D0,0.0D0)
C
      Z=ZZ
      T=DCMPLX(0.0D0,0.0D0)
2     R=SQRT(DREAL(Z)**2+DIMAG(Z)**2)
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
C
      ONE=DCMPLX(1.0D0,0.0D0)
      TWO=ONE*2.0D0
C
      Z=ZZ
      T=DCMPLX(0.0D0,0.0D0)
2     R=SQRT(DREAL(Z)**2+DIMAG(Z)**2)
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
C
      ONE=DCMPLX(1.0D0,0.0D0)
      SIX=ONE*6.0D0
C
      Z=ZZ
      T=DCMPLX(0.0D0,0.0D0)
2     R=SQRT(DREAL(Z)**2+DIMAG(Z)**2)
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
      RES=(V1-V2)/8.0D0
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION SUM1(K,N)
C     -----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  SINGLE (ALTERNATING) HARMONIC SUM
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      K1 = ABS(K)
      S1 = 0.0D0
C
      DO 10 L=1,N
         IF ((K. LE. 0) .AND. (MOD(L,2) .NE. 0)) THEN
            S1 = S1 - 1.0D0/DBLE(L)**K1
         ELSE
            S1 = S1 + 1.0D0/DBLE(L)**K1
         END IF
 10   CONTINUE
C
      SUM1 = S1
      RETURN
      END
               
      DOUBLE PRECISION FUNCTION SUM2(J,K,N)
C     -------------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  DOUBLE (ALTERNATING) HARMONIC SUM
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      J1 = ABS(J)
      K1 = ABS(K)
C
      S1 = 0.0D0
      S2 = 0.0D0
C
      DO 10 L=1,N
         IF ((K .LE. 0) .AND. (MOD(L,2) .NE. 0)) THEN
            S1 = S1 - 1.0D0/DBLE(L)**K1
         ELSE
            S1 = S1 + 1.0D0/DBLE(L)**K1
         END IF
         IF ((J .LE. 0) .AND. (MOD(L,2) .NE. 0)) THEN
            S2 = S2 - S1/DBLE(L)**J1
         ELSE
            S2 = S2 + S1/DBLE(L)**J1
         END IF
 10   CONTINUE
C
      SUM2 = S2
      RETURN
      END
      DOUBLE PRECISION FUNCTION SUM3(I,J,K,N)
C     ---------------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  TRIPLE (ALTERNATING) HARMONIC SUM
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      I1 = ABS(I)
      J1 = ABS(J)
      K1 = ABS(K)
C
      S1 = 0.0D0
      S2 = 0.0D0
      S3 = 0.0D0
C
      DO 10 L=1,N
         IF ((K .LE. 0) .AND. (MOD(L,2) .NE. 0)) THEN
            S1 = S1 - 1.0D0/DBLE(L)**K1
         ELSE
            S1 = S1 + 1.0D0/DBLE(L)**K1
         END IF
         IF ((J .LE. 0) .AND. (MOD(L,2) .NE. 0)) THEN
            S2 = S2 - S1/DBLE(L)**J1
         ELSE
            S2 = S2 + S1/DBLE(L)**J1
         END IF
         IF ((I .LE. 0) .AND. (MOD(L,2) .NE. 0)) THEN
            S3 = S3 - S2/DBLE(L)**I1
         ELSE
            S3 = S3 + S2/DBLE(L)**I1
         END IF
 10   CONTINUE
C
      SUM3 = S3
      RETURN
      END
      DOUBLE PRECISION FUNCTION SUM4(I,J,K,L,N)
C     -----------------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  QUADRUPLE (ALTERNATING) HARMONIC SUM
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      I1 = ABS(I)
      J1 = ABS(J)
      K1 = ABS(K)
      L1 = ABS(L)
C
      T1 = 0.0D0
      T2 = 0.0D0
      T3 = 0.0D0
      T4 = 0.0D0
C
      DO 10 M=1,N
         IF ((L .LE. 0) .AND. (MOD(M,2) .NE. 0)) THEN
            T1 = T1 - 1.0D0/DBLE(M)**L1
         ELSE
            T1 = T1 + 1.0D0/DBLE(M)**L1
         END IF
         IF ((K .LE. 0) .AND. (MOD(M,2) .NE. 0)) THEN
            T2 = T2 - T1/DBLE(M)**K1
         ELSE
            T2 = T2 + T1/DBLE(M)**K1
         END IF
         IF ((J .LE. 0) .AND. (MOD(M,2) .NE. 0)) THEN
            T3 = T3 - T2/DBLE(M)**J1
         ELSE
            T3 = T3 + T2/DBLE(M)**J1
         END IF
         IF ((I .LE. 0) .AND. (MOD(M,2) .NE. 0)) THEN
            T4 = T4 - T3/DBLE(M)**I1
         ELSE
            T4 = T4 + T3/DBLE(M)**I1
         END IF
 10   CONTINUE
C
      SUM4 = T4
      RETURN
      END
      DOUBLE PRECISION FUNCTION SUM5(I,J,K,L,I10,N)
C     ---------------------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  FIVEFOLD  (ALTERNATING) HARMONIC SUM
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      I1 = ABS(I)
      J1 = ABS(J)
      K1 = ABS(K)
      L1 = ABS(L)
      II1= ABS(I10)
C
      T1 = 0.0D0
      T2 = 0.0D0
      T3 = 0.0D0
      T4 = 0.0D0
      T5 = 0.0D0
C
      DO 10 M=1,N
         IF ((I10 .LE. 0) .AND. (MOD(M,2) .NE. 0)) THEN
            T1 = T1 - 1.0D0/DBLE(M)**II1
         ELSE
            T1 = T1 + 1.0D0/DBLE(M)**II1
         END IF
         IF ((L .LE. 0) .AND. (MOD(M,2) .NE. 0)) THEN
            T2 = T2 - T1/DBLE(M)**L1
         ELSE
            T2 = T2 + T1/DBLE(M)**L1
         END IF
         IF ((K .LE. 0) .AND. (MOD(M,2) .NE. 0)) THEN
            T3 = T3 - T2/DBLE(M)**K1
         ELSE
            T3 = T3 + T2/DBLE(M)**K1
         END IF
         IF ((J .LE. 0) .AND. (MOD(M,2) .NE. 0)) THEN
            T4 = T4 - T3/DBLE(M)**J1
         ELSE
            T4 = T4 + T3/DBLE(M)**J1
         END IF
         IF ((I .LE. 0) .AND. (MOD(M,2) .NE. 0)) THEN
            T5 = T5 - T4/DBLE(M)**I1
         ELSE
            T5 = T5 + T4/DBLE(M)**I1
         END IF
 10   CONTINUE
C
      SUM5 = T5
      RETURN
      END
      DOUBLE PRECISION FUNCTION FLI3(X)
C     -----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  LI3(X) FOR -1. LE . X . LE .+1
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      A=1D0
      F=0D0
      AN=0D0
      TCH=1D-16
1     AN=AN+1D0
      A=A*X
      B=A/AN**3
      F=F+B
      IF(ABS(B)-TCH)2,2,1
2     FLI3=F
      END
      DOUBLE PRECISION FUNCTION FLI4(X)
C     -----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  LI4(X) FOR -1. LE . X . LE .+1
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      A=1D0
      F=0D0
      AN=0D0
      TCH=1D-16
1     AN=AN+1D0
      A=A*X
      B=A/AN**4
      F=F+B
      IF(ABS(B)-TCH)2,2,1
2     FLI4 = F
      END
      DOUBLE PRECISION FUNCTION S12(Z)
C     --------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  NIELSEN INTEGRAL: S12(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      EXTERNAL FINT
C
      EPS = 1.0D-8
      KEY = 2
      MAX = 10000
C
      F = DAIND1(0.0D0,Z,FINT,EPS,KEY,MAX,KOU,EST)
      S12 = F/2.0D0
      RETURN
      END
      DOUBLE PRECISION FUNCTION FINT(Y)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  SUBSIDIARY ROUTINE FOR THE NIELSEN INTEGRAL: S12(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      FINT = LOG(1.0D0-Y)**2/Y
      RETURN
      END
      DOUBLE PRECISION FUNCTION YI1(T)
C     --------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  I1(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      EXTERNAL YIF
C
      YI1 = DAIND1(0.0D0,T,YIF,EPS,KEY,MAX,KOU,EST)
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION YIF(R)
C     --------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  INTEGRAND OF I1(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C
      YIF = LOG(1.0D0+R)*LOG(1.0D0-R)/R
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION DAIND(A,B,FUN,EPS,KEY,MAX,KOUNT,EST)
C     --------------------------------------------------------------
C************************************************************************
C
C---  INTEGRATION ROUTINE:
C     CF. R. PIESSENS, ANGEW. INFORMATIK, VOL. 9 (1973) 399.
C
C************************************************************************
C
C  INPUTPARAMETERS
C  A,B      LIMITS OF THE INTEGRATION INTERVAL
C  FUN      FUNCTION TO BE INTEGRATED (TO BE DECLARED EXTERNAL IN THE MAIN PR.)
C  EPS      ABSOLUTE OR RELATIVE TOLERANCE,DEPENDING OF THE VALUE OF 'KEY'
C  KEY      =1 THEN 'EPS' DENOTES AN ABSOLUTE, =2 THEN A RELATIVE TOLERANCE
C  MAX      UPPER BOUND ON THE NUMBERS OF INTEGRAND EVALUATIONS (MAX.LE.10000)
C
C  OUTPUTPARAMETERS
C  KOUNT    NUMBER OF INTEGRAND EVALUATIONS
C  EST      ESTIMATION OF THE ABSOLUTE ERROR OF THE APPROXIMATION
      IMPLICIT REAL*8 (A-H,O-Z)
      DOUBLE PRECISION MAXIM,MINIM,MODUL1,MODUL2
      INTEGER RANG(130)
      DIMENSION
     &AINIT(250),END(250),EPSIL(250),PART(250),W1(5),W2(5),W3(6),
     &     X1(5),X2(5)
      DATA X1/0.973906528517D+0,0.865063366689D+0,0.679409568299D+0,
     *         0.433395394129D+0,0.148874338981D+0/
      DATA X2/0.995657163026D+0,0.930157491356D+0,0.780817726586D+0,
     *         0.562757134669D+0,0.294392862701D+0/
      DATA W1/0.666713443087D-1,0.149451349151D+0,0.219086362516D+0,
     *         0.269266719310D+0,0.295524224715D+0/
      DATA W2/0.325581623080D-1,0.750396748109D-1,0.109387158802D+0,
     *         0.134709217311D+0,0.147739104901D+0/
      DATA W3/0.116946388674D-1,0.547558965744D-1,0.931254545837D-1,
     *         0.123491976262D+0,0.142775938577D+0,0.149445554003D+0/
      DATA TOL/0.23D-15/
      EXTERNAL FUN
      MAX1 = (MAX+21)/42+1
      MAX2 = MAX1/2+2
      ALFA = A
      BETA = B
      MAAT = 1
C EVALUATION OF GAUSSIAN AND KRONROD FORMULAS
   10 S = 0.5D+0*(BETA-ALFA)
      U = 0.5D+0*(BETA+ALFA)
      RES1 = 0.0D+0
      RES2 = W3(6)*FUN(U)
      DO 20 K = 1,5
        C = S*X1(K)
        C = FUN(C+U)+FUN(U-C)
        RES1 = RES1+W1(K)*C
        RES2 = RES2+W2(K)*C
        C = S*X2(K)
   20   RES2 =RES2+W3(K)*(FUN(C+U)+FUN(U-C))
      PAT = RES2*S
      MODUL2 = ABS(PAT-RES1*S)
      IF(MAAT.GT.1) GOTO 50
      EST = MODUL2
      BINT = PAT
      KOUNT =21
      PART(1) = BINT
      GOTO 90
   30 RANG(1) = 1
      AINIT(1) = A
      END(1) = B
      EPSIL(1) = EST
   40 NR = RANG(1)
      BINT = BINT-PART(NR)
      EST =EST-EPSIL(NR)
C THE SUBINTERVAL WITH LARGEST ERROR IS SPLIT UP INTO TWO EQUAL PARTS
      ALFA = AINIT(NR)
      BETA = (AINIT(NR)+END(NR))*0.5
      JJ = 1
      MAAT = MAAT+1
      GOTO 10
   50 EST = EST+MODUL2
      BINT = BINT+PAT
      IF(JJ.EQ.0) GOTO 60
      MODUL1 = MODUL2
      PAT1 = PAT
      ALFA = BETA
      BETA = END(NR)
      JJ = 0
      GOTO 10
   60 MA = MAAT
      IF(MAAT.GT.MAX2) MA = MAX1+3-MAAT
      IF(MODUL1.GT.MODUL2) GOTO 70
      EPSIL(NR) = MODUL2
      EPSIL(MAAT) = MODUL1
      AINIT(MAAT) = AINIT(NR)
      AINIT(NR) = ALFA
      END(MAAT) = ALFA
      MAXIM = MODUL2
      MINIM = MODUL1
      PART(NR) = PAT
      PART(MAAT) = PAT1
      GOTO 80
   70 EPSIL(NR) = MODUL1
      EPSIL(MAAT) = MODUL2
      END(MAAT) = BETA
      END(NR) = ALFA
      AINIT(MAAT) = ALFA
      MAXIM = MODUL1
      MINIM = MODUL2
      PART(NR) = PAT1
      PART(MAAT) = PAT
   80 KOUNT = KOUNT+42
C TEST ON THE NUMBER OF FUNCTION EVALUATIONS
      IF(KOUNT.GE.MAX) GOTO 190
   90 GOTO (100,110),KEY
C TEST ON ABSOLUTE ACCURACY
  100 IF(EST.LE.EPS) GOTO 190
      GOTO 120
C TEST ON RELATIVE ACCURACY
  110 IF(ABS(EPS*BINT).LE.TOL) GOTO 100
      IF(EST.LE.ABS(EPS*BINT)) GOTO 190
  120 IF(MAAT.EQ.1) GOTO 30
      IF(MAAT.GT.2) GOTO 130
      RANG(2) = 2
      GOTO 40
  130 MB = MA-1
C SEARCH FOR THE SUBINTERVAL WITH LARGEST ERROR
      DO 140 I = 2,MB
        IR = RANG(I)
        IF(MAXIM.GE.EPSIL(IR)) GOTO 150
  140   RANG(I-1) = RANG(I)
      RANG(MB) = NR
      RANG(MA) = MAAT
      GOTO 40
  150 RANG(I-1) = NR
      DO 160 K = I,MB
        IR = RANG(K)
        IF(MINIM.GE.EPSIL(IR)) GOTO 170
  160   CONTINUE
      RANG(MA) = MAAT
      GOTO 40
  170 DO 180 I = K,MB
        KK = MB-I+K
  180   RANG(KK+1) = RANG(KK)
      RANG(K) = MAAT
      GOTO 40
C CALCULATION OF THE INTEGRAL
  190 AIND1 = 0.0D+0
      DO 200 K = 1,MAAT
  200   AIND1 = AIND1+PART(K)
      IF(AIND1.EQ.0.0D+0)
     &        WRITE(6,*) '**** AIND=0.**** EST NOT CALCULATED'
      IF(AIND1.NE.0.0D+0) EST=EST/AIND1
      DAIND=AIND1
      RETURN
      END

      DOUBLE PRECISION FUNCTION DAIND1(A,B,FUN,EPS,KEY,MAX,KOUNT,EST)
C     --------------------------------------------------------------
C************************************************************************
C
C---  INTEGRATION ROUTINE:
C     CF. R. PIESSENS, ANGEW. INFORMATIK, VOL. 9 (1973) 399.
C
C************************************************************************
C  INPUTPARAMETERS
C  A,B      LIMITS OF THE INTEGRATION INTERVAL
C  FUN      FUNCTION TO BE INTEGRATED (TO BE DECLARED EXTERNAL IN THE MAIN PR.)
C  EPS      ABSOLUTE OR RELATIVE TOLERANCE,DEPENDING OF THE VALUE OF 'KEY'
C  KEY      =1 THEN 'EPS' DENOTES AN ABSOLUTE, =2 THEN A RELATIVE TOLERANCE
C  MAX      UPPER BOUND ON THE NUMBERS OF INTEGRAND EVALUATIONS (MAX.LE.10000)
C
C  OUTPUTPARAMETERS
C  KOUNT    NUMBER OF INTEGRAND EVALUATIONS
C  EST      ESTIMATION OF THE ABSOLUTE ERROR OF THE APPROXIMATION
      IMPLICIT REAL*8 (A-H,O-Z)
      DOUBLE PRECISION MAXIM,MINIM,MODUL1,MODUL2
      INTEGER RANG(130)
      DIMENSION
     &AINIT(250),END(250),EPSIL(250),PART(250),W1(5),W2(5),W3(6),
     &     X1(5),X2(5)
      DATA X1/0.973906528517D+0,0.865063366689D+0,0.679409568299D+0,
     *         0.433395394129D+0,0.148874338981D+0/
      DATA X2/0.995657163026D+0,0.930157491356D+0,0.780817726586D+0,
     *         0.562757134669D+0,0.294392862701D+0/
      DATA W1/0.666713443087D-1,0.149451349151D+0,0.219086362516D+0,
     *         0.269266719310D+0,0.295524224715D+0/
      DATA W2/0.325581623080D-1,0.750396748109D-1,0.109387158802D+0,
     *         0.134709217311D+0,0.147739104901D+0/
      DATA W3/0.116946388674D-1,0.547558965744D-1,0.931254545837D-1,
     *         0.123491976262D+0,0.142775938577D+0,0.149445554003D+0/
      DATA TOL/0.23D-15/
      EXTERNAL FUN
      MAX1 = (MAX+21)/42+1
      MAX2 = MAX1/2+2
      ALFA = A
      BETA = B
      MAAT = 1
C EVALUATION OF GAUSSIAN AND KRONROD FORMULAS
   10 S = 0.5D+0*(BETA-ALFA)
      U = 0.5D+0*(BETA+ALFA)
      RES1 = 0.0D+0
      RES2 = W3(6)*FUN(U)
      DO 20 K = 1,5
        C = S*X1(K)
        C = FUN(C+U)+FUN(U-C)
        RES1 = RES1+W1(K)*C
        RES2 = RES2+W2(K)*C
        C = S*X2(K)
   20   RES2 =RES2+W3(K)*(FUN(C+U)+FUN(U-C))
      PAT = RES2*S
      MODUL2 = ABS(PAT-RES1*S)
      IF(MAAT.GT.1) GOTO 50
      EST = MODUL2
      BINT = PAT
      KOUNT =21
      PART(1) = BINT
      GOTO 90
   30 RANG(1) = 1
      AINIT(1) = A
      END(1) = B
      EPSIL(1) = EST
   40 NR = RANG(1)
      BINT = BINT-PART(NR)
      EST =EST-EPSIL(NR)
C THE SUBINTERVAL WITH LARGEST ERROR IS SPLIT UP INTO TWO EQUAL PARTS
      ALFA = AINIT(NR)
      BETA = (AINIT(NR)+END(NR))*0.5
      JJ = 1
      MAAT = MAAT+1
      GOTO 10
   50 EST = EST+MODUL2
      BINT = BINT+PAT
      IF(JJ.EQ.0) GOTO 60
      MODUL1 = MODUL2
      PAT1 = PAT
      ALFA = BETA
      BETA = END(NR)
      JJ = 0
      GOTO 10
   60 MA = MAAT
      IF(MAAT.GT.MAX2) MA = MAX1+3-MAAT
      IF(MODUL1.GT.MODUL2) GOTO 70
      EPSIL(NR) = MODUL2
      EPSIL(MAAT) = MODUL1
      AINIT(MAAT) = AINIT(NR)
      AINIT(NR) = ALFA
      END(MAAT) = ALFA
      MAXIM = MODUL2
      MINIM = MODUL1
      PART(NR) = PAT
      PART(MAAT) = PAT1
      GOTO 80
   70 EPSIL(NR) = MODUL1
      EPSIL(MAAT) = MODUL2
      END(MAAT) = BETA
      END(NR) = ALFA
      AINIT(MAAT) = ALFA
      MAXIM = MODUL1
      MINIM = MODUL2
      PART(NR) = PAT1
      PART(MAAT) = PAT
   80 KOUNT = KOUNT+42
C TEST ON THE NUMBER OF FUNCTION EVALUATIONS
      IF(KOUNT.GE.MAX) GOTO 190
   90 GOTO (100,110),KEY
C TEST ON ABSOLUTE ACCURACY
  100 IF(EST.LE.EPS) GOTO 190
      GOTO 120
C TEST ON RELATIVE ACCURACY
  110 IF(ABS(EPS*BINT).LE.TOL) GOTO 100
      IF(EST.LE.ABS(EPS*BINT)) GOTO 190
  120 IF(MAAT.EQ.1) GOTO 30
      IF(MAAT.GT.2) GOTO 130
      RANG(2) = 2
      GOTO 40
  130 MB = MA-1
C SEARCH FOR THE SUBINTERVAL WITH LARGEST ERROR
      DO 140 I = 2,MB
        IR = RANG(I)
        IF(MAXIM.GE.EPSIL(IR)) GOTO 150
  140   RANG(I-1) = RANG(I)
      RANG(MB) = NR
      RANG(MA) = MAAT
      GOTO 40
  150 RANG(I-1) = NR
      DO 160 K = I,MB
        IR = RANG(K)
        IF(MINIM.GE.EPSIL(IR)) GOTO 170
  160   CONTINUE
      RANG(MA) = MAAT
      GOTO 40
  170 DO 180 I = K,MB
        KK = MB-I+K
  180   RANG(KK+1) = RANG(KK)
      RANG(K) = MAAT
      GOTO 40
C CALCULATION OF THE INTEGRAL
  190 AIND1 = 0.0D+0
      DO 200 K = 1,MAAT
  200   AIND1 = AIND1+PART(K)
      IF(AIND1.EQ.0.0D+0)
     &        WRITE(6,*) '**** AIND=0.**** EST NOT CALCULATED'
      IF(AIND1.NE.0.0D+0) EST=EST/AIND1
      DAIND1=AIND1
      RETURN
      END
      DOUBLE PRECISION FUNCTION FCT1(N)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN MOMENT FOR POSITIVE INTEGER ARGUMENT OF FKT1(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ NN
C
      EXTERNAL FKT1
C
      NN=N
      F = DAIND(0.0D0,1.0D0,FKT1,EPS,KEY,MAX,KOU,EST)
C
      FCT1=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKT1(X)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR:  LOG(1+X)/(1+X)
C---  INTEGRAND FOR
C---  NUMERICAL INTEGRAL FOR POSITIVE INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ N
C
C
      F = LOG(1.0D0+X)/(1.0D0+X)*X**N
C
      FKT1=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FCT2(N)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN MOMENT FOR POSITIVE INTEGER ARGUMENT OF FKT2(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ NN
C
      EXTERNAL FKT2
C
      NN=N
      F = DAIND(0.0D0,1.0D0,FKT2,EPS,KEY,MAX,KOU,EST)
C
      FCT2=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKT2(X)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR:  LOG^2(1+X)/(1+X)
C---  INTEGRAND FOR
C---  NUMERICAL INTEGRAL FOR POSITIVE INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ N
C
C
      T = LOG(1.0D0+X)
      F = T**2*X**N/(1.0D0+X)
C
      FKT2=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FCT3(N)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN MOMENT FOR POSITIVE INTEGER ARGUMENT OF FKT3(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ NN
C
      EXTERNAL FKT3
C
      NN=N
      F = DAIND(0.0D0,1.0D0,FKT3,EPS,KEY,MAX,KOU,EST)
C
      FCT3=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKT3(X)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR:  LI2(X)/(1+X)
C---  INTEGRAND FOR
C---  NUMERICAL INTEGRAL FOR POSITIVE INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ N
C
C
      F = FLI2(X)/(1.0D0+X)*X**(N)
C
      FKT3=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FCT4(N)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN MOMENT FOR POSITIVE INTEGER ARGUMENT OF FKT4(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ NN
C
      EXTERNAL FKT4
C
      NN=N
      F = DAIND(0.0D0,1.0D0,FKT4,EPS,KEY,MAX,KOU,EST)
C
      FCT4=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKT4(X)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR:  LI2(-X)/(1+X)
C---  INTEGRAND FOR
C---  NUMERICAL INTEGRAL FOR POSITIVE INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ N
C
C
      F = FLI2(-X)/(1.0D0+X)*X**(N)
C
      FKT4=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FCT5(N)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN MOMENT FOR POSITIVE INTEGER ARGUMENT OF FKT5(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ NN
C
      EXTERNAL FKT5
C
      NN=N
      F = DAIND(0.0D0,1.0D0,FKT5,EPS,KEY,MAX,KOU,EST)
C
      FCT5=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKT5(X)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR:  LOG(X)*LI2(X)/(1+X)
C---  INTEGRAND FOR
C---  NUMERICAL INTEGRAL FOR POSITIVE INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ N
C
C
      F = LOG(X)*FLI2(X)/(1.0D0+X)*X**(N)
C
      FKT5=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FCT6(N)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN MOMENT FOR POSITIVE INTEGER ARGUMENT OF FKT6(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ NN
C
      EXTERNAL FKT6
C
      NN=N
      F = DAIND(0.0D0,1.0D0,FKT6,EPS,KEY,MAX,KOU,EST)
C
      FCT6=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKT6(X)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR:  LI3(X)/(1+X)
C---  INTEGRAND FOR
C---  NUMERICAL INTEGRAL FOR POSITIVE INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ N
C
C
      F = FLI3(X)/(1.0D0+X)*X**(N)
C
      FKT6=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FCT7(N)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN MOMENT FOR POSITIVE INTEGER ARGUMENT OF FKT7(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ NN
C
      EXTERNAL FKT7
C
      NN=N
      F = DAIND(0.0D0,1.0D0,FKT7,EPS,KEY,MAX,KOU,EST)
C
      FCT7=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKT7(X)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR:  LI3(-X)/(1+X)
C---  INTEGRAND FOR
C---  NUMERICAL INTEGRAL FOR POSITIVE INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ N
C
C
      F = FLI3(-X)/(1.0D0+X)*X**(N)
C
      FKT7=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FCT8(N)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN MOMENT FOR POSITIVE INTEGER ARGUMENT OF FKT8(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ NN
C
      EXTERNAL FKT8
C
      NN=N
      F = DAIND(0.0D0,1.0D0,FKT8,EPS,KEY,MAX,KOU,EST)
C
      FCT8=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKT8(X)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR:  S12(X)/(1+X)
C---  INTEGRAND FOR
C---  NUMERICAL INTEGRAL FOR POSITIVE INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ N
C
C
      F = S12(X)/(1.0D0+X)*X**(N)
C
      FKT8=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FCT9(N)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN MOMENT FOR POSITIVE INTEGER ARGUMENT OF FKT9(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ NN
C
      EXTERNAL FKT9
C
      NN=N
      F = DAIND(0.0D0,1.0D0,FKT9,EPS,KEY,MAX,KOU,EST)
C
      FCT9=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKT9(X)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR:  S12(-X)/(1+X)
C---  INTEGRAND FOR
C---  NUMERICAL INTEGRAL FOR POSITIVE INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ N
C
C
      F = S12(-X)/(1.0D0+X)*X**(N)
C
      FKT9=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FCT10(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN MOMENT FOR POSITIVE INTEGER ARGUMENT OF FKT10(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ NN
C
      EXTERNAL FKT10
C
      NN=N
      F = DAIND(0.0D0,1.0D0,FKT10,EPS,KEY,MAX,KOU,EST)
C
      FCT10=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKT10(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR:  I1(X)/(1+X)
C---  INTEGRAND FOR
C---  NUMERICAL INTEGRAL FOR POSITIVE INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ N
      EXTERNAL YI1
C
C
      T1=X**N
      F = YI1(X)*T1/(1.0D0+X)         
C
      FKT10=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FCT11(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN MOMENT FOR POSITIVE INTEGER ARGUMENT OF FKT11(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ NN
C
      EXTERNAL FKT11
C
      NN=N
      F = DAIND(0.0D0,1.0D0,FKT11,EPS,KEY,MAX,KOU,EST)
C
      FCT11=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKT11(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR:  LOG(1-X)*LI2(X)/(1+X)
C---  INTEGRAND FOR
C---  NUMERICAL INTEGRAL FOR POSITIVE INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ N
C
C
      F = LOG(1.0D0-X)*FLI2(X)/(1.0D0+X)*X**(N)
C
      FKT11=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FCT12(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN MOMENT FOR POSITIVE INTEGER ARGUMENT OF FKT12(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ NN
C
      EXTERNAL FKT12
C
      NN=N
      F = DAIND(0.0D0,1.0D0,FKT12,EPS,KEY,MAX,KOU,EST)
C
      FCT12=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKT12(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR:  LOG(1-X)*LI2(-X)/(1+X)
C---  INTEGRAND FOR
C---  NUMERICAL INTEGRAL FOR POSITIVE INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ N
C
C
      X2=X**2
      F = LOG(1.0D0-X)*FLI2(-X)/(1.0D0+X)*X**(N)
C
      FKT12=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FCT13(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN MOMENT FOR POSITIVE INTEGER ARGUMENT OF FKT13(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ NN
C
      EXTERNAL FKT13
C
      NN=N
      F = DAIND(0.0D0,1.0D0,FKT13,EPS,KEY,MAX,KOU,EST)
C
      FCT13=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKT13(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR:  LOG(1+X)*LI2(-X)/(1+X)
C---  INTEGRAND FOR
C---  NUMERICAL INTEGRAL FOR POSITIVE INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ N
C
      F = LOG(1.0D0+X)*FLI2(-X)/(1.0D0+X)*X**(N)
C
      FKT13=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FCT14(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN MOMENT FOR POSITIVE INTEGER ARGUMENT OF FKT14(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ NN
C
      EXTERNAL FKT14
C
      NN=N
      F = DAIND(0.0D0,1.0D0,FKT14,EPS,KEY,MAX,KOU,EST)
C
      FCT14=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKT14(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR:  (LOG^2(1+X)-LOG^2(2))/(X-1)
C---  INTEGRAND FOR
C---  NUMERICAL INTEGRAL FOR POSITIVE INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ N
C
      T = LOG(1.0D0+X)
      F = (T**2             -D2**2)/(X-1.0D0)*X**(N)
C
      FKT14=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FCT15(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN MOMENT FOR POSITIVE INTEGER ARGUMENT OF FKT15(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ NN
C
      EXTERNAL FKT15
C
      NN=N
      F = DAIND(0.0D0,1.0D0,FKT15,EPS,KEY,MAX,KOU,EST)
C
      FCT15=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKT15(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR:  (LOG(1+X)-LOG(2))/(X-1)*LI2(X)
C---  INTEGRAND FOR
C---  NUMERICAL INTEGRAL FOR POSITIVE INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ N
C
C
      F = (LOG(1.0D0+X)-D2)/(X-1.0D0)*FLI2(X)*X**(N)
C
      FKT15=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FCT16(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN MOMENT FOR POSITIVE INTEGER ARGUMENT OF FKT16(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ NN
C
      EXTERNAL FKT16
C
      NN=N
      F = DAIND(0.0D0,1.0D0,FKT16,EPS,KEY,MAX,KOU,EST)
C
      FCT16=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKT16(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR:  (LOG(1+X)-LOG(2))/(X-1)*LI2(-X)
C---  INTEGRAND FOR
C---  NUMERICAL INTEGRAL FOR POSITIVE INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ N
C
C
      F =  (LOG(1.0D0+X)-D2)/(X-1.0D0)*X**(N)*FLI2(-X)
C
      FKT16=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FCT17(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN MOMENT FOR POSITIVE INTEGER ARGUMENT OF FKT17(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ NN
C
      EXTERNAL FKT17
C
      NN=N
      F = DAIND(0.0D0,1.0D0,FKT17,EPS,KEY,MAX,KOU,EST)
C
      FCT17=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKT17(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR:  LOG(X)*LOG^2(1+X)/(X-1)
C---  INTEGRAND FOR
C---  NUMERICAL INTEGRAL FOR POSITIVE INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ N
C
C
      T = LOG(1.0D0+X)
      F = LOG(X)* T**2*(X**N      )/(X-1.0D0)
C
      FKT17=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FCT18(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN MOMENT FOR POSITIVE INTEGER ARGUMENT OF FKT18(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ NN
C
      EXTERNAL FKT18
C
      NN=N
      F = DAIND(0.0D0,0.99999999999999D0,FKT18,EPS,KEY,MAX,KOU,EST)
C
      FCT18=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKT18(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR:  (LI2(X)-ZETA2)/(X-1)
C---  INTEGRAND FOR
C---  NUMERICAL INTEGRAL FOR POSITIVE INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ N
C
C
      F = (FLI2(X)-ZETA2)/(X-1.0D0)*X**N
C
      FKT18=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FCT19(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN MOMENT FOR POSITIVE INTEGER ARGUMENT OF FKT19(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ NN
C
      EXTERNAL FKT19
C
      NN=N
      F = DAIND(0.0D0,1.0D0,FKT19,EPS,KEY,MAX,KOU,EST)
C
      FCT19=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKT19(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR:  (LI2(-X)+ZETA2/2)/(X-1)
C---  INTEGRAND FOR
C---  NUMERICAL INTEGRAL FOR POSITIVE INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ N
C
C
      F = (FLI2(-X)+ZETA2/2.0D0)/(X-1.0D0)*X**N
C
      FKT19=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FCT20(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN MOMENT FOR POSITIVE INTEGER ARGUMENT OF FKT20(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ NN
C
      EXTERNAL FKT20
C
      NN=N
      F = DAIND(0.0D0,1.0D0,FKT20,EPS,KEY,MAX,KOU,EST)
C
      FCT20=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKT20(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR:  (LI3(X)-ZETA3)/(X-1)
C---  INTEGRAND FOR
C---  NUMERICAL INTEGRAL FOR POSITIVE INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ N
C
C
      F = (FLI3(X)-ZETA3)/(X-1.0D0)*X**N
C
      FKT20=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FCT21(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN MOMENT FOR POSITIVE INTEGER ARGUMENT OF FKT21(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ NN
C
      EXTERNAL FKT21
C
      NN=N
      F = DAIND(0.0D0,0.99999999999D0,FKT21,EPS,KEY,MAX,KOU,EST)
C
      FCT21=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKT21(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR:  (S12(X)-ZETA3)/(X-1)
C---  INTEGRAND FOR
C---  NUMERICAL INTEGRAL FOR POSITIVE INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ N
C
C
      F = (S12(X)-ZETA3)/(X-1.0D0)*X**N
C
      FKT21=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FCT22(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN MOMENT FOR POSITIVE INTEGER ARGUMENT OF FKT22(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ NN
C
      EXTERNAL FKT22
C
      NN=N
      F = DAIND(0.0D0,1.0D0,FKT22,EPS,KEY,MAX,KOU,EST)
C
      FCT22=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKT22(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR:  LOG(X)*LI2(X)/(X-1)
C---  INTEGRAND FOR
C---  NUMERICAL INTEGRAL FOR POSITIVE INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ N
C
C
      F = LOG(X)*FLI2(X)/(X-1.0D0)*X**N
C
      FKT22=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FCT23(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN MOMENT FOR POSITIVE INTEGER ARGUMENT OF FKT23(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ NN
C
      EXTERNAL FKT23
C
      NN=N
      F = DAIND(0.0D0,1.0D0,FKT23,EPS,KEY,MAX,KOU,EST)
C
      FCT23=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKT23(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR:  (LI3(-X)+3/4*ZETA2)/(X-1)
C---  INTEGRAND FOR
C---  NUMERICAL INTEGRAL FOR POSITIVE INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ N
C
C
      F = (FLI3(-X)+3.0D0/4.0D0*ZETA3)/(X-1.0D0)*X**N
C
      FKT23=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FCT24(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN MOMENT FOR POSITIVE INTEGER ARGUMENT OF FKT24(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ NN
C
      EXTERNAL FKT24
C
      NN=N
      F = DAIND(0.0D0,1.0D0,FKT24,EPS,KEY,MAX,KOU,EST)
C
      FCT24=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKT24(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR:  (I1(X)+5/8*ZETA3)/(X-1)
C---  INTEGRAND FOR
C---  NUMERICAL INTEGRAL FOR POSITIVE INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ N
C
      F = (YI1(X)+5.0D0/8.0D0*ZETA3)/(X-1)*X**N
C
      FKT24=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FCT25(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN MOMENT FOR POSITIVE INTEGER ARGUMENT OF FKT25(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ NN
C
      EXTERNAL FKT25
C
      NN=N
      F = DAIND(0.0D0,1.0D0,FKT25,EPS,KEY,MAX,KOU,EST)
C
      FCT25=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKT25(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR:  (S12(X)-ZETA3)/(X-1)
C---  INTEGRAND FOR
C---  NUMERICAL INTEGRAL FOR POSITIVE INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ N
C
C
      F = (S12(-X)-ZETA3/8.0D0)/(X-1.0D0)*X**N
C
      FKT25=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FCT26(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN MOMENT FOR POSITIVE INTEGER ARGUMENT OF FKT25(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ NN
C
      EXTERNAL FKT26
C
      NN=N
      F = DAIND(0.0D0,1.0D0,FKT26,EPS,KEY,MAX,KOU,EST)
C
      FCT26=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKT26(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR:  LOG^3(1-X)/(1+X)
C---  INTEGRAND FOR
C---  NUMERICAL INTEGRAL FOR POSITIVE INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
      COMMON /EP/EPS
      COMMON /KM/KEY,MAX
      COMMON/EXPO/ N
C
      F=(LOG(1.0D0-X))**3/(1.0D0+X)*X**N
C
      FKT26=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION XCG1(N)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(1+X)/(1+X)
C---  SUM REPRESENTATION FOR INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
C
      DL=D2
      N1=N
C
      IF(MOD(N1,2).NE.1) THEN
         IFA= 1
         ELSE
         IFA=-1
      ENDIF
C
      T1=SUM2(-1,1,N1)
      T2=SUM1(1,N1)
      T3=SUM1(-1,N1)
      T4=SUM1(-2,N1)
C
      T=(DL*DL/2.0D0+T1-(T2-T3)*DL-T2*T3-T4)*IFA
C
      XCG1=T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION XCG2(N)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG^2(1+X)/(1+X)
C---  SUM REPRESENTATION FOR INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
C
      DL=D2
      N1=N
      IF(MOD(N1,2).NE.1) THEN
         IFA= 1
         ELSE
         IFA=-1
      ENDIF
      T1=SUM3(1,1,-1,N1)
      T2=SUM2(1,-1,N1)
      T3=SUM2(1,1,N1)
      T4=SUM1(1,N1)
      T5=SUM1(-1,N1)
C
      T=(T1-DL*(T2-T3)-DL**2/2.0D0*(T4-T5)+DL**3/6.0D0)*IFA*2.0D0
C
      XCG2=T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION XCG3(N)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LI2(X)/(1+X)
C---  SUM REPRESENTATION FOR INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
C
      DL=D2
      N1=N
      IF(MOD(N1,2).NE.1) THEN
         IFA= 1
         ELSE
         IFA=-1
      ENDIF
      T1=SUM2(-2,1,N1)
      T2=SUM1(-1,N1)
C
      IFA=-IFA
      T= (T1-ZETA2*T2+5.0D0/8.0D0*ZETA3-ZETA2*DL)*IFA
C
      XCG3=T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION XCG4(N)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LI2(-X)/(1+X)
C---  SUM REPRESENTATION FOR INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
C
      DL=D2
      N1=N
      IF(MOD(N1,2).NE.1) THEN
         IFA= 1
         ELSE
         IFA=-1
      ENDIF
      T1=SUM2(2,-1,N1)
      T2=SUM1(-1,N1)
      T3=SUM1(2,N1)
      T4=SUM1(-2,N1)
C
      IFA=-IFA
      T= (T1+DL*(T3-T4)+ZETA2*T2/2.0D0-ZETA3/4.0D0+ZETA2*DL/2.0D0)*IFA
C
      XCG4=T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION XCG5(N)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(X)*LI2(X)/(1+X)
C---  SUM REPRESENTATION FOR INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
C
      DL=D2
      XLI4=ZLI4
      N1=N
      IF(MOD(N1,2).NE.1) THEN
         IFA= 1
         ELSE
         IFA=-1
      ENDIF
      T0=SUM2(-2,2,N1)
      T1=SUM2(-3,1,N1)
      T2=SUM1(-2,N1)
      T=(T0+2.0D0*T1-2.0D0*ZETA2*T2-ZETA2**2*3.0D0/40.0D0)*IFA
C
      XCG5=T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION XCG6(N)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LI3(X)/(1+X)
C---  SUM REPRESENTATION FOR INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
C
      DL=D2
      XLI4=ZLI4
      N1=N
      IF(MOD(N1,2).NE.1) THEN
         IFA= 1
         ELSE
         IFA=-1
      ENDIF
      T1=SUM2(-3,1,N1)
      T2=SUM1(-2,N1)
      T3=SUM1(-1,N1)
C
      T=(T1-ZETA2*T2+ZETA3*T3+3.0D0/5.0D0*ZETA2**2-2.0D0*XLI4
     &  -3.0D0/4.0D0*ZETA3*DL+ZETA2*DL**2/2.0D0-DL**4/12.0D0)*IFA
C
      XCG6=T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION XCG7(N)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LI3(-X)/(1+X)
C---  SUM REPRESENTATION FOR INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
C
      DL=D2
      XLI4=ZLI4
      N1=N
      IF(MOD(N1,2).NE.1) THEN
         IFA= 1
         ELSE
         IFA=-1
      ENDIF
      T1=SUM2(3,-1,N1)
      T2=SUM1(3,N1)
      T3=SUM1(-3,N1)
      T4=SUM1(-2,N1)
      T5=SUM1(-1,N1)
C
      T=(T1+DL*(T2-T3)+ZETA2/2.0D0*T4-3.0D0/4.0D0*ZETA3*T5
     &   +ZETA2**2/8.0D0
     &   -3.0D0/4.0D0*ZETA3*DL)*IFA
C
      XCG7=T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION XCG8(N)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: S12(X)/(1+X)
C---  SUM REPRESENTATION FOR INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
C
      DL=D2
      XLI4=ZLI4
      N1=N
      IF(MOD(N1,2).NE.1) THEN
         IFA= 1
         ELSE
         IFA=-1
      ENDIF
      T1=SUM3(-2,1,1,N1)
      T2=SUM1(-1,N1)
C
      T=-IFA*(T1+XLI4-ZETA3*T2-ZETA2**2/8.0D0-ZETA3*DL/8.0D0
     &  -ZETA2*DL**2/4.0D0+DL**4/24.0D0)
C
      XCG8=T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION XCG9(N)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: S12(-X)/(1+X)
C---  SUM REPRESENTATION FOR INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
C
      DL=D2
      XLI4=ZLI4
C
      N1=N
      IF(MOD(N1,2).NE.1) THEN
         IFA= 1
         ELSE
         IFA=-1
      ENDIF
      T1=SUM3(2,1,-1,N1)
      T2=SUM2(2,1,N1)
      T3=SUM2(2,-1,N1)
      T4=SUM1(2,N1)
      T5=SUM1(-2,N1)
      T6=SUM1(-1,N1)
C
      T=-(T1+DL*(T2-T3)-DL*DL/2.0D0*(T4-T5)-ZETA3/8.0D0*T6-3.0D0*XLI4
     &  +6.0D0/5.0D0*ZETA2**2-11.0D0/4.0D0*ZETA3*DL
     &  +3.0D0/4.0D0*ZETA2*DL*DL-DL**4/8.0D0)*IFA
C
      IFA=-IFA
      XCG9=T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION XCG10(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: I1(X)/(1+X)
C---  SUM REPRESENTATION FOR INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
C
      DL=D2
      XLI4=ZLI4
C
      N1=N
      IF(MOD(N1,2).NE.1) THEN
         IFA= 1
         ELSE
         IFA=-1
      ENDIF
C
      T1=SUM3(-2,-1,-1,N1)
      T2=SUM3(2,-1,1,N1)
      T3=SUM2(-2,1,N1)
      T4=SUM2(-2,-1,N1)
      T5=SUM1(2,N1)
      T6=SUM1(-2,N1)
      T7=SUM1(-1,N1)
C
      T= -(T1+T2-DL*(T3-T4)+(T5-T6)/2.0D0*(ZETA2-DL**2)
     &    +5.0D0/8.0D0*ZETA3*(T7+DL)-3.0D0/20.0D0*ZETA2**2)*IFA
C
      IFA=-IFA
      XCG10=T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION XCG11(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(1-X)*LI2(X)/(1+X)
C---  SUM REPRESENTATION FOR INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
C
      DL=D2
      XLI4=ZLI4
C
      N1=N
      IF(MOD(N1,2).NE.1) THEN
         IFA= 1
         ELSE
         IFA=-1
      ENDIF
      T1=SUM3(-1,2,1,N1)
      T2=SUM2(-1,1,N1)
      T3=SUM3(-2,1,1,N1)
      T4=SUM1(-1,N1)
C
      S12R=-IFA*(T3+XLI4-ZETA3*T4-ZETA2**2/8.0D0-ZETA3*DL/8.0D0
     &  -ZETA2*DL**2/4.0D0+DL**4/24.0D0)
C
      T=(T1-ZETA2*T2-19.0D0/40.0D0*ZETA2**2+XLI4+ZETA3*DL/4.0D0
     &  +DL**2*ZETA2/4.0D0+DL**4/24.0D0)*IFA-2.0D0*S12R
C
      XCG11=T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION XCG12(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(1-X)*LI2(-X)/(1+X)
C---  SUM REPRESENTATION FOR INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
C
      DL=D2
      XLI4=ZLI4
C
      N1=N
      IF(MOD(N1,2).NE.1) THEN
         IFA= 1
         ELSE
         IFA=-1
      ENDIF
C
      T1=SUM3(2,-1,1,N1)
      T2=SUM3(-2,-1,-1,N1)
      T3=SUM3(-1,-2,-1,N1)
      T4=SUM2(-2, 1,N1)
      T5=SUM2(-1,2,N1)
      T6=SUM2(-1,1,N1)
      T7=SUM1(-1,N1)
      T8=SUM1(-2,N1)
      T9=SUM1(3,N1)
      T10=SUM1(2,N1)
      T11=SUM1(-1,N1)
C
      T=(T1+T2+T3-(T4+T5)*DL+T6*ZETA2/2.0D0+T7*T8*DL+T9*DL
     &+(ZETA2-DL**2)/2.0D0*(T10-T8)+5.0D0/8.0D0*ZETA3*T7-4.0D0*XLI4
     &+3.0D0/2.0D0*ZETA2**2-21.0D0/8.0D0*ZETA3*DL
     &+3.0D0/4.0D0*ZETA2*DL**2-DL**4/6.0D0)*IFA
C
      XCG12=T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION XCG13(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(1+X)*LI2(-X)/(1+X)
C---  SUM REPRESENTATION FOR INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
C
      DL=D2
      XLI4=ZLI4
      N1=N
      IF(MOD(N1,2).NE.1) THEN
         IFA= 1
         ELSE
         IFA=-1
      ENDIF

C
      T1=SUM3(1,2,-1,N1)
      T2=SUM3(2,1,-1,N1)
      T3=SUM2(2,1,N1)
      T4=SUM2(1,-2,N1)
      T5=SUM2(2,-1,N1)
      T6=SUM1(1,N1)
      T7=SUM1(2,N1)
      T8=SUM1(3,N1)
      T9=SUM1(-1,N1)
      T10=SUM2(1,-1,N1)
      T11=SUM1(-2,N1)
      T12=SUM1(1 ,N1)
C
      T=(T1+2.0D0*T2+(T3-T4-2*T5+T6*T7+T8-ZETA2/2.0D0*T9)*DL
     &  +T10*ZETA2/2.0D0-(T7-T11)*DL**2
     &  -(ZETA3/4.0D0-ZETA2*DL/2.0D0)*T12
     &  -3.0D0*XLI4+6.0D0/5.0D0*ZETA2**2-21.0D0/8.0D0*ZETA3*DL
     &  +ZETA2/2.0D0*DL**2-DL**4/8.0D0)*IFA
C
      XCG13=T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION XCG14(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: (LOG^2(1+X)-LOG^2(2))/(X-1)
C---  SUM REPRESENTATION FOR INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
C
      DL=D2
      XLI4=ZLI4
      N1=N
C
      T1=SUM3(-1,1,-1,N1)
      T2=SUM2(-1,-1,N1)
      T3=SUM2(-1,1,N1)
      T4=SUM1(1,N1)
      T5=SUM1(-1,N1)
C
      T=(T1-DL*(T2-T3)+DL*DL/2.0D0*(T4-T5))*2.0D0
      T=T-DL**2*T4-ZETA3/4.0D0+ZETA2*DL-2.0D0/3.0D0*DL**3
C
      XCG14=T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION XCG15(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: (LOG(1+X)-LOG(2))/(X-1)*LI2(X)
C---  SUM REPRESENTATION FOR INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
C
      DL=D2
      XLI4=ZLI4
      N1=N
C
      T1=SUM3(-1,-2,1,N1)
      T2=SUM3(2,-1,-1,N1)
      T3=SUM3(-2,-1,1,N1)
      T4=SUM2(-1,-1,N1)
      T5=SUM1(1,N1)
      T6=SUM1(-1,N1)
      T7=SUM2(2,1,N1)
      T8=SUM2(2,-1,N1)
      T9=SUM1(2,N1)
      T10=SUM1(-2,N1)
C
      T=T1+T2+T3-ZETA2*T4-(5.0D0/8.0D0*ZETA3-ZETA2*DL)*(T5-T6)
     &  +5.0D0/8.0D0*ZETA3*T5-DL*(T7-T8)
     &  -(ZETA2-DL**2)/2.0D0*(T9-T10)
C >>
     &  -DL*(   ZETA2*T5)
     &  +19.0D0/40.0D0*ZETA2**2-XLI4+7.0D0/4.0D0*ZETA3*DL
     &  -ZETA2*DL**2/4.0D0-DL**4/24.0D0
     &  +DL*(SUM2(2,1,N1)-ZETA3*2.0D0)
C
      XCG15=T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION XCG16(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: (LOG(1+x)-LOG(2))/(X-1)*LI2(-X)
C---  SUM REPRESENTATION FOR INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
C
      DL=D2
      XLI4=ZLI4
      N1=N
C
      T= 2.0D0*SUM3(-2,1,-1,N)+SUM3(-1,2,-1,N)
     & +DL*(2.0D0*(-SUM2(-2,-1,N)+SUM2(-2,1,N))
     &            -SUM2(-1,-2,N)+SUM2(-1,2,N))
     & +ZETA2/2.0D0*SUM2(-1,-1,N)-DL**2*(SUM1(-2,N)-SUM1(2,N))
     & -ZETA3/4.0D0*SUM1(1,N)
     & -(ZETA3/4.0D0-ZETA2*DL/2.0D0)*(SUM1(-1,N)-SUM1(1,N))
C>>
     & +DL*(SUM2(-2,-1,N)-DL*(SUM1(2,N)-SUM1(-2,N))
     & +ZETA2/2.0D0*SUM1(1,N))
     & -33.0D0/20.0D0*ZETA2**2+4.0D0*XLI4
     & +13.0D0/4.0D0*ZETA3*DL-3.0D0/4.0D0*ZETA2*DL**2+DL**4/6.0D0
C
      XCG16=T
C
      RETURN


      END
      DOUBLE PRECISION FUNCTION XCG17(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(X)*LOG^2(1+X)/(X-1)
C---  SUM REPRESENTATION FOR INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
C
      N1=N
      DL=D2
      XLI4=ZLI4
C
      T1= SUM2(-1,-2,N)                                                  "
      T2= SUM2(-1, 2,N)
      T3= SUM2(-1,-1,N)
      T4= SUM2(-1, 1,N)
      T5= SUM1(-1,N)
      T6= SUM1( 1,N)
      T7= SUM2(-2, 1,N)
      T8= SUM2(-2,-1,N)
      T9= SUM1(-2,   N)
      T10=SUM1( 2,   N)
C
      U1=DL*(T1-T2)+(ZETA3/8.0D0)*(T5-T6)
      U2=-ZETA2/2.0D0*T4
      U3=-(T7-T8)*DL+(T9-T10)*DL**2/2.0D0 +ZETA3/8.0D0*T6
C
      TEST1B=2.0D0*(SUM3(-1,2,-1,N)+SUM3(-1,1,-2,N)+SUM3(-2,1,-1,N))
      T     =  (U1+U2+U3)*2.0D0-TEST1B
      CO    =  7.0D0/4.0D0*ZETA2**2-4.0D0*ZLI4-21.0D0/4.0D0*ZETA3*D2
     &        +5.0D0/2.0D0*ZETA2*D2**2-D2**4/6.0D0
C
      XCG17= T+CO
      RETURN
      END
      DOUBLE PRECISION FUNCTION XCG18(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: (LI2(X)-ZETA2)/(X-1)
C---  SUM REPRESENTATION FOR INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
C
      N1=N
      DL=D2
      XLI4=ZLI4
C
      N1=N
C
      T1=SUM2(2,1,N1)
C
      T=-T1+2.0D0*ZETA3
C
      XCG18= T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION XCG19(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: (LI2(-X)+ZETA2/2)/(X-1)
C---  SUM REPRESENTATION FOR INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
C
      N1=N
      DL=D2
      XLI4=ZLI4
C
      T1=SUM2(-2,-1,N1)
      T2=SUM1(2,N1)
      T3=SUM1(-2,N1)
      T4=SUM1(1,N1)
C
      T=-T1+DL*(T2-T3)-5.0D0/8.0D0*ZETA3
C
      XCG19= T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION XCG20(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: (LI3(X)-ZETA3)/(X-1)
C---  SUM REPRESENTATION FOR INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
C
      N1=N
      DL=D2
      XLI4=ZLI4
C
      T1=SUM1(1,N1)
      T2=SUM1(2,N1)
      T3=SUM2(3,1,N1)
C
      T=         -ZETA2*T2+T3+ZETA2**2/2.0D0
C
      XCG20= T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION XCG21(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: (S12(X)-ZETA3)/(X-1)
C---  SUM REPRESENTATION FOR INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
C
      N1=N
      DL=D2
      XLI4=ZLI4
C
      T1=SUM3(2,1,1,N1)
      T2=SUM1(1,N1)
C
      T=-T1+6.0D0/5.0D0*ZETA2**2
C
      XCG21= T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION XCG22(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(X)*LI2(X)/(X-1)
C---  SUM REPRESENTATION FOR INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
C
      DL=D2
      XLI4=ZLI4
C
      N1=N
C
      T1=SUM1(2,N1)
      T2=SUM1(4,N1)
      T3=SUM2(3,1,N1)
C
      T=-2.0D0*ZETA2*T1+T1**2/2.0D0+T2/2.0D0+2.0D0*T3
     &  +3.0D0/10.0D0*ZETA2**2
C
      XCG22=T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION XCG23(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: (LI3(-X)+3*ZETA3/4)/(X-1)
C---  SUM REPRESENTATION FOR INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
C
      DL=D2
      XLI4=ZLI4
C
      N1=N
C
      T1=SUM2(-3,-1,N1)
      T2=SUM1(3,N1)
      T3=SUM1(-3,N1)
      T4=SUM1(2,N1)
      T5=SUM1(1,N1)
C
      T= T1-(T2-T3)*DL+ZETA2/2.0D0*T4-3.0D0/4.0D0*ZETA3*T5
     &  +2.0D0*XLI4-11.0D0/10.0D0*ZETA2**2+7.0D0/4.0D0*ZETA3*DL
     &  -ZETA2/2.0D0*DL**2+DL**4/12.0D0+3.0D0/4.0D0*ZETA3*T5
C
      XCG23= T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION XCG24(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: (I1(X)+5/8*ZETA3)/(X-1)
C---  SUM REPRESENTATION FOR INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
C
      DL=D2
      XLI4=ZLI4
      N1=N
C
      T1=SUM3(2,-1,-1,N1)
      T2=SUM3(-2,-1, 1,N1)
      T3=SUM2(2, 1,   N1)
      T4=SUM2(2,-1,   N1)
      T5=SUM1(2,N1)
      T6=SUM1(-2,N1)
      T7=SUM1(1,N1)
C
      T= 
C        -5.0D0/8.0D0*ZETA3*T7
     &  -T1-T2+DL*(T3-T4)+(ZETA2-DL**2)/2.0D0*(T5-T6)
     &  +ZETA2**2/4.0D0-2.0D0*XLI4-7.0D0/4.0D0*ZETA3*DL
     &   +ZETA2*DL**2/2.0D0        -DL**4/12.0D0
C
      XCG24=T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION XCG25(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: (S12(-X)-ZETA(3)/8)/(X-1)
C---  SUM REPRESENTATION FOR INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
C
      DL=D2
      XLI4=ZLI4
C
      N1=N
C
      T1=SUM3(-2,1,-1,N1)
      T2=SUM2(-2,1,N1)
      T3=SUM2(-2,-1,N1)
      T4=SUM1(-2,N1)
      T5=SUM1(2,N1)
      T6=SUM1(1,N1)
C
      T= -T1-(T2-T3)*DL+(T4-T5)*DL**2/2.0D0+3.0D0/40.0D0*ZETA2**2
C
      XCG25=T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION XCG26(N)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR:  LOG^3(1-X)/(1+X)
C---  SUM REPRESENTATION FOR INTEGER N
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /PIZ/PI,ZETA2,ZETA3,ZLI4,D2
C
      DL=D2
      XLI4=ZLI4
C
      T1=SUM4(-1,1,1,1,N)
      IF(MOD(N,2).EQ.1) IFA=1
      IF(MOD(N,2).EQ.0) IFA=-1
C
      T=(T1+XLI4)*6*IFA
C
      XCG26=T
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKN1(X)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  FUNCTION               LOG(1+X)/(1+X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      F = LOG(1.0D0+X)/(1.0D0+X)
C
      FKN1=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKN2(X)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  FUNCTION               LOG^2(1+X)/(1+X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C
      T = LOG(1.0D0+X)
      F = T**2/(1.0D0+X)
C
      FKN2=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKN3(X)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  FUNCTION               LI2(X)/(1+X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      F = FLI2(X)/(1.0D0+X)
C
      FKN3=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKN4(X)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  FUNCTION               LI2(-X)/(1+X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C
      F = FLI2(-X)/(1.0D0+X)
C
      FKN4=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKN5(X)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  FUNCTION               LOG(X)*LI2(X)/(1+X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      F = LOG(X)*FLI2(X)/(1.0D0+X)
C
      FKN5=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKN6(X)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  FUNCTION               LI3(X)/(1+X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      F = FLI3(X)/(1.0D0+X)
C
      FKN6=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKN7(X)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  FUNCTION               LI3(-X)/(1+X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      F = FLI3(-X)/(1.0D0+X)
C
      FKN7=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKN8(X)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  FUNCTION               S12(X))/(1+X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      F = S12(X)/(1.0D0+X)
C
      FKN8=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKN9(X)
C     ---------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  FUNCTION               S12(-X))/(1+X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      F = S12(-X)/(1.0D0+X)
C
      FKN9=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKN10(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  FUNCTION               I1(X)/(1+X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      EXTERNAL YI1
C
      F = YI1(X)/(1.0D0+X)         
C
      FKN10=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKN11(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  FUNCTION               LOG(1-X)*LI2(X)/(1+X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      F = LOG(1.0D0-X)*FLI2(X)/(1.0D0+X)
C
      FKN11=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKN12(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  FUNCTION               LOG(1-X)*LI2(-X)/(1+X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      F = LOG(1.0D0-X)*FLI2(-X)/(1.0D0+X)
C
      FKN12=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKN13(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  FUNCTION               LOG(1+X)*LI2(-X)/(1+X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      F = LOG(1.0D0+X)*FLI2(-X)/(1.0D0+X)
C
      FKN13=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKN14(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  FUNCTION               (LOG^2(1+X)-LOG^2(2))/(X-1)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      T = LOG(1.0D0+X)
      F = (T**2 -(LOG(2.0D0))**2)/(X-1.0D0)
C
      FKN14=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKN15(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  FUNCTION               (LOG(1+X)-LOG(2))/(X-1)*LI2(X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      F = (LOG(1.0D0+X)-LOG(2.0D0))/(X-1.0D0)*FLI2(X)
C
      FKN15=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKN16(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  FUNCTION               (LOG(1+X)-LOG(2))/(X-1)*LI2(-X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      F =  (LOG(1.0D0+X)-LOG(2.0D0))/(X-1.0D0)*FLI2(-X)
C
      FKN16=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKN17(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  FUNCTION                LOG(X)*LOG^2(1+X)/(X-1)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      T = LOG(1.0D0+X)
      F = LOG(X)* T**2/(X-1.0D0)
C
      FKN17=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKN18(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  FUNCTION                (LI2(X)-ZETA2)/(X-1) 
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /ZET2/ ZETA2
C
      F = (FLI2(X)-ZETA2)/(X-1.0D0)
C
      FKN18=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKN19(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  FUNCTION                (LI2(-X)+ZETA2/2)/(X-1)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /ZET2/ ZETA2
C
      F = (FLI2(-X)+ZETA2/2.0D0)/(X-1.0D0)        
C
      FKN19=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKN20(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  FUNCTION                (LI3(X)-ZETA3)/(X-1)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /ZET3/ ZETA3
C
      F = (FLI3(X)-ZETA3)/(X-1.0D0)
C
      FKN20=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKN21(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  FUNCTION                (S12(X)-ZETA3)/(X-1)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /ZET3/ ZETA3
C
      F = (S12(X)-ZETA3)/(X-1.0D0)       
C
      FKN21=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKN22(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  FUNCTION                LOG(X)*LI2(X)/(X-1)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      F = LOG(X)*FLI2(X)/(X-1.0D0)
C
      FKN22=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKN23(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  FUNCTION                (LI3(-X)+3/4*ZETA3)/(X-1)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /ZET3/ ZETA3
C
      F = (FLI3(-X)+3.0D0/4.0D0*ZETA3)/(X-1.0D0)
C
      FKN23=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKN24(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  FUNCTION                (I1(X)+5/8*ZETA3)/(X-1)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /ZET3/ ZETA3
C
      F = (YI1(X)+5.0D0/8.0D0*ZETA3)/(X-1)
C
      FKN24=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKN25(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  FUNCTION                (S12(-X)-ZETA3/8)/(X-1)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /ZET3/ ZETA3
C
      F = (S12(-X)-ZETA3/8.0D0)/(X-1.0D0)
C
      FKN25=F
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FKN26(X)
C     ----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  FUNCTION                LOG^3(1-X)/(1+X)
C
C************************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      F=(LOG(1.0D0-X))**3/(1.0D0+X)
C
      FKN26=F
C
      RETURN
      END
      SUBROUTINE UINIT
C     ----------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  USER INITIALIZATION OF RUNNING FLAGS
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C---  CHANGE THE DEFAULT VALUES OF FLAGS AND RUNNING PARAMETERS
C
      COMMON/IAPP  / IAPP
      COMMON/RUN / IRUN
      COMMON /EP/EPS
      COMMON/TEST/ ITEST1,ITEST2,ITEST3
      COMMON/MOMPA/ NMIN,NMAX
      COMMON/FUNPA/ IMIN,IMAX
C
      IRUN = 0
      EPS = 1.0D-9
      ITEST1=1
      ITEST2=1
      ITEST3=1
      IAPP  =1
      IMIN=1
      IMAX=26
      NMIN=1
      NMAX=20
C
      RETURN
      END
      SUBROUTINE URUN
C     ----------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  USER RUNNING:
C---  HERE THE FUNCTIONS & SUBROUTINES
C---  FCTi(N), ACGi(Z) and XCGi(N) may be accessed and combined
C---  to other structures
C
C************************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16
     &         ACG1,ACG2,ACG3,ACG4,ACG5,ACG6,ACG7,ACG8,ACG9,ACG10,
     &         ACG11,ACG12,ACG13,ACG14,ACG15,ACG16,ACG17,ACG18,ACG19,
     &         ACG20,ACG21,ACG22,ACG23,ACG24,ACG25,ACG26
C
      EXTERNAL FCT1,FCT2,FCT3,FCT4,FCT5,FCT6,FCT7,FCT8,FCT9,FCT10,
     &         FCT11,FCT12,FCT13,FCT14,FCT15,FCT16,FCT17,FCT18,FCT19,
     &         FCT20,FCT21,FCT22,FCT23,FCT24,FCT25,FCT26
      EXTERNAL XCG1,XCG2,XCG3,XCG4,XCG5,XCG6,XCG7,XCG8,XCG9,XCG10,
     &         XCG11,XCG12,XCG13,XCG14,XCG15,XCG16,XCG17,XCG18,XCG19,
     &         XCG20,XCG21,XCG22,XCG23,XCG24,XCG25,XCG26
      EXTERNAL ACG1,ACG2,ACG3,ACG4,ACG5,ACG6,ACG7,ACG8,ACG9,ACG10,
     &         ACG11,ACG12,ACG13,ACG14,ACG15,ACG16,ACG17,ACG18,ACG19,
     &         ACG20,ACG21,ACG22,ACG23,ACG24,ACG25,ACG26
C
      RETURN
      END
      SUBROUTINE UOUT
C     ----------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  USER OUTPUT
C
C************************************************************************
C
      RETURN
      END
