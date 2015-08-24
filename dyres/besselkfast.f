c.....Function BK(n,z)
c.....BK(n,z) is the n-derivative of BesselK[nu,z]
c.....with respect to nu in nu=1

c.....Itilde defined as in the paper

      function Itilde(m)
      implicit none
      double precision zbesselk0,zbesselk1,zbesselk2,zbesselk3
      double precision argum,Itilde
      double precision Eulergamma,b0,z2,z3,logx
      double precision xmio 
      integer m
      common/xmio/xmio


      Eulergamma=0.577215664902d0
      z2=1.64493406685d0
      z3=1.20205690316d0
      b0=2*dexp(-Eulergamma)

C

      argum=b0*xmio
      logx=dlog(xmio)
      
      if (m.eq.1) then
         Itilde=-zbesselk0(argum)/xmio**2
      elseif (m.eq.2) then
         Itilde=2d0/xmio**2*(zbesselk0(argum)*logx-zbesselk1(argum))
      elseif (m.eq.3) then
         Itilde=-3/xmio**2*(zbesselk0(argum)*(logx**2-z2)
     &          -2*zbesselk1(argum)*logx+zbesselk2(argum))
      elseif (m.eq.4) then
         Itilde=4/xmio**2*(zbesselk0(argum)*(logx**3-3*z2*logx+2*z3)
     &         -3*zbesselk1(argum)*(logx**2-z2)+3*zbesselk2(argum)*logx
     &         -zbesselk3(argum))
      endif


      return
      end

C     n-derivative of the function BesselK[nu,z] 
C     with respect to nu for nu=1
C     NOTE: IT IS MULTIPLIED by z


      function zbesselk0(z)
      implicit none
      double precision zbesselk0,z,zm,loz,egamma,pi,r(0:11),zbk0
      data r(0)/ 6.03844076705d2/
      data r(1)/ -1.21597891877d2/
      data r(2)/ 2.72488273113d1/
      data r(3)/ -6.88391426811d0/
      data r(4)/ 1.99353173375d0/
      data r(5)/ -0.676592588425d0/
      data r(6)/ 0.277576446533d0/
      data r(7)/ -0.144195556641d0/
      data r(8)/ 0.102539062500d0/
      data r(9)/ -0.1171875d0/
      data r(10)/ 0.375d0/
      data r(11)/ 1d0/

      egamma=0.577215664902d0
      pi=3.14159265359d0
      if(z.lt.1.5d0) then
      zm=z/2d0
      loz=dlog(zm)

      zbesselk0=1d0+z*zm*(loz-0.5d0*(1-2*egamma))
     &         +zm**4*(loz-0.5d0*(2.5d0-2*egamma))
     &         +zm**6/6d0*(loz-0.5d0*(10d0/3-2*egamma))
     &         +zm**8/72d0*(loz-0.5d0*(47d0/12-2*egamma))
     &         +zm**10/1440d0*(loz-0.5d0*(131d0/30-2*egamma))
      elseif(z.ge.1.5d0.and.z.lt.4d0) then
      zbesselk0=zbk0(z)
      else
      zbesselk0=dsqrt(Pi/2d0)*(z)**(-10.5d0)*dexp(-z)*
     &    (r(0)+r(1)*z+r(2)*z**2+r(3)*z**3+r(4)*z**4+r(5)*z**5+
     &     r(6)*z**6+r(7)*z**7+r(8)*z**8+r(9)*z**9+r(10)*z**10+
     &     r(11)*z**11)
      endif
      return
      end


      function zbesselk1(z)
      implicit none
      double precision zbesselk1,z,zm,loz,egamma,pi,r(0:11),zbk1
      data r(0) /-5.51335896122d2/
      data r(1) /1.10017140269d2/
      data r(2) /-2.43805296996d1/
      data r(3) /6.07404200127d0/
      data r(4) /-1.72772750258d0/
      data r(5) /0.572501420975d0/
      data r(6) /-0.227108001709d0/
      data r(7) /0.112152099609d0/
      data r(8) /-7.32421875d-2/
      data r(9) /7.03125d-2/
      data r(10) /-0.125d0/
      data r(11) /1d0/

      egamma=0.577215664902d0
      pi=3.14159265359d0
      if(z.lt.1.5d0) then
      zm=z/2d0
      loz=dlog(zm)
      zbesselk1=-(loz+egamma)-zm**2*(loz-1+egamma)
     &         -0.25d0*zm**4*(loz-1.5d0+egamma)
     &         -zm**6/36d0*(loz-11d0/6+egamma)
     &         -zm**8/576d0*(loz-25d0/12+egamma)
      elseif(z.ge.1.5d0.and.z.lt.4d0) then
      zbesselk1=zbk1(z)
      else
      zbesselk1=dsqrt(Pi/2d0)*(z)**(-11.5d0)*dexp(-z)*
     &    (r(0)+r(1)*z+r(2)*z**2+r(3)*z**3+r(4)*z**4+r(5)*z**5+
     &     r(6)*z**6+r(7)*z**7+r(8)*z**8+r(9)*z**9+r(10)*z**10+
     &     r(11)*z**11)
      endif
      return
      end


      function zbesselk2(z)
      implicit none
      double precision zbesselk2,z,a(0:13),loz,zm,pi,r(0:11),zbk2
      data a(0) / 1.15443132980306572d0/
      data a(1) / 1.97811199065594511d0/
      data a(2) / 0.154431329803065721d0/
      data a(3) / 4.801792651508824500d0/
      data a(4) / 0.806235643470665767d0/
      data a(5) /-0.672784335098467139d0/
      data a(6) / 3.285072828402112960d0/
      data a(7) /-1.945338757678943440d0/
      data a(8) /-0.181575166960855634d0/
      data a(9) / 0.694195147571435559d0/
      data a(10)/-0.607655744858515573d0/
      data a(11)/-0.019182189839330562d0/
      data a(12)/ 0.068894530444636532d0/
      data a(13)/-0.070514317816328185d0/

      data r(0) /3.19461756880d4/
      data r(1) /-5.82903207466d3/
      data r(2) /1.17069096329d3/
      data r(3) /-2.61456867712d2/
      data r(4) /6.57620334072d1/
      data r(5) /-18.9305966582d0/
      data r(6) /6.37010269165d0/
      data r(7) /-2.57905883789d0/
      data r(8) /1.30957031250d0/
      data r(9) /-0.888020833333d0/
      data r(10) /0.875d0/
      data r(11) /1d0/

      pi=3.14159265359d0
      if(z.lt.1.5d0) then
      zm=z/2
      loz=dlog(zm)
      
      zbesselk2=loz**2+a(0)*loz+a(1)
     &       +zm**2*(2*loz**3/3d0+a(2)*loz**2+a(3)*loz+a(4))
     &       +zm**4*(loz**3/3d0+a(5)*loz**2+a(6)*loz+a(7))
     &       +zm**6*(loz**3/18d0+a(8)*loz**2+a(9)*loz+a(10))
     &       +zm**8*(loz**3/216d0+a(11)*loz**2+a(12)*loz+a(13))
      elseif(z.ge.1.5d0.and.z.lt.4d0) then
      zbesselk2=zbk2(z)
      else
      zbesselk2=dsqrt(Pi/2d0)*(z)**(-11.5d0)*dexp(-z)*
     &    (r(0)+r(1)*z+r(2)*z**2+r(3)*z**3+r(4)*z**4+r(5)*z**5+
     &     r(6)*z**6+r(7)*z**7+r(8)*z**8+r(9)*z**9+r(10)*z**10+
     &     r(11)*z**11)
      endif
      return
      end


      function zbesselk3(z)
      implicit none
      double precision zbesselk3,z,b(0:14),loz,zm,pi,r(0:9),zbk3

      data b(0) / 1.731646994704598580d0/
      data b(1) / 5.934335971967835330d0/
      data b(2) / 5.444874456485317730d0/
      data b(3) /-1.268353005295401420d0/ 
      data b(4) / 8.471041982558638170d0/
      data b(5) /-3.026167526073320430d0/
      data b(6) /-0.692088251323850355d0/ 
      data b(7) / 2.809848746963509900d0/
      data b(8) /-2.161466255000085060d0/
      data b(9) /-0.104676472369316706d0/
      data b(10)/ 0.381989731242156681d0/
      data b(11)/-0.367492827636283900d0/
      data b(12)/-0.007844362856415627d0/
      data b(13)/ 0.027796539630842606d0/
      data b(14)/-0.029917436634978395d0/

      data r(0)/-3.19152148877d3/
      data r(1)/7.05641513542d2/
      data r(2)/-1.75295543138d2/
      data r(3)/4.96775524480d1/
      data r(4)/-1.63798988342d1/
      data r(5)/6.45276489258d0/
      data r(6)/-3.15332031250d0/
      data r(7)/2.0234375d0/
      data r(8)/-1.875d0/
      data r(9)/3d0/

      pi=3.14159265359d0
      if(z.lt.1.5d0) then
       zm=z/2
       loz=dlog(zm)
      
       zbesselk3=loz**3+b(0)*loz**2+b(1)*loz+b(2)
     &       +zm**2*(loz**3+b(3)*loz**2+b(4)*loz+b(5))
     &       +zm**4*(loz**3/4d0+b(6)*loz**2+b(7)*loz+b(8))
     &       +zm**6*(loz**3/36d0+b(9)*loz**2+b(10)*loz+b(11))
     &       +zm**8*(loz**3/576d0+b(12)*loz**2+b(13)*loz+b(14))
       zbesselk3=-zbesselk3
      elseif(z.ge.1.5d0.and.z.lt.4d0) then
      zbesselk3=zbk3(z)
      else
       zbesselk3=dsqrt(Pi/2d0)*(z)**(-10.5d0)*dexp(-z)*
     &    (r(0)+r(1)*z+r(2)*z**2+r(3)*z**3+r(4)*z**4+r(5)*z**5+
     &     r(6)*z**6+r(7)*z**7+r(8)*z**8+r(9)*z**9)
      endif
      return
      end

C     n-derivative of the function BesselK[nu,z] 
C     with respect to nu for nu=1, multiplied by z

C     Interpolated form from z=1 to z=5

      function zbk0(x)
      implicit none
      integer i,j
      real*8 xa(1:21),ya(1:21),xx(1:5),yy(1:5),zbk0,x,y,dy   
	   DATA XA/1d0,1.2d0,1.4d0,1.6d0,1.8d0,2d0,2.2d0,2.4d0,
     &             2.6d0,2.8d0,3d0,3.2d0,3.4d0,3.6d0,3.8d0,4d0,
     &             4.2d0,4.4d0,4.6d0,4.8d0,5d0/
           DATA YA/0.601907230197d0,0.521510869273d0,0.449170263122d0,
     &             0.385014258172d0,0.328721579643d0,0.279731763633d0,
     &             0.237372982262d0,0.200939613012d0,0.169738517152d0,
     &             1.431155197d-1,1.20469293385d-1,1.01257264676d-1,
     &             8.49965460188d-2,7.12618632710d-2,5.96817704982d-2,
     &             4.99339955491d-2,4.17404598909d-2,3.48623147904d-2,
     &             2.90952007636d-2,2.42648469315d-2,2.02230672273d-2/
  
         call dylocate(xa,21,x,j)
         if (j.lt.2) j=2
         if (j.gt.20) j=20
         do i=1,4
         xx(i)=xa(j+i-2)
         yy(i)=ya(j+i-2)
         enddo
         call xpolint(xx,yy,4,x,y,dy)
	    zbk0=y
         if(x.lt.1d0) zbk0=0d0

       end

      function zbk1(x)
      implicit none
      integer i,j
      real*8 xa(1:21),ya(1:21),xx(1:5),yy(1:5),zbk1,x,y,dy   
	   DATA XA/1d0,1.2d0,1.4d0,1.6d0,1.8d0,2d0,2.2d0,2.4d0,
     &             2.6d0,2.8d0,3d0,3.2d0,3.4d0,3.6d0,3.8d0,4d0,
     &             4.2d0,4.4d0,4.6d0,4.8d0,5d0/
           DATA YA/0.421024438241d0,0.318508220287d0,0.243655061182d0,
     &             0.187954751969d0,0.14593140049d0,0.11389387275d0,
     &             8.92690056716d-2,7.02173415434d-2,5.53983032863d-2,
     &             4.3819981975d-2,3.47395043863d-2,2.75949976751d-2,
     &             2.19580188068d-2,1.74996410181d-2,1.39658845342d-2,
     &             1.11596760859d-2,8.92745154154d-3,7.14911062331d-3,
     &             5.73042291729d-3,4.59724631672d-3,3.69109833404d-3/
  
         call dylocate(xa,21,x,j)
         if (j.lt.2) j=2
         if (j.gt.20) j=20
         do i=1,4
         xx(i)=xa(j+i-2)
         yy(i)=ya(j+i-2)
         enddo
         call xpolint(xx,yy,4,x,y,dy)
	    zbk1=y
         if(x.lt.1d0) zbk1=0d0

       end

      function zbk2(x)
      implicit none
      integer i,j
      real*8 xa(1:21),ya(1:21),xx(1:5),yy(1:5),zbk2,x,y,dy   
	   DATA XA/1d0,1.2d0,1.4d0,1.6d0,1.8d0,2d0,2.2d0,2.4d0,
     &             2.6d0,2.8d0,3d0,3.2d0,3.4d0,3.6d0,3.8d0,4d0,
     &             4.2d0,4.4d0,4.6d0,4.8d0,5d0/

           DATA YA/0.680674889238d0,0.491594990820d0,0.362195189493d0,
     &             0.270816279222d0,0.204792327117d0,0.156253340447d0,
     &             0.120082673664d0,9.28359137653d-2,7.21300592003d-2,
     &             5.62799660668d-2,4.40726122853d-2,3.46219598130d-2,
     &             2.72728851066d-2,2.15360151576d-2,1.70426018226d-2,
     &             1.35127158151d-2,1.07324867217d-2,8.53760219545d-3,
     &             6.80120891232d-3,5.42495329856d-3,4.33228989306d-3/ 

         call dylocate(xa,21,x,j)
         if (j.lt.2) j=2
         if (j.gt.20) j=20
         do i=1,4
         xx(i)=xa(j+i-2)
         yy(i)=ya(j+i-2)
         enddo
         call xpolint(xx,yy,4,x,y,dy)
	    zbk2=y
         if(x.lt.1d0) zbk2=0d0

       end

          function zbk3(x)
      implicit none
      integer i,j
      real*8 xa(1:21),ya(1:21),xx(1:5),yy(1:5),zbk3,x,y,dy   
	   DATA XA/1d0,1.2d0,1.4d0,1.6d0,1.8d0,2d0,2.2d0,2.4d0,
     &             2.6d0,2.8d0,3d0,3.2d0,3.4d0,3.6d0,3.8d0,4d0,
     &             4.2d0,4.4d0,4.6d0,4.8d0,5d0/

           DATA YA/0.923433129276d0,0.604688213722d0,0.408298815558d0,
     &             0.282121092418d0,0.198476215508d0,0.141664802216d0,
     &             1.02324263482d-1,7.46480617096d-2,5.49203947626d-2,
     &             4.07017495249d-2,3.03561460004d-2,2.27666860711d-2,
     &             1.71591652247d-2,1.29898581976d-2,9.87254609222d-3,
     &             7.53015239679d-3,5.76215402177d-3,4.42230076800d-3,
     &             3.40318600534d-3,2.62543960976d-3,2.03008127036d-3/

         call dylocate(xa,21,x,j)
         if (j.lt.2) j=2
         if (j.gt.20) j=20
         do i=1,4
         xx(i)=xa(j+i-2)
         yy(i)=ya(j+i-2)
         enddo
         call xpolint(xx,yy,4,x,y,dy)
	    zbk3=y
         if(x.lt.1d0) zbk3=0d0

       end




       SUBROUTINE DYLOCATE(XX,N,X,J)
       INTEGER J,N,JL,JM,JU
       DOUBLE PRECISION X,XX(N)
       JL=0
       JU=N+1
 10    IF(JU-JL.GT.1)THEN
              JM=(JU+JL)/2
              IF((XX(N).GE.XX(1)).EQV.(X.GE.XX(JM))) then
                  JL=JM
              ELSE
                 JU=JM
              ENDIF
       GOTO 10
       ENDIF
       IF(X.EQ.XX(1))THEN
           J=1
       ELSE IF(X.EQ.XX(N)) THEN
           J=N-1 
       ELSE 
              J=JL
       ENDIF
       RETURN
       END


       SUBROUTINE XPOLINT (XA,YA,N,X,Y,DY)
       IMPLICIT DOUBLE PRECISION (A-H, O-Z)
       PARAMETER (NMAX=10)
       DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
       NS=1
       DIF=ABS(X-XA(1))
       DO 11 I=1,N
       DIFT=ABS(X-XA(I))
       IF (DIFT.LT.DIF) THEN
         NS=I
         DIF=DIFT
       ENDIF
       C(I)=YA(I)
       D(I)=YA(I)
 11    CONTINUE
       Y=YA(NS)
       NS=NS-1
       DO 13 M=1,N-1
       DO 12 I=1,N-M
         HO=XA(I)-X
         HP=XA(I+M)-X
         W=C(I+1)-D(I)
         DEN=HO-HP
         IF(DEN.EQ.0.)PAUSE
         DEN=W/DEN
         D(I)=HP*DEN
         C(I)=HO*DEN
 12    CONTINUE
       IF (2*NS.LT.N-M)THEN
         DY=C(NS+1)
       ELSE
         DY=D(NS)
         NS=NS-1
       ENDIF
       Y=Y+DY
 13    CONTINUE
       RETURN
       END
