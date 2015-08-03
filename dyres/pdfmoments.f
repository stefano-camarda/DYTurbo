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

c     nodes for gaussian quadrature
      double precision xxx24(24),www24(24)
      data xxx24/-0.0640568928626056,0.0640568928626056,
     .     -0.1911188674736163,0.1911188674736163,
     .     -0.3150426796961634,0.3150426796961634,
     .     -0.4337935076260451,0.4337935076260451,
     .     -0.5454214713888396,0.5454214713888396,
     .     -0.6480936519369755,0.6480936519369755,
     .     -0.7401241915785544,0.7401241915785544,
     .     -0.8200019859739029,0.8200019859739029,
     .     -0.8864155270044011,0.8864155270044011,
     .     -0.9382745520027328,0.9382745520027328,
     .     -0.9747285559713095,0.9747285559713095,
     .     -0.9951872199970213,0.9951872199970213/
      
      data www24/0.1279381953467522,0.1279381953467522,	
     .     0.1258374563468283,0.1258374563468283,
     .     0.1216704729278034,0.1216704729278034,
     .     0.1155056680537256,0.1155056680537256,
     .     0.1074442701159656,0.1074442701159656,
     .     0.0976186521041139,0.0976186521041139,
     .     0.0861901615319533,0.0861901615319533,
     .     0.0733464814110803,0.0733464814110803,
     .     0.0592985849154368,0.0592985849154368,
     .     0.0442774388174198,0.0442774388174198,
     .     0.0285313886289337,0.0285313886289337,
     .     0.0123412297999872,0.0123412297999872/

      integer intervals,i,j
      double precision muf

c     select beam
      if (beam.eq.1) then
         hadron = ih1
      elseif (beam.eq.2) then
         hadron = ih1
      endif
c     factorization scale
      muf=91.1876d0
c     muf=2D0

c     boundaries of integration      
      xa = 1D-10
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
      

      intervals=100
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

