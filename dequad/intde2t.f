!     test of intde2.f
!
      include 'intde2.f'
!
      program main
      implicit real*8 (a - i, o - z)
      parameter (lenaw = 8000)
      parameter (tiny = 1.0d-307)
      dimension aw(0 : lenaw - 1)
      external f1, f2, f3, f4, f5, f6
      common  / nfunc / nn
!     
      call intdeini(lenaw, tiny, 1.0d-15, aw)
      nn = 0
      call intde(f1, 0.0d0, 1.0d0, aw, i, err)
      write (*, *) 'I_1=int_0^1 1/sqrt(x) dx'
      write (*, *) '  I_1=', i, ',  err=', err, ',  n=', nn
      nn = 0
      call intde(f2, 0.0d0, 2.0d0, aw, i, err)
      write (*, *) 'I_2=int_0^2 sqrt(4-x*x) dx'
      write (*, *) '  I_2=', i, ',  err=', err, ',  n=', nn
!     
      call intdeiini(lenaw, tiny, 1.0d-15, aw)
      nn = 0
      call intdei(f3, 0.0d0, aw, i, err)
      write (*, *) 'I_3=int_0^infty 1/(1+x*x) dx'
      write (*, *) '  I_3=', i, ',  err=', err, ',  n=', nn
      nn = 0
      call intdei(f4, 0.0d0, aw, i, err)
      write (*, *) 'I_4=int_0^infty exp(-x)/sqrt(x) dx'
      write (*, *) '  I_4=', i, ',  err=', err, ',  n=', nn
!     
      call intdeoini(lenaw, tiny, 1.0d-15, aw)
      nn = 0
      call intdeo(f5, 0.0d0, 1.0d0, aw, i, err)
      write (*, *) 'I_5=int_0^infty sin(x)/x dx'
      write (*, *) '  I_5=', i, ',  err=', err, ',  n=', nn
      nn = 0
      call intdeo(f6, 0.0d0, 1.0d0, aw, i, err)
      write (*, *) 'I_6=int_0^infty cos(x)/sqrt(x) dx'
      write (*, *) '  I_6=', i, ',  err=', err, ',  n=', nn
      end
!
!
      function f1(x)
      implicit real*8 (a - i, o - z)
      common  / nfunc / n
      n = n + 1
      f1 = 1 / sqrt(x)
      end
!
!
      function f2(x)
      implicit real*8 (a - i, o - z)
      common  / nfunc / n
      n = n + 1
      f2 = sqrt(4 - x * x)
      end
!
!
      function f3(x)
      implicit real*8 (a - i, o - z)
      common  / nfunc / n
      n = n + 1
      f3 = 1 / (1 + x * x)
      end
!
!
      function f4(x)
      implicit real*8 (a - i, o - z)
      common  / nfunc / n
      n = n + 1
      f4 = exp(-x) / sqrt(x)
      end
!
!
      function f5(x)
      implicit real*8 (a - i, o - z)
      common  / nfunc / n
      n = n + 1
      f5 = sin(x) / x
      end
!
!
      function f6(x)
      implicit real*8 (a - i, o - z)
      common  / nfunc / n
      n = n + 1
      f6 = cos(x) / sqrt(x)
      end
!
!
