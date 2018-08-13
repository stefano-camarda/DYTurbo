      subroutine jyndd ( n, x, bjn, djn, fjn, byn, dyn, fyn )

c*********************************************************************72
c
cc JYNDD computes Bessel functions Jn(x) and Yn(x), first and second derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    07 April 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the order.
c
c    Input, double precision X, the argument.
c
c    Output, double precision BJN, DJN, FJN, BYN, DYN, FYN, the values of
c    Jn(x), Jn'(x), Jn"(x), Yn(x), Yn'(x), Yn"(x).
c
      implicit none

      double precision bj(102)
      double precision bjn
      double precision byn
      double precision bs
      double precision by(102)
      double precision djn
      double precision dyn
      double precision e0
      double precision ec
      double precision f
      double precision f0
      double precision f1
      double precision fjn
      double precision fyn
      integer k
      integer m
      integer mt
      integer n
      integer nt
      double precision s1
      double precision su
      double precision x

      do nt = 1, 900
        mt = int ( 0.5D+00 * log10 ( 6.28D+00 * nt )
     &    - nt * log10 ( 1.36D+00 * abs ( x ) / nt ) )
        if ( 20 .lt. mt ) then
          go to 10
        end if
      end do

10    continue

      m = nt
      bs = 0.0D+00
      f0 = 0.0D+00
      f1 = 1.0D-35
      su = 0.0D+00
      do k = m, 0, -1
        f = 2.0D+00 * ( k + 1.0D+00 ) * f1 / x - f0
        if ( k .le. n + 1 ) then
          bj(k+1) = f
        end if
        if ( k .eq. 2 * int ( k / 2 ) ) then
          bs = bs + 2.0D+00 * f
          if ( k .ne. 0 ) then
            su = su + ( -1.0D+00 ) ** ( k / 2 ) * f / k
          end if
        end if
        f0 = f1
        f1 = f
      end do

      do k = 0, n + 1
        bj(k+1) = bj(k+1) / ( bs - f )
      end do

      bjn = bj(n+1)
      ec = 0.5772156649015329D+00
      e0 = 0.3183098861837907D+00
      s1 = 2.0D+00 * e0 * ( log ( x / 2.0D+00 ) + ec ) * bj(1)
      f0 = s1 - 8.0D+00 * e0 * su / ( bs - f )
      f1 = ( bj(2) * f0 - 2.0D+00 * e0 / x ) / bj(1)

      by(1) = f0
      by(2) = f1
      do k = 2, n + 1 
        f = 2.0D+00 * ( k - 1.0D+00 ) * f1 / x - f0
        by(k+1) = f
        f0 = f1
        f1 = f
      end do

      byn = by(n+1)
      djn = - bj(n+2) + n * bj(n+1) / x
      dyn = - by(n+2) + n * by(n+1) / x
      fjn = ( n * n / ( x * x ) - 1.0D+00 ) * bjn - djn / x
      fyn = ( n * n / ( x * x ) - 1.0D+00 ) * byn - dyn / x

      return
      end
      subroutine jyzo ( n, nt, rj0, rj1, ry0, ry1 )

c*********************************************************************72
c
cc JYZO computes the zeros of Bessel functions Jn(x), Yn(x) and derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    28 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the order of the Bessel functions.
c
c    Input, integer NT, the number of zeros.
c
c    Output, double precision RJ0(NT), RJ1(NT), RY0(NT), RY1(NT), the zeros 
c    of Jn(x), Jn'(x), Yn(x), Yn'(x).
c
      implicit none

      integer nt

      double precision bjn
      double precision byn
      double precision djn
      double precision dyn
      double precision fjn
      double precision fyn
      integer l
      integer n
      double precision rj0(nt)
      double precision rj1(nt)
      double precision ry0(nt)
      double precision ry1(nt)
      double precision x
      double precision x0

      if ( n .le. 20 ) then
        x = 2.82141D+00 + 1.15859D+00 * dble ( n ) 
      else
        x = n + 1.85576D+00 * dble ( n ) ** 0.33333D+00 
     &    + 1.03315D+00 / dble ( n ) ** 0.33333D+00
      end if

      l = 0

10    continue

      x0 = x
      call jyndd ( n, x, bjn, djn, fjn, byn, dyn, fyn )
      x = x - bjn / djn
      if ( 1.0D-09 .lt. abs ( x - x0 ) ) then
        go to 10
      end if
      l = l + 1
      rj0(l) = x
      x = x + 3.1416D+00 + ( 0.0972D+00 + 0.0679D+00 * dble ( n ) 
     &  - 0.000354D+00 * dble ( n ) ** 2 ) / l

      if ( l .lt. nt ) then
        go to 10
      end if

      if ( n .le. 20 ) then
        x = 0.961587D+00 + 1.07703D+00 * dble ( n ) 
      else
        x =  dble ( n ) + 0.80861D+00 * dble ( n ) ** 0.33333D+00
     &    + 0.07249D+00 / dble ( n ) ** 0.33333D+00
      end if

      if ( n .eq. 0 ) then
        x = 3.8317D+00
      end if

      l = 0

20    continue

      x0 = x
      call jyndd ( n, x, bjn, djn, fjn, byn, dyn, fyn )
      x = x - djn / fjn
      if ( 1.0D-09 .lt. abs ( x - x0 ) ) then
        go to 20
      end if
      l = l + 1
      rj1(l) = x
      x = x + 3.1416D+00 + ( 0.4955D+00 + 0.0915D+00 * dble ( n ) 
     &  - 0.000435D+00 * dble ( n ) ** 2 ) / l

      if ( l .lt. nt ) then
        go to 20
      end if

      if ( n .le. 20 ) then
        x = 1.19477D+00 + 1.08933D+00 * dble ( n ) 
      else
        x =  dble ( n ) + 0.93158D+00 * dble ( n ) ** 0.33333D+00
     &    + 0.26035D+00 / dble ( n ) ** 0.33333D+00
      end if
 
      l = 0

30    continue

      x0 = x
      call jyndd ( n, x, bjn, djn, fjn, byn, dyn, fyn )
      x = x - byn / dyn

      if ( 1.0D-09 .lt. abs ( x - x0 ) ) then
        go to 30
      end if

      l = l + 1
      ry0(l) = x
      x = x + 3.1416D+00 + ( 0.312D+00 + 0.0852D+00 * dble ( n ) 
     &  - 0.000403D+00 * dble ( n ) ** 2 ) / l

      if ( l .lt. nt ) then
        go to 30
      end if

      if ( n .le. 20 ) then
        x = 2.67257D+00 + 1.16099D+00 * dble ( n ) 
      else
        x =  dble ( n ) + 1.8211D+00 * dble ( n ) ** 0.33333D+00
     &    + 0.94001D+00 / dble ( n ) ** 0.33333D+00
      end if
  
      l = 0

40    continue

      x0 = x
      call jyndd ( n, x, bjn, djn, fjn, byn, dyn, fyn )
      x = x - dyn / fyn
      if ( 1.0D-09 .lt. abs ( x - x0 ) ) then
        go to 40
      end if
      l = l + 1
      ry1(l) = x
      x = x + 3.1416D+00 + ( 0.197D+00 + 0.0643D+00 * dble ( n ) 
     &  -0.000286D+00 * dble ( n ) ** 2 ) / l 

      if ( l .lt. nt ) then
        go to 40
      end if

      return
      end
