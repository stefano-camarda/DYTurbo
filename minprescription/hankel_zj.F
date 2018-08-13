      subroutine ch12n ( n, z, nm, chf1, chd1, chf2, chd2 )

c*********************************************************************72
c
cc CH12N computes Hankel functions of first and second kinds, complex argument.
c
c  Discussion:
c
c    Both the Hankel functions and their derivatives are computed.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    26 July 2012
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
c    Input, integer N, the order of the functions.
c
c    Input, complex*16 Z, the argument.
c
c    Output, integer NM, the highest order computed.
c
c    Output, complex*16 CHF1(0:n), CHD1(0:n), CHF2(0:n), CHD2(0:n), the
c    values of Hn(1)(z), Hn(1)'(z), Hn(2)(z), Hn(2)'(z).
c
      implicit none

      integer n

      complex*16 cbi(0:250)
      complex*16 cbj(0:250)
      complex*16 cbk(0:250)
      complex*16 cby(0:250)
      complex*16 cdi(0:250)
      complex*16 cdj(0:250)
      complex*16 cdk(0:250)
      complex*16 cdy(0:250)
      complex*16 chd1(0:n)
      complex*16 chd2(0:n)
      complex*16 cf1
      complex*16 cfac
      complex*16 chf1(0:n)
      complex*16 chf2(0:n)
      complex*16 ci
      integer k
      integer nm
      double precision pi
      complex*16 z
      complex*16 zi

      ci = cmplx ( 0.0D+00, 1.0D+00 )
      pi = 3.141592653589793D+00

      if ( dimag ( z ) .lt. 0.0D+00 ) then

        call cjynb ( n, z, nm, cbj, cdj, cby, cdy )

        do k = 0, nm
          chf1(k) = cbj(k) + ci * cby(k)
          chd1(k) = cdj(k) + ci * cdy(k)
        end do

        zi = ci * z
        call ciknb ( n, zi, nm, cbi, cdi, cbk, cdk )
        cfac = -2.0D+00 / ( pi * ci )

        do k = 0, nm
          chf2(k) = cfac * cbk(k)
          chd2(k) = cfac * ci * cdk(k)
          cfac = cfac * ci
        end do

      else if ( 0.0D+00 .lt. dimag ( z ) ) then

        zi = - ci * z
        call ciknb ( n, zi, nm, cbi, cdi, cbk, cdk )
        cf1 = -ci
        cfac = 2.0D+00 / ( pi * ci )

        do k = 0, nm
          chf1(k) = cfac * cbk(k)
          chd1(k) = -cfac * ci * cdk(k)
          cfac = cfac * cf1
        end do

        call cjynb ( n, z, nm, cbj, cdj, cby, cdy )

        do k = 0, nm
          chf2(k) = cbj(k) - ci * cby(k)
          chd2(k) = cdj(k) - ci * cdy(k)
        end do

      else

        call cjynb ( n, z, nm, cbj, cdj, cby, cdy )

        do k = 0, nm
          chf1(k) = cbj(k) + ci * cby(k)
          chd1(k) = cdj(k) + ci * cdy(k)
          chf2(k) = cbj(k) - ci * cby(k)
          chd2(k) = cdj(k) - ci * cdy(k)
        end do

      end if

      return
      end
      subroutine cjynb ( n, z, nm, cbj, cdj, cby, cdy )

c*********************************************************************72
c
cc CJYNB: Bessel functions, derivatives, Jn(z) and Yn(z) of complex argument.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    03 August 2012
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
c    Input, integer N, the order of Jn(z) and Yn(z).
c
c    Input, complex*16 Z, the argument of Jn(z) and Yn(z).
c
c    Output, integer NM, the highest order computed.
c
c    Output, complex*16 CBJ(0:N), CDJ(0:N), CBY(0:N), CDY(0:N), 
c    the values of Jn(z), Jn'(z), Yn(z), Yn'(z).
c
      implicit none

      integer n

      double precision a(4)
      double precision a0
      double precision a1(4)
      double precision b(4)
      double precision b1(4)
      complex*16 cbj(0:n)
      complex*16 cbj0
      complex*16 cbj1
      complex*16 cbjk
      complex*16 cbs
      complex*16 cby(0:n)
      complex*16 cby0
      complex*16 cby1
      complex*16 cdj(0:n)
      complex*16 cdy(0:n)
      complex*16 ce
      complex*16 cf
      complex*16 cf1
      complex*16 cf2
      complex*16 cp0
      complex*16 cp1
      complex*16 cq0
      complex*16 cq1
      complex*16 cs0
      complex*16 csu
      complex*16 csv
      complex*16 ct1
      complex*16 ct2
      complex*16 cu
      complex*16 cyy
      double precision el
      integer k
      integer m
      integer msta1
      integer msta2
      integer nm
      double precision pi
      double precision r2p
      double precision y0
      complex*16 z

      save a
      save a1
      save b
      save b1

      data a    / -0.7031250000000000D-01, 0.1121520996093750D+00,
     &            -0.5725014209747314D+00, 0.6074042001273483D+01/

      data a1    / 0.1171875000000000D+00,-0.1441955566406250D+00,
     &             0.6765925884246826D+00,-0.6883914268109947D+01/

      data b     / 0.7324218750000000D-01,-0.2271080017089844D+00,
     &             0.1727727502584457D+01,-0.2438052969955606D+02/

      data b1    / -0.1025390625000000D+00,0.2775764465332031D+00,
     &             -0.1993531733751297D+01,0.2724882731126854D+02/

      el = 0.5772156649015329D+00
      pi = 3.141592653589793D+00
      r2p = 0.63661977236758D+00
      y0 = abs ( dimag ( z ) )
      a0 = cdabs ( z )
      nm = n

      if ( a0 .lt. 1.0D-100 ) then
        do k = 0, n
          cbj(k) = cmplx ( 0.0D+00, 0.0D+00 )
          cdj(k) = cmplx ( 0.0D+00, 0.0D+00 )
          cby(k) = - cmplx ( 1.0D+30, 0.0D+00 )
          cdy(k) = cmplx ( 1.0D+30, 0.0D+00 )
        end do
        cbj(0) = cmplx ( 1.0D+00, 0.0D+00 )
        cdj(1) = cmplx ( 0.5D+00, 0.0D+00 )
        return
      end if

      if ( a0 .le. 300.0D+00 .or. 80 .lt. n ) then

        if ( n .eq. 0 ) then
          nm = 1
        end if
        m = msta1 ( a0, 200 )
        if ( m .lt. nm ) then
          nm = m
        else
          m = msta2 ( a0, nm, 15 )
        end if

        cbs = cmplx ( 0.0D+00, 0.0D+00 )
        csu = cmplx ( 0.0D+00, 0.0D+00 )
        csv = cmplx ( 0.0D+00, 0.0D+00 )
        cf2 = cmplx ( 0.0D+00, 0.0D+00 )
        cf1 = cmplx ( 1.0D-30, 0.0D+00 )

        do k = m, 0, -1
          cf = 2.0D+00 * ( k + 1.0D+00 ) / z * cf1 - cf2
          if ( k .le. nm ) then
            cbj(k) = cf
          end if
          if ( k .eq. 2 * int ( k / 2 ) .and. k .ne. 0 ) then
            if ( y0 .le. 1.0D+00 ) then
              cbs = cbs + 2.0D+00 * cf
            else
              cbs = cbs + ( -1.0D+00 ) ** ( k / 2 ) * 2.0D+00 * cf
            end if
            csu = csu + ( -1.0D+00 ) ** ( k / 2 ) * cf / k
          else if ( 1 .lt. k ) then
            csv = csv + ( -1.0D+00 ) ** ( k / 2 ) * k 
     &        / ( k * k - 1.0D+00 ) * cf
          end if
          cf2 = cf1
          cf1 = cf
        end do

        if ( y0 .le. 1.0D+00 ) then
          cs0 = cbs + cf
        else
          cs0 = ( cbs + cf ) / cdcos ( z )
        end if

        do k = 0, nm
          cbj(k) = cbj(k) / cs0
        end do

        ce = cdlog ( z / 2.0D+00 ) + el
        cby(0) = r2p * ( ce * cbj(0) - 4.0D+00 * csu / cs0 )
        cby(1) = r2p * ( - cbj(0) / z + ( ce - 1.0D+00 ) * cbj(1) 
     &    - 4.0D+00 * csv / cs0 )

      else

        ct1 = z - 0.25D+00 * pi
        cp0 = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, 4
          cp0 = cp0 + a(k) * z ** ( - 2 * k )
        end do
        cq0 = -0.125D+00 / z
        do k = 1, 4
          cq0 = cq0 + b(k) * z ** ( - 2 * k - 1 )
        end do
        cu = cdsqrt ( r2p / z )
        cbj0 = cu * ( cp0 * cdcos ( ct1 ) - cq0 * cdsin ( ct1 ) )
        cby0 = cu * ( cp0 * cdsin ( ct1 ) + cq0 * cdcos ( ct1 ) )
        cbj(0) = cbj0
        cby(0) = cby0
        ct2 = z - 0.75D+00 * pi
        cp1 = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, 4
          cp1 = cp1 + a1(k) * z ** ( - 2 * k )
        end do
        cq1 = 0.375D+00 / z
        do k = 1, 4
          cq1 = cq1 + b1(k) * z ** ( - 2 * k - 1 )
        end do
        cbj1 = cu * ( cp1 * cdcos ( ct2 ) - cq1 * cdsin ( ct2 ) )
        cby1 = cu * ( cp1 * cdsin ( ct2 ) + cq1 * cdcos ( ct2 ) )
        cbj(1) = cbj1
        cby(1) = cby1
        do k = 2, nm
          cbjk = 2.0D+00 * ( k - 1.0D+00 ) / z * cbj1 - cbj0
          cbj(k) = cbjk
          cbj0 = cbj1
          cbj1 = cbjk
        end do
      end if

      cdj(0) = -cbj(1)
      do k = 1, nm
        cdj(k) = cbj(k-1) - k / z * cbj(k)
      end do

      if ( 1.0D+00 .lt. cdabs ( cbj(0) ) ) then
        cby(1) = ( cbj(1) * cby(0) - 2.0D+00 / ( pi * z ) ) / cbj(0)
      end if

      do k = 2, nm
        if ( cdabs ( cbj(k-2) ) .le. cdabs ( cbj(k-1) ) ) then
          cyy = ( cbj(k) * cby(k-1) - 2.0D+00 / ( pi * z ) ) / cbj(k-1)
        else
          cyy = ( cbj(k) * cby(k-2) - 4.0D+00 * ( k - 1.0D+00 )
     &      / ( pi * z * z ) ) / cbj(k-2)
        end if
        cby(k) = cyy
      end do

      cdy(0) = -cby(1)
      do k = 1, nm
        cdy(k) = cby(k-1) - k / z * cby(k)
      end do

      return
      end
      subroutine ciknb ( n, z, nm, cbi, cdi, cbk, cdk )

c*********************************************************************72
c
cc CIKNB computes modified Bessel functions In(z) and Kn(z) for complex argument.
c
c  Discussion:
c
c    This procedure also evaluates the derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    30 July 2012
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
c    Input, integer N, the order of In(z) and Kn(z).
c
c    Input, complex*16 Z, the argument.
c
c    Output, integer NM, the highest order computed.
c
c    Output, complex*16 CBI((0:N), CDI(0:N), CBK(0:N), CDK(0:N), the values of
c    In(z), In'(z), Kn(z), Kn'(z).
c
      implicit none

      integer n

      double precision a0
      complex*16 c
      complex*16 ca0
      complex*16 cbi(0:n)
      complex*16 cbkl
      complex*16 cbs
      complex*16 cdi(0:n)
      complex*16 cbk(0:n)
      complex*16 cdk(0:n)
      complex*16 cf
      complex*16 cf0
      complex*16 cf1
      complex*16 cg
      complex*16 cg0
      complex*16 cg1
      complex*16 ci
      complex*16 cr
      complex*16 cs0
      complex*16 csk0
      double precision el
      double precision fac
      integer k
      integer k0
      integer l
      integer m
      integer msta1
      integer msta2
      integer nm
      double precision pi
      double precision vt
      complex*16 z
      complex*16 z1

      pi = 3.141592653589793D+00
      el = 0.57721566490153D+00
      a0 = cdabs ( z )
      nm = n

      if ( a0 .lt. 1.0D-100 ) then
        do k = 0, n
          cbi(k) = cmplx ( 0.0D+00, 0.0D+00 )
          cbk(k) = cmplx ( 1.0D+30, 0.0D+00 )
          cdi(k) = cmplx ( 0.0D+00, 0.0D+00 )
          cdk(k) = - cmplx ( 1.0D+30, 0.0D+00 )
        end do
        cbi(0) = cmplx ( 1.0D+00, 0.0D+00 )
        cdi(1) = cmplx ( 0.5D+00, 0.0D+00 ) 
        return
      end if

      ci = cmplx ( 0.0D+00, 1.0D+00 )

      if ( real ( z ) .lt. 0.0D+00 ) then
        z1 = -z
      else
        z1 = z
      end if

      if ( n .eq. 0 ) then
        nm = 1
      end if

      m = msta1 ( a0, 200 )

      if ( m .lt. nm ) then
        nm = m
      else
        m = msta2 ( a0, nm, 15 )
      end if

      cbs = 0.0D+00
      csk0 = 0.0D+00
      cf0 = 0.0D+00
      cf1 = 1.0D-100

      do k = m, 0, -1
        cf = 2.0D+00 * ( k + 1.0D+00 ) * cf1 / z1 + cf0
        if ( k .le. nm ) then
          cbi(k) = cf
        end if
        if ( k .ne. 0 .and. k .eq. 2 * int ( k / 2 ) ) then
          csk0 = csk0 + 4.0D+00 * cf / k
        end if
        cbs = cbs + 2.0D+00 * cf
        cf0 = cf1
        cf1 = cf
      end do

      cs0 = cdexp ( z1 ) / ( cbs - cf )

      do k = 0, nm
        cbi(k) = cs0 * cbi(k)
      end do

      if ( a0 .le. 9.0D+00 ) then

        cbk(0) = - ( cdlog ( 0.5D+00 * z1 ) + el ) * cbi(0) + cs0 * csk0
        cbk(1) = ( 1.0D+00 / z1 - cbi(1) * cbk(0) ) / cbi(0)

      else

        ca0 = cdsqrt ( pi / ( 2.0D+00 * z1 ) ) * cdexp ( -z1 )

        if ( a0 .lt. 25.0D+00 ) then
          k0 = 16
        else if ( a0 .lt. 80.0D+00 ) then
          k0 = 10
        else if ( a0 .lt. 200.0D+00 ) then
          k0 = 8
        else
          k0 = 6
        end if

        do l = 0, 1
          cbkl = 1.0D+00
          vt = 4.0D+00 * l
          cr = cmplx ( 1.0D+00, 0.0D+00 )
          do k = 1, k0
            cr = 0.125D+00 * cr 
     &        * ( vt - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) / ( k * z1 )
            cbkl = cbkl + cr
          end do
          cbk(l) = ca0 * cbkl
        end do
      end if

      cg0 = cbk(0)
      cg1 = cbk(1)
      do k = 2, nm
        cg = 2.0D+00 * ( k - 1.0D+00 ) / z1 * cg1 + cg0
        cbk(k) = cg
        cg0 = cg1
        cg1 = cg
      end do

      if ( real ( z ) .lt. 0.0D+00 ) then
        fac = 1.0D+00
        do k = 0, nm
          if ( dimag ( z ) .lt. 0.0D+00 ) then
            cbk(k) = fac * cbk(k) + ci * pi * cbi(k)
          else
            cbk(k) = fac * cbk(k) - ci * pi * cbi(k)
          end if
          cbi(k) = fac * cbi(k)
          fac = - fac
        end do
      end if

      cdi(0) = cbi(1)
      cdk(0) = -cbk(1)
      do k = 1, nm
        cdi(k) = cbi(k-1) - k / z * cbi(k)
        cdk(k) = - cbk(k-1) - k / z * cbk(k)
      end do

      return
      end
      function msta1 ( x, mp )

c*********************************************************************72
c
cc MSTA1 determines a backward recurrence starting point for Jn(x).
c
c  Discussion:
c
c    This procedure determines the starting point for backward  
c    recurrence such that the magnitude of    
c    Jn(x) at that point is about 10^(-MP).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    08 July 2012
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
c    Input, double precision X, the argument.
c
c    Input, integer MP, the negative logarithm of the desired magnitude.
c
c    Output, integer MSTA1, the starting point.
c
      implicit none

      double precision a0
      double precision envj
      double precision f
      double precision f0
      double precision f1
      integer it
      integer mp
      integer msta1
      integer n0
      integer n1
      integer nn
      double precision x

      a0 = abs ( x )
      n0 = int ( 1.1D+00 * a0 ) + 1
      f0 = envj ( n0, a0 ) - mp
      n1 = n0 + 5
      f1 = envj ( n1, a0 ) - mp
      do it = 1, 20       
        nn = n1 - ( n1 - n0 ) / ( 1.0D+00 - f0 / f1 )                  
        f = envj ( nn, a0 ) - mp
        if ( abs ( nn - n1 ) .lt. 1 ) then
          go to 10
        end if
        n0 = n1
        f0 = f1
        n1 = nn
        f1 = f
      end do

 10   continue

      msta1 = nn

      return
      end
      function msta2 ( x, n, mp )

c*********************************************************************72
c
cc MSTA2 determines a backward recurrence starting point for Jn(x).
c
c  Discussion:
c
c    This procedure determines the starting point for a backward
c    recurrence such that Jn(x) has MP significant digits.
c
c    Jianming Jin supplied a modification to this code on 12 January 2016.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    14 January 2016
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
c    Input, double precision X, the argument of Jn(x).
c
c    Input, integer N, the order of Jn(x).
c
c    Input, integer MP, the number of significant digits.
c
c    Output, integer MSTA2, the starting point.
c
      implicit none

      double precision a0
      double precision ejn
      double precision envj
      double precision f
      double precision f0
      double precision f1
      double precision hmp
      integer it
      integer mp
      integer msta2
      integer n
      integer n0
      integer n1
      integer nn
      double precision obj
      double precision x

      a0 = abs ( x )
      hmp = 0.5D+00 * mp
      ejn = envj ( n, a0 )

      if ( ejn .le. hmp ) then
        obj = mp
c
c  Original code:
c
c       n0 = int ( 1.1D+00 * a0 )
c
c  Updated code:
c
        n0 = int ( 1.1D+00 * a0 ) + 1
      else
        obj = hmp + ejn
        n0 = n
      end if

      f0 = envj ( n0, a0 ) - obj
      n1 = n0 + 5
      f1 = envj ( n1, a0 ) - obj

      do it = 1, 20
        nn = n1 - ( n1 - n0 ) / ( 1.0D+00 - f0 / f1 )
        f = envj ( nn, a0 ) - obj
        if ( abs ( nn - n1 ) .lt. 1 ) then
          go to 10
        end if
        n0 = n1
        f0 = f1
        n1 = nn
        f1 = f
      end do

10    continue

      msta2 = nn + 10

      return
      end
      function envj ( n, x )

c*********************************************************************72
c
cc ENVJ is a utility function used by MSTA1 and MSTA2.
c
c  Discussion:
c
c    ENVJ estimates -log(Jn(x)) from the estimate
c    Jn(x) approx 1/sqrt(2*pi*n) * ( e*x/(2*n))^n
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    14 January 2016
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c    Modifications suggested by Vincent Lafage, 11 January 2016.
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
c    Input, integer N, the order of the Bessel function.
c
c    Input, double precision X, the absolute value of the argument.
c
c    Output, double precision ENVJ, the value.
c
      implicit none

      double precision envj
      double precision logten
      integer n
      double precision n_r8
      double precision r8_gamma_log
      double precision x
c
c  Original code
c
      if ( .true. ) then

        envj = 0.5D+00 * log10 ( 6.28D+00 * n ) 
     &    - n * log10 ( 1.36D+00 * x / n )
c
c  Modification suggested by Vincent Lafage.
c
      else

        n_r8 = dble ( n )
        logten = log ( 10.0D+00 )

        envj = r8_gamma_log ( n_r8 + 1.0D+00 ) / logten 
     &    - n_r8 * log10 ( x )

      end if

      return
      end
