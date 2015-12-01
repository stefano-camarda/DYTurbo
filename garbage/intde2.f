! DE-Quadrature
! Numerical Automatic Integrator for Improper Integral
!     method    : Double Exponential (DE) Transformation
!     dimension : one
!     table     : use
! subroutines
!     intde  : integrator of f(x) over (a,b).
!     intdei : integrator of f(x) over (a,infinity), 
!                  f(x) is non oscillatory function.
!     intdeo : integrator of f(x) over (a,infinity), 
!                  f(x) is oscillatory function.
!
!
! intde
!     [description]
!         I = integral of f(x) over (a,b)
!     [declaration]
!         external f
!     [usage]
!         call intdeini(lenaw, tiny, eps, aw)  ! initialization of aw
!         ...
!         call intde(f, a, b, aw, i, err)
!     [parameters]
!         lenaw     : length of aw (integer)
!         tiny      : minimum value that 1/tiny does not 
!                     overflow (real*8)
!         eps       : relative error requested (real*8)
!         aw        : points and weights of the quadrature 
!                     formula, aw(0...lenaw-1) (real*8)
!         f         : integrand f(x) (real*8 function)
!         a         : lower limit of integration (real*8)
!         b         : upper limit of integration (real*8)
!         i         : approximation to the integral (real*8)
!         err       : estimate of the absolute error (real*8)
!     [remarks]
!         initial parameters
!             lenaw > 1000, 
!             IEEE double :
!                 lenaw = 8000
!                 tiny = 1.0d-307
!         function
!             f(x) needs to be analytic over (a,b).
!         relative error
!             eps is relative error requested excluding 
!             cancellation of significant digits.
!             i.e. eps means : (absolute error) / 
!                              (integral_a^b |f(x)| dx).
!             eps does not mean : (absolute error) / I.
!         error message
!             err >= 0 : normal termination.
!             err < 0  : abnormal termination.
!                        i.e. convergent error is detected :
!                            1. f(x) or (d/dx)^n f(x) has 
!                               discontinuous points or sharp 
!                               peaks over (a,b).
!                               you must divide the interval 
!                               (a,b) at this points.
!                            2. relative error of f(x) is 
!                               greater than eps.
!                            3. f(x) has oscillatory factor 
!                               and frequency of the oscillation 
!                               is very high.
!
!
! intdei
!     [description]
!         I = integral of f(x) over (a,infinity), 
!             f(x) has not oscillatory factor.
!     [declaration]
!         external f
!     [usage]
!         call intdeiini(lenaw, tiny, eps, aw)  ! initialization of aw
!         ...
!         call intdei(f, a, aw, i, err)
!     [parameters]
!         lenaw     : length of aw (integer)
!         tiny      : minimum value that 1/tiny does not 
!                     overflow (real*8)
!         eps       : relative error requested (real*8)
!         aw        : points and weights of the quadrature 
!                     formula, aw(0...lenaw-1) (real*8)
!         f         : integrand f(x) (real*8 function)
!         a         : lower limit of integration (real*8)
!         i         : approximation to the integral (real*8)
!         err       : estimate of the absolute error (real*8)
!     [remarks]
!         initial parameters
!             lenaw > 1000, 
!             IEEE double :
!                 lenaw = 8000
!                 tiny = 1.0d-307
!         function
!             f(x) needs to be analytic over (a,infinity).
!         relative error
!             eps is relative error requested excluding 
!             cancellation of significant digits.
!             i.e. eps means : (absolute error) / 
!                              (integral_a^infinity |f(x)| dx).
!             eps does not mean : (absolute error) / I.
!         error message
!             err >= 0 : normal termination.
!             err < 0  : abnormal termination.
!                        i.e. convergent error is detected :
!                            1. f(x) or (d/dx)^n f(x) has 
!                               discontinuous points or sharp 
!                               peaks over (a,infinity).
!                               you must divide the interval 
!                               (a,infinity) at this points.
!                            2. relative error of f(x) is 
!                               greater than eps.
!                            3. f(x) has oscillatory factor 
!                               and decay of f(x) is very slow 
!                               as x -> infinity.
!
!
! intdeo
!     [description]
!         I = integral of f(x) over (a,infinity), 
!             f(x) has oscillatory factor :
!             f(x) = g(x) * sin(omega * x + theta) as x -> infinity.
!     [declaration]
!         external f
!     [usage]
!         call intdeoini(lenaw, tiny, eps, aw)  ! initialization of aw
!         ...
!         call intdeo(f, a, omega, aw, i, err)
!     [parameters]
!         lenaw     : length of aw (integer)
!         tiny      : minimum value that 1/tiny does not 
!                     overflow (real*8)
!         eps       : relative error requested (real*8)
!         aw        : points and weights of the quadrature 
!                     formula, aw(0...lenaw-1) (real*8)
!         f         : integrand f(x) (real*8 function)
!         a         : lower limit of integration (real*8)
!         omega     : frequency of oscillation (real*8)
!         i         : approximation to the integral (real*8)
!         err       : estimate of the absolute error (real*8)
!     [remarks]
!         initial parameters
!             lenaw > 1000, 
!             IEEE double :
!                 lenaw = 8000
!                 tiny = 1.0d-307
!         function
!             f(x) needs to be analytic over (a,infinity).
!         relative error
!             eps is relative error requested excluding 
!             cancellation of significant digits.
!             i.e. eps means : (absolute error) / 
!                              (integral_a^R |f(x)| dx).
!             eps does not mean : (absolute error) / I.
!         error message
!             err >= 0 : normal termination.
!             err < 0  : abnormal termination.
!                        i.e. convergent error is detected :
!                            1. f(x) or (d/dx)^n f(x) has 
!                               discontinuous points or sharp 
!                               peaks over (a,infinity).
!                               you must divide the interval 
!                               (a,infinity) at this points.
!                            2. relative error of f(x) is 
!                               greater than eps.
!
!
      subroutine intdeini(lenaw, tiny, eps, aw)
      integer lenaw
      real*8 tiny, eps, aw(0 : lenaw - 1)
      real*8 efs, hoff
      integer noff, nk, k, j
      real*8 pi2, tinyln, epsln, h0, ehp, ehm, h, t, ep, em, xw, wg
! ---- adjustable parameter ----
      efs = 0.1d0
      hoff = 8.5d0
! ------------------------------
      pi2 = 2 * atan(1.0d0)
      tinyln = -log(tiny)
      epsln = 1 - log(efs * eps)
      h0 = hoff / epsln
      ehp = exp(h0)
      ehm = 1 / ehp
      aw(2) = eps
      aw(3) = exp(-ehm * epsln)
      aw(4) = sqrt(efs * eps)
      noff = 5
      aw(noff) = 0.5d0
      aw(noff + 1) = h0
      aw(noff + 2) = pi2 * h0 * 0.5d0
      h = 2
      nk = 0
      k = noff + 3
   10 continue
          t = h * 0.5d0
   20     continue
              em = exp(h0 * t)
              ep = pi2 * em
              em = pi2 / em
              j = k
   30         continue
                  xw = 1 / (1 + exp(ep - em))
                  wg = xw * (1 - xw) * h0
                  aw(j) = xw
                  aw(j + 1) = wg * 4
                  aw(j + 2) = wg * (ep + em)
                  ep = ep * ehp
                  em = em * ehm
                  j = j + 3
              if (ep .lt. tinyln .and. j .le. lenaw - 3) goto 30
              t = t + h
              k = k + nk
          if (t .lt. 1) goto 20
          h = h * 0.5d0
          if (nk .eq. 0) then
              if (j .gt. lenaw - 6) j = j - 3
              nk = j - noff
              k = k + nk
              aw(1) = nk
          end if
      if (2 * k - noff - 3 .le. lenaw) goto 10
      aw(0) = k - 3
      end
!
      subroutine intde(f, a, b, aw, i, err)
      real*8 f, a, b, aw(0 : *), i, err
      integer noff, lenawm, nk, k, j, jtmp, jm, m, klim
      real*8 epsh, ba, ir, xa, fa, fb, errt, errh, errd, h, iback, 
     &    irback
      noff = 5
      lenawm = int(aw(0) + 0.5d0)
      nk = int(aw(1) + 0.5d0)
      epsh = aw(4)
      ba = b - a
      i = f((a + b) * aw(noff))
      ir = i * aw(noff + 1)
      i = i * aw(noff + 2)
      err = abs(i)
      k = nk + noff
      j = noff
   10 continue
          j = j + 3
          xa = ba * aw(j)
          fa = f(a + xa)
          fb = f(b - xa)
          ir = ir + (fa + fb) * aw(j + 1)
          fa = fa * aw(j + 2)
          fb = fb * aw(j + 2)
          i = i + (fa + fb)
          err = err + (abs(fa) + abs(fb))
      if (aw(j) .gt. epsh .and. j .lt. k) goto 10
      errt = err * aw(3)
      errh = err * epsh
      errd = 1 + 2 * errh
      jtmp = j
      do while (abs(fa) .gt. errt .and. j .lt. k)
          j = j + 3
          fa = f(a + ba * aw(j))
          ir = ir + fa * aw(j + 1)
          fa = fa * aw(j + 2)
          i = i + fa
      end do
      jm = j
      j = jtmp
      do while (abs(fb) .gt. errt .and. j .lt. k)
          j = j + 3
          fb = f(b - ba * aw(j))
          ir = ir + fb * aw(j + 1)
          fb = fb * aw(j + 2)
          i = i + fb
      end do
      if (j .lt. jm) jm = j
      jm = jm - (noff + 3)
      h = 1
      m = 1
      klim = k + nk
      do while (errd .gt. errh .and. klim .le. lenawm)
          iback = i
          irback = ir
   20     continue
              jtmp = k + jm
              do j = k + 3, jtmp, 3
                  xa = ba * aw(j)
                  fa = f(a + xa)
                  fb = f(b - xa)
                  ir = ir + (fa + fb) * aw(j + 1)
                  i = i + (fa + fb) * aw(j + 2)
              end do
              k = k + nk
              j = jtmp
   30         continue
                  j = j + 3
                  fa = f(a + ba * aw(j))
                  ir = ir + fa * aw(j + 1)
                  fa = fa * aw(j + 2)
                  i = i + fa
              if (abs(fa) .gt. errt .and. j .lt. k) goto 30
              j = jtmp
   40         continue
                  j = j + 3
                  fb = f(b - ba * aw(j))
                  ir = ir + fb * aw(j + 1)
                  fb = fb * aw(j + 2)
                  i = i + fb
              if (abs(fb) .gt. errt .and. j .lt. k) goto 40
          if (k .lt. klim) goto 20
          errd = h * (abs(i - 2 * iback) + abs(ir - 2 * irback))
          h = h * 0.5d0
          m = m * 2
          klim = 2 * klim - noff
      end do
      i = i * (h * ba)
      if (errd .gt. errh) then
          err = -errd * (m * abs(ba))
      else
          err = err * aw(2) * (m * abs(ba))
      end if
      end
!
!
      subroutine intdeiini(lenaw, tiny, eps, aw)
      integer lenaw
      real*8 tiny, eps, aw(0 : lenaw - 1)
      real*8 efs, hoff
      integer noff, nk, k, j
      real*8 pi4, tinyln, epsln, h0, ehp, ehm, h, t, ep, em, xp, 
     &    xm, wp, wm
! ---- adjustable parameter ----
      efs = 0.1d0
      hoff = 11.0d0
! ------------------------------
      pi4 = atan(1.0d0)
      tinyln = -log(tiny)
      epsln = 1 - log(efs * eps)
      h0 = hoff / epsln
      ehp = exp(h0)
      ehm = 1 / ehp
      aw(2) = eps
      aw(3) = exp(-ehm * epsln)
      aw(4) = sqrt(efs * eps)
      noff = 5
      aw(noff) = 1
      aw(noff + 1) = 4 * h0
      aw(noff + 2) = 2 * pi4 * h0
      h = 2
      nk = 0
      k = noff + 6
   10 continue
          t = h * 0.5d0
   20     continue
              em = exp(h0 * t)
              ep = pi4 * em
              em = pi4 / em
              j = k
   30         continue
                  xp = exp(ep - em)
                  xm = 1 / xp
                  wp = xp * ((ep + em) * h0)
                  wm = xm * ((ep + em) * h0)
                  aw(j) = xm
                  aw(j + 1) = xp
                  aw(j + 2) = xm * (4 * h0)
                  aw(j + 3) = xp * (4 * h0)
                  aw(j + 4) = wm
                  aw(j + 5) = wp
                  ep = ep * ehp
                  em = em * ehm
                  j = j + 6
              if (ep .lt. tinyln .and. j .le. lenaw - 6) goto 30
              t = t + h
              k = k + nk
          if (t .lt. 1) goto 20
          h = h * 0.5d0
          if (nk .eq. 0) then
              if (j .gt. lenaw - 12) j = j - 6
              nk = j - noff
              k = k + nk
              aw(1) = nk
          end if
      if (2 * k - noff - 6 .le. lenaw) goto 10
      aw(0) = k - 6
      end
!
      subroutine intdei(f, a, aw, i, err)
      real*8 f, a, aw(0 : *), i, err
      integer noff, lenawm, nk, k, j, jtmp, jm, m, klim
      real*8 epsh, ir, fp, fm, errt, errh, errd, h, iback, irback
      noff = 5
      lenawm = int(aw(0) + 0.5d0)
      nk = int(aw(1) + 0.5d0)
      epsh = aw(4)
      i = f(a + aw(noff))
      ir = i * aw(noff + 1)
      i = i * aw(noff + 2)
      err = abs(i)
      k = nk + noff
      j = noff
   10 continue
          j = j + 6
          fm = f(a + aw(j))
          fp = f(a + aw(j + 1))
          ir = ir + (fm * aw(j + 2) + fp * aw(j + 3))
          fm = fm * aw(j + 4)
          fp = fp * aw(j + 5)
          i = i + (fm + fp)
          err = err + (abs(fm) + abs(fp))
      if (aw(j) .gt. epsh .and. j .lt. k) goto 10
      errt = err * aw(3)
      errh = err * epsh
      errd = 1 + 2 * errh
      jtmp = j
      do while (abs(fm) .gt. errt .and. j .lt. k)
          j = j + 6
          fm = f(a + aw(j))
          ir = ir + fm * aw(j + 2)
          fm = fm * aw(j + 4)
          i = i + fm
      end do
      jm = j
      j = jtmp
      do while (abs(fp) .gt. errt .and. j .lt. k)
          j = j + 6
          fp = f(a + aw(j + 1))
          ir = ir + fp * aw(j + 3)
          fp = fp * aw(j + 5)
          i = i + fp
      end do
      if (j .lt. jm) jm = j
      jm = jm - (noff + 6)
      h = 1
      m = 1
      klim = k + nk
      do while (errd .gt. errh .and. klim .le. lenawm)
          iback = i
          irback = ir
   20     continue
              jtmp = k + jm
              do j = k + 6, jtmp, 6
                  fm = f(a + aw(j))
                  fp = f(a + aw(j + 1))
                  ir = ir + (fm * aw(j + 2) + fp * aw(j + 3))
                  i = i + (fm * aw(j + 4) + fp * aw(j + 5))
              end do
              k = k + nk
              j = jtmp
   30         continue
                  j = j + 6
                  fm = f(a + aw(j))
                  ir = ir + fm * aw(j + 2)
                  fm = fm * aw(j + 4)
                  i = i + fm
              if (abs(fm) .gt. errt .and. j .lt. k) goto 30
              j = jtmp
   40         continue
                  j = j + 6
                  fp = f(a + aw(j + 1))
                  ir = ir + fp * aw(j + 3)
                  fp = fp * aw(j + 5)
                  i = i + fp
              if (abs(fp) .gt. errt .and. j .lt. k) goto 40
          if (k .lt. klim) goto 20
          errd = h * (abs(i - 2 * iback) + abs(ir - 2 * irback))
          h = h * 0.5d0
          m = m * 2
          klim = 2 * klim - noff
      end do
      i = i * h
      if (errd .gt. errh) then
          err = -errd * m
      else
          err = err * (aw(2) * m)
      end if
      end
!
!
      subroutine intdeoini(lenaw, tiny, eps, aw)
      integer lenaw
      real*8 tiny, eps, aw(0 : lenaw - 1)
      integer lmax
      real*8 efs, enoff, pqoff, ppoff
      integer noff0, nk0, noff, k, nk, j
      real*8 pi4, tinyln, epsln, frq4, per2, pp, pq, ehp, ehm, h, 
     &    t, ep, em, tk, xw, wg, xa
! ---- adjustable parameter ----
      lmax = 5
      efs = 0.1d0
      enoff = 0.40d0
      pqoff = 2.9d0
      ppoff = -0.72d0
! ------------------------------
      pi4 = atan(1.0d0)
      tinyln = -log(tiny)
      epsln = 1 - log(efs * eps)
      frq4 = 1 / (2 * pi4)
      per2 = 4 * pi4
      pq = pqoff / epsln
      pp = ppoff - log(pq * pq * frq4)
      ehp = exp(2 * pq)
      ehm = 1 / ehp
      aw(3) = lmax
      aw(4) = eps
      aw(5) = sqrt(efs * eps)
      noff0 = 6
      nk0 = 1 + int(enoff * epsln)
      aw(1) = nk0
      noff = 2 * nk0 + noff0
      wg = 0
      xw = 1
      do k = 1, nk0
          wg = wg + xw
          aw(noff - 2 * k) = wg
          aw(noff - 2 * k + 1) = xw
          xw = xw * (nk0 - k) / k
      end do
      wg = per2 / wg
      do k = noff0, noff - 2, 2
          aw(k) = aw(k) * wg
          aw(k + 1) = aw(k + 1) * wg
      end do
      xw = exp(pp - 2 * pi4)
      aw(noff) = sqrt(xw * (per2 * 0.5d0))
      aw(noff + 1) = xw * pq
      aw(noff + 2) = per2 * 0.5d0
      h = 2
      nk = 0
      k = noff + 3
   10 continue
          t = h * 0.5d0
   20     continue
              em = exp(2 * pq * t)
              ep = pi4 * em
              em = pi4 / em
              tk = t
              j = k
   30         continue
                  xw = exp(pp - ep - em)
                  wg = sqrt(frq4 * xw + tk * tk)
                  xa = xw / (tk + wg)
                  wg = (pq * xw * (ep - em) + xa) / wg
                  aw(j) = xa
                  aw(j + 1) = xw * pq
                  aw(j + 2) = wg
                  ep = ep * ehp
                  em = em * ehm
                  tk = tk + 1
                  j = j + 3
              if (ep .lt. tinyln .and. j .le. lenaw - 3) goto 30
              t = t + h
              k = k + nk
          if (t .lt. 1) goto 20
          h = h * 0.5d0
          if (nk .eq. 0) then
              if (j .gt. lenaw - 6) j = j - 3
              nk = j - noff
              k = k + nk
              aw(2) = nk
          end if
      if (2 * k - noff - 3 .le. lenaw) goto 10
      aw(0) = k - 3
      end
!
      subroutine intdeo(f, a, omega, aw, i, err)
      real*8 f, a, omega, aw(0 : *), i, err
      integer lenawm, nk0, noff0, nk, noff, lmax, m, k, j, jm, l
      real*8 eps, per, perw, w02, ir, h, iback, irback, t, tk, 
     &    xa, fm, fp, errh, s0, s1, s2, errd
      lenawm = int(aw(0) + 0.5d0)
      nk0 = int(aw(1) + 0.5d0)
      noff0 = 6
      nk = int(aw(2) + 0.5d0)
      noff = 2 * nk0 + noff0
      lmax = int(aw(3) + 0.5d0)
      eps = aw(4)
      per = 1 / abs(omega)
      w02 = 2 * aw(noff + 2)
      perw = per * w02
      i = f(a + aw(noff) * per)
      ir = i * aw(noff + 1)
      i = i * aw(noff + 2)
      err = abs(i)
      h = 2
      m = 1
      k = noff
   10 continue
          iback = i
          irback = ir
          t = h * 0.5d0
   20     continue
              if (k .eq. noff) then
                  tk = 1
                  k = k + nk
                  j = noff
   30             continue
                      j = j + 3
                      xa = per * aw(j)
                      fm = f(a + xa)
                      fp = f(a + xa + perw * tk)
                      ir = ir + (fm + fp) * aw(j + 1)
                      fm = fm * aw(j + 2)
                      fp = fp * (w02 - aw(j + 2))
                      i = i + (fm + fp)
                      err = err + (abs(fm) + abs(fp))
                      tk = tk + 1
                  if (aw(j) .gt. eps .and. j .lt. k) goto 30
                  errh = err * aw(5)
                  err = err * eps
                  jm = j - noff
              else
                  tk = t
                  do j = k + 3, k + jm, 3
                      xa = per * aw(j)
                      fm = f(a + xa)
                      fp = f(a + xa + perw * tk)
                      ir = ir + (fm + fp) * aw(j + 1)
                      fm = fm * aw(j + 2)
                      fp = fp * (w02 - aw(j + 2))
                      i = i + (fm + fp)
                      tk = tk + 1
                  end do
                  j = k + jm
                  k = k + nk
              end if
              do while (abs(fm) .gt. err .and. j .lt. k)
                  j = j + 3
                  fm = f(a + per * aw(j))
                  ir = ir + fm * aw(j + 1)
                  fm = fm * aw(j + 2)
                  i = i + fm
              end do
              fm = f(a + perw * tk)
              s2 = w02 * fm
              i = i + s2
              if (abs(fp) .gt. err .or. abs(s2) .gt. err) then
                  l = 0
   40             continue
                      l = l + 1
                      s0 = 0
                      s1 = 0
                      s2 = fm * aw(noff0 + 1)
                      do j = noff0 + 2, noff - 2, 2
                          tk = tk + 1
                          fm = f(a + perw * tk)
                          s0 = s0 + fm
                          s1 = s1 + fm * aw(j)
                          s2 = s2 + fm * aw(j + 1)
                      end do
                      if (s2 .le. err .or. l .ge. lmax) goto 50
                      i = i + w02 * s0
                  goto 40
   50             continue
                  i = i + s1
                  if (s2 .gt. err) err = s2
              end if
              t = t + h
          if (t .lt. 1) goto 20
          if (m .eq. 1) then
              errd = 1 + 2 * errh
          else
              errd = h * (abs(i - 2 * iback) + abs(ir - 2 * irback))
          end if
          h = h * 0.5d0
          m = m * 2
      if (errd .gt. errh .and. 2 * k - noff .le. lenawm) goto 10
      i = i * (h * per)
      if (errd .gt. errh) then
          err = -errd * per
      else
          err = err * (per * m * 0.5d0)
      end if
      end
!
