      function Li2(x)
      implicit none
      double precision x,li2,ddilog,fli2
      external ddilog,fli2

      Li2 = ddilog(x)
      Li2 = fli2(x)
      
      return
      end

c     Interfaces to chaplin harmonic polylogarithmic functions
      
!c     Dilogarithmic function Li2(x)
!      function Li2(x)
!      implicit none
!      double precision x,li2
!      double complex z,HPL2
!      external HPL2
!
!      z = x
!      Li2 = real(hpl2(0,1,z))
!      
!      return
!      end

!c     Polylogarithmic function Li3(x)
!      function Li3(x)
!      implicit none
!      double precision x,li3
!      double complex z,HPL3
!      external HPL3
!
!      z = x
!      Li3 = real(hpl3(0,0,1,x))
!      
!      return
!      end
