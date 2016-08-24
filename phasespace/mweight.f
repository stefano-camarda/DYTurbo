      subroutine mweight_breitw(x,mmin2,mmax2,rmass,rwidth,m2,jac)
      implicit none
c     Given a number 0 < x < 1 generates a breit-wigner-unweighted mass-squared
c     distribution m2 and multiplies jac by the corresponding jacobian
c     of the change of variable
      double precision x,mmin2,mmax2,rmass,rwidth,m2,jac
      double precision almin,almax,al,tanal

      almin=atan((mmin2-rmass**2)/rmass/rwidth)
      almax=atan((mmax2-rmass**2)/rmass/rwidth)
      al=(almax-almin)*x+almin
      tanal=tan(al)

      m2=rmass**2+rmass*rwidth*tanal
c---- bw=(1d0+tanal**2)*rmass**2*rwidth**2
      jac=jac*(almax-almin)*rmass*rwidth*(1d0+tanal**2)
      return
      end

      subroutine mweight_flat(x,mmin,mmax,m,jac)
      implicit none
      double precision x,mmin,mmax,m,jac

      m=mmin+(mmax-mmin)*x
      jac=jac*(mmax-mmin)
      return
      end
