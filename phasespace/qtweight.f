      subroutine qtweight(x,qtmin,qtmax,qt,jac)
      implicit none
      double precision x,qtmin,qtmax,qt,jac
      double precision qtmin2,qtmax2,qt2
      double precision tiny
      parameter (tiny=1D-5)
      double precision a,b,xx

      qtmin2 = qtmin**2
      qtmax2 = qtmax**2;
      a = 1d0/(1+log(qtmax2/tiny))
      b = 1d0/(1+log(qtmin2/tiny));
      xx = a + (b-a) * x;
      jac = jac*(b-a);
      qt2=tiny*exp(1d0/xx - 1d0);
      jac=jac*qt2/xx**2
      qt=sqrt(qt2)
      jac=jac/qt/2.;
      return
      end

      subroutine qt2weight(x,qtmin2,qtmax2,qt2,jac)
      implicit none
      double precision x,qtmin2,qtmax2,qt2,jac
      double precision tiny
      parameter (tiny=1D-5)
      double precision a,b,xx

      a = 1d0/(1+log(qtmax2/tiny))
      b = 1d0/(1+log(qtmin2/tiny));
      xx = a + (b-a) * x;
      jac = jac*(b-a);
      qt2=tiny*exp(1d0/xx - 1d0);
      jac=jac*qt2/xx**2
      return
      end

      subroutine qtweight_flat(x,qtmin,qtmax,qt,jac)
      implicit none
      double precision x,qtmin,qtmax,qt,jac
      
      qt=qtmin+(qtmax-qtmin)*x
      jac=jac*(qtmax-qtmin)
      return
      end
