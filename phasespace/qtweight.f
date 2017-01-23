      subroutine qtweight(x,qtmin,qtmax,qt,jac)
      implicit none
      double precision x,qtmin,qtmax,qt,jac
      double precision qtmin2,qtmax2,qt2
      double precision tiny
      parameter (tiny=1d-5) !be carefull, if qtmin^2 < tiny the phase space generation is screwed up
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
      parameter (tiny=1d-5)
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

      subroutine qtweight_res(x,qtmin,qtmax,qt,jac)
      implicit none
      double precision x,qtmin,qtmax,qt,jac
      double precision qtmin2,qtmax2,qt2
      double precision tiny,esp
      parameter (tiny=1d-3,esp=0.5d0)
      double precision a,b,xx
      a = log(qtmin/tiny)**esp
      b = log(qtmax/tiny)**esp
      xx= a + (b-a)*x
      jac=jac*(b-a)
      qt = tiny*exp(xx**(1d0/esp));
      jac=jac*qt*(xx**(1d0/esp-1d0))/esp;
      return
      end
      
c     other options
c  double tiny = 1E-5;
c  double a = 1./(1+log(qtmx/tiny));
c  double b = 1./(1+log(qtmn/tiny));
c  double x2 = a + (b-a) * x[1];
c  jac = jac * (b-a);
c  double qt=tiny*exp(1./x2 - 1);
c  jac=jac*qt/pow(x2,2);

c  double exp = 1./3.a;
c  double a = pow(qtmn,exp);
c  double b = pow(qtmx,exp);
c  double x2=a+(b-a)*x[1];
c  jac=jac*(b-a);
c  double qt = pow(x2,1./exp);
c  jac=jac*pow(x2,1./exp-1)/exp;

c  double a = log(qtmn);
c  double b = log(qtmx);
c  double x2=a+(b-a)*x[1];
c  jac=jac*(b-a);
c  double qt = exp(x2);
c  jac=jac*qt;


c  double tiny = 1E-3;
c  double a = log(log(qtmn/tiny));
c  double b = log(log(qtmx/tiny));
c  double x2=a+(b-a)*x[1];
c  jac=jac*(b-a);
c  double qt = tiny*exp(exp(x2));
c  jac=jac*qt*exp(x2);
  
c  double base = 100000.;
c  double a = log(qtmn)/log(base);
c  double b = log(qtmx)/log(base);
c  double x2=a+(b-a)*x[1];
c  jac=jac*(b-a);
c  double qt = exp(x2*log(base));
c  jac=jac*qt*log(base);
      
