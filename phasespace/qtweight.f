      subroutine qtweight_lo(x,qtmin,qtmax,q2,qt,jac)
      implicit none
      double precision x,qtmin,qtmax,qt,jac
      double precision qtmin2,qtmax2,qt2,q2
      double precision lna,lnb
      
c     The behaviour of V+jet and CT ds^2/dqt^2
c     in the qt -> 0 limit at order n is:
c     Sum_k=1^2n Q^2/qt^2 * (ln(Q^2/qt^2))^(k-1)
c     See Eq. 4 of PRL 98, 222002 (2007)

c     The master formula for deriving a change of variable x -> t
c     for the integration of f(x) which makes the integrand flat is given here:
c     https://www.ippp.dur.ac.uk/~krauss/Lectures/QuarksLeptons/Basics/PS_3.html#TrivialUnweighting
c     The formula is:
c     x -> t = Int_x0^x f(x')dx' / Int_x0^x1 f(x')dx'
c     At O(alphas):
c     qt2 -> t = ln(q2/qt2)*(2+ln(q2/qt2)) - lna / (lnb-lna)
c     with lna = ln(q2/qtmin2)*(2+ln(q2/qtmin2))
c     and lnb =  ln(q2/qtmax2)*(2+ln(q2/qtmax2))
c     qt2 = q2 * exp(1-sqrt(1+lna+t(lnb-lna)))
c     dqt2 = q2*(lnb-lna)/(2*sqrt(1+lna+t(lnb-lna))) * exp(1-sqrt(1+lna+t(lnb-lna))) * dt
c     dqt2 = (lnb-lna)/(2*sqrt(1+lna+t(lnb-lna))) * qt2 * dt
      
c      qtmin2 = qtmin**2
c      qtmax2 = qtmax**2
c      lna = log(q2/qtmin2)*(2d0+log(q2/qtmin2))
c      lnb = log(q2/qtmax2)*(2d0+log(q2/qtmax2))
c      qt2 = q2 * exp(1d0-sqrt(1d0+lna+x*(lnb-lna)))
c      jac=-jac*(lnb-lna)/(2*sqrt(1d0+lna+x*(lnb-lna)))*qt2
c      qt=sqrt(qt2)
c      jac=jac/qt/2d0

c     Actually the leading behaviour for qt->0 is ~Q^2/qt^2
      qt = qtmax*qtmin/((qtmin-qtmax)*x+qtmax)
      jac=-jac*(qtmin-qtmax)*qtmin*qtmax/((qtmin-qtmax)*x+qtmax)**2

c     However none of the above seems to be more efficient than the change of variable qt2=tiny*exp(1d0/xx - 1d0)
      
      return
      end


      subroutine qtweight(x,qtmin,qtmax,qt,jac)
      implicit none
      double precision x,qtmin,qtmax,qt,jac
      double precision qtmin2,qtmax2,qt2
      double precision tiny
      parameter (tiny=1d-5) !be carefull, if qtmin^2 < tiny the phase space generation is screwed up
      double precision a,b,xx

      qtmin2 = qtmin**2
      qtmax2 = qtmax**2
      a = 1d0/(1+log(qtmax2/tiny))
      b = 1d0/(1+log(qtmin2/tiny))
      xx = a + (b-a) * x
      jac = jac*(b-a)
      qt2=tiny*exp(1d0/xx - 1d0)
      jac=jac*qt2/xx**2
      qt=sqrt(qt2)
      jac=jac/qt/2d0
      return
      end

      subroutine qt2weight(x,qtmin2,qtmax2,qt2,jac)
      implicit none
      double precision x,qtmin2,qtmax2,qt2,jac
      double precision tiny
      parameter (tiny=1d-5)
      double precision a,b,xx

      a = 1d0/(1+log(qtmax2/tiny))
      b = 1d0/(1+log(qtmin2/tiny))
      xx = a + (b-a) * x
      jac = jac*(b-a)
      qt2=tiny*exp(1d0/xx - 1d0)
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
      qt = tiny*exp(xx**(1d0/esp))
      jac=jac*qt*(xx**(1d0/esp-1d0))/esp
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
      
