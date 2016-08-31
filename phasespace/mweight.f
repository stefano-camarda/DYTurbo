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

c  //split the mass in 3 pieces,
c  //and unweight the breit wigner only from -5 width to +5 width
c  //no sense, I have to split the mass before entering the integration
c  //and set breit wigner unweighting or not
c  double bwmn = opts.rmass -5*opts.rwidth;
c  double bwmx = opts.rmass +5*opts.rwidth;
c  double bwmnsq = pow(bwmn,2);
c  double bwmxsq = pow(bwmx,2);
c  double m,m2,wt;
c  double xl = 0.25;
c  double xu = 0.75;
c  if (x[0] < xl)
c    {
c      double x1=x[0]/(xl-0.);
c      m=phasespace::mmin+(bwmn-phasespace::mmin)*x[0];
c      jac=jac*(bwmn-phasespace::mmin)/(xl-0.);
c    }
c  else if (x[0] > xu)
c    {
c      double x1=(x[0]-xu)/(1.-xu);
c      m=bwmx+(phasespace::mmax-bwmx)*x1;
c      jac=jac*(phasespace::mmax-bwmx)/(1.-xu);
c    }
c  else
c    {
c      double xbw=(x[0]-xl)/(xu-xl);
c      breitw_(xbw,bwmnsq,bwmxsq,opts.rmass,opts.rwidth,m2,wt);
c      m=sqrt(m2);
c      wt=wt/(xu-xl);
c      jac=jac*wt;
c    }
      
