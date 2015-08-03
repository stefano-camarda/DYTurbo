      subroutine breitw(x1,mminsq,mmaxsq,rmass,rwidth,msq,wt)       
      implicit none
      include 'constants.f'
c---- Given a number 0<x<1 generate a mass-squared msq and a weight wt 
c---- such that mminsq<msq<mmaxsq
c---- points are generated around resonance position rmass, but 
c---- breit-wigner should still be included in the matrix element
c     wt is the jacobian between integration in msq and integration in x1
      double precision x1,mminsq,mmaxsq,rmass,rwidth,msq,wt
      double precision almin,almax,al,tanal
      include 'zerowidth.f'


      if (zerowidth) then
          tanal=0d0
          almax=+pi/two
          almin=-pi/two
      else
          almin=datan((mminsq-rmass**2)/rmass/rwidth)
          almax=datan((mmaxsq-rmass**2)/rmass/rwidth)
          al=(almax-almin)*x1+almin
          tanal=dtan(al)
      endif

      msq=rmass**2+rmass*rwidth*tanal
c---- bw=(1d0+tanal**2)*rmass**2*rwidth**2
      wt=(almax-almin)*rmass*rwidth*(1d0+tanal**2)
      return
      end

