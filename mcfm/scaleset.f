      subroutine scaleset(q2)
      implicit none
      include 'constants.f'
      include 'scale.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'nwz.f'
      include 'facscale.f'
      include 'nlooprun.f'

      double precision q2,scalemax,amz,dyalphas
      common/couple/amz
      
      scale=dsqrt(q2)
      facscale=dsqrt(q2)

      scalemax=3000d0

c--- catch absurdly large scales      
      if  (scale.gt.scalemax) then
       scale=scalemax
       facscale=scalemax
      endif
       

c--- run alpha_s
      as=dyalphas(scale,amz,nlooprun)
        
      ason2pi=as/twopi
      ason4pi=as/fourpi
      gsq=fourpi*as
      musq=scale**2
      
      return
            
      end
