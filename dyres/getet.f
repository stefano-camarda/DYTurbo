      double precision function getet(E,px,py,pz)
c--- given (E,px,py,pz) for a four-vector, calculates the corresponding
c--- Et or Pt, depending on the parameter that is set in this routine
c--- Note that, at present, this routine is only called by
c--- genclust_cone.f
      implicit none
      double precision E,px,py,pz,etsq
      logical useEt
      parameter (useEt=.true.)
      
      if (useEt) then
c--- this is the formula for Et
        etsq=px**2+py**2
        getet=dsqrt(etsq)*E/dsqrt(etsq+pz**2)
      else
c--- this is the formula for pt
        getet=dsqrt(px**2+py**2)
      endif
      
      return
      end
      
