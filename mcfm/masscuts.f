      subroutine masscuts(s,*)
      implicit none
      include 'constants.f'
      include 'limits.f'
      double precision s(mxpart,mxpart)
      
      if (  (s(3,4) .lt. wsqmin) 
     . .or. (s(3,4) .gt. wsqmax))
     .  return 1
      
      return
      end

