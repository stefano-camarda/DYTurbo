CC NEW: cuts only on QCD partons

      subroutine smalls(s,npart,*)
c    cut if radiated parton too close
      implicit none
      include 'constants.f'
      include 'cutoff.f'
      integer npart
      double precision s(mxpart,mxpart)

CC    Case H->34
   
      if (npart .eq. 3) then
      if ( 
     .      (-s(1,5) .lt. cutoff)
     . .or. (-s(2,5) .lt. cutoff)
     . ) return 1

      elseif (npart .eq. 4) then
      if ( 
     .      (-s(1,5) .lt. cutoff)
     . .or. (-s(2,5) .lt. cutoff)
     . .or. (-s(1,6) .lt. cutoff)
     . .or. (-s(2,6) .lt. cutoff)
     . .or. (+s(5,6) .lt. cutoff)
     . ) return 1
     
CC    Case H->3456

      elseif (npart .eq. 5) then
      if ( 
     .      (-s(1,7) .lt. cutoff)
     . .or. (-s(2,7) .lt. cutoff)
     . ) return 1

      elseif (npart .eq. 6) then
      if ( 
     .      (-s(1,7) .lt. cutoff)
     . .or. (-s(2,7) .lt. cutoff)
     . .or. (-s(1,8) .lt. cutoff)
     . .or. (-s(2,8) .lt. cutoff)
     . .or. (+s(7,8) .lt. cutoff)
     . ) return 1
      else
      write(*,*)'ERROR in SMALLS'
      stop
      endif      

      return
      end
