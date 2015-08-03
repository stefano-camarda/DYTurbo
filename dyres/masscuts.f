      subroutine masscuts(s,*)
      implicit none
      include 'constants.f'
      include 'npart.f'
      include 'limits.f'
      logical first
      double precision s(mxpart,mxpart)
      integer nqcdjets,nqcdstart
      common/nqcdjets/nqcdjets,nqcdstart
      data first/.true./
      save first

      
      if (  (s(3,4) .lt. wsqmin) 
     . .or. (s(3,4) .gt. wsqmax))
     .  return 1
      
c      if ((npart .gt. 3) .and. (nqcdjets .lt. 2)) then
c        if (  (s(5,6) .lt. bbsqmin) 
c     .   .or. (s(5,6) .gt. bbsqmax))
c     .    return 1
c      endif
     
   98 format(' *      ',f8.2,'  <   ',a12,'  < ',f8.2,'      *')
   99 format(' *          ',f8.2,'  <   ',a3,'  < ',f8.2,'           *')
     
      return
      end

