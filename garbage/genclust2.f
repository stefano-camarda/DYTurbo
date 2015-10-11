      subroutine genclust2(q,R,qfinal,isub)
c--- this is a wrapper routine for the jet clustering algorithm
c--- either re-route to:
c---  genclust_kt.f     for kt clustering
c---  genclust_cone.f   for cone algorithm
      implicit none
      include 'constants.f'
      include 'clustering.f'
      double precision q(mxpart,4),qfinal(mxpart,4),R
      integer nqcdjets,nqcdstart,isub
      logical first
      character*4 part
      common/part/part
      common/nqcdjets/nqcdjets,nqcdstart
      data first/.true./
      save first

      if (algorithm .eq. 'ktal') then
        call genclust_kt(q,R,qfinal,isub)
      else
        call genclust_cone(q,R,qfinal,isub)
      endif

      return
      end
      
