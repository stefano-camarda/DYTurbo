      subroutine zeromsq(msq,msqv)
      implicit none
      include 'constants.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf)

      msq=0d0
      msqv=0d0
      return
      
      do j=-nf,nf
      do k=-nf,nf
        msq(j,k)=0d0
        msqv(j,k)=0d0
      enddo      
      enddo      

      return
      end
