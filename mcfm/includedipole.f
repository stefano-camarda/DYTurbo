      logical function includedipole(nd,ptrans)
      implicit none
      include 'constants.f'
      include 'npart.f'
      include 'ptilde.f'
      double precision ptrans(mxpart,4),pjet(mxpart,4),rcut,pt34
      double precision qt2,xqtcut,Q2
      integer i,j,nd,nqcdjets,nqcdstart,notag,isub
      logical cuts,failedcuts,makecuts

      common/nqcdjets/nqcdjets,nqcdstart
      common/rcut/rcut
      common/makecuts/makecuts
      common/notag/notag
      common/qtcut/xqtcut
      
      logical binner
      external binner

      includedipole=.true.

CC    Compute qt2,Q2

      qt2=(ptrans(3,1)+ptrans(4,1))**2+(ptrans(3,2)+ptrans(4,2))**2      

      Q2=(ptrans(3,4)+ptrans(4,4))**2
     #    -(ptrans(3,1)+ptrans(4,1))**2      
     #    -(ptrans(3,2)+ptrans(4,2))**2      
     #    -(ptrans(3,3)+ptrans(4,3))**2 

c      print *,nd, q2, qt2

      if(dsqrt(qt2/Q2).lt.xqtcut) then
         includedipole=.false.
         return
      endif

      do j=1,4
         do i=1,npart+2 !i=3,4
            ptildejet(nd,i,j)=ptrans(i,j)
            pjet(i,j)=ptrans(i,j)
         enddo
      enddo

c--- check the lepton cuts
      failedcuts=(cuts(pjet,0))
      if (failedcuts) includedipole=.false.
      if (.not.binner(pjet(3,:),pjet(4,:))) includedipole=.false.
      
      return
      end
            
