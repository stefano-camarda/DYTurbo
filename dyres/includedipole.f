      logical function includedipole(nd,ptrans)
      implicit none
      include 'constants.f'
      include 'clustering.f'
      include 'npart.f'
      include 'ptilde.f'
      include 'jetlabel.f'
      double precision ptrans(mxpart,4),pjet(mxpart,4),rcut,pt34
      double precision qt2,xqtcut,Q2
      integer i,j,nd,nqcdjets,nqcdstart,notag,isub
      logical cuts,failedcuts,makecuts,isolation,isol

      common/nqcdjets/nqcdjets,nqcdstart
      common/rcut/rcut
      common/makecuts/makecuts
      common/notag/notag
      common/qtcut/xqtcut
      common/isol/isol
      
      logical binner
      external binner

      includedipole=.true.

      if (nd .gt. 0) then
        isub=1
      else
        isub=0
      endif

CC    Compute qt2,Q2

       qt2=(ptrans(3,1)+ptrans(4,1))**2+(ptrans(3,2)+ptrans(4,2))**2      

       Q2=(ptrans(3,4)+ptrans(4,4))**2
     #    -(ptrans(3,1)+ptrans(4,1))**2      
     #    -(ptrans(3,2)+ptrans(4,2))**2      
     #    -(ptrans(3,3)+ptrans(4,3))**2 


      if(dsqrt(qt2/Q2).lt.xqtcut) then
       includedipole=.false.
      else

       call genclust2(ptrans,rcut,pjet,isub)
       do j=1,4
        do i=1,npart+2
        ptildejet(nd,i,j)=pjet(i,j)
        enddo
       enddo


CC    Insert here isolation cut

      if(isol) then
       if(isolation(ptrans).eqv..false.) includedipole=.false.
      endif 
 
c--- check the lepton cuts

        if (makecuts) then
           failedcuts=(cuts(pjet,jets).or.
     .          (.not.binner(pjet(3,:),pjet(4,:))))
          if (failedcuts) includedipole=.false.
        endif
 
      endif
      
      return
      end
            
