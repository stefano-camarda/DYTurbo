      subroutine evtwriter(p,wgt,ii)
      implicit none
      include 'clustering.f'
      include 'constants.f'
      include 'cutoff.f'
      include 'jetlabel.f'
      include 'npart.f'
      include 'mxdim.f'
      include 'process.f'
      include 'removebr.f'
      include 'masses.f'


      character*2 plabel(mxpart)
      integer it,ch,ii,pr
      integer n,switch,nplotmax,i,j,myswitch 
      integer itmx1,ncall1,itmx2,ncall2
      common/plabel/plabel
      common/myswitch/myswitch
      common/it/it
      common/iterat/itmx1,ncall1,itmx2,ncall2
      common/pr/pr

      character tag*4

      double precision wgt,wt,p(mxpart,4)

      integer eventpart,nqcdjets,nqcdstart

      integer mysw 
      common/mysw/mysw


      if(ii.eq.0) ch=120
      if(ii.ge.1) ch=121

!      if(ii.eq.12) ch=122
!      if(ii.eq.13) ch=123

!           write(*,*) mysw



      if ((it-itmx2+pr).ge.1) then
 
      write(ch,"(A3,2x,I2,2x,I10)")  '<e>',it-itmx2+pr,mysw
      write(ch,"(1x,I2,1x,E23.15)")  ii,wgt*dfloat(itmx2)/dfloat(pr)
      do i=3,4  
      write(ch,"(1x,A2,2x,4E15.7)")  ! leptons 4-momenta  
     &  plabel(i),p(i,1),p(i,2),p(i,3),p(i,4)
      enddo
        write(ch,"(A4)")  '</e>' 

        else
      mysw=0
      endif

      return 
      end
      

      
