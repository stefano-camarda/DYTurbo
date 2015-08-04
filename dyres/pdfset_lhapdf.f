*****************
* LHAPDF version*
*****************
      subroutine pdfini
      implicit none
      include 'masses.f'
      include 'lhapdf.f'
      include 'PDFerrors.f'
      include 'pdlabel.f'
      double precision amz,alphasPDF
      logical validPDF
      character*30 oldPDFname
      integer i
      logical lhapdfs
      common/lhapdfs/lhapdfs

      common/couple/amz

      lhapdfs=.true.
      
c      if (newinput .eqv. .false.) then
c        open(unit=21,file='lhapdf.DAT',status='old',err=999)
c        call checkversion(21,'lhapdf.DAT')
c        read(21,*) PDFname
c        read(21,*) PDFmember            
c        close(21)
c      endif
      
      oldPDFname=PDFname
      validPDF=.true.
      i=0
   20 continue
      i=i+1    
      if ((oldPDFname(i:i) .eq. '.') .or.
     .    (oldPDFname(i:i) .eq. ' ') .or.
     .    (oldPDFname(i:i) .eq. '[')) then
        validPDF=.true.
        if (oldPDFname(i:i+6) .eq. '.LHgrid') then        
          PDFname=oldPDFname(1:i-1)//'.LHgrid'
        else
          PDFname=oldPDFname(1:i-1)//'.LHpdf'
        endif
      endif  
      if ((i .lt. 20) .and. (validPDF .eqv. .false.)) goto 20
      
      if (validPDF .eqv. .false.) then
        write(6,*) 'Problem with PDFname'
        write(6,*)
        stop
      endif
      

      write(6,*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
      write(6,*) 'C                                                  C'
      write(6,*) 'C            DYRES  now calling LHAPDF             C'
      write(6,*) 'C                                                  C'
      write(6,98) 'PDFname',PDFname(1:20)
      write(6,99) 'PDFmember',PDFmember
      write(6,*) 'C                                                  C'
      write(6,*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
      write(6,*)

      call InitPDFset('PDFsets/'//trim(PDFname))
      
      if (PDFmember .lt. 0) then
        PDFerrors=.true.
        call numberPDF(maxPDFsets)
        if (maxPDFsets .gt. 50) then
          write(6,*) 'ERROR: Max. number of error sets is 50!'
          stop
        endif
        write(6,*)
        write(6,*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'        
        write(6,*) 'C        Calculating errors using      C'
        write(6,*) 'C        ',maxPDFsets,' sets of error PDFs        C'
        write(6,*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
        call InitPDF(0)
        amz=alphasPDF(zmass)
        currentPDF=0
      else  
        call InitPDF(PDFmember)
        amz=alphasPDF(zmass)
      endif

c--- rename pdlabel to get sensible output name
      pdlabel=PDFname(1:7)

      return
 
   98 format(' C            ',a7,' ',a20,'          C')
   99 format(' C                ',a10,i3,'                     C')

  999 write(6,*) 'Error reading lhapdf.DAT'
      call flush(6)
      stop

      end
 

