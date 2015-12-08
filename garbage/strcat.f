      subroutine strcat(str1,str2,str)
c concatenates str1 and str2 into str. Ignores trailing blanks of str1,str2
      character *(*) str1,str2,str
      l1=istrl(str1)
      l2=istrl(str2)
      l =len(str)
      if(l.lt.l1+l2) then
          write(*,*) 'error: l1+l2>l in strcat'
          write(*,*) 'l1=',l1,' str1=',str1
          write(*,*) 'l2=',l2,' str2=',str2
          write(*,*) 'l=',l
          stop
      endif
      if(l1.ne.0) str(1:l1)=str1(1:l1)
      if(l2.ne.0) str(l1+1:l1+l2)=str2(1:l2)
      if(l1+l2+1.le.l) str(l1+l2+1:l)= ' '
      end

      function istrl(string)
c returns the position of the last non-blank character in string
      character * (*) string
      i = len(string)
      dowhile(i.gt.0.and.string(i:i).eq.' ')
         i=i-1
      enddo
      istrl = i
      end


      subroutine cstring(string)
c     center a string in its length (MG July 2007)
      implicit none
      character * (*) string
      integer i,j,l,l1,istrl
      l=len(string)
      l1=istrl(string)
      i=(l-l1)/2


      do j=0,l1-1
      string(l1+i-j:l1+i-j)=string(l1-j:l1-j)
      enddo

      do j=1,i
      string(j:j)=' '
      enddo

      return
      end

