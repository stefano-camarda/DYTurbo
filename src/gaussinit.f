      subroutine gaussinit
      include 'gauss.f'
      integer i,j

      do i = 1, 64
         do j = 1, 64
            xxx(i,j) = 0
            www(i,j) = 0
         enddo
      enddo

      do i = 1, 64
         if (i.le.4) xxx(4,i) = xxx4(i)
         if (i.le.8) xxx(8,i) = xxx8(i)
         if (i.le.10) xxx(10,i) = xxx10(i)
         if (i.le.20) xxx(20,i) = xxx20(i)
         if (i.le.24) xxx(24,i) = xxx24(i)
         if (i.le.40) xxx(40,i) = xxx40(i)
         if (i.le.50) xxx(50,i) = xxx50(i)
         if (i.le.64) xxx(64,i) = xxx64(i)

         if (i.le.4) www(4,i) = www4(i)
         if (i.le.8) www(8,i) = www8(i)
         if (i.le.10) www(10,i) = www10(i)
         if (i.le.20) www(20,i) = www20(i)
         if (i.le.24) www(24,i) = www24(i)
         if (i.le.40) www(40,i) = www40(i)
         if (i.le.50) www(50,i) = www50(i)
         if (i.le.64) www(64,i) = www64(i)
      enddo

      return
      end
