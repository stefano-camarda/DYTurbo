      subroutine iniflavreduce
      implicit none
      integer nproc
      common/nproc/nproc
      include 'flred.f'

c     reduce flavour double loop to single loop
      if (nproc.eq.3) then
         nfme = 10

         jjm(1)=1
         jjm(2)=-1
         jjm(3)=2
         jjm(4)=-2
         jjm(5)=3
         jjm(6)=-3
         jjm(7)=4
         jjm(8)=-4
         jjm(9)=5
         jjm(10)=-5

         kkm(1)=-1
         kkm(2)=1
         kkm(3)=-2
         kkm(4)=2
         kkm(5)=-3
         kkm(6)=3
         kkm(7)=-4
         kkm(8)=4
         kkm(9)=-5
         kkm(10)=5

         jjs(1)=2	
         jjs(2)=-2	
         jjs(3)=1	
         jjs(4)=-1	
         jjs(5)=3	
         jjs(6)=-3	
         jjs(7)=4	
         jjs(8)=-4	
         jjs(9)=5	
         jjs(10)=-5

         kks(1)=-2
         kks(2)=2 
         kks(3)=-1
         kks(4)=1 
         kks(5)=-3
         kks(6)=3 
         kks(7)=-4
         kks(8)=4 
         kks(9)=-5
         kks(10)=5 

      elseif(nproc.eq.1) then
         nfme = 12

         jjm(1)=2	
         jjm(2)=-1	
         jjm(3)=2	
         jjm(4)=-3	
         jjm(5)=2	
         jjm(6)=-5	
         jjm(7)=4	
         jjm(8)=-3	
         jjm(9)=4	
         jjm(10)=-1	
         jjm(11)=4	
         jjm(12)=-5	

         kkm(1)=-1
         kkm(2)=2 
         kkm(3)=-3
         kkm(4)=2 
         kkm(5)=-5
         kkm(6)=2 
         kkm(7)=-3
         kkm(8)=4 
         kkm(9)=-1
         kkm(10)=4 
         kkm(11)=-5
         kkm(12)=4 

         jjs(1)=1	
         jjs(2)=-2	
         jjs(3)=1	
         jjs(4)=-3	
         jjs(5)=1	
         jjs(6)=-5	
         jjs(7)=4	
         jjs(8)=-3	
         jjs(9)=4	
         jjs(10)=-2	
         jjs(11)=4	          
         jjs(12)=-5	

         kks(1)=-2
         kks(2)=1 
         kks(3)=-3
         kks(4)=1 
         kks(5)=-5
         kks(6)=1 
         kks(7)=-3
         kks(8)=4 
         kks(9)=-2
         kks(10)=4 
         kks(11)=-5
         kks(12)=4 
         
      elseif(nproc.eq.2) then
         nfme = 12
         
         jjm(1)=1	
         jjm(2)=-2	
         jjm(3)=3	
         jjm(4)=-2	
         jjm(5)=5	
         jjm(6)=-2	
         jjm(7)=3	
         jjm(8)=-4	
         jjm(9)=1	
         jjm(10)=-4	
         jjm(11)=5	
         jjm(12)=-4

         kkm(1)=-2	
         kkm(2)=1	
         kkm(3)=-2	
         kkm(4)=3	
         kkm(5)=-2	
         kkm(6)=5	
         kkm(7)=-4	
         kkm(8)=3	
         kkm(9)=-4	
         kkm(10)=1	
         kkm(11)=-4	
         kkm(12)=5	

         jjs(1)=2	
         jjs(2)=-1	
         jjs(3)=3	
         jjs(4)=-1	
         jjs(5)=5	
         jjs(6)=-1	
         jjs(7)=3	
         jjs(8)=-4	
         jjs(9)=2	
         jjs(10)=-4	
         jjs(11)=5	
         jjs(12)=-4	


         kks(1)=-1        
         kks(2)=2         
         kks(3)=-1        
         kks(4)=3         
         kks(5)=-1        
         kks(6)=5         
         kks(7)=-4        
         kks(8)=3         
         kks(9)=-4        
         kks(10)=2         
         kks(11)=-4        
         kks(12)=5         
      endif

      return
      end
