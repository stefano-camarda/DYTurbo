c********************************
c.....Dilogarithmic function Li2 
c********************************
     
c.....Function Li2(x) for arguments x < = 1.0                             
      function dyLI2(x)                                                     
      implicit none                                                           
      real*8 X,Y,T,S,A,PI3,PI6,ZERO,ONE,HALF,MALF,MONE,MTWO  
      real*8 C(0:18),H,ALFA,B0,B1,B2,LI2OLD                 
      real*8 dyLi2                                             
      integer i                                                       
      
      data ZERO /0.0d0/, ONE /1.0d0/                               
      data HALF /0.5d0/, MALF /-0.5d0/                             
      data MONE /-1.0d0/, MTWO /-2.0d0/                            
      data PI3 /3.289868133696453d0/, PI6 /1.644934066848226d0/               
      data C( 0) / 0.4299669356081370d0/
      data C( 1) / 0.4097598753307711d0/                              
      data C( 2) /-0.0185884366501460d0/                              
      data C( 3) / 0.0014575108406227d0/                              
      data C( 4) /-0.0001430418444234d0/                              
      data C( 5) / 0.0000158841554188d0/                              
      data C( 6) /-0.0000019078495939d0/                              
      data C( 7) / 0.0000002419518085d0/                              
      data C( 8) /-0.0000000319334127d0/                              
      data C( 9) / 0.0000000043454506d0/                              
      data C(10) /-0.0000000006057848d0/                              
      data C(11) / 0.0000000000861210d0/                              
      data C(12) /-0.0000000000124433d0/                              
      data C(13) / 0.0000000000018226d0/                              
      data C(14) /-0.0000000000002701d0/                              
      data C(15) / 0.0000000000000404d0/                              
      data C(16) /-0.0000000000000061d0/                              
      data C(17) / 0.0000000000000009d0/                              
      data C(18) /-0.0000000000000001d0/                              
      
      if(x .gt. 1.00000000001d0) then                                    
         write(6,*)'problems in dyLI2'
         write(6,*)'x=',x 
         stop                                               
      elseif(x .gt. 1.0d0) then                                          
         x = 1.d0                                                      
      endif                                                              
      if(X .eq. ONE) then                                                
         LI2OLD=PI6
         dyLI2=LI2OLD                                                       
         return                                                            
      else if(X .eq. MONE) then                                          
         LI2OLD=MALF*PI6
         dyLI2=LI2OLD                                                  
         return                                                            
      end if                                                             
      T=-X                                                               
      if(T .le. MTWO) then                                               
         Y=MONE/(ONE+T)                                                    
         S=ONE                                                             
         A=-PI3+HALF*(LOG(-T)**2-LOG(ONE+ONE/T)**2)                        
      else if(T .lt. MONE) then                                          
         Y=MONE-T                                                          
         S=MONE                                                            
         A=LOG(-T)                                                         
         A=-PI6+A*(A+LOG(ONE+ONE/T))                                       
      else if(T .le. MALF) then                                          
         Y=(MONE-T)/T                                                      
         S=ONE                                                             
         A=LOG(-T)                                                         
         A=-PI6+A*(MALF*A+LOG(ONE+T))                                      
      else if(T .lt. ZERO) then                                          
         Y=-T/(ONE+T)                                                      
         S=MONE                                                            
         A=HALF*LOG(ONE+T)**2                                              
      else if(T .le. ONE) then                                           
         Y=T                                                               
         S=ONE                                                             
         A=ZERO                                                            
      else                                                               
         Y=ONE/T                                                           
         S=MONE                                                            
         A=PI6+HALF*LOG(T)**2                                              
      end if                                                             
      
      H=Y+Y-ONE                                                          
      ALFA=H+H                                                           
      B1=ZERO                                                            
      B2=ZERO                                                            
      do I = 18,0,-1                                                    
         B0=C(I)+ALFA*B1-B2                                               
         B2=B1                                                            
         B1=B0                                                            
      enddo                                                              
      LI2OLD=-(S*(B0-H*B2)+A) 
      dyLI2=LI2OLD
      end     
