      SUBROUTINE partons(xmu,X,FX,NF,ISET,IH)
      implicit none
      integer Iprtn,ih,Irt,NF,ISET
      double precision fx(-NF:NF),x,xmu,sqxmu,fPDF(-6:6),temp
c---  ih1=+1 proton 
c---  ih1=-1 pbar 

C---set to zero if x out of range
      if (x .gt. 1d0) then
          do Iprtn=-5,5
             fx(Iprtn)=0d0
          enddo
          return
      endif
      sqxmu=sqrt(xmu)
      call evolvePDF(x,sqxmu,fPDF)

      if (ih.eq.1) then
        do Iprtn=-5,5
          fx(Iprtn)=fPDF(Iprtn)/x
        enddo
      elseif(ih.eq.-1) then
        do Iprtn=-5,5
          fx(+Iprtn)=fPDF(-Iprtn)/x
       enddo
      endif
! Change u and d quarks
      temp=fx(1)
      fx(1)=fx(2)
      fx(2)=temp
      temp=fx(-1)
      fx(-1)=fx(-2)
      fx(-2)=temp

      fx(-5)=0d0
      fx(-4)=0d0
      fx(-3)=0d0
      fx(-2)=0d0
      fx(-1)=0d0
      fx(0)=0d0
      fx(1)=0d0
c      fx(2)=0d0
      fx(3)=0d0
      fx(4)=0d0
      fx(5)=0d0
      
      return
      end

      subroutine flavour
      implicit real *8 (a-h,o-z)
      real *8 fh1(-5:5),fh2(-5:5),sumckm(5)
      integer nf,ih1,ih2,ic,isetproton,nloop,ord,prodflag
      common/scales2/xmur,xmuf,xmur2,xmuf2
      common/isetproton/isetproton
      common/pdf/ih1,ih2
      common/nf/nf
      common/fractions/x1,x2
      common/quarks/eq(5),alq(5),arq(5),ckm(6,6),delta(5,5),tau3(5,5)
      common/couplings/xw,cw,sw,alpha0
      common/prodflag/prodflag
      common/luminosities/xlumgg,xlumqg,xlumgq,xlumqgtr,xlumgqtr,
     /     xlumqqb,xlumqqbtr,xlumqqbdbb,xlumqqbdbc,xlumqqbdcc,
     /     xlumqqbLL,xlumqqbLR,xlumqq,xlumqqeaa,xlumqqead,
     /     xlumqqLL,xlumqqLR
      common/flagfit/flagfit
      real *8 siggamma,sigint,sigz,sigw
      common/sigs/siggamma,sigint,sigz,sigw
            call partons(xmuf2,x1,fh1,nf,isetproton,ih1)
            call partons(xmuf2,x2,fh2,nf,isetproton,ih2)
c      print *,x1,sqrt(xmuf2),fh1(-1),fh1(-2),fh1(-3),fh1(-4),fh1(-5)
c            call fdist(ih1,x1,sqrt(xmuf2),fh1)
c            call fdist(ih2,x2,sqrt(xmuf2),fh2)
c      print *,x1,sqrt(xmuf2),fh1(-2),fh1(-1),fh1(-3),fh1(-4),fh1(-5)

      sumckm(1)=ckm(1,2)**2+ckm(4,2)**2
      sumckm(2)=ckm(1,2)**2+ckm(1,3)**2+ckm(1,5)**2
      sumckm(3)=ckm(4,2)**2+ckm(4,3)**2+ckm(4,5)**2
      sumckm(4)=ckm(1,3)**2+ckm(4,3)**2
      sumckm(5)=ckm(6,2)**2+ckm(6,3)**2+ckm(6,5)**2

c.....gg luminosity
      xlumgg=0d0
      pdf=fh1(0)*fh2(0)
      do i=1,nf
         do j=1,nf
            flav=0   
            if (prodflag.eq.1) then
               flav=2*eq(j)**2*delta(i,j)*siggamma
            elseif (prodflag.eq.2) then
               flav=1/(2*xw)*ckm(i,j)**2*sigw
            elseif (prodflag.eq.21.or.prodflag.eq.22) then
               flav=1/(2*xw)*ckm(i,j)**2/2d0*sigw
            elseif (prodflag.eq.3) then
               flav=((1/(2*sw*cw)*tau3(i,j)-delta(i,j)*eq(j)*sw/cw)**2+
     /              (-delta(i,j)*eq(j)*sw/cw)**2)*sigz
            elseif (prodflag.eq.4) then
               flav=(((1/(2*sw*cw)*tau3(i,j)-delta(i,j)*eq(j)*sw/cw)+
     /              (-delta(i,j)*eq(j)*sw/cw))*2d0*eq(j))*sigint
            elseif (prodflag.eq.5) then
               flav=2*eq(j)**2*delta(i,j)*siggamma+
     / ((1/(2*sw*cw)*tau3(i,j)-delta(i,j)*eq(j)*sw/cw)**2+
     /              (-delta(i,j)*eq(j)*sw/cw)**2)*sigz +
     / (((1/(2*sw*cw)*tau3(i,j)-delta(i,j)*eq(j)*sw/cw)+
     /              (-delta(i,j)*eq(j)*sw/cw))*2d0*eq(j))*sigint
            endif
            xlumgg=xlumgg+pdf*flav
         enddo
      enddo
c.....qg luminosity
      xlumqg=0d0
      do i=1,nf
         pdf=0
         pdf=(fh1(i)+fh1(-i))*fh2(0)
         do j=1,nf
            flav=0
            if (prodflag.eq.1) then
               flav=2*eq(j)**2*delta(i,j)*siggamma
            elseif (prodflag.eq.2) then
               flav=1/(2*xw)*delta(i,j)*sumckm(i)*sigw
            elseif (prodflag.eq.21) then
               pdf=(fh1(i)*sumckm(i)*
     /              (delta(i,1)+delta(i,4))+
     /              fh1(-j)*sumckm(j)*
     /              (delta(j,2)+delta(j,3)+delta(j,5)))*fh2(0) 
               flav=1/(2*xw)*delta(i,j)*sigw
            elseif (prodflag.eq.22) then
               pdf=(fh1(-i)*sumckm(i)*
     /              (delta(i,1)+delta(i,4))+
     /              fh1(j)*sumckm(j)*
     /              (delta(j,2)+delta(j,3)+delta(j,5)))*fh2(0) 
               flav=1/(2*xw)*delta(i,j)*sigw
            elseif (prodflag.eq.3) then
               flav=((1/(2*sw*cw)*tau3(i,j)-delta(i,j)*eq(j)*sw/cw)**2+
     /            (-delta(i,j)*eq(j)*sw/cw)**2)*sigz
            elseif (prodflag.eq.4) then
               flav=(((1/(2*sw*cw)*tau3(i,j)-delta(i,j)*eq(j)*sw/cw)+
     /              (-delta(i,j)*eq(j)*sw/cw))*2d0*eq(j))*sigint
            elseif (prodflag.eq.5) then
               flav=2*eq(j)**2*delta(i,j)*siggamma+
     /((1/(2*sw*cw)*tau3(i,j)-delta(i,j)*eq(j)*sw/cw)**2+
     /            (-delta(i,j)*eq(j)*sw/cw)**2)*sigz +
     / (((1/(2*sw*cw)*tau3(i,j)-delta(i,j)*eq(j)*sw/cw)+
     /              (-delta(i,j)*eq(j)*sw/cw))*2d0*eq(j))*sigint
            endif
            xlumqg=xlumqg+pdf*flav
         enddo
      enddo

c.....gq luminosity
      xlumgq=0d0
      do i=1,nf
         pdf=0
         pdf=(fh2(i)+fh2(-i))*fh1(0)
         do j=1,nf
            flav=0
            if (prodflag.eq.1) then
               flav=2*eq(j)**2*delta(i,j)*siggamma
            elseif (prodflag.eq.2) then
               flav=1/(2*xw)*delta(i,j)*sumckm(i)*sigw
            elseif (prodflag.eq.21) then
               pdf=(fh2(i)*sumckm(i)*
     /              (delta(i,1)+delta(i,4))+
     /              fh2(-j)*sumckm(j)*
     /              (delta(j,2)+delta(j,3)+delta(j,5)))*fh1(0) 
               flav=1/(2*xw)*delta(i,j)*sigw
            elseif (prodflag.eq.22) then
               pdf=(fh2(-i)*sumckm(i)*
     /              (delta(i,1)+delta(i,4))+
     /              fh2(j)*sumckm(j)*
     /              (delta(j,2)+delta(j,3)+delta(j,5)))*fh1(0) 
               flav=1/(2*xw)*delta(i,j)*sigw
            elseif (prodflag.eq.3) then
               flav=((1/(2*sw*cw)*tau3(i,j)-delta(i,j)*eq(j)*sw/cw)**2+
     /              (-delta(i,j)*eq(j)*sw/cw)**2)*sigz
            elseif (prodflag.eq.4) then
               flav=(((1/(2*sw*cw)*tau3(i,j)-delta(i,j)*eq(j)*sw/cw)+
     /              (-delta(i,j)*eq(j)*sw/cw))*2d0*eq(j))*sigint
            elseif (prodflag.eq.5) then
               flav=2*eq(j)**2*delta(i,j)*siggamma+
     /((1/(2*sw*cw)*tau3(i,j)-delta(i,j)*eq(j)*sw/cw)**2+
     /              (-delta(i,j)*eq(j)*sw/cw)**2)*sigz+
     /(((1/(2*sw*cw)*tau3(i,j)-delta(i,j)*eq(j)*sw/cw)+
     /              (-delta(i,j)*eq(j)*sw/cw))*2d0*eq(j))*sigint
            endif
            xlumgq=xlumgq+pdf*flav
         enddo
      enddo      

c.....triangular loop qg luminosity (only for Z production)
      xlumqgtr=0d0
      do i=1,nf 
         pdf=0
         pdf=(fh1(i)+fh1(-i))*fh2(0)*tau3(i,i)
         flav=0d0
         if (prodflag.eq.1) then
            flav=0d0
         elseif (prodflag.eq.2) then
            flav=0d0
         elseif (prodflag.eq.21) then
            flav=0d0
         elseif (prodflag.eq.22) then
            flav=0d0
         elseif (prodflag.eq.4) then
            flav=0d0
         elseif (prodflag.eq.3) then
            flav=-(1/(2*sw*cw))*(1/(2*sw*cw))*sigz
         elseif (prodflag.eq.5) then
            flav=-(1/(2*sw*cw))*(1/(2*sw*cw))*sigz
         endif
         xlumqgtr=xlumqgtr+pdf*flav
      enddo

c.....triangular loop gq luminosity (only for Z production)
      xlumgqtr=0d0
      do i=1,nf
         pdf=0
         pdf=(fh2(i)+fh2(-i))*fh1(0)*tau3(i,i)
         flav=0d0
         if (prodflag.eq.1) then
            flav=0d0
         elseif (prodflag.eq.2) then
            flav=0d0
         elseif (prodflag.eq.21) then
            flav=0d0
         elseif (prodflag.eq.22) then
            flav=0d0
         elseif (prodflag.eq.4) then
            flav=0d0
         elseif (prodflag.eq.3) then
            flav=-(1/(2*sw*cw))*(1/(2*sw*cw))*sigz
         elseif (prodflag.eq.5) then
            flav=-(1/(2*sw*cw))*(1/(2*sw*cw))*sigz
         endif
         xlumgqtr=xlumgqtr+pdf*flav
      enddo      
      
c.....qqb luminosity
      xlumqqb=0d0
      do i=1,nf
         do j=1,nf
            pdf=0
            pdf=fh1(i)*fh2(-j)+fh1(-i)*fh2(j)
            flav=0
            if (prodflag.eq.1) then
               flav=2*eq(j)**2*delta(i,j)*siggamma
            elseif (prodflag.eq.2) then
               flav=1/(2*xw)*ckm(i,j)**2*sigw
            elseif (prodflag.eq.21) then
               pdf=0
               pdf=(fh1(i)*fh2(-j)+fh2(i)*fh1(-j))*
     /              (delta(i,1)+delta(i,4))*
     /              (delta(j,2)+delta(j,3)+delta(j,5)) 
               flav=1/(2*xw)*ckm(i,j)**2*sigw
            elseif (prodflag.eq.22) then
               pdf=0
               pdf=(fh1(-i)*fh2(j)+fh2(-i)*fh1(j))*
     /              (delta(i,1)+delta(i,4))*
     /              (delta(j,2)+delta(j,3)+delta(j,5)) 
               flav=1/(2*xw)*ckm(i,j)**2*sigw
            elseif (prodflag.eq.3) then
               flav=((1/(2*sw*cw)*tau3(i,j)-delta(i,j)*eq(j)*sw/cw)**2+
     /              (-delta(i,j)*eq(j)*sw/cw)**2)*sigz
            elseif (prodflag.eq.4) then
               flav=(((1/(2*sw*cw)*tau3(i,j)-delta(i,j)*eq(j)*sw/cw)+
     /              (-delta(i,j)*eq(j)*sw/cw))*2d0*eq(j))*sigint
            elseif (prodflag.eq.5) then
               flav=2*eq(j)**2*delta(i,j)*siggamma+
     /((1/(2*sw*cw)*tau3(i,j)-delta(i,j)*eq(j)*sw/cw)**2+
     /              (-delta(i,j)*eq(j)*sw/cw)**2)*sigz+
     /(((1/(2*sw*cw)*tau3(i,j)-delta(i,j)*eq(j)*sw/cw)+
     /              (-delta(i,j)*eq(j)*sw/cw))*2d0*eq(j))*sigint
            endif
            xlumqqb=xlumqqb+pdf*flav
         enddo
      enddo      

c.....triangular loop and Dab qqb luminosity (only for Z production)
      xlumqqbtr=0d0
      do i=1,nf
         pdf=0
         pdf=(fh1(i)*fh2(-i)+fh1(-i)*fh2(i))*tau3(i,i)
         flav=0
         if (prodflag.eq.1) then
            flav=0d0
         elseif (prodflag.eq.2) then
            flav=0d0
         elseif (prodflag.eq.21) then
            flav=0d0
         elseif (prodflag.eq.22) then
            flav=0d0
         elseif (prodflag.eq.4) then
            flav=0d0
         elseif (prodflag.eq.3) then
            flav=-(1/(2*sw*cw))*(1/(2*sw*cw))*sigz
         elseif (prodflag.eq.5) then
            flav=-(1/(2*sw*cw))*(1/(2*sw*cw))*sigz
         endif
         xlumqqbtr=xlumqqbtr+pdf*flav
      enddo      
      
c.....Dbb qqb luminosity
      xlumqqbdbb=0d0
      do m=1,nf
         pdf=0
         pdf=fh1(m)*fh2(-m)+fh1(-m)*fh2(m)
         do i=1,nf
           do j=1,nf
              flav=0
              if (prodflag.eq.1) then
                 flav=2*eq(j)**2*delta(i,j)*siggamma
              elseif (prodflag.eq.2) then
                 flav=1/(2*xw)*ckm(i,j)**2*sigw
              elseif (prodflag.eq.21.or.prodflag.eq.22) then
                 flav=1/(2*xw)*ckm(i,j)**2/2d0*sigw
              elseif (prodflag.eq.3) then
                flav=((1/(2*sw*cw)*tau3(i,j)-delta(i,j)*eq(j)*sw/cw)**2+
     /                (-delta(i,j)*eq(j)*sw/cw)**2)*sigz
              elseif (prodflag.eq.5) then
                 flav=2*eq(j)**2*delta(i,j)*siggamma+
     /((1/(2*sw*cw)*tau3(i,j)-delta(i,j)*eq(j)*sw/cw)**2+
     /                (-delta(i,j)*eq(j)*sw/cw)**2)*sigz
              endif
              xlumqqbdbb=xlumqqbdbb+pdf*flav
           enddo
        enddo
      enddo

c.....Dbc and Dbd qqb luminosity
      xlumqqbdbc=0
      do i=1,nf
         do j=1,nf
            pdf=0
            pdf=fh1(i)*fh2(-j)+fh1(-i)*fh2(j)
            flav=0
            if (prodflag.eq.1) then
               flav=2*eq(j)**2*delta(i,j)*siggamma
            elseif (prodflag.eq.2) then
               flav=1/(2*xw)*delta(i,j)*sumckm(i)*sigw
            elseif (prodflag.eq.21) then
               pdf=(fh1(-j)*fh2(j)*sumckm(j)*
     /              (delta(j,2)+delta(j,3)+delta(j,5))+ 
     /              fh2(-i)*fh1(i)*sumckm(i)*
     /              (delta(i,1)+delta(i,4)))
               flav=1/(2*xw)*delta(i,j)*sigw
            elseif (prodflag.eq.22) then
               pdf=(fh1(j)*fh2(-j)*sumckm(j)*
     /              (delta(j,2)+delta(j,3)+delta(j,5))+ 
     /              fh2(i)*fh1(-i)*sumckm(i)*
     /              (delta(i,1)+delta(i,4)))
               flav=1/(2*xw)*delta(i,j)*sigw
            elseif (prodflag.eq.3) then
              flav=((1/(2*sw*cw)*tau3(i,j)-delta(i,j)*eq(j)*sw/cw)**2+
     /              (-delta(i,j)*eq(j)*sw/cw)**2)*sigz
            elseif (prodflag.eq.5) then
               flav=2*eq(j)**2*delta(i,j)*siggamma+
     /((1/(2*sw*cw)*tau3(i,j)-delta(i,j)*eq(j)*sw/cw)**2+
     /              (-delta(i,j)*eq(j)*sw/cw)**2)*sigz
            endif
            xlumqqbdbc=xlumqqbdbc+pdf*flav
         enddo
      enddo      

c.....Dcc and Ddd qqb luminosity
      xlumqqbdcc=0d0
      do i=1,nf
         do j=1,nf
            pdf=0
            pdf=(fh1(i)+fh1(-i))*(fh2(j)+fh2(-j))
            flav=0
            if (prodflag.eq.1) then
               flav=2*eq(i)**2*delta(i,i)*siggamma
            elseif (prodflag.eq.2) then
               flav=1/(2*xw)*sumckm(i)*sigw
            elseif (prodflag.eq.21) then
               pdf=(fh1(i)*sumckm(i)*
     /              (delta(i,1)+delta(i,4))+
     /              fh1(-i)*sumckm(i)*
     /              (delta(i,2)+delta(i,3)+delta(i,5)))* 
     /              (fh2(j)+fh2(-j)) 
               flav=1/(2*xw)*sigw
            elseif (prodflag.eq.22) then
               pdf=(fh1(-i)*sumckm(i)*
     /              (delta(i,1)+delta(i,4))+
     /              fh1(i)*sumckm(i)*
     /              (delta(i,2)+delta(i,3)+delta(i,5)))* 
     /              (fh2(j)+fh2(-j)) 
               flav=1/(2*xw)*sigw
            elseif (prodflag.eq.3) then
               flav=((1/(2*sw*cw)*tau3(i,i)-delta(i,i)*eq(i)*sw/cw)**2+
     /              (-delta(i,i)*eq(i)*sw/cw)**2)*sigz
            elseif (prodflag.eq.5) then
               flav=2*eq(i)**2*delta(i,i)*siggamma+
     /((1/(2*sw*cw)*tau3(i,i)-delta(i,i)*eq(i)*sw/cw)**2+
     /              (-delta(i,i)*eq(i)*sw/cw)**2)*sigz
            endif
            xlumqqbdcc=xlumqqbdcc+pdf*flav
         enddo
      enddo      

c.....LL qqb luminosity
      xlumqqbLL=0
      do i=1,nf        
         do j=1,nf
            flav=0
            if (prodflag.eq.1) then
               flav=2*(eq(i)*(fh1(i)-fh1(-i)))*(eq(j)*(fh2(-j)-fh2(j)))
     /*siggamma
            elseif (prodflag.eq.2) then
               flav=0d0
            elseif (prodflag.eq.3) then
               flav=(((1/(2*sw*cw)*tau3(i,i)-eq(i)*sw/cw)*fh1(i)+
     /              (eq(i)*sw/cw)*fh1(-i))*
     /              ((1/(2*sw*cw)*tau3(j,j)-eq(j)*sw/cw)*fh2(-j)+
     /              (eq(j)*sw/cw)*fh2(j))+
     /              ((1/(2*sw*cw)*tau3(i,i)-eq(i)*sw/cw)*fh1(-i)+
     /              (eq(i)*sw/cw)*fh1(i))*
     /              ((1/(2*sw*cw)*tau3(j,j)-eq(j)*sw/cw)*fh2(j)+
     /              (eq(j)*sw/cw)*fh2(-j)))*sigz
            elseif (prodflag.eq.5) then
               flav=2*(eq(i)*(fh1(i)-fh1(-i)))*(eq(j)*(fh2(-j)-fh2(j)))
     /*siggamma+
     /              (((1/(2*sw*cw)*tau3(i,i)-eq(i)*sw/cw)*fh1(i)+
     /              (eq(i)*sw/cw)*fh1(-i))*
     /              ((1/(2*sw*cw)*tau3(j,j)-eq(j)*sw/cw)*fh2(-j)+
     /              (eq(j)*sw/cw)*fh2(j))+
     /              ((1/(2*sw*cw)*tau3(i,i)-eq(i)*sw/cw)*fh1(-i)+
     /              (eq(i)*sw/cw)*fh1(i))*
     /              ((1/(2*sw*cw)*tau3(j,j)-eq(j)*sw/cw)*fh2(j)+
     /              (eq(j)*sw/cw)*fh2(-j)))*sigz
            endif
            xlumqqbLL=xlumqqbLL+flav
         enddo
      enddo      

c.....LR qqb luminosity
      xlumqqbLR=0
      do i=1,nf        
         do j=1,nf
            flav=0
            if (prodflag.eq.1) then
               flav=2*(eq(i)*(fh1(i)-fh1(-i)))*(eq(j)*(fh2(-j)-fh2(j)))
     /*siggamma
            elseif (prodflag.eq.2) then
               flav=0d0
            elseif (prodflag.eq.3) then
               flav=(((1/(2*sw*cw)*tau3(i,i)-eq(i)*sw/cw)*fh1(i)+
     /              (eq(i)*sw/cw)*fh1(-i))*
     /              (-(1/(2*sw*cw)*tau3(j,j)-eq(j)*sw/cw)*fh2(j)+
     /              (-eq(j)*sw/cw)*fh2(-j))+
     /              ((1/(2*sw*cw)*tau3(i,i)-eq(i)*sw/cw)*fh1(-i)+
     /              (eq(i)*sw/cw)*fh1(i))*
     /              (-(1/(2*sw*cw)*tau3(j,j)-eq(j)*sw/cw)*fh2(-j)+
     /              (-eq(j)*sw/cw)*fh2(j)))*sigz
            elseif (prodflag.eq.5) then
               flav=2*(eq(i)*(fh1(i)-fh1(-i)))*(eq(j)*(fh2(-j)-fh2(j)))
     /*siggamma+
     /              (((1/(2*sw*cw)*tau3(i,i)-eq(i)*sw/cw)*fh1(i)+
     /              (eq(i)*sw/cw)*fh1(-i))*
     /              (-(1/(2*sw*cw)*tau3(j,j)-eq(j)*sw/cw)*fh2(j)+
     /              (-eq(j)*sw/cw)*fh2(-j))+
     /              ((1/(2*sw*cw)*tau3(i,i)-eq(i)*sw/cw)*fh1(-i)+
     /              (eq(i)*sw/cw)*fh1(i))*
     /              (-(1/(2*sw*cw)*tau3(j,j)-eq(j)*sw/cw)*fh2(-j)+
     /              (-eq(j)*sw/cw)*fh2(j)))*sigz
            endif
            xlumqqbLR=xlumqqbLR+flav
         enddo
      enddo      

c.....qq luminosity
      xlumqq=0
      do i=1,nf
         do j=1,nf
            pdf=0
            pdf=fh1(i)*fh2(j)+fh1(-i)*fh2(-j)   
            flav=0
            if (prodflag.eq.1) then
               flav=2*eq(j)**2*delta(i,j)*siggamma
            elseif (prodflag.eq.2) then
               flav=1/(2*xw)*ckm(i,j)**2*sigw
            elseif (prodflag.eq.21) then
               pdf=0
               pdf=(fh1(i)*fh2(j)+fh1(-j)*fh2(-i))*
     /              (delta(i,1)+delta(i,4))*
     /              (delta(j,2)+delta(j,3)+delta(j,5)) 
               flav=1/(2*xw)*ckm(i,j)**2*sigw
            elseif (prodflag.eq.22) then
               pdf=0
               pdf=(fh1(j)*fh2(i)+fh1(-i)*fh2(-j))*
     /              (delta(i,1)+delta(i,4))*
     /              (delta(j,2)+delta(j,3)+delta(j,5)) 
               flav=1/(2*xw)*ckm(i,j)**2*sigw
            elseif (prodflag.eq.3) then
               flav=((1/(2*sw*cw)*tau3(i,j)-delta(i,j)*eq(j)*sw/cw)**2+
     /              (-delta(i,j)*eq(j)*sw/cw)**2)*sigz
            elseif (prodflag.eq.5) then
               flav=2*eq(j)**2*delta(i,j)*siggamma+
     /((1/(2*sw*cw)*tau3(i,j)-delta(i,j)*eq(j)*sw/cw)**2+
     /              (-delta(i,j)*eq(j)*sw/cw)**2)*sigz
            endif
            xlumqq=xlumqq+pdf*flav
         enddo
      enddo

c.....Ead and Ebc qq luminosity (here q = q')
      xlumqqead=0
      do i=1,nf
         pdf=0
         pdf=fh1(i)*fh2(i)+fh1(-i)*fh2(-i)
         do j=1,nf
            flav=0
            if (prodflag.eq.1) then
               flav=2*eq(j)**2*delta(i,j)*siggamma
            elseif (prodflag.eq.2) then
               flav=1/(2*xw)*delta(i,j)*sumckm(i)*sigw
            elseif (prodflag.eq.21) then
               pdf=0
               pdf=(fh1(i)*fh2(i)*sumckm(i)*
     /              (delta(i,1)+delta(i,4))+
     /              fh2(-j)*fh1(-j)*sumckm(j)*
     /              (delta(j,2)+delta(j,3)+delta(j,5))) 
               flav=1/(2*xw)*delta(i,j)*sigw
            elseif (prodflag.eq.22) then
               pdf=0
               pdf=(fh1(j)*fh2(j)*sumckm(j)*
     /              (delta(j,2)+delta(j,3)+delta(j,5))+
     /              fh2(-i)*fh1(-i)*sumckm(i)*
     /              (delta(i,1)+delta(i,4))) 
               flav=1/(2*xw)*delta(i,j)*sigw
            elseif (prodflag.eq.3) then
               flav=((1/(2*sw*cw)*tau3(i,j)-delta(i,j)*eq(j)*sw/cw)**2+
     /              (-delta(i,j)*eq(j)*sw/cw)**2)*sigz
            elseif (prodflag.eq.5) then
               flav=2*eq(j)**2*delta(i,j)*siggamma+
     /((1/(2*sw*cw)*tau3(i,j)-delta(i,j)*eq(j)*sw/cw)**2+
     /              (-delta(i,j)*eq(j)*sw/cw)**2)*sigz
            endif
            xlumqqead=xlumqqead+pdf*flav
         enddo
      enddo      

c     Need to add missing luminosities of Eq. (2.21) of [Gonsalves, Pawlowsky, Wai]
c.....Eaa and Ecc qq luminosity (to be implemented)
      xlumqqeaa=0d0
c      do i=1,nf
c         do j=1,nf
c            pdf=0
c            pdf=(fh1(i)+fh1(i))*(fh2(j)+fh2(j))
c            flav=0
c            if (prodflag.eq.1) then
c               flav=2*eq(j)**2*delta(i,j)*siggamma
c            elseif (prodflag.eq.2) then
c               flav=1/(2*xw)*delta(i,j)*sumckm(i)*sigw
c            elseif (prodflag.eq.21) then
c               pdf=0
c               pdf=(fh1(i)*fh2(i)*sumckm(i)*
c     /              (delta(i,1)+delta(i,4))+
c     /              fh2(-j)*fh1(-j)*sumckm(j)*
c     /              (delta(j,2)+delta(j,3)+delta(j,5))) 
c               flav=1/(2*xw)*delta(i,j)*sigw
c            elseif (prodflag.eq.22) then
c               pdf=0
c               pdf=(fh1(j)*fh2(j)*sumckm(j)*
c     /              (delta(j,2)+delta(j,3)+delta(j,5))+
c     /              fh2(-i)*fh1(-i)*sumckm(i)*
c     /              (delta(i,1)+delta(i,4))) 
c               flav=1/(2*xw)*delta(i,j)*sigw
c            elseif (prodflag.eq.3) then
c               flav=((1/(2*sw*cw)*tau3(i,j)-delta(i,j)*eq(j)*sw/cw)**2+
c     /              (-delta(i,j)*eq(j)*sw/cw)**2)*sigz
c            elseif (prodflag.eq.5) then
c               flav=2*eq(j)**2*delta(i,j)*siggamma+
c     /((1/(2*sw*cw)*tau3(i,j)-delta(i,j)*eq(j)*sw/cw)**2+
c     /              (-delta(i,j)*eq(j)*sw/cw)**2)*sigz
c            endif
c            xlumqqeaa=xlumqqeaa+pdf*flav
c         enddo
c      enddo      

      xlumqqLL=0d0
      xlumqqLR=0d0

      return
      end
