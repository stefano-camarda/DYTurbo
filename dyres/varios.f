       subroutine evaluatekin(eta,qt,q2)
       IMPLICIT DOUBLE PRECISION (A-I,L-Z)
       INTEGER IFIT,imod
       COMMON/transverse/qtt
       COMMON/ IFIT/ IFIT
!       COMMON/logs/xlog0,xlog1,xlog2,xlog3
       COMMON/modified/imod
       qtt=qt
       qt2=qt**2
C       
       IF (eta.le.-4.50001d0) then
            IFIT=14
       ELSEIF (eta.le.-4.0001d0) then
            IFIT=13
       ELSEIF (eta.le.-3.50001d0) then
            IFIT=12
       ELSEIF (eta.le.-3.0001d0) then
            IFIT=11
       ELSEIF (eta.le.-2.0001d0) then
            IFIT=10
       ELSEIF (eta.le.-1.0001d0) then
            IFIT=9
       ELSEIF (eta.lt.0.0001d0) then
            IFIT=8
       ELSEIF (eta.lt.1.0001) then
            IFIT=1
       ELSEIF (eta.lt.2.0001) then
            IFIT=2
       ELSEIF (eta.lt.3.0001) then
            IFIT=3
       ELSEIF (eta.lt.3.50001d0) then
            IFIT=4
       ELSEIF (eta.lt.4.0001) then
            IFIT=5
       ELSEIF (eta.lt.4.50001d0) then
            IFIT=6
       ELSEIF (eta.lt.5.0001d0) then
            IFIT=7
       ELSE
!        write(6,*) 'ETA too large', eta
!        STOP
           IFIT=7
       ENDIF

!      if (imod.eq.1) then
!            xlog0=((-qt2/2)* IK(1) )*4/qt2
!            xlog1=((-qt2/4) * Ik(2) )*4/qt2
!            xlog2=((-qt2/6) * IK(3) )*4/qt2
!            xlog3=( -(qt2/8)* IK(4) )*4/qt2 ! + 4*1.2020569d0*xlog0 no va!!   
!      elseif (imod.eq.0) then
!            xlog0=1d0 *4/qt2
!            xlog1=DLOG(q2/qt2) *4/qt2
!            xlog2=(DLOG(q2/qt2))**2 *4/qt2
!            xlog3=(DLOG(q2/qt2))**3 *4/qt2
!      else
!           write(6,*)'wrong imod= ',imod
!           stop
!      endif       
       
c       IFIT=1
       RETURN
       END
       

      subroutine fiteador(xtauf,muf2,etam)
       IMPLICIT DOUBLE PRECISION (A-I,L-Z)
c
       dimension A1UV(30),A2UV(30),A3UV(30),A4UV(30),
     .       A5UV(30),A6UV(30),A7UV(30),A8UV(30)
       dimension A1DV(30),A2DV(30),A3DV(30),A4DV(30),
     .       A5DV(30),A6DV(30),A7DV(30),A8DV(30)
       dimension A1US(30),A2US(30),A3US(30),A4US(30),
     .       A5US(30),A6US(30),A7US(30),A8US(30)
       dimension A1DS(30),A2DS(30),A3DS(30),A4DS(30),
     .       A5DS(30),A6DS(30),A7DS(30),A8DS(30)
       dimension A1SS(30),A2SS(30),A3SS(30),A4SS(30),
     .       A5SS(30),A6SS(30),A7SS(30),A8SS(30)
       dimension A1GL(30),A2GL(30),A3GL(30),A4GL(30),
     .       A5GL(30),A6GL(30),A7GL(30),A8GL(30)
       dimension A1CH(30),A2CH(30),A3CH(30),A4CH(30),
     .       A5CH(30),A6CH(30),A7CH(30),A8CH(30)
       dimension A1BO(30),A2BO(30),A3BO(30),A4BO(30),
     .       A5BO(30),A6BO(30),A7BO(30),A8BO(30)
c
      dimension A1UVp(30),A2UVp(30),A3UVp(30),A4UVp(30),
     .       A5UVp(30),A6UVp(30),A7UVp(30),A8UVp(30)
       dimension A1DVp(30),A2DVp(30),A3DVp(30),A4DVp(30),
     .       A5DVp(30),A6DVp(30),A7DVp(30),A8DVp(30)
       dimension A1USp(30),A2USp(30),A3USp(30),A4USp(30),
     .       A5USp(30),A6USp(30),A7USp(30),A8USp(30)
       dimension A1DSp(30),A2DSp(30),A3DSp(30),A4DSp(30),
     .       A5DSp(30),A6DSp(30),A7DSp(30),A8DSp(30)
       dimension A1SSp(30),A2SSp(30),A3SSp(30),A4SSp(30),
     .       A5SSp(30),A6SSp(30),A7SSp(30),A8SSp(30)
       dimension A1GLp(30),A2GLp(30),A3GLp(30),A4GLp(30),
     .       A5GLp(30),A6GLp(30),A7GLp(30),A8GLp(30)
       dimension A1CHp(30),A2CHp(30),A3CHp(30),A4CHp(30),
     .       A5CHp(30),A6CHp(30),A7CHp(30),A8CHp(30)
       dimension A1BOp(30),A2BOp(30),A3BOp(30),A4BOp(30),
     .       A5BOp(30),A6BOp(30),A7BOp(30),A8BOp(30)


       DIMENSION yav(30)
       INTEGER NFITMAX
       integer j

c
c beam 1
      common/ CUV1/ A1UV,A2UV,A3UV,A4UV,A5UV,A6UV,A7UV,A8UV
      common/ CDV1/ A1DV,A2DV,A3DV,A4DV,A5DV,A6DV,A7DV,A8DV
      common/ CUS1/ A1US,A2US,A3US,A4US,A5US,A6US,A7US,A8US
      common/ CDS1/ A1DS,A2DS,A3DS,A4DS,A5DS,A6DS,A7DS,A8DS
      common/ CSS1/ A1SS,A2SS,A3SS,A4SS,A5SS,A6SS,A7SS,A8SS
      common/ CGL1/ A1GL,A2GL,A3GL,A4GL,A5GL,A6GL,A7GL,A8GL
      common/ CCH1/ A1CH,A2CH,A3CH,A4CH,A5CH,A6CH,A7CH,A8CH
      common/ CBO1/ A1BO,A2BO,A3BO,A4BO,A5BO,A6BO,A7BO,A8BO
c beam 2
      common/ CUV2/ A1UVp,A2UVp,A3UVp,A4UVp,A5UVp,A6UVp,A7UVp,A8UVp
      common/ CDV2/ A1DVp,A2DVp,A3DVp,A4DVp,A5DVp,A6DVp,A7DVp,A8DVp
      common/ CUS2/ A1USp,A2USp,A3USp,A4USp,A5USp,A6USp,A7USp,A8USp
      common/ CDS2/ A1DSp,A2DSp,A3DSp,A4DSp,A5DSp,A6DSp,A7DSp,A8DSp
      common/ CSS2/ A1SSp,A2SSp,A3SSp,A4SSp,A5SSp,A6SSp,A7SSp,A8SSp
      common/ CGL2/ A1GLp,A2GLp,A3GLp,A4GLp,A5GLp,A6GLp,A7GLp,A8GLp
      common/ CCH2/ A1CHp,A2CHp,A3CHp,A4CHp,A5CHp,A6CHp,A7CHp,A8CHp
      common/ CBO2/ A1BOp,A2BOp,A3BOp,A4BOp,A5BOp,A6BOp,A7BOp,A8BOp
      common/expp/aa
      common/NFITMAX/NFITMAX

      LOGICAL fileexist


!CM       AA=2.5d0
       AA=3d0
       
       yav(1)=0.5d0 
       yav(2)=1.5d0 
       yav(3)=2.5d0 
       yav(4)=3.25d0 
       yav(5)=3.75d0 
       yav(6)=4.d0 
       yav(7)=4.5d0 
       if (etam.lt.1.0001d0) yav(1)=0d0
       if (etam.lt.2.0001d0) yav(2)=yav(1)
       if (etam.lt.3.0001d0) yav(3)=yav(2)
       if (etam.lt.3.50001d0) yav(4)=yav(3)
       if (etam.lt.4.0001d0) yav(5)=yav(4)
       if (etam.lt.4.50001d0) yav(6)=yav(5)
       if (etam.lt.5.0001d0) yav(7)=yav(6)

       yav(8)=-0.5d0 
       yav(9)=-1.5d0 
       yav(10)=-2.5d0 
       yav(11)=-3.25d0 
       yav(12)=-3.75d0 
       yav(13)=-4.d0 
       yav(14)=-4.5d0 
       if (etam.le.1.0001d0) yav(8)=0d0
       if (etam.le.2.0001d0) yav(9)=yav(8)
       if (etam.le.3.0001d0) yav(10)=yav(9)
       if (etam.le.3.50001d0) yav(11)=yav(10)
       if (etam.le.4.0001d0) yav(12)=yav(11)
       if (etam.le.4.50001d0) yav(13)=yav(12)
       if (etam.le.5.0001d0) yav(14)=yav(13)

        NFITMAX=14

c        do ji=1,14
c        yav(ji)=0d0
c        enddo

c        inquire( FILE='pdffit.grid', EXIST=fileexist ) 
c        if ( fileexist ) then
c           open(unit=101,file='pdffit.grid',status='unknown')
c           do kk=1,NFITMAX
c              read(101,*) A1UV(kk),A2UV(kk),A3UV(kk),
c     .             A4UV(kk),A5UV(kk),A6UV(kk),A7UV(kk),A8UV(kk)
c              read(101,*) A1DV(kk),A2DV(kk),A3DV(kk),
c     .             A4DV(kk),A5DV(kk),A6DV(kk),A7DV(kk),A8DV(kk)
c              read(101,*) A1US(kk),A2US(kk),A3US(kk),
c     .             A4US(kk),A5US(kk),A6US(kk),A7US(kk),A8US(kk)
c              read(101,*) A1DS(kk),A2DS(kk),A3DS(kk),
c     .             A4DS(kk),A5DS(kk),A6DS(kk),A7DS(kk),A8DS(kk)
c              read(101,*) A1SS(kk),A2SS(kk),A3SS(kk),
c     .             A4SS(kk),A5SS(kk),A6SS(kk),A7SS(kk),A8SS(kk)
c              read(101,*) A1GL(kk),A2GL(kk),A3GL(kk),
c     .             A4GL(kk),A5GL(kk),A6GL(kk),A7GL(kk),A8GL(kk)
c              read(101,*) A1CH(kk),A2CH(kk),A3CH(kk),
c     .             A4CH(kk),A5CH(kk),A6CH(kk),A7CH(kk),A8CH(kk)
c              read(101,*) A1BO(kk),A2BO(kk),A3BO(kk),
c     .             A4BO(kk),A5BO(kk),A6BO(kk),A7BO(kk),A8BO(kk)
c
c              read(101,*) A1UVp(kk),A2UVp(kk),A3UVp(kk),
c     .             A4UVp(kk),A5UVp(kk),A6UVp(kk),A7UVp(kk),A8UVp(kk)
c              read(101,*) A1DVp(kk),A2DVp(kk),A3DVp(kk),
c     .             A4DVp(kk),A5DVp(kk),A6DVp(kk),A7DVp(kk),A8DVp(kk)
c              read(101,*) A1USp(kk),A2USp(kk),A3USp(kk),
c     .             A4USp(kk),A5USp(kk),A6USp(kk),A7USp(kk),A8USp(kk)
c              read(101,*) A1DSp(kk),A2DSp(kk),A3DSp(kk),
c     .             A4DSp(kk),A5DSp(kk),A6DSp(kk),A7DSp(kk),A8DSp(kk)
c              read(101,*) A1SSp(kk),A2SSp(kk),A3SSp(kk),
c     .             A4SSp(kk),A5SSp(kk),A6SSp(kk),A7SSp(kk),A8SSp(kk)
c              read(101,*) A1GLp(kk),A2GLp(kk),A3GLp(kk),
c     .             A4GLp(kk),A5GLp(kk),A6GLp(kk),A7GLp(kk),A8GLp(kk)
c              read(101,*) A1CHp(kk),A2CHp(kk),A3CHp(kk),
c     .             A4CHp(kk),A5CHp(kk),A6CHp(kk),A7CHp(kk),A8CHp(kk)
c              read(101,*) A1BOp(kk),A2BOp(kk),A3BOp(kk),
c     .             A4BOp(kk),A5BOp(kk),A6BOp(kk),A7BOp(kk),A8BOp(kk)
c           enddo
c           close(101)
c           return
c        endif

         write(*,*)'Waiting for PDF fit ...'
         
        do kk=1,NFITMAX 
c        
        xtau1=xtauf**0.5* dexp(yav(kk))
        
        if (xtau1.gt.1d0) then
          write(6,*)xtauf**0.5,yav(kk),kk,etamax
          stop
          endif
c        
         call fiter(1,xtau1,muf2,A1UV(kk),A2UV(kk),A3UV(kk),
     .      A4UV(kk),A5UV(kk),A6UV(kk),A7UV(kk),A8UV(kk))
         call fiter(2,xtau1,muf2,A1DV(kk),A2DV(kk),A3DV(kk),
     .      A4DV(kk),A5DV(kk),A6DV(kk),A7DV(kk),A8DV(kk))
         call fiter(3,xtau1,muf2,A1US(kk),A2US(kk),A3US(kk),
     .      A4US(kk),A5US(kk),A6US(kk),A7US(kk),A8US(kk))
         call fiter(4,xtau1,muf2,A1DS(kk),A2DS(kk),A3DS(kk),
     .      A4DS(kk),A5DS(kk),A6DS(kk),A7DS(kk),A8DS(kk))
         call fiter(5,xtau1,muf2,A1SS(kk),A2SS(kk),A3SS(kk),
     .      A4SS(kk),A5SS(kk),A6SS(kk),A7SS(kk),A8SS(kk))
         call fiter(6,xtau1,muf2,A1GL(kk),A2GL(kk),A3GL(kk),
     .      A4GL(kk),A5GL(kk),A6GL(kk),A7GL(kk),A8GL(kk))
         call fiter(7,xtau1,muf2,A1CH(kk),A2CH(kk),A3CH(kk),
     .      A4CH(kk),A5CH(kk),A6CH(kk),A7CH(kk),A8CH(kk))
         call fiter(8,xtau1,muf2,A1BO(kk),A2BO(kk),A3BO(kk),
     .      A4BO(kk),A5BO(kk),A6BO(kk),A7BO(kk),A8BO(kk))
c
        xtau2=xtauf**0.5* dexp(-yav(kk))

        if (xtau2.gt.1d0) then
          write(6,*)xtauf**0.5,-yav(kk),kk,etamax
          stop
          endif          
c
         call fiter(1,xtau2,muf2,A1UVp(kk),A2UVp(kk),A3UVp(kk),
     .      A4UVp(kk),A5UVp(kk),A6UVp(kk),A7UVp(kk),A8UVp(kk))
         call fiter(2,xtau2,muf2,A1DVp(kk),A2DVp(kk),A3DVp(kk),
     .      A4DVp(kk),A5DVp(kk),A6DVp(kk),A7DVp(kk),A8DVp(kk))
         call fiter(3,xtau2,muf2,A1USp(kk),A2USp(kk),A3USp(kk),
     .      A4USp(kk),A5USp(kk),A6USp(kk),A7USp(kk),A8USp(kk))
         call fiter(4,xtau2,muf2,A1DSp(kk),A2DSp(kk),A3DSp(kk),
     .      A4DSp(kk),A5DSp(kk),A6DSp(kk),A7DSp(kk),A8DSp(kk))
         call fiter(5,xtau2,muf2,A1SSp(kk),A2SSp(kk),A3SSp(kk),
     .      A4SSp(kk),A5SSp(kk),A6SSp(kk),A7SSp(kk),A8SSp(kk))
         call fiter(6,xtau2,muf2,A1GLp(kk),A2GLp(kk),A3GLp(kk),
     .      A4GLp(kk),A5GLp(kk),A6GLp(kk),A7GLp(kk),A8GLp(kk))
         call fiter(7,xtau2,muf2,A1CHp(kk),A2CHp(kk),A3CHp(kk),
     .      A4CHp(kk),A5CHp(kk),A6CHp(kk),A7CHp(kk),A8CHp(kk))
         call fiter(8,xtau2,muf2,A1BOp(kk),A2BOp(kk),A3BOp(kk),
     .      A4BOp(kk),A5BOp(kk),A6BOp(kk),A7BOp(kk),A8BOp(kk))
c
        enddo 
        write(*,*)'PDF fit ended'
c     Dump PDF fit
        open(unit=101,file='pdffit.grid',status='unknown')
        do kk=1,NFITMAX
           write(101,*) A1UV(kk),A2UV(kk),A3UV(kk),
     .          A4UV(kk),A5UV(kk),A6UV(kk),A7UV(kk),A8UV(kk)
           write(101,*) A1DV(kk),A2DV(kk),A3DV(kk),
     .          A4DV(kk),A5DV(kk),A6DV(kk),A7DV(kk),A8DV(kk)
           write(101,*) A1US(kk),A2US(kk),A3US(kk),
     .          A4US(kk),A5US(kk),A6US(kk),A7US(kk),A8US(kk)
           write(101,*) A1DS(kk),A2DS(kk),A3DS(kk),
     .          A4DS(kk),A5DS(kk),A6DS(kk),A7DS(kk),A8DS(kk)
           write(101,*) A1SS(kk),A2SS(kk),A3SS(kk),
     .          A4SS(kk),A5SS(kk),A6SS(kk),A7SS(kk),A8SS(kk)
           write(101,*) A1GL(kk),A2GL(kk),A3GL(kk),
     .          A4GL(kk),A5GL(kk),A6GL(kk),A7GL(kk),A8GL(kk)
           write(101,*) A1CH(kk),A2CH(kk),A3CH(kk),
     .          A4CH(kk),A5CH(kk),A6CH(kk),A7CH(kk),A8CH(kk)
           write(101,*) A1BO(kk),A2BO(kk),A3BO(kk),
     .          A4BO(kk),A5BO(kk),A6BO(kk),A7BO(kk),A8BO(kk)

           write(101,*) A1UVp(kk),A2UVp(kk),A3UVp(kk),
     .          A4UVp(kk),A5UVp(kk),A6UVp(kk),A7UVp(kk),A8UVp(kk)
           write(101,*) A1DVp(kk),A2DVp(kk),A3DVp(kk),
     .          A4DVp(kk),A5DVp(kk),A6DVp(kk),A7DVp(kk),A8DVp(kk)
           write(101,*) A1USp(kk),A2USp(kk),A3USp(kk),
     .          A4USp(kk),A5USp(kk),A6USp(kk),A7USp(kk),A8USp(kk)
           write(101,*) A1DSp(kk),A2DSp(kk),A3DSp(kk),
     .          A4DSp(kk),A5DSp(kk),A6DSp(kk),A7DSp(kk),A8DSp(kk)
           write(101,*) A1SSp(kk),A2SSp(kk),A3SSp(kk),
     .          A4SSp(kk),A5SSp(kk),A6SSp(kk),A7SSp(kk),A8SSp(kk)
           write(101,*) A1GLp(kk),A2GLp(kk),A3GLp(kk),
     .          A4GLp(kk),A5GLp(kk),A6GLp(kk),A7GLp(kk),A8GLp(kk)
           write(101,*) A1CHp(kk),A2CHp(kk),A3CHp(kk),
     .          A4CHp(kk),A5CHp(kk),A6CHp(kk),A7CHp(kk),A8CHp(kk)
           write(101,*) A1BOp(kk),A2BOp(kk),A3BOp(kk),
     .          A4BOp(kk),A5BOp(kk),A6BOp(kk),A7BOp(kk),A8BOp(kk)
        enddo
        close(101)

c beam 1
c      common/ CUV1/ A1UV,A2UV,A3UV,A4UV,A5UV,A6UV,A7UV,A8UV
c      common/ CDV1/ A1DV,A2DV,A3DV,A4DV,A5DV,A6DV,A7DV,A8DV
c      common/ CUS1/ A1US,A2US,A3US,A4US,A5US,A6US,A7US,A8US
c      common/ CDS1/ A1DS,A2DS,A3DS,A4DS,A5DS,A6DS,A7DS,A8DS
c      common/ CSS1/ A1SS,A2SS,A3SS,A4SS,A5SS,A6SS,A7SS,A8SS
c      common/ CGL1/ A1GL,A2GL,A3GL,A4GL,A5GL,A6GL,A7GL,A8GL
c      common/ CCH1/ A1CH,A2CH,A3CH,A4CH,A5CH,A6CH,A7CH,A8CH
c      common/ CBO1/ A1BO,A2BO,A3BO,A4BO,A5BO,A6BO,A7BO,A8BO
cc beam 2
c      common/ CUV2/ A1UVp,A2UVp,A3UVp,A4UVp,A5UVp,A6UVp,A7UVp,A8UVp
c      common/ CDV2/ A1DVp,A2DVp,A3DVp,A4DVp,A5DVp,A6DVp,A7DVp,A8DVp
c      common/ CUS2/ A1USp,A2USp,A3USp,A4USp,A5USp,A6USp,A7USp,A8USp
c      common/ CDS2/ A1DSp,A2DSp,A3DSp,A4DSp,A5DSp,A6DSp,A7DSp,A8DSp
c      common/ CSS2/ A1SSp,A2SSp,A3SSp,A4SSp,A5SSp,A6SSp,A7SSp,A8SSp
c      common/ CGL2/ A1GLp,A2GLp,A3GLp,A4GLp,A5GLp,A6GLp,A7GLp,A8GLp
c      common/ CCH2/ A1CHp,A2CHp,A3CHp,A4CHp,A5CHp,A6CHp,A7CHp,A8CHp
c      common/ CBO2/ A1BOp,A2BOp,A3BOp,A4BOp,A5BOp,A6BOp,A7BOp,A8BOp


       
      return
      end


      subroutine fiter(idist,x,q2,a1,a2,a3,a4,a5,a6,a7,a8)
      implicit double precision (A-H,O-Z)
      external fcng1,chi2
      dimension nprm(8), vstrt(8),stp(8),bl(8),bu(8),arglis(16) 
      dimension parsal(8)
      character*10 pnam(8) 
      data nprm /1,2,3,4,5,6,7,8/  
      data pnam  / 'A1','A2','A3','A4','A5', 'A6', 'A7', 'A8'/
      data vstrt /  1D0, -.1D0, 5D0, 3D0, -1D0, 0D0, 0D0, 0D0/
      data stp   /  1D0, 1D0, 1D0, 10D0, 15D0, 25D0, 25D0, 25D0/
c      data bl    /0D0, -0.8D0, 4D0, -20D0, -40D0, -80D0,-80D0,-80D0/
c      data bu    /0D0, 1.5D0, 19D0, 600D0, 40D0, 80D0, 80D0, 80D0/

c      data bl    /0D0, -0.8D0, 5D0, -5D0, -30D0, -30D0,-30D0,-30D0/
c      data bu    /0D0, 1.2D0, 12D0, 50D0, 30D0, 30D0, 30D0, 30D0/

      data bl    /0D0, -0.8D0, 4D0, -20D0, -60D0, -90D0,-90D0,-90D0/
      data bu    /0D0, 1.5D0, 12D0, 600D0, 60D0, 90D0, 90D0, 90D0/

      common/sal/parsal
      common/distribucion/id
      common/xq2/rx,rq2
      id=idist
      rx=x
      rq2=q2
c.....initialization :
      call mninit(5,1,1)
c     call mninit(5,6,6)
c.....definitions of the parameters :
      do 11 i = 1, 8
         call mnparm (nprm(i),pnam(i),VSTRT(i),STP(i),BL(i),BU(i),
     .        ierflg,chi2)
         if (ierflg .ne. 0) then
            write (6,*) ' unable to define parameter no.', i
            stop
         end if
 11   continue
c.....output  
      arglis(1) = 0.
      call mnexcm(fcng1,'set print',arglis,1,ierflg,chi2)
c.....first call :
      arglis(1) = 1.            !   IFLAG = 1 
      call mnexcm(fcng1,'CALL FCN',arglis,1,ierflg,chi2)
c.....simplex fit :
   
c      goto 100

      arglis(1)=5000.
      call mnexcm(fcng1,'SIMPLEX',arglis,1,ierflg,chi2)
      call mnexcm(fcng1,'MIGRAD',arglis,1,ierflg,chi2)
      
      arglis(1)=20000.
      call mnexcm(fcng1,'SIMPLEX',arglis,1,ierflg,chi2)
      call mnexcm(fcng1,'MIGRAD',arglis,1,ierflg,chi2)
      arglis(1)=28000.
      call mnexcm(fcng1,'SIMPLEX',arglis,1,ierflg,chi2)
      call mnexcm(fcng1,'MIGRAD',arglis,1,ierflg,chi2)
 100  continue     

c      arglis(1)=500.
c      call mnexcm(fcng1,'SIMPLEX',arglis,1,ierflg,chi2)
c      call mnexcm(fcng1,'MIGRAD',arglis,1,ierflg,chi2)

cc.....last call :
      arglis(1) = 3             !   iflag = 3
      call mnexcm (fcng1, 'call fcn', arglis, 1, ierflg,chi2)
c.....stop :
      call mnexcm (fcng1,'stop',arglis,1,ierflg,chi2)
      a1=parsal(1)
      a2=parsal(2)
      a3=parsal(3)
      a4=parsal(4)
      a5=parsal(5)
      a6=parsal(6)
      a7=parsal(7)
      a8=parsal(8)
 1200 format(4F8.4)
      return
      end
      
      subroutine fcng1 (npar, g, f, x, iflag ,chi2)
      implicit double precision (a-h, o-z)
      dimension x(*), g(*)
      external chi2 
      f = chi2 (x)
      return
      end
            
      double precision function chi2(param)
      implicit double precision (a-h,o-z)
      dimension param(8),xx1(74),xx2(228)
      dimension parsal(8)
      common/sal/parsal
      common/distribucion/id
      common/xq2/rx,rq2
      common/expp/aa
      data npoints1/74/
      data npoints2/228/

      data xx1/1.d-5,2.d-5,4.d-5,6.d-5,8.d-5,9D-5,
     .     1.D-4,2.D-4,3.D-4,4.D-4,5.D-4,6.D-4,8.D-4,
     .     1.D-3,2.D-3,3D-3,4.D-3,5.D-3,6.D-3,7D-3,8d-3,9d-3,
     .     1.D-2,2.D-2,3D-2,4.D-2,5D-2,
     .     6.D-2,6.5d-2,7D-2,7.5d-2,8.D-2,8.5d-2,9d-2,9.5d-2,
     .     .1D0,.11d0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
     .     .5D0,.525D0,.55D0,.575D0,.6D0,
     .     .625d0,.65d0,.675d0,.7d0,.725d0,.75d0,.775d0,.8d0,
     .     .825d0,.85d0,.875d0,.9d0,.92d0,.94d0,.96d0,.98d0,1d0/

      data xx2/1.d-5,2.d-5,3.d-5,4.d-5,5.d-5,6.d-5,7.d-5,8.d-5,9D-5,
     . 1.D-4,1.5D-4,2.D-4,2.5D-4,3.D-4,3.5D-4,4.D-4,4.5D-4,5.D-4,
     . 5.5D-4,6.D-4,6.5D-4,7.D-4,7.5D-4,8.D-4,8.5D-4,9.D-4,9.5D-4,
     .0.001,0.0015,0.002,0.0025,0.003,0.0035,0.004,0.0045,0.005,
     .0.0055,0.006,0.0065,0.007,0.0075,0.008,0.0085,0.009,0.0095,
     .0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.050,
     .0.055,0.06,0.065,0.07,0.075,0.08,0.085,0.09,0.095,
     .0.1,0.105,0.11,0.115,0.12,0.125,0.13,0.135,    
     .0.14,0.145,0.15,0.155,0.16,0.165,0.17,0.175,0.18,   !80
     .0.185,0.19,0.195,0.20,0.205,    
     .0.21,0.215,0.22,0.225,0.23,0.235,    
     .0.24,0.245,0.25,0.255,0.26,0.265,0.27,0.275,0.28,0.285,0.29,0.295,    
     .0.30,0.305,0.31,0.315,0.32,0.325,0.33,0.335,0.34,0.345,0.35,0.355,
     .0.36,0.365,0.37,0.375,0.38,0.385,0.39,    
     .0.395,0.40,0.405,0.41,0.415,0.42,0.425,0.43,0.435,0.44,.445,0.45,  !134
     .0.455,0.46,0.465,0.47,0.475,0.48,0.485,0.49,0.495,0.50,0.505,    
     .0.51,0.515,0.52,0.525,0.53,0.535,0.54,0.545,0.55,0.555,0.56,
     .0.565,0.57,0.575,0.58,0.585,0.59,0.595,0.60,0.605,0.61,0.615,0.62,
     .0.625,0.63,0.635,0.64,0.645,0.65,0.655,0.66,0.665,0.67,0.675,0.68,  !180
     .0.685,0.69,0.70,0.705,0.71,0.715,0.72,0.725,0.73,0.735,0.74,
     .0.745,0.75,0.755,0.76,0.765,0.77,0.78,.785,
     .0.79,0.795,0.80,0.805,0.81,
     .0.815,0.82,0.825,0.83,0.835,
     .0.84,0.845,0.85,
     .0.855,0.86,0.865,0.87,0.875,
     .0.88,0.885,0.89,0.895,0.90,
     .0.91,0.93,0.95,0.97,0.99,1.0/! 228 

      a1=(param(1))
      a2=(param(2))
      a3=(param(3))
      a4=(param(4))
      a5=(param(5))
      a6=(param(6))
      a7=(param(7))
      a8=(param(8))
      chi2=0.d0
      
c.....aa can be changed to try to improve the fit
c.....the common changes it at the same time in the subroutine f0moments
!      aa=2.5
!      aa=3.5
      aa=3
      npoints=npoints2
c      aextra=1
c      if (id.eq.6) then
c          npoints=npoints1
c          aextra=30d0
c      endif
      
c.....take less/more points to improve speed/accuracy
      do i=1,npoints,1
         x=rx**( (1-xx2(i)))
         
c      if (id.eq.6) x=rx**(1-xx2(i))

         if (x.gt.1d0) goto 133
         
         f=a1*x**a2*(1-x)**a3*(1+a4*x+a5*x**(0.5)+a6*x**(1.5)+A7*X**2
     .        +A8*X**(aa))

         call distrit(x,sqrt(rq2),unv,dnv,us,ds,str,chm,bot,glu,1)
         if(id.eq.1) then
            f2=unv
         elseif(id.eq.2) then
            f2=dnv
         elseif(id.eq.3) then
            f2=us
         elseif(id.eq.4) then
            f2=ds
         elseif(id.eq.5) then
            f2=str
         elseif(id.eq.6) then
            f2=glu
         elseif(id.eq.7) then
            f2=chm
         elseif(id.eq.8) then
            f2=bot
         endif
       aaa=1d0    !ADDED
         if (x.lt.rx)   aaa=30d0  !ADDED
c      aextra=1
c      if (id.eq.6) then
c          npoints=npoints1                                                                                                  c           aextra=30d0
c      endif
         if(id.gt.2) then
         chi2=chi2+(f2-f)**2*aaa*x**1.2   !*aextra
         elseif(id.eq.1) then
         chi2=chi2+(f2-f)**2*300*aaa
         elseif(id.eq.2) then
         chi2=chi2+(f2-f)**2*300*aaa
         endif


 133   continue        
      enddo
      
      do i=1,npoints,1
         x=rx**( (1+xx2(i)))     

c      if (id.eq.6) x=rx**(1+xx2(i))

         if (x.lt.1d-5) goto 134
         f=a1*x**a2*(1-x)**a3*(1+a4*x+a5*x**(0.5)+a6*x**(1.5)+A7*X**2
     .        +A8*X**(aa))
         call distrit(x,sqrt(rq2),unv,dnv,us,ds,str,chm,bot,glu,1)
         if(id.eq.1) then
            f2=unv
         elseif(id.eq.2) then
            f2=dnv
         elseif(id.eq.3) then
            f2=us
         elseif(id.eq.4) then
            f2=ds
         elseif(id.eq.5) then
            f2=str
         elseif(id.eq.6) then
            f2=glu
         elseif(id.eq.7) then
            f2=chm
         elseif(id.eq.8) then
            f2=bot
         endif
        aaaa=2d0
         if (x.gt.rx)  aaaa=35d0

c      aextra=1
c      if (id.eq.6) then
c          npoints=npoints1
c          aextra=1d0
c      endif


         if(id.gt.2) then
         chi2=chi2+(f2-f)**2*aaaa*x**1.2  !*aextra
         elseif(id.eq.1) then
         chi2=chi2+(f2-f)**2*aaaa  !*300
         elseif(id.eq.2) then
         chi2=chi2+(f2-f)**2*aaaa  !*300  
         endif
 134   continue        
        enddo
      
C.....return the parameters to save the last set
      do i=1,8
         parsal(i)=param(i)
      enddo
      return
      end  
     
      
      subroutine distrit(x,q,upv,dnv,usea,dsea,str,chm,bot,glu,ippbar)
C.....here I call the new  sets !!!!! (in program prog_pdf.f)
C.....gives always x*distribution!!!!!
      real*8 upv,dnv,usea,dsea,str,chm,bot,glu,x,q
      common/isetproton/isetproton
      real*8 FX(-5:5),sq2,sx
      sq2=(q*q)
      sx=(x)
!      call partons(sq2,sx,fx,5,isetproton,ippbar)
      call fdist(ippbar,sx,q,fx)
!CM      call fittedpartons(sq2,sx,fx,5,isetproton,ippbar)
! CHANGED DEFINITION U <--> D
!      usea=(fx(-1))*x  
      usea=(fx(-2))*x  
!      dsea=(fx(-2))*x 
      dsea=(fx(-1))*x 
      str=(fx(-3))*x  
      chm=(fx(-4))*x  
      bot=(fx(-5))*x  
      glu=(fx(0))*x
!      upv=((fx(1))*x-usea )  
!      dnv=((fx(2))*x-dsea )  
      upv=((fx(2))*x-usea )  
      dnv=((fx(1))*x-dsea )  
      return
      end
       
      
       SUBROUTINE INITOFIT   
       IMPLICIT DOUBLE PRECISION (A - Z)
       INTEGER k,ik,NFITMAX
       COMPLEX*16 CCp,CCm, Np(136),Nm(136),XN
       COMPLEX*16 UV,DV,US,DS,SS,GL,CH,BO
       COMPLEX*16 uval,dval,usea,dsea,ssea,glu,charm,bot
      dimension UV(30),DV(30),US(30),DS(30),SS(30),GL(30),CH(30),BO(30)
      COMPLEX*16 MellinH2qq,MellinH2gg,MellinH2gq
       COMPLEX*16 UVP(136,30),DVP(136,30),USP(136,30),DSP(136,30),
     .         SSP(136,30),GLP(136,30),CHP(136,30),BOP(136,30)
       COMPLEX*16 UVM(136,30),DVM(136,30),USM(136,30),DSM(136,30),
     .           SSM(136,30),GLM(136,30),CHM(136,30),BOM(136,30)
       COMPLEX*16 UVP2(136,30),DVP2(136,30),USP2(136,30),DSP2(136,30),
     .           SSP2(136,30),GLP2(136,30),CHP2(136,30),BOP2(136,30)
       COMPLEX*16 UVM2(136,30),DVM2(136,30),USM2(136,30),DSM2(136,30),
     .           SSM2(136,30),GLM2(136,30),CHM2(136,30),BOM2(136,30)
     
       COMPLEX*16 QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F,
     1              QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F
     
       COMPLEX*16 QQIP(136),QGFP(136), GQIP(136), GGIP(136), GGFP(136),
     1            NS1MIP(136), NS1PIP(136), NS1FP(136),QQ1FP(136), 
     2       QG1FP(136), GQ1IP(136), GQ1FP(136), GG1IP(136), GG1FP(136)
 
       COMPLEX*16 QQIM(136),QGFM(136), GQIM(136), GGIM(136), GGFM(136),
     1           NS1MIM(136), NS1PIM(136), NS1FM(136),QQ1FM(136), 
     2       QG1FM(136), GQ1IM(136), GQ1FM(136), GG1IM(136), GG1FM(136)

c       COMPLEX*16 H2QQp(136),H2QGp(136),H2GQp(136),H2GGp(136)
c       COMPLEX*16 H2QQm(136),H2QGm(136),H2GQm(136),H2GGm(136)
       COMPLEX*16 C2qgMp(136),C2NSqqMp(136),C2SqqbMp(136),
     .            C2NSqqbMp(136)
       COMPLEX*16 C2qgMm(136),C2NSqqMm(136),C2SqqbMm(136),
     .            C2NSqqbMm(136)
       COMPLEX*16 C2qg,C2NSqq,C2Sqqb,C2NSqqb
     
       COMMON / ANOMP/QQIp, QGFp, GQIp, GGIp, GGFp, NS1MIp, NS1PIp, 
     1          NS1Fp, QQ1Fp, QG1Fp, GQ1Ip, GQ1Fp, GG1Ip, GG1Fp
       COMMON / ANOMM/QQIm, QGFm, GQIm, GGIm, GGFm, NS1MIm, NS1PIm, 
     1          NS1Fm, QQ1Fm, QG1Fm, GQ1Im, GQ1Fm, GG1Im, GG1Fm
       
       COMMON / DISTP1/ UVP,DVP,USP,DSP,SSP,GLP,CHP,BOP
       COMMON / DISTM1/ UVM,DVM,USM,DSM,SSM,GLM,CHM,BOM
       COMMON / DISTP2/ UVP2,DVP2,USP2,DSP2,SSP2,GLP2,CHP2,BOP2
       COMMON / DISTM2/ UVM2,DVM2,USM2,DSM2,SSM2,GLM2,CHM2,BOM2
       COMMON / MOMS2    / Np,Nm,CCP,CCm
       common/NFITMAX/NFITMAX

       COMMON / H2COEF /C2qgMp,C2NSqqMp,C2SqqbMp,C2NSqqbMp,
     .                  C2qgMm,C2NSqqMm,C2SqqbMm,C2NSqqbMm

c       COMMON / H2COEF /H2QQp,H2QGp,H2GQp,H2GGp,H2QQm,H2QGm,H2GQm,H2GGm
       include 'quadrules.f'

      Write(6,*)'Start initialization'

c Beam 1        
       do k=1,136
C POSITIVE BRANCH
       XN=Np(k)  
       call F0MOMENTS1(XN,UV,DV,US,DS,SS,GL,CH,BO)
       do ik=1,NFITMAX
        UVp(k,ik)= UV(ik)
        DVp(k,ik)= DV(ik)
        USp(k,ik)= US(ik)
        DSp(k,ik)= DS(ik)
        SSp(k,ik)= SS(ik)
        GLp(k,ik)= GL(ik)
        CHp(k,ik)= CH(ik)
        BOp(k,ik)= BO(ik)
       enddo 
        CALL ANCALC (QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F,
     1              QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, XN)
        QQIp(k) = QQI 
        QGFp(k) = QGF 
        GQIp(k) = GQI 
        GGIp(k) = GGI 
        GGFp(k) = GGF 
        NS1MIp(k) = NS1MI  
        NS1PIp(k) = NS1PI 
        NS1Fp(k) = NS1F            
        QQ1Fp(k) = QQ1F 
        QG1Fp(k) = QG1F 
        GQ1Ip(k) = GQ1I 
        GQ1Fp(k) = GQ1F 
        GG1Ip(k) = GG1I 
        GG1Fp(k) = GG1F 
        
C Compute H2 coefficients      
c        H2QQp(k)=MellinH2qq(XN) 
c        H2QGp(k)=0d0
c        H2GQp(k)=MellinH2gq(XN)
c        H2GGp(k)=MellinH2gg(XN) 

        call H2calc(C2qg,C2NSqqb,C2NSqq,C2Sqqb,XN)

        C2qgMp(k)= C2qg
        C2NSqqMp(k)= C2NSqq
        C2SqqbMp(k)= C2Sqqb
        C2NSqqbMp(k)=  C2NSqqb      
       
C NEGATIVE BRANCH
         XN=Nm(k)
       call F0MOMENTS1(XN,UV,DV,US,DS,SS,GL,CH,BO)
       do ik=1,NFITMAX
        UVm(k,ik)= UV(ik)
        DVm(k,ik)= DV(ik)
        USm(k,ik)= US(ik)
        DSm(k,ik)= DS(ik)
        SSm(k,ik)= SS(ik)
        GLm(k,ik)= GL(ik)
        CHm(k,ik)= CH(ik)
        BOm(k,ik)= BO(ik)
        enddo
        CALL ANCALC (QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F,
     1              QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, XN)
        QQIm(k) = QQI 
        QGFm(k) = QGF 
        GQIm(k) = GQI 
        GGIm(k) = GGI 
        GGFm(k) = GGF
        NS1MIm(k) = NS1MI  
        NS1PIm(k) = NS1PI 
        NS1Fm(k) = NS1F            
        QQ1Fm(k) = QQ1F 
        QG1Fm(k) = QG1F 
        GQ1Im(k) = GQ1I 
        GQ1Fm(k) = GQ1F 
        GG1Im(k) = GG1I 
        GG1Fm(k) = GG1F 

C Compute H2 coefficients      
c        H2QQm(k)=MellinH2qq(XN) 
c        H2QGm(k)=0d0 
c        H2GQm(k)=MellinH2gq(XN) 
c        H2GGm(k)=MellinH2gg(XN) 

       call H2calc(C2qg,C2NSqqb,C2NSqq,C2Sqqb,XN)

        C2qgMm(k)= C2qg
        C2NSqqMm(k)= C2NSqq
        C2SqqbMm(k)= C2Sqqb
        C2NSqqbMm(k)=  C2NSqqb              
        enddo

c Beam 2
       do k=1,136
C POSITIVE BRANCH
       XN=Np(k)  
       call F0MOMENTS2(XN,UV,DV,US,DS,SS,GL,CH,BO)
       do ik=1,NFITMAX
        UVp2(k,ik)= UV(ik)
        DVp2(k,ik)= DV(ik)
        USp2(k,ik)= US(ik)
        DSp2(k,ik)= DS(ik)
        SSp2(k,ik)= SS(ik)
        GLp2(k,ik)= GL(ik)
        CHp2(k,ik)= CH(ik)
        BOp2(k,ik)= BO(ik)
      enddo
C NEGATIVE BRANCH
      XN=Nm(k)
      call F0MOMENTS2(XN,UV,DV,US,DS,SS,GL,CH,BO)
       do ik=1,NFITMAX
        UVm2(k,ik)= UV(ik)
        DVm2(k,ik)= DV(ik)
        USm2(k,ik)= US(ik)
        DSm2(k,ik)= DS(ik)
        SSm2(k,ik)= SS(ik)
        GLm2(k,ik)= GL(ik)
        CHm2(k,ik)= CH(ik)
        BOm2(k,ik)= BO(ik)
        enddo
        enddo

        call CACHEANOM
        Write(6,*)'End initialization'
        return
        end

        
c n=8 gaussian quadrature
      SUBROUTINE INITO   
      IMPLICIT DOUBLE PRECISION (A - Z)
      DOUBLE PRECISION  WN(136),WZ(8),Zs(8),DOWN(17),UP(17)
      COMPLEX*16 CCp,CCm, Np(136),Nm(136)
      INTEGER NMAX, I1, I2
      INTEGER K, I3
      COMMON / WEIGHTS2 / WN
      COMMON / MOMS2    / Np,Nm,CCP,CCm
                                        
       DATA WZ
     1  / 0.10122 85362 90376,  0.22238 10344 53374, 
     2    0.31370 66458 77887,  0.36268 37833 78362, 
     3    0.36268 37833 78362,  0.31370 66458 77887,
     4    0.22238 10344 53374,  0.10122 85362 90376/
       DATA ZS
     1  /-0.96028 98564 97536, -0.79666 64774 13627, 
     2   -0.52553 24099 16329, -0.18343 46424 95650, 
     3    0.18343 46424 95650,  0.52553 24099 16329,
     4    0.79666 64774 13627,  0.96028 98564 97536/
*...INTEGRATION CONTOUR PARAMETERS :
       DATA DOWN / 0.D0, 0.5D0, 1.D0, 2.D0, 3.D0, 4.D0, 6.D0, 8.D0,
     1     1.D1, 1.2D1, 1.5D1, 1.8D1, 2.1D1, 2.4D1, 2.7D1, 3.D1, 3.3D1/
c       DATA DOWN/0d0,0.5d0,1d0,2d0,3d0,4d0,
c     .      6d0,8d0,10d0,14d0,18D0,22d0,
c     .      36d0,42d0,48d0,54d0,62d0/
       include 'quadrules.f'
C     
c     C is the starting point on the real axis for the positive and negative part of the integration path
c     D is a scaling factor for the total length of the integration path
c     PHI is the angle in the complex plane of the positive part of the integration path
       
c     C = 2.3d0       
c**********************************
       if (approxpdf.eq.1) then
c     Original settings
          C = 1d0
          D = 1d0
          PHI = 3.141592654 * 3./4.
       else
c     Modified settings to allow numerical integration of PDFs melling moment (need real part of moments always > 0)
          C = 1d0
          D = 1.5d0             !sqrt(2d0)
          PHI = 3.141592654 * 1./2.
       endif
c**********************************
       
       CO = D*DCOS (PHI)
       SI = D*DSIN (PHI)
       CCp = CMPLX (CO, SI)
       DO 31 I1 = 1, 16
         UP(I1) = DOWN(I1+1)
  31    CONTINUE
      UP(17) = 36.D0
*...SUPPORT POINTS AND WEIGHTS FOR THE GAUSS INTEGRATION : 
*    (THE FACTOR (UP-DOWN)/2 IS INCLUDED IN THE WEIGHTS)
       K = 0
       DO 2 I2 = 1, 17
         SUM  = UP(I2) + DOWN(I2) 
         DIFF = UP(I2) - DOWN(I2) 
       DO 3 I3 = 1, 8
         K = K + 1
         Z = 0.5 * (SUM + DIFF * ZS(I3))
         WN(K) = DIFF / 2.* WZ(I3) 
         Np(K)  = CMPLX (C+CO*Z+1.,SI*Z)
  3    CONTINUE
  2    CONTINUE 
       CO = D*DCOS (PHI)
       SI = -D*DSIN (PHI)
       CCm = CMPLX (CO, SI)
       DO 1 I1 = 1, 16
         UP(I1) = DOWN(I1+1)
  1    CONTINUE
      UP(17) = 36.D0
*...SUPPORT POINTS AND WEIGHTS FOR THE GAUSS INTEGRATION : 
*    (THE FACTOR (UP-DOWN)/2 IS INCLUDED IN THE WEIGHTS)
       K = 0
       DO 20 I2 = 1, 17
         SUM  = UP(I2) + DOWN(I2) 
         DIFF = UP(I2) - DOWN(I2) 
       DO 30 I3 = 1, 8
         K = K + 1
         Z = 0.5 * (SUM + DIFF * ZS(I3))
         Nm(K)  = CMPLX (C+CO*Z+1.,SI*Z)         
  30    CONTINUE
  20    CONTINUE 
        RETURN                           
        END

        
!      FUNCTION ADPINT (F, A, B, AERR, RERR, ERREST, IER)
!      IMPLICIT DOUBLE PRECISION (A-H, O-Z)     
!c.....Integral of F(X) from A to B, with error
!c.....less than ABS(AERR) + ABS(RERR*INTEGRAL)
!c.....Best estimate of error returned in ERREST.
!c.....Error code is IER: zero if OK, non-zero if in trouble.
!      INTEGER OLDINT
!      EXTERNAL F
!      PARAMETER (MAXINT = 500)
!     
!c     Work space:
!      COMMON / ADPWRK2 / NUMINT
!      COMMON / ADPWRK / U(MAXINT), V(MAXINT), FU(MAXINT),
!     >     FV(MAXINT), FW(MAXINT), ERR(MAXINT), RESULT(MAXINT)
!
!      SAVE / ADPWRK /
!      SAVE / ADPWRK2 /
!      
!      IER = 0
!      NUMINT = 5
!      DX = (B-A)/ NUMINT
!      DO 1  I = 1, NUMINT
!         IF (I .EQ. 1)  THEN
!            U(I) = A
!            FU(I) = F(U(I))
!         ELSE
!            U(I) = V(I-1)
!            FU(I) = FV(I-1)
!         ENDIF
!         IF (I .EQ. NUMINT) THEN
!            V(I) = B
!         ELSE
!            V(I) = A + DX * I
!         ENDIF
!         FV(I) = F(V(I))
!         CALL ADPCAL(F,I)
! 1    CONTINUE
!      
! 2    CONTINUE
!      
!c.....Error estimate:
!      
!      ADPINT = 0.
!      ERREST = 0.
!      DO 3  I = 1, NUMINT
!         ADPINT = ADPINT + RESULT(I)
!         ERREST = ERREST + ERR(I)
! 3    CONTINUE
!      TARGET = ABS(AERR) + ABS(RERR * ADPINT)
!      IF (ERREST .GT. TARGET)  THEN
!         OLDINT = NUMINT
!         DO 4 I = 1, OLDINT
!            IF (ERR(I)*2*OLDINT .GT. TARGET)  CALL ADPSPL(F,I,IER)
! 4       CONTINUE
!         IF (IER .EQ. 0)  GOTO 2
!      ENDIF
!      RETURN
!      END
!      
!      FUNCTION INTUSE ()
!      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
!
!c.....Return number of intervals used last call to ADPINT
!
!      PARAMETER (MAXINT = 500)
!      COMMON / ADPWRK2 / NUMINT
!      COMMON / ADPWRK / U(MAXINT), V(MAXINT), FU(MAXINT),
!     >     FV(MAXINT), FW(MAXINT), ERR(MAXINT), RESULT(MAXINT)
!
!      INTUSE = NUMINT
!      RETURN
!      END
!      
!      SUBROUTINE ADPCAL (F,I)
!
!c.....Fill in details of interval I given endpoints
!      
!      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
!      EXTERNAL F
!      PARAMETER (MAXINT = 500)
!      COMMON / ADPWRK2 / NUMINT
!      COMMON / ADPWRK / U(MAXINT), V(MAXINT), FU(MAXINT),
!     >     FV(MAXINT), FW(MAXINT), ERR(MAXINT), RESULT(MAXINT)
!      
!      FW(I) = F( (U(I) + V(I)) /2.)
!      DX = V(I) - U(I)
!      RESULT(I) = DX * (FU(I) + 4. * FW(I) + FV(I)) / 6.
!      ERR(I) = ABS(DX * (FU(I) - 2. * FW(I) + FV(I)) / 12.)
!      RETURN
!      END
!      
!      SUBROUTINE ADPSPL (F, I, IER)
!
!c.....Split interval I
!
!      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
!      EXTERNAL F
!      PARAMETER (MAXINT = 500)
!      COMMON / ADPWRK2 / NUMINT
!      COMMON / ADPWRK / U(MAXINT), V(MAXINT), FU(MAXINT),
!     >     FV(MAXINT), FW(MAXINT), ERR(MAXINT), RESULT(MAXINT)
!      
!      IF (NUMINT .GE. MAXINT)  THEN
!         IER = 1
!         RETURN
!      ENDIF
!      NUMINT = NUMINT + 1
!      V(NUMINT) = V(I)
!      U(NUMINT) = (U(I) + V(I)) / 2.
!      V(I) = U(NUMINT)
!      FV(NUMINT) = FV(I)
!      FU(NUMINT) = FW(I)
!      FV(I) = FW(I)
!      CALL ADPCAL (F, I)
!      CALL ADPCAL (F, NUMINT)
!      RETURN
!      END

c.....Function BK(n,z)
c.....BK(n,z) is the n-derivative of BesselK[nu,z]
c.....with respect to nu in nu=1

c     NOW b0->b0p INCLUDED

      function IK(m)
      implicit none
      external BK
      real *8 BK,IK,argum,dloqt
      integer m
      include 'scales_inc.f'
      include 'sudakov_inc.f'
      include 'const.h'

      argum=b0p*qt/q
      dloqt=DLOG(a_param*qt/q)
      
      if (m.eq.1) then
         IK=-(2*b0p/q/qt)*BK(0,argum)
      elseif (m.eq.2) then
         IK=-(4*b0p/q/qt)*(BK(1,argum)-BK(0,argum)*DLOqt)
      elseif (m.eq.3) then
         IK=(b0p/q/qt)*(-6*BK(2,argum)+12*BK(1,argum)*DLOqt+
     /                 BK(0,argum)*(pi**2-6*(DLOqt**2)))
      elseif (m.eq.4) then
         IK=-(4*b0p/q/qt)*(2*BK(3,argum)-6*DLOqt*BK(2,argum)+
     /                    BK(1,argum)*(6*(DLOqt)**2-pi**2)+
     /                    BK(0,argum)*(pi**2*DLOqt-
     /                                 2*(DLOqt)**3-4*Z3))
      endif
      return
      end

      function BK(n,z)
      implicit none
      external fb
      real *8 bk,fb,errest,z,zz,max,adpint
      integer n,nn,ifail
      common/nuorder/nn
      common/zz/zz
      nn=n
      zz=z
      max=10d0
      bk=adpint(fb,0d0,max,1d-10,1d-5,errest,ifail)
      return
      end
      
      
      function fb(t)
      implicit none
      integer nn,nu
      real *8 fb,t,zz
      common/nuorder/nn
      common/zz/zz
      nu=1
      if(nn.eq.0) then
         fb=dexp(-zz*dcosh(t))*dcosh(nu*t)
      elseif(nn.eq.1) then
         fb=dexp(-zz*dcosh(t))*t*dsinh(nu*t)
      elseif(nn.eq.2) then
         fb=dexp(-zz*dcosh(t))*t*t*dcosh(nu*t)
      elseif(nn.eq.3) then
         fb=dexp(-zz*dcosh(t))*t*t*t*dsinh(nu*t)
      endif      
      return
      end
      




