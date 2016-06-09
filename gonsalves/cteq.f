*     CTEQ Parton Densities
      SUBROUTINE cteq(nset,x,qsq,pden)
      IMPLICIT NONE
      INTEGER nset,mode,i
      DOUBLE PRECISION x,qsq,pden(-6:6)
      DOUBLE PRECISION q,upv,dnv,usea,dsea,str,chm,bot,glu
      DOUBLE PRECISION Ctq6Pdf
      EXTERNAL Ctq6Pdf

      q=dsqrt(qsq)
      IF (nset.EQ.1) THEN
          mode=1
          CALL SetCtq6(mode)
          DO i=-5,5
              pden(i) = Ctq6Pdf(i,x,q)
          ENDDO
          pden(6) = 0d0
          pden(-6) = 0d0
      ELSE
          PRINT *,' bad CTEQ mode number',nset
          STOP
      ENDIF

      RETURN     
      END
