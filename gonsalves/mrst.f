*     Martin Roberts Stirling Thorne Parton Densities
      SUBROUTINE mrst(nset,x,qsq,pden)
      IMPLICIT NONE
      INTEGER nset,mode
      DOUBLE PRECISION x,qsq,pden(-6:6)
      DOUBLE PRECISION q,upv,dnv,usea,dsea,str,chm,bot,glu

      q=dsqrt(qsq)
      IF (nset.EQ.0) THEN
          mode = 1
          CALL mrstlo(x,q,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
      ELSEIF (nset.EQ.1.OR.nset.EQ.2) THEN
          mode=nset
          CALL mrst2002(x,q,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
      ELSE
          PRINT *,' bad MRST mode number',nset
          STOP
      ENDIF
      pden(0) = glu/x
      pden(1) = (upv+usea)/x
      pden(-1) = usea/x
      pden(2) = (dnv+dsea)/x
      pden(-2) = dsea/x
      pden(3) = str/x
      pden(-3) = pden(3)
      pden(4) = chm/x
      pden(-4) = pden(4)
      pden(5) = bot/x
      pden(-5) = pden(5)
      pden(6) = 0d0
      pden(-6) = pden(6)

      RETURN     
      END
