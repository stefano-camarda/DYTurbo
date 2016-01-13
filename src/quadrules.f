      include 'options.f'
      integer pdfrule
      parameter (pdfrule=64)

c     quadrature rules for mellin inversions
!!!!! mdim=mellinrule*mellintervales has to be lower than 136
      integer mellinintervals,mellinrule,mdim
      parameter (mellinintervals=1)
      parameter (mellinrule=64)
      parameter (mdim=mellinrule*mellinintervals)

c     rapidity integral in resummed 2d integration
      integer yintervals,yrule
      parameter (yintervals=1)
      parameter (yrule=20)

c     alfa beta integral in counter term
      integer ctintervals,ctrule,ctdim
      parameter (ctintervals=1)
      parameter (ctrule=64)
      parameter (ctdim=ctrule*ctintervals)
c     common block initialised in ctquadinit
      double precision ctx(ctdim)
      double precision ctw(ctdim)
      common/ctweights/ctx,ctw

c     qt integral in counterterm
      integer qtintervals,qtrule
      parameter (qtintervals=1)
      parameter (qtrule=20)
