      complex *16 cfx1(-5:5,512), cfx2p(-5:5,512), cfx2m(-5:5,512)
      common/creno/cfx1,cfx2p,cfx2m
!$OMP THREADPRIVATE(/creno/)
