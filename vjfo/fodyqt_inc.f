      double precision as
      common/asNEW/as
!$OMP THREADPRIVATE(/asNEW/)

      double precision siggamma,sigint,sigz,sigw
      common/sigs/siggamma,sigint,sigz,sigw
!$OMP THREADPRIVATE(/sigs/)

      double precision yv,expyp,expym
      common/yv/yv,expyp,expym
!$OMP THREADPRIVATE(/yv/)
