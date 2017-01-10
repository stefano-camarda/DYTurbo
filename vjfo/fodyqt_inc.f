      double precision as
      common/asNEW/as
!$OMP THREADPRIVATE(/asNEW/)

      double precision asp
      common/asp/asp
!$OMP THREADPRIVATE(/asp/)
      
      double precision siggamma,sigint,sigz,sigw
      common/sigs/siggamma,sigint,sigz,sigw
!$OMP THREADPRIVATE(/sigs/)

      double precision vud,vus,vub,vcd,vcs,vcb,vtd,vts,vtb
      common/ckm/vud,vus,vub,vcd,vcs,vcb,vtd,vts,vtb

      double precision eq(5),alq(5),arq(5),ckm(6,6),delta(5,5),tau3(5,5)
      common/quarks/eq,alq,arq,ckm,delta,tau3
