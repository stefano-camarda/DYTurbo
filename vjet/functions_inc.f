      double precision dt,du,ds,dst,dsu,dtu,la
      common/denominators/dt,du,ds,dst,dsu,dtu,la
!$OMP THREADPRIVATE(/denominators/)

      double precision fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
      common/transcendentals/fs,ft,fu,fs2,fm2,fmu2,fa,fst,fsu,
     /     fla,fstu,ftu,flat,flau,f1t,f2t,f1u,f2u
!$OMP THREADPRIVATE(/transcendentals/)
      
      double precision s2
      common/recmass/s2
!$OMP THREADPRIVATE(/recmass/)
