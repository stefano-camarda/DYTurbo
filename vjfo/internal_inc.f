      double precision qt,q,q2
      common/internal/qt,q,q2
!$OMP THREADPRIVATE(/internal/)

      double precision tm
      common/tm/tm
!$OMP THREADPRIVATE(/tm/)

      double precision sh,th,uh,s2
      common/mand/sh,th,uh,s2
!$OMP THREADPRIVATE(/mand/)

      double precision x1,x2
      common/fractions/x1,x2
!$OMP THREADPRIVATE(/fractions/)

      double precision asp
      common/asp/asp
!$OMP THREADPRIVATE(/asp/)
      
