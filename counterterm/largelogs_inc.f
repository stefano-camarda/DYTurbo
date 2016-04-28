c     Large logs calculated in ctqtint
      double precision LL1jk(-nf:nf,-nf:nf),LL2jk(-nf:nf,-nf:nf)
      double precision LL3jk(-nf:nf,-nf:nf),LL4jk(-nf:nf,-nf:nf)
      common/largelogs/LL1jk,LL2jk,LL3jk,LL4jk
!$OMP THREADPRIVATE(/largelogs/)
