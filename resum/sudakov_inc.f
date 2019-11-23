      integer flagrealcomplex,imod,iord
      common/iorder/iord
      common/flagrealcomplex/flagrealcomplex
      common/modified/imod

      double precision aass
      common/aass/aass
!$OMP THREADPRIVATE(/aass/)

      double precision a_param,b0p
      common/a_param/a_param,b0p
!$OMP THREADPRIVATE(/a_param/)

      double complex aexp,aexpB,aexpC,aexpD
      common/exponent/aexp,aexpB,aexpC,aexpD
!$OMP THREADPRIVATE(/exponent/)

      double precision rloga,rlogq2mur2
      common/rlogs/rloga,rlogq2mur2
!$OMP THREADPRIVATE(/rlogs/)

      double precision rblim
      complex *16 cblim
      common/blimit/rblim,cblim
!$OMP THREADPRIVATE(/blimit/)
      
      complex *16 log1y
      common/csud/log1y
!$OMP THREADPRIVATE(/csud/)
