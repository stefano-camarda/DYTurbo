       DOUBLE COMPLEX P0NS,P0SG
       DOUBLE COMPLEX P1NS,P1SG
       DOUBLE COMPLEX P2NS,P2SG
       DOUBLE COMPLEX A2NS,A2SG
       DOUBLE COMPLEX SSCHLP,SSTR2P,SSTR3P
       DOUBLE COMPLEX U1,R1,U1H
       DOUBLE COMPLEX U2,R2,U2H,UNS2
       DOUBLE COMPLEX R,E
       include 'dimensions.f'
*
* ---------------------------------------------------------------------
*
* Singlet and non-singlet splitting functions at LO,NLO,NNLO
*
       COMMON / PNS0   / P0NS (NDIM, NFMIN:NFMAX)
!$OMP THREADPRIVATE(/PNS0/)
       COMMON / PSG0   / P0SG (NDIM, NFMIN:NFMAX, 2, 2)
!$OMP THREADPRIVATE(/PSG0/)
       COMMON / PNS1   / P1NS (NDIM, NFMIN:NFMAX, 3)
!$OMP THREADPRIVATE(/PNS1/)
       COMMON / PSG1   / P1SG (NDIM, NFMIN:NFMAX, 2, 2)
!$OMP THREADPRIVATE(/PSG1/)
       COMMON / PNS2   / P2NS (NDIM, NFMIN:NFMAX, 3)
!$OMP THREADPRIVATE(/PNS2/)
       COMMON / PSG2   / P2SG (NDIM, NFMIN:NFMAX, 2, 2)
!$OMP THREADPRIVATE(/PSG2/)
       COMMON / SPSUMS / SSCHLP(NDIM), SSTR2P(NDIM), SSTR3P(NDIM)
!$OMP THREADPRIVATE(/SPSUMS/)
       COMMON / ANS2   / A2NS (NDIM)
!$OMP THREADPRIVATE(/ANS2/)
       COMMON / ASG2   / A2SG (NDIM, 2, 2)
!$OMP THREADPRIVATE(/ASG2/)

       COMMON / U1SG   / U1(NDIM, NFMIN:NFMAX, 2, 2)
!$OMP THREADPRIVATE(/U1SG/)
       COMMON / R1SG   / R1(NDIM, NFMIN:NFMAX, 2, 2)
!$OMP THREADPRIVATE(/R1SG/)
       COMMON / U1HSG  / U1H(NUMAX, NDIM, NFMIN:NFMAX, 2, 2)
!$OMP THREADPRIVATE(/U1HSG/)
       COMMON / U2NS   / UNS2(NUMAX, NDIM, NFMIN:NFMAX, 3)
!$OMP THREADPRIVATE(/U2NS/)
       COMMON / U2SG   / U2(NDIM, NFMIN:NFMAX, 2, 2)
!$OMP THREADPRIVATE(/U2SG/)
       COMMON / R2SG   / R2(NDIM, NFMIN:NFMAX, 2, 2)
!$OMP THREADPRIVATE(/R2SG/)
       COMMON / U2HSG  / U2H(NUMAX, NDIM, NFMIN:NFMAX, 2, 2)
!$OMP THREADPRIVATE(/U2HSG/)

       COMMON / LSG    / R(NDIM, NFMIN:NFMAX, 2),
     1                   E(NDIM, NFMIN:NFMAX, 2, 2, 2)
!$OMP THREADPRIVATE(/LSG/)
       
