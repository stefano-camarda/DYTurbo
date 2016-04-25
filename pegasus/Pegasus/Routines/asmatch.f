*
* ..File: asmatch.f  
*
*
* ..The threshold matching of the QCD coupling in the MS(bar) scheme,  
*    a_s = alpha_s(mu_r^2)/(4 pi),  for  NF -> NF + 1  active flavours 
*    up to order a_s^4 (NNNLO).
*
* ..The value  ASNF  of a_s for NF flavours at the matching scale, the 
*    logarithm  LOGRH = ln (mu_r^2/m_H^2) -- where m_H is the pole mass
*    of the heavy quark -- and  NF  are passed as arguments to the 
*    function  ASNF1.  The order of the expansion  NAORD  (defined as 
*    the 'n' in N^nLO) is provided by the common-block  ASPAR.
*
* ..The matching coefficients are inverted from Chetyrkin, Kniehl and
*    Steinhauser, Phys. Rev. Lett. 79 (1997) 2184. The QCD colour
*    factors have been hard-wired in these results. The lowest integer 
*    values of the Zeta function are given by the common-block  RZETA.
*
* =====================================================================
*
*
       FUNCTION ASNF1 (ASNF, LOGRH, NF)
*
       IMPLICIT DOUBLE PRECISION (A - Z)
       INTEGER NF, NAORD, NASTPS, PRVCLL, K1, K2
       DIMENSION CMC(3,0:3)
*
* ---------------------------------------------------------------------
*
* ..Input common-blocks 
*
       COMMON / ASPAR  / NAORD, NASTPS
       COMMON / RZETA  / ZETA(6)
*
* ..Variables to be saved for the next call
*
       SAVE CMC, CMCI30, CMCF30, CMCF31, CMCI31, PRVCLL
*
* ---------------------------------------------------------------------
*
* ..The coupling-constant matching coefficients (CMC's) up to NNNLO 
*   (calculated and saved in the first call of this routine)
*
       IF (PRVCLL .NE. 1) THEN
*
         CMC(1,0) =  0.D0
         CMC(1,1) =  2./3.D0
*
         CMC(2,0) = 14./3.D0
         CMC(2,1) = 38./3.D0
         CMC(2,2) =  4./9.D0  
*
         CMCI30 = + 80507./432.D0 * ZETA(3) + 58933./1944.D0 
     1            + 128./3.D0 * ZETA(2) * (1.+ DLOG(2.D0)/3.D0)
         CMCF30 = - 64./9.D0 * (ZETA(2) + 2479./3456.D0)
         CMCI31 =   8941./27.D0
         CMCF31 = - 409./27.D0
         CMC(3,2) = 511./9.D0
         CMC(3,3) = 8./27.D0
*
         PRVCLL = 1
*
       END IF
*
* ---------------------------------------------------------------------
*
* ..The N_f dependend CMC's, and the alpha_s matching at order NAORD 
*
       CMC(3,0) = CMCI30 + NF * CMCF30
       CMC(3,1) = CMCI31 + NF * CMCF31
*
       ASNF1 = ASNF
       IF (NAORD .EQ. 0) GO TO 1
       ASP   = ASNF
*
       DO 11 K1 = 1, NAORD 
         ASP = ASP * ASNF
         LRHP = 1.D0
*
       DO 12 K2 = 0, K1
         ASNF1 = ASNF1 + ASP * CMC(K1,K2) * LRHP
         LRHP = LRHP * LOGRH
*
  12   CONTINUE
  11   CONTINUE
*
* ---------------------------------------------------------------------
*
   1   RETURN
       END
*
* =================================================================av==
