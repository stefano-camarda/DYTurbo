*
* ..File: pnsg0mom.f
*
*
* ..The LO (alpha_s^1) non-singlet (subroutine  PNS0MOM)  and singlet 
*    (subroutine  PSG0MOM)  QCD splitting function  P0NS  and  P0SG  
*    in N-space for  NF = 3...6  (parameters NFMIN and NFMAX) massless 
*    quark flavours.
*    The coupling constant is normalized as  a_s = alpha_s/(4*pi). 
*
* ..These quantities are determined on an external NDIM-dimensional 
*    array  NA  of complex Mellin moments provided by the common-block 
*    MOMS.  The non-singlet and singlet results are written to the 
*    common-blocks  PNS0  and  PSG0, respectively. The last two singlet
*    array arguments refer to the quark-gluon mixing matrix with  1 = q 
*    and  2 = g.
*
* ..The SU(N_colours=3) colour factors  CF, CA and TF  are taken from 
*    the common-block  COLOUR.  The harmonic sums S_1(N) are provided 
*    by the common-block  HSUMS.
*
* =====================================================================
*
*
       SUBROUTINE PNS0MOM 
*
       IMPLICIT DOUBLE COMPLEX (A - Z)
       INTEGER NMAX, NDIM, NFMIN, NFMAX, KN, NF
       PARAMETER (NDIM = 144, NFMIN = 3, NFMAX = 6)
       DOUBLE PRECISION CF, CA, TR
*
* ---------------------------------------------------------------------
*
* ..Input common-blocks 
*
       COMMON / MOMS   / NA (NDIM)
       COMMON / NNUSED / NMAX
       COMMON / HSUMS  / S(NDIM,6)
       COMMON / COLOUR / CF, CA, TR
*
* ..Output common-block 
*
       COMMON / PNS0   / P0NS (NDIM, NFMIN:NFMAX)
*
* ---------------------------------------------------------------------
*
* ..Mellin-N loop
*
       DO 1 KN = 1, NMAX
*
       N  = NA(KN)
       PNSA = 3. - 4.* S(KN,1) + 2./(N * (N+1.))
*
* ..Output to the array
*
       DO 2 NF = NFMIN, NFMAX
*
       P0NS(KN,NF) = CF * PNSA
*
* ---------------------------------------------------------------------
*
  2    CONTINUE
  1    CONTINUE
*
       RETURN
       END
*
* =====================================================================
*
*
       SUBROUTINE PSG0MOM 
*
       IMPLICIT DOUBLE COMPLEX (A - Z)
       INTEGER NMAX, NDIM, NFMIN, NFMAX, KN, NF
       PARAMETER (NDIM = 144, NFMIN = 3, NFMAX = 6)
       DOUBLE PRECISION CF, CA, TR
*
* ---------------------------------------------------------------------
*
* ..Input common-blocks  
*
       COMMON / MOMS   / NA (NDIM)
       COMMON / NNUSED / NMAX
       COMMON / HSUMS  / S(NDIM,6)
       COMMON / COLOUR / CF, CA, TR
*
* ..Output common-block  
*
       COMMON / PSG0   / P0SG (NDIM, NFMIN:NFMAX, 2, 2)
*
* ---------------------------------------------------------------------
*
* ..Mellin-N loop
*
       DO 1 KN = 1, NMAX
*
       N  = NA(KN)
       NS = N * N
       N1 = N + 1.
       N2 = N + 2.
       NM = N - 1.
*
       PQQA = 3. - 4.* S(KN,1) + 2./(N * N1)
       PQGA = 4.* (NS + N + 2.) / (N * N1 * N2)
       PGQA = 2.* (NS + N + 2.) / (N * N1 * NM)
       PGGA = 11./3.D0 - 4.* S(KN,1) + 4./(N * NM) + 4./(N1 * N2) 
       PGGB = - 4./3.D0
*
* ..Flavour-number loop and output to the array
*
       DO 2 NF = NFMIN, NFMAX
*
       P0SG(KN,NF,1,1) = CF * PQQA
       P0SG(KN,NF,1,2) = TR * NF * PQGA
       P0SG(KN,NF,2,1) = CF * PGQA
       P0SG(KN,NF,2,2) = CA * PGGA + TR * NF * PGGB
*
* ---------------------------------------------------------------------
*
  2    CONTINUE
  1    CONTINUE
*
       RETURN
       END
*
* =================================================================av==
