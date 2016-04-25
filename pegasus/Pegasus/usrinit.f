*
* ..File: usrinit.f
*
*
* ..Initialization parameters for speed, the heavy-flavour treatment,
*    the mode and order of the evolution, and the value of  FR2 = 
*    mu_f^2/mu_r^2  used by INITEVOL (IPAR) if called with  IPAR = 2.
*    See section 4.1 of the manual.
*
* =====================================================================
*
*
       SUBROUTINE USRINIT (IFAST, IVFNS, NFF, IMODEV, NPORD, FR2)
* 
       IMPLICIT DOUBLE PRECISION (A - Z)
       INTEGER IFAST, IVFNS, NFF, IMODEV, NPORD
*
       IFAST  = 0
       IVFNS  = 1
       NFF    = 4
       IMODEV = 1
       NPORD  = 1
       FR2    = 1.D0
*
       RETURN
       END
*
* =================================================================av==
