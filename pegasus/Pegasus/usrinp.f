* ..File: usrinp.f  
*
*
* ..Input values for initial scale and alpha_s, the heavy-quark masses 
*    and the parton-distribution shapes used by  INITINP (IPAR)  if 
*    called with  IPAR = 2.  See section 4.2 of the manual.
*                   
* =====================================================================
*
*
       SUBROUTINE USRINP (PUV, PDV, PLM, PLP, PSM, PSP, PGL, M20, ASI,
     1                    MC2, MB2, MT2, NFORM, IMOMIN, ISSIMP)
* 
       IMPLICIT DOUBLE PRECISION (A - Z)
       INTEGER IMOMIN, ISSIMP, NFORM 
       DIMENSION PUV(6), PDV(6), PLM(6), PLP(6), PSM(6), PSP(6), PGL(6)
*
       M20 = 2.D0
       ASI = 0.35000
*
* ..The otherwise irrelevant shift of m_c ensures that the input is
*   returned for m_c = mu_0 also in the VFNS. Useful for input checks.
*
       MC2 = 2.00001D0
       MB2 = 20.25D0
       MT2 = 3.0625D4
*
       NFORM  = 1
       IMOMIN = 1
       ISSIMP = 0
*
       PUV(1) =  2.0D0
       PUV(2) =  0.8D0
       PUV(3) =  3.0D0
       PUV(4) =  0.5D0
       PUV(5) =  0.0D0
       PUV(6) =  0.0D0
*
       PDV(1) =  1.0D0
       PDV(2) =  0.8D0
       PDV(3) =  4.0D0
       PDV(4) =  0.5D0
       PDV(5) =  0.0D0
       PDV(6) =  0.0D0
*
       PLM(1) =  0.193987D0
       PLM(2) =  0.9D0
       PLM(3) =  6.0D0
       PLM(4) =  0.5D0
       PLM(5) =  0.0D0
       PLM(6) =  0.0D0
*           
       PLP(1) =  0.136565D0
       PLP(2) = -0.1D0
       PLP(3) =  6.0D0
       PLP(4) =  0.5D0
       PLP(5) =  0.0D0
       PLP(6) = -0.5D0
*          
       PSM(1) =  0.0d0
       PSM(2) =  0.9D0
       PSM(3) =  6.0D0
       PSM(4) =  0.5D0
       PSM(5) =  0.0D0
       PSM(6) =  0.0D0
*          
       PSP(1) =  0.027313d0
       PSP(2) = -0.1D0
       PSP(3) =  6.0D0
       PSP(4) =  0.5D0
       PSP(5) =  0.0D0
       PSP(6) = -0.5D0
*
       PGL(1) =  1.0D0
       PGL(2) = -0.1D0
       PGL(3) =  5.0D0
       PGL(4) =  0.5D0
       PGL(5) =  0.0D0
       PGL(6) =  0.0D0
*
       RETURN
       END
*
* =================================================================av==
