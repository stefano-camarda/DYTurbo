* integration.f: adaptative integration code, probably from John Collins and Wu Ki Tung (http://users.phys.psu.edu/~collins/software/)
It is called from BK(n,z) (derivative of bessel function), which is called from IK(n), which is called by evaluatekin, but commented, so actually not used
* mellinh2.f: Bluemlein functions PSI0(ZZ,RES), PSI1(ZZ,RES), PSI2(ZZ,RES), PSI3(ZZ,RES), BET(ZZ,RES), BET1(ZZ,RES), BET2(ZZ,RES), BET3(ZZ,RES) (from https://www-zeuthen.desy.de/~blumlein/CODE/ANCONT/ancont.f Comput.Phys.Commun. 133 (2000) 76-104)
* resinit.f: initialisation for resummed integrand (integrands/main2.f). Part of the initialisation should be moved to src/rescinit.f
* besselkfast.f: Bessel functions. !!! LOCATE and POLINT functions are from "Numerical Recipes in Fortran 77" Be careful, this code has a copyright !!! (-> locate renamed to dylocate to avoid conflicts with QCDNUM, polint renamed to xpolint probably to avoid conflict with CTEQ PDFs)
* const.h: DYRES constants (some repetition with MCFM constants, to clean up)
* countDYnew-RES.f: counterterm integrand
* lowintHst-RES.f: resummed integrand
* csection.f: ratio of original/approximated PDF (used only for the original integration)
* gen4MIO.f: 4 body phase space generation, slightly modified from MCFM gen4.f
* genBORN2qT.f: phase space generation for the counterterm, used only for the original counterterm
* myli2.f: DYRES version of ddilog.f from CERNLIB
* rescoeff.f: resummation coefficients
* scales.h: DYRES common variables for scales
* varios.f: PDF mellin moments fit
* anom.f: anomalous dimensions and wilson coefficients in Mellin space (extracted from main2.f)
