* integration.f adaptative integration code, probably from John Collins and Wu Ki Tung (http://users.phys.psu.edu/~collins/software/)
It is called from BK(n,z) (derivative of bessel function), which is called from IK(n), which is called by evaluatekin, but commented, so actually not used
* vegas_common.f not used, need to check if any variable is still accessed before removing
* virtonly.f not used, need to check if any variable is still accessed before removing
* realonly.f not used, need to check if any variable is still accessed before removing
* mellinh2.f Bluemlein functions PSI0(ZZ,RES), PSI1(ZZ,RES), PSI2(ZZ,RES), PSI3(ZZ,RES), BET(ZZ,RES), BET1(ZZ,RES), BET2(ZZ,RES), BET3(ZZ,RES) (from https://www-zeuthen.desy.de/~blumlein/CODE/ANCONT/ancont.f)
