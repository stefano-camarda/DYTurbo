list of files which have been modified from the "original" mcfm implementation
- dipolesub.f (speed up by applying kinematic cuts on dipoles counterevents before calculating matrix elements)
- alfamz_lhapdf.f (direct call to lhapdf alphas mz routine, changed name to dyalphasmz to avoid fortran memory conflicts)
- ddilog.F dilogarithm function from CERNLIB (http://cernlib.web.cern.ch/cernlib/download/2006_source/src/mclibs/isajet/utils/cern_lib/ddilog.F)
- phi1_2m.f modified: switched off jbranch switching -> add a setting to enable or disable jbranch switching
- phase3.f could try to uncomment smin=wsqmin ! for more efficient generation
- ran1.f random number generator (used only by lowinstHst-RES for phi angles generation)
