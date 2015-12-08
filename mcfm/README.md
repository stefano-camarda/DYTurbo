list of files which have been modified from the "original" mcfm implementation
- dipolesub.f: (speed up by applying kinematic cuts on dipoles counterevents before calculating matrix elements)
- alfamz_lhapdf.f: (direct call to lhapdf alphas mz routine, changed name to dyalphasmz to avoid fortran memory conflicts)
- masses.f: (changed name to dymasses to avoid fortran memory conflicts)
- includedipole.f: apply cuts on each dipole, modified for binning, removed jet clustering, xqtcut

Description of other files in the package
- ddilog.F dilogarithm function from CERNLIB (http://cernlib.web.cern.ch/cernlib/download/2006_source/src/mclibs/isajet/utils/cern_lib/ddilog.F)
- phi1_2m.f modified: switched off jbranch switching -> add a setting to enable or disable jbranch switching
- phase3.f could try to uncomment smin=wsqmin ! for more efficient generation
- ran1.f random number generator (used only by lowinstHst-RES for phi angles generation)
- vegas_common.f: not used, but some routines are still accessing some of the variables, need to cleanup before removing
