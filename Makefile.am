ACLOCAL_AMFLAGS=-I m4
#do not ditribute test directories
# add cuba to SUBDIRS to ship cuba
DIST_SUBDIRS = histo user cernlib ancont pegasus cubature dequad bequad sandia smolpack besche bessel minprescription specialfunctions dyres dynnlo mcfm integrands resum counterterm v2j born phasespace vjet src
SUBDIRS      = ${DIST_SUBDIRS}

# add test directories if enabled
if HAVE_TEST
SUBDIRS += test
endif
if HAVE_DYRES
SUBDIRS += origdyres
endif


dist_bin_SCRIPTS = install-cuba install-lhapdf install-chaplin install-cubature
EXTRA_DIST = README.md

bin_SCRIPTS = dyturbo-config
CLEANFILES = $(bin_SCRIPTS)

sharedir = ${datadir}/${PACKAGE}
dist_share_DATA = ${srcdir}/input/*.in \
		  ${srcdir}/pdfsets.list \
		  ${srcdir}/scripts/compile_grid.sh \
		  ${srcdir}/scripts/run_grid.sh \
		  ${srcdir}/scripts/run_on_batch_tmpl.sh \
		  ${srcdir}/scripts/run_prun.sh \
		  ${srcdir}/scripts/submit_jobs.sh

#distribute the cuba library
#lib_LIBRARIES = libcuba.a
#libcuba_a_SOURCES =
#libcuba_a_LIBADD  = Cuba-4.2/libcuba.a

installcheck-local:
	${bindir}/dyturbo ${sharedir}/test.in
	echo 's_qt->Integral()' | root -l -b results.root
	rm -rf ./tmp_results_main.dat
	rm -rf ./results.txt

distclean-local:
	rm -rf Cuba-$(CUBA_VERSION)
	rm -rf lhapdf6
	rm -rf chaplin-$(CHAPLIN_VERSION)
	rm -rf cubature/*.c cubature/*.h cubature/NEWS cubature/COPYING cubature/README cubature/ChangeLog
	rm -rf origdyres/*.f origdyres/*.h origdyres/DYRes-v1.0
