ACLOCAL_AMFLAGS=-I m4
#do not ditribute test directories
# add cuba to SUBDIRS to ship cuba
DIST_SUBDIRS = histo user cernlib ancont pegasus cubature dequad sandia smolpack dyres dynnlo mcfm integrands resum counterterm v2j eba phasespace vjfo src
SUBDIRS      = ${DIST_SUBDIRS}

# add test directories if enabled
if HAVE_TEST
SUBDIRS += test
endif
if HAVE_DYRES
SUBDIRS += origdyres
endif


dist_bin_SCRIPTS = install-cuba install-lhapdf install-chaplin
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

distclean-local:
	rm -rf Cuba-$(CUBA_VERSION)
	rm -rf lhapdf6
	rm -rf chaplin-$(CHAPLIN_VERSION)
	rm -rf origdyres/*.f origdyres/*.h origdyres/DYRes-v1.0
