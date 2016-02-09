ACLOCAL_AMFLAGS=-I m4
#add cuba to SUBDIRS to ship cuba
SUBDIRS = cernlib dequad dyres dynnlo mcfm integrands resum src

#add test directories if enabled
if HAVE_TEST
SUBDIRS += test
endif
if HAVE_DYRES
SUBDIRS += origdyres
endif

#do not ditribute test directories
DIST_SUBDIRS = cernlib dequad dyres dynnlo mcfm integrands resum src

dist_bin_SCRIPTS = install-cuba install-lhapdf
EXTRA_DIST = README.md

bin_SCRIPTS = dyturbo-config
CLEANFILES = $(bin_SCRIPTS)

sharedir = ${datadir}/${PACKAGE}
dist_share_DATA = ${srcdir}/input/*.in ${srcdir}/pdfsets.list ${srcdir}/scripts/submit_DYTURBO.sh ${srcdir}/scripts/DYTURBO_*.in ${srcdir}/scripts/run_DYTURBO_Array_TMPL.sh

#distribute the cuba library
#lib_LIBRARIES = libcuba.a
#libcuba_a_SOURCES =
#libcuba_a_LIBADD  = Cuba-4.2/libcuba.a

distclean-local:
	rm -rf Cuba-$(CUBA_VERSION)
	rm -rf origdyres/*.f origdyres/*.h origdyres/DYRes-v1.0
