ACLOCAL_AMFLAGS=-I m4
#add cuba to SUBDIRS to ship cuba
SUBDIRS = cernlib dequad dyres dynnlo mcfm integrands resum src #origdyres test

#nodist_SUBDIRS = origdyres #test

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

clean-local:
	rm -rf Cuba-$(CUBA_VERSION)
