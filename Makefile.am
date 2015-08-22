ACLOCAL_AMFLAGS=-I m4
SUBDIRS = cernlib dyres src
EXTRA_DIST= install-cuba install-lhapdf pdfsets.list README.md
bin_SCRIPTS = dyturbo-config

sharedir = ${datadir}/${PACKAGE}
dist_share_DATA = input/*.in

#distribute the cuba library
lib_LIBRARIES = libcuba.a
libcuba_a_SOURCES =
libcuba_a_LIBADD  = Cuba-4.2/libcuba.a

pkginclude_HEADERS = Cuba-4.2/cuba.h
