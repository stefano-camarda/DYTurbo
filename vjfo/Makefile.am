AUTOMAKE_OPTIONS = foreign
#  -I$(srcdir)/../dequad
#	$(CUBA_CPPFLAGS) \
#	-DSHAREDIR=\"${datadir}/${PACKAGE}\"
AM_CPPFLAGS = -fpermissive -Wl,Bstatic -I$(srcdir)/../src -I$(srcdir)/../resum $(LHAPDF_CPPFLAGS)

AM_FFLAGS = -I$(srcdir)/../mcfm #-I$(srcdir)/../dyres 

noinst_LTLIBRARIES = libvjfo.la
libvjfo_la_SOURCES = vjfo.f vjint.C

pkginclude_HEADERS = vjint.h
