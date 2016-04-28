AUTOMAKE_OPTIONS = foreign
#  -I$(srcdir)/../dequad
#	$(CUBA_CPPFLAGS) \
#	-DSHAREDIR=\"${datadir}/${PACKAGE}\"
AM_CPPFLAGS = -fpermissive -Wl,Bstatic #-I$(srcdir)/../src -I$(srcdir)/../dequad $(LHAPDF_CPPFLAGS)

AM_FFLAGS =  -I$(srcdir)/../dyres -I$(srcdir)/../mcfm -I$(srcdir)/../resum -I$(srcdir)/../src

noinst_LTLIBRARIES = libcounterterm.la
libcounterterm_la_SOURCES = ctint.f

pkginclude_HEADERS = xmio_inc.f largelogs_inc.f
