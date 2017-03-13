AUTOMAKE_OPTIONS = foreign
#  -I$(srcdir)/../dequad
#	$(CUBA_CPPFLAGS) \
#	-DSHAREDIR=\"${datadir}/${PACKAGE}\"
AM_CPPFLAGS = -fpermissive -Wl,Bstatic -I$(srcdir)/../resum -I$(srcdir)/../src -I$(srcdir)/../mcfm -I$(srcdir)/../dyres -I$(srcdir)/../phasespace $(CUBA_CPPFLAGS) #-I$(srcdir)/../dequad $(LHAPDF_CPPFLAGS)

AM_FFLAGS =  -I$(srcdir)/../dyres -I$(srcdir)/../mcfm -I$(srcdir)/../resum -I$(srcdir)/../src

noinst_LTLIBRARIES = libcounterterm.la
libcounterterm_la_SOURCES = countint.f ctint.C qtint.C

pkginclude_HEADERS = xmio_inc.f largelogs_inc.f ctint.h qtint.h
