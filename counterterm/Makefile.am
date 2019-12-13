AUTOMAKE_OPTIONS = foreign
#	$(CUBA_CPPFLAGS) \
#	-DSHAREDIR=\"${datadir}/${PACKAGE}\"
AM_CPPFLAGS = -fpermissive -Wl,Bstatic \
	-I$(srcdir)/../resum \
	-I$(srcdir)/../src \
	-I$(srcdir)/../mcfm \
	-I$(srcdir)/../dyres \
	-I$(srcdir)/../phasespace \
	-I$(srcdir)/../dequad \
	-I$(srcdir)/../bessel \
	$(CUBA_CPPFLAGS) $(LHAPDF_CPPFLAGS)

AM_FFLAGS =  -I$(srcdir)/../dyres -I$(srcdir)/../mcfm -I$(srcdir)/../resum -I$(srcdir)/../src

noinst_LTLIBRARIES = libcounterterm.la
libcounterterm_la_SOURCES = countint.f ctint.C qtint.C itilde.C ctmellin.C

pkginclude_HEADERS = xmio_inc.f largelogs_inc.f ctint.h qtint.h itilde.h ctmellin.h
