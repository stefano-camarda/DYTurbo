AUTOMAKE_OPTIONS = foreign
#  -I$(srcdir)/../dequad
#	$(CUBA_CPPFLAGS) \
#	-DSHAREDIR=\"${datadir}/${PACKAGE}\"
AM_CPPFLAGS = -fpermissive -Wl,Bstatic -I$(srcdir)/../src -I$(srcdir)/../resum -I$(srcdir)/../cubature -I$(srcdir)/../mcfm -I$(srcdir)/../phasespace -I$(srcdir)/../resum -I$(srcdir)/../histo -I$(srcdir)/../user $(LHAPDF_CPPFLAGS) $(CUBA_CPPFLAGS)

AM_FFLAGS = -I$(srcdir)/../mcfm #-I$(srcdir)/../dyres 

noinst_LTLIBRARIES = libvjfo.la
libvjfo_la_SOURCES = vjfo.f vjint.C vjloint.C sing.C

pkginclude_HEADERS = vjint.h fodyqt_inc.f functions_inc.f internal_inc.f luminosities_inc.f scales2_inc.f vjloint.h
