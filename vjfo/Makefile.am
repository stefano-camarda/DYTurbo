AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = -fpermissive -Wl,Bstatic \
	-I$(srcdir)/../src \
	-I$(srcdir)/../resum \
	-I$(builddir)/../cubature \
	-I$(srcdir)/../mcfm \
	-I$(srcdir)/../phasespace \
	-I$(srcdir)/../resum \
	-I$(srcdir)/../histo \
	-I$(srcdir)/../user \
	$(LHAPDF_CPPFLAGS) $(CUBA_CPPFLAGS)

AM_FFLAGS = -I$(srcdir)/../mcfm

noinst_LTLIBRARIES = libvjfo.la
libvjfo_la_SOURCES = vjint.C vjloint.C sing.C delta.C luminosity.C

libvjfo_la_SOURCES += utils.f abcdefunc.f

pkginclude_HEADERS = vjint.h fodyqt_inc.f functions_inc.f luminosities_inc.f scales2_inc.f vjloint.h luminosity.h dyqcd.f
