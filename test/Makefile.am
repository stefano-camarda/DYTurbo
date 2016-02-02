AUTOMAKE_OPTIONS = foreign

if HAVE_TEST
AM_CPPFLAGS = -fpermissive -I$(srcdir)/../include -I$(srcdir)/../dequad -I$(srcdir)/../resum -I$(srcdir)/../src \
	-Wl,Bstatic \
	$(CUBA_CPPFLAGS) $(LHAPDF_CPPFLAGS) \
	-DSHAREDIR=\"${datadir}/${PACKAGE}\"

AM_FFLAGS =  -I$(srcdir)/../dyres -I$(srcdir)/../mcfm -I$(srcdir)/../src

AM_CPPFLAGS += $(ROOT_CFLAGS)

bin_PROGRAMS = dytests

dytests_SOURCES = dytests.C
dytests_LDADD = ../src/libdyturbo.la \
	../cernlib/libcern.la ../dequad/libdequad.la \
	-lgfortran \
	$(CUBA_LDFLAGS) $(LHAPDF_LDFLAGS) -lgsl -lgslcblas $(ROOT_LDFLAGS)
endif
