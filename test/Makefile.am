AUTOMAKE_OPTIONS = foreign

if HAVE_TEST
AM_CPPFLAGS = -fpermissive -I$(srcdir)/../include -I$(srcdir)/../histo -I$(srcdir)/../user -I$(srcdir)/../dequad -I$(srcdir)/../resum -I$(srcdir)/../src -I$(srcdir)/../mcfm -I$(srcdir)/../phasespace \
	      -Wl,Bstatic \
	      $(CUBA_CPPFLAGS) $(LHAPDF_CPPFLAGS) \
	      -DSHAREDIR=\"${datadir}/${PACKAGE}\"
AM_CPPFLAGS += $(GMOCK_CPPFLAGS)

AM_FFLAGS =  -I$(srcdir)/../dyres -I$(srcdir)/../mcfm -I$(srcdir)/../src

if HAVE_ROOT
AM_CPPFLAGS += $(ROOT_CFLAGS)
endif

bin_PROGRAMS = points TurboHist_unittest

points_SOURCES = points.C finitemapping.C
points_LDADD = ../src/libdyturbo.la \
	       ../cernlib/libcern.la ../dequad/libdequad.la \
	       -lgfortran \
	       $(CUBA_LDFLAGS) $(LHAPDF_LDFLAGS) -lgsl -lgslcblas $(ROOT_LDFLAGS)

TurboHist_unittest_SOURCES = TurboHist_unittest.cxx
TurboHist_unittest_LDADD   = $(GMOCK_LDFLAGS) ../histo/libturbohist.la


test:  TurboHist_unittest
	./TurboHist_unittest

endif

#
#bin_PROGRAMS += dytests
#
#dytests_SOURCES = dytests.C
#dytests_LDADD = ../src/libdyturbo.la \
#	../cernlib/libcern.la ../dequad/libdequad.la \
#	-lgfortran \
#	$(CUBA_LDFLAGS) $(LHAPDF_LDFLAGS) -lgsl -lgslcblas $(ROOT_LDFLAGS)
#
#endif

