AUTOMAKE_OPTIONS = foreign

if HAVE_TEST
AM_CPPFLAGS = -fpermissive -I$(srcdir)/../include -I$(srcdir)/../histo -I$(srcdir)/../user -I$(srcdir)/../dequad -I$(srcdir)/../resum -I$(srcdir)/../src -I$(srcdir)/../mcfm -I$(srcdir)/../phasespace \
	      -Wl,Bstatic \
	      $(CUBA_CPPFLAGS) $(LHAPDF_CPPFLAGS) \
	      -DSHAREDIR=\"${datadir}/${PACKAGE}\"
AM_CPPFLAGS += $(GMOCK_CPPFLAGS)

AM_FFLAGS =  -I$(srcdir)/../dyres -I$(srcdir)/../mcfm -I$(srcdir)/../src


bin_PROGRAMS = points TurboHist_unittest Kinematics_unittest HistoHandlerSTL_unittest HistoHandlerROOT_unittest Cuts_unittest

points_SOURCES = points.C finitemapping.C
points_LDADD = ../src/libdyturbo.la \
	       ../cernlib/libcern.la ../dequad/libdequad.la \
	       -lgfortran \
	       $(CUBA_LDFLAGS) $(LHAPDF_LDFLAGS) -lgsl -lgslcblas $(ROOT_LDFLAGS)

# unit tests
TurboHist_unittest_SOURCES = TurboHist_unittest.cxx
TurboHist_unittest_LDADD   = $(GMOCK_LDFLAGS) ../histo/libturbohist.la

Kinematics_unittest_SOURCES = Kinematics_unittest.cxx
Kinematics_unittest_LDADD   = $(GMOCK_LDFLAGS) ../histo/libkin.la

HistoHandlerSTL_unittest_SOURCES = HistoHandlerSTL_unittest.cxx
HistoHandlerSTL_unittest_LDADD   = $(GMOCK_LDFLAGS) ../histo/libhist.la

HistoHandlerROOT_unittest_SOURCES = HistoHandlerROOT_unittest.cxx
HistoHandlerROOT_unittest_LDADD   = $(GMOCK_LDFLAGS) ../histo/libhist.la

Cuts_unittest_SOURCES = Cuts_unittest.cxx
Cuts_unittest_LDADD   = $(GMOCK_LDFLAGS) ../histo/libkin.la




if HAVE_ROOT
AM_CPPFLAGS += $(ROOT_CFLAGS)
TurboHist_unittest_LDADD += $(ROOT_LDFLAGS)
HistoHandlerROOT_unittest_LDADD += $(ROOT_LDFLAGS)

endif

check:  TurboHist_unittest
	./TurboHist_unittest
	./Kinematics_unittest
	./HistoHandlerSTL_unittest
	./HistoHandlerROOT_unittest
	rm -f *.dat *.root

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

