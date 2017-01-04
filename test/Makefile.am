AUTOMAKE_OPTIONS = foreign

if HAVE_TEST
AM_CPPFLAGS = -fpermissive \
	      -I$(srcdir)/../src \
	      -I$(srcdir)/../include \
	      -I$(srcdir)/../dequad \
	      -I$(srcdir)/../mcfm \
	      -I$(srcdir)/../dyres \
	      -I$(srcdir)/../resum \
	      -I$(srcdir)/../vjfo \
	      -I$(srcdir)/../counterterm \
	      -I$(srcdir)/../eba \
	      -I$(srcdir)/../cubature \
	      -I$(srcdir)/../phasespace \
	      -I$(srcdir)/../histo \
	      -I$(srcdir)/../user \
	      -Wl,Bstatic \
	      $(CUBA_CPPFLAGS) $(LHAPDF_CPPFLAGS) \
	      -DSHAREDIR=\"${datadir}/${PACKAGE}\"
AM_CPPFLAGS += $(GMOCK_CPPFLAGS)

DYTURBOLIBS= ../src/libdyturbo.la ../cernlib/libcern.la \
	     -lgfortran $(CUBA_LDFLAGS) $(LHAPDF_LDFLAGS)

AM_FFLAGS =  -I$(srcdir)/../dyres -I$(srcdir)/../mcfm -I$(srcdir)/../src

bin_PROGRAMS = points TurboHist_unittest Kinematics_unittest  HistoHandler_unittest Cuts_unittest DYTurbo_unittest

points_SOURCES = points.C finitemapping.C
points_LDADD = ../src/libdyturbo.la \
	       ../cernlib/libcern.la ../dequad/libdequad.la \
	       -lgfortran \
	       $(CUBA_LDFLAGS) $(LHAPDF_LDFLAGS) -lgsl -lgslcblas $(ROOT_LDFLAGS)


# unit tests
TurboHist_unittest_SOURCES = TurboHist_unittest.cxx
TurboHist_unittest_LDADD   = $(GMOCK_LDFLAGS) ../histo/libhisto.la

Kinematics_unittest_SOURCES = Kinematics_unittest.cxx
Kinematics_unittest_LDADD   = $(GMOCK_LDFLAGS) $(DYTURBOLIBS)

HistoHandler_unittest_SOURCES = HistoHandler_unittest.cxx
HistoHandler_unittest_LDADD   = $(GMOCK_LDFLAGS) $(DYTURBOLIBS)

Cuts_unittest_SOURCES = Cuts_unittest.cxx
Cuts_unittest_LDADD   = $(GMOCK_LDFLAGS) $(DYTURBOLIBS)

DYTurbo_unittest_SOURCES = DYTurbo_unittest.cxx
DYTurbo_unittest_LDADD   = $(GMOCK_LDFLAGS) $(DYTURBOLIBS)

if HAVE_ROOT
AM_CPPFLAGS += $(ROOT_CFLAGS)
Kinematics_unittest_LDADD += $(ROOT_LDFLAGS)
TurboHist_unittest_LDADD += $(ROOT_LDFLAGS)
HistoHandler_unittest_LDADD += $(ROOT_LDFLAGS)
DYTurbo_unittest_LDADD += $(ROOT_LDFLAGS)
Cuts_unittest_LDADD   += $(ROOT_LDFLAGS)
endif #ROOT


check:  TurboHist_unittest
	./TurboHist_unittest
	./Kinematics_unittest
	./Cuts_unittest
	./HistoHandler_unittest
	./DYTurbo_unittest
	rm -f *.dat *.root
	@echo ALL OK

endif #TEST

