AUTOMAKE_OPTIONS = foreign
AM_CPPFLAGS = -fpermissive -Wl,Bstatic -I$(srcdir)/../user

if HAVE_ROOT
AM_CPPFLAGS +=  $(ROOT_CFLAGS)
endif

noinst_LTLIBRARIES = libhisto.la
libhisto_la_SOURCES = Kinematics.cxx HistoHandler.cxx HistoObjects.cxx TurboHist.cxx

pkginclude_HEADERS  = Kinematics.h KinematicDefinitions.h \
		      HistoHandler.h HistoObjects.h HistoObjectsROOT.h \
		      HistoObjectsSTL.h HistoBook.h \
		      TurboHist.h
