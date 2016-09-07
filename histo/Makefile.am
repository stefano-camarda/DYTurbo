AUTOMAKE_OPTIONS = foreign
AM_CPPFLAGS = -fpermissive -Wl,Bstatic -I$(srcdir)/../user

if HAVE_ROOT
AM_CPPFLAGS +=  $(ROOT_CFLAGS)
endif


noinst_LTLIBRARIES = libkin.la libhist.la libturbohist.la
libkin_la_SOURCES = Kinematics.cxx
libhist_la_SOURCES = HistoHandler.cxx HistoObjects.cxx
libturbohist_la_SOURCES = TurboHist.cxx


pkginclude_HEADERS  = Kinematics.h KinematicDefinitions.h
pkginclude_HEADERS += HistoHandler.h HistoObjects.h HistoObjectsROOT.h HistoObjectsSTL.h HistoBook.h
pkginclude_HEADERS += TurboHist.h
