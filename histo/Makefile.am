AUTOMAKE_OPTIONS = foreign
AM_CPPFLAGS = -fpermissive -Wl,Bstatic -I$(srcdir)/../user

if HAVE_ROOT
AM_CPPFLAGS +=  $(ROOT_CFLAGS)
endif

noinst_LTLIBRARIES = libhisto.la
libhisto_la_SOURCES = Kinematics.cxx HistoHandler.cxx HistoObjects.cxx TurboHist.cxx

pkginclude_HEADERS  = Kinematics.h KinematicDefinitions.h KinematicCuts.h KinUtils.h\
		      HistoHandler.h HistoObjects.h HistoBase.h HistoSpecialization.h \
		      HistoObjectsROOT.h HistoObjectsSTL.h HistoBook.h \
		      TurboHist.h TurboHist_Binning.h TurboHist_Counter.h TurboHist_File.h \
		      TurboHist_HBase.h TurboHist_H1.h TurboHist_H2.h TurboHist_H3.h TurboHist_P1.h TurboHist_P2.h
