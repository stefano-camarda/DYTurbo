AUTOMAKE_OPTIONS = foreign
AM_CPPFLAGS = -fpermissive -Wl,Bstatic -I$(srcdir)/../user


noinst_LTLIBRARIES = libkin.la libhist.la libturbohist.la
libkin_la_SOURCES = Kinematics.cxx
libhist_la_SOURCES = HistoHadler.cxx
libturbohist_la_SOURCES = TurboHist.cxx

pkginclude_HEADERS = Kinematics.h HistoHadler.h TurboHist.h
