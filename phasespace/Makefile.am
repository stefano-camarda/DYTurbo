AUTOMAKE_OPTIONS = foreign
AM_CPPFLAGS = -fpermissive -Wl,Bstatic \
		-I$(srcdir)/../mcfm \
		-I$(srcdir)/../src #-I$(srcdir)/../resum $(CUBA_CPPFLAGS) #$(LHAPDF_CPPFLAGS)

AM_FFLAGS = -I$(srcdir)/../mcfm

noinst_LTLIBRARIES = libphasespace.la
libphasespace_la_SOURCES = psset.C psgen.C psbin.C gen4p.C

libphasespace_la_SOURCES += dyboost.f boostV.f qtweight.f mweight.f genp.f dyspinoru.f rotate.f dygen4.f dygen3.f

pkginclude_HEADERS = phasespace.h
