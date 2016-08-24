AUTOMAKE_OPTIONS = foreign
AM_CPPFLAGS = -fpermissive -Wl,Bstatic -I$(srcdir)/../resum -I$(srcdir)/../src -I$(srcdir)/../mcfm -I$(srcdir)/../phasespace $(CUBA_CPPFLAGS) #$(LHAPDF_CPPFLAGS)

noinst_LTLIBRARIES = libeba.la
libeba_la_SOURCES = loint.C

pkginclude_HEADERS = loint.h
