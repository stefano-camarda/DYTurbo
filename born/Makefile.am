AUTOMAKE_OPTIONS = foreign
AM_CPPFLAGS = -fpermissive -Wl,Bstatic -I$(srcdir)/../resum -I$(srcdir)/../src -I$(srcdir)/../dyres -I$(srcdir)/../dynnlo -I$(srcdir)/../mcfm -I$(srcdir)/../phasespace $(CUBA_CPPFLAGS) $(LHAPDF_CPPFLAGS)

noinst_LTLIBRARIES = libborn.la
libborn_la_SOURCES = loint.C lomellin.C

pkginclude_HEADERS = loint.h lomellin.h
