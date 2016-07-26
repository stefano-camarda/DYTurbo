AUTOMAKE_OPTIONS = foreign
AM_CPPFLAGS = -fpermissive -Wl,Bstatic #-I$(srcdir)/../src -I$(srcdir)/../dequad -I$(srcdir)/../cubature -I$(srcdir)/../mcfm $(LHAPDF_CPPFLAGS)

AM_FFLAGS = -I$(srcdir)/../mcfm -I$(srcdir)/../src

noinst_LTLIBRARIES = libv2j.la
libv2j_la_SOURCES = v2jint.f

#pkginclude_HEADERS = v2jint.h
