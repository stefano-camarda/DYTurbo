AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = -fpermissive -Wl,Bstatic \
                -I$(srcdir)/../specialfunctions \
                -I$(srcdir)/../src \
                -I$(srcdir)/../mcfm \
                -I$(srcdir)/../resum

noinst_LTLIBRARIES = libbessel.la
libbessel_la_SOURCES = hankel.C 588.f

pkginclude_HEADERS = hankel.h 588_interface.h

