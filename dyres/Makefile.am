AUTOMAKE_OPTIONS = foreign

AM_FFLAGS = -I$(srcdir)/../mcfm -I$(srcdir)/../src -I$(srcdir)/../resum -I$(srcdir)/../counterterm
AM_CPPFLAGS = -I$(srcdir)/../mcfm -I$(srcdir)/../src

noinst_LTLIBRARIES = libdyres.la
libdyres_la_SOURCES = varios.f csection.f integration.f resinit.f besselkfast.f countDYnew-RES.f gen4MIO.f genBORN2qT.f anom.f dyres_interface.C

pkginclude_HEADERS = const.h rescoeff.f scales.h dyres_interface.h
