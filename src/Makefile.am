AM_CPPFLAGS = -fpermissive -I$(srcdir)/../include -I$(srcdir)/../Cuba-4.2 -I$(srcdir)/../dequad $(LHAPDF_CPPFLAGS)
bin_PROGRAMS = dyturbo
dyturbo_SOURCES = dyturbo.C integr.C settings.C init.C finintegr.C finitemapping.C cubacall.C
dyturbo_LDADD = ../dyres/libdyres.a ../lib/libcern.a -lgfortran $(CUBA_LDFLAGS) $(LHAPDF_LDFLAGS) -lgsl -lgslcblas #-lifcore -limf -lifport -lintlc -lquadmath

pkginclude_HEADERS = finintegr.h  integr.h  interface.h  settings.h cubacall.h
