AM_CPPFLAGS = -fpermissive -I$(srcdir)/../include -I$(srcdir)/../Cuba-4.2 -I$(srcdir)/../dequad $(LHAPDF_CPPFLAGS)
bin_PROGRAMS = dyfast dyfin
dyfast_SOURCES = dyfast.C integr.C settings.C init.C
dyfast_LDADD = ../dyres/libdyres.a ../lib/libcern.a -lgfortran $(CUBA_LDFLAGS) $(LHAPDF_LDFLAGS) -lgsl -lgslcblas #-lifcore -limf -lifport -lintlc -lquadmath

dyfin_SOURCES = dyfin.C dyct.C finitemapping.C integr.C settings.C init.C
dyfin_LDADD = ../dyres/libdyres.a ../lib/libcern.a -lgfortran -L$(srcdir)/../Cuba-4.2 -lcuba $(LHAPDF_LDFLAGS) -lgsl -lgslcblas

pkginclude_HEADERS = dyct.h  integr.h  interface.h  settings.h
