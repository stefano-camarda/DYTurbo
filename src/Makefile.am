AM_CPPFLAGS = -ggdb -fpermissive -I$(srcdir)/../include -I$(srcdir)/../Cuba-4.2 -I$(srcdir)/../dequad $(LHAPDF_CPPFLAGS) $(ROOT_CFLAGS)
BASE_SRC =  integr.C settings.C init.C plotter.C
bin_PROGRAMS = dyfast dyfin
dyfast_SOURCES = dyfast.C $(BASE_SRC)
dyfast_LDADD = ../dyres/libdyres.a ../lib/libcern.a -lgfortran $(CUBA_LDFLAGS) $(LHAPDF_LDFLAGS) -lgsl -lgslcblas $(ROOT_LDFLAGS) $(ROOT_LIBS) #-lifcore -limf -lifport -lintlc -lquadmath

dyfin_SOURCES = dyfin.C dyct.C dyrealvirt.C $(BASE_SRC)
dyfin_LDADD = ../dyres/libdyres.a ../lib/libcern.a -lgfortran -L$(srcdir)/../Cuba-4.2 -lcuba $(LHAPDF_LDFLAGS) -lgsl -lgslcblas $(ROOT_LDFLAGS) $(ROOT_LIBS)

pkginclude_HEADERS = dyct.h  dyfast.h  dyfin.h  dyrealvirt.h  integr.h  interface.h  settings.h plotter.h
