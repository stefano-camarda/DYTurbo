AUTOMAKE_OPTIONS = foreign
AM_CPPFLAGS = -fpermissive -I$(srcdir)/../include -I$(srcdir)/../Cuba-4.2 -I$(srcdir)/../dequad $(LHAPDF_CPPFLAGS) $(ROOT_CFLAGS)

lib_LTLIBRARIES = libdyturbo.la
libdyturbo_la_SOURCES = integr.C settings.C init.C finintegr.C finitemapping.C cubacall.C plotter.C

bin_PROGRAMS = dyturbo
dyturbo_SOURCES = dyturbo.C #integr.C settings.C init.C finintegr.C finitemapping.C cubacall.C plotter.C
dyturbo_LDADD = libdyturbo.la ../dyres/libdyres.la ../cernlib/libcern.la -lgfortran \
	$(CUBA_LDFLAGS) $(LHAPDF_LDFLAGS) -lgsl -lgslcblas $(ROOT_LDFLAGS) $(ROOT_LIBS) #-lifcore -limf -lifport -lintlc -lquadmath

pkginclude_HEADERS = finintegr.h integr.h interface.h settings.h cubacall.h plotter.h finitemapping.h lines.C
