AUTOMAKE_OPTIONS = foreign
AM_CPPFLAGS = -fpermissive -I$(srcdir)/../include -I$(srcdir)/../dequad $(CUBA_CPPFLAGS) $(LHAPDF_CPPFLAGS) $(ROOT_CFLAGS) \
	-DSHAREDIR=\"${datadir}/${PACKAGE}\"

lib_LTLIBRARIES = libdyturbo.la
libdyturbo_la_SOURCES = resintegr.C integr.C settings.C init.C finintegr.C finitemapping.C cubacall.C plotter.C switch.C
libdyturbo_la_LIBADD = ../dyres/libdyres.la ../dynnlo/libdynnlo.la ../mcfm/libmcfm.la

bin_PROGRAMS = dyturbo dytests merger
dyturbo_SOURCES = dyturbo.C
dyturbo_LDADD = libdyturbo.la ../Cuba-4.2/libcuba.a ../cernlib/libcern.la -lgfortran \
	$(CUBA_LDFLAGS) $(LHAPDF_LDFLAGS) -lgsl -lgslcblas $(ROOT_LDFLAGS) $(ROOT_LIBS) #-lifcore -limf -lifport -lintlc -lquadmath

dytests_SOURCES = dytests.C
dytests_LDADD = libdyturbo.la ../Cuba-4.2/libcuba.a ../cernlib/libcern.la -lgfortran \
	$(CUBA_LDFLAGS) $(LHAPDF_LDFLAGS) -lgsl -lgslcblas $(ROOT_LDFLAGS) $(ROOT_LIBS) #-lifcore -limf -lifport -lintlc -lquadmath

dytests_SOURCES = merger.C
dytests_LDADD = $(ROOT_LDFLAGS) $(ROOT_LIBS)

pkginclude_HEADERS = resintegr.h finintegr.h integr.h interface.h settings.h cubacall.h plotter.h finitemapping.h lines.C init.h switch.h options.f
