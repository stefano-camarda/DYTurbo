AUTOMAKE_OPTIONS = foreign
AM_CPPFLAGS = -fpermissive -I$(srcdir)/../include -I$(srcdir)/../dequad -I$(srcdir)/../resum \
	-Wl,Bstatic \
	$(CUBA_CPPFLAGS) $(LHAPDF_CPPFLAGS) \
	-DSHAREDIR=\"${datadir}/${PACKAGE}\"

AM_FFLAGS =  -I$(srcdir)/../dyres -I$(srcdir)/../mcfm #-I$(srcdir)/../src

noinst_LTLIBRARIES = libdyturbobase.la
libdyturbobase_la_SOURCES = resintegr.C integr.C settings.C init.C finintegr.C finitemapping.C cubacall.C plotter.C switch.C cuts.C usercuts.C pdf.C printsettings.C coupling.C
libdyturbobase_la_SOURCES += flavreduce.f gaussinit.f rescinit.f pdfmoments.f

lib_LTLIBRARIES = libdyturbo.la
libdyturbo_la_SOURCES =
libdyturbo_la_LIBADD =  libdyturbobase.la ../integrands/libintegrands.la ../dyres/libdyres.la ../dynnlo/libdynnlo.la ../mcfm/libmcfm.la ../resum/libresum.la

bin_PROGRAMS = dyturbo
dyturbo_SOURCES = dyturbo.C
dyturbo_LDADD = libdyturbo.la \
	../cernlib/libcern.la ../dequad/libdequad.la \
	-lgfortran \
	$(CUBA_LDFLAGS) $(LHAPDF_LDFLAGS) -lgsl -lgslcblas #-lifcore -limf -lifport -lintlc -lquadmath

if HAVE_ROOT
dyturbo_LDADD += $(ROOT_LDFLAGS)
AM_CPPFLAGS += $(ROOT_CFLAGS)

bin_PROGRAMS += merger

merger_SOURCES = merger.C
merger_LDADD = $(ROOT_LDFLAGS)
endif

pkginclude_HEADERS = resintegr.h finintegr.h integr.h interface.h settings.h cubacall.h plotter.h finitemapping.h init.h switch.h options.f flred.f gauss.f quadrules.f cuts.h coupling.h pdf.h printsettings.h
