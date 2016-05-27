AUTOMAKE_OPTIONS = foreign
AM_CPPFLAGS = -fpermissive \
	-I$(srcdir)/../include \
	-I$(srcdir)/../dequad \
	-I$(srcdir)/../mcfm \
	-I$(srcdir)/../resum \
	-I$(srcdir)/../vjfo \
	-I$(srcdir)/../cubature \
	-Wl,Bstatic \
	$(CUBA_CPPFLAGS) $(LHAPDF_CPPFLAGS) \
	-DSHAREDIR=\"${datadir}/${PACKAGE}\"

AM_FFLAGS =  -I$(srcdir)/../dyres \
	-I$(srcdir)/../mcfm

noinst_LTLIBRARIES = libdyturbobase.la
libdyturbobase_la_SOURCES = resintegr.C ctintegr.C finintegr.C integr.C phasespace.C settings.C init.C cubacall.C plotter.C switch.C cuts.C usercuts.C pdf.C printsettings.C coupling.C clock_real.C banner.C
libdyturbobase_la_SOURCES += flavreduce.f gaussinit.f rescinit.f pdfmoments.f
libdyturbobase_la_SOURCES += AiMoments.cxx TLVUtils.cxx

lib_LTLIBRARIES = libdyturbo.la
libdyturbo_la_SOURCES =
libdyturbo_la_LIBADD = libdyturbobase.la \
	../integrands/libintegrands.la \
	../resum/libresum.la \
	../counterterm/libcounterterm.la \
	../vjfo/libvjfo.la \
	../dynnlo/libdynnlo.la \
	../dyres/libdyres.la \
	../mcfm/libmcfm.la \
	../dyqt/libdyqt.la \
	../pegasus/libpegasus.la \
	../dequad/libdequad.la \
	../cubature/libcubature.la

bin_PROGRAMS = dyturbo
dyturbo_SOURCES = dyturbo.C
dyturbo_LDADD = libdyturbo.la \
	../cernlib/libcern.la \
	-lgfortran \
	$(CUBA_LDFLAGS) $(LHAPDF_LDFLAGS) #-lgsl -lgslcblas -lifcore -limf -lifport -lintlc -lquadmath

if HAVE_ROOT
dyturbo_LDADD += $(ROOT_LDFLAGS)
AM_CPPFLAGS += $(ROOT_CFLAGS)

bin_PROGRAMS += merger

merger_SOURCES = merger.C
merger_LDADD = $(ROOT_LDFLAGS)
endif

pkginclude_HEADERS = resintegr.h ctintegr.h finintegr.h integr.h phasespace.h interface.h settings.h cubacall.h plotter.h init.h switch.h options.f flred.f gauss.f quadrules.f cuts.h coupling.h pdf.h printsettings.h clock_real.h fcomplex.h banner.h
