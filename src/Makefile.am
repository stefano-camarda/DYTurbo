AUTOMAKE_OPTIONS = foreign
AM_CPPFLAGS = -fpermissive \
	-I$(srcdir)/../include \
	-I$(srcdir)/../dequad \
	-I$(srcdir)/../sandia \
	-I$(srcdir)/../smolpack \
	-I$(srcdir)/../dyres \
	-I$(srcdir)/../mcfm \
	-I$(srcdir)/../resum \
	-I$(srcdir)/../counterterm \
	-I$(srcdir)/../vjet \
	-I$(srcdir)/../born \
	-I$(builddir)/../cubature \
	-I$(srcdir)/../phasespace \
	-I$(srcdir)/../histo \
	-I$(srcdir)/../user \
	-Wl,Bstatic \
	$(CUBA_CPPFLAGS) $(LHAPDF_CPPFLAGS) \
	-DSHAREDIR=\"${datadir}/${PACKAGE}\"

AM_FFLAGS =  -I$(srcdir)/../dyres \
	-I$(srcdir)/../mcfm

noinst_LTLIBRARIES = libdyturbobase.la
libdyturbobase_la_SOURCES = dyturbo.C resintegr.C ctintegr.C finintegr.C bornintegr.C \
			    omegaintegr.C interface.C settings.C init.C cubacall.C switch.C  \
			    pdf.C coupling.C clock_real.C banner.C print_table.C abint.C clenshawcurtisrules.C gaussrules.C numbers.C chebyshev.C codata.C #solvew.C
libdyturbobase_la_SOURCES += flavreduce.f gaussinit.f rescinit.f pdfmoments.f chaplin.f
libdyturbobase_la_SOURCES += AiMoments.cxx TLVUtils.cxx

lib_LTLIBRARIES = libdyturbo.la
libdyturbo_la_SOURCES =
libdyturbo_la_LIBADD = libdyturbobase.la \
		       ../histo/libhisto.la \
		       ../phasespace/libphasespace.la \
		       ../integrands/libintegrands.la \
		       ../resum/libresum.la \
		       ../counterterm/libcounterterm.la \
		       ../vjet/libvjet.la \
		       ../v2j/libv2j.la \
		       ../born/libborn.la \
		       ../dynnlo/libdynnlo.la \
		       ../dyres/libdyres.la \
		       ../mcfm/libmcfm.la \
		       ../pegasus/libpegasus.la \
		       ../dequad/libdequad.la \
		       ../bessel/libbessel.la \
		       ../sandia/libsandia.la \
		       ../smolpack/libsmolpack.la \
		       ../besche/libbesche.la \
		       ../minprescription/libminprescription.la \
		       ../specialfunctions/libspecialfunctions.la \
		       ../cubature/libcubature.la \
		       ../cernlib/libcern.la \
		       -lgfortran -lblas \
		       $(LHAPDF_LDFLAGS) #-lgsl -lgslcblas -lifcore -limf -lifport -lintlc -lquadmath


bin_PROGRAMS = dyturbo
dyturbo_SOURCES = main.C
dyturbo_LDADD = libdyturbo.la $(CUBA_LDFLAGS) $(CHAPLIN_LDFLAGS)

if HAVE_ROOT
dyturbo_LDADD += $(ROOT_LDFLAGS)
AM_CPPFLAGS += $(ROOT_CFLAGS)

bin_PROGRAMS += dyturbo-merger

dyturbo_merger_SOURCES = merger.C
dyturbo_merger_LDADD = $(ROOT_LDFLAGS)
endif

pkginclude_HEADERS = dyturbo.h resintegr.h ctintegr.h finintegr.h bornintegr.h omegaintegr.h \
		     interface.h settings.h cubacall.h switch.h options.f flred.f \
		     gauss.f quadrules.f coupling.h pdf.h clock_real.h \
		     fcomplex.h banner.h cxxopts.hpp isnan.h handy_typdefs.h abint.h \
		     clenshawcurtisrules.h gaussrules.h numbers.h parton.h AiMoments.h TLVUtils.h chebyshev.h codata.h #solvew.h
