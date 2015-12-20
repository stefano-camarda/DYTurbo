AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = -fpermissive $(LHAPDF_CPPFLAGS)

noinst_LTLIBRARIES = liborigintegrands.la

liborigintegrands_la_SOURCES = main2.f intde2.f alfamz.f varios.f csection.f newton1.f integration.f mellinh2.f pdfset_lhapdf.f pdf_lhapdf.f mdata.f dyinit.f setup-RES.f coupling.f strcat.f ran1.f ckmfill.f plotter.f couplz.f branch.f mbook.f evtwriter.f ptyrap.f

bin_PROGRAMS = dyrestest
dyrestest_SOURCES = dyrestest.C
dyrestest_LDADD = liborigintegrands.la ../cernlib/libcern.la \
	-lgfortran $(LHAPDF_LDFLAGS)

