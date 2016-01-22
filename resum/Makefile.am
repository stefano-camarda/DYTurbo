AUTOMAKE_OPTIONS = foreign
#  -I$(srcdir)/../dequad
#	$(CUBA_CPPFLAGS) $(LHAPDF_CPPFLAGS) \
#	-DSHAREDIR=\"${datadir}/${PACKAGE}\"
AM_CPPFLAGS = -fpermissive -Wl,Bstatic -I$(srcdir)/../src

#AM_FFLAGS =  -I$(srcdir)/../dyres -I$(srcdir)/../mcfm #-I$(srcdir)/../src

noinst_LTLIBRARIES = libresum.la
libresum_la_SOURCES = gaussrules.C mellinint.C mesq.C rapint.C anomalous.C resconst.C hcoefficients.C pdfevol.C

pkginclude_HEADERS = gaussrules.h mellinint.h mesq.h rapint.h anomalous.h resconst.h hcoefficients.h pdfevol.h
