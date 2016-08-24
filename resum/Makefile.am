AUTOMAKE_OPTIONS = foreign
#  -I$(srcdir)/../dequad
#	$(CUBA_CPPFLAGS) \
#	-DSHAREDIR=\"${datadir}/${PACKAGE}\"
AM_CPPFLAGS = -fpermissive -Wl,Bstatic -I$(srcdir)/../src -I$(srcdir)/../dequad -I$(srcdir)/../cubature -I$(srcdir)/../mcfm -I$(srcdir)/../phasespace $(LHAPDF_CPPFLAGS) $(CUBA_CPPFLAGS)

AM_FFLAGS =  -I$(srcdir)/../dyres -I$(srcdir)/../phasespace #-I$(srcdir)/../mcfm #-I$(srcdir)/../src

noinst_LTLIBRARIES = libresum.la
libresum_la_SOURCES = gaussrules.C mellinint.C mesq.C rapint.C anomalous.C resconst.C hcoefficients.C pdfevol.C besselint.C bessel.f resint.C pegasus.C sudakov.f

pkginclude_HEADERS = gaussrules.h mellinint.h mesq.h rapint.h anomalous.h resconst.h hcoefficients.h pdfevol.h besselint.h resint.h pegasus.h parton.h sudakov_inc.f pdfn_inc.f scales_inc.f sigmaij_inc.f
