
AUTOMAKE_OPTIONS = foreign
#  -I$(srcdir)/../dequad
#	$(CUBA_CPPFLAGS) \
#	-DSHAREDIR=\"${datadir}/${PACKAGE}\"
AM_CPPFLAGS = -fpermissive -Wl,Bstatic \
		-I$(srcdir)/../src \
		-I$(srcdir)/../dequad \
		-I$(srcdir)/../bequad \
		-I$(srcdir)/../besche \
		-I$(srcdir)/../minprescription \
		-I$(srcdir)/../specialfunctions \
		-I$(srcdir)/../bessel \
		-I$(builddir)/../cubature \
		-I$(srcdir)/../dyres \
		-I$(srcdir)/../mcfm \
		-I$(srcdir)/../phasespace \
		$(LHAPDF_CPPFLAGS) $(CUBA_CPPFLAGS)

AM_FFLAGS =  -I$(srcdir)/../dyres \
	     -I$(srcdir)/../phasespace -I$(srcdir)/../mcfm -I$(srcdir)/../src

noinst_LTLIBRARIES = libresum.la
libresum_la_SOURCES = mellinint.C mesq.C rapint.C anomalous.C resconst.C hcoefficients.C hcoeff.C pdfevol.C besselint.C bessel.f resint.C pegasus.C sudakov.f sudakovff.C npff.C

pkginclude_HEADERS = mellinint.h mesq.h rapint.h anomalous.h resconst.h hcoefficients.h hcoeff.h pdfevol.h besselint.h resint.h pegasus.h sudakov_inc.f pdfn_inc.f scales_inc.f sigmaij_inc.f sudakovff.h npff.h
