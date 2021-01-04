
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
		-I$(srcdir)/../dynnlo \
		-I$(srcdir)/../phasespace \
		-I$(srcdir)/../harmonicsums \
		-I$(srcdir)/../polylog \
		$(LHAPDF_CPPFLAGS) $(CUBA_CPPFLAGS) $(CERES_CPPFLAGS) $(EIGEN_CPPFLAGS) $(POLYLOG_CPPFLAGS) \
		-DSHAREDIR=\"${datadir}/${PACKAGE}\"

AM_FFLAGS =  -I$(srcdir)/../dyres \
	     -I$(srcdir)/../phasespace \
	     -I$(srcdir)/../mcfm \
             -I$(srcdir)/../src

noinst_LTLIBRARIES = libresum.la
libresum_la_SOURCES = mellinint.C mesq.C rapint.C anomalous.C resconst.C hcoefficients.C hcoeff.C expc.C alphas.C pdfevol.C evolnative.C evolnum.C besselint.C bessel.f resint.C pegasus.C sudakov.f sudakovff.C npff.C blim.C ccoeff.C icoeff.C ifunc.C gint.C pmom.C muf.C hcoeff_check.C cmom.C

pkginclude_HEADERS = mellinint.h mesq.h rapint.h anomalous.h resconst.h hcoefficients.h hcoeff.h expc.h alphas.h pdfevol.h evolnative.h evolnum.h besselint.h resint.h pegasus.h sudakov_inc.f pdfn_inc.f scales_inc.f sigmaij_inc.f sudakovff.h npff.h blim.h ccoeff.h icoeff.h ifunc.h gint.h pmom.h muf.h hcoeff_check.h cmom.h
