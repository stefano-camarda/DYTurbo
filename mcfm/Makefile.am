AUTOMAKE_OPTIONS = foreign

AM_FFLAGS = -I$(srcdir)/../dyres

noinst_LTLIBRARIES = libmcfm.la

#original files
libmcfm_la_SOURCES = qqb_w1jet_gs.f  qqb_w1jet_v.f  qqb_w1jet_z.f  qqb_w2jet.f  qqb_z1jet.f  qqb_z1jet_gs.f  qqb_z1jet_v.f  qqb_z1jet_z.f  qqb_z2jet.f z2jetsq.f w2jetsq.f A51.f A52.f A5NLO.f ampqqb_qqb.f aqqb_zbb.f boost.f branch.f breitw.f ckmfill.f lfunctions.f dclaus.f i3m.f transform.f ddilog.f lnrat.f

#modified files
libmcfm_la_SOURCES +=  alfamz_lhapdf.f

pkginclude_HEADERS = agq.f alfacut.f b0.f ckm.f ckm1.f clustering.f
