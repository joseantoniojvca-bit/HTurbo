AUTOMAKE_OPTIONS = foreign
#  -I$(srcdir)/../dequad
#	$(CUBA_CPPFLAGS) \
#	-DSHAREDIR=\"${datadir}/${PACKAGE}\"
AM_CPPFLAGS = -fpermissive \
		-I$(srcdir)/../src \
		-I$(srcdir)/../dequad \
		-I$(srcdir)/../bequad \
		-I$(srcdir)/../besche \
		-I$(srcdir)/../minprescription \
		-I$(srcdir)/../specialfunctions \
		-I$(srcdir)/../bessel \
		-I$(srcdir)/../fbt \
		-I$(builddir)/../cubature \
		-I$(srcdir)/../dyres \
		-I$(srcdir)/../mcfm \
		-I$(srcdir)/../dynnlo \
		-I$(srcdir)/../phasespace \
		-I$(srcdir)/../harmonicsums \
		-I$(srcdir)/../polylog \
		-I$(srcdir)/../crundec \
		$(HELLX_CPPFLAGS) $(HELL_CPPFLAGS) $(LHAPDF_CPPFLAGS) $(NEOPDF_CPPFLAGS) $(CUBA_CPPFLAGS) $(CERES_CPPFLAGS) $(EIGEN_CPPFLAGS) $(POLYLOG_CPPFLAGS) \
		-DSHAREDIR=\"${datadir}/${PACKAGE}\"

AM_FFLAGS =  -I$(srcdir)/../dyres \
	     -I$(srcdir)/../phasespace \
	     -I$(srcdir)/../mcfm \
             -I$(srcdir)/../src

noinst_LTLIBRARIES = libresum.la

libresum_la_SOURCES = hcoeff_check.C mellinint.C mesq.C rapint.C anomalous.C resconst.C hcoefficients.C hcoeff.C expc.C alphas.C alphaqed.C pdfevol.C evolnative.C evolnum.C hell.C besselint.C bessel.f resint.C pegasus.C sudakov.f sudakovff.C npff.C gj.C gk.C blim.C ccoeff.C icoeff.C ifunc.C gint.C pmom.C muf.C mures.C cmom.C melfun.C c2gqmom.C c3gqmom.C dymellinh2.f gamma4ns.C gamma4ns_NF.C gamma4nsp.C gamma4nsm.C gamma4nsv.C gamma4pqqps.C gamma4pqgps.C gamma4pgqps.C gamma4pggps.C gamma4sg.f commute.C cevol.C evoltrunc.C bstar.C btilde.C muffast1d.C muffast2d.C muf1d.C muf2d.C mufeqqb.C mufeqg.C mufqqb.C mufqg.C mufqq.C mufqqp.C mufqqbp.C mufgg.C mufqbq.C sgfun.C gamma1.C gamma2.C rad.C cusp.C mellinh2Higgs.f H1_Int.f Cl2.f90  Cl3.f90 Cl4.f90 Cl5.f90 Cl6.f90 Li2.f90 Li3.f90 Li4.f90 Li5.f90 Li6.f90 log.f90

#for the Higgs
#libresum_la_SOURCES += c2ggmom.C c3gqmom.C c3ggmom.C

pkginclude_HEADERS = mellinint.h mesq.h rapint.h anomalous.h resconst.h hcoefficients.h hcoeff.h expc.h alphas.h alphaqed.h pdfevol.h evolnative.h evolnum.h hell.h besselint.h resint.h pegasus.h sudakov_inc.f pdfn_inc.f scales_inc.f sigmaij_inc.f sudakovff.h npff.h blim.h ccoeff.h icoeff.h ifunc.h gint.h pmom.h muf.h mures.h hcoeff_check.h cmom.h melfun.h gamma4ns.h gamma4ns_NF.h gamma4pqqps.h gamma4pqgps.h gamma4pgqps.h gamma4pggps.h gamma4sg.h commute.h cevol.h evoltrunc.h bstar.h btilde.h sgfun.h gamma1.h gamma2.h dytmatrix.h rad.h cusp.h
