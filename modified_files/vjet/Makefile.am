AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = -fpermissive -Wl,Bstatic \
	-I$(srcdir)/../src \
	-I$(srcdir)/../resum \
	-I$(builddir)/../cubature \
	-I$(srcdir)/../mcfm \
	-I$(srcdir)/../phasespace \
	-I$(srcdir)/../resum \
	-I$(srcdir)/../histo \
	-I$(srcdir)/../user \
	$(LHAPDF_CPPFLAGS) $(CUBA_CPPFLAGS) $(NEOPDF_CPPFLAGS)

AM_FFLAGS = -I$(srcdir)/../mcfm

noinst_LTLIBRARIES = libvjet.la
libvjet_la_SOURCES = vjint.C vjloint.C sing.C delta.C luminosity.C params.cc mi-rand-t.cc user_2fot.cc threevector.cc lorentzvector.cc phasespace_n2i2f0.cc dipole_emission.cc Evento.cc hhc2ph1jet.cc ampq2g1p2.cc ampq2g2p2.cc ampq2g3.cc ampbase.cc li2.cc cl2.cc color.cc ampq4p2.cc split.cc weight.cc hhc2phjet.cc flkern.cc mstwpdf_cpp.cc process.cc H1.cpp ampBORN.cpp qt_splitt.cc  Itilde.cpp wadp.cpp InteItilde.cpp H2f.cpp integral-hhc2ph.cc Ampi.cc

libvjet_la_SOURCES += utils.f abcdefunc.f abcdefunc_qed.f

pkginclude_HEADERS = vjint.h fodyqt_inc.f functions_inc.f luminosities_inc.f scales2_inc.f vjloint.h luminosity.h dyqcd.f
