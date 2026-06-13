AUTOMAKE_OPTIONS = foreign
AM_CPPFLAGS = -fpermissive \
	-I$(srcdir)/../include \
	-I$(srcdir)/../dequad \
	-I$(srcdir)/../bequad \
	-I$(srcdir)/../zeta \
	-I$(srcdir)/../sandia \
	-I$(srcdir)/../sparse \
	-I$(srcdir)/../smolpack \
	-I$(srcdir)/../dyres \
	-I$(srcdir)/../mcfm \
	-I$(srcdir)/../as3vjet \
	-I$(srcdir)/../resum \
	-I$(srcdir)/../counterterm \
	-I$(srcdir)/../qt \
	-I$(srcdir)/../hpoly \
	-I$(srcdir)/../crundec \
	-I$(srcdir)/../vjet \
	-I$(srcdir)/../born \
	-I$(builddir)/../cubature \
	-I$(srcdir)/../phasespace \
	-I$(srcdir)/../histo \
	-I$(srcdir)/../user \
	-I$(srcdir)/../polylog \
	$(HELLX_CPPFLAGS) $(HELL_CPPFLAGS) $(CUBA_CPPFLAGS) $(LHAPDF_CPPFLAGS) $(CERES_CPPFLAGS) $(EIGEN_CPPFLAGS) $(NEOPDF_CPPFLAGS) \
	-DSHAREDIR=\"${datadir}/${PACKAGE}\"

AM_FFLAGS =  -I$(srcdir)/../dyres \
	-I$(srcdir)/../mcfm

noinst_LTLIBRARIES = libdyturbobase.la
libdyturbobase_la_SOURCES = dyturbo.C resintegr.C ctintegr.C fpcintegr.C finintegr.C bornintegr.C fointegr.C \
			    omegaintegr.C interface.C settings.C init.C cubacall.C switch.C  \
			    pdf.C coupling.C clock_real.C banner.C print_table.C abint.C clenshawcurtisrules.C gaussrules.C sparserules.C \
			    propagator.C numbers.C chebyshev.C codata.C scales.C taucut.C mellinpdf.C constants.C psi.C dilog.C isnan.C #solvew.C
if HAVE_CERES
libdyturbobase_la_SOURCES += cerespdf.C
endif
libdyturbobase_la_SOURCES += flavreduce.f gaussinit.f rescinit.f pdfmoments.f mjj.f fortranflush.f
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
		       ../qt/libqt.la \
		       ../hpoly/libhpoly.la \
		       ../mcfm/libmcfm.la \
		       ../O3as3vj/libO3as3vj.la \
		       ../as3vjet/libas3vjet.la \
		       ../pegasus/libpegasus.la \
		       ../dequad/libdequad.la \
		       ../bequad/libbequad.la \
		       ../polylog/libpolylog.la \
		       ../crundec/libcrundec.la \
		       ../bessel/libbessel.la \
		       ../fbt/libfbt.la \
		       ../sandia/libsandia.la \
		       ../sparse/libsparse.la \
		       ../smolpack/libsmolpack.la \
		       ../besche/libbesche.la \
		       ../minprescription/libminprescription.la \
		       ../specialfunctions/libspecialfunctions.la \
		       ../cubature/libcubature.la \
		       -lgfortran \
		       $(LHAPDF_LDFLAGS) -lLHAPDF \
		       $(NEOPDF_LDFLAGS) \
		       $(CUBA_LDFLAGS) \
		       $(QCDLOOP_LDFLAGS) \
		       $(CHAPLIN_LDFLAGS) \
		       $(CERES_LDFLAGS) \
		       $(POLYLOG_LDFLAGS) \
		       $(HELLX_LDFLAGS) \
		       $(HELL_LDFLAGS) \
		       -lquadmath \
	               -lm
#		       ../cernlib/libcern.la
#		       ../zeta/libzeta.la
#		       ../harmonicsums/libhs.la
#-lgsl -lgslcblas -lifcore -limf -lifport -lintlc -lquadmath
#-llapack -lblas -lgslcblas

#libdyturbo_la_LDFLAGS = -all-static -Wl,-Bstatic

if HAVE_ROOT
libdyturbo_la_LIBADD += $(ROOT_LDFLAGS)
endif

if HAVE_CERES
libdyturbo_la_LIBADD += -lgsl -lgslcblas
endif

if HAVE_HPL
libdyturbo_la_LIBADD += -lginac
endif


bin_PROGRAMS = dyturbo
dyturbo_SOURCES = main.C
#dyturbo_LDADD = -Wl,-Bstatic -ldyturbo -Wl,-Bdynamic
dyturbo_LDADD = libdyturbo.la

dyturbo_LDFLAGS = -static $(NEOPDF_LDFLAGS) #-all-static #-Wl,-Bstatic

if HAVE_ROOT
#dyturbo_LDADD += $(ROOT_LDFLAGS)
AM_CPPFLAGS += $(ROOT_CFLAGS)

bin_PROGRAMS += dyturbo-merger

dyturbo_merger_SOURCES = merger.C
dyturbo_merger_LDADD = $(ROOT_LDFLAGS)

bin_PROGRAMS += dyturbo-csk

dyturbo_csk_SOURCES = csk.C
dyturbo_csk_LDADD = libdyturbo.la $(ROOT_LDFLAGS)
endif

pkginclude_HEADERS = dyturbo.h resintegr.h ctintegr.h fpcintegr.h finintegr.h bornintegr.h fointegr.h omegaintegr.h \
		     interface.h settings.h cubacall.h switch.h options.f flred.f \
		     gauss.f quadrules.f coupling.h pdf.h clock_real.h \
		     fcomplex.h banner.h cxxopts.hpp isnan.h handy_typdefs.h abint.h \
		     clenshawcurtisrules.h gaussrules.h sparserules.h numbers.h parton.h AiMoments.h TLVUtils.h chebyshev.h codata.h \
		     propagator.h scales.h taucut.h mellinpdf.h constants.h cerespdf.h psi.h accel.h flush.h dilog.h  #solvew.h
