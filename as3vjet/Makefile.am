AUTOMAKE_OPTIONS = foreign

AM_FFLAGS = -ffixed-line-length-none -I$(srcdir)/../src -I$(srcdir)/../O3as3vj
AM_FCFLAGS = -cpp -ffixed-line-length-none -std=legacy -I$(srcdir)/../O3as3vj
AM_CPPFLAGS = -I$(srcdir)/../src -I$(srcdir)/../phasespace -I$(srcdir)/../resum

#disable Ofast as it gives unpredictable results for some functions
#FFLAGS   += -O3 #-ffast-math -fallow-store-data-races -fstack-arrays -fno-protect-parens #-funroll-loops
#FCFLAGS  += -O3 #-ffast-math -fallow-store-data-races -fstack-arrays -fno-protect-parens #-funroll-loops
#CXXFLAGS += -O3 #-ffast-math -fallow-store-data-races -fstack-arrays -fno-protect-parens #-funroll-loops

noinst_LTLIBRARIES = libas3vjetpar.la libas3vjet.la

.NOTPARALLEL: libas3vjet.la

libas3vjetpar_la_SOURCES = avh_olo.f90 \
	tdhpl.f \
	hplog.f \
	helamp.C \
	Li3.f \
	WGPLG.f \
	ptyrap.f \
	r.f \
	genparton.f \
	ampqqbgll.f \
	I1qq.f \
	I1gg.f \
	I2qq.f \
	I2gg.f \
	xbeam1bis.f \
	xbeam2bis.f \
	ap.f \
	auxfunctions.f \
	plus.f \
	soft1.f90 \
	as3vjet_dotem.f \
	as3vjet_includedipole.f \
	fillcoeff2dhpl210.f \
	fillcoeff2dhpl310.f \
	fillcoeff2dhpl320.f \
	fillcoeff2dhpl320e.f \
	fillcoeff2dhpl321.f \
	fillcoeff2dhpl321e.f \
	fillcoeff2dhpl321u.f \
	fillcoeff2dhplaux.f \
	spinorz.f \
	spinork.f \
	dipolesubx.f \
	gen_njets.f \
	genpt.f \
	phi3m.f \
	storedip.f \
	dipolesubx_new.f \
	checkndotp.f \
	dipolesub.f \
	transform.f \
	storeptilde.f \
	qqb_wp2jetx_new.f \
	qqb_wm2jetx_new.f \
	qqb_w2jet_v.f \
	qqb_w2jet_z.f \
	xwqqgg_v.f \
	qqbw2j_loop.f \
	a61LLL.f \
	a61LRL.f \
	atrLLL.f \
	atrLRL.f \
	a6treeg.f \
	a61g.f \
	fax.f \
	faxsl.f \
	a6g.f \
	a6.f \
	a6routine.f \
	t.f \
	vvg.f \
	fcc.f \
	fsc.f \
	fpp.f \
	fsl.f \
	fpm.f \
	vv.f \
	as3vjet_w2jetsq.f \
	amp_q_qbq_qb.f \
	qqb_w2jet_gs_new.f \
	qqb_w2jet_g.f \
	qqb_wp2jet_g.f \
	qqb_wm2jet_g.f \
	qqb_wm2jet_gs_new.f \
	qqb_wp2jet_gs_new.f \
	xwqqggg.f \
	qqb_w2jet_gvec.f \
	qqb_wm2jet_gvecx_new.f \
	qqb_wp2jet_gvecx_new.f \
	qqb_wp2jet_gs.f \
	w2jetnx.f \
	storecsv_px.f \
	storecsv_qx.f \
	qqb_w2jetx.f \
	ampwqqqqg.f \
	qqb_w2jet_gvecx.f \
	qqb_w2jet.f \
	ggZZmassamp_new.f \
	Acalc.f \
	LRcalc.f \
	ZZmassivebox.f \
	ZZmassiveboxtri.f \
	ZZmassivebub.f \
	ZZmassivetri.f \
	ggZZparsecheck.f \
	ZZbox1LL.f \
	ZZbox2LL.f \
	ZZD062x1x34LLmp.f \
	ZZD02x1x34LLmp.f \
	ZZC01x2LLmp.f \
	ZZC01x34LLmp.f \
	ZZC012x34LLmp.f \
	ggZZcapture.f \
	ZZmbc.f \
	ZZtri1_2LL.f \
	ZZtri1_34LL.f \
	ZZtri12_34LL.f \
	Tri3masscoeff.f \
	qqb_z2jetx_new.f \
	qqb_z2jet_z.f \
	atreez.f \
	a61.f \
	a62.f \
	a63.f \
	atree.f \
	a61z.f \
	a62z.f \
	a63z.f \
	xzqqgg_v.f \
	msq_z2jetx.f \
	xzqqgg.f \
	a6ax.f \
	fzip.f \
	fmt.f \
	BDKqqbggAxAmp.f \
	qqbZggtree.f \
	qqbggAxslCoeffs.f \
	qqbggAxbox3x12x4.f \
	qqbggAxbox3x4x12.f \
	qqbggAxtri12x3x456.f \
	qqbggAxtri123x4x56.f \
	A6axBDK.f \
	Ltfunctions.f \
	as3vjet_z2jetsq.f \
	as3vjet_subqcd.f \
	ampqqb_qqb.f \
	aqqb_zbb.f \
	qqb_z2jet_gs_new.f \
	qqb_z2jet_g.f \
	msq_ZqqQQg.f \
	xzqqggg.f \
	qqb_z2jet_gvecx_new.f \
	qqb_z2jet_gvec.f \
	z2jetsqn.f \
	subqcdn.f \
	amp_qqggg.f \
	makemb.f \
	makem.f \
	nagyqqqqg.f \
	qqb_z2jet.f \
	qqb_z1jet.f \
	virt5.f \
	A5NLO.f \
	A53.f \
	fvs.f \
	fvf.f \
	FPFPcc.f \
	FMPFcc.f \
	FPMFcc.f \
	FFPPcc.f \
	FFPMcc.f \
	FFMPcc.f \
	FPFMcc.f \
	FPPFcc.f \
	FMPFsc.f \
	FFPPsc.f \
	FFPMsc.f \
	FFMPsc.f \
	FPPFsc.f \
	FPMFsc.f \
	FPFPsc.f \
	FPFMsc.f \
	FFPMccT.f \
	FFPMscT.f \
	FFPMscTtilde.f \
	FFPMccTtilde.f \
	FPFMscT.f \
	FPFMccTtilde.f \
	Master1.f \
	Master2.f \
	Master3.f \
	Master2a.f \
	Master3a.f \
	genclust2.f \
	genclust_cone.f \
	jetreorder.f \
	genclust_hqrk.f \
	miscclust.f \
	is_functions.f \
	deltarj.f \
	checkjets.f \
	getbs.f \
	couplz.f \
	couplz_cms.f \
	as3vjet_dipoles.f \
	as3vjet_qqb_z1jet_v.f \
	as3vjet_interface.C

libas3vjet_la_SOURCES =  mod_ptveto.f90 \
	assemblejet.f90 \
	mod_Multichannel.f90 \
	smallnew.f \
	as3vjet_smalls.f \
	jet.f90 \
	mod_SCET.f90 \
	gen_lops.f \
	mod_nnlo_w1jet.f90 \
	mod_nnlo_z1jet.f90 \
	scetint.f \
	scet_xmsq.f \
	jettagging.f \
	maketaucut.f \
	swapjet.f \
	genclust_kt.f \
	mod_qcdloop_c.f \
	qqb_z2jet_v.f \
	loopI2_generic.f \
	loopI3_generic.f \
	loopI4_generic.f \
	F1anom.f \
	A6texact.f \
	ZZintegraleval.f \
	Bdiff.f \
	qqbZgg_floop.f \
	Ftexact.f \
	fmtfull.f \
	mod_Integration.f90 \
	as3vjet_virtint.f \
	as3vjet_realint.f


libas3vjet_la_LIBADD = libas3vjetpar.la

#vvconfig_m.f 
#mod_Scalevar.f90 

#nnlo_z1jet.mod: mod_nnlo_z1jet.f90
#	gfortran mod_nnlo_z1jet.f90 -o mod_nnlo_z1jet.o
#
#types.mod: types_mod.f90
#	gfortran types_mod.f90 -o typs_mod.o
#
#avh_olo.mod: avh_olo.f90
#	gfortran avh_olo.f90 -o avh_olo.o
#
#ptveto.mod: mod_ptveto.f90
#	gfortran mod_ptveto.f90 -o mod_ptveto.o
#
#mod_PDFerrors.o: types.mod
#jettagging.o: types.mod
#jettagging.o: types.mod
#mod_Integration.o: types.mod
#constants_mod.o: types.mod
#xspenz.o: types.mod
#mod_SCET.o: types.mod
#mod_Multichannel.o: vvconfig_m.o
#mod_nnlo_z1jet_hfun.o: constants_mod.o
#mod_nnlo_z1jet.o: multichannel.mod scet.mod scet_jet.mod nnlo_z1jet_hfun.mod # mod_SCET.o jet.o
#mod_nnlo_w1jet.o: scet_jet.mod types.mod mod_nnlo_z1jet_hfun.mod scet.mod
#gen_lops.o: multichannel.mod
#scet_xmsq.o: nnlo_z1jet.mod
#scetint.o: mod_SCET.o types.mod
#gen3taucut.o: ptveto.mod
#as3vjet_virtint.o: mod_Scalevar.o
#qqb_z2jet_v.o: mod_qcdloop_c.o
#qqbZgg_floop.o: loopI2_generic.mod loopI3_generic.mod loopI4_generic.mod
#loopI2_generic.o: avh_olo.mod
#loopI3_generic.o: avh_olo.mod
#loopI4_generic.o: avh_olo.mod
#genclust_kt.o: jettagging.o
#ZZintegraleval.o: loopI2_generic.mod

#constants.mod
#multichannel.mod
#nnlo_z1jet_hfun.mod
#pdferrors.mod
#ptveto.mod
#scet_jet.mod
#scet.mod
#types.mod
#vvconfig_m.mod

CLEANFILES = *.mod


#mod_nnlo_h1jet.f90 
#mod_MCFMStorage.f90 lhapdf_fortran.f90 

pkginclude_HEADERS = tiny.f useet.f first.f reset.f nodecay.f new_pspace.f zerowidth.f breit.f dynamicscale.f toploops.f cutoff.f mpicommon.f Cabibbo.f epinv2.f epinv.f ckm.f facscale.f scale.f zprods_decl.f sprods_com.f zprods_com.f as3vjet_ewcharge.f as3vjet_zcouple.f as3vjet_ewcouple.f as3vjet_qcdcouple.f as3vjet_masses.f nwz.f npart.f ptilde.f hdecaymode.f frag.f maxd.f cplx.h hpls.f plabel.f nqcdjets.f beamtype.f ipsgen.f kpart.f mxdim.f as3vjet_constants.f nflav.f scet_beta.f zeta.f scet_const.f nf.f limits.f mxpart.f kprocess.f types.f nproc.f xmin.f x1x2.f energy.f taucut.f scheme.f zcouple_cms.f blha.f noglue.f agq.f PR_new.f PR_cs_new.f PR_h2j.f PR_twojet.f PR_stop.f PR_mix.f msq_cs.f msq_struc.f msq_mix.f lc.f maxwt.f b0.f wts_bypart.f nores.f stopbmass.f flags.f outputoptions.f runstring.f bypart.f initialscales.f ppmax.f hbbparams.f ewcorr.f couple.f nlooprun.f badpoint.f anomcoup.f TRbadpoint.f dipolescale.f incldip.f jetlabel.f pp.f ppwp2j.f ppwm2j.f heldefs.f alfacut.f mmsq_cs.f qqbggintnames.f docheck.f loopI2_inc.f loopI2c_inc.f scalarselect.f loopI3_inc.f loopI3c_inc.f loopI4_inc.f loopI4c_inc.f ggZZcomputemp.f ZZdlabels.f ZZclabels.f blabels.f ggZZintegrals.f tri123x4x56coeffs.f clustering.f bbproc.f notag.f is_functions_com.f userap.f removebr.f qqgg.f msqv_cs.f mmsqv_cs.f colstruc.f vegas_common.f as3vjet_interface.h helamp.h
