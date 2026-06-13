#include "ctint.h"

#include "qtint.h"
#include "mesq.h"
#include "parton.h"
#include "settings.h"
#include "omegaintegr.h"
#include "pdf.h"
#include "scales.h"
#include "phasespace.h"
#include "abint.h"
#include "resconst.h"
#include "constants.h"
#include "dyres_interface.h"
#include "dynnlo_interface.h"
#include "qt_interface.h"
#include "hpoly_interface.h"
#include "switch.h"
#include "isnan.h"
#include "alphaqed.h"

#include <iostream>
#include <string.h>
#include <math.h>

using namespace std;
using namespace resconst;
using namespace constants;

double ctint::intfuncqt(double fxip, double fxi0, double xxi0, double z, double Preg, double PD0, double PD1, double PD2, double Pdelta)
{
  //This function can be used if PD0,PD1 and PD2 don't depend on z
  double commu,commufuori; //,D0int,D1int,D2int;
  double IntFuncqT,IntPgggggg;
  //external D0int,D1int,D2int

  commu=-log(xxi0)*(fxip-z*fxi0)/(1.0-z);
  IntPgggggg=commu*(PD0 + PD1*log(1.0-z) + PD2*pow(log(1.0-z),2));
  commufuori=fxi0*(Pdelta - PD0*d0int_(xxi0) - PD1*d1int_(xxi0));
  commufuori=commufuori - fxi0*(PD2*d2int_(xxi0));
  IntPgggggg=IntPgggggg + commufuori;
  IntPgggggg=IntPgggggg - log(xxi0)*Preg*fxip;
       
  return IntPgggggg;
}


//Counterterm to be subtracted from V+j to get a finite cross section at qt->0
void ctint::calc(double costh, double m, double qt, double y, int mode, double f[], int order, bool fpc)
{

	//cout << "ctint.C" << endl;
	if (opts.nproc == 4) {

		for (int npdf = 0; npdf < opts.totpdf; npdf++)
			f[npdf] = 0.;

		/*
		 * //limit the integration to qt = m
		 * double jac = 1;
		 * double qtp = qt;
		 * if (qtp >= m*0.999)
		 *	return 0.;
		 * else
		 *	qt = qtp/sqrt(1-pow(qtp/m,2));
		 * jac = pow(m/sqrt(m*m-qtp*qtp),3);
		 */

		///////////////////////////////////////////////////////
		double m2 = m*m;
		double qt2 = qt*qt;
		double exppy = exp(y);
		double expmy = 1./exppy;
		double tau = sqrt(m2/pow(opts.sroot,2));

		//amplitudes are set in qtint::calc(), no need to calculate them
		/*
		 * mesq::setpropagators(m);
		 * double cthmom0, cthmom1, cthmom2;
		 * if (mode == 0)
		 *	{
		 *	//apply here lepton cuts
		 *	cthmom0 = 1.;
		 *	cthmom1 = costh;
		 *	cthmom2 = pow(costh,2);
	}
	else if (mode == 1)
	{
	//phasespace::set_mqtyphi(m, 0., y);
	phasespace::set_qt(0.);
	omegaintegr::genV4p();
	omegaintegr::cthmoments(cthmom0,cthmom1,cthmom2);
	}

	mesq::setmesq(cthmom0, cthmom1, cthmom2);
	//cout << cthmom0 << "  " << cthmom1 << "  " << cthmom2 << endl;
	*/

		//calculate Bjorken x1 x2
		double x1 = tau*exppy;
		double x2 = tau*expmy;

		if (x1 >= 1 || x2 >= 1)
			return;

		//Set scales
		scales::set(m);
		scales::mcfm();
		double muf = scales::fac;
		double mur = scales::ren;

		/*
		 * //Set factorization scale
		 * double muf, mur;
		 * if (opts.dynamicscale)
		 *	{
		 *	muf = m*opts.kmufac;
		 *	mur = m*opts.kmuren;
		 *	double mur2 = mur*mur;
		 *	scaleset_(mur2); //set renormalization and factorization scales, and calculate ason2pi and ason4pi
	}
	else
	{
	muf = opts.rmass*opts.kmufac;
	mur = opts.rmass*opts.kmuren;
	}
	*/

		//a-parameter of the resummation scale

		//dynamic case
		if (opts.fmures > 0)
			a_param_.a_param_ = 1./opts.kmures;

		//for fixed resummation scale need to recompute a_param
		else
			a_param_.a_param_ = m/scales::res;

		//PDFs
		double fx1[2*MAXNF+1],fx2[2*MAXNF+1];
		fdist_(opts.ih1,x1,muf,fx1);
		fdist_(opts.ih2,x2,muf,fx2);
		//////////////////////////////////////////////////////////

		double LR, LF, LQ;
		if (opts.order >= 2)
			LR = log(m2/pow(mur,2));
		LF = log(m2/pow(muf,2));
		LQ = 2.*log(a_param_.a_param_);

		/*
		 * //Calculate large logs LL1-LL4
		 * double LL1,LL2,LL3,LL4;
		 * if (mode == 0 || mode == 1)
		 *	{
		 *	//xmio is used in besselkfast for Itilde
		 *	xmio_.xmio_ = sqrt(qt2/(m2/pow(a_param_.a_param_,2)));
		 *
		 *	//LL1,LL2,LL3,LL4: large log (squared) corresponding to eq. (136)
		 *	//In this way normalization is fixed to dsigma/dqt2
		 *
		 *	//depends on xmio, which depends on qt2 and m2
		 *	LL1 = itilde_(one)/pow(m2,2)*pow(a_param_.a_param_,2);
		 *	LL2 = itilde_(two)/pow(m2,2)*pow(a_param_.a_param_,2);
		 *	LL3 = itilde_(three)/pow(m2,2)*pow(a_param_.a_param_,2);
		 *	LL4 = itilde_(four)/pow(m2,2)*pow(a_param_.a_param_,2);
	}
	*/

		//Scaled momentum fractions
		double lx1, lx2;
		double pqqintx1, pqqintx2;
		double d0intx1, d0intx2;
		double d1intx1, d1intx2;
		if (opts.order >= 1)
		{
			lx1 = log(x1);
			lx2 = log(x2);
			pqqintx1 = pqqint_(x1);
			pqqintx2 = pqqint_(x2);
			// Higgs production
			d0intx1 = d0int_(x1);
			d0intx2 = d0int_(x2);
		}
		if (opts.order >= 2)
		{
			d1intx1 = d1int_(x1);
			d1intx2 = d1int_(x2);
		}

		//Preliminary loop for caching
		double cz1[abint::abdim];
		double cz2[abint::abdim];
		double oz1[abint::abdim];
		double oz2[abint::abdim];
		// --> cache the logs, and pass them as parameters to the CxxPxx, Pxxxx functions
		double lz1[abint::abdim];
		double lz2[abint::abdim];
		double l1z1[abint::abdim];
		double l1z2[abint::abdim];
		double lzoz1[abint::abdim];
		double lzoz2[abint::abdim];
		if (opts.order >= 1)
			for (int ab = 0; ab < abint::abdim; ab++)
			{
				cz1[ab] = pow(x1,abint::abx[ab]);
				cz2[ab] = pow(x2,abint::abx[ab]);
				oz1[ab] = 1./(1.-cz1[ab]);
				oz2[ab] = 1./(1.-cz2[ab]);
				lz1[ab] = log(cz1[ab]);
				lz2[ab] = log(cz2[ab]);
				l1z1[ab] = log(1.-cz1[ab]);
				l1z2[ab] = log(1.-cz2[ab]);
				lzoz1[ab] = log(1.-cz1[ab])/(1.-cz1[ab]);
				lzoz2[ab] = log(1.-cz2[ab])/(1.-cz2[ab]);
			}

			// skip PDF loop in the preconditioning phase
			int maxpdf=0;
			if (dofill_.doFill_ != 0) maxpdf = opts.totpdf;

			// start PDF loop
			for (int npdf = 0; npdf < maxpdf; npdf++) {

				dysetpdf_(npdf);

				//cache scaled PDFs
				double fx1p[abint::abdim][2*MAXNF+1],fx2p[abint::abdim][2*MAXNF+1];
				if (opts.order >= 1)
					for (int ab = 0; ab < abint::abdim; ab++)
					{
						double fx1temp[2*MAXNF+1],fx2temp[2*MAXNF+1];
						double xx1 = pow(x1,(1-abint::abx[ab]));
						double xx2 = pow(x2,(1-abint::abx[ab]));
						fdist_(opts.ih1,xx1,muf,fx1temp);
						fdist_(opts.ih2,xx2,muf,fx2temp);
						memcpy(fx1p[ab], fx1temp, (2*MAXNF+1)*sizeof(double));
						memcpy(fx2p[ab], fx2temp, (2*MAXNF+1)*sizeof(double));
					}

					double sumfx1p[abint::abdim];
					double sumfx2p[abint::abdim];
					if (opts.order >= 2)
						for (int ab = 0; ab < abint::abdim; ab++)
						{
							sumfx1p[ab] = 0.;
							sumfx2p[ab] = 0.;
							for (int f = 0; f < MAXNF; f++)
							{
								sumfx1p[ab] += fx1p[ab][f]+fx1p[ab][parton::charge_conj(parton::pdgid(f))];
								sumfx2p[ab] += fx2p[ab][f]+fx2p[ab][parton::charge_conj(parton::pdgid(f))];
							}
						}

						// Start calculation
						double asopi = qcdcouple_.ason2pi_*2.;

						//loop on born subprocesses, i.e. born incoming partons ij
						//double lumi[mesq::totpch];
						double sig11[mesq::totpch];
						double sig12[mesq::totpch];
						double sig21[mesq::totpch];
						double sig22[mesq::totpch];
						double sig24[mesq::totpch];
						double sig23[mesq::totpch];
						for (int sp = 0; sp < mesq::totpch; sp++)
						{
							//simplify notation
							//double bornmesqij = real(mesq::mesqij[sp]); //born level amplitudes
							parton::pdgid i = mesq::pid1[sp];         //parton 1
							parton::pdgid j = mesq::pid2[sp];         //parton 2
							parton::pdgid g = parton::G;              //gluon
							parton::pdgid im = parton::charge_conj(i);
							parton::pdgid jm = parton::charge_conj(j);

							//LO term (there is no counterterm at LO...)..................................................................................................

							//Simplest term without convolutions
							double tdelta = fx1[5]*fx2[5];

							//NLO terms
							double th1st = 0; //this piece is used only at NNLO?
							double th1stF = 0;
							double th1stQ = 0; //this piece is used only at NNLO

							//H1st delta term
							th1st += 2*resconst::C1ggdelta*tdelta;

							//add resummation scale dependence
							th1stQ += -(resconst::B1g+resconst::A1g/2.*LQ)*LQ*tdelta; //this piece is used only at NNLO

							//alfa loop (first leg)
							for (int a = 0; a < abint::abdim; a++) {

								if (cz1[a] >= 1) continue;

								//H1st muf dependence, gammaqq and gammaqg:
								// th1stF += (-lx1*((fx1p[a][i]-fx1[i]*cz1[a])*pqq_(cz1[a])+fx1p[a][g]*dypqg_(cz1[a])))*fx2[j] * abint::abw[a]; // Drell-Yan
								th1stF += (-lx1*((fx1p[a][g]-fx1[5]*cz1[a])*3/(1.-cz1[a])+fx1p[a][g]*pggreg_(cz1[a]))
								)*fx2[5] * abint::abw[a];

								// Loop over quark flavours
								for (int pidq = 1; pidq <= MAXNF; pidq++) {

									// cout << "\npidq = "<< pidq << endl;
									// cout << "f(j) = "<< fx1p[0][MAXNF+pidq] << endl;
									// cout << "f(-j) = "<< fx1p[0][MAXNF-pidq] << endl;

									//H1st non delta terms
									// th1st += (fx1p[a][i]*cqq_(cz1[a])+fx1p[a][g]*cqg_(cz1[a]))*(-lx1)*fx2[j] * abint::abw[a]; // Drell-Yan
									th1st += (fx1p[a][MAXNF+pidq]+fx1p[a][MAXNF-pidq])*cgq_(cz1[a])*(-lx1)*fx2[5] * abint::abw[a]; // Honly

									// H1st muf dependence: Pgq, second leg
									th1stF += (-lx1*((fx1p[a][MAXNF+pidq]+fx1p[a][MAXNF-pidq])*pgq_(cz1[a])))*fx2[5] * abint::abw[a]; // FonlyPgq

								}

							}

							th1stF += -3. * d0intx1 * tdelta;

							//beta loop (second leg)
							for (int b = 0; b < abint::abdim; b++) {

								if (cz2[b] >= 1) continue;

								//H1st muf dependence, gammaqq and gammaqg:
								// th1stF += (-lx2*((fx2p[b][j]-fx2[j]*cz2[b])*pqq_(cz2[b])+fx2p[b][g]*dypqg_(cz2[b])))*fx1[i] * abint::abw[b]; // Drell-Yan
								th1stF += (-lx2*((fx2p[b][g]-fx2[5]*cz2[b])*3/(1.-cz2[b])+fx2p[b][g]*pggreg_(cz2[b]))
								)*fx1[5] * abint::abw[b];

								// Loop over quark flavours
								for (int pidq = 1; pidq <= MAXNF; pidq++) {

									// cout << "\npidq = "<< pidq << endl;
									// cout << "f(j) = "<< fx1p[b][MAXNF+pidq] << endl;
									// cout << "f(-j) = "<< fx1p[b][MAXNF-pidq] << endl;

									//H1st non delta terms
									// th1st += (fx2p[b][j]*cqq_(cz2[b])+fx2p[b][g]*cqg_(cz2[b]))*(-lx2)*fx1[i] * abint::abw[b]; // Drell-Yan
									th1st += (fx2p[b][MAXNF+pidq]+fx2p[b][MAXNF-pidq])*cgq_(cz2[b])*(-lx2)*fx1[5] * abint::abw[b]; // Honly

									// H1st muf dependence: Pgq, second leg
									th1stF += (-lx2*((fx2p[b][MAXNF+pidq]+fx2p[b][MAXNF-pidq])*pgq_(cz2[b])))*fx1[5] * abint::abw[b]; // FonlyPgq

								}

							}

							th1stF += -3. * d0intx2 * tdelta;

							th1stF = th1stF + 2*resconst::beta0 * tdelta;

							// NLL sigma contributions https://arxiv.org/abs/hep-ph/0508068v1 eqs. (63), (64)

							sig12[sp] = - 0.5 * resconst::A1g * tdelta;

							sig11[sp] = - (resconst::B1g + resconst::A1g * LQ) * tdelta - th1stF;


							//                        cout << " sig11[sp] ; " << sig11[sp] << endl;

							if (opts.order == 1) continue;


							cout << " NNLO  ; "  << endl;

							//end NLO

							//NNLO terms..................................................................................................................................
							double tcga = 0;
							double tgamma2 = 0;
							double tgaga = 0;

							//alfa loop
							double diff1 = 0;
							double diff10 = 0;
							double diff1f = 0;
							double diffg1f = 0;
							double diffg10 = 0;
							double diffc1f = 0;
							double diffc10 = 0;
							for (int a = 0; a < abint::abdim; a++) {

								if (cz1[a] >= 1) continue;

								//(gamma+gamma)*(gamma+gamma) term

								//First part: one gamma for each leg: gluon channel
								diff10 += 3*(-lx1*(fx1p[a][g]-fx1[5]*cz1[a])/(1.-cz1[a])-d0intx1*fx1[5]) * abint::abw[a];
								diff10 += (resconst::beta0*fx1[5]-lx1*fx1p[a][g]*pggreg_(cz1[a])) * abint::abw[a];

								//(C+C)*(gamma+gamma) term: diagonal part - gamma first leg, C second
								diffg10 += (-lx1*((fx1p[a][g]-fx1[5]*cz1[a])*3./(1.-cz1[a])+pggreg_(cz1[a])*fx1p[a][g])
								+fx1[5]*(resconst::beta0-3.*d0intx1)) * abint::abw[a];
								diffc10 += resconst::C1ggdelta*fx1[5] * abint::abw[a];

								//Second part: gamma*gamma terms
								//Pij * Pjk = D1ijjk (log(1-z)/(1-z))_+ + D0ijjk/(1-z)_+
								//          + Pijjk(z) + Deltaijjk delta(1-z)
								//First leg
								diff1 += (-lx1*((fx1p[a][g]-fx1[5]*cz1[a])
								*(resconst::D0gggg*oz1[a]+resconst::D1gggg*lzoz1[a])
								+fx1p[a][g]*pggggreg_(cz1[a])+fx1p[a][g]*pgqqg_(cz1[a]))
								+(resconst::Deltagggg-resconst::D0gggg*d0intx1-resconst::D1gggg*d1intx1)*fx1[5]) * abint::abw[a];
								//End of (gamma+gamma)*(gamma+gamma) term

								//C*gamma: first leg (ignore delta term in Cgg: taken into account with th1stF)
								// tcga += (fx1p[a][i]*cqqpqq_(cz1[a])+fx1p[a][g]*(cqqpqg_(cz1[a])+cqgpgg_(cz1[a])))*(-lx1)*fx2[j] * abint::abw[a];
								tcga += cgqpqg_(cz1[a])*(-lx1)*fx1p[a][g]*fx2[5] * abint::abw[a];

								//P2qq non-singlet: 1/(1-z)_+
								// tgamma2 += 2./3.*resconst::Kappa*(-lx1*(fx1p[a][i]-fx1[i]*cz1[a])*oz1[a]-d0intx1*fx1[i])*fx2[j] * abint::abw[a];
								tgamma2 += (1.5*resconst::Kappa*(-lx1*(fx1p[a][g]-fx1[5]*cz1[a])*oz1[a]-d0intx1*fx1[5])
								-lx1*p2gg_(cz1[a])*fx1p[a][g])*fx2[5] * abint::abw[a];

								// Loop over quark flavours
								for (int pidq = 1; pidq <= MAXNF; pidq++) {

									// Include Pgqqq and Pgggq
									diff1 += -lx1*(pgqqq_(cz1[a])+pgggq_(cz1[a]))*(fx1p[a][MAXNF+pidq]+fx1p[a][MAXNF-pidq]) * abint::abw[a];
									diff1f += -lx1*pgq_(cz1[a])*(fx1p[a][MAXNF+pidq]+fx1p[a][MAXNF-pidq]) * abint::abw[a];

									// Gamma first leg
									diffg1f += (-lx1*(fx1p[a][MAXNF+pidq]-fx1p[a][MAXNF-pidq])*pgq_(cz1[a])) * abint::abw[a];

									// C first leg
									diffc1f += (-lx1*(fx1p[a][MAXNF+pidq]-fx1p[a][MAXNF-pidq])*cgq_(cz1[a])) * abint::abw[a];

									// C*gamma: first leg (ignore delta term in Cgg: taken into account with th1stF)
									tcga += cgqpqq_(cz1[a])*(-lx1)*(fx1p[a][MAXNF+pidq]-fx1p[a][MAXNF-pidq])*fx2[5] * abint::abw[a];

									// Start 2-loop AP - gamma2: gq channel
									tgamma2 += -lx1*p2gq_(cz1[a])*(fx1p[a][MAXNF+pidq]+fx1p[a][MAXNF-pidq])*fx2[5] * abint::abw[a];

								}

							}

							//beta loop
							double diff2 = 0;
							double diff20 = 0;
							double diff2f = 0;
							double diffg2f = 0;
							double diffg20 = 0;
							double diffc2f = 0;
							double diffc20 = 0;
							for (int b = 0; b < abint::abdim; b++) {

								if (cz2[b] >= 1) continue;

								//(gamma+gamma)*(gamma+gamma) term

								//First part: one gamma for each leg
								diff20 += 3.*(-lx2*(fx2p[b][g]-fx2[5]*cz2[b])/(1.-cz2[b])-d0intx2*fx2[5]) * abint::abw[b];
								diff20 += (resconst::beta0*fx2[5]-lx2*fx2p[b][g]*pggreg_(cz2[b])) * abint::abw[b];

								//(C+C)*(gamma+gamma) term: diagonal part - gamma second leg, C first
								diffg20 += (-lx2*((fx2p[b][g]-fx2[5]*cz2[b])*3./(1.-cz2[b])+pggreg_(cz2[b])*fx2p[b][g])
								+fx2[5]*(resconst::beta0-3.*d0intx2)) * abint::abw[b];
								diffc20 += resconst::C1ggdelta*fx2[5] * abint::abw[b];

								//Second part: gamma*gamma terms
								//Pij * Pjk = D1ijjk (log(1-z)/(1-z))_+ + D0ijjk/(1-z)_+
								//          + Pijjk(z) + Deltaijjk delta(1-z)
								//Second leg
								diff2 += (-lx2*((fx2p[b][g]-fx2[5]*cz2[b])
								*(resconst::D0gggg*oz2[b]+resconst::D1gggg*lzoz2[b])
								+fx2p[b][g]*pggggreg_(cz2[b])+fx2p[b][g]*pgqqg_(cz2[b]))
								+(resconst::Deltagggg-resconst::D0gggg*d0intx2-resconst::D1gggg*d1intx2)*fx2[5]) * abint::abw[b];
								//End of (gamma+gamma)*(gamma+gamma) term

								//C*gamma: second leg (ignore delta term in Cqq: taken into account with th1stF)
								// tcga += (fx2p[b][j]*cqqpqq_(cz2[b])+fx2p[b][g]*(cqqpqg_(cz2[b])+cqgpgg_(cz2[b])))*(-lx2)*fx1[i] * abint::abw[b];
								tcga += cgqpqg_(cz2[b])*(-lx2)*fx2p[b][g]*fx1[5] * abint::abw[b];

								//Start 2-loop AP
								// Gluon + pure singlet
								//f == gluon piece
								// tgamma2 += fx2p[b][g]*p2qg_(cz2[b])*(-lx2)*fx1[i] * abint::abw[b];

								//P2qq non-singlet: 1/(1-z)_+
								// tgamma2 += 2./3.*resconst::Kappa*(-lx2*(fx2p[b][j]-fx2[j]*cz2[b])*oz2[b]-d0intx2*fx2[j])*fx1[i] * abint::abw[b];
								tgamma2 += (1.5*resconst::Kappa*(-lx2*(fx2p[b][g]-fx2[5]*cz2[b])*oz2[b]-d0intx2*fx2[5])
								-lx2*p2gg_(cz2[b])*fx2p[b][g])*fx1[5] * abint::abw[b];

								// Loop over quark flavours
								for (int pidq = 1; pidq <= MAXNF; pidq++) {

									// Include Pgqqq and Pgggq
									diff2 += -lx2*(pgqqq_(cz2[b])+pgggq_(cz2[b]))*(fx2p[b][MAXNF+pidq]+fx2p[b][MAXNF-pidq]) * abint::abw[b];
									diff2f += -lx2*pgq_(cz2[b])*(fx2p[b][MAXNF+pidq]+fx2p[b][MAXNF-pidq]) * abint::abw[b];

									// Gamma second leg
									diffg2f += (-lx2*(fx2p[b][MAXNF+pidq]-fx2p[b][MAXNF-pidq])*pgq_(cz2[b])) * abint::abw[b];

									// Gamma second leg
									diffc2f += (-lx2*(fx2p[b][MAXNF+pidq]-fx2p[b][MAXNF-pidq])*cgq_(cz2[b])) * abint::abw[b];

									// C*gamma: first leg (ignore delta term in Cgg: taken into account with th1stF)
									tcga += cgqpqq_(cz2[b])*(-lx2)*(fx2p[b][MAXNF+pidq]-fx2p[b][MAXNF-pidq])*fx1[5] * abint::abw[b];

									// gamma2: gq channel
									tgamma2 += -lx2*p2gq_(cz2[b])*(fx2p[b][MAXNF+pidq]+fx2p[b][MAXNF-pidq])*fx1[5] * abint::abw[b];

								}

							}

							tgaga += diff1*fx2[5];
							tgaga += diff2*fx1[5];
							tgaga=tgaga+2*(diff10*diff20+diff1f*diff2f+diff10*diff2f+diff1f*diff20);

							tcga += (diffc10*diffg20+diffc1f*diffg2f+diffc10*diffg2f+diffc1f*diffg20);
							tcga += (diffg10*diffc20+diffg1f*diffc2f+diffg10*diffc2f+diffg1f*diffc20);
							// tcga += (diffc10*diffg20+diffc1f*diffg2f);
							// tcga += (diffg10*diffc20+diffg1f*diffc2f);
							tcga = 0.;

							// NNLL sigma contributions - https://arxiv.org/abs/hep-ph/0508068v1 eqs. (66), (67), (68), (69)

							sig24[sp] = pow(resconst::A1g, 2) / 8 * tdelta;
							// sig24[sp] = 0.;

							sig23[sp] = -resconst::beta0*resconst::A1g/3.*tdelta-0.5*resconst::A1g*sig11[sp];
							// sig23[sp] = 0.;

							sig22[sp] = 0.5*(resconst::beta0*resconst::A1g*(LR-LQ)-resconst::A2g)*tdelta
							-0.5*resconst::A1g*(th1st+(LF-LQ)*th1stF)
							-0.5*(resconst::B1g-resconst::beta0)*sig11[sp]
							+0.5*resconst::B1g*th1stF
							+0.5*tgaga
							// Add mur dependence from H1st
							+resconst::beta0*resconst::A1g*LR*tdelta
							// New Q scale dependence
							+0.5*resconst::A1g*LQ * ( - sig11[sp] + 2*th1stF
							+tdelta * ( - resconst::beta0 + resconst::B1g + 0.5*resconst::A1g*LQ ) );
							// sig22[sp] = 0.;

							sig21[sp] = -resconst::beta0*LR*sig11[sp]-resconst::B1g*(th1st+LF*th1stF)
							-LF*tgaga-resconst::B2g*tdelta+resconst::beta0*th1st-tcga-tgamma2
							// Include missing delta term from C*gamma
							-resconst::C1ggdelta*th1stF
							// // Include missng term from contact term in 2 loop AP
							-2.*resconst::Delta2gg*tdelta;
							// Include mur dependence from H1st
							+2*resconst::beta0*LR*(resconst::B1g*tdelta+th1stF)
							// New Q scale dependence
							+ LQ * ( sig11[sp]*(resconst::beta0-resconst::B1g-0.5*resconst::A1g*LQ)
							-resconst::A2g*tdelta
							+resconst::A1g*(2.*resconst::beta0*LR*tdelta-th1st+(LQ-LF)*th1stF)
							+resconst::B1g*th1stF+tgaga );
							// sig21[sp] = 0.;

						}

						double xmsq = 0.;

						for (int sp = 0; sp < mesq::totpch; sp++)
						{

							//			  cout << "loopc  : "  << endl;
							//as/pi factor
							double sig1 = (sig12[sp]*qtint::LL2_mesqij[sp]+sig11[sp]*qtint::LL1_mesqij[sp])*asopi;
							xmsq += -sig1;


							//				cout << " qtint::LL2_mesqij[sp] : " << qtint::LL2_mesqij[sp] << endl;
							//	cout << " qtint::LL1_mesqij[sp] : " << qtint::LL1_mesqij[sp] << endl;

							if (opts.order == 1) continue;

							//(as/pi)^2 factor
							double sig2 = (sig24[sp]*qtint::LL4_mesqij[sp]+sig23[sp]*qtint::LL3_mesqij[sp]+sig22[sp]*qtint::LL2_mesqij[sp]+sig21[sp]*qtint::LL1_mesqij[sp])*pow(asopi,2);

							// sum O(as) and O(as^2) contributions
							xmsq += -sig2;
						}

						double shad = pow(opts.sroot,2);

						//double fbGeV2=0.38937966e12;
						//double fac = M_PI*6./fbGeV2*(2*m2);
						//double flux = fbGeV2/(2.*x1*x2*shad);
						//double ps = 1./shad;
						//double norm = 1./16./M_PI / 2. / M_PI;
						//xmsq = xmsq * ps * fac * flux * norm;

						//                        cout << " xmsq  " << xmsq << endl;

						if(opts.calc_type == "OLD")
						{
							xmsq = xmsq/shad * 3./8. /2./M_PI;
							//xmsq = xmsq/shad * 3./8. /2./M_PI/2./M_PI;
						}
						else
						{
							xmsq = xmsq/shad/16.;
						}
						xmsq = xmsq * m2; //Is this a Jacobian for dsigma/dq2 -> dsigma/dq?
						xmsq = xmsq * 2*M_PI; //phiV integration
						if (isnan_ofast(xmsq))
							cout << m << " " << y << " " << costh << "  " << xmsq << endl;

				//if (swtch < 0.01) return 0.;// do not apply this cut to avoid discontinuities. Instead the phase space is limited to qt and m switching limits
				xmsq = xmsq;//*swtch; //switching function is inside qtint

				//xmsq = xmsq*jac;//jacobian for the change of variable qt=qtp/sqrt(1-qtp^2/m^2)

				f[npdf] = xmsq;

			} //end loop on pdf

			return;

	}
	//Drell-Yan process
	else{

		for (int npdf = 0; npdf < opts.totpdf; npdf++)
			f[npdf] = 0.;

		/*
		//limit the integration to qt = m
		double jac = 1;
		double qtp = qt;
		if (qtp >= m*0.999)
			return 0.;
		else
			qt = qtp/sqrt(1-pow(qtp/m,2));
		jac = pow(m/sqrt(m*m-qtp*qtp),3);
		*/

		///////////////////////////////////////////////////////
		double m2 = m*m;
		double qt2 = qt*qt;
		double exppy = exp(y);
		double expmy = 1./exppy;
		double tau = sqrt(m2/pow(opts.sroot,2));

		//amplitudes are set in qtint::calc(), no need to calculate them
		/*
		mesq::setpropagators(m);
		double cthmom0, cthmom1, cthmom2;
		if (mode == 0)
			{
			//apply here lepton cuts
			cthmom0 = 1.;
			cthmom1 = costh;
			cthmom2 = pow(costh,2);
			}
		else if (mode == 1)
			{
			//phasespace::set_mqtyphi(m, 0., y);
			phasespace::set_qt(0.);
			omegaintegr::genV4p();
			omegaintegr::cthmoments(cthmom0,cthmom1,cthmom2);
			}

		mesq::setmesq(cthmom0, cthmom1, cthmom2);
		//cout << cthmom0 << "  " << cthmom1 << "  " << cthmom2 << endl;
		*/

		//calculate Bjorken x1 x2
		double x1 = tau*exppy;
		double x2 = tau*expmy;

		if (x1 >= 1 || x2 >= 1)
			return;

		//kinematical limits on y
		double ylim = fabs(log(tau));
		double ay = fabs(y);
		if (fabs(y) > ylim)
			return;
		if (fabs(*(long*)& ay - *(long*)& ylim) < 2 ) //check also equality
			return;

		//Set scales
		scales::set(m);
		scales::mcfm();
		scales::dyres(m);
		double muf = scales::fac;
		double mur = scales::ren;

		/*
		//Set factorization scale
		double muf, mur;
		if (opts.dynamicscale)
			{
			muf = m*opts.kmufac;
			mur = m*opts.kmuren;
			double mur2 = mur*mur;
			scaleset_(mur2); //set renormalization and factorization scales, and calculate ason2pi and ason4pi
			}
		else
			{
			muf = opts.rmass*opts.kmufac;
			mur = opts.rmass*opts.kmuren;
			}
		*/

		//a-parameter of the resummation scale

		//dynamic case
		//if (opts.fmures > 0)
		//a_param_.a_param_ = 1./opts.kmures;

		//for fixed resummation scale need to recompute a_param
		//else
		//a_param_.a_param_ = m/scales::res;

		//PDFs
		double fx1[2*MAXNF+1],fx2[2*MAXNF+1];
		fdist_(opts.ih1,x1,muf,fx1);
		fdist_(opts.ih2,x2,muf,fx2);
		pdf::fdist(1,x1,muf,fx1);
		pdf::fdist(2,x2,muf,fx2);
		//////////////////////////////////////////////////////////

		double LR, LF, LQ;
		if (order >= 2)
			LR = log(m2/pow(mur,2));
		LF = log(m2/pow(muf,2));
		LQ = 2.*log(a_param_.a_param_);

		double LQ2 = pow(LQ,2);
		double LQ3 = pow(LQ,3);
		double LQ4 = pow(LQ,4);
		double LQ5 = pow(LQ,5);

		double LR2 = pow(LR,2);

		double LF2 = pow(LF,2);

		/*
		//Calculate large logs LL1-LL4
		double LL1,LL2,LL3,LL4;
		if (mode == 0 || mode == 1)
			{
			//xmio is used in besselkfast for Itilde
			xmio_.xmio_ = sqrt(qt2/(m2/pow(a_param_.a_param_,2)));

			//LL1,LL2,LL3,LL4: large log (squared) corresponding to eq. (136)
			//In this way normalization is fixed to dsigma/dqt2

			//depends on xmio, which depends on qt2 and m2
			LL1 = itilde_(one)/pow(m2,2)*pow(a_param_.a_param_,2);
			LL2 = itilde_(two)/pow(m2,2)*pow(a_param_.a_param_,2);
			LL3 = itilde_(three)/pow(m2,2)*pow(a_param_.a_param_,2);
			LL4 = itilde_(four)/pow(m2,2)*pow(a_param_.a_param_,2);
			}
		*/

		//Scaled momentum fractions
		double lx1, lx2;
		double pqqintx1, pqqintx2;
		double d0intx1, d0intx2;
		double d1intx1, d1intx2;
		if (order >= 1)
			{
			lx1 = log(x1);
			lx2 = log(x2);
			pqqintx1 = pqqint_(x1);
			pqqintx2 = pqqint_(x2);
			}
		if (order >= 2)
			{
			d0intx1 = d0int_(x1);
			d0intx2 = d0int_(x2);
			d1intx1 = d1int_(x1);
			d1intx2 = d1int_(x2);
			}

		//Preliminary loop for caching
		double cz1[abint::abdim];
		double cz2[abint::abdim];
		double oz1[abint::abdim];
		double oz2[abint::abdim];
		// --> cache the logs, and pass them as parameters to the CxxPxx, Pxxxx functions
		double lz1[abint::abdim];
		double lz2[abint::abdim];
		double l1z1[abint::abdim];
		double l1z2[abint::abdim];
		double lzoz1[abint::abdim];
		double lzoz2[abint::abdim];
		if (order >= 1)
			for (int ab = 0; ab < abint::abdim; ab++)
			{
			cz1[ab] = pow(x1,abint::abx[ab]);
			cz2[ab] = pow(x2,abint::abx[ab]);
			oz1[ab] = 1./(1.-cz1[ab]);
			oz2[ab] = 1./(1.-cz2[ab]);
			lz1[ab] = log(cz1[ab]);
			lz2[ab] = log(cz2[ab]);
			l1z1[ab] = log(1.-cz1[ab]);
			l1z2[ab] = log(1.-cz2[ab]);
			lzoz1[ab] = log(1.-cz1[ab])/(1.-cz1[ab]);
			lzoz2[ab] = log(1.-cz2[ab])/(1.-cz2[ab]);
			}

		//tc2ga1
		double p1qqc2qqhsreg1[abint::abdim];
		double p1qqc2qqhsd01[abint::abdim];
		double p1qqc2qqhsd11[abint::abdim];
		double p1qqc2qqhsdelta1[abint::abdim];
		double p1qqc2qqb1[abint::abdim];
		double p1qqc2qqp1[abint::abdim];
		double p1gqc2qg1[abint::abdim];
		double p1qgc2qq1[abint::abdim];
		double p1qgc2qqb1[abint::abdim];
		double p1qgc2qqp1[abint::abdim];
		double c2qqreghschm1[abint::abdim];
		double c2qqb1[abint::abdim];
		double c2qqp1[abint::abdim];
		double c2qghschm1[abint::abdim];
		double p1ggc2qg1[abint::abdim];

		double p1qqc2qqhsreg2[abint::abdim];
		double p1qqc2qqhsd02[abint::abdim];
		double p1qqc2qqhsd12[abint::abdim];
		double p1qqc2qqhsdelta2[abint::abdim];
		double p1qqc2qqb2[abint::abdim];
		double p1qqc2qqp2[abint::abdim];
		double p1gqc2qg2[abint::abdim];
		double p1qgc2qq2[abint::abdim];
		double p1qgc2qqb2[abint::abdim];
		double p1qgc2qqp2[abint::abdim];
		double c2qqreghschm2[abint::abdim];
		double c2qqb2[abint::abdim];
		double c2qqp2[abint::abdim];
		double c2qghschm2[abint::abdim];
		double p1ggc2qg2[abint::abdim];

		//tgamma3
		double p3qqv1[abint::abdim];
		double p3qqvd01[abint::abdim];
		double p3qqvdelta1[abint::abdim];
		double p3qqbv1[abint::abdim];
		double p3sqq1[abint::abdim];
		double p3sqqb1[abint::abdim];
		double p3qg1[abint::abdim];
		double p3qqv2[abint::abdim];
		double p3qqvd02[abint::abdim];
		double p3qqvdelta2[abint::abdim];
		double p3qqbv2[abint::abdim];
		double p3sqq2[abint::abdim];
		double p3sqqb2[abint::abdim];
		double p3qg2[abint::abdim];
		if (order >= 3)
			for (int ab = 0; ab < abint::abdim; ab++)
			{
			//if (cz1[ab] >= 1) continue;
			cacheuv_(cz1[ab]);
			p1qqc2qqhsreg1[ab]   =  p1qqc2qqhsreg_(cz1[ab]);
			p1qqc2qqhsd01[ab]    =  p1qqc2qqhsd0_(cz1[ab]);
			p1qqc2qqhsd11[ab]    =  p1qqc2qqhsd1_(cz1[ab]);
			p1qqc2qqhsdelta1[ab] =  p1qqc2qqhsdelta_(cz1[ab]);
			p1qqc2qqb1[ab]	     =  p1qqc2qqb_(cz1[ab]);
			p1qqc2qqp1[ab]	     =  p1qqc2qqp_(cz1[ab]);
			p1gqc2qg1[ab]	     =  p1gqc2qg_(cz1[ab]);
			p1qgc2qq1[ab]	     =  p1qgc2qq_(cz1[ab]);
			p1qgc2qqb1[ab]	     =  p1qgc2qqb_(cz1[ab]);
			p1qgc2qqp1[ab]	     =  p1qgc2qqp_(cz1[ab]);
			c2qqreghschm1[ab]    =  c2qqreghschm_(cz1[ab]);
			c2qqb1[ab]	     =  c2qqb_(cz1[ab]);
			c2qqp1[ab]	     =  c2qqp_(cz1[ab]);
			c2qghschm1[ab]       =  c2qghschm_(cz1[ab]);
			p1ggc2qg1[ab]        =  p1ggc2qg_(cz1[ab]);
			p3qqv1[ab]      = p3qqv_(cz1[ab]);
			p3qqvd01[ab]	= p3qqvd0_(cz1[ab]);
			p3qqvdelta1[ab]	= p3qqvdelta_(cz1[ab]);
			p3qqbv1[ab]	= p3qqbv_(cz1[ab]);
			p3sqq1[ab]	= p3sqq_(cz1[ab]);
			p3sqqb1[ab]	= p3sqqb_(cz1[ab]);
			p3qg1[ab]	= p3qg_(cz1[ab]);
			}
		if (order >= 3)
			for (int ab = 0; ab < abint::abdim; ab++)
			{
			//if (cz2[ab] >= 1) continue;
			cacheuv_(cz2[ab]);
			p1qqc2qqhsreg2[ab]   =  p1qqc2qqhsreg_(cz2[ab]);
			p1qqc2qqhsd02[ab]    =  p1qqc2qqhsd0_(cz2[ab]);
			p1qqc2qqhsd12[ab]    =  p1qqc2qqhsd1_(cz2[ab]);
			p1qqc2qqhsdelta2[ab] =  p1qqc2qqhsdelta_(cz2[ab]);
			p1qqc2qqb2[ab]	     =  p1qqc2qqb_(cz2[ab]);
			p1qqc2qqp2[ab]	     =  p1qqc2qqp_(cz2[ab]);
			p1gqc2qg2[ab]	     =  p1gqc2qg_(cz2[ab]);
			p1qgc2qq2[ab]	     =  p1qgc2qq_(cz2[ab]);
			p1qgc2qqb2[ab]	     =  p1qgc2qqb_(cz2[ab]);
			p1qgc2qqp2[ab]	     =  p1qgc2qqp_(cz2[ab]);
			c2qqreghschm2[ab]    =  c2qqreghschm_(cz2[ab]);
			c2qqb2[ab]	     =  c2qqb_(cz2[ab]);
			c2qqp2[ab]	     =  c2qqp_(cz2[ab]);
			c2qghschm2[ab]       =  c2qghschm_(cz2[ab]);
			p1ggc2qg2[ab]        =  p1ggc2qg_(cz2[ab]);
			p3qqv2[ab]      = p3qqv_(cz2[ab]);
			p3qqvd02[ab]	= p3qqvd0_(cz2[ab]);
			p3qqvdelta2[ab]	= p3qqvdelta_(cz2[ab]);
			p3qqbv2[ab]	= p3qqbv_(cz2[ab]);
			p3sqq2[ab]	= p3sqq_(cz2[ab]);
			p3sqqb2[ab]	= p3sqqb_(cz2[ab]);
			p3qg2[ab]	= p3qg_(cz2[ab]);
			}

		// skip PDF loop in the preconditioning phase
		int maxpdf=0;
		if (dofill_.doFill_ != 0) maxpdf = opts.totpdf;

		// start PDF loop
		for (int npdf = 0; npdf < maxpdf; npdf++)
			{
			dysetpdf_(npdf); // Inicializa los members de las PDFs y su alpha_s asociada en src/pdf.C
			// Si hay un solo member (e.g npfd=0) calcula el miembro central solo

			//cache scaled PDFs
			double fx1p[abint::abdim][2*MAXNF+1],fx2p[abint::abdim][2*MAXNF+1];
			double f1Ap[abint::abdim],f2Ap[abint::abdim];
			if (order >= 1)
			for (int ab = 0; ab < abint::abdim; ab++)
			{
				double fx1temp[2*MAXNF+1],fx2temp[2*MAXNF+1];
				double xx1 = pow(x1,(1-abint::abx[ab]));
				double xx2 = pow(x2,(1-abint::abx[ab]));
				fdist_(opts.ih1,xx1,muf,fx1temp);
				fdist_(opts.ih2,xx2,muf,fx2temp);
				pdf::fdist(1,xx1,muf,fx1temp);
				pdf::fdist(2,xx2,muf,fx2temp);
				memcpy(fx1p[ab], fx1temp, (2*MAXNF+1)*sizeof(double));
				memcpy(fx2p[ab], fx2temp, (2*MAXNF+1)*sizeof(double));
				pdf::fdistphot(1,xx1,muf,f1Ap[ab]);
				pdf::fdistphot(2,xx2,muf,f2Ap[ab]);
			}

			double sumfx1p[abint::abdim];
			double sumfx2p[abint::abdim];
			if (order >= 2)
			for (int ab = 0; ab < abint::abdim; ab++)
			{
				sumfx1p[ab] = 0.;
				sumfx2p[ab] = 0.;
				for (int f = 0; f < MAXNF; f++)
				{
				sumfx1p[ab] += fx1p[ab][f]+fx1p[ab][parton::charge_conj(parton::pdgid(f))];
				sumfx2p[ab] += fx2p[ab][f]+fx2p[ab][parton::charge_conj(parton::pdgid(f))];
				}
			}

			double fx1pafx1pmanoti[2*MAXNF+1][abint::abdim];
			double fx2pafx2pmanotj[2*MAXNF+1][abint::abdim];
			memset(fx1pafx1pmanoti, 0, sizeof(fx1pafx1pmanoti));
			memset(fx2pafx2pmanotj, 0, sizeof(fx2pafx2pmanotj));
			if (order >= 2)
			for (int ab = 0; ab < abint::abdim; ab++)
			for (int f = 0; f < MAXNF; f++)
				{
				fx1pafx1pmanoti[f][ab] = fx1pafx1pmanoti[parton::charge_conj(parton::pdgid(f))][ab] = sumfx1p[ab] - (fx1p[ab][f]+fx1p[ab][parton::charge_conj(parton::pdgid(f))]);
				fx2pafx2pmanotj[f][ab] = fx2pafx2pmanotj[parton::charge_conj(parton::pdgid(f))][ab] = sumfx2p[ab] - (fx2p[ab][f]+fx2p[ab][parton::charge_conj(parton::pdgid(f))]);
				}

			double sumposfx1p[abint::abdim];
			double sumnegfx1p[abint::abdim];
			double sumposfx2p[abint::abdim];
			double sumnegfx2p[abint::abdim];
			memset(sumposfx1p, 0, sizeof(sumposfx1p));
			memset(sumnegfx1p, 0, sizeof(sumnegfx1p));
			memset(sumposfx2p, 0, sizeof(sumposfx2p));
			memset(sumnegfx2p, 0, sizeof(sumnegfx2p));
			if (order >= 3)
			for (int ab = 0; ab < abint::abdim; ab++)
			{
				for (int f = 0; f < MAXNF; f++)
				{
				sumposfx1p[ab] += fx1p[ab][f];
				sumnegfx1p[ab] += fx1p[ab][parton::charge_conj(parton::pdgid(f))];
				sumposfx2p[ab] += fx2p[ab][f];
				sumnegfx2p[ab] += fx2p[ab][parton::charge_conj(parton::pdgid(f))];
				}
			}
			// Start calculation
			double asopi;
			if (opts.alphaslha)
			asopi = pdf::alphas(scales::ren)/M_PI;
			else
			asopi = pdf::rgktalphas(scales::ren)/M_PI;

			//Use alpha(mur) QED for initial state photons
			double facqedis = (alphaqed::calc(pow(scales::ren,2),opts.qedorder+1)*4.*M_PI)/resconst::a0qed;

			//loop on born subprocesses, i.e. born incoming partons ij
			//double lumi[mesq::totpch];
			double sig11[mesq::totpch];
			double sig12[mesq::totpch];
			double sig21[mesq::totpch];
			double sig22[mesq::totpch];
			double sig23[mesq::totpch];
			double sig24[mesq::totpch];
			double sig31[mesq::totpch];
			double sig32[mesq::totpch];
			double sig33[mesq::totpch];
			double sig34[mesq::totpch];
			double sig35[mesq::totpch];
			double sig36[mesq::totpch];
			//QED
			double sig11_qed[mesq::totpch] = {0.};
			double sig12_qed[mesq::totpch] = {0.};
			for (int sp = 0; sp < mesq::totpch; sp++)
			{
			//simplify notation
			//double bornmesqij = real(mesq::mesqij[sp]); //born level amplitudes
			parton::pdgid i = mesq::pid1[sp];         //parton 1
			parton::pdgid j = mesq::pid2[sp];         //parton 2
			parton::pdgid g = parton::G;              //gluon
			parton::pdgid im = parton::charge_conj(i);
			parton::pdgid jm = parton::charge_conj(j);

			//LO term (there is no counterterm at LO...)
			//Simplest term without convolutions
			double tdelta = fx1[i]*fx2[j];

			//NLO terms
			double th1st = 0; //this piece is used only at NNLO?
			double th1stF = 0;
			double th1stQ = 0; //this piece is used only at NNLO

			double th1stF_qed = 0;

			//H1st delta term
			th1st += 2*resconst::C1qqdelta*tdelta;

			//add resummation scale dependence
			th1stQ += -(resconst::B1q+resconst::A1q/2.*LQ)*LQ*tdelta; //this piece is used only at NNLO

			//alfa loop (first leg)
			for (int a = 0; a < abint::abdim; a++)
				{
				if (cz1[a] >= 1) continue;

				//H1st non delta terms
				th1st += (fx1p[a][i]*cqq_(cz1[a])+fx1p[a][g]*cqg_(cz1[a]))*(-lx1)*fx2[j] * abint::abw[a];

				//H1st muf dependence, gammaqq and gammaqg:
				th1stF     += (-lx1*((fx1p[a][i]-fx1[i]*cz1[a])*pqq_(cz1[a])+fx1p[a][g]*dypqg_(cz1[a])))*fx2[j] * abint::abw[a];
				th1stF_qed += (-lx1*((fx1p[a][i]-fx1[i]*cz1[a])*3./4.*pqq_(cz1[a])))*fx2[j] * abint::abw[a]*(parton::chsq(i)+parton::chsq(j))/2.;
				th1stF_qed += (-lx1*( +f1Ap[a]*6.*dypqg_(cz1[a])   ))*fx2[j] * abint::abw[a]*parton::chsq(j) * facqedis; //photon PDF
				}
			th1stF     += -pqqintx1*tdelta;
			th1stF_qed += -3./4.*pqqintx1*tdelta*(parton::chsq(i)+parton::chsq(j))/2.;

			//beta loop (second leg)
			for (int b = 0; b < abint::abdim; b++)
				{
				if (cz2[b] >= 1) continue;

				//H1st non delta terms
				th1st += (fx2p[b][j]*cqq_(cz2[b])+fx2p[b][g]*cqg_(cz2[b]))*(-lx2)*fx1[i] * abint::abw[b];

				//H1st muf dependence, gammaqq and gammaqg:
				th1stF     += (-lx2*((fx2p[b][j]-fx2[j]*cz2[b])*pqq_(cz2[b])+fx2p[b][g]*dypqg_(cz2[b])))*fx1[i] * abint::abw[b];
				th1stF_qed += (-lx2*((fx2p[b][j]-fx2[j]*cz2[b])*3./4.*pqq_(cz2[b])))*fx1[i] * abint::abw[b]*(parton::chsq(i)+parton::chsq(j))/2.;
				th1stF_qed += (-lx2*(f2Ap[b]*6.*dypqg_(cz2[b]) ))*fx1[i] * abint::abw[b]*parton::chsq(i) * facqedis; //photon PDF
				}
			th1stF     += -pqqintx2*tdelta;
			th1stF_qed += -3./4.*pqqintx2*tdelta*(parton::chsq(i)+parton::chsq(j))/2.;


			sig12[sp] = -0.5*resconst::A1q*tdelta;
			sig11[sp] = -(resconst::B1q+resconst::A1q*LQ)*tdelta - th1stF;

			//QED
			if (opts.qed && opts.qedorder_ct >= 1)
				{
				sig12_qed[sp] = -0.5*resconst::A1qp*tdelta*(parton::chsq(i)+parton::chsq(j))/2.;
				sig11_qed[sp] = (-(resconst::B1qp+resconst::A1qp*LQ)*tdelta)*(parton::chsq(i)+parton::chsq(j))/2.;
				sig11_qed[sp] += (-th1stF_qed); //include photon PDF
				}

			if (order == 1) continue;
			//end NLO

			//NNLO terms
			double tcga = 0;
			double tgamma2 = 0;
			double tgaga = 0;

			//alfa loop
			double diffg1f = 0;
			double diffg10 = 0;
			double diff1 = 0;
			double diffc1f = 0;
			double diffc10 = 0;

			for (int a = 0; a < abint::abdim; a++)
				{
				if (cz1[a] >= 1) continue;

				//(gamma+gamma)*(gamma+gamma) term

				//First part: one gamma for each leg
				diffg1f += (-lx1*(fx1p[a][i]-fx1[i]*cz1[a])*pqq_(cz1[a]) - pqqintx1*fx1[i]) * abint::abw[a];
				diffg10 += -lx1*fx1p[a][g]*dypqg_(cz1[a]) * abint::abw[a];

				//Second part: gamma*gamma terms
				//Pij * Pjk = D1ijjk (log(1-z)/(1-z))_+ + D0ijjk/(1-z)_+
				//          + Pijjk(z) + Deltaijjk delta(1-z)
				//First leg
				diff1 += (-lx1*((fx1p[a][i]-fx1[i]*cz1[a])
						*(resconst::D0qqqq*oz1[a]+resconst::D1qqqq*lzoz1[a])
						+fx1p[a][i]*pqqqq_(cz1[a])+fx1p[a][g]*(pqqqg_(cz1[a])+pqggg_(cz1[a])))
					+(resconst::Deltaqqqq-resconst::D0qqqq*d0intx1-resconst::D1qqqq*d1intx1)
					*fx1[i]) * abint::abw[a];

				//Include Pqggq
				diff1 += -lx1*sumfx1p[a]*pqggq_(cz1[a]) * abint::abw[a];
				//End of (gamma+gamma)*(gamma+gamma) term

				//Start  (C+C)*(gamma+gamma) term
				//C first leg, gamma second leg
				diffc1f += (-lx1*fx1p[a][i]*cqq_(cz1[a])+resconst::C1qqdelta*fx1[i]) * abint::abw[a];
				diffc10 += -lx1*fx1p[a][g]*cqg_(cz1[a]) * abint::abw[a];

				//C*gamma: first leg (ignore delta term in Cqq: taken into account with th1stF)
				tcga += (fx1p[a][i]*cqqpqq_(cz1[a])+fx1p[a][g]*(cqqpqg_(cz1[a])+cqgpgg_(cz1[a])))*(-lx1)*fx2[j] * abint::abw[a];

				//Add Cqg*Pgq contribution
				tcga += sumfx1p[a]*cqgpgq_(cz1[a])*(-lx1)*fx2[j] * abint::abw[a];

				//Start 2-loop AP
				// Gluon + pure singlet
				//f == gluon piece
				tgamma2 += fx1p[a][g]*p2qg_(cz1[a])*(-lx1)*fx2[j] * abint::abw[a];

				//f != gluon piece
				tgamma2 += sumfx1p[a]*p2qqs_(cz1[a])*(-lx1)*fx2[j] * abint::abw[a];

				//P2qq non-singlet: regular part
				tgamma2 += fx1p[a][i]*p2qqv_(cz1[a])*(-lx1)*fx2[j] * abint::abw[a];

				//P2qq non-singlet: 1/(1-z)_+
				tgamma2 += 2./3.*resconst::Kappa*(-lx1*(fx1p[a][i]-fx1[i]*cz1[a])*oz1[a]-d0intx1*fx1[i])*fx2[j] * abint::abw[a];

				//P2qqb non singlet
				tgamma2 += fx1p[a][im]*p2qqbv_(cz1[a])*(-lx1)*fx2[j] * abint::abw[a];
				}

			//beta loop
			double diffg2f = 0;
			double diffg20 = 0;
			double diff2 = 0;
			double diffc2f = 0;
			double diffc20 = 0;
			for (int b = 0; b < abint::abdim; b++)
				{
				if (cz2[b] >= 1) continue;

				//(gamma+gamma)*(gamma+gamma) term

				//First part: one gamma for each leg
				diffg2f += (-lx2*(fx2p[b][j]-fx2[j]*cz2[b])*pqq_(cz2[b]) - pqqintx2*fx2[j]) * abint::abw[b];
				diffg20 += -lx2*fx2p[b][g]*dypqg_(cz2[b]) * abint::abw[b];

				//Second part: gamma*gamma terms
				//Pij * Pjk = D1ijjk (log(1-z)/(1-z))_+ + D0ijjk/(1-z)_+
				//          + Pijjk(z) + Deltaijjk delta(1-z)
				//Second leg
				diff2 += (-lx2*((fx2p[b][j]-fx2[j]*cz2[b])
						*(resconst::D0qqqq*oz2[b]+resconst::D1qqqq*lzoz2[b])
						+fx2p[b][j]*pqqqq_(cz2[b])+fx2p[b][g]*(pqqqg_(cz2[b])+pqggg_(cz2[b])))
					+(resconst::Deltaqqqq-resconst::D0qqqq*d0intx2-resconst::D1qqqq*d1intx2)
					*fx2[j])* abint::abw[b];
				//Include Pqggq
				diff2 += -lx2*sumfx2p[b]*pqggq_(cz2[b]) * abint::abw[b];
				//End of (gamma+gamma)*(gamma+gamma) term

				//Start  (C+C)*(gamma+gamma) term
				//gamma first leg, C second leg
				diffc2f += (-lx2*fx2p[b][j]*cqq_(cz2[b])+resconst::C1qqdelta*fx2[j]) * abint::abw[b];
				diffc20 += -lx2*fx2p[b][g]*cqg_(cz2[b]) * abint::abw[b];

				//C*gamma: second leg (ignore delta term in Cqq: taken into account with th1stF)
				tcga += (fx2p[b][j]*cqqpqq_(cz2[b])+fx2p[b][g]*(cqqpqg_(cz2[b])+cqgpgg_(cz2[b])))*(-lx2)*fx1[i] * abint::abw[b];

				//Add Cqg*Pgq contribution
				tcga += sumfx2p[b]*cqgpgq_(cz2[b])*(-lx2)*fx1[i] * abint::abw[b];

				//Start 2-loop AP
				// Gluon + pure singlet
				//f == gluon piece
				tgamma2 += fx2p[b][g]*p2qg_(cz2[b])*(-lx2)*fx1[i] * abint::abw[b];

				//f != gluon piece
				tgamma2 += sumfx2p[b]*p2qqs_(cz2[b])*(-lx2)*fx1[i] * abint::abw[b];

				//P2qq non-singlet: regular part
				tgamma2 += fx2p[b][j]*p2qqv_(cz2[b])*(-lx2)*fx1[i] * abint::abw[b];

				//P2qq non-singlet: 1/(1-z)_+
				tgamma2 += 2./3.*resconst::Kappa*(-lx2*(fx2p[b][j]-fx2[j]*cz2[b])*oz2[b]-d0intx2*fx2[j])*fx1[i] * abint::abw[b];

				//P2qqb non singlet
				tgamma2 += fx2p[b][jm]*p2qqbv_(cz2[b])*(-lx2)*fx1[i] * abint::abw[b];
				}

			tgamma2 += +2.*tdelta*resconst::Delta2qq; // Delta2qq beelongs to tgamma2. tgamma2 nees to be defined in a closed way in order to be used at N3LO

			tgaga=tgaga+2*(diffg10*diffg20+diffg1f*diffg2f+diffg10*diffg2f+diffg1f*diffg20);
			tgaga += diff1*fx2[j];
			tgaga += diff2*fx1[i];

			tcga += (diffc10*diffg20+diffc1f*diffg2f+diffc10*diffg2f+diffc1f*diffg20);
			tcga += (diffg10*diffc20+diffg1f*diffc2f+diffg10*diffc2f+diffg1f*diffc20);

			sig24[sp] = pow(resconst::A1q,2) / 8. * tdelta;
			sig23[sp] = -resconst::beta0*resconst::A1q/3.*tdelta-0.5*resconst::A1q*sig11[sp];
			sig22[sp] = 0.5*(resconst::beta0*resconst::A1q*(LR-LQ)-resconst::A2q)*tdelta
				-0.5*resconst::A1q*(th1st+th1stQ+(LF-LQ)*th1stF)
				-0.5*(resconst::B1q+resconst::A1q*LQ-resconst::beta0)*sig11[sp]
				+0.5*(resconst::B1q+resconst::A1q*LQ)*th1stF
				+0.5*tgaga;
			sig21[sp] =
				-resconst::beta0*(LR-LQ)*sig11[sp]
				-(resconst::B1q+resconst::A1q*LQ)*(th1st+th1stQ+(LF-LQ)*th1stF)
				-(LF-LQ)*tgaga
				-(resconst::B2q+resconst::A2q*LQ)*tdelta
				+resconst::beta0*th1st
				+(resconst::B1q+0.5*resconst::A1q*LQ)*LQ*th1stF
				-tcga-tgamma2;

			//include missing delta term from C*gamma (no factor 2 here)
			sig21[sp] += -resconst::C1qqdelta*th1stF;
			//include missing term from contact term in 2 loop AP
			// sig21[sp] += -2*resconst::Delta2qq*tdelta; // This part belongs to tgamma2. The N3LO contribution uses the complete and correct tgamma2 14.10.2021

					//end NNLO
			if (order == 2) continue;

			// ----------------------------------------------------------------------------------------------------------------------------------------
				//NNNLO terms written in the Hard Scheme
			// ----------------------------------------------------------------------------------------------------------------------------------------
			double tgagaga = 0;
			double tcgah = 0;
			double tcgaga = 0;
			double tga1ga2 = 0;
			double tc2ga1 = 0;
			double tc1ga2 = 0;
			double tccga = 0;
			double tgamma3 = 0;

			double C1caC1cbarb = 0;
			double C1cbarbpiuC1ca = 0;
			double tH2st = 0;
			// tgagaga definitions
			double tgagaqL1 = 0;
			double tgagaqL2 = 0;
			double diffqq1  = 0;
			double diffqq2 = 0;
			double diffqg1 = 0;
			double diffqg2 = 0;
			double tgagagL1 = 0;
			double tgagagL2 = 0;
			diffg10 = 0;// Probabilmente sia possibile usare quella NNLO. Occhio con lo schema CSS
			diffg20 = 0;// Probabilmente sia possibile usare quella NNLO. Occhio con lo schema CSS
			diffg1f = 0;// Probabilmente sia possibile usare quella NNLO. Occhio con lo schema CSS
			diffg2f = 0;// Probabilmente sia possibile usare quella NNLO. Occhio con lo schema CSS
			// tcgaga definitions
			double tcgahqL1 = 0;
			double tcgahgL1 = 0;
			double tcgahqL2 = 0;
			double tcgahgL2 = 0;
			diffc1f = 0; //Si deve definire di nuovo perche' il CT[NNLO] e' scritto nello schema CSS
			diffc2f = 0; //Si deve definire di nuovo perche' il CT[NNLO] e' scritto nello schema CSS
			// tc2ga1 definitions
			diff1 = 0;
			diff2 = 0;
			double C2qqL1 = 0;
			double C2qqL2 = 0;
			double C2qgL1 = 0;
			double C2qgL2 = 0;
			// tc2ga1 definitions
			double gamma2qqL1 = 0;
			double gamma2qqL2 = 0;
			double gamma2qgL1 = 0;
			double gamma2qgL2 = 0;
			//Conversions:
			//IntFuncqT -> intfuncqt

			//z1 -> cz1[a]
			//z2 -> cz2[b]

			//fx1p(0) -> fx1p[a][g]
			//fx2p(0) -> fx2p[b][g]
			//fx1p(j) -> fx1p[a][i]
			//fx2p(k) -> fx2p[b][j]
			//fx10(j) -> fx1[i]
			//fx20(k) -> fx2[j]

			//fx1p(-j) -> fx1p[a][im]
			//fx2p(-k) -> fx2p[b][jm]

			//xx10 -> x1
			//xx20 -> x2

			//fx1plfx1pml -> sumfx1p[a]
			//fx2plfx2pml -> sumfx2p[b]

			//msqc(j,k) -> (remove)
			//fl* -> (remove)

			//
			// ----------------------------------------------------------------------------------------------------------------------------------------
			// tgagaga
			// ----------------------------------------------------------------------------------------------------------------------------------------
			//Leg 1
			for (int a = 0; a < abint::abdim; a++)
				{
				cacheuv_(cz1[a]);
				//if (cz1[a] >= 1) continue;
				// qqb Channel
				tgagaga  += intfuncqt(fx1p[a][i],fx1[i],x1,cz1[a],p1qqp1qqp1qqreg_(cz1[a]),p1qqp1qqp1qqd0_(cz1[a]),p1qqp1qqp1qqd1_(cz1[a]),p1qqp1qqp1qqd2_(cz1[a]),p1qqp1qqp1qqdelta_(cz1[a]))*fx2[j] *abint::abw[a]; //PqqPqqPqq
				tgagaga  += - 2.*lx1*(sumfx1p[a])*p1gqp1qqp1qg_(cz1[a])*fx2[j] *abint::abw[a]; //PqqPqgPgq
				// PqqPqgPgq The Contribution is the same as the previous one ->  we simply put a factor 2
				tgagaga  += - lx1*(sumfx1p[a])*p1qgp1ggp1gq_(cz1[a])*fx2[j] *abint::abw[a];		                                                                                                            //PqgPggPgq
				// qg Channel
				tgagaga += - lx1*fx1p[a][g]*p1qqp1qqp1qg_(cz1[a])*fx2[j]*abint::abw[a]; //PqqPqqPqg
				tgagaga += - lx1*fx1p[a][g]*p1qqp1qgp1gg_(cz1[a])*fx2[j]*abint::abw[a]; //PqqPqgPgg
				tgagaga += - lx1*fx1p[a][g]*p1qgp1gqp1qg_(cz1[a])*fx2[j]*abint::abw[a]*2.*NF;//PqgPgqPqg ! Candidate for Nf factor?!!!!
				tgagaga += - lx1*fx1p[a][g]*p1qgp1ggp1gg_(cz1[a])*fx2[j]*abint::abw[a]; // PqgPggPgg
				// Start contributions used for simple products contributions
				diffqq1  += (-lx1*(fx1p[a][i]-fx1[i]*cz1[a])*pqq_(cz1[a]) - pqqintx1*fx1[i])*abint::abw[a]; // porque no se puede usar la de NNLO? si se usa fuera del loop?
				tgagaqL1 += (-lx1*((fx1p[a][i]-fx1[i]*cz1[a]) *(D0qqqq/(1.-cz1[a])+D1qqqq*log(1.-cz1[a])/(1.-cz1[a]))
						+fx1p[a][i]*pqqqq_(cz1[a]))+(Deltaqqqq-D0qqqq*d0intx1-D1qqqq*d1intx1)*fx1[i])*abint::abw[a];
				tgagaqL1 += -lx1*(sumfx1p[a])*pqggq_(cz1[a])*abint::abw[a];
				diffg10  += -lx1*fx1p[a][g]*dypqg_(cz1[a])*abint::abw[a];
				tgagagL1 += -lx1*(fx1p[a][g]*(pqqqg_(cz1[a])+pqggg_(cz1[a])))*abint::abw[a];
				}
			//Leg 2
			for (int b = 0; b < abint::abdim; b++)
				{
				cacheuv_(cz2[b]);
				//if (cz2[b] >= 1) continue;
				// qqb Channel
				tgagaga += intfuncqt(fx2p[b][j],fx2[j],x2,cz2[b],p1qqp1qqp1qqreg_(cz2[b]),p1qqp1qqp1qqd0_(cz2[b]),p1qqp1qqp1qqd1_(cz2[b]),p1qqp1qqp1qqd2_(cz2[b]),p1qqp1qqp1qqdelta_(cz2[b]))*fx1[i] *abint::abw[b]; //PqqPqqPqq
				tgagaga += - 2.*lx2*(sumfx2p[b])*p1gqp1qqp1qg_(cz2[b])*fx1[i]*abint::abw[b]; //PqqPqgPgq
				// PqqPqgPgq The Contribution is the same as the previous one ->  we simply put a factor 2
				tgagaga += - lx2*(sumfx2p[b])*p1qgp1ggp1gq_(cz2[b])*fx1[i] *abint::abw[b]; //PqgPggPgq
				// qg Channel
				tgagaga += - lx2*fx2p[b][g]*p1qqp1qqp1qg_(cz2[b])*fx1[i]*abint::abw[b]; //PqqPqqPqg
				tgagaga += - lx2*fx2p[b][g]*p1qqp1qgp1gg_(cz2[b])*fx1[i]*abint::abw[b]; //PqqPqgPgg
				tgagaga += - lx2*fx2p[b][g]*p1qgp1gqp1qg_(cz2[b])*fx1[i]*abint::abw[b]*2.*NF;//PqgPgqPqg  ! Candidate for Nf factor?!!!!
				tgagaga += - lx2*fx2p[b][g]*p1qgp1ggp1gg_(cz2[b])*fx1[i]*abint::abw[b]; //PqgPggPgg
				// Start contributions for simple products contributions ! Checked 09.05.2020
				diffqq2  += (-lx2*(fx2p[b][j]-fx2[j]*cz2[b])*pqq_(cz2[b])- pqqint_(x2)*fx2[j])*abint::abw[b];
				diffg20  += -lx2*fx2p[b][g]*dypqg_(cz2[b])*abint::abw[b];
				tgagaqL2 += (-lx2*((fx2p[b][j]-fx2[j]*cz2[b]) *(D0qqqq/(1.-cz2[b])+D1qqqq*log(1.-cz2[b])/(1.-cz2[b]))
						+fx2p[b][j]*pqqqq_(cz2[b]))+(Deltaqqqq-D0qqqq*d0intx2-D1qqqq*d1intx2)*fx2[j])*abint::abw[b];
				tgagaqL2 += -lx2*(sumfx2p[b])*pqggq_(cz2[b])*abint::abw[b];
				tgagagL2 += -lx2*(fx2p[b][g]*(pqqqg_(cz2[b])+pqggg_(cz2[b])))*abint::abw[b];
				}
			diffqg1 = diffg10;
			diffqg2 = diffg20;
			diffg1f = diffqq1;
			diffg2f = diffqq2;
			// Start simple products contributions ! Checked 09.05.2020
			// qqb Channel
			tgagaga  += +3.*tgagaqL1*diffqq2;
			tgagaga  += +3.*tgagaqL2*diffqq1;
			// End qqb Channel
			// qg Channel
			tgagaga  += + 3.*tgagaqL2*diffqg1;
			tgagaga  += + 3.*tgagaqL1*diffqg2;
			tgagaga  += + 3.*tgagagL2*diffqq1;
			tgagaga  += + 3.*tgagagL1*diffqq2;
			// End qg Channel
			// gg Channel
			tgagaga  += + 3.*tgagagL1*diffg20;
			tgagaga  += + 3.*tgagagL2*diffg10;
			// End gg Channel
			// End tgagaga Contribution
			// ----------------------------------------------------------------------------------------------------------------------------------------

			//		sig21[sp] = tgagaga;
			//       	if (order == 2) continue;

			// ----------------------------------------------------------------------------------------------------------------------------------------
			// tcgaga
			// ----------------------------------------------------------------------------------------------------------------------------------------
			for (int a = 0; a < abint::abdim; a++)
				{
				cacheuv_(cz1[a]);
				//if (cz1[a] >= 1) continue;
				// qqb Channel
				// CqqPqqPqq Leg 1 Typo solved on 09.05.2020
				tcgaga += - lx1*fx1p[a][i]*c1qqp1qqp1qq_(cz1[a])*fx2[j]*abint::abw[a];
				// CqqPqgPgq Leg 1 Checked 09.05.2020
				tcgaga += - lx1*(sumfx1p[a])*c1qqp1qgp1gq_(cz1[a])*fx2[j]*abint::abw[a];
				// CqgPgqPqq Leg 1 Checked 09.05.2020
				tcgaga += - lx1*(sumfx1p[a])*c1qgp1gqp1qq_(cz1[a])*fx2[j]*abint::abw[a];
				// CqgPggPgq Leg 1 Checked 09.05.2020
				tcgaga += - lx1*(sumfx1p[a])*c1qgp1ggp1gq_(cz1[a])*fx2[j]*abint::abw[a];
				// End qqb Channel
				//       ! qg Channel Checked 09.05.2020
				// CqqPqqPqg Leg 1
				tcgaga += - lx1*fx1p[a][g]*c1qqp1qqp1qg_(cz1[a])*fx2[j]*abint::abw[a];
				// CqqPqgPgg Leg 1
				tcgaga += - lx1*fx1p[a][g]*c1qqp1qgp1gg_(cz1[a])*fx2[j]*abint::abw[a];
				// CqgPgqPqg Leg 1
				tcgaga += - lx1*fx1p[a][g]*c1qgp1gqp1qg_(cz1[a])*fx2[j]*2.*NF*abint::abw[a]; // Candidate for NF factor?!!!!
				// CqgPggPgg Leg 1
				tcgaga += - lx1*fx1p[a][g]*c1qgp1ggp1gg_(cz1[a])*fx2[j]*abint::abw[a];
				// Start contributions for simple products contributions Checked 09.05.2020
				tcgahqL1 += (fx1p[a][i]*cqqpqq_(cz1[a]))*(-lx1)*abint::abw[a];
				tcgahqL1 += + sumfx1p[a]*cqgpgq_(cz1[a])*(-lx1)*abint::abw[a];
				tcgahgL1 += fx1p[a][g]*(cqqpqg_(cz1[a])+cqgpgg_(cz1[a]))*(-lx1)*abint::abw[a];
				diffc1f += (-lx1*fx1p[a][i]*cqq_(cz1[a])+0.0*fx1[i]) * abint::abw[a];
				}
			for (int b = 0; b < abint::abdim; b++)
				{
				cacheuv_(cz2[b]);
				//if (cz2[b] >= 1) continue;
				// qqb Channel
				// CqqPqqPqq Leg 2 Typo solved on 09.05.2020
				tcgaga += - lx2*fx2p[b][j]*c1qqp1qqp1qq_(cz2[b])*fx1[i]*abint::abw[b];
				// CqqPqgPgq Leg 2 Checked 09.05.2020
				tcgaga += - lx2*(sumfx2p[b])*c1qqp1qgp1gq_(cz2[b])*fx1[i]*abint::abw[b];
				// CqgPgqPqq Leg 2 Checked 09.05.2020
				tcgaga += - lx2*(sumfx2p[b])*c1qgp1gqp1qq_(cz2[b])*fx1[i]*abint::abw[b];
				// CqgPggPgq Leg 2 Checked 09.05.2020
				tcgaga += - lx2*(sumfx2p[b])*c1qgp1ggp1gq_(cz2[b])*fx1[i]*abint::abw[b];
				//       ! qg Channel Checked 09.05.2020
				// CqqPqqPqg Leg 2 Typo solved 09.05.2020
				tcgaga += - lx2*fx2p[b][g]*c1qqp1qqp1qg_(cz2[b])*fx1[i]*abint::abw[b];
				// CqqPqgPgg Leg 2 Typo solved 09.05.2020
				tcgaga += - lx2*fx2p[b][g]*c1qqp1qgp1gg_(cz2[b])*fx1[i]*abint::abw[b];
				// CqgPgqPqg Leg 2 Typo solved 09.05.2020
				tcgaga += - lx2*fx2p[b][g]*c1qgp1gqp1qg_(cz2[b])*fx1[i]*2.*NF*abint::abw[b]; // Candidate for NF factor?!!!!
				// CqgPggPgg Leg 2 Typo solved 09.05.2020
				tcgaga += - lx2*fx2p[b][g]*c1qgp1ggp1gg_(cz2[b])*fx1[i]*abint::abw[b];
				// Start contributions for simple products contributions Checked 09.05.2020
				tcgahqL2 += (fx2p[b][j]*cqqpqq_(cz2[b]))*(-lx2)*abint::abw[b];
				tcgahqL2 += + sumfx2p[b]*cqgpgq_(cz2[b])*(-lx2)*abint::abw[b];
				tcgahgL2 += fx2p[b][g]*(cqqpqg_(cz2[b])+cqgpgg_(cz2[b]))*(-lx2)*abint::abw[b];
				diffc2f += (-lx2*fx2p[b][j]*cqq_(cz2[b])+0.0*fx2[j]) * abint::abw[b];
				}
			// Start simple products contributions ! Checked 09.05.2020
			// qqb Channel
			tcgaga += + diffc2f*tgagaqL1;
			tcgaga += + diffc1f*tgagaqL2;
			tcgaga += + 2.*diffg1f*tcgahqL2;
			tcgaga += + 2.*diffg2f*tcgahqL1;
			// End qqb Channel
			// qg Channel
			tcgaga += + diffc20*tgagaqL1;
			tcgaga += + diffc10*tgagaqL2;
			tcgaga += + diffc1f*tgagagL2;
			tcgaga += + diffc2f*tgagagL1;
			tcgaga += + 2.*tcgahqL1*diffg20;
			tcgaga += + 2.*tcgahqL2*diffg10;
			tcgaga += + 2.*tcgahgL1*diffg2f;
			tcgaga += + 2.*tcgahgL2*diffg1f;
			// End qg Channel
			// gg Channel Checked 09.05.2020
			tcgaga += + diffc20*tgagagL1;
			tcgaga += + diffc10*tgagagL2;
			tcgaga += + 2.*diffg20*tcgahgL1;
			tcgaga += + 2.*diffg10*tcgahgL2;
			// End gg Channel
			// ----------------------------------------------------------------------------------------------------------------------------------------
			//		sig21[sp] = tcgaga;
			//       	if (order == 2) continue;
			// ----------------------------------------------------------------------------------------------------------------------------------------
			// tc2ga1
			// ----------------------------------------------------------------------------------------------------------------------------------------
			//Leg 1
			for (int a = 0; a < abint::abdim; a++)
				{
				//if (cz1[a] >= 1) continue;
				// qqb Channel
				tc2ga1 += + intfuncqt(fx1p[a][i],fx1[i],x1,cz1[a],p1qqc2qqhsreg1[a],p1qqc2qqhsd01[a],p1qqc2qqhsd11[a],0.,p1qqc2qqhsdelta1[a])*fx2[j]*abint::abw[a]; // P1qqC2qq
				// qq
				tc2ga1 += + fx1p[a][im]*p1qqc2qqb1[a]*(-lx1)*fx2[j]*abint::abw[a]; // P1qqC2qqb
				// qqp
				tc2ga1 += + (-lx1)*(fx1pafx1pmanoti[i][a])*p1qqc2qqp1[a]*fx2[j]*abint::abw[a]; // P1qqC2qqp
				tc2ga1 += + (-lx1)*(sumfx1p[a])*p1gqc2qg1[a]*fx2[j]*abint::abw[a]; // P1gqC2qg
				// qg Channel
				// General C2qq*P1qg
				tc2ga1 += + (-lx1)*fx1p[a][g]*p1qgc2qq1[a]*fx2[j]*abint::abw[a];
				tc2ga1 += + (-lx1)*fx1p[a][g]*p1qgc2qqb1[a]*fx2[j]*abint::abw[a];
				tc2ga1 += + (-lx1)*fx1p[a][g]*p1qgc2qqp1[a]*fx2[j]*2.*(NF-1.)*abint::abw[a];//OK
				tc2ga1 += + (-lx1)*fx1p[a][g]*p1ggc2qg1[a]*fx2[j]*abint::abw[a]; // P1ggC2qg
				// Start contributions for simple products contributions
				diff1  = -lx1*(fx1p[a][i]-fx1[i]*cz1[a])*H2qqD0/(1.-cz1[a])*abint::abw[a];
				C2qqL1 +=  0.5*diff1;
				C2qqL1 += -0.5*H2qqD0*d0intx1*fx1[i]*abint::abw[a];
				C2qqL1 += + fx1p[a][i]*c2qqreghschm1[a]*(-lx1)*abint::abw[a];
				C2qqL1 += + fx1p[a][im]*c2qqb1[a]*(-lx1)*abint::abw[a];
				C2qqL1 += + fx1pafx1pmanoti[i][a]*c2qqp1[a]*(-lx1)*abint::abw[a]; // Problem in fx1pafx1pmanoti[i][a] FIXED 12.10.2021
				C2qgL1 += fx1p[a][g]*c2qghschm1[a]*(-lx1)*abint::abw[a];
				}
			//Leg 2
			for (int b = 0; b < abint::abdim; b++)
				{
				//if (cz2[b] >= 1) continue;
				// qqb Channel
				tc2ga1 += + intfuncqt(fx2p[b][j],fx2[j],x2,cz2[b],p1qqc2qqhsreg2[b],p1qqc2qqhsd02[b],p1qqc2qqhsd12[b],0.,p1qqc2qqhsdelta2[b])*fx1[i]*abint::abw[b]; // P1qqC2qq
				// qq
				tc2ga1 += + fx2p[b][jm]*p1qqc2qqb2[b]*(-lx2)*fx1[i]*abint::abw[b]; // P1qqC2qqb
				// qqp
				tc2ga1 += + (-lx2)*(fx2pafx2pmanotj[j][b])*p1qqc2qqp2[b]*fx1[i]*abint::abw[b]; // P1qqC2qqp
				tc2ga1 += + (-lx2)*(sumfx2p[b])*p1gqc2qg2[b]*fx1[i]*abint::abw[b]; // P1gqC2qg
				// qg Channel
				// General C2qq*P1qg
				tc2ga1 += + (-lx2)*fx2p[b][g]*p1qgc2qq2[b]*fx1[i]*abint::abw[b];
				tc2ga1 += + (-lx2)*fx2p[b][g]*p1qgc2qqb2[b]*fx1[i]*abint::abw[b];
				tc2ga1 += + (-lx2)*fx2p[b][g]*p1qgc2qqp2[b]*fx1[i]*2.*(NF-1.)*abint::abw[b]; // OK
				tc2ga1 += + (-lx2)*fx2p[b][g]*p1ggc2qg2[b]*fx1[i]*abint::abw[b]; // P1ggC2qg
				// Start contributions for simple products contributions
				diff2  = -lx2*(fx2p[b][j]-fx2[j]*cz2[b])*H2qqD0/(1.-cz2[b])*abint::abw[b];
				C2qqL2 += 0.5*diff2;
				C2qqL2 += -0.5*H2qqD0*d0intx2*fx2[j]*abint::abw[b];
				C2qqL2 += + fx2p[b][j]*c2qqreghschm2[b]*(-lx2)*abint::abw[b];
				C2qqL2 += + fx2p[b][jm]*c2qqb2[b]*(-lx2)*abint::abw[b];
				C2qqL2 += + fx2pafx2pmanotj[j][b]*c2qqp2[b]*(-lx2)*abint::abw[b]; // Problen in fx1pafx1pmanoti[i][a] FIXED 12.10.2021
				C2qgL2 += fx2p[b][g]*c2qghschm2[b]*(-lx2)*abint::abw[b];
				}
			// Start simple products contributions
			// qqb Channel
			tc2ga1 += + diffg2f*C2qqL1; // Leg 1     OK
			tc2ga1 += + diffg1f*C2qqL2; // Leg 2     OK
			// qg Channel
			tc2ga1 += + diffg1f*C2qgL2;//OK
			tc2ga1 += + diffg2f*C2qgL1;//OK
			tc2ga1 += + diffg10*C2qqL2;//OK
			tc2ga1 += + diffg20*C2qqL1;//OK
			// gg Channel
			tc2ga1 += + diffg10*C2qgL2;//OK
			tc2ga1 += + diffg20*C2qgL1;//OK
			// ----------------------------------------------------------------------------------------------------------------------------------------
			// sig21[sp] = tc2ga1;
			// if (order == 2) continue;
			// ----------------------------------------------------------------------------------------------------------------------------------------
			// tc1ga2
			// ----------------------------------------------------------------------------------------------------------------------------------------
			for (int a = 0; a < abint::abdim; a++)
			{
				cacheuv_(cz1[a]);
				//if (cz1[a] >= 1) continue;
				// tc1ga2 Checked 10.05.2020
				// qq Channel
				// Singlet case C1qqP2qqS
				//  Leg 1
				tc1ga2 += + sumfx1p[a] * c1qqp2qqs_(cz1[a])*(-lx1)*fx2[j]*abint::abw[a];
				// c------ Non-singlet case C1qqP2qqV
				// // Leg 1
				tc1ga2 += + fx1p[a][i]*c1qqp2qqv_(cz1[a])*(-lx1)*fx2[j]*abint::abw[a];
				// // Non-singlet case C1qqP2qqbV
				// // Leg 1
				tc1ga2 += + fx1p[a][im]*c1qqp2qqbv_(cz1[a])*(-lx1)*fx2[j]*abint::abw[a];
				// C1qgP2gq
				// Leg 1
				tc1ga2 += + (-lx1)*(sumfx1p[a])*c1qgp2gq_(cz1[a])*fx2[j]*abint::abw[a];
				// qg Channel
				// Leg 1
				tc1ga2 += + fx1p[a][g]*(c1qqp2qg_(cz1[a])+c1qgp2gg_(cz1[a]))*(-lx1)*fx2[j]*abint::abw[a];
				// Start contributions for simple products contributions
				gamma2qqL1 += sumfx1p[a]*p2qqs_(cz1[a])*(-lx1)*abint::abw[a];
				gamma2qqL1 += + fx1p[a][i]*p2qqv_(cz1[a])*(-lx1)*abint::abw[a];
				gamma2qqL1 += fx1[i]*resconst::Delta2qq*abint::abw[a];
				gamma2qqL1 += + 2./3.*Kappa*(-lx1*(fx1p[a][i]-fx1[i]*cz1[a])/(1.-cz1[a])-d0intx1*fx1[i])*abint::abw[a];
				gamma2qqL1 += + fx1p[a][im]*p2qqbv_(cz1[a])*(-lx1)*abint::abw[a];
				gamma2qgL1 += fx1p[a][g]*p2qg_(cz1[a])*(-lx1)*abint::abw[a];
			}
			for (int b = 0; b < abint::abdim; b++)
				{
				cacheuv_(cz2[b]);
				//if (cz2[b] >= 1) continue;
				// qqb Channel
				// Leg 2
				tc1ga2 += + sumfx2p[b] * c1qqp2qqs_(cz2[b])*(-lx2)*fx1[i]*abint::abw[b];
				// c------ Non-singlet case C1qqP2qqV
				// Leg 2
				tc1ga2 += + fx2p[b][j]*c1qqp2qqv_(cz2[b])*(-lx2)*fx1[i]*abint::abw[b];
				// Non-singlet case C1qqP2qqbV
				// Leg 2
				tc1ga2 += + fx2p[b][jm]*c1qqp2qqbv_(cz2[b])*(-lx2)*fx1[i]*abint::abw[b];
				// C1qgP2gq
				// Leg 2
				tc1ga2 += + (-lx2)*(sumfx2p[b])*c1qgp2gq_(cz2[b])*fx1[i]*abint::abw[b];
				// qg Channel
				// Leg 2
				tc1ga2 += + fx2p[b][g]*(c1qqp2qg_(cz2[b])+c1qgp2gg_(cz2[b]))*(-lx2)*fx1[i]*abint::abw[b];
				// Start contributions for simple products contributions
				gamma2qqL2 +=  sumfx2p[b]*p2qqs_(cz2[b])*(-lx2)*abint::abw[b];
				gamma2qqL2 += + fx2p[b][j]*p2qqv_(cz2[b])*(-lx2)*abint::abw[b];
				gamma2qqL2 +=   fx2[j]*resconst::Delta2qq*abint::abw[b];
				gamma2qqL2 += + 2./3.*Kappa*(-lx2*(fx2p[b][j]-fx2[j]*cz2[b])/(1.-cz2[b])-d0intx2*fx2[j])*abint::abw[b];
				gamma2qqL2 += + fx2p[b][jm]*p2qqbv_(cz2[b])*(-lx2)*abint::abw[b];
				gamma2qgL2 += fx2p[b][g]*p2qg_(cz2[b])*(-lx2)*abint::abw[b];
				}
			// Start simple products contributions
			// Leg 1
			tc1ga2 += + diffc2f * gamma2qqL1;
			// // Leg 2
			tc1ga2 += + diffc1f * gamma2qqL2;
			// // End qqb Channel
			// // qg Channel
			// // Leg 1
			tc1ga2 += + diffc20*gamma2qqL1;
			tc1ga2 += + diffc2f*gamma2qgL1;
			// // Leg 2
			tc1ga2 += + diffc10*gamma2qqL2;
			tc1ga2 += + diffc1f*gamma2qgL2;
			// End qg Channel
			// gg Channel
			// Leg 1
			tc1ga2 += + diffc20*gamma2qgL1;
			// Leg 2
			tc1ga2 += + diffc10*gamma2qgL2;
			// End gg Channel
			// ----------------------------------------------------------------------------------------------------------------------------------------
			//sig21[sp] = tc1ga2;
			//if (order == 2) continue;
			// ----------------------------------------------------------------------------------------------------------------------------------------
			// tga1ga2
			// ----------------------------------------------------------------------------------------------------------------------------------------
			for (int a = 0; a < abint::abdim; a++)
				{
				cacheuv_(cz1[a]);
				//if (cz1[a] >= 1) continue;
				// tga1ga2 Checked 10.05.2020
				// Usage IntFuncqT(fxip,fxi0,xxi0,z,Preg,PD0,PD1,PD2,Pdelta)
				// qq Channel
				// Leg 1
				tga1ga2 += + 2.*0.5*intfuncqt(fx1p[a][i],fx1[i],x1,cz1[a],p1qqp2qqvreg_(cz1[a]),p1qqp2qqvd0_(cz1[a]),p1qqp2qqvd1_(cz1[a]),0.,p1qqp2qqvdelta_(cz1[a]))*fx2[j]*abint::abw[a];
				// Singlet case P2qqS
				// Leg 1
				tga1ga2 += + 2.*0.5*sumfx1p[a] *p1qqp2qqs_(cz1[a])*(-lx1)*fx2[j]*abint::abw[a];
				// P2qqb non singlet P2qqbV
				// Leg 1
				tga1ga2 += + 2.*0.5*fx1p[a][im]*p1qqp2qqbv_(cz1[a])*(-lx1)*fx2[j]*abint::abw[a];
				// P1qgP2gq +P1gqP2qg
				// Leg 1
				tga1ga2 += + 0.5*(-lx1)*(sumfx1p[a])*(p1qgp2gq_(cz1[a])+p1gqp2qg_(cz1[a]))*fx2[j]*abint::abw[a];
				// End qqb Channel
				// qg Channel
				// P1qqP2qg + P1qgP2gg + P1ggP2qg
				// Leg 1
				tga1ga2 += + 0.5*(-lx1)*fx1p[a][g]*(p1qqp2qg_(cz1[a])+p1qgp2gg_(cz1[a])+p1ggp2qg_(cz1[a]))*fx2[j]*abint::abw[a];
				// P1qgP2qq (all contributions inside)
				// Leg 1
				tga1ga2 += + 0.5*(-lx1)*fx1p[a][g]*p1qgp2qq_(cz1[a])*fx2[j]*abint::abw[a];
				}
			for (int b = 0; b < abint::abdim; b++)
				{
				cacheuv_(cz2[b]);
				//if (cz2[b] >= 1) continue;
				// qqb Channel
				// Leg 2
				tga1ga2 += + 2.*0.5*intfuncqt(fx2p[b][j],fx2[j],x2,cz2[b],p1qqp2qqvreg_(cz2[b]),p1qqp2qqvd0_(cz2[b]),p1qqp2qqvd1_(cz2[b]),0.,p1qqp2qqvdelta_(cz2[b]))*fx1[i]*abint::abw[b];
				// Singlet case P2qqS
				// Leg 2
				tga1ga2 += + 2.*0.5*sumfx2p[b] *p1qqp2qqs_(cz2[b])*(-lx2)*fx1[i]*abint::abw[b];
				// P2qqb non singlet P2qqbV
				// Leg 2
				tga1ga2 += + 2.*0.5*fx2p[b][jm]*p1qqp2qqbv_(cz2[b])*(-lx2)*fx1[i]*abint::abw[b];
				// P1qgP2gq +P1gqP2qg
				// Leg 2
				tga1ga2 += + 0.5*(-lx2)*(sumfx2p[b])*(p1qgp2gq_(cz2[b])+p1gqp2qg_(cz2[b]))*fx1[i]*abint::abw[b];
				// End qqb Channel
				// qg Channel
				// P1qqP2qg + P1qgP2gg + P1ggP2qg
				// Leg 2
				tga1ga2 += + 0.5*(-lx2)*fx2p[b][g]*(p1qqp2qg_(cz2[b])+p1qgp2gg_(cz2[b])+p1ggp2qg_(cz2[b]))*fx1[i]*abint::abw[b];
				// P1qgP2qq (all contributions inside)
				// Leg 2
				tga1ga2 += + 0.5*(-lx2)*fx2p[b][g]*p1qgp2qq_(cz2[b])*fx1[i]*abint::abw[b];
				}
			// Start simple products contributions
			// Leg 1
			tga1ga2 += + gamma2qqL2*diffqq1;
			// Leg 2
			tga1ga2 += + gamma2qqL1*diffqq2;
			// end qqb Channel
			// Start simple products contributions
			// qg Channel
			// Leg 1
			tga1ga2 += + gamma2qgL2*diffqq1;
			// Leg 2
			tga1ga2 += + gamma2qgL1*diffqq2;
			// Leg 1
			tga1ga2 += + gamma2qqL2*diffqg1;
			// Leg 2
			tga1ga2 += + gamma2qqL1*diffqg2;
			// End qg Channe
			// gg Channel
			// Leg 1
			tga1ga2 += + diffg10*gamma2qgL2; // Corrected 10.05.2020
			// Typo found : tga1ga2 += + diffqg1*gamma2qgL1
			// Leg 2
			tga1ga2 += + diffg20*gamma2qgL1; // Corrected 10.05.2020
			// Typo found : tga1ga2 += + diffqg2*gamma2qgL2
			// End gg Channel
			//   // ----------------------------------------------------------------------------------------------------------------------------------------
			//   sig21[sp] = tga1ga2;
			//   if (order == 2) continue;
			// ----------------------------------------------------------------------------------------------------------------------------------------
			// tccga
			// ----------------------------------------------------------------------------------------------------------------------------------------
			// tccga !Checked 09.05.2020
			// qq Channel
			// Start simple products contributions
			// Leg 1
			tccga += + tcgahqL1*diffc2f; //Checked 09.05.2020
			// Leg 2
			tccga += + tcgahqL2*diffc1f; //Checked 09.05.2020
			// End qqb Channel
			// qg Channel
			// Start simple products contributions
			// Leg 1
			tccga += + tcgahqL1*diffc20; //Checked 09.05.2020
			// Leg 2
			tccga += + tcgahqL2*diffc10; //Checked 09.05.2020
			// Leg 1
			tccga += + tcgahgL1*diffc2f; //Checked 09.05.2020
			// Leg 2
			tccga += + tcgahgL2*diffc1f; //Checked 09.05.2020
			// End qg Channel
			// gg Channel
			// Leg 1
			tccga += + tcgahgL1*diffc20; //Checked 09.05.2020
			// Leg 2
			tccga += + tcgahgL2*diffc10; //Checked 09.05.2020
				// end gg Channel
			// // ----------------------------------------------------------------------------------------------------------------------------------------
			// sig21[sp] = tccga;
			// if (order == 2) continue;
			// ----------------------------------------------------------------------------------------------------------------------------------------
			// tgamma3
			// ----------------------------------------------------------------------------------------------------------------------------------------
			//Leg 1
			for (int a = 0; a < abint::abdim; a++)
				{
				cacheuv_(cz1[a]);
				//if (cz1[a] >= 1) continue;
				//qq Channel
				tgamma3 += + intfuncqt(fx1p[a][i],fx1[i],x1,cz1[a],p3qqv1[a],p3qqvd01[a],0.,0.,p3qqvdelta1[a])*fx2[j]*abint::abw[a]; //non singlet P3qqV
				tgamma3 += + fx1p[a][im]*p3qqbv1[a]*(-lx1)*fx2[j]*abint::abw[a];		                                           //non singlet P3qqbV
				//singlet
				if (j > g)
				{
				tgamma3 += + (sumnegfx1p[a])*p3sqq1[a]*(-lx1)*fx2[j]*abint::abw[a];
				tgamma3 += + (sumposfx1p[a])*p3sqqb1[a]*(-lx1)*fx2[j]*abint::abw[a];
				}
				else
				{
				tgamma3 += + (sumposfx1p[a])*p3sqq1[a]*(-lx1)*fx2[j]*abint::abw[a];
				tgamma3 += + (sumnegfx1p[a])*p3sqqb1[a]*(-lx1)*fx2[j]*abint::abw[a];
				}
				//qg Channel
				tgamma3+=+fx1p[a][g]*p3qg1[a]*(-lx1)*fx2[j]*abint::abw[a];
				}
			//Leg 2
			for (int b = 0; b < abint::abdim; b++)
				{
				cacheuv_(cz2[b]);
				//if (cz2[b] >= 1) continue;
				//qq Channel
				tgamma3 += + intfuncqt(fx2p[b][j],fx2[j],x2,cz2[b],p3qqv2[b],p3qqvd02[b],0.,0.,p3qqvdelta2[b])*fx1[i]*abint::abw[b]; //non singlet P3qqV
				tgamma3 += + fx2p[b][jm]*p3qqbv2[b]*(-lx2)*fx1[i]*abint::abw[b];                                                     //non singlet P3qqbV
				//singlet
				if (i > g)
				{
				tgamma3 += + (sumnegfx2p[b])*p3sqq2[b]*(-lx2)*fx1[i]*abint::abw[b];
				tgamma3 += + (sumposfx2p[b])*p3sqqb2[b]*(-lx2)*fx1[i]*abint::abw[b];
				}
				else
				{
				tgamma3 += + (sumposfx2p[b])*p3sqq2[b]*(-lx2)*fx1[i]*abint::abw[b];
				tgamma3 += + (sumnegfx2p[b])*p3sqqb2[b]*(-lx2)*fx1[i]*abint::abw[b];
				}
				//qg Channel
				tgamma3+=+fx2p[b][g]*p3qg2[b]*(-lx2)*fx1[i]*abint::abw[b];
				}
			// ----------------------------------------------------------------------------------------------------------------------------------------
			//sig21[sp] = tgamma3;
			//if (order == 2) continue;
			// ----------------------------------------------------------------------------------------------------------------------------------------
			// C1caC1cbarb and tH2st
			// ----------------------------------------------------------------------------------------------------------------------------------------
			// Leg 1
			for (int a = 0; a < abint::abdim; a++)
				{
				cacheuv_(cz1[a]);
				if (cz1[a] >= 1) continue;
				C1cbarbpiuC1ca += +(fx1p[a][i]*cqq_(cz1[a])+fx1p[a][g]*cqg_(cz1[a]))*(-lx1)*fx2[j]*abint::abw[a];
				tH2st += +fx1p[a][i]*cqq_(cz1[a])*(-lx1)*fx2[j]*C1qqdelta*0.0*abint::abw[a];// Since in the Hard Scheme C1qqdelta is ZERO
				tH2st += +fx1p[a][g]*cqg_(cz1[a])*(-lx1)*fx2[j]*C1qqdelta*0.0*abint::abw[a];// Since in the Hard Scheme C1qqdelta is ZERO
				tH2st += +0.5*(-lx1*(fx1p[a][i]-fx1[i]*cz1[a])*H2qqD0/(1.-cz1[a]))*fx2[j]*abint::abw[a];
				tH2st += +fx1p[a][i]*c2qqreghschm_(cz1[a])*(-lx1)*fx2[j] *abint::abw[a];
				tH2st += +fx1p[a][g]*c2qghschm_(cz1[a])*(-lx1)*fx2[j] *abint::abw[a];
				tH2st += +fx1p[a][im]*c2qqb_(cz1[a])*(-lx1)*fx2[j] *abint::abw[a];
				tH2st += +fx1pafx1pmanoti[i][a]*c2qqp_(cz1[a])*(-lx1)*fx2[j] *abint::abw[a];
				tH2st += -0.5*H2qqD0*d0intx1*fx1[i]*fx2[j]*abint::abw[a];
				}
			// Leg 2
			for (int b = 0; b < abint::abdim; b++)
				{
				cacheuv_(cz2[b]);
				if (cz2[b] >= 1) continue;
				C1cbarbpiuC1ca += +(fx2p[b][j]*cqq_(cz2[b])+fx2p[b][g]*cqg_(cz2[b]))*(-lx2)*fx1[i]*abint::abw[b];
				tH2st += +fx2p[b][j]*cqq_(cz2[b])*(-lx2)*fx1[i]*C1qqdelta*0.0 *abint::abw[b];// Since in the Hard Scheme C1qqdelta is ZERO
				tH2st += +fx2p[b][g]*cqg_(cz2[b])*(-lx2)*fx1[i]*C1qqdelta*0.0 *abint::abw[b];// Since in the Hard Scheme C1qqdelta is ZERO
				tH2st += +0.5*(-lx2*(fx2p[b][j]-fx2[j]*cz2[b])*H2qqD0/(1.-cz2[b]))*fx1[i]*abint::abw[b];
				tH2st += +fx2p[b][j]*c2qqreghschm_(cz2[b])*(-lx2)*fx1[i]*abint::abw[b];
				tH2st += +fx2p[b][g]*c2qghschm_(cz2[b])*(-lx2)*fx1[i] *abint::abw[b];
				tH2st += +fx2p[b][jm]*c2qqb_(cz2[b])*(-lx2)*fx1[i] *abint::abw[b];
				tH2st += +fx2pafx2pmanotj[j][b]*c2qqp_(cz2[b])*(-lx2)*fx1[i]*abint::abw[b];
				tH2st += -0.5*H2qqD0*d0intx2*fx1[i]*fx2[j]*abint::abw[b];
				}
			// Start Simple product contributions
			C1caC1cbarb += diffc10*diffc20;
			C1caC1cbarb += diffc1f*diffc2f;
			C1caC1cbarb += diffc10*diffc2f;
			C1caC1cbarb += diffc1f*diffc20;
			// sig21[sp] = tH2st+C1caC1cbarb;
			// if (order == 2) continue;
			// ----------------------------------------------------------------------------------------------------------------------------------------
			// ----------------------------------------------------------------------------------------------------------------------------------------
			// tcgah
			// ----------------------------------------------------------------------------------------------------------------------------------------
			for (int a = 0; a < abint::abdim; a++)
				{
				if (cz1[a] >= 1) continue;
				tcgah += ( fx1p[a][i] *cqqpqq_(cz1[a])+ fx1p[a][g] *(cqqpqg_(cz1[a])+cqgpgg_(cz1[a])))*(-lx1)*fx2[j]*abint::abw[a];
				tcgah += ( sumfx1p[a] )*cqgpgq_(cz1[a])*(-lx1)*fx2[j]*abint::abw[a];

				}
			for (int b = 0; b < abint::abdim; b++)
				{
				if (cz2[b] >= 1) continue;
				tcgah += ( fx2p[b][j] *cqqpqq_(cz2[b])+ fx2p[b][g] *(cqqpqg_(cz2[b])+cqgpgg_(cz2[b])))*(-lx2)*fx1[i]*abint::abw[b];
				tcgah += ( sumfx2p[b] )*cqgpgq_(cz2[b])*(-lx2)*fx1[i]*abint::abw[b];
				}
			// Start simple product contributions
			tcgah += diffg10*diffc20+ diffg1f*diffc2f +diffg10*diffc2f+ diffg1f*diffc20;
				tcgah += diffc10*diffg20+ diffc1f*diffg2f +diffc10*diffg2f+ diffc1f*diffg20;
			// ----------------------------------------------------------------------------------------------------------------------------------------
			// Coefficients for Third order
			double H1q = 2.*(pi2-8.)/3.;
			double H2q = -2561./144.+ 127.*NF/72. + 3.*pi2/2. -19.*NF*pi2/81.+ 49.*pi4/324. + 58.*zeta3/9. + 8.*NF*zeta3/27.;
			double CFH = 4./3.;
			double B2qHSCHM = (1.0/48.0)*CFH*(2.0*NF*(1.0
									+ 4.0*zeta2) +
								9.0*CFH*(-1.0 + 8.0*zeta2 - 16.0*zeta3)
								+    3.0*(-17.0 - 44.0*zeta2 + 72.0*zeta3));

			double B3qHSCHM = (-4623*CF)/32.
				+ (20149*CA*CA*CF)/1152. - (29*CF*CF)/8. +
				(11*CA*CF*CF)/384. - (29*CF*CF*CF)/64. + (763*CF*NF)/48. -
				(853*CA*CF*NF)/144. + (23*CF*CF*NF)/32. + (17*CF*NF*NF)/288. +
				(1295*CF*pi2)/72. - (353*CA*CA*CF*pi2)/162. +
				(131*CF*CF*pi2)/24. - (319*CA*CF*CF*pi2)/288. -
				(3*CF*CF*CF*pi2)/32. - (469*CF*NF*pi2)/216. +
				(983*CA*CF*NF*pi2)/1296. - (CF*CF*NF*pi2)/18. -
				(251*CF*pi4)/480. + (11*CA*CA*CF*pi4)/1440. +
				(CF*CF*pi4)/288. + (121*CA*CF*CF*pi4)/2160. -
				(CF*CF*CF*pi4)/20. + (11*CF*NF*pi4)/360. -
				(CA*CF*NF*pi4)/720. - (29*CF*CF*NF*pi4)/2160. +
				(2825*CF*zeta3)/24. - (649*CA*CA*CF*zeta3)/108. - 47*CF*CF*zeta3 +
				(55*CA*CF*CF*zeta3)/8. - (17*CF*CF*CF*zeta3)/8. - (215*CF*NF*zeta3)/36. +
				(37*CA*CF*NF*zeta3)/54. + (17*CF*CF*NF*zeta3)/12. +
				(CF*NF*NF*zeta3)/9. - (CF*CF*pi2*zeta3)/4. + (CF*CF*CF*pi2*zeta3)/6. -
				(45*CF*zeta5)/4. - (45*CF*CF*zeta5)/4. + (15*CF*CF*CF*zeta5)/2.;


			double pcF = 0.;

			double ThetaA=0.5;
			double ThetaB=B1q - (4.*beta0)/3. + A1q*LQ;
			double ThetaC=-(A1q*beta0)/3. + (A1q*(B1q + A1q*LQ))/2.;
			double ThetaD=0.5;
			double ThetaF=pow(beta0,2)/3. + pow(B1q + A1q*LQ,2)/2. -
				(beta0*(5*B1q + A1q*(2.*LQ + 3.*LR)))/6.;
			double Thetadelta=(-2*A1q*beta1 - 3.*A2q*(B1q + A1q*LQ)
						- pow(B1q + A1q*LQ,3) + beta0*(B1q + A1q*LQ)*(B1q +
											A1q*(-2.*LQ + 3*LR)))/6.;

			double CH1=(-A2q - pow(B1q + A1q*LQ,2) + beta0*(B1q + A1q*LR))/2.;

			double alfa32=(2.*H1q - 2.*LQ*(B1q + (A1q*LQ)/2.) +
					2*beta0*(LF + LQ - LR*(2. + pcF)))/2.;

			double Cdelta32=(-A3q - (2.*A2q*beta0 + A1q*(beta1 +
									pow(beta0,2)*(LQ - LR)))*(LQ - LR))/2.;

			double CL1=   (4.*B2qHSCHM + 4.*A2q*LQ + beta0*
					(2.*H1q + A1q*LQ*(3.*LQ - 4*LR) + 2.*B1q*(LQ - 2.*LR) -
					2.*beta0*LR*pcF))/4.;

			double C1L1=(-2.*B2qHSCHM + 2.*beta1 - 2.*A2q*LQ - 2.*pow(B1q,2)*LQ - 2.*B1q
					*beta0*LQ +
					4.*pow(beta0,2)*LQ - 3.*A1q*B1q*pow(LQ,2) - 2.*A1q*beta0*pow(LQ,2) -
					pow(A1q,2)*LQ3 + 2.*B1q*beta0*LR - 4.*pow(beta0,2)*LR +
					2.*A1q*beta0*LQ*LR - 2.*pow(beta0,2)*LR*pcF)/2.;

			double C1L2=pow(H1q,2) - H2q - beta1*LF + B2qHSCHM*LQ - 2.*pow(beta0,2)*LF*LQ+
				(A2q*LQ2)/2. + (pow(B1q,2)*LQ2)/2. + (B1q*beta0*LQ2)/2. +
				pow(beta0,2)*LQ2 + (A1q*B1q*LQ3)/2. + (A1q*beta0*LQ3)/3. +
				(pow(A1q,2)*LQ4)/8. + beta1*LR + 2.*pow(beta0,2)*LF*LR -
				B1q*beta0*LQ*LR - (A1q*beta0*LQ2*LR)/2.
				- pow(beta0,2)*LR2 +
				H1q*(B1q*LF - 2.*B1q*LQ + A1q*LF*LQ - (3.*A1q*LQ2)/2. +
				beta0*LR) + beta1*LR*pcF + pow(beta0,2)*LF*LR*pcF -
				(3*pow(beta0,2)*LR2*pcF)/2. - (pow(beta0,2)*LR2*pow(pcF,2))/2.;

			double C1L3= -B3qHSCHM - beta1*H1q + B1q*H2q + B1q*beta1*LF -
				pow(B1q,2)*H1q*LF - A3q*LQ +
				A2q*H1q*LQ + pow(B1q,2)*H1q*LQ + B1q*beta0*H1q*LQ -
				2.*pow(beta0,2)*H1q*LQ + A1q*H2q*LQ + 2.*B1q*pow(beta0,2)*LF*LQ +
				A1q*beta1*LF*LQ - 2.*A1q*B1q*H1q*LF*LQ -
				(3.*A2q*B1q*LQ2)/2. -
				(pow(B1q,3)*LQ2)/2. - 2.*A2q*beta0*LQ2 -
				(3.*pow(B1q,2)*beta0*LQ2)/2. - (A1q*beta1*LQ2)/2. +
				2.*A1q*B1q*H1q*LQ2 + A1q*beta0*H1q*LQ2 +
				2.*A1q*pow(beta0,2)*LF*LQ2 - pow(A1q,2)*H1q*LF*LQ2 - A1q*A2q
				*LQ3 -
				A1q*pow(B1q,2)*LQ3 - (7.*A1q*B1q*beta0*LQ3)/3. -
				A1q*pow(beta0,2)*LQ3 + pow(A1q,2)*H1q*LQ3 -
				(5.*pow(A1q,2)*B1q*LQ4)/8. -
				(5.*pow(A1q,2)*beta0*LQ4)/6. - (pow(A1q,3)*LQ5)/8. -
				2.*B1q*beta0*H1q*LR +2.*pow(beta0,2)*H1q*LR -
				2.*B1q*pow(beta0,2)*LF*LR +
				2.*A2q*beta0*LQ*LR + 2.*pow(B1q,2)*beta0*LQ*LR -
				2.*A1q*beta0*H1q*LQ*LR - 2.*A1q*pow(beta0,2)*LF*LQ*LR +
				3.*A1q*B1q*beta0*LQ2*LR + A1q*pow(beta0,2)*LQ2*LR +
				pow(A1q,2)*beta0*LQ3*LR +
				B2qHSCHM*(H1q - 2.*B1q*LQ - 2.*beta0*LQ - (3.*A1q*LQ2)/2. +
					2.*beta0*LR) - B1q*beta1*LR*pcF + beta0*beta1*LR*pcF -
				B1q*beta0*H1q*LR*pcF + pow(beta0,2)*H1q*LR*pcF -
				B1q*pow(beta0,2)*LF*LR*pcF - B1q*pow(beta0,2)*LQ*LR*pcF +
				2.*pow(beta0,3)*LQ*LR*pcF - A1q*beta1*LQ*LR*pcF -
				A1q*beta0*H1q*LQ*LR*pcF - A1q*pow(beta0,2)*LF*LQ*LR*pcF -
				(A1q*pow(beta0,2)*LQ2*LR*pcF)/2. +
				(3.*B1q*pow(beta0,2)*LR2*pcF)/2. - 2*pow(beta0,3)*LR2*pcF +
				(3.*A1q*pow(beta0,2)*LQ*LR2*pcF)/2. +
				(B1q*pow(beta0,2)*LR2*pow(pcF,2))/2. - pow(beta0,3)*LR2*pow(pcF,2) +
				(A1q*pow(beta0,2)*LQ*LR2*pow(pcF,2))/2.;

			double Hache1h =  C1cbarbpiuC1ca + H1q*tdelta -
				((2.*B1q*LQ + A1q*LQ*LQ + 2.*beta0*LR*pcF)*tdelta)/2. +
				LF*th1stF - LQ*th1stF;

			double Hache2h = tdelta*H2q + ((LQ*(-8.*B2qHSCHM - 4.*A2q*LQ +
								(2.*B1q + A1q*LQ)*(-4.*H1q
											+ LQ*(2.*B1q + A1q*LQ))))/8. -
							beta1*LR*pcF + (pow(beta0,2)*LR2*pcF*(1. + pcF))/2. +
							(beta0*(-6.*H1q*LR*(1. + pcF)-
								3.*B1q*LQ*(LQ - 2.*LR*(1. + pcF))+
								A1q*LQ2*(-2.*LQ + 3.*LR*(1. + pcF))))/6.)*tdelta;


			// ACTIVATE FOR FULL HARD SCHEME!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			Hache2h += C1caC1cbarb + C1cbarbpiuC1ca*H1q;
			// Now this two contributions are in the tH2st contribution
			// since the single coefficients C(2)_ij contain this information
			Hache2h += -
				B1q*C1cbarbpiuC1ca*LQ +
				beta0*C1cbarbpiuC1ca*LQ - (A1q*C1cbarbpiuC1ca*LQ2)/2. -
				beta0*C1cbarbpiuC1ca*LR - beta0*C1cbarbpiuC1ca*LR*pcF +
				LF*tcgah - LQ*tcgah + (LF2*tgaga)/2. - LF*LQ*tgaga +
				(LQ2*tgaga)/2. + LF*tgamma2 - LQ*tgamma2 +
				H1q*LF*th1stF + (beta0*LF2*th1stF)/2.
				- H1q*LQ*th1stF -
				B1q*LF*LQ*th1stF + B1q*LQ2*th1stF -
				(beta0*LQ2*th1stF)/2. - (A1q*LF*LQ2*th1stF)/2. +
				(A1q*LQ3*th1stF)/2. - beta0*LF*LR*th1stF +
				beta0*LQ*LR*th1stF - beta0*LF*LR*pcF*th1stF +
				beta0*LQ*LR*pcF*th1stF + tH2st;

			// Sigmas for Third order
			sig36[sp]=-(pow(A1q,3)*tdelta)/48.;

			sig35[sp]=-(A1q*((A1q*sig11[sp])/4. + sig23[sp]))/2.;

			sig34[sp]=-(pow(A1q,2)*Hache1h)/8. - (A1q*sig22[sp])/2. - (A1q*sig11[sp]*
											beta0)/3. - (A1q*pow(beta0,2)*tdelta)/4.;

			sig33[sp]=-(beta0*tgaga)/6. - tgagaga/6. - A1q*sig21[sp]*ThetaA -
				sig22[sp]*ThetaB - Hache1h*ThetaC - A2q*sig11[sp]*ThetaD +
				tdelta*Thetadelta - sig11[sp]*ThetaF;


			sig32[sp]=CH1*Hache1h - (A1q*Hache2h)/2. + (beta1*sig11[sp])/2. -
				(B1q - beta0 + A1q*LQ)*sig21[sp] - (beta0*tcgah)/2. + tcgaga/2. +
				Cdelta32*tdelta + tga1ga2 + (alfa32*tgaga)/2. +
				((LF - LQ)*tgagaga)/2. + CL1*th1stF;

			sig31[sp]=2.*beta0*C1caC1cbarb + C1L1*Hache1h +
				C1cbarbpiuC1ca*H1q*(B1q + A1q*LQ) - Hache2h*(B1q + A1q*LQ) -
				C1L2*sig11[sp] + (H1q - LQ*(B1q + (A1q*LQ)/2.))*sig21[sp] - tc1ga2 -
				tc2ga1 - tccga + beta0*(LF - 3.*LQ + LR*(2. + pcF))*tcgah +
				(-LF + LQ)*tcgaga + C1L3*tdelta + 2.*(-LF + LQ)*tga1ga2 +
				(beta0*(-LF + LQ)*(LF + 3.*LQ - 2.*LR*(2. + pcF))*tgaga)/2. -
				(pow(-LF + LQ,2)*tgagaga)/2. +
				beta0*(-2.*LQ + LR*(2. + pcF))*tgamma2 - tgamma3 + 2.*beta0*tH2st;



			}

			double xmsq = 0.;
			for (int sp = 0; sp < mesq::totpch; sp++)
			{
			//as/pi factor
			double sig1 = (sig12[sp]*qtint::LL2_mesqij[sp]+sig11[sp]*qtint::LL1_mesqij[sp])*asopi + (sig12_qed[sp]*qtint::LL2_mesqij[sp]+sig11_qed[sp]*qtint::LL1_mesqij[sp])*resconst::a0qed/M_PI;
			if (order < 3 || fpc) xmsq += -sig1;
			if (order == 1) continue;

			//(as/pi)^2 factor
			double sig2 = ((sig24[sp]*qtint::LL4_mesqij[sp]+sig23[sp]*qtint::LL3_mesqij[sp]+sig22[sp]*qtint::LL2_mesqij[sp])+sig21[sp]*qtint::LL1_mesqij[sp])*pow(asopi,2);
			// sig2 =  sig21[sp]*qtint::LL3_mesqij[sp]*pow(asopi,3);
			if (order < 3 || fpc) xmsq += -sig2;
			if (order == 2) continue;

			//(as/pi)^3 factor
			double sig3 = ((sig36[sp]*qtint::LL6_mesqij[sp]+sig35[sp]*qtint::LL5_mesqij[sp]+sig34[sp]*qtint::LL4_mesqij[sp]+sig33[sp]*qtint::LL3_mesqij[sp])+sig32[sp]*qtint::LL2_mesqij[sp]+(sig31[sp]*qtint::LL1_mesqij[sp]))*pow(asopi,3);
			xmsq += -sig3;
			}

			double shad = pow(opts.sroot,2);
			//double fbGeV2=0.38937966e12;
			//double fac = M_PI*6./fbGeV2*(2*m2);
			//double flux = fbGeV2/(2.*x1*x2*shad);
			//double ps = 1./shad;
			//double norm = 1./16./M_PI / 2. / M_PI;
			//xmsq = xmsq * ps * fac * flux * norm;

			xmsq = xmsq/shad * 3./8. /2./M_PI;
			//xmsq = xmsq/shad * 3./8. /2./M_PI/2./M_PI;

			xmsq = xmsq * m2; //Is this a Jacobian for dsigma/dq2 -> dsigma/dq?
			xmsq = xmsq * 2*M_PI; //phiV integration
			if (isnan_ofast(xmsq))
			cout << m << " " << y << " " << costh << "  " << xmsq << endl;

			//switching --> switching function is inside qtint, do not apply
			//double swtch;
			//if (mode == 0 || mode == 1)
			//  swtch = switching::swtch(qt, m);
			//else if (mode == 2)
			//  swtch=1.; // qt integration already performed

			//if (swtch < 0.01) return 0.;// do not apply this cut to avoid discontinuities. Instead the phase space is limited to qt and m switching limits
			xmsq = xmsq;//*swtch; //switching function is inside qtint

			//xmsq = xmsq*jac;//jacobian for the change of variable qt=qtp/sqrt(1-qtp^2/m^2)

			f[npdf] = xmsq;
			} //end loop on pdf

		return;
	}
}
