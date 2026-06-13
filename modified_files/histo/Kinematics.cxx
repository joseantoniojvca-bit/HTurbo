#ifndef Kinematics_CXX
#define Kinematics_CXX
/**
 * @file Kinematics.cxx
 * Description of this macro
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-08-26
 */

#include "Kinematics.h"
#include "KinematicDefinitions.h"
#include "KinematicCuts.h"
#include <iostream>

namespace Kinematics{

	// -- photons
	double exp2phycut;
	double exp2crack_eta_max;
	double exp2crack_eta_min;
	double fTmin12;
	double fTmin22;
	// -- leptons
	double exp2lycut;
	double lpt2cut;
	double mt2cut;
	double et2misscut;

	void init()
	{
		// -- Photon cuts
		exp2phycut = exp(2*opts.phycut);
		fTmin12 = opts.phpt_cut_h*opts.phpt_cut_h;
		fTmin22 = opts.phpt_cut_s*opts.phpt_cut_s;
		exp2crack_eta_max = exp(2*opts.crack_eta_max);
		exp2crack_eta_min = exp(2*opts.crack_eta_min);
		// -- lepton cuts
		exp2lycut = exp(2*opts.lycut);
		lpt2cut = pow(opts.lptcut,2);
		mt2cut = pow(opts.mtcut,2);
		et2misscut = pow(opts.etmisscut,2);
	}

	// input values definitions
	double p3[4];
	double p4[4];
	double event_weight;

	// All flags are stored here
	ObservableFlags flags;

	void SetKinematics(double _p3[4], double _p4[4], double wgt=1.){
		isIntegratorMode=false;

		p3[0]=_p3[0];
		p3[1]=_p3[1];
		p3[2]=_p3[2];
		p3[3]=_p3[3];
		p4[0]=_p4[0];
		p4[1]=_p4[1];
		p4[2]=_p4[2];
		p4[3]=_p4[3];

		event_weight=wgt;
		// Set all variables to re-calculate
		//for (auto ptr_isCalculated : flags) *ptr_isCalculated=false;
	}

	void SetMiddlePoint(){
		isIntegratorMode=true;
		for (auto ptr_isCalculated : flags) *ptr_isCalculated=false;
	}

	bool isIntegratorMode=false;

	namespace Cuts {

		const bool KeepEvent=true;
		const bool SkipEvent=false;

		struct StandardCuts : public CutBase {
			// define here all observables and variables
			LepCh      ch1;
			LepPT      pt1;
			LepAbsRap  absy1;
			ALpCh      ch2;
			ALpPT      pt2;
			ALpAbsRap  absy2;
			BosMT      mt;
			MET        etmiss;
			//absolute-rapidity-ordered leptons
			LepCPT lcpt;
			LepFPT lfpt;
			LepCAbsEta lcy;
			LepFAbsEta lfy;
			//CS frame
			CosThCS cthCS;
			//PhiCS phiCS;

			// define here cut decision

			bool operator()(){

				/*
				 *	      if (opts.lptcut > 0){
				 *		if (ch1()!=0 && pt1() < opts.lptcut) return SkipEvent;
				 *		if (ch2()!=0 && pt2() < opts.lptcut) return SkipEvent;
			}
			if (opts.lycut < 100){
				if (ch1()!=0 && absy1() > opts.lycut) return SkipEvent;
				if (ch2()!=0 && absy2() > opts.lycut) return SkipEvent;
			}
			if (opts.mtcut     >0   && mt()     < opts.mtcut     ) return SkipEvent;
			if (opts.etmisscut >0   && etmiss() < opts.etmisscut ) return SkipEvent;
			*/

				if (opts.lepptcut   >0   && pt1()    < opts.lepptcut   ) return SkipEvent;
				if (opts.alpptcut   >0   && pt2()    < opts.alpptcut   ) return SkipEvent;
				if (opts.lepycut    <100 && absy1()  > opts.lepycut    ) return SkipEvent;
				if (opts.alpycut    <100 && absy2()  > opts.alpycut    ) return SkipEvent;
				if ((opts.lcymax    <100 || opts.lcymin    > 0) && (lcy()  > opts.lcymax || lcy()  < opts.lcymin)    ) return SkipEvent;
				if ((opts.lfymax    <100 || opts.lfymin    > 0) && (lfy()  > opts.lfymax || lfy()  < opts.lfymin)    ) return SkipEvent;
				if (opts.lcptcut   >0   && lcpt()    < opts.lcptcut   ) return SkipEvent;
				if (opts.lfptcut   >0   && lfpt()    < opts.lfptcut   ) return SkipEvent;

				if ((opts.cthCSmin  > -1 || opts.cthCSmax  < 1) && (cthCS() < opts.cthCSmin || cthCS() > opts.cthCSmax)) return SkipEvent;


				return KeepEvent;
			}
		} standard_cuts;

		UserCuts user_cuts;

		bool KeepThisEvent(double _p3[4], double _p4[4]){
			//SetKinematics(_p3,_p4);

			isIntegratorMode=false;
			if (!opts.makecuts) return KeepEvent;
			event_weight=1.;

			//rewritten cuts in plain form since they are faster
			p3[0]=_p3[0];
			p3[1]=_p3[1];
			p3[2]=_p3[2];
			p3[3]=_p3[3];
			p4[0]=_p4[0];
			p4[1]=_p4[1];
			p4[2]=_p4[2];
			p4[3]=_p4[3];

			double ptsq;
			double ptsqh;
			double ptsqs;
			double Q2;

			double exp2y;
			// -- Photon cuts
			if (opts.nproc == 4 && opts.phpt_cut_h > 0.)
			{
				ptsqh = pow(p3[0],2)+pow(p3[1],2);
				ptsqs = pow(p4[0],2)+pow(p4[1],2);
				if (ptsqh  <  ptsqs ) swap(fTmin12,fTmin22);
				if (ptsqh  <  fTmin12 || ptsqs < fTmin22) return  SkipEvent;
			}

			if (opts.nproc == 4 && opts.phycut > 0.)
			{
				// -- first photon
				exp2y = (p3[3] + fabs(p3[2])) / (p3[3] - fabs(p3[2]));
				if (exp2y > exp2phycut) return SkipEvent;
				// -- Crack detector
				if(opts.apply_crack)
					if (exp2y < exp2crack_eta_max && exp2y > exp2crack_eta_min) return SkipEvent;
					// -- second photon
					exp2y = (p4[3] + fabs(p4[2])) / (p4[3] - fabs(p4[2]));
				if (exp2y > exp2phycut) return SkipEvent;
				// -- Crack detector
				if(opts.apply_crack)
					if (exp2y < exp2crack_eta_max && exp2y > exp2crack_eta_min) return SkipEvent;
					Q2 = 2.*(p3[0]*p4[0]+p3[1]*p4[1]+p3[2]*p4[2]-p3[3]*p4[3]);
				// cout << " Q from Kinematics : " << sqrt(Q2) << endl;
				// if( sqrt(Q2) <  phasespace::mmin || sqrt(Q2) >  phasespace::mmax ) return SkipEvent;
			}

			// -- lepton cuts
			if (opts.lptcut > 0 && opts.nproc != 4)
				if (opts.nproc == 3)
				{
					ptsq = pow(p3[0],2)+pow(p3[1],2);
					if (ptsq < lpt2cut) return SkipEvent;
					ptsq = pow(p4[0],2)+pow(p4[1],2);
					if (ptsq < lpt2cut)  return SkipEvent;
				}

				else {
					ptsq = (opts.nproc==1 ? pow(p4[0],2)+pow(p4[1],2) : pow(p3[0],2)+pow(p3[1],2));
					if (ptsq < lpt2cut)  return SkipEvent;
				}

				if (opts.lycut < 100)
					if (opts.nproc == 3 && opts.nproc != 4)
					{
						exp2y = (p3[3] + fabs(p3[2])) / (p3[3] - fabs(p3[2]));
						if (exp2y > exp2lycut) return SkipEvent;
						exp2y = (p4[3] + fabs(p4[2])) / (p4[3] - fabs(p4[2]));
						if (exp2y > exp2lycut) return SkipEvent;
					}
					else
					{
						exp2y =  (opts.nproc==1 ? (p4[3] + fabs(p4[2])) / (p4[3] - fabs(p4[2])) : (p3[3] + fabs(p3[2])) / (p3[3] - fabs(p3[2])));
						if (exp2y > exp2lycut) return SkipEvent;
					}


					if (opts.nproc != 3 && opts.mtcut > 0 && opts.nproc != 4)
					{
						double mt2 = 2 * (sqrt(p3[0]*p3[0]+p3[1]*p3[1])*sqrt(p4[0]*p4[0]+p4[1]*p4[1]) - (p3[0]*p4[0]+p3[1]*p4[1]));
						if (mt2 < mt2cut) return SkipEvent;
					}

					if (opts.nproc != 3 && opts.etmisscut > 0 && opts.nproc != 4)
					{
						ptsq = (opts.nproc==1 ? pow(p3[0],2)+pow(p3[1],2) : pow(p4[0],2)+pow(p4[1],2));
						if (ptsq < et2misscut) return SkipEvent;
					}


					if (opts.makecuts){
						if ( standard_cuts() == SkipEvent ) return SkipEvent ;
						if ( user_cuts() == SkipEvent ) return SkipEvent ;
					}
					return KeepEvent;
		}
	}
}


#endif /* ifndef Kinematics_CXX */

