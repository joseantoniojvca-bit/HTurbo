#include "ctmellin.h"

#include "qtint.h"
#include "mesq.h"
#include "anomalous.h"
#include "pegasus.h"
#include "pmom.h"
#include "ccoeff.h"
#include "pdfevol.h"
#include "evolnative.h"
#include "parton.h"
#include "settings.h"
#include "omegaintegr.h"
#include "rapint.h"
#include "pdf.h"
#include "scales.h"
#include "alphaqed.h"
#include "phasespace.h"
#include "abint.h"
#include "resconst.h"
#include "dyres_interface.h"
#include "switch.h"
#include "isnan.h"

#include <iostream>
#include <string.h>
#include <math.h>

using namespace std;

//using namespace anomalous;
using namespace pmom;


//Counterterm to be subtracted from V+j to get a finite cross section at qt->0
void ctmellin::calc(double m, double qt, double y, int mode, double f[])
{
  //mode 0: differential in m,y,pt,costh
  //mode 1: integrated in costh
  //mode 2: integrated in costh and pt
  //mode 3: integrated in costh, y, and pt

  for (int npdf = 0; npdf < opts.totpdf; npdf++)
    f[npdf] = 0.;
      

  ///////////////////////////////////////////////////////
  double m2 = m*m;
  //  double qt2 = qt*qt;
  //  double exppy = exp(y);
  //  double expmy = 1./exppy;

  //Set scales
  scales::set(m);
  scales::mcfm();
  scales::dyres(m);
  double muren = scales::ren;
  double mufac = scales::fac;
  double mures = scales::res;
  //double muf = scales::fac;
  //double mur = scales::ren;

  //update PDFs in Mellin space at the starting scale, if the factorisation scale is proportional to mll
  //!!! create a special function/module to deal with non-evolved PDFs
  if (opts.melup >= 1)
    {
      mellinint::allocate();
      ccoeff::allocate();
      //pegasus::allocate();
      pmom::allocate();

      mellinint::updategauss();
      ccoeff::calc();
      pegasus::calc_mellin();
      pmom::calc();
    }
  evolnative::allocate();
  evolnative::update();
  pdfevol::allocate_fx();
  
  if (!opts.mellin1d && mode == 3)
    {
      double ylim = 0.5*log(pow(opts.sroot,2)/phasespace::m2);
      double ymn = min(max(-ylim, phasespace::ymin),ylim);
      double ymx = max(min(ylim, phasespace::ymax),-ylim);
      rapint::allocate();
	
      if (opts.makecuts)
	{
	  //there is a potential issue here, when lepton cuts are applied
	  //the rapidity dependent exponential are cached assuming integration between ymin and ymax
	  //for consistency, has to keep the integration between ymin and ymax
	  //rapint::integrate(phasespace::ymin,phasespace::ymax,phasespace::m);
	  rapint::numint(ymn,ymx,phasespace::m);
	}
      else
	rapint::integrate(ymn,ymx,phasespace::m);
    }
  

  /*
  //Set factorization scale
  double muf, mur;
  if (opts.dynamicscale)
    {
      //Cannot use dynamic scale, throw an error!!!
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

  //a-parameter of the resummation scale, set it for the dynamic case
  //if (opts.fmures > 0)
  //a_param_.a_param_ = 1./opts.kmures;
  
  //for fixed resummation scale need to recompute a_param
  //else
  //a_param_.a_param_ = m/scales::res;
  //////////////////////////////////////////////////////////

  double LR, LF, LQ;
  if (opts.order_ct >= 2)
    LR = log(m2/pow(muren,2));
  LF = log(m2/pow(mufac,2));
  //LQ = 2.*log(a_param_.a_param_);
  LQ = log(m2/pow(mures,2));

  // skip PDF loop in the preconditioning phase
  int maxpdf=0;
  if (dofill_.doFill_ != 0) maxpdf = opts.totpdf;
      
  // Start calculation
  double H1q = 0.;

  //factors to compensate the as/2/pi normalization (a factor of 2 for each power of alphas)
  double beta0N = resconst::beta0*2.;
  double beta1N = resconst::beta1*4.;
  double beta2N = resconst::beta2*8.;
  double A1gN = resconst::A1g*2.;
  double A2gN = resconst::A2g*4.;
  double A3gN = resconst::A3g*8.;
  double B1gN = resconst::B1g*2.;
  double B2gN = resconst::B2g*4.;
  double A1qN = resconst::A1q*2.;
  double A2qN = resconst::A2q*4.;
  double A3qN = resconst::A3q*8.;
  double B1qN = resconst::B1q*2.;
  double B2qN = resconst::B2q*4.;

  //QED
  double A1qpN = resconst::A1qp*2.;
  double A2qpN = resconst::A2qp*4.;
  double B1qpN = resconst::B1qp*2.;

  double asopi = pdf::alphas(scales::ren)/M_PI;

  double qedfacfs = resconst::a0qed/M_PI;
  double qedfacis = (alphaqed::calc(pow(scales::ren,2),opts.qedorder+1)*4.*M_PI)/M_PI; //Use alpha(mur) QED for initial state photons

  double bjx = m2/pow(opts.sroot,2);
  double ax = log(bjx);
  double ax1 = (ax+2*y)/2.; 
  double ax2 = (ax-2*y)/2.; 
  
  //double lumi[mesq::totpch];
  double sig11[mesq::totpch];
  double sig12[mesq::totpch];
  double sig21[mesq::totpch];
  double sig22[mesq::totpch];
  double sig23[mesq::totpch];
  double sig24[mesq::totpch];
  //QED
  double sig11_qed[mesq::totpch];
  double sig12_qed[mesq::totpch];
  fill(sig11,sig11+mesq::totpch, 0.);
  fill(sig12,sig12+mesq::totpch, 0.);
  fill(sig21,sig21+mesq::totpch, 0.);
  fill(sig22,sig22+mesq::totpch, 0.);
  fill(sig23,sig23+mesq::totpch, 0.);
  fill(sig24,sig24+mesq::totpch, 0.);
  fill(sig11_qed,sig11_qed+mesq::totpch, 0.);
  fill(sig12_qed,sig12_qed+mesq::totpch, 0.);
  
  //Start mellin loop
  if (opts.mellin1d && mode == 3)
    for (int n = 0; n < mellinint::mdim; n++)
      {
	//retrieve PDFs in Mellin space at the factorisation scale
	pdfevol::retrievemuf_1d(n);
	complex<double>* fn1 = pdfevol::fn1;
	complex<double>* fn2 = pdfevol::fn2;

	complex<double> fn1A = pdfevol::fotn1;
	complex<double> fn2A = pdfevol::fotn2;
      
	int ni = anomalous::index(n,mesq::positive);

	//Mellin transform
	complex <double> cexp = exp(-mellinint::Np[n] * ax)/M_PI * mellinint::CCp/complex <double>(0.,1);

	//cout << mellinint::mdim << "  " << n << "  " << mellinint::Np[n] << "  " << fn1[4] << "  " << fn2[4] << endl;
      
	//loop on born subprocesses, i.e. born incoming partons ij
	complex<double> s11,s12,s21,s22,s23,s24;
	complex<double> s11_qed,s12_qed;
	  for (int sp = 0; sp < mesq::totpch; sp++)
	  {
	    //simplify notation
	    //double bornmesqij = real(mesq::mesqij[sp]); //born level amplitudes
	    parton::partid i = mesq::pidn1[sp];         //parton 1
	    parton::partid j = mesq::pidn2[sp];         //parton 2
	    parton::partid g = parton::g;              //gluon
	    parton::partid im = parton::charge_conjn(i);
	    parton::partid jm = parton::charge_conjn(j);

	    //LO term (there is no counterterm at LO...)
	    //Simplest term without convolutions
	    complex <double> tdelta = fn1[i]*fn2[j];
	  
	    complex <double> th1stF = 2.*(2.*pmom::gamma1qq[ni])*tdelta + (2.*pmom::gamma1qg[ni])*(fn1[i]*fn2[g]+fn1[g]*fn2[j]);

	    s12 = -0.5*A1qN*tdelta;
	    s11 = -(B1qN+A1qN*LQ)*tdelta - th1stF;

	    //QED
	    if (opts.qed && opts.qedorder_ct >= 1)
	      {
		s12_qed =  qedfacfs*(-0.5*A1qpN*tdelta)*(parton::chsq(i)+parton::chsq(j))/2.;
		s11_qed =  qedfacfs*(-(B1qpN+A1qpN*LQ+2.*2.*pmom::gamma1qqqed[ni])*tdelta)*(parton::chsq(i)+parton::chsq(j))/2.;
		s11_qed += qedfacis*(-2.*pmom::gamma1qgamma[ni]*(fn1[i]*fn2A*parton::chsq(i)+fn1A*fn2[j]*parton::chsq(j))); //photon PDF
	      }
	  
	    sig11[sp] += real(s11 * cexp * mellinint::wn[n]);
	    sig12[sp] += real(s12 * cexp * mellinint::wn[n]);
	    //sig11[sp] += real(tdelta * cexp * mellinint::wn[n])*real(mesq::mesqij[sp]);
	    sig11_qed[sp] += real(s11_qed * cexp * mellinint::wn[n]);
	    sig12_qed[sp] += real(s12_qed * cexp * mellinint::wn[n]);

	  
	    //end NLO
	    if (opts.order_ct == 1) continue;
	  
	  
	    complex <double> th1stQ = -(B1qN+A1qN*LQ/2.)*LQ*tdelta;
	    complex <double> th1st = 2.*2.*ccoeff::C1qq[ni]*tdelta + 2.*ccoeff::C1qg[ni]*(fn1[i]*fn2[g]+fn1[g]*fn2[j]);
	  
	    s24 = pow(A1qN,2)/8.*tdelta;
	    s23 = -beta0N*A1qN/3.*tdelta-0.5*A1qN*s11;
	   
	    complex <double> tgaga = 
	      +4.*pow(2.*pmom::gamma1qq[ni],2)*tdelta
	      +3.*2.*pmom::gamma1qq[ni]*2.*pmom::gamma1qg[ni]*(fn1[i]*fn2[g]+fn1[g]*fn2[j])
	      +2.*pow(2.*pmom::gamma1qg[ni],2)*(fn1[g]*fn2[g])
	      +2.*pmom::gamma1qg[ni]*2.*pmom::gamma1gg[ni]*(fn1[i]*fn2[g]+fn1[g]*fn2[j]);

	    for (int k = 0; k < 2*MAXNF+1; k++)
	      if (k != g)
		tgaga += 2.*pmom::gamma1qg[ni]*2.*pmom::gamma1gq[ni]*(fn1[i]*fn2[k]+fn1[k]*fn2[j]);

	    s22 =
	      0.5*(beta0N*A1qN*(LR-LQ)-A2qN)*tdelta
	      -0.5*A1qN*(H1q*tdelta+th1st+th1stQ+(LF-LQ)*th1stF)
	      -0.5*(B1qN+A1qN*LQ-beta0N)*s11
	      +0.5*(B1qN+A1qN*LQ)*th1stF
	      +0.5*tgaga;

	    //	  //add this piece to tgaga
	    //	  for (int k = 0; k < 2*MAXNF+1; k++)
	    //	    if (k != g)
	    //	      s22 += 0.5*2.*pmom::gamma1qg[ni]*2.*pmom::gamma1gq[ni]*(fn1[i]*fn2[k]+fn1[k]*fn2[j]);

	    s21 =
	      -beta0N*(LR-LQ)*s11
	      -(B1qN+A1qN*LQ)*(H1q*tdelta+th1stQ+th1st+(LF-LQ)*th1stF)
	      -(LF-LQ)*tgaga
	      -(B2qN+A2qN*LQ)*tdelta
	      +beta0N*th1st
	      +(B1qN+A1qN*LQ/2.)*LQ*th1stF
	      - (2.*ccoeff::C1qg[ni]*(2.*pmom::gamma1qq[ni]+2.*pmom::gamma1gg[ni])+(H1q+2.*2.*ccoeff::C1qq[ni])*2.*pmom::gamma1qg[ni])*(fn1[i]*fn2[g]+fn1[g]*fn2[j])  //tcga
	      -(2.*(4.*pmom::gamma2qq[ni]))*tdelta //tgamma2
	      -(4.*pmom::gamma2qg[ni])*(fn1[i]*fn2[g]+fn1[g]*fn2[j])
	      -(4.*pmom::gamma2qqb[ni])*(fn1[i]*fn2[jm]+fn1[im]*fn2[j]) //tgamma2
	      -2.*2.*ccoeff::C1qg[ni]*2.*pmom::gamma1qg[ni]*(fn1[g]*fn2[g]) //tgamma2
	      -2.*H1q*2.*pmom::gamma1qq[ni]*tdelta
	      -4.*2.*ccoeff::C1qq[ni]*2.*pmom::gamma1qq[ni]*tdelta;

	    for (int k = 0; k < 2*MAXNF+1; k++)
	      {
		if (k != g && k != j && k != jm)
		  s21 -= 4.*pmom::gamma2qqp[ni]*(fn1[i]*fn2[k]);
		if (k != g && k != i && k != im)
		  s21 -= 4.*pmom::gamma2qqp[ni]*(fn1[k]*fn2[j]);
	      }


	    for (int k = 0; k < 2*MAXNF+1; k++)
	      if (k != g)
		s21 -= 2.*ccoeff::C1qg[ni]*2.*pmom::gamma1gq[ni]*(fn1[i]*fn2[k]+fn1[k]*fn2[j]);
	  
	    sig21[sp] += real(s21 * cexp * mellinint::wn[n]);
	    sig22[sp] += real(s22 * cexp * mellinint::wn[n]);
	    sig24[sp] += real(s24 * cexp * mellinint::wn[n]);
	    sig23[sp] += real(s23 * cexp * mellinint::wn[n]);
	  }
      }
  else
    for (int sign = mesq::positive; sign <= mesq::negative; sign++)
      for (int i1 = 0; i1 < mellinint::mdim; i1++)
	for (int i2 = 0; i2 < mellinint::mdim; i2++)
	  {
	    int ii1 = anomalous::index(i1,mesq::positive);
	    int ii2 = anomalous::index(i2,sign);

	    //retrieve PDFs in Mellin space at the factorisation scale
	    pdfevol::retrievemuf(i1,i2,sign);

	    complex<double>* fn1 = pdfevol::fn1;
	    complex<double>* fn2 = pdfevol::fn2;

	    complex<double> fn1A = pdfevol::fotn1;
	    complex<double> fn2A = pdfevol::fotn2;

	    complex <double> cexp;

	    if (mode == 3)
	      if ( sign == mesq::positive)
		cexp = -rapint::Ith0p[mellinint::index(i1,i2)]/2.;
	      else
		cexp = rapint::Ith0m[mellinint::index(i1,i2)]/2.;
	    else
	      {
		complex <double> cex1, cex2;
		cex1 = exp(-mellinint::Np_1[i1] * ax1) / M_PI * mellinint::CCp;
		if ( sign == mesq::positive)
		  cex2 = exp(-mellinint::Np_2[i2] * ax2) / M_PI * mellinint::CCp;
		else
		  cex2 = -exp(-mellinint::Nm_2[i2] * ax2) / M_PI * mellinint::CCm; //Is this the complex conjugate of the above?
		cexp = -cex1*cex2*mellinint::wn_1[i1]*conj(mellinint::wn_2[i2]);
	      }
	    
	    complex<double> s11,s12,s21,s22,s23,s24 = 0.;
	    complex<double> s11_qed,s12_qed = 0.;

	    for (int sp = 0; sp < mesq::totpch; sp++)
	      {
		parton::partid i = mesq::pidn1[sp];         //parton 1
		parton::partid j = mesq::pidn2[sp];         //parton 2
		parton::partid g = parton::g;              //gluon
		parton::partid im = parton::charge_conjn(i);
		parton::partid jm = parton::charge_conjn(j);

		//LO term (there is no counterterm at LO...)
		//Simplest term without convolutions
		complex <double> tdelta = fn1[i]*fn2[j];
		complex <double> th1stF = 2.*(pmom::gamma1qq_1[ii1]+pmom::gamma1qq_2[ii2])*tdelta
		  + 2.*(pmom::gamma1qg_2[ii2]*fn1[i]*fn2[g])
		  + 2.*(pmom::gamma1qg_1[ii1]*fn1[g]*fn2[j]);

		s12 = -0.5*A1qN*tdelta;
		s11 = -(B1qN+A1qN*LQ)*tdelta - th1stF;

		//QED
		if (opts.qed && opts.qedorder_ct >= 1)
		  {
		    s12_qed = qedfacfs*(-0.5*A1qpN*tdelta)*(parton::chsq(i)+parton::chsq(j))/2.;
		    s11_qed = qedfacfs*(-(B1qpN+A1qpN*LQ+2.*pmom::gamma1qqqed_1[ii1])*tdelta)*(parton::chsq(i)+parton::chsq(j))/2.;
		    s11_qed += qedfacfs*(-(B1qpN+A1qpN*LQ+2.*pmom::gamma1qqqed_2[ii2])*tdelta)*(parton::chsq(i)+parton::chsq(j))/2.;
		    s11_qed += qedfacis*(-2.*pmom::gamma1qgamma_1[ii1]*(fn1A*fn2[j]*parton::chsq(j))); //photon PDF
		    s11_qed += qedfacis*(-2.*pmom::gamma1qgamma_2[ii2]*(fn1[i]*fn2A*parton::chsq(i))); //photon PDF
		  }
		
		sig11[sp] += real(s11 * cexp)/2.;
		sig12[sp] += real(s12 * cexp)/2.;
		sig11_qed[sp] += real(s11_qed * cexp)/2.;
		sig12_qed[sp] += real(s12_qed * cexp)/2.;
		
		//end NLO
		if (opts.order_ct == 1) continue;

		complex <double> th1stQ = -(B1qN+A1qN*LQ/2.)*LQ*tdelta;
		complex <double> th1st = 2.*(ccoeff::C1qq_1[ii1]+ccoeff::C1qq_2[ii2])*tdelta
		  + 2.*ccoeff::C1qg_1[ii1]*fn1[g]*fn2[j]
		  + 2.*ccoeff::C1qg_2[ii2]*fn1[i]*fn2[g];
	  
		s24 = pow(A1qN,2)/8.*tdelta;
		s23 = -beta0N*A1qN/3.*tdelta-0.5*A1qN*s11;

		complex <double> tgaga = 
 		  +4.*pmom::gamma1qq_1[ii1]*pmom::gamma1qq_1[ii1]*tdelta
		  +4.*pmom::gamma1qq_1[ii1]*pmom::gamma1qq_2[ii2]*tdelta
		  +4.*pmom::gamma1qq_2[ii2]*pmom::gamma1qq_1[ii1]*tdelta
		  +4.*pmom::gamma1qq_2[ii2]*pmom::gamma1qq_2[ii2]*tdelta
		  +4.*2.*pmom::gamma1qq_1[ii1]*pmom::gamma1qg_2[ii2]*(fn1[i]*fn2[g])
		  +2.*2.*pmom::gamma1qq_2[ii2]*pmom::gamma1qg_2[ii2]*(fn1[i]*fn2[g])
		  +2.*2.*pmom::gamma1qq_1[ii1]*pmom::gamma1qg_1[ii1]*(fn1[g]*fn2[j])
		  +4.*2.*pmom::gamma1qq_2[ii2]*pmom::gamma1qg_1[ii1]*(fn1[g]*fn2[j])
		  +2.*4.*pmom::gamma1qg_1[ii1]*pmom::gamma1qg_2[ii2]*(fn1[g]*fn2[g])
		  +2.*pmom::gamma1qg_2[ii2]*2.*pmom::gamma1gg_2[ii2]*(fn1[i]*fn2[g])
		  +2.*pmom::gamma1qg_1[ii1]*2.*pmom::gamma1gg_1[ii1]*(fn1[g]*fn2[j])
		  ;

		for (int k = 0; k < 2*MAXNF+1; k++)
		  if (k != g)
		    tgaga +=
		      +2.*pmom::gamma1qg_2[ii2]*2.*pmom::gamma1gq_2[ii2]*(fn1[i]*fn2[k])
		      +2.*pmom::gamma1qg_1[ii1]*2.*pmom::gamma1gq_1[ii1]*(fn1[k]*fn2[j]);

		s22 =
		  0.5*(beta0N*A1qN*(LR-LQ)-A2qN)*tdelta
		  -0.5*A1qN*(H1q*tdelta+th1st+th1stQ+(LF-LQ)*th1stF)
		  -0.5*(B1qN+A1qN*LQ-beta0N)*s11
		  +0.5*(B1qN+A1qN*LQ)*th1stF
		  +0.5*tgaga;

		complex <double> tgamma2, tcga;
		tgamma2 = 
		  +(2.*(2.*(pmom::gamma2qq_1[ii1]+pmom::gamma2qq_2[ii2])))*tdelta
		  +(4.*pmom::gamma2qg_2[ii2])*(fn1[i]*fn2[g])
		  +(4.*pmom::gamma2qg_1[ii1])*(fn1[g]*fn2[j])
		  +(4.*pmom::gamma2qqb_2[ii2])*(fn1[i]*fn2[jm])
		  +(4.*pmom::gamma2qqb_1[ii1])*(fn1[im]*fn2[j]);

		for (int k = 0; k < 2*MAXNF+1; k++)
		  {
		    if (k != g && k != j && k != jm)
		      tgamma2 += 4.*pmom::gamma2qqp_2[ii2]*(fn1[i]*fn2[k]);
		    if (k != g && k != i && k != im)
		      tgamma2 += 4.*pmom::gamma2qqp_1[ii1]*(fn1[k]*fn2[j]);
		  }

		tcga = 
		  +2.*ccoeff::C1qg_2[ii2]*2.*pmom::gamma1gg_2[ii2]*(fn1[i]*fn2[g])
		  +2.*ccoeff::C1qg_1[ii1]*2.*pmom::gamma1gg_1[ii1]*(fn1[g]*fn2[j])
		  +2.*ccoeff::C1qq_2[ii2]*2.*pmom::gamma1qg_2[ii2]*(fn1[i]*fn2[g])
		  +2.*ccoeff::C1qq_1[ii1]*2.*pmom::gamma1qg_1[ii1]*(fn1[g]*fn2[j])
		  +4.*ccoeff::C1qq_1[ii1]*pmom::gamma1qq_1[ii1]*tdelta
		  +4.*ccoeff::C1qq_2[ii2]*pmom::gamma1qq_2[ii2]*tdelta
		  +4.*ccoeff::C1qg_2[ii2]*pmom::gamma1qq_1[ii1]*(fn1[i]*fn2[g])
		  +4.*ccoeff::C1qg_1[ii1]*pmom::gamma1qq_2[ii2]*(fn1[g]*fn2[j])
		  +4.*ccoeff::C1qq_1[ii1]*pmom::gamma1qg_2[ii2]*(fn1[i]*fn2[g])
		  +4.*ccoeff::C1qq_2[ii2]*pmom::gamma1qg_1[ii1]*(fn1[g]*fn2[j])
		  +2.*ccoeff::C1qg_1[ii1]*2.*pmom::gamma1qg_2[ii2]*(fn1[g]*fn2[g])
		  +2.*ccoeff::C1qg_2[ii2]*2.*pmom::gamma1qg_1[ii1]*(fn1[g]*fn2[g])
		  +4.*ccoeff::C1qq_1[ii1]*pmom::gamma1qq_2[ii2]*tdelta
		  +4.*ccoeff::C1qq_2[ii2]*pmom::gamma1qq_1[ii1]*tdelta;

		for (int k = 0; k < 2*MAXNF+1; k++)
		  if (k != g)
		    tcga +=
		      +2.*ccoeff::C1qg_2[ii2]*2.*pmom::gamma1gq_2[ii2]*(fn1[i]*fn2[k])
		      +2.*ccoeff::C1qg_1[ii1]*2.*pmom::gamma1gq_1[ii1]*(fn1[k]*fn2[j]);

		s21 =
		  -beta0N*(LR-LQ)*s11
		  -(B1qN+A1qN*LQ)*(H1q*tdelta+th1stQ+th1st+(LF-LQ)*th1stF)
		  -(LF-LQ)*tgaga
		  -(B2qN+A2qN*LQ)*tdelta
		  +beta0N*th1st
		  +(B1qN+A1qN*LQ/2.)*LQ*th1stF
		  -tgamma2
		  -tcga;
		
		sig21[sp] += real(s21 * cexp)/2.;
		sig22[sp] += real(s22 * cexp)/2.;
		sig24[sp] += real(s24 * cexp)/2.;
		sig23[sp] += real(s23 * cexp)/2.;
	      }
	  }
  
  double xmsq = 0.;
  for (int sp = 0; sp < mesq::totpch; sp++)
    {
      //as/pi factor
      double sig1 = (sig12[sp]*qtint::LL2_mesqij[sp]+sig11[sp]*qtint::LL1_mesqij[sp])*asopi/2. + (sig12_qed[sp]*qtint::LL2_mesqij[sp]+sig11_qed[sp]*qtint::LL1_mesqij[sp])/2.;
      xmsq += -sig1;
      
      if (opts.order_ct == 1) continue;
      //(as/pi)^2 factor
      double sig2 = (sig24[sp]*qtint::LL4_mesqij[sp]+sig23[sp]*qtint::LL3_mesqij[sp]+sig22[sp]*qtint::LL2_mesqij[sp]+sig21[sp]*qtint::LL1_mesqij[sp])*pow(asopi/2.,2);

      if (opts.cmwfo) sig2 += -resconst::K1*sig1*asopi;
      
      // sum O(as) and O(as^2) contributions	
      xmsq += -sig2;
    }	  

  //Do not fully understand this factor
  xmsq = xmsq * 3./8. /2./M_PI;

  //m^2/s factor
  xmsq = xmsq * m2/pow(opts.sroot,2);

  //phiV integration  
  xmsq = xmsq * 2*M_PI;

  if (isnan_ofast(xmsq))
    cout << m << " " << xmsq << endl;
  
  //switching function is inside qtint, do not apply
  
  f[0] = xmsq;

  //cout << "xmsq " << xmsq << "  " << asopi << endl;

  if (opts.melup == 2)
    {
      mellinint::free();
      ccoeff::free();
      //pegasus::free();
      pmom::free();
    }
  evolnative::free();
  
  pdfevol::free_fx();

  if (!opts.mellin1d && mode == 3)
    rapint::free();
  
  return;
}
