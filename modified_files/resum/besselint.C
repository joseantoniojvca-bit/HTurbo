#include "resconst.h"
#include "settings.h"
#include "besselint.h"
#include "pdfevol.h"
#include "pegasus.h"
#include "mesq.h"
#include "hcoefficients.h"
#include "hcoeff.h"
#include "expc.h"
#include "cevol.h"
#include "muf.h"
#include "phasespace.h"
#include "gaussrules.h"
#include "resint.h"
#include "isnan.h"
#include "pdf.h"
#include "sudakovff.h"
#include "npff.h"
#include "ccoeff.h"
#include "gint.h"
#include "scales.h"

#include <LHAPDF/LHAPDF.h>

double C1 = 1;
double C3 = 1;

//fortran interface
void besselint_(double &b, double &qt, double &q2)
{
  resint::_qt = qt;
  resint::_m = sqrt(q2);
  besselint::bint(b);
};

complex <double> besselint::bint(complex <double> b)
{
  if (b == 0.)
    //if (opts.modlog && abs(b) < 1e-8)
    return 1.; //Should be correct, since J_0(0) = 1, and S(0) = 1 (Sudakov). May be different in the case that the L~ = L+1 matching is not used.
  
  double qt = resint::_qt; //better take this from phasespace::qt

  /***************************************************** Moved to scales (begin)
  //pdfevol::bscale is used for
  //alpq (C1)
  //pdfevol::alpr (C3)
  //pdfevol::XL (C3)
  //aexp (C1)
  //it also used in evolmode 1 and 4 as factorisation scale
  

  //The soft scale of the resummation is b0/b. The factor a = mll/mures is introduced because the b scale is used in alphasl
  //Introduce two scales:
  //mub = C1*b0/b is the lower integration limit of the Sudakov integral
  complex <double> mub = resint::a * resconst::b0/b * C1;
  //pdfevol::bscale = C3*b0/b is the "soft" factorisation scale at which the Wilson coefficient functions are evaluated (see arxiv:1309.1393 for details)
  pdfevol::bscale = resint::a * resconst::b0/b * C3;
  
  //convert to fortran complex number
  fcomplex fscale2_mub = fcx(pow(mub,2));
  fcomplex fscale2_mufb = fcx(pow(pdfevol::bscale,2));

  //freeze PDF evolution below a certain scale
  //  if (pdfevol::bscale < 5.)
  //    pdfevol::bscale = 5.;

  //pdfevol::bstarscale and pdfevol::bstartilde are not used (be careful, they do no have the C3 scale)

  //bstarscale = a*b0/bstar (final scale used for the PDF evolution)
  double bstar = real(b / sqrt(1.+(b*b)/(blimit_.rblim_*blimit_.rblim_)));
  //double bstar = real(b / sqrt(1.+(b*b)/(1.1229190*1.1229190)));  //// -->hard-coded!! --> should allow a different blim in the PDF evolution as a setting
  pdfevol::bstarscale = resconst::b0*resint::a/bstar;

  //qbstar is used in evolmode 3
  //qbstar = b0/bstar, it corresponds to pdfevol::bsstarcale but without a_param
  pdfevol::qbstar = resconst::b0/bstar*C3;

  //simulate pythia ISR factorisation scale, which is muf^2 = qt^2 (i.e. avoid the bstar prescription to freeze the factorisation scale)
  //pdfevol::qbstar = resconst::b0/b*C3;
  
  //bstartilde is bstarscale (qbstar) with the modification L -> L~, which freezes the scale at muf
  //bstartilde is used in evolmode 2, for the direct mellin transfrom at each scale
  //pdfevol::bstartilde = pdfevol::bstarscale * resint::mufac / sqrt((pow(pdfevol::bstarscale,2) + resint::mufac2));
  if (opts.modlog)
    pdfevol::bstartilde = pdfevol::qbstar * resint::mures / sqrt((pow(pdfevol::qbstar,2) + resint::mures2));
  else
    pdfevol::bstartilde = pdfevol::qbstar;


  //complex PDF scale to be used for the minimal prescription
  pdfevol::bcomplex = resconst::b0/b;
  
  //  cout << b << "  " << bstar << "  " << blimit_.rblim_ << endl;
  *****************************************************/ // Moved to pdfevol::scales (end)

  
  //The integration from b to qt space is done with the bstar prescription (real axis in the b space), and use the bessel function
  //The (complex) integration in the minimal prescription would require hankel functions
  //********************
  //qt and b dependence (bessel function) (is xj0 a jacobian? Probably yes, the Jacobian for the change of variable from cartesian to radial coordinate, which translates from Fourier to Hankel transform)
  complex <double> xj0;
  double qtb = qt*real(b);
  if (resint::_mode == 3 || resint::_mode == 4) //qt-integrated mode --> Needs to be implemented also for the other prescriptions
    {
      //bstar prescription
      if (opts.bprescription == 0)
	xj0 = 2.*fort_besj1_(qtb)/b;
      //Integrate up to blim with besche
      else if (opts.bprescription == 1)
	xj0 = 2./b;
      //Minimal prescription
      else if (opts.bprescription == 2 || opts.bprescription == 3 || opts.bprescription == 5)
	xj0 = 2./b;
      //Local bstar prescription
      else if (opts.bprescription == 4)
	xj0 = 2.*fort_besj1_(qtb)/b;
    }
  else
    {
      //bstar prescription
      if (opts.bprescription == 0)
	xj0 = 2.*fort_besj0_(qtb);
	//xj0 = 2.*fort_besj1_(qtb)/real(b); //--> this line is to check the qt-integral as a function of qt used in mode 3
      //Integrate up to blim with besche
      else if (opts.bprescription == 1)
	xj0 = 2.;
      //Minimal prescription
      else if (opts.bprescription == 2 || opts.bprescription == 3 || opts.bprescription == 5)
	xj0 = 2.;
      //Local bstar prescription
      else if (opts.bprescription == 4)
	xj0 = 2.*fort_besj0_(qtb);
    }
  //********************
  //(Eq. 3.8 of https://arxiv.org/pdf/1805.05916.pdf)
  //double Q2 = pow(scales::res,2);
  //complex <double> b2 = pow(b,2);
  //double b02 = pow(resconst::b0,2);
  //if (opts.modlog && opts.p == 1)
  //  xj0 *= sqrt(Q2*b2/b02/(Q2*b2/b02+1.));
  //else if (opts.modlog && opts.p > 1)
  //  xj0 *= sqrt(pow(Q2*b2/b02,opts.p)/(1.+pow(Q2*b2/b02,opts.p)));

  //if (opts.modlog && opts.p == 1)
  //  xj0 *= scales::res/qt/(1.+scales::res/qt);
  //else if (opts.modlog && opts.p > 1)
  //  xj0 *= pow(scales::res/qt,opts.p)/(1.+pow(scales::res/qt,opts.p));
  //********************
  
  //The Sudakov is m and b dependent

  //numerical Sudakov
  if (opts.numsud || opts.numexpc)
    gint::calc(b);
  
  //fortran
  //fcomplex fb = fcx(b);
  //complex <double> sudak=cx(s_(fb));
  //complex <double> invres = complex<double>(0.,0.);

  //Set scales for the PDF evolution (b-dependence)
  //pdfevol::scales(b);
  //pdfevol::alphasl(b); //alphasl gives the LL/NLL/NNLL evolution of alphas from Qres=Q/a_param to q2=b0^2/b^2
  //pdfevol::alphaqedl(b); //alphasl gives the LL/NLL/NNLL evolution of alphas from Qres=Q/a_param to q2=b0^2/b^2

  //C++
  //Calculate Sudakovs
  complex <double> sudak=sudakov::sff(b);
  if (sudak == 0.)
    return 0.;

  //non-perturbative form factor
  //sudak *= exp(-opts.g_param*pow(b,2));
  complex <double> ff = npff::S(b,resint::_m,resint::x1,resint::x2);
  if (!opts.flavour_kt)
    sudak *= ff;

  //********************
  //b, qt and mass dependence
  complex <double> factorfin = b*xj0*sudak;
  //cout << "where is the nan " << b << "  " << xj0 << "  " << sudak << endl;
  //********************

  //Do not need the Mellin transform for the LL case, because the HN coefficient is 1 --> Use PDFs in x space
  if (opts.order == 0 && opts.xspace && (opts.evolmode == 0 || opts.evolmode == 1)) //if (opts.order == 0 && opts.xspace && (opts.evolmode == 1))
    return factorfin;
  /*
  //Perform PDF evolution
  if (!(opts.order == 0 && opts.xspace))
    pdfevol::evolution();

  // Cache the positive and negative branch of coefficients which depend only on one I index
  if (opts.order_expc > 0)
    {
      expc::reset();
      if (opts.expc < 8)
	expc::calc(b);
      if (opts.expc == 8)
	{
	  cevol::calc(b);
	  hcoeff::reset();
	  hcoeff::calcb();
	}
      if (opts.expc == 9)
	{
	  cevol::css(b);
	  hcoeff::reset();
	  hcoeff::calcb();
	}
      if (opts.expc == 10)
	{
	  cevol::ctilde(b);
	  hcoeff::reset();
	  hcoeff::calcb();
	}
    }
    */
    // Set scales for evolution in pdfevol
    //alphasl gives the LL/NLL/NNLL evolution of alphas from Qres=Q/a_param to q2=b0^2/b^2
    pdfevol::alphasl(b);
    pdfevol::scales(b);
      
    //Perform PDF evolution
    pdfevol::evolution();

    // Cache the positive and negative branch of coefficients which depend only on one I index
    expc::reset();
    expc::calc(b);
    muf::reset();
    if (opts.mufvar
      && (opts.kmufac != opts.kmures))     //Save expensive muf b-dependent scale variations when not needed
    muf::calc(b);

    //Apply flavour dependent non-perturbative form factors
    if (opts.flavour_kt)
      pdfevol::flavour_kt();
      
  complex <double> fun = 0.;
      
  mellinint::reset();
      
  //Do not need the Mellin transform for the LL case, because the HN coefficient is 1 --> Use PDFs in x space
  if (opts.order == 0 && opts.xspace && (opts.evolmode == 2 || opts.evolmode == 3 || opts.evolmode == 4))
    //if (opts.order == 0 && opts.xspace && (opts.evolmode == 0 || opts.evolmode == 2 || opts.evolmode == 3 || opts.evolmode == 4))
    {
      double muf;
      
      if (opts.evolmode == 0 || opts.evolmode == 1)
	muf = resint::mufac;
      else
	muf = real(pdfevol::mubstartilde);

	if (resint::_mode < 2 || resint::_mode == 4) //rapidity differential
	  fun = mesq::loxs(resint::x1, resint::x2, muf);
	else                   //rapidity integrated
	  fun = mesq::loxs(resint::tau, muf);      
    }
  else
    {
      if (opts.qed)
	{
	  if (opts.nproc == 3)
	    {
	      //1d mellin
	      if (opts.mellin1d)
		{
		  fun = mellinint::calc1d(resconst::equ2,resconst::eqd2);
		}
	      else
		{
		  fun = mellinint::calc2d(resconst::equ2,resconst::eqd2);
		}
	    }
	  else if (opts.nproc == 2 || opts.nproc == 1)
	    {
	      //1d mellin
	      if (opts.mellin1d)
		// Here no up down separation! Calling of the overloading with (equ2,eqd2)
		fun   = mellinint::calc1d(resconst::equ2,resconst::eqd2);	  
	      // Here the same We want the overloading of the function which is taking both of the squared charges
	      else
		fun = mellinint::calc2d(resconst::equ2,resconst::eqd2);
	    }
	}
      else
	{
      //1d mellin
      if (opts.mellin1d)
        if (opts.mufvar)
          fun = mellinint::calc1d_muf();
      else
        fun = mellinint::calc1d();
      else
        if (opts.mufvar)
          fun = mellinint::calc2d_muf();
      else
        fun = mellinint::calc2d();
	}
    }

    complex <double> invres;

    invres = fun*factorfin;
  
  if (opts.deltaismb)
    invres += mellinint::bbfun*factorfin*(sudakov::S_bb-1.);
      
  if (opts.deltaismc)
    invres += mellinint::ccfun*factorfin*(sudakov::S_cc-1.);
  
  //complex <double> invres = fun*real(factorfin); //with bstar prescription factorfin is real
  //complex <double> invres = factorfin/resint::_m*(8./3.)*pow(opts.sroot,2)/2.; // --> Check unitarity of Sudakov integral
  //complex <double> invres = factorfin;
      
  //  if (factorfin == 0.)
  //    invres = 0.;
      
  //cout << "invres = " << invres << ", qt = " << qt << ", b = "  << b << ", pdf*mesq = " << fun << ", S " << sudakov::S << " bj0 = " << b*xj0 << endl;
  //cout << "invres = " << invres << ", qt = " << qt << ", b = "  << b << ", pdf = " << pdfevol::fx1[6] << ", mesq = " << mesq::mesqij_expy[mesq::index(0,0,0,mesq::positive)] << endl;
  //cout << setprecision(16) << "C++ " << b << "  " << invres << "  " << fun << "  " << factorfin << endl;
  //cout << setprecision(16) << " b " << b << " bstar " << bstar << " besselint " << invres << " J0 " << fort_besj0_(qtb) << " sud " << sudak << " fun " << fun << endl;
  
  //cout << "invres " << b << "  " << fun << "  " << factorfin << "  " << invres << endl;
  if ( isnan_ofast(real(invres)) || isnan_ofast(imag(invres)))
    {
      //cout << "Warning, invres = " << invres << ", qt = " << qt << ", b = "  << b << ", pdf*mesq = " << fun << ", S*bj0 = " << factorfin << endl;
      //cout << "Warning, invres = " << invres << ", qt = " << qt << ", b = "  << b << ", pdf = " << pdfevol::fx1[6] << ", mesq = " << mesq::mesqij_expy[mesq::index(0,0,0,mesq::positive)] << endl;
      invres = 0;
    }

    //cout << "  invres : " <<  invres << endl;
  
  return invres;
}
