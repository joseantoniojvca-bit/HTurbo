#include "npff.h"
#include "settings.h"
#include "resconst.h"
#include "constants.h"
#include "pdfevol.h"
#include "resint.h"
//#include <LHAPDF/LHAPDF.h>
#include "pdf.h"
#include "blim.h"
#include "gaussrules.h"
#include "specialfunctions_interface.h"
#include "scales.h"
#include "bstar.h"
#include "alphas.h"
#include "btilde.h"
#include <iostream>

using namespace constants;
using namespace resconst;

complex <double> npff::uvff;
complex <double> npff::usff;
complex <double> npff::dvff;
complex <double> npff::dsff;
complex <double> npff::ssff;
complex <double> npff::chff;
complex <double> npff::boff;
complex <double> npff::glff;

void npff::init()
{
  //calculate g0(blim) according to Eq. (81) of https://arxiv.org/pdf/1412.3820.pdf
  int rule = 100;

  double q0 = resconst::b0/opts.bref;
  double q1 = resconst::b0/opts.blim;
  double cc = 0.5*(q1+q0);
  double mm = 0.5*(q1-q0);
  double res = 0;
  for (int i = 0; i < rule; i++)
    {
      double mu = cc+mm*gr::xxx[rule-1][i];
      double fac = 2.*constants::CF/M_PI*gr::www[rule-1][i]*mm;
      res += fac*pdf::alphas(mu)/mu;
      //cout << mu << "  " << fac << "  " << pdf::alphas(mu,opts.cmw) << endl;
    }

  opts.g0 += res;
  //cout << " g0 is now " << opts.g0 << endl;
}


//TODO: Separate the gk (Collins-Soper kernel) and gj (TMD) terms


complex <double> npff::S(complex <double> b, double m, double x1, double x2)
{
  complex <double> ff = 1;
  complex <double> gkk = 0.;
  complex <double> gjj = 0.;

  //General formula with TMD gj and CSK gk
  if (opts.npff <= 2)
    {
      gkk = gk(b);
      gjj = gj(b,x1,x2);
      if (!opts.gjxnspace || opts.evolmode != 2 || opts.xspace)
	gjj += gjx(b,x1)+gjx(b,x2);
      ff = exp(-gjj -gkk*log(pow(m/opts.Q0,2)));
    }
  
  //  if (opts.npff == 0)
  //{
  //Gaussian form
  //BLNY (blim = 0.5)
  //g1 = 0.21;  //g2 = 0.68;  //g3 = -0.13;  //Q0 = 3.2;
  //Konychev Nadolsky (blim = 1.5)
  //g1 = 0.20;  //g2 = 0.19;  //g3 = -0.03;  //Q0 = 3.2;  
  //M. Hirai, H. Kawamura, K. Tanaka, https://inspirehep.net/record/1229722 (minimal prescription, MSTW08)
  //g1 = 0.330;  //g2 = 0.066;  //g3 = 0.;  //Q0 = 2.*1.3;
  //g1 = 0.50;   //g2 = 2.0;  //g3 = 0.;  //Q0 = 50;
  //}
  
  //Exponential form
  //if (opts.npff == 1)
  //{
  //}

  if (opts.npff == 2)
    {
      //Logaritmic term as in Eq. (13) of https://arxiv.org/pdf/1406.3073.pdf (See also Eq.(67) of https://arxiv.org/pdf/1412.3820.pdf and Eq. (10a) of https://inspirehep.net/literature/217570)
      if (opts.g2 != 0)
	ff *= exp(-opts.g2*log(1.+pow(real(b)/blim::sudakov,2))*log(m/opts.Q0));
    }

  //Dispersive approach (DMW) https://arxiv.org/abs/hep-ph/9512336  https://iopscience.iop.org/article/10.1088/1126-6708/1999/07/012/pdf
  //# 3: Dokshitzer-Marchesini-Webber    npff(b) = exp(-(a2 * (log(m)+1.) - a2p) * b^2/2);
  //# Dokshitzer-Marchesini-Webber options (npff = 3)
  //a2 = 0.2
  //a2p = 0
  //Also the N-independent part is included in expc
  if (opts.npff == 3)
    ff = exp(-(opts.a2 * (2.*log(m)+1.) - opts.a2p) * pow(real(b),2)/2. * 2.); //one factor for each leg -> multiply by 2
  //ff = 1.;

  //Pavia19 fit
  if (opts.npff == 4)
    {
      //From Table 5 of https://arxiv.org/pdf/1912.07550.pdf
      double lambda = opts.lambda; //0.580;
      double g2 = opts.g2A; //0.036;   //ln(m)-dependent quadratic (Gaussian) term
      double g2B = opts.g2B; //0.012;  //ln(m)-dependent quartic term
      double N1 = opts.NA; //0.625;
      double sigma = opts.sigmaA; //0.370;
      double alpha = opts.alphaA; //0.205;
      double N1B = opts.NB; //0.044;
      double sigmaB = opts.sigmaB; //0.356;
      double alphaB = opts.alphaB; //0.069;
      double Q0 = 1;
      
      //https://arxiv.org/pdf/1912.07550.pdf Eq. (2.39) and (2.40)
      double g1_1  = N1  / (x1*sigma)  * exp(-1./(2.*pow(sigma,2))  * pow(log(x1/alpha),2) );
      double g1B_1 = N1B / (x1*sigmaB) * exp(-1./(2.*pow(sigmaB,2)) * pow(log(x1/alphaB),2) );

      double g1_2  = N1  / (x2*sigma)  * exp(-1./(2.*pow(sigma,2))  * pow(log(x2/alpha),2) );
      double g1B_2 = N1B / (x2*sigmaB) * exp(-1./(2.*pow(sigmaB,2)) * pow(log(x2/alphaB),2) );

      double four = 4.;
      complex <double> ff_1 = (1.-lambda)/(1.+g1_1*pow(b,2)/four) + lambda*exp(-g1B_1*pow(b,2)/four);
      complex <double> ff_2 = (1.-lambda)/(1.+g1_2*pow(b,2)/four) + lambda*exp(-g1B_2*pow(b,2)/four);
      //ff = ff_1*ff_2*pow(exp(-(g2+g2B*pow(b,2)) *log(pow(m/Q0,2)*pow(b,2)/4.)),2);
      ff = ff_1*ff_2*pow(exp(-    (g2+g2B*4.*pow(abs(b),2)/four) * log(pow(m/Q0,2))* pow(abs(b),2)/four),2);
      //ff =     exp(- 2.*(g2+g2B*pow(abs(b),2)) * log(m*m/Q0/Q0  )* pow(abs(b),2)/4.);

    }

  //SV19 fit
  if (opts.npff == 5)
    {
      
      /*
      //x1 and x2 defined as in https://arxiv.org/pdf/1912.06532.pdf Eq.(2.42)
      double exppy = exp(phasespace::y);
      double expmy = 1./exppy;
      tau = sqrt(m*m+pow(phasespace::qt,2))/opts.sroot;
      double xx1 = tau*exppy;
      double xx2 = tau*expmy;
      */
      double xx1 = x1;
      double xx2 = x2;

      //https://arxiv.org/pdf/1912.06532.pdf Table 9
      double lambda1 = 0.224;
      double lambda2 = 9.24;
      double lambda3 = 375.;
      double lambda4 = 2.15;
      double lambda5 = -4.97;
      
      //https://arxiv.org/pdf/1912.06532.pdf Eq.(2.86)
      complex <double> ff_1 = exp((-lambda1*(1.-xx1) + lambda2*xx1+ xx1*(1.-xx1)*lambda5) / sqrt(1.+lambda3*pow(xx1,lambda4)*pow(b,2)) * pow(b,2));
      complex <double> ff_2 = exp((-lambda1*(1.-xx2) + lambda2*xx2+ xx2*(1.-xx2)*lambda5) / sqrt(1.+lambda3*pow(xx2,lambda4)*pow(b,2)) * pow(b,2));
      ff = ff_1*ff_2;
    }

  if (opts.npff == 6)
    ff = exp( -opts.g1/sqrt(1.+max(0.,opts.lambda)*pow(b,2)) * pow(b,2));

  //Reduce the dependence on bmax using https://arxiv.org/pdf/1412.3820.pdf Eq. (67) and Nucl.Phys. B284, 253 (1987), Eq. (10a)
  if (opts.npff == 7)
    {
      //ff = exp(opts.lambda * CF*pdf::alphas(real(pdfevol::mubstar),opts.cmw)/M_PI*log(1.+pow(real(b)/blim::sudakov,2)));
      //ff = exp(-8.*CF*pdf::alphas(real(pdfevol::mubstar),opts.cmw)/M_PI*( log(1.+pow(real(b)/blim::sudakov,2)) - log(1.+pow(real(b)/3.0,2)) ) );
      //ff = exp(-2.* CF*pdf::alphas(real(pdfevol::mubstar),opts.cmw)/M_PI*log(1.+pow(real(b)/blim::sudakov,2)));
      ff = exp(opts.lambda*log(1.+pow(real(b)/blim::sudakov,2)));

      double gk = opts.g0 * (1. - exp(- (CF*pdf::alphas(real(pdfevol::mubstar))*pow(real(b),2))/(M_PI*opts.g0*pow(blim::sudakov,2))));
      ff *= exp(-gk*log(pow(m/opts.Q0,2)));

      //Gaussian term
      if (opts.g1 > 0)
	ff *= exp(-opts.g1 * pow(abs(b),2));

	
      //Quartic term
      if (opts.q > 0)
	ff *= exp(-opts.q * pow(abs(b),4));

      //cout << "Collins-Rogers " << " mubstar " << pdfevol::mubstar << " alphas " << pdf::alphas(real(pdfevol::mubstar),opts.cmw) << " lha alphas " << pdf::extalphas(real(pdfevol::mubstar)) << endl;
    }

  //gk in the CdFG formalism
  if (opts.npff == 8)
    {
      //Perturbative correction S(b)-S(b*) expanded in b^2/blim^2
      const int p = 1;
      
      double x = pow(real(b)/blim::sudakov,2);
      double x2 = pow(x,2);
      double x3 = pow(x,3);
      double as = resint::aass;
      double as2 = pow(as,2);
      double as3 = pow(as,3);
      double as4 = pow(as,4);
      double as5 = pow(as,5);
      double Q = scales::res;
      double LL = real(btilde::blog(b,Q)); //do not use bstar here, but b
      double LL2 = pow(LL,2);
      double LL3 = pow(LL,3);
      double LL4 = pow(LL,4);
      double LL5 = pow(LL,5);
      //double LL = 1.;      

      double f2k = 0.;
      double f4k = 0.;
      double f6k = 0.;

      //LL
      f2k += resconst::A1q*LL*x*as;
      f2k += resconst::A1q*resconst::beta0*LL2*x*as2;
      f2k += resconst::A1q*pow(resconst::beta0,2)*LL3*x*as3;
      f2k += resconst::A1q*pow(resconst::beta0,3)*LL4*x*as4;
      f2k += resconst::A1q*pow(resconst::beta0,4)*LL5*x*as5;

      //NLL
      if (opts.order_sudak >= 1)
	{
	  f2k += resconst::B1q*x*as;
	  f2k += resconst::B1q*resconst::beta0*LL*x*as2;
	  f2k += resconst::B1q*pow(resconst::beta0,2)*LL2*x*as3;
	  f2k += resconst::B1q*pow(resconst::beta0,3)*LL3*x*as4;
	  f2k += resconst::B1q*pow(resconst::beta0,4)*LL4*x*as5;
	  f2k += resconst::A2q*LL*x*as2;
	  f2k += 2.*resconst::A2q*resconst::beta0*LL2*x*as3;
	  f2k += 3.*resconst::A2q*pow(resconst::beta0,2)*LL3*x*as4;
	  f2k += 4.*resconst::A2q*pow(resconst::beta0,3)*LL4*x*as5;
	  f2k += resconst::A1q*resconst::beta1*LL2*x*as3;
	  f2k += 5./2.*resconst::A1q*resconst::beta0*resconst::beta1*LL3*x*as4;
	  f2k += 13./3.*resconst::A1q*pow(resconst::beta0,2)*resconst::beta1*LL4*x*as5;
	}

      //NNLL	 
      if (opts.order_sudak >= 2)
	{
	  f2k += resconst::B2q*x*as2;
	  f2k += 2.*resconst::B2q*resconst::beta0*LL*x*as3;
	  f2k += 3.*resconst::B2q*pow(resconst::beta0,2)*LL2*x*as4;
	  f2k += 4.*resconst::B2q*pow(resconst::beta0,3)*LL3*x*as5;
	  f2k += resconst::B1q*resconst::beta1*LL*x*as3;
	  f2k += 5./2.*resconst::B1q*resconst::beta0*resconst::beta1*LL2*x*as4;
	  f2k += 13./3.*resconst::B1q*pow(resconst::beta0,2)*resconst::beta1*LL3*x*as5;
	  f2k += resconst::A3q*LL*x*as3;
	  f2k += 3.*resconst::A3q*resconst::beta0*LL2*x*as4;
	  f2k += 6.*resconst::A3q*pow(resconst::beta0,2)*LL3*x*as5;
	  f2k += 2.*resconst::A2q*resconst::beta1*LL2*x*as4;
	  f2k += 7.*resconst::A2q*resconst::beta0*resconst::beta1*LL3*x*as5;
	  f2k += resconst::A1q*resconst::beta2*LL2*x*as4;
	  f2k += 3./2.*resconst::A1q*pow(resconst::beta1,2)*LL3*x*as5;
	  f2k += 3.*resconst::A1q*resconst::beta0*resconst::beta2*LL3*x*as5;
	}

      //N3LL
      if (opts.order_sudak >= 3)
	{
	  f2k += resconst::B3q*x*as3;
	  f2k += 3.*resconst::B3q*resconst::beta0*LL*x*as4;
	  f2k += 6.*resconst::B3q*pow(resconst::beta0,2)*LL2*x*as5;
	  f2k += 2.*resconst::B2q*resconst::beta1*LL*x*as4;
	  f2k += 7.*resconst::B2q*resconst::beta0*resconst::beta1*LL2*x*as5;
	  f2k += resconst::B1q*resconst::beta2*LL*x*as4;
	  f2k += 3./2.*resconst::B1q*pow(resconst::beta1,2)*LL2*x*as5;
	  f2k += 3.*resconst::B1q*resconst::beta0*resconst::beta2*LL2*x*as5;
	  f2k += resconst::A4q*LL*x*as4;
	  f2k += 4.*resconst::A4q*resconst::beta0*LL2*x*as5;
	  f2k += 3.*resconst::A3q*resconst::beta1*LL2*x*as5;
	  f2k += resconst::A1q*resconst::beta3*LL2*x*as5;
	}
      
      //N4LL
      if (opts.order_sudak >= 4)
	{
	  f2k += resconst::B4q*x*as4;
	  f2k += 4.*resconst::B4q*resconst::beta0*LL*x*as5;
	  f2k += 3.*resconst::B3q*resconst::beta1*LL*x*as5;
	  f2k += - 2.*resconst::B2q*resconst::beta2*LL*x*as5;
	  f2k += resconst::B1q*resconst::beta3*LL*x*as5;
	  f2k += resconst::A5q*LL*x*as5;
	}

      ff *= exp(-f2k);

      //double g2k = opts.g0 * pow(1. - exp(- pow(f2k/opts.g0,p)),1./p);
      //double g2k = opts.g0 * pow(exp(pow(f2k/opts.g0,p))-1.,1./p);
      //ff *= exp(-g2k);
      
      //b^4
      //LL
      f4k += - 1./2.*resconst::A1q*as*x2;
      f4k += - 1./2.*resconst::A1q*LL*as*x2;
      f4k += - resconst::A1q*resconst::beta0*LL*as2*x2;
      f4k += - 1./2.*resconst::A1q*resconst::beta0*LL2*as2*x2;
      f4k += - 3./2.*resconst::A1q*pow(resconst::beta0,2)*LL2*as3*x2;
      f4k += - 1./2.*resconst::A1q*pow(resconst::beta0,2)*LL3*as3*x2;
      f4k += - 2.*resconst::A1q*pow(resconst::beta0,3)*LL3*as4*x2;
      f4k += - 1./2.*resconst::A1q*pow(resconst::beta0,3)*LL4*as4*x2;
      f4k += - 5./2.*resconst::A1q*pow(resconst::beta0,4)*LL4*as5*x2;
      f4k += - 1./2.*resconst::A1q*pow(resconst::beta0,4)*LL5*as5*x2;

      //NLL
      if (opts.order_sudak >= 1)
	{
	  f4k += - 1./2.*resconst::B1q*as*x2;
	  f4k += - 1./2.*resconst::B1q*resconst::beta0*as2*x2;
	  f4k += - 1./2.*resconst::B1q*resconst::beta0*LL*as2*x2;
	  f4k += - resconst::B1q*pow(resconst::beta0,2)*LL*as3*x2;
	  f4k += - 1./2.*resconst::B1q*pow(resconst::beta0,2)*LL2*as3*x2;
	  f4k += - 3./2.*resconst::B1q*pow(resconst::beta0,3)*LL2*as4*x2;
	  f4k += - 1./2.*resconst::B1q*pow(resconst::beta0,3)*LL3*as4*x2;
	  f4k += - 2.*resconst::B1q*pow(resconst::beta0,4)*LL3*as5*x2;
	  f4k += - 1./2.*resconst::B1q*pow(resconst::beta0,4)*LL4*as5*x2;
	  f4k += - 1./2.*resconst::A2q*as2*x2;
	  f4k += - 1./2.*resconst::A2q*LL*as2*x2;
	  f4k += - 2.*resconst::A2q*resconst::beta0*LL*as3*x2;
	  f4k += - resconst::A2q*resconst::beta0*LL2*as3*x2;
	  f4k += - 9./2.*resconst::A2q*pow(resconst::beta0,2)*LL2*as4*x2;
	  f4k += - 3./2.*resconst::A2q*pow(resconst::beta0,2)*LL3*as4*x2;
	  f4k += - 8.*resconst::A2q*pow(resconst::beta0,3)*LL3*as5*x2;
	  f4k += - 2.*resconst::A2q*pow(resconst::beta0,3)*LL4*as5*x2;
	  f4k += - resconst::A1q*resconst::beta1*LL*as3*x2;
	  f4k += - 1./2.*resconst::A1q*resconst::beta1*LL2*as3*x2;
	  f4k += - 15./4.*resconst::A1q*resconst::beta0*resconst::beta1*LL2*as4*x2;
	  f4k += - 5./4.*resconst::A1q*resconst::beta0*resconst::beta1*LL3*as4*x2;
	  f4k += - 26./3.*resconst::A1q*pow(resconst::beta0,2)*resconst::beta1*LL3*as5*x2;
	  f4k += - 13./6.*resconst::A1q*pow(resconst::beta0,2)*resconst::beta1*LL4*as5*x2;
	 }

      //NNLL	 
      if (opts.order_sudak >= 2)
	{
	  f4k += - 1./2.*resconst::B2q*as2*x2;
	  f4k += - resconst::B2q*resconst::beta0*as3*x2;
	  f4k += - resconst::B2q*resconst::beta0*LL*as3*x2;
	  f4k += - 3./2.*resconst::B2q*pow(resconst::beta0,2)*as4*x2*LL2;
	  f4k += - 3*resconst::B2q*pow(resconst::beta0,2)*LL*as4*x2;
	  f4k += - 2*resconst::B2q*pow(resconst::beta0,3)*as5*x2*LL3;
	  f4k += - 6*resconst::B2q*pow(resconst::beta0,3)*as5*x2*LL2;
	  f4k += - 1./2.*resconst::B1q*resconst::beta1*as3*x2;
	  f4k += - 1./2.*resconst::B1q*resconst::beta1*LL*as3*x2;
	  f4k += - 5./4.*resconst::B1q*resconst::beta0*resconst::beta1*as4*x2*LL2;
	  f4k += - 5./2.*resconst::B1q*resconst::beta0*resconst::beta1*LL*as4*x2;
	  f4k += - 13./6.*resconst::B1q*pow(resconst::beta0,2)*resconst::beta1*as5*x2*LL3;
	  f4k += - 13./2.*resconst::B1q*pow(resconst::beta0,2)*resconst::beta1*as5*x2*LL2;
	  f4k += - 1./2.*resconst::A3q*as3*x2;
	  f4k += - 1./2.*resconst::A3q*LL*as3*x2;
	  f4k += - 3./2.*resconst::A3q*resconst::beta0*as4*x2*LL2;
	  f4k += - 3*resconst::A3q*resconst::beta0*LL*as4*x2;
	  f4k += - 3*resconst::A3q*pow(resconst::beta0,2)*as5*x2*LL3;
	  f4k += - 9*resconst::A3q*pow(resconst::beta0,2)*as5*x2*LL2;
	  f4k += - resconst::A2q*resconst::beta1*as4*x2*LL2;
	  f4k += - 2*resconst::A2q*resconst::beta1*LL*as4*x2;
	  f4k += - 7./2.*resconst::A2q*resconst::beta0*resconst::beta1*as5*x2*LL3;
	  f4k += - 21./2.*resconst::A2q*resconst::beta0*resconst::beta1*as5*x2*LL2;
	  f4k += - 1./2.*resconst::A1q*resconst::beta2*as4*x2*LL2;
	  f4k += - resconst::A1q*resconst::beta2*LL*as4*x2;
	  f4k += - 3./4.*resconst::A1q*pow(resconst::beta1,2)*as5*x2*LL3;
	  f4k += - 9./4.*resconst::A1q*pow(resconst::beta1,2)*as5*x2*LL2;
	  f4k += - 3./2.*resconst::A1q*resconst::beta0*resconst::beta2*as5*x2*LL3;
	  f4k += - 9./2.*resconst::A1q*resconst::beta0*resconst::beta2*as5*x2*LL2;
	}


      //N3LL
      if (opts.order_sudak >= 3)
	{
	  f4k += - 1./2.*resconst::B3q*as3*x2                                           ;
	  f4k += - 3./2.*resconst::B3q*resconst::beta0*as4*x2		    		  ;
	  f4k += - 3./2.*resconst::B3q*resconst::beta0*LL*as4*x2	    		  ;
	  f4k += - 3*resconst::B3q*pow(resconst::beta0,2)*as5*x2*LL2	    		  ;
	  f4k += - 6*resconst::B3q*pow(resconst::beta0,2)*LL*as5*x2	    		  ;
	  f4k += - resconst::B2q*resconst::beta1*as4*x2		    		  ;
	  f4k += - resconst::B2q*resconst::beta1*LL*as4*x2		    		  ;
	  f4k += - 7./2.*resconst::B2q*resconst::beta0*resconst::beta1*as5*x2*LL2    	  ;
	  f4k += - 7*resconst::B2q*resconst::beta0*resconst::beta1*LL*as5*x2	    	  ;
	  f4k += - 1./2.*resconst::B1q*resconst::beta2*as4*x2		    		  ;
	  f4k += - 1./2.*resconst::B1q*resconst::beta2*LL*as4*x2	    		  ;
	  f4k += - 3./4.*resconst::B1q*pow(resconst::beta1,2)*as5*x2*LL2   		  ;
	  f4k += - 3./2.*resconst::B1q*pow(resconst::beta1,2)*LL*as5*x2    		  ;
	  f4k += - 3./2.*resconst::B1q*resconst::beta0*resconst::beta2*as5*x2*LL2    	  ;
	  f4k += - 3*resconst::B1q*resconst::beta0*resconst::beta2*LL*as5*x2	    	  ;
	  f4k += - 1./2.*resconst::A4q*as4*x2		    				  ;
	  f4k += - 1./2.*resconst::A4q*LL*as4*x2		    			  ;
	  f4k += - 2*resconst::A4q*resconst::beta0*as5*x2*LL2		    		  ;
	  f4k += - 4*resconst::A4q*resconst::beta0*LL*as5*x2		    		  ;
	  f4k += - 3./2.*resconst::A3q*resconst::beta1*as5*x2*LL2	    		  ;
	  f4k += - 3*resconst::A3q*resconst::beta1*LL*as5*x2		    		  ;
	  f4k += - 1./2.*resconst::A1q*resconst::beta3*as5*x2*LL2	    		  ;
	  f4k += - resconst::A1q*resconst::beta3*LL*as5*x2                              ;
	}

      //N4LL
      if (opts.order_sudak >= 4)
	{
	  f4k += - 1./2.*B3q*as3*x2                     ;
	  f4k += - 3./2.*B3q*beta0*as4*x2		    ;
	  f4k += - 3./2.*B3q*beta0*LL*as4*x2	    ;
	  f4k += - 3*B3q*pow(beta0,2)*as5*x2*LL2	    ;
	  f4k += - 6*B3q*pow(beta0,2)*LL*as5*x2	    ;
	  f4k += - B2q*beta1*as4*x2			    ;
	  f4k += - B2q*beta1*LL*as4*x2		    ;
	  f4k += - 7./2.*B2q*beta0*beta1*as5*x2*LL2	    ;
	  f4k += - 7*B2q*beta0*beta1*LL*as5*x2	    ;
	  f4k += - 1./2.*B1q*beta2*as4*x2		    ;
	  f4k += - 1./2.*B1q*beta2*LL*as4*x2	    ;
	  f4k += - 3./4.*B1q*pow(beta1,2)*as5*x2*LL2    ;
	  f4k += - 3./2.*B1q*pow(beta1,2)*LL*as5*x2	    ;
	  f4k += - 3./2.*B1q*beta0*beta2*as5*x2*LL2	    ;
	  f4k += - 3*B1q*beta0*beta2*LL*as5*x2	    ;
	  f4k += - 1./2.*A4q*as4*x2			    ;
	  f4k += - 1./2.*A4q*LL*as4*x2		    ;
	  f4k += - 2*A4q*beta0*as5*x2*LL2		    ;
	  f4k += - 4*A4q*beta0*LL*as5*x2		    ;
	  f4k += - 3./2.*A3q*beta1*as5*x2*LL2	    ;
	  f4k += - 3*A3q*beta1*LL*as5*x2		    ;
	  f4k += - 1./2.*A1q*beta3*as5*x2*LL2	    ;
	  f4k += - A1q*beta3*LL*as5*x2                  ;
	}

      f6k += + 1./2.*A1q*as*x3;
      f6k += + 1./3.*A1q*LL*as*x3;
      f6k += + 1./3.*A1q*beta0*as2*x3;
      f6k += + 1./3.*A1q*beta0*as2*x3*LL2;
      f6k += + A1q*beta0*LL*as2*x3;
      f6k += + 1./3.*A1q*pow(beta0,2)*as3*x3*LL3;
      f6k += + 3./2.*A1q*pow(beta0,2)*as3*x3*LL2;
      f6k += + A1q*pow(beta0,2)*LL*as3*x3;
      f6k += + 1./3.*A1q*pow(beta0,3)*as4*x3*LL4;
      f6k += + 2*A1q*pow(beta0,3)*as4*x3*LL3;
      f6k += + 2*A1q*pow(beta0,3)*as4*x3*LL2;
      f6k += + 1./3.*A1q*pow(beta0,4)*as5*x3*LL5;
      f6k += + 5./2.*A1q*pow(beta0,4)*as5*x3*LL4;
      f6k += + 10./3.*A1q*pow(beta0,4)*as5*x3*LL3;
      
      //ff *= exp(-f4k);
      
      //double g4k = opts.g0 * pow(1. - exp(- pow(f4k/opts.g0,p)),1./p);
      //double g4k = opts.g0 * pow(exp(pow(f4k/opts.g0,p))-1.,1./p);
      //ff *= exp(-g4k);

      //cout << "Collins-Rogers " << " mubstar " << pdfevol::mubstar << " alphas " << pdf::alphas(real(pdfevol::mubstar),opts.cmw) << " lha alphas " << pdf::extalphas(real(pdfevol::mubstar)) << endl;


      //keep only leading terms
      //f2k = resconst::A1q*LL*x*as;
      //f4k = -1./2.*resconst::A1q*LL*as*x2;

      //double nesting of gk terms
      double g0 = opts.g0;
      double g0q = 1.0;
      double gk = g0  * (1. - exp(- f2k/g0) + 1./2.* pow(f2k/g0,2) - g0q * (1. - exp(f4k/g0q)));
      //double g0s = 1.0;
      //double gk = g0  * (1. - exp(- f2k/g0) + 1./2.* pow(f2k/g0,2) - g0q * (1. - exp(f4k/g0q)) - 1./6.* pow(f2k/g0,3) + g0s * (1. - exp(-f6k/g0s)));
      ff = exp(-gk);

      //ff = exp(-f2k);
      //ff = exp(-f2k-f4k);

      //cout << endl;
      
      //cout << "LL " << LL << " mubstar " << real(pdfevol::mubstar) << endl;
      //cout << "npff 8 " << "A1q " << resconst::A1q << " LL*as/log(m2)  " << (LL*as+beta0*LL2*as2+pow(beta0,2)*LL3*as3+pow(beta0,3)*LL4*as4+pow(beta0,4)*LL5*as5)/log(pow(m,2)) << " x     " << x << endl;
      //cout << "npff 2 " << "CF    " << CF            << " as(mub)/pi     " << pdf::alphas (real(pdfevol::mubstar),opts.cmw)/M_PI << " b/blim " << pow(real(b),2)/pow(blim::sudakov,2) << endl;
      //cout << " f2k " << f2k << endl;
      //cout << " gk " << (CF*pdf::alphas(real(pdfevol::mubstar),opts.cmw)/M_PI*pow(real(b),2))/(pow(blim::sudakov,2))*log(pow(m/opts.Q0,2)) << endl;

      //Gaussian term
      if (opts.g1 != 0)
        ff *= exp(-opts.g1/sqrt(1.+max(0.,opts.lambda)*pow(b,2)) * pow(abs(b),2));

      //Quartic term
      double lambdaq = 1;
      if (opts.q != 0)
      	ff *= exp(- (opts.q/abs(opts.q))* pow(lambdaq,4)  * (1. - exp(- abs(opts.q)*pow(real(b),4)/pow(lambdaq,4))));

      //Sextic term
      if (opts.sx > 0)
	ff *= exp(-opts.sx * pow(abs(b),6));
    }
  

  if (opts.npff == 9)
    {
      double Linf = opts.Linf;
      double L2   = opts.L2;
      double L4   = opts.L4;
      double linf = opts.linf;
      double l2   = opts.l2;
      double l4   = opts.l4;

      double bb = real(b);
      double bb2 = pow(bb,2);
      double bb3 = pow(bb,3);
      double bb4 = pow(bb,4);
      
      double gj = Linf*bb*2.*tanh(L2/Linf*bb +L4/Linf*bb3);
      double gk = linf      *tanh(l2/linf*bb2+l4/linf*bb4);

      ff *= exp(-gj -gk*log(pow(m/opts.Q0,2)));
    }


  //Calculate perturbative correction to gk from S(b)-S(b*), exactly without any expansion in b^2/blim^2
  if (opts.npff == 10)
    {
      double as = resint::aass;
      double Q = scales::res;
      double LL = real(btilde::blog(b,Q)); //do not use bstar here, but b
      double LLst = real(btilde::blog(bstar::calc(b,blim::sudakov),Q));
      double y = beta0*resint::aass*LL;
      double yst = beta0*resint::aass*LLst;

      double fk = 0.;

      //LL
      fk +=
	+ A1q/beta0*LLst
	- A1q/beta0*LL
	+ log(1 - yst)*A1q/beta0*pow(yst,-1)*LLst
	- log(1 - y)*A1q/beta0*pow(y,-1)*LL;

      //NLL
      if (opts.order_sudak >= 1)
	{
	  fk +=
	    - A2q*pow(beta0,-2)*yst*pow((1-yst),-1)
	    + A2q*pow(beta0,-2)*y*pow((1-y),-1)
	    + A1q*pow(beta0,-3)*beta1*yst*pow((1-yst),-1)
	    - A1q*pow(beta0,-3)*beta1*y*pow((1-y),-1)
	    + log(1 - yst)*B1q*pow(beta0,-1)
	    - log(1 - yst)*A2q*pow(beta0,-2)
	    + log(1 - yst)*A1q*pow(beta0,-3)*beta1*pow((1-yst),-1)
	    + 1./2.*pow(log(1 - yst),2)*A1q*pow(beta0,-3)*beta1
	    - log(1 - y)*B1q*pow(beta0,-1)
	    + log(1 - y)*A2q*pow(beta0,-2)
	    - log(1 - y)*A1q*pow(beta0,-3)*beta1*pow((1-y),-1)
	    - 1./2.*pow(log(1 - y),2)*A1q*pow(beta0,-3)*beta1;
	}	

      //NNLL
      if (opts.order_sudak >= 2)
	{
	  fk +=
	    - B2q*pow(beta0,-1)*yst*as*pow((1-yst),-1)
	    + B2q*pow(beta0,-1)*y*as*pow((1-y),-1)
	    + B1q*pow(beta0,-2)*beta1*yst*as*pow((1-yst),-1)
	    - B1q*pow(beta0,-2)*beta1*y*as*pow((1-y),-1)
	    - 1./2.*A3q*pow(beta0,-2)*pow(yst,2)*as*pow((1-yst),-2)
	    + 1./2.*A3q*pow(beta0,-2)*pow(y,2)*as*pow((1-y),-2)
	    - A2q*pow(beta0,-3)*beta1*yst*as*pow((1-yst),-2)
	    + 3./2.*A2q*pow(beta0,-3)*beta1*pow(yst,2)*as*pow((1-yst),-2)
	    + A2q*pow(beta0,-3)*beta1*y*as*pow((1-y),-2)
	    - 3./2.*A2q*pow(beta0,-3)*beta1*pow(y,2)*as*pow((1-y),-2)
	    + 1./2.*A1q*pow(beta0,-4)*pow(beta1,2)*pow(yst,2)*as*pow(
								     (1-yst),-2)
	    - 1./2.*A1q*pow(beta0,-4)*pow(beta1,2)*pow(y,2)*as*pow((1-y),-2)
	    + A1q*pow(beta0,-3)*beta2*yst*as*pow((1-yst),-2)
	    - 3./2.*A1q*pow(beta0,-3)*beta2*pow(yst,2)*as*pow(
							      (1-yst),-2)
	    - A1q*pow(beta0,-3)*beta2*y*as*pow((1-y),-2)
	    + 3./2.*A1q*pow(beta0,-3)*beta2*pow(y,2)*as*pow((1-y),-2)
	    + log(1 - yst)*B1q*pow(beta0,-2)*beta1*as*pow((1-yst),-1)
	    - log(1 - yst)*A2q*pow(beta0,-3)*beta1*as*pow((1-yst),-2)
	    + 2*log(1 - yst)*A2q*pow(beta0,-3)*beta1*yst*as*pow((1-yst),-2)
	    + log(1 - yst)*A1q*pow(beta0,-4)*pow(beta1,2)*yst*as*pow(
								     (1-yst),-2)
	    - log(1 - yst)*A1q*pow(beta0,-4)*pow(beta1,2)*pow(yst,2)*as*pow(
									    (1-yst),-2)
	    + log(1 - yst)*A1q*pow(beta0,-3)*beta2*as*pow((1-yst),-2)
	    - 2*log(1 - yst)*A1q*pow(beta0,-3)*beta2*yst*as*pow((1-yst),-2)
	    + log(1 - yst)*A1q*pow(beta0,-3)*beta2*pow(yst,2)*as*pow(
								     (1-yst),-2)
	    + 1./2.*pow(log(1 - yst),2)*A1q*pow(beta0,-4)*pow(beta1,2)
	    *as*pow((1-yst),-2)
	    - pow(log(1 - yst),2)*A1q*pow(beta0,-4)*pow(beta1,2)*yst*as*pow(
									    (1-yst),-2)
	    - log(1 - y)*B1q*pow(beta0,-2)*beta1*as*pow((1-y),-1)
	    + log(1 - y)*A2q*pow(beta0,-3)*beta1*as*pow((1-y),-2)
	    - 2*log(1 - y)*A2q*pow(beta0,-3)*beta1*y*as*pow((1-y),-2)
	    - log(1 - y)*A1q*pow(beta0,-4)*pow(beta1,2)*y*as*pow((1-y),-2)
	    + log(1 - y)*A1q*pow(beta0,-4)*pow(beta1,2)*pow(y,2)*as*pow(
									(1-y),-2)
	    - log(1 - y)*A1q*pow(beta0,-3)*beta2*as*pow((1-y),-2)
	    + 2*log(1 - y)*A1q*pow(beta0,-3)*beta2*y*as*pow((1-y),-2)
	    - log(1 - y)*A1q*pow(beta0,-3)*beta2*pow(y,2)*as*pow((1-y),-2)
	    - 1./2.*pow(log(1 - y),2)*A1q*pow(beta0,-4)*pow(beta1,2)*as*pow(
									    (1-y),-2)
	    + pow(log(1 - y),2)*A1q*pow(beta0,-4)*pow(beta1,2)*y*as*
	    pow((1-y),-2);
	}


      //NNNLL
      if (opts.order_sudak >= 3)
	{
	  fk +=
	    - B3q*pow(beta0,-1)*yst*pow(as,2)*pow((1-yst),-2)
	    + 1./2.*B3q*pow(beta0,-1)*pow(yst,2)*pow(as,2)*pow((1-yst),-2)
	    + B3q*pow(beta0,-1)*y*pow(as,2)*pow((1-y),-2)
	    - 1./2.*B3q*pow(beta0,-1)*pow(y,2)*pow(as,2)*pow((1-y),-2)
	    + B2q*pow(beta0,-2)*beta1*yst*pow(as,2)*pow((1-yst),-2)
	    - 1./2.*B2q*pow(beta0,-2)*beta1*pow(yst,2)*pow(as,2)*pow(
								     (1-yst),-2)
	    - B2q*pow(beta0,-2)*beta1*y*pow(as,2)*pow((1-y),-2)
	    + 1./2.*B2q*pow(beta0,-2)*beta1*pow(y,2)*pow(as,2)*pow((1-y),-2)
	    + 1./2.*B1q*pow(beta0,-3)*pow(beta1,2)*pow(yst,2)*pow(as,2)*pow(
									    (1-yst),-2)
	    - 1./2.*B1q*pow(beta0,-3)*pow(beta1,2)*pow(y,2)*pow(as,2)*pow(
									  (1-y),-2)
	    - 1./2.*B1q*pow(beta0,-2)*beta2*pow(yst,2)*pow(as,2)*pow(
								     (1-yst),-2)
	    + 1./2.*B1q*pow(beta0,-2)*beta2*pow(y,2)*pow(as,2)*pow(
								   (1-y),-2)
	    - 1./2.*A4q*pow(beta0,-2)*pow(yst,2)*pow(as,2)*pow((1-yst),-3)
	    + 1./6.*A4q*pow(beta0,-2)*pow(yst,3)*pow(as,2)*pow((1-yst),-3)
	    + 1./2.*A4q*pow(beta0,-2)*pow(y,2)*pow(as,2)*pow((1-y),-3)
	    - 1./6.*A4q*pow(beta0,-2)*pow(y,3)*pow(as,2)*pow((1-y),-3)
	    - 1./2.*A3q*pow(beta0,-3)*beta1*yst*pow(as,2)*pow((1-yst),-3)
	    + 5./4.*A3q*pow(beta0,-3)*beta1*pow(yst,2)*pow(as,2)*pow(
								     (1-yst),-3)
	    - 5./12.*A3q*pow(beta0,-3)*beta1*pow(yst,3)*pow(as,2)*pow(
								      (1-yst),-3)
	    + 1./2.*A3q*pow(beta0,-3)*beta1*y*pow(as,2)*pow((1-y),-3)
	    - 5./4.*A3q*pow(beta0,-3)*beta1*pow(y,2)*pow(as,2)*pow((1-y),-3)
	    + 5./12.*A3q*pow(beta0,-3)*beta1*pow(y,3)*pow(as,2)*pow(
								    (1-y),-3)
	    + 5./18.*A2q*pow(beta0,-4)*pow(beta1,2)*yst*pow(as,2)*pow(
								      (1-yst),-3)
	    - 19./36.*A2q*pow(beta0,-4)*pow(beta1,2)*pow(yst,2)*pow(as,2)*
	    pow((1-yst),-3)
	    + 19./108.*A2q*pow(beta0,-4)*pow(beta1,2)*pow(yst,3)*pow(as,2)*
	    pow((1-yst),-3)
	    - 5./18.*A2q*pow(beta0,-4)*pow(beta1,2)*y*pow(as,2)*pow(
								    (1-y),-3)
	    + 19./36.*A2q*pow(beta0,-4)*pow(beta1,2)*pow(y,2)*pow(as,2)*pow(
									    (1-y),-3)
	    - 19./108.*A2q*pow(beta0,-4)*pow(beta1,2)*pow(y,3)*pow(as,2)*
	    pow((1-y),-3)
	    + 5./18.*A1q*pow(beta0,-5)*pow(beta1,3)*yst*pow(as,2)*pow(
								      (1-yst),-3)
	    - 19./36.*A1q*pow(beta0,-5)*pow(beta1,3)*pow(yst,2)*pow(as,2)*
	    pow((1-yst),-3)
	    + 1./108.*A1q*pow(beta0,-5)*pow(beta1,3)*pow(yst,3)*pow(
								    as,2)*pow((1-yst),-3)
	    - 5./18.*A1q*pow(beta0,-5)*pow(beta1,3)*y*pow(as,2)*pow(
								    (1-y),-3)
	    + 19./36.*A1q*pow(beta0,-5)*pow(beta1,3)*pow(y,2)*pow(as,2)*pow(
									    (1-y),-3)
	    - 1./108.*A1q*pow(beta0,-5)*pow(beta1,3)*pow(y,3)*pow(as,2)*pow(
									    (1-y),-3)
	    + 13./18.*A1q*pow(beta0,-4)*beta1*beta2*yst*pow(as,2)*pow(
								      (1-yst),-3)
	    - 53./36.*A1q*pow(beta0,-4)*beta1*beta2*pow(yst,2)*pow(as,2)*
	    pow((1-yst),-3)
	    + 71./108.*A1q*pow(beta0,-4)*beta1*beta2*pow(yst,3)*pow(as,2)*
	    pow((1-yst),-3)
	    - 13./18.*A1q*pow(beta0,-4)*beta1*beta2*y*pow(as,2)*pow(
								    (1-y),-3)
	    + 53./36.*A1q*pow(beta0,-4)*beta1*beta2*pow(y,2)*pow(as,2)
	    *pow((1-y),-3)
	    - 71./108.*A1q*pow(beta0,-4)*beta1*beta2*pow(y,3)*pow(as,2)*pow(
									    (1-y),-3)
	    + 1./2.*A1q*pow(beta0,-3)*beta3*yst*pow(as,2)*pow((1-yst),-3)
	    - 5./4.*A1q*pow(beta0,-3)*beta3*pow(yst,2)*pow(as,2)*pow(
								     (1-yst),-3)
	    + 7./12.*A1q*pow(beta0,-3)*beta3*pow(yst,3)*pow(as,2)*pow(
								      (1-yst),-3)
	    - 1./2.*A1q*pow(beta0,-3)*beta3*y*pow(as,2)*pow((1-y),-3)
	    + 5./4.*A1q*pow(beta0,-3)*beta3*pow(y,2)*pow(as,2)*pow((1-y),-3)
	    - 7./12.*A1q*pow(beta0,-3)*beta3*pow(y,3)*pow(as,2)*pow(
								    (1-y),-3)
	    + log(1 - yst)*B2q*pow(beta0,-2)*beta1*pow(as,2)*pow((1-yst),-2)
	    - 1./2.*log(1 - yst)*A3q*pow(beta0,-3)*beta1*pow(as,2)*pow(
								       (1-yst),-3)
	    + 3./2.*log(1 - yst)*A3q*pow(beta0,-3)*beta1*yst*pow(as,2)
	    *pow((1-yst),-3)
	    + 5./18.*log(1 - yst)*A2q*pow(beta0,-4)*pow(beta1,2)*pow(as,2)*
	    pow((1-yst),-3)
	    - 1./2.*log(1 - yst)*A2q*pow(beta0,-4)*pow(beta1,2)*yst*pow(
									as,2)*pow((1-yst),-3)
	    + 5./18.*log(1 - yst)*A1q*pow(beta0,-5)*pow(beta1,3)*pow(as,2)*
	    pow((1-yst),-3)
	    - 1./2.*log(1 - yst)*A1q*pow(beta0,-5)*pow(beta1,3)*yst*pow(
									as,2)*pow((1-yst),-3)
	    - 1./2.*log(1 - yst)*A1q*pow(beta0,-5)*pow(beta1,3)*pow(yst,2)*
	    pow(as,2)*pow((1-yst),-3)
	    - 1./2.*log(1 - yst)*A1q*pow(beta0,-5)*pow(beta1,3)*pow(yst,3)*
	    pow(as,2)*pow((1-yst),-3)
	    + 13./18.*log(1 - yst)*A1q*pow(beta0,-4)*beta1*beta2*pow(as,2)*
	    pow((1-yst),-3)
	    + 1./2.*log(1 - yst)*A1q*pow(beta0,-4)*beta1*beta2*yst*
	    pow(as,2)*pow((1-yst),-3)
	    - 2*log(1 - yst)*A1q*pow(beta0,-4)*beta1*beta2*yst*pow(as,2)*
	    pow((1-yst),-2)
	    - 3./2.*log(1 - yst)*A1q*pow(beta0,-4)*beta1*beta2*pow(yst,2)*
	    pow(as,2)*pow((1-yst),-3)
	    + 1./2.*log(1 - yst)*A1q*pow(beta0,-4)*beta1*beta2*pow(yst,3)*
	    pow(as,2)*pow((1-yst),-3)
	    + 1./2.*log(1 - yst)*A1q*pow(beta0,-3)*beta3*pow(as,2)
	    - 1./2.*pow(log(1 - yst),2)*B1q*pow(beta0,-3)*pow(beta1,2)*pow(
									   as,2)*pow((1-yst),-2)
	    + 1./6.*pow(log(1 - yst),2)*A2q*pow(beta0,-4)*pow(beta1,2)*pow(
									   as,2)*pow((1-yst),-3)
	    - 1./2.*pow(log(1 - yst),2)*A2q*pow(beta0,-4)*pow(beta1,2)*yst*
	    pow(as,2)*pow((1-yst),-3)
	    + 1./6.*pow(log(1 - yst),2)*A1q*pow(beta0,-5)*pow(beta1,3)
	    *pow(as,2)*pow((1-yst),-3)
	    - pow(log(1 - yst),2)*A1q*pow(beta0,-5)*pow(beta1,3)*yst*pow(
									 as,2)*pow((1-yst),-3)
	    + 1./3.*pow(log(1 - yst),2)*A1q*pow(beta0,-4)*beta1*beta2*pow(
									  as,2)*pow((1-yst),-3)
	    - pow(log(1 - yst),2)*A1q*pow(beta0,-4)*beta1*beta2*yst*pow(
									as,2)*pow((1-yst),-3)
	    - 1./6.*pow(log(1 - yst),3)*A1q*pow(beta0,-5)*pow(beta1,3)*pow(
									   as,2)*pow((1-yst),-3)
	    + 1./2.*pow(log(1 - yst),3)*A1q*pow(beta0,-5)*pow(beta1,3)*yst*
	    pow(as,2)*pow((1-yst),-3)
	    - log(1 - y)*B2q*pow(beta0,-2)*beta1*pow(as,2)*pow((1-y),-2)
	    + 1./2.*log(1 - y)*A3q*pow(beta0,-3)*beta1*pow(as,2)*pow(
								     (1-y),-3)
	    - 3./2.*log(1 - y)*A3q*pow(beta0,-3)*beta1*y*pow(as,2)*
	    pow((1-y),-3)
	    - 5./18.*log(1 - y)*A2q*pow(beta0,-4)*pow(beta1,2)*pow(as,2)*
	    pow((1-y),-3)
	    + 1./2.*log(1 - y)*A2q*pow(beta0,-4)*pow(beta1,2)*y*pow(as,2)*
	    pow((1-y),-3)
	    - 5./18.*log(1 - y)*A1q*pow(beta0,-5)*pow(beta1,3)*pow(as,2)*
	    pow((1-y),-3)
	    + 1./2.*log(1 - y)*A1q*pow(beta0,-5)*pow(beta1,3)*y*pow(as,2)*
	    pow((1-y),-3)
	    + 1./2.*log(1 - y)*A1q*pow(beta0,-5)*pow(beta1,3)*pow(y,2)*pow(
									   as,2)*pow((1-y),-3)
	    + 1./2.*log(1 - y)*A1q*pow(beta0,-5)*pow(beta1,3)*pow(y,3)*pow(
									   as,2)*pow((1-y),-3)
	    - 13./18.*log(1 - y)*A1q*pow(beta0,-4)*beta1*beta2*pow(as,2)*
	    pow((1-y),-3)
	    - 1./2.*log(1 - y)*A1q*pow(beta0,-4)*beta1*beta2*y*pow(
								   as,2)*pow((1-y),-3)
	    + 2*log(1 - y)*A1q*pow(beta0,-4)*beta1*beta2*y*pow(as,2)*pow(
									 (1-y),-2)
	    + 3./2.*log(1 - y)*A1q*pow(beta0,-4)*beta1*beta2*pow(y,2)*pow(
									  as,2)*pow((1-y),-3)
	    - 1./2.*log(1 - y)*A1q*pow(beta0,-4)*beta1*beta2*pow(y,3)*pow(
									  as,2)*pow((1-y),-3)
	    - 1./2.*log(1 - y)*A1q*pow(beta0,-3)*beta3*pow(as,2)
	    + 1./2.*pow(log(1 - y),2)*B1q*pow(beta0,-3)*pow(beta1,2)*pow(
									 as,2)*pow((1-y),-2)
	    - 1./6.*pow(log(1 - y),2)*A2q*pow(beta0,-4)*pow(beta1,2)*pow(
									 as,2)*pow((1-y),-3)
	    + 1./2.*pow(log(1 - y),2)*A2q*pow(beta0,-4)*pow(beta1,2)*y*pow(
									   as,2)*pow((1-y),-3)
	    - 1./6.*pow(log(1 - y),2)*A1q*pow(beta0,-5)*pow(beta1,3)*
	    pow(as,2)*pow((1-y),-3)
	    + pow(log(1 - y),2)*A1q*pow(beta0,-5)*pow(beta1,3)*y*pow(as,2)*
	    pow((1-y),-3)
	    - 1./3.*pow(log(1 - y),2)*A1q*pow(beta0,-4)*beta1*beta2*pow(
									as,2)*pow((1-y),-3)
	    + pow(log(1 - y),2)*A1q*pow(beta0,-4)*beta1*beta2*y*pow(as,2)*
	    pow((1-y),-3)
	    + 1./6.*pow(log(1 - y),3)*A1q*pow(beta0,-5)*pow(beta1,3)*pow(
									 as,2)*pow((1-y),-3)
	    - 1./2.*pow(log(1 - y),3)*A1q*pow(beta0,-5)*pow(beta1,3)*y*pow(
									   as,2)*pow((1-y),-3)
	    ;
	}

//      //NNNNLL
//      if (opts.order_sudak >= 4)
//	{
//	  fk +=
//	    - B4q*pow(beta0,-1)*yst*pow(as,3)*pow((1-yst),-3)
//	    + B4q*pow(beta0,-1)*pow(yst,2)*pow(as,3)*pow((1-yst),-3)
//	    - 1./3.*B4q*pow(beta0,-1)*pow(yst,3)*pow(as,3)*pow((1-yst),-3)
//	    + B4q*pow(beta0,-1)*y*pow(as,3)*pow((1-y),-3)
//	    - B4q*pow(beta0,-1)*pow(y,2)*pow(as,3)*pow((1-y),-3)
//	    + 1./3.*B4q*pow(beta0,-1)*pow(y,3)*pow(as,3)*pow((1-y),-3)
//	    + B3q*pow(beta0,-2)*beta1*yst*pow(as,3)*pow((1-yst),-3)
//	    - B3q*pow(beta0,-2)*beta1*pow(yst,2)*pow(as,3)*pow((1-yst),-3)
//	    + 1./3.*B3q*pow(beta0,-2)*beta1*pow(yst,3)*pow(as,3)*pow(
//								     (1-yst),-3)
//	    - B3q*pow(beta0,-2)*beta1*y*pow(as,3)*pow((1-y),-3)
//	    + B3q*pow(beta0,-2)*beta1*pow(y,2)*pow(as,3)*pow((1-y),-3)
//	    - 1./3.*B3q*pow(beta0,-2)*beta1*pow(y,3)*pow(as,3)*pow((1-y),-3)
//	    - B2q*pow(beta0,-3)*pow(beta1,2)*pow(yst,2)*pow(as,3)*pow(
//								      (1-yst),-3)
//	    + 1./3.*B2q*pow(beta0,-3)*pow(beta1,2)*pow(yst,3)*pow(
//								  as,3)*pow((1-yst),-3)
//	    + B2q*pow(beta0,-3)*pow(beta1,2)*pow(y,2)*pow(as,3)*pow(
//								    (1-y),-3)
//	    - 1./3.*B2q*pow(beta0,-3)*pow(beta1,2)*pow(y,3)*pow(as,3)*pow(
//									  (1-y),-3)
//	    + B2q*pow(beta0,-2)*beta2*pow(yst,2)*pow(as,3)*pow((1-yst),-3)
//	    - 1./3.*B2q*pow(beta0,-2)*beta2*pow(yst,3)*pow(as,3)*pow(
//								     (1-yst),-3)
//	    - B2q*pow(beta0,-2)*beta2*pow(y,2)*pow(as,3)*pow((1-y),-3)
//	    + 1./3.*B2q*pow(beta0,-2)*beta2*pow(y,3)*pow(as,3)*pow((1-y),-3)
//	    - 1./2.*B1q*pow(beta0,-4)*pow(beta1,3)*pow(yst,2)*pow(as,3)*pow(
//									    (1-yst),-3)
//	    + 1./3.*B1q*pow(beta0,-4)*pow(beta1,3)*pow(yst,3)*pow(as,3)*pow(
//									    (1-yst),-3)
//	    + 1./2.*B1q*pow(beta0,-4)*pow(beta1,3)*pow(y,2)*pow(as,3)*
//	    pow((1-y),-3)
//	    - 1./3.*B1q*pow(beta0,-4)*pow(beta1,3)*pow(y,3)*pow(as,3)*pow(
//									  (1-y),-3)
//	    + B1q*pow(beta0,-3)*beta1*beta2*pow(yst,2)*pow(as,3)*pow(
//								     (1-yst),-3)
//	    - 2./3.*B1q*pow(beta0,-3)*beta1*beta2*pow(yst,3)*pow(as,3)*pow(
//									   (1-yst),-3)
//	    - B1q*pow(beta0,-3)*beta1*beta2*pow(y,2)*pow(as,3)*pow((1-y),-3)
//	    + 2./3.*B1q*pow(beta0,-3)*beta1*beta2*pow(y,3)*pow(as,3)*pow(
//									 (1-y),-3)
//	    - 1./2.*B1q*pow(beta0,-2)*beta3*pow(yst,2)*pow(as,3)*pow(
//								     (1-yst),-3)
//	    + 1./3.*B1q*pow(beta0,-2)*beta3*pow(yst,3)*pow(as,3)*pow(
//								     (1-yst),-3)
//	    + 1./2.*B1q*pow(beta0,-2)*beta3*pow(y,2)*pow(as,3)*pow(
//								   (1-y),-3)
//	    - 1./3.*B1q*pow(beta0,-2)*beta3*pow(y,3)*pow(as,3)*pow((1-y),-3)
//	    - 1./2.*A5q*pow(beta0,-2)*pow(yst,2)*pow(as,3)*pow((1-yst),-4)
//	    + 1./3.*A5q*pow(beta0,-2)*pow(yst,3)*pow(as,3)*pow((1-yst),-4)
//	    - 1./12.*A5q*pow(beta0,-2)*pow(yst,4)*pow(as,3)*pow((1-yst),-4)
//	    + 1./2.*A5q*pow(beta0,-2)*pow(y,2)*pow(as,3)*pow((1-y),-4)
//	    - 1./3.*A5q*pow(beta0,-2)*pow(y,3)*pow(as,3)*pow((1-y),-4)
//	    + 1./12.*A5q*pow(beta0,-2)*pow(y,4)*pow(as,3)*pow((1-y),-4)
//	    - 1./3.*A4q*pow(beta0,-3)*beta1*yst*pow(as,3)*pow((1-yst),-4)
//	    + 7./6.*A4q*pow(beta0,-3)*beta1*pow(yst,2)*pow(as,3)*pow(
//								     (1-yst),-4)
//	    - 7./9.*A4q*pow(beta0,-3)*beta1*pow(yst,3)*pow(as,3)*pow(
//								     (1-yst),-4)
//	    + 7./36.*A4q*pow(beta0,-3)*beta1*pow(yst,4)*pow(as,3)*pow(
//								      (1-yst),-4)
//	    + 1./3.*A4q*pow(beta0,-3)*beta1*y*pow(as,3)*pow((1-y),-4)
//	    - 7./6.*A4q*pow(beta0,-3)*beta1*pow(y,2)*pow(as,3)*pow((1-y),-4)
//	    + 7./9.*A4q*pow(beta0,-3)*beta1*pow(y,3)*pow(as,3)*pow((1-y),-4)
//	    - 7./36.*A4q*pow(beta0,-3)*beta1*pow(y,4)*pow(as,3)*pow(
//								    (1-y),-4)
//	    + 1./3.*A3q*pow(beta0,-4)*pow(beta1,2)*yst*pow(as,3)*pow(
//								     (1-yst),-4)
//	    - 2./3.*A3q*pow(beta0,-4)*pow(beta1,2)*pow(yst,2)*pow(as,3)*pow(
//									    (1-yst),-4)
//	    + 13./9.*A3q*pow(beta0,-4)*pow(beta1,2)*pow(yst,3)*pow(as,3)*
//	    pow((1-yst),-4)
//	    - 13./36.*A3q*pow(beta0,-4)*pow(beta1,2)*pow(yst,4)*pow(as,3)*
//	    pow((1-yst),-4)
//	    - 1./3.*A3q*pow(beta0,-4)*pow(beta1,2)*y*pow(as,3)*pow((1-y),-4)
//	    + 2./3.*A3q*pow(beta0,-4)*pow(beta1,2)*pow(y,2)*pow(as,3)*pow(
//									  (1-y),-4)
//	    - 13./9.*A3q*pow(beta0,-4)*pow(beta1,2)*pow(y,3)*pow(as,3)
//	    *pow((1-y),-4)
//	    + 13./36.*A3q*pow(beta0,-4)*pow(beta1,2)*pow(y,4)*pow(as,3)*pow(
//									    (1-y),-4)
//	    - A3q*pow(beta0,-3)*beta2*pow(yst,3)*pow(as,3)*pow((1-yst),-4)
//	    + 1./4.*A3q*pow(beta0,-3)*beta2*pow(yst,4)*pow(as,3)*pow(
//								     (1-yst),-4)
//	    + A3q*pow(beta0,-3)*beta2*pow(y,3)*pow(as,3)*pow((1-y),-4)
//	    - 1./4.*A3q*pow(beta0,-3)*beta2*pow(y,4)*pow(as,3)*pow((1-y),-4)
//	    - 1./3.*A2q*pow(beta0,-5)*pow(beta1,3)*yst*pow(as,3)*pow(
//								     (1-yst),-4)
//	    + 7./6.*A2q*pow(beta0,-5)*pow(beta1,3)*pow(yst,2)*pow(as,3)*pow(
//									    (1-yst),-4)
//	    - 16./9.*A2q*pow(beta0,-5)*pow(beta1,3)*pow(yst,3)*pow(as,3)*
//	    pow((1-yst),-4)
//	    + 25./36.*A2q*pow(beta0,-5)*pow(beta1,3)*pow(yst,4)*pow(
//								    as,3)*pow((1-yst),-4)
//	    + 1./3.*A2q*pow(beta0,-5)*pow(beta1,3)*y*pow(as,3)*pow((1-y),-4)
//	    - 7./6.*A2q*pow(beta0,-5)*pow(beta1,3)*pow(y,2)*pow(as,3)*pow(
//									  (1-y),-4)
//	    + 16./9.*A2q*pow(beta0,-5)*pow(beta1,3)*pow(y,3)*pow(as,3)*pow(
//									   (1-y),-4)
//	    - 25./36.*A2q*pow(beta0,-5)*pow(beta1,3)*pow(y,4)*pow(as,3)*pow(
//									    (1-y),-4)
//	    + 1./3.*A2q*pow(beta0,-4)*beta1*beta2*yst*pow(as,3)*pow(
//								    (1-yst),-4)
//	    - 7./6.*A2q*pow(beta0,-4)*beta1*beta2*pow(yst,2)*pow(as,3)*pow(
//									   (1-yst),-4)
//	    + 22./9.*A2q*pow(beta0,-4)*beta1*beta2*pow(yst,3)*pow(as,3)*pow(
//									    (1-yst),-4)
//	    - 10./9.*A2q*pow(beta0,-4)*beta1*beta2*pow(yst,4)*pow(
//								  as,3)*pow((1-yst),-4)
//	    - 1./3.*A2q*pow(beta0,-4)*beta1*beta2*y*pow(as,3)*pow((1-y),-4)
//	    + 7./6.*A2q*pow(beta0,-4)*beta1*beta2*pow(y,2)*pow(as,3)*pow(
//									 (1-y),-4)
//	    - 22./9.*A2q*pow(beta0,-4)*beta1*beta2*pow(y,3)*pow(as,3)*pow(
//									  (1-y),-4)
//	    + 10./9.*A2q*pow(beta0,-4)*beta1*beta2*pow(y,4)*pow(as,3)*pow(
//									  (1-y),-4)
//	    - 2./3.*A2q*pow(beta0,-3)*beta3*pow(yst,3)*pow(as,3)*pow(
//								     (1-yst),-4)
//	    + 5./12.*A2q*pow(beta0,-3)*beta3*pow(yst,4)*pow(as,3)*pow(
//								      (1-yst),-4)
//	    + 2./3.*A2q*pow(beta0,-3)*beta3*pow(y,3)*pow(as,3)*pow((1-y),-4)
//	    - 5./12.*A2q*pow(beta0,-3)*beta3*pow(y,4)*pow(as,3)*pow(
//								    (1-y),-4)
//	    - 1./6.*A1q*pow(beta0,-6)*pow(beta1,4)*pow(yst,3)*pow(
//								  as,3)*pow((1-yst),-4)
//	    - 1./4.*A1q*pow(beta0,-6)*pow(beta1,4)*pow(yst,4)*pow(as,3)*pow(
//									    (1-yst),-4)
//	    + 1./6.*A1q*pow(beta0,-6)*pow(beta1,4)*pow(y,3)*pow(as,3)*pow(
//									  (1-y),-4)
//	    + 1./4.*A1q*pow(beta0,-6)*pow(beta1,4)*pow(y,4)*pow(as,3)*pow(
//									  (1-y),-4)
//	    + 1./3.*A1q*pow(beta0,-5)*pow(beta1,2)*beta2*yst*pow(as,3)*pow(
//									   (1-yst),-4)
//	    - 7./6.*A1q*pow(beta0,-5)*pow(beta1,2)*beta2*pow(yst,2)*pow(
//									as,3)*pow((1-yst),-4)
//	    + 13./9.*A1q*pow(beta0,-5)*pow(beta1,2)*beta2*pow(yst,3)*pow(
//									 as,3)*pow((1-yst),-4)
//	    + 5./36.*A1q*pow(beta0,-5)*pow(beta1,2)*beta2*pow(yst,4)*pow(
//									 as,3)*pow((1-yst),-4)
//	    - 1./3.*A1q*pow(beta0,-5)*pow(beta1,2)*beta2*y*pow(as,3)*
//	    pow((1-y),-4)
//	    + 7./6.*A1q*pow(beta0,-5)*pow(beta1,2)*beta2*pow(y,2)*pow(as,3)*
//	    pow((1-y),-4)
//	    - 13./9.*A1q*pow(beta0,-5)*pow(beta1,2)*beta2*pow(y,3)*pow(as,3)
//	    *pow((1-y),-4)
//	    - 5./36.*A1q*pow(beta0,-5)*pow(beta1,2)*beta2*pow(y,4)*pow(as,3)
//	    *pow((1-y),-4)
//	    - 1./3.*A1q*pow(beta0,-4)*pow(beta2,2)*yst*pow(as,3)*pow(
//								     (1-yst),-4)
//	    + 7./6.*A1q*pow(beta0,-4)*pow(beta2,2)*pow(yst,2)*pow(as,3)*pow(
//									    (1-yst),-4)
//	    - 13./9.*A1q*pow(beta0,-4)*pow(beta2,2)*pow(yst,3)*pow(as,3)*
//	    pow((1-yst),-4)
//	    + 7./36.*A1q*pow(beta0,-4)*pow(beta2,2)*pow(yst,4)*pow(as,3)*
//	    pow((1-yst),-4)
//	    + 1./3.*A1q*pow(beta0,-4)*pow(beta2,2)*y*pow(as,3)*pow(
//								   (1-y),-4)
//	    - 7./6.*A1q*pow(beta0,-4)*pow(beta2,2)*pow(y,2)*pow(as,3)*pow(
//									  (1-y),-4)
//	    + 13./9.*A1q*pow(beta0,-4)*pow(beta2,2)*pow(y,3)*pow(as,3)*pow(
//									   (1-y),-4)
//	    - 7./36.*A1q*pow(beta0,-4)*pow(beta2,2)*pow(y,4)*pow(as,3)*pow(
//									   (1-y),-4)
//	    - 1./3.*A1q*pow(beta0,-4)*beta1*beta3*yst*pow(as,3)*pow(
//								    (1-yst),-4)
//	    + 7./6.*A1q*pow(beta0,-4)*beta1*beta3*pow(yst,2)*pow(as,3)*pow(
//									   (1-yst),-4)
//	    - 17./18.*A1q*pow(beta0,-4)*beta1*beta3*pow(yst,3)*pow(as,3)*
//	    pow((1-yst),-4)
//	    + 5./18.*A1q*pow(beta0,-4)*beta1*beta3*pow(yst,4)*pow(as,3)*pow(
//									    (1-yst),-4)
//	    + 1./3.*A1q*pow(beta0,-4)*beta1*beta3*y*pow(as,3)*pow(
//								  (1-y),-4)
//	    - 7./6.*A1q*pow(beta0,-4)*beta1*beta3*pow(y,2)*pow(as,3)*pow(
//									 (1-y),-4)
//	    + 17./18.*A1q*pow(beta0,-4)*beta1*beta3*pow(y,3)*pow(as,3)*pow(
//									   (1-y),-4)
//	    - 5./18.*A1q*pow(beta0,-4)*beta1*beta3*pow(y,4)*pow(as,3)*pow(
//									  (1-y),-4)
//	    + 1./3.*A1q*pow(beta0,-3)*beta4*yst*pow(as,3)*pow((1-yst),-4)
//	    - 7./6.*A1q*pow(beta0,-3)*beta4*pow(yst,2)*pow(as,3)*pow(
//								     (1-yst),-4)
//	    + 10./9.*A1q*pow(beta0,-3)*beta4*pow(yst,3)*pow(as,3)*pow(
//								      (1-yst),-4)
//	    - 13./36.*A1q*pow(beta0,-3)*beta4*pow(yst,4)*pow(as,3)*pow(
//								       (1-yst),-4)
//	    - 1./3.*A1q*pow(beta0,-3)*beta4*y*pow(as,3)*pow((1-y),-4)
//	    + 7./6.*A1q*pow(beta0,-3)*beta4*pow(y,2)*pow(as,3)*pow((1-y),-4)
//	    - 10./9.*A1q*pow(beta0,-3)*beta4*pow(y,3)*pow(as,3)*pow(
//								    (1-y),-4)
//	    + 13./36.*A1q*pow(beta0,-3)*beta4*pow(y,4)*pow(as,3)*pow(
//								     (1-y),-4)
//	    + log(1 - yst)*B3q*pow(beta0,-2)*beta1*pow(as,3)*pow((1-yst),-3)
//	    - log(1 - yst)*B1q*pow(beta0,-4)*pow(beta1,3)*yst*pow(as,3)*pow(
//									    (1-yst),-3)
//	    + log(1 - yst)*B1q*pow(beta0,-3)*beta1*beta2*yst*pow(as,3)*pow(
//									   (1-yst),-3)
//	    - 1./3.*log(1 - yst)*A4q*pow(beta0,-3)*beta1*pow(as,3)*pow(
//								       (1-yst),-4)
//	    + 4./3.*log(1 - yst)*A4q*pow(beta0,-3)*beta1*yst*pow(as,3)*pow(
//									   (1-yst),-4)
//	    + 1./3.*log(1 - yst)*A3q*pow(beta0,-4)*pow(beta1,2)*pow(
//								    as,3)*pow((1-yst),-4)
//	    - 1./3.*log(1 - yst)*A3q*pow(beta0,-4)*pow(beta1,2)*yst*pow(
//									as,3)*pow((1-yst),-4)
//	    - 1./3.*log(1 - yst)*A2q*pow(beta0,-5)*pow(beta1,3)*pow(as,3)*
//	    pow((1-yst),-4)
//	    + 4./3.*log(1 - yst)*A2q*pow(beta0,-5)*pow(beta1,3)*yst*pow(
//									as,3)*pow((1-yst),-4)
//	    - 3*log(1 - yst)*A2q*pow(beta0,-5)*pow(beta1,3)*pow(yst,2)*pow(
//									   as,3)*pow((1-yst),-4)
//	    + 1./3.*log(1 - yst)*A2q*pow(beta0,-4)*beta1*beta2*pow(as,3)*
//	    pow((1-yst),-4)
//	    - 4./3.*log(1 - yst)*A2q*pow(beta0,-4)*beta1*beta2*yst*pow(as,3)
//	    *pow((1-yst),-4)
//	    + 3*log(1 - yst)*A2q*pow(beta0,-4)*beta1*beta2*pow(yst,2)*pow(
//									  as,3)*pow((1-yst),-4)
//	    - 1./2.*log(1 - yst)*A1q*pow(beta0,-6)*pow(beta1,4)*pow(
//								    yst,2)*pow(as,3)*pow((1-yst),-4)
//	    + 1./3.*log(1 - yst)*A1q*pow(beta0,-6)*pow(beta1,4)*pow(yst,3)*
//	    pow(as,3)*pow((1-yst),-4)
//	    - 1./3.*log(1 - yst)*A1q*pow(beta0,-6)*pow(beta1,4)*pow(yst,4)*
//	    pow(as,3)*pow((1-yst),-4)
//	    + 1./3.*log(1 - yst)*A1q*pow(beta0,-5)*pow(beta1,2)*beta2*pow(
//									  as,3)*pow((1-yst),-4)
//	    - 4./3.*log(1 - yst)*A1q*pow(beta0,-5)*pow(beta1,2)*beta2*yst*
//	    pow(as,3)*pow((1-yst),-4)
//	    + 2*log(1 - yst)*A1q*pow(beta0,-5)*pow(beta1,2)*beta2*pow(yst,2)
//	    *pow(as,3)*pow((1-yst),-4)
//	    - 2*log(1 - yst)*A1q*pow(beta0,-5)*pow(beta1,2)*beta2*pow(yst,3)
//	    *pow(as,3)*pow((1-yst),-4)
//	    + log(1 - yst)*A1q*pow(beta0,-5)*pow(beta1,2)*beta2*pow(yst,4)*
//	    pow(as,3)*pow((1-yst),-4)
//	    - 1./3.*log(1 - yst)*A1q*pow(beta0,-4)*pow(beta2,2)*pow(
//								    as,3)*pow((1-yst),-4)
//	    + 4./3.*log(1 - yst)*A1q*pow(beta0,-4)*pow(beta2,2)*yst*pow(
//									as,3)*pow((1-yst),-4)
//	    - 2*log(1 - yst)*A1q*pow(beta0,-4)*pow(beta2,2)*pow(yst,2)*pow(
//									   as,3)*pow((1-yst),-4)
//	    + 4./3.*log(1 - yst)*A1q*pow(beta0,-4)*pow(beta2,2)*pow(yst,3)*
//	    pow(as,3)*pow((1-yst),-4)
//	    - 1./3.*log(1 - yst)*A1q*pow(beta0,-4)*pow(beta2,2)*pow(yst,4)*
//	    pow(as,3)*pow((1-yst),-4)
//	    - 1./3.*log(1 - yst)*A1q*pow(beta0,-4)*beta1*beta3*pow(as,3)*
//	    pow((1-yst),-4)
//	    + 4./3.*log(1 - yst)*A1q*pow(beta0,-4)*beta1*beta3*yst*pow(as,3)
//	    *pow((1-yst),-4)
//	    - 3./2.*log(1 - yst)*A1q*pow(beta0,-4)*beta1*beta3*pow(yst,2)*
//	    pow(as,3)*pow((1-yst),-4)
//	    + 5./3.*log(1 - yst)*A1q*pow(beta0,-4)*beta1*beta3*pow(
//								   yst,3)*pow(as,3)*pow((1-yst),-4)
//	    - 2./3.*log(1 - yst)*A1q*pow(beta0,-4)*beta1*beta3*pow(yst,4)*
//	    pow(as,3)*pow((1-yst),-4)
//	    + 1./3.*log(1 - yst)*A1q*pow(beta0,-3)*beta4*pow(as,3)*pow(
//								       (1-yst),-4)
//	    - 4./3.*log(1 - yst)*A1q*pow(beta0,-3)*beta4*yst*pow(as,3)*pow(
//									   (1-yst),-4)
//	    + 2*log(1 - yst)*A1q*pow(beta0,-3)*beta4*pow(yst,2)*pow(as,3)*
//	    pow((1-yst),-4)
//	    - 4./3.*log(1 - yst)*A1q*pow(beta0,-3)*beta4*pow(yst,3)*pow(
//									as,3)*pow((1-yst),-4)
//	    + 1./3.*log(1 - yst)*A1q*pow(beta0,-3)*beta4*pow(yst,4)*pow(
//									as,3)*pow((1-yst),-4)
//	    + pow(log(1 - yst),2)*B2q*pow(beta0,-3)*pow(beta1,2)*pow(as,3)*
//	    pow((1-yst),-3)
//	    - 1./2.*pow(log(1 - yst),2)*B1q*pow(beta0,-4)*pow(beta1,3)
//	    *pow(as,3)*pow((1-yst),-3)
//	    + 1./2.*pow(log(1 - yst),2)*A3q*pow(beta0,-4)*pow(beta1,2)*pow(
//									   as,3)*pow((1-yst),-4)
//	    - 2*pow(log(1 - yst),2)*A3q*pow(beta0,-4)*pow(beta1,2)*yst*pow(
//									   as,3)*pow((1-yst),-4)
//	    - pow(log(1 - yst),2)*A2q*pow(beta0,-5)*pow(beta1,3)*yst*pow(
//									 as,3)*pow((1-yst),-4)
//	    - 1./2.*pow(log(1 - yst),2)*A1q*pow(beta0,-6)*pow(beta1,4)*yst*
//	    pow(as,3)*pow((1-yst),-4)
//	    + 3./2.*pow(log(1 - yst),2)*A1q*pow(beta0,-6)*pow(beta1,4)*pow(
//									   yst,2)*pow(as,3)*pow((1-yst),-4)
//	    - 3./2.*pow(log(1 - yst),2)*A1q*pow(beta0,-5)*pow(beta1,2)*beta2
//	    *pow(yst,2)*pow(as,3)*pow((1-yst),-4)
//	    + 1./3.*pow(log(1 - yst),3)*B1q*pow(beta0,-4)*pow(beta1,3)*pow(
//									   as,3)*pow((1-yst),-3)
//	    - 1./3.*pow(log(1 - yst),3)*A2q*pow(beta0,-5)*pow(beta1,3)
//	    *pow(as,3)*pow((1-yst),-4)
//	    + 4./3.*pow(log(1 - yst),3)*A2q*pow(beta0,-5)*pow(beta1,3)*yst*
//	    pow(as,3)*pow((1-yst),-4)
//	    - 1./6.*pow(log(1 - yst),3)*A1q*pow(beta0,-6)*pow(beta1,4)*pow(
//									   as,3)*pow((1-yst),-4)
//	    + pow(log(1 - yst),3)*A1q*pow(beta0,-6)*pow(beta1,4)*yst*pow(
//									 as,3)*pow((1-yst),-4)
//	    + 1./12.*pow(log(1 - yst),4)*A1q*pow(beta0,-6)*pow(beta1,4)*pow(
//									    as,3)*pow((1-yst),-4)
//	    - 1./3.*pow(log(1 - yst),4)*A1q*pow(beta0,-6)*pow(beta1,4)*yst*
//	    pow(as,3)*pow((1-yst),-4)
//	    - log(1 - y)*B3q*pow(beta0,-2)*beta1*pow(as,3)*pow((1-y),-3)
//	    + log(1 - y)*B1q*pow(beta0,-4)*pow(beta1,3)*y*pow(as,3)*pow(
//									(1-y),-3)
//	    - log(1 - y)*B1q*pow(beta0,-3)*beta1*beta2*y*pow(as,3)*
//	    pow((1-y),-3)
//	    + 1./3.*log(1 - y)*A4q*pow(beta0,-3)*beta1*pow(as,3)*pow(
//								     (1-y),-4)
//	    - 4./3.*log(1 - y)*A4q*pow(beta0,-3)*beta1*y*pow(as,3)*pow(
//								       (1-y),-4)
//	    - 1./3.*log(1 - y)*A3q*pow(beta0,-4)*pow(beta1,2)*pow(as,3)*pow(
//									    (1-y),-4)
//	    + 1./3.*log(1 - y)*A3q*pow(beta0,-4)*pow(beta1,2)*y*pow(as,3)*
//	    pow((1-y),-4)
//	    + 1./3.*log(1 - y)*A2q*pow(beta0,-5)*pow(beta1,3)*pow(as,3)*pow(
//									    (1-y),-4)
//	    - 4./3.*log(1 - y)*A2q*pow(beta0,-5)*pow(beta1,3)*y*pow(as,3)*
//	    pow((1-y),-4)
//	    + 3*log(1 - y)*A2q*pow(beta0,-5)*pow(beta1,3)*pow(y,2)*pow(as,3)
//	    *pow((1-y),-4)
//	    - 1./3.*log(1 - y)*A2q*pow(beta0,-4)*beta1*beta2*pow(as,3)
//	    *pow((1-y),-4)
//	    + 4./3.*log(1 - y)*A2q*pow(beta0,-4)*beta1*beta2*y*pow(as,3)*
//	    pow((1-y),-4)
//	    - 3*log(1 - y)*A2q*pow(beta0,-4)*beta1*beta2*pow(y,2)*pow(as,3)*
//	    pow((1-y),-4)
//	    + 1./2.*log(1 - y)*A1q*pow(beta0,-6)*pow(beta1,4)*pow(y,2)*pow(
//									   as,3)*pow((1-y),-4)
//	    - 1./3.*log(1 - y)*A1q*pow(beta0,-6)*pow(beta1,4)*pow(y,3)*pow(
//									   as,3)*pow((1-y),-4)
//	    + 1./3.*log(1 - y)*A1q*pow(beta0,-6)*pow(beta1,4)*pow(y,4)*pow(
//									   as,3)*pow((1-y),-4)
//	    - 1./3.*log(1 - y)*A1q*pow(beta0,-5)*pow(beta1,2)*beta2*pow(
//									as,3)*pow((1-y),-4)
//	    + 4./3.*log(1 - y)*A1q*pow(beta0,-5)*pow(beta1,2)*beta2*y*pow(
//									  as,3)*pow((1-y),-4)
//	    - 2*log(1 - y)*A1q*pow(beta0,-5)*pow(beta1,2)*beta2*pow(
//								    y,2)*pow(as,3)*pow((1-y),-4)
//	    + 2*log(1 - y)*A1q*pow(beta0,-5)*pow(beta1,2)*beta2*pow(y,3)*
//	    pow(as,3)*pow((1-y),-4)
//	    - log(1 - y)*A1q*pow(beta0,-5)*pow(beta1,2)*beta2*pow(y,4)*pow(
//									   as,3)*pow((1-y),-4)
//	    + 1./3.*log(1 - y)*A1q*pow(beta0,-4)*pow(beta2,2)*pow(as,3)*pow(
//									    (1-y),-4)
//	    - 4./3.*log(1 - y)*A1q*pow(beta0,-4)*pow(beta2,2)*y*pow(as,3)*
//	    pow((1-y),-4)
//	    + 2*log(1 - y)*A1q*pow(beta0,-4)*pow(beta2,2)*pow(y,2)*pow(as,3)
//	    *pow((1-y),-4)
//	    - 4./3.*log(1 - y)*A1q*pow(beta0,-4)*pow(beta2,2)*pow(y,3)*pow(
//									   as,3)*pow((1-y),-4)
//	    + 1./3.*log(1 - y)*A1q*pow(beta0,-4)*pow(beta2,2)*pow(y,4)*pow(
//									   as,3)*pow((1-y),-4)
//	    + 1./3.*log(1 - y)*A1q*pow(beta0,-4)*beta1*beta3*pow(as,3)
//	    *pow((1-y),-4)
//	    - 4./3.*log(1 - y)*A1q*pow(beta0,-4)*beta1*beta3*y*pow(as,3)*
//	    pow((1-y),-4)
//	    + 3./2.*log(1 - y)*A1q*pow(beta0,-4)*beta1*beta3*pow(y,2)*pow(
//									  as,3)*pow((1-y),-4)
//	    - 5./3.*log(1 - y)*A1q*pow(beta0,-4)*beta1*beta3*pow(y,3)*pow(
//									  as,3)*pow((1-y),-4)
//	    + 2./3.*log(1 - y)*A1q*pow(beta0,-4)*beta1*beta3*pow(y,4)*pow(
//									  as,3)*pow((1-y),-4)
//	    - 1./3.*log(1 - y)*A1q*pow(beta0,-3)*beta4*pow(as,3)*pow(
//								     (1-y),-4)
//	    + 4./3.*log(1 - y)*A1q*pow(beta0,-3)*beta4*y*pow(as,3)*pow(
//								       (1-y),-4)
//	    - 2*log(1 - y)*A1q*pow(beta0,-3)*beta4*pow(y,2)*pow(as,3)*pow(
//									  (1-y),-4)
//	    + 4./3.*log(1 - y)*A1q*pow(beta0,-3)*beta4*pow(y,3)*pow(
//								    as,3)*pow((1-y),-4)
//	    - 1./3.*log(1 - y)*A1q*pow(beta0,-3)*beta4*pow(y,4)*pow(as,3)*
//	    pow((1-y),-4)
//	    - pow(log(1 - y),2)*B2q*pow(beta0,-3)*pow(beta1,2)*pow(as,3)*
//	    pow((1-y),-3)
//	    + 1./2.*pow(log(1 - y),2)*B1q*pow(beta0,-4)*pow(beta1,3)*pow(
//									 as,3)*pow((1-y),-3)
//	    - 1./2.*pow(log(1 - y),2)*A3q*pow(beta0,-4)*pow(beta1,2)*pow(
//									 as,3)*pow((1-y),-4)
//	    + 2*pow(log(1 - y),2)*A3q*pow(beta0,-4)*pow(beta1,2)*y*pow(as,3)
//	    *pow((1-y),-4)
//	    + pow(log(1 - y),2)*A2q*pow(beta0,-5)*pow(beta1,3)*y*pow(as,3)*
//	    pow((1-y),-4)
//	    + 1./2.*pow(log(1 - y),2)*A1q*pow(beta0,-6)*pow(beta1,4)*y*pow(
//									   as,3)*pow((1-y),-4)
//	    - 3./2.*pow(log(1 - y),2)*A1q*pow(beta0,-6)*pow(beta1,4)*
//	    pow(y,2)*pow(as,3)*pow((1-y),-4)
//	    + 3./2.*pow(log(1 - y),2)*A1q*pow(beta0,-5)*pow(beta1,2)*beta2*
//	    pow(y,2)*pow(as,3)*pow((1-y),-4)
//	    - 1./3.*pow(log(1 - y),3)*B1q*pow(beta0,-4)*pow(beta1,3)*pow(
//									 as,3)*pow((1-y),-3)
//	    + 1./3.*pow(log(1 - y),3)*A2q*pow(beta0,-5)*pow(beta1,3)*pow(
//									 as,3)*pow((1-y),-4)
//	    - 4./3.*pow(log(1 - y),3)*A2q*pow(beta0,-5)*pow(beta1,3)*y*pow(
//									   as,3)*pow((1-y),-4)
//	    + 1./6.*pow(log(1 - y),3)*A1q*pow(beta0,-6)*pow(beta1,4)*pow(
//									 as,3)*pow((1-y),-4)
//	    - pow(log(1 - y),3)*A1q*pow(beta0,-6)*pow(beta1,4)*y*pow(as,3)*
//	    pow((1-y),-4)
//	    - 1./12.*pow(log(1 - y),4)*A1q*pow(beta0,-6)*pow(beta1,4)*pow(
//									  as,3)*pow((1-y),-4)
//	    + 1./3.*pow(log(1 - y),4)*A1q*pow(beta0,-6)*pow(beta1,4)*y
//	    *pow(as,3)*pow((1-y),-4)
//	    ;
//	}
      

      //remove log(Q^2) dependence
      fk = fk/log(pow(Q,2));

      //cout << "exact fk " << fk << endl;
      double gk = opts.g0  * (1. - exp(- fk/opts.g0));
      //gk = fk; --> use this to benchmark blim variations
      
      ff = exp(-gk*log(pow(m/opts.Q0,2)));

      //Gaussian term
      if (opts.g1 != 0)
        ff *= exp(-opts.g1/sqrt(1.+max(0.,opts.lambda)*pow(b,2)) * pow(abs(b),2));

      //Quartic term
      double lambdaq = 1;
      if (opts.q != 0)
      	ff *= exp(- (opts.q/abs(opts.q))* pow(lambdaq,4)  * (1. - exp(- abs(opts.q)*pow(real(b),4)/pow(lambdaq,4))));

      //Sextic term
      if (opts.sx > 0)
	ff *= exp(-opts.sx * pow(abs(b),6));
    }
  
  //Additional flavor dependent form factors
  if (opts.flavour_kt)
    {
      uvff = exp(-opts.g1_uv*pow(b,2));// * ff;
      usff = exp(-opts.g1_us*pow(b,2));// * ff;
      dvff = exp(-opts.g1_dv*pow(b,2));// * ff;
      dsff = exp(-opts.g1_ds*pow(b,2));// * ff;
      ssff = exp(-opts.g1_ss*pow(b,2));// * ff;
      chff = exp(-opts.g1_ch*pow(b,2));// * ff;
      boff = exp(-opts.g1_bo*pow(b,2));// * ff;
      glff = exp(-opts.g1_gl*pow(b,2));// * ff;
    }

  //universal form factor
  return ff;
}
