#include "npff.h"
#include "settings.h"
#include "constants.h"
#include "pdf.h"
#include "pdfevol.h"
#include "blim.h"

using namespace constants;

complex <double> npff::gk(complex <double> b)
{
  double fgk = 0.;

  double bb;
  //bb = real(b);
  bb = abs(b);
  
  if (opts.npff == 0)
    {
      double g = max(0.,opts.g2);
      g = g/2.; //Divide by a factor of 2 because the BLNY formula is g2*ln(Q/Q0) while we use g2*ln(Q^2/Q0^2)
      fgk = g * pow(bb,2);
    }

  //Exponential form
  if (opts.npff == 1)
    {
      fgk = 0.;
    }
  
  //Collins, Rogers, large bT form (https://arxiv.org/abs/1412.3820 Eq.(39), (79))
  //g0 = 0.3; //Q0 = 1.6;
  else if (opts.npff == 2)
    {
      //double gk = opts.g0 * (1. - exp(- (CF*pdf::extalphas(real(pdfevol::mubstar))*pow(bb,2))/(M_PI*opts.g0*pow(blim::sudakov,2))));
      //double gk = opts.g0 * (1. - exp(- (CF*pdf::alphas(real(pdfevol::mubstar),opts.cmw)*pow(bb,2))/(M_PI*opts.g0*pow(blim::sudakov,2))));
      //cout << "Collins-Rogers " << " mubstar " << pdfevol::mubstar << " alphas " << pdf::alphas(real(pdfevol::mubstar),opts.cmw) << " lha alphas " << pdf::extalphas(real(pdfevol::mubstar)) << endl;

      fgk = opts.g0 * (1. - exp(- (CF*pdf::alphas(real(pdfevol::mubstar))*pow(bb,2))/(M_PI*opts.g0*pow(blim::sudakov,2)) +
				(opts.qCS/abs(opts.qCS)) * (1 - exp(-abs(opts.qCS)*pow(bb,4)))));

      /*
      double as = pdf::alphas(real(pdfevol::mubstar),opts.cmw)/M_PI;
      double x = pow(real(b/blim::sudakov),2);
      //double gk = opts.g0 * (1. - exp(- resconst::A1q*as*x/opts.g0));
      double fk = 0.;
      double gk = opts.g0 * (1. - exp(- (resconst::A1q*as+resconst::A2q*as*as)*x/opts.g0));
      //double fk = opts.g0 * (1. - exp(- resconst::B1q*as*x/opts.g0));
      ff = exp(-gk*log(pow(m/opts.Q0,2)) - fk);
      */
    }
  
  return fgk;
}
