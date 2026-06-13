#include "sudakovff.h"
#include "resconst.h"
#include "constants.h"
#include "resint.h"
#include "blim.h"
#include "scales.h"
#include "besselint.h"
#include "pegasus.h"
#include "npff.h"
#include "isnan.h"
#include "gint.h"
#include "alphas.h"
#include "bstar.h"
#include "btilde.h"
#include "pdf.h"
#include "Li2.hpp"
#include <iostream>
#include <iomanip>

using namespace std;
using namespace resconst;
using namespace resint;
using namespace polylogarithm;

const double pi2 = M_PI*M_PI;
const double pi3 = pi2*M_PI;
const double pi4 = pi2*pi2;

//int sudakov::nf;
//double sudakov::beta0;
//double sudakov::beta1;
//double sudakov::beta2;
complex <double> sudakov::log1y;
complex <double> sudakov::log1yp;
//double sudakov::ry;
complex <double> sudakov::logS;
complex <double> sudakov::S;

complex <double> sudakov::S_bb;
complex <double> sudakov::S_cc;

//void sudakov::setnf(int nff)
//{
//  nf = nff;
//  beta0=(33.-2.*nf)/12.;
//  beta1=(153.-19.*nf)/24.;
//  beta2=2857./128.-5033.*nf/1152.+325.*nf*nf/3456.;
//}


complex <double> sudakov::f0(complex <double> y)
{
  //FFN beta
  //return (resconst::A1q/beta0(fabs(b)))*(y+log1y)/(y);
  //VFN beta
  //cout << "f0" << "  " << A1q << "  " << beta0 << "  " << y << " " << log1y << "  " << (A1q/beta0)*(y+log1y)/(y) << endl;

  return (A1q/beta0)*(y+log1y)/(y);
}

complex <double> sudakov::f1(complex <double> y) //Eq.23 of arXiv:hep-ph/0508068
{
  complex <double> F1 = 
    ((A1q*beta1)/(pow(beta0,3)))*((1./2.)*log1y*log1y +
					    (y)/(1.-y)+log1y/(1.-y)) -
    (A2q/(pow(beta0,2)))*(log1y+(y)/(1.-y)) + 
    (B1q/beta0)*log1y +
    (A1q/beta0)*(y/(1.-y)+log1y)*resint::rlogq2mur2; // !!! should be rlogq2mur2 -> log(mures^2/mur^2) = (logmur2q2-2.*loga) (is this a bug in DYRES?)
    //    a dependence      
  F1=F1-2.*resint::rloga*A1q/beta0*y/(1.-y);
  return F1;
}
complex <double> sudakov::f2(complex <double> y)
{
  complex <double> F2 = 
    ((A2q*beta1)/(pow(beta0,3)))*((y/2.)*((3.*y-2.)/pow(1.-y,2))-
				  ((1.-2.*y)*log1y/(1.-y)/(1.-y)))
    -(B2q/beta0)*((y)/(1.-y))
    +(B1q*beta1/pow(beta0,2))*((y)/(1.-y)+log1y/(1.-y))
    -(A3q/2./pow(beta0,2))*(y)*(y)/(1.-y)/(1.-y)
    +A1q*((pow(beta1,2)/2./pow(beta0,4))*(1.-2.*y)/(1.-y)
	  /(1.-y)*log1y*log1y +
	  log1y*((beta0*beta2-pow(beta1,2))/(pow(beta0,4))+
		 pow(beta1,2)/pow(beta0,4)/(1.-y)) +
	  // pow(beta1,4)/pow(beta0,4)/(1.-y)) +
	  (y)/(2.*pow(beta0,4)*(1.-y)*(1.-y))*
	  (beta0*beta2*(2.-3.*y)+pow(beta1,2)*y)) -
    (A1q/2.)*(y)*(y)/(1.-y)/(1.-y)*resint::rlogq2mur2*resint::rlogq2mur2 + //should be rlogq2mur2 -> log(mures^2/mur^2) = (logmur2q2-2.*loga) (is this a bug in DYRES?)
    resint::rlogq2mur2*(B1q*y/(1.-y)+A2q/beta0*y*y/(1.-y)/(1.-y)+
			A1q*beta1/pow(beta0,2)*(y/(1.-y)+(1.-2.*y)/(1.-y)/(1.-y)*log1y))
    +2.*C1qqn*((y)/(1.-y));
    
  //    a dependence (now without constant term)
  F2 += 2*A1q*y*(y-2.)/pow(1.-y,2)*pow(resint::rloga,2)-resint::rloga
    *(2.*B1q*y/(1.-y)+2.*y/beta0*A2q/pow(1.-y,2)
      -2.*A1q*beta1/pow(beta0,2)*y*log1y/pow(1.-y,2))
    + A1q*resint::rloga*resint::rlogq2mur2*y*2./pow(1.-y,2);
    
    return F2;
}

//Delta IS
complex <double> sudakov::g1deltaismb(complex <double> y, complex <double> b)
{
  complex <double> deltais = 0.;  
  if (pdf::mb > b0/real(b))
    {
      double lamR = beta0*aass*log(pow(scales::res/pdf::mb,2));
      if (opts.deltaismb)
	deltais = -4./3.*A1q/beta0/y*(y-lamR+(1.-lamR)*log( (1.-y)/(1.-lamR) ) );
    }
  return deltais;
}

complex <double> sudakov::g1deltaismc(complex <double> y, complex <double> b)
{
  complex <double> deltais = 0.;  
  if (pdf::mc > b0/real(b))
    {
      double lamR = beta0*aass*log(pow(scales::res/pdf::mc,2));
      if (opts.deltaismc)
	deltais = -4./3.*A1q/beta0/y*(y-lamR+(1.-lamR)*log( (1.-y)/(1.-lamR) ) );
    }
  return deltais;
}

//Delta FS
complex <double> sudakov::g1deltafs(complex <double> b)
{
  complex <double> deltafs = 0.;  
  double beta0_nf5 = alphas::bet0[5-alphas::NFMIN];
  double beta0_nf4 = alphas::bet0[4-alphas::NFMIN];
  double beta0_nf3 = alphas::bet0[3-alphas::NFMIN];
      
  if (opts.deltafsmb)
    if (pdf::mb > b0/real(b))
      {
	double asmb = pdf::alphas(pdf::mb);
	complex <double> blog_mb = log(pow(pdf::mb*b/b0,2));
	complex <double> y5_mb = beta0_nf5*asmb/M_PI*blog_mb;
	complex <double> y4_mb = beta0_nf4*asmb/M_PI*blog_mb;
	complex <double> g15_mb = (A1q*(y5_mb + log(1.-y5_mb)))/(beta0_nf5*y5_mb);
	complex <double> g14_mb = (A1q*(y4_mb + log(1.-y4_mb)))/(beta0_nf4*y4_mb);
	deltafs += (blog_mb*(g14_mb - g15_mb) - log(pow(pdf::mb/scales::res,2))*A1q/beta0_nf5*(beta0_nf5/beta0_nf4*log(1.-y4_mb)-log(1.-y5_mb)));
      }
  
  if (opts.deltafsmc)
    if (pdf::mc > b0/real(b))
      {
	double asmc = pdf::alphas(pdf::mc);
	complex <double> blog_mc = log(pow(pdf::mc*b/b0,2));
	complex <double> y4_mc = beta0_nf4*asmc/M_PI*blog_mc;
	complex <double> y3_mc = beta0_nf3*asmc/M_PI*blog_mc;
	complex <double> g14_mc = (A1q*(y4_mc + log(1.-y4_mc)))/(beta0_nf4*y4_mc);
	complex <double> g13_mc = (A1q*(y3_mc + log(1.-y3_mc)))/(beta0_nf3*y3_mc);
	deltafs += (blog_mc*(g13_mc - g14_mc) - log(pow(pdf::mc/scales::res,2))*A1q/beta0_nf4*(beta0_nf4/beta0_nf3*log(1.-y3_mc)-log(1.-y4_mc)));
      }
  
  return deltafs;
}


complex <double> sudakov::g1(complex <double> y)
{
  // Higgs production
  if (opts.nproc == 4) {
    A1q = A1g;
  }

  if (abs(y) < 1e-16)
    return A1q/beta0;
  
  return (A1q*(y + log1y))/(beta0*y);
  //return (A1q*(1. + log1y/y))/(beta0);
  
  //return (A1q*(ry + log1y))/(beta0*ry);
}
complex <double> sudakov::g1p(complex <double> yp)
{

  if (abs(yp) < 1e-16)
    return A1qp/beta0p;
  
  return (A1qp*(yp + log1yp))/(beta0p*yp);
  //return (A1q*(1. + log1y/y))/(beta0);
  
  //return (A1q*(ry + log1y))/(beta0*ry);
}
complex <double> sudakov::g2(complex <double> y)
{
  // Higgs production
  if (opts.nproc == 4) {

    A1q = A1g;
    A2q = A2g;
    B1q = B1g;

  }

  double B1qbar = B1q + A1q*LQ;
  //double LQR = 0.;
  double LQR = LR-LQ;
  return (A2q*(-(y/(1. - y)) - log1y))/pow(beta0,2) + (B1qbar*log1y)/beta0 + (A1q*beta1*(y/(1. - y) + log1y/(1. - y) + pow(log1y,2)/2.))/pow(beta0,3)
    + (A1q*(-y + (-1. + y)*log1y)*(LQR))/(beta0*(-1. + y));

  //return
  //  + (A2q*(-(ry/(1. - y)) - log1y))/pow(beta0,2)
  //  + (B1qbar*log1y)/beta0
  //  + (A1q*beta1*(ry/(1. - y) + log1y/(1. - y) + pow(log1y,2)/2.))/pow(beta0,3)
  //  + (A1q*(-ry + (-1. + ry)*log1y)*(LQR))/(beta0*(-1. + y));

}
complex <double> sudakov::g2p(complex <double> yp)
{
  double B1qbarp = B1qp + A1qp*LQ;
  //double LQR = 0.;
  double LQR = LR-LQ;
  return (A2qp*(-(yp/(1. - yp)) - log1yp))/pow(beta0p,2) + (B1qbarp*log1yp)/beta0p + (A1qp*beta1p*(yp/(1. - yp) + log1yp/(1. - yp) + pow(log1yp,2)/2.))/pow(beta0p,3)
    + (A1qp*(-yp + (-1. + yp)*log1yp)*(LQR))/(beta0p*(-1. + yp));
}

// delta function for W boson production . in dyqt deltaWQED was called deltaf0Wqed
complex <double> sudakov::deltaWQED(complex <double> yp)
{
	return D1w_QED/beta0p*log1yp;	
}


// dyqtg11 is the first mixed QCD+QED contribution at LL 
complex <double> sudakov::g11(complex <double> y, complex <double> yp, double eq2)
{
  double LQR = LR-LQ;
  complex <double> R1=y/yp;
  complex <double> R2=yp/y;
  // El coeficiente A11=0 en QCD-QED mixed                                                                                                  
  return (A1q*beta01_QCD)/(beta0*beta0*beta0p)*((log(R2/(R2-1.))+log(1.-y))*log(1.-R2*y)+(R2*(y-1.)*log(1.-y)+(1.-R2*y)*log(1.-R2*y))/((R2-1.)*(y-1.))-Li2(1./(1.-R2))+Li2((R2*y-1.)/(R2-1.)))+(A1qp*beta01_QED*eq2)/(beta0*beta0p*beta0p)*((log(R1/(R1-1.))+log(1.-yp))*log(1.-R1*yp)+(R1*(yp-1.)*log(1.-yp)+(1.-R1*yp)*log(1.-R1*yp))/((R1-1.)*(yp-1.))-Li2(1./(1.-R1))+Li2((R1*yp-1.)/(R1-1.)));
}
complex <double> sudakov::g3(complex <double> y)
{
  // Higgs procution
  if (opts.nproc == 4) {

    A1q = A1g;
    A2q = A2g;
    A3q = A3g;
    B1q = B1g;
    B2q = B2g;

    double B1qbar = B1q + A1q*LQ;
    double B2qbar = B2q + A2q*LQ;
    //double LQR = 0.;
    double LQR = LR-LQ;
    return (A1q*beta2*y)/(pow(beta0,3)*pow(1. - y,2)) - (3.*A1q*beta2*pow(y,2))/(2.*pow(beta0,3)*pow(1. - y,2)) - (A3q*pow(y,2))/(2.*pow(beta0,2)*pow(-1. + y,2)) + (A2q*beta1*y*(-2. + 3.*y))/(2.*pow(beta0,3)*pow(1. - y,2)) - (B2qbar*y)/(beta0 - beta0*y) + (A1q*beta2*log1y)/(pow(beta0,3)*pow(1. - y,2)) - (2.*A1q*beta2*y*log1y)/(pow(beta0,3)*pow(1. - y,2)) + (A1q*beta2*pow(y,2)*log1y)/(pow(beta0,3)*pow(1. - y,2)) + (A2q*beta1*(-1. + 2.*y)*log1y)/(pow(beta0,3)*pow(1. - y,2)) + (B1qbar*beta1*(y + log1y))/(pow(beta0,2)*(1. - y)) - (A1q*pow(beta1,2)*(y + log1y)*(-y + (-1. + 2.*y)*log1y))/(2.*pow(beta0,4)*pow(-1. + y,2))
    - ((LQR)*(2.*A1q*beta1*(-1. + 2.*y)*log1y + y*(2.*B1qbar*pow(beta0,2)*(-1. + y) + 2.*A1q*beta1*(-1. + y) - 2.*A2q*beta0*y + A1q*pow(beta0,2)*y*(LQR))))/(2.*pow(beta0,2)*pow(-1. + y,2))
    + 2.*C1ggn*((y)/(1.-y));

  }

  double B1qbar = B1q + A1q*LQ;
  double B2qbar = B2q + A2q*LQ;
  //double LQR = 0.;
  double LQR = LR-LQ;
  complex <double> G3 = (A1q*beta2*y)/(pow(beta0,3)*pow(1. - y,2)) - (3.*A1q*beta2*pow(y,2))/(2.*pow(beta0,3)*pow(1. - y,2)) - (A3q*pow(y,2))/(2.*pow(beta0,2)*pow(-1. + y,2)) + (A2q*beta1*y*(-2. + 3.*y))/(2.*pow(beta0,3)*pow(1. - y,2)) - (B2qbar*y)/(beta0 - beta0*y) + (A1q*beta2*log1y)/(pow(beta0,3)*pow(1. - y,2)) - (2.*A1q*beta2*y*log1y)/(pow(beta0,3)*pow(1. - y,2)) + (A1q*beta2*pow(y,2)*log1y)/(pow(beta0,3)*pow(1. - y,2)) + (A2q*beta1*(-1. + 2.*y)*log1y)/(pow(beta0,3)*pow(1. - y,2)) + (B1qbar*beta1*(y + log1y))/(pow(beta0,2)*(1. - y)) - (A1q*pow(beta1,2)*(y + log1y)*(-y + (-1. + 2.*y)*log1y))/(2.*pow(beta0,4)*pow(-1. + y,2))
    - ((LQR)*(2.*A1q*beta1*(-1. + 2.*y)*log1y + y*(2.*B1qbar*pow(beta0,2)*(-1. + y) + 2.*A1q*beta1*(-1. + y) - 2.*A2q*beta0*y + A1q*pow(beta0,2)*y*(LQR))))/(2.*pow(beta0,2)*pow(-1. + y,2));

  //double bound = log(2.)/resint::aass;
  ////if (opts.gbound)
  //if (real(G3) > bound)
  //  G3 = complex <double> (bound, imag(G3));

  return G3;
  
  //return
  //  +(A1q*beta2*ry)/(pow(beta0,3)*pow(1. - y,2))
  //  -(3.*A1q*beta2*pow(ry,2))/(2.*pow(beta0,3)*pow(1. - y,2))
  //  -(A3q*pow(ry,2))/(2.*pow(beta0,2)*pow(-1. + y,2))
  //  +(A2q*beta1*ry*(-2. + 3.*ry))/(2.*pow(beta0,3)*pow(1. - y,2))
  //  -(B2qbar*ry)/(beta0 - beta0*y)
  //  +(A1q*beta2*log1y)/(pow(beta0,3)*pow(1. - y,2))
  //  -(2.*A1q*beta2*ry*log1y)/(pow(beta0,3)*pow(1. - y,2))
  //  +(A1q*beta2*pow(ry,2)*log1y)/(pow(beta0,3)*pow(1. - y,2))
  //  +(A2q*beta1*(-1. + 2.*ry)*log1y)/(pow(beta0,3)*pow(1. - y,2))
  //  +(B1qbar*beta1*(ry + log1y))/(pow(beta0,2)*(1. - y))
  //  -(A1q*pow(beta1,2)*(ry + log1y)*(-ry + (-1. + 2.*ry)*log1y))/(2.*pow(beta0,4)*pow(-1. + y,2))
  //  
  //  - (
  //     (LQR)*(
  //	      + 2.*A1q*beta1*(-1. + 2.*ry)*log1y
  //	      + ry*(
  //		   + 2.*B1qbar*pow(beta0,2)*(-1. + y)
  //		   + 2.*A1q*beta1*(-1. + y)
  //		   - 2.*A2q*beta0*ry
  //		   + A1q*pow(beta0,2)*ry*(LQR)
  //		   )
  //	      )
  //     )
  //  /(2.*pow(beta0,2)*pow(-1. + y,2));
  
}
complex <double> sudakov::g4(complex <double> y)
{
  // CJU - N3LL
  // Higgs procution
  if (opts.nproc == 4) {

    A1q = A1g;
    A2q = A2g;
    A3q = A3g;
    A4q = A4g;
    B1q = B1g;
    B2q = B2g;
    B3q = B3g;

    double B1qbar = B1q + A1q*LQ;
    double B2qbar = B2q + A2q*LQ;
    double B3qbar = B3q + A3q*LQ;

    //double LQR = 0.;
    double LQR = LR-LQ;

    bool term = true;  //switch off one term which gives divergent Sudakov
    return (B3qbar*(-2. + y)*y)/(2.*beta0*pow(1. - y,2)) + (A4q*(-3. + y)*pow(y,2))/(6.*pow(beta0,2)*pow(1. - y,3)) - (A3q*beta1*y*(6. - 15.*y + 5.*pow(y,2)))/(12.*pow(beta0,3)*pow(1. - y,3)) + (A3q*beta1*(-1. + 3.*y)*log1y)/(2.*pow(beta0,3)*pow(1. - y,3)) + (B2qbar*beta1*((2. - y)*y + 2.*log1y))/(2.*pow(beta0,2)*pow(1. - y,2)) + (B1qbar*((pow(beta1,2) - beta0*beta2)*pow(y,2) - pow(beta1,2)*pow(log1y,2)))/(2.*pow(beta0,3)*pow(1. - y,2)) - (A2q*(y*(8.*beta0*beta2*pow(y,2) + pow(beta1,2)*(-6. + (9. - 11.*y)*y)) + 6.*pow(beta1,2)*log1y*(-1. + y + (-1. + 3.*y)*log1y)))/(12.*pow(beta0,4)*pow(1. - y,3)) - (A1q*(y*(2.*pow(beta1,3)*pow(y,2) + pow(beta0,2)*beta3*(-6. + (15. - 7.*y)*y) + beta0*beta1*beta2*(6. + 5.*(-3. + y)*y)) + 2.*log1y*(3.*(-(pow(beta0,2)*beta3*pow(1. - y,3)) + pow(beta1,3)*pow(y,2)*(1. + y) + beta0*beta1*beta2*(1. + y*(-3. + 2.*(1. - y)*y))) + pow(beta1,3)*log1y*(3.*y + (1. - 3.*y)*log1y))))/(12.*pow(beta0,5)*pow(1. - y,3))
    + ((LQR)*(3.*y*(-2.*B2qbar*pow(beta0,3)*(-2. + y)*(-1. + y) - A2q*beta0*beta1*(-2. + y)*(-1. + y) + y*(A3q*pow(beta0,2)*(-3. + y) + A1q*(pow(beta1,2) - beta0*beta2)*(1. + y))) + 3.*beta1*log1y*(2.*(-(B1qbar*pow(beta0,2)*(-1. + y)) + A1q*beta1*y + A2q*beta0*(-1. + 3.*y)) + term*A1q*beta1*(1. - 3.*y)*log1y) + 3.*pow(beta0,2)*(y*(A1q*beta1 + B1qbar*pow(beta0,2)*(-2. + y)*(-1. + y) - A2q*beta0*(-3. + y)*y) + A1q*beta1*(1. - 3.*y)*log1y)*(LQR) + A1q*pow(beta0,4)*(-3. + y)*pow(y,2)*pow(LQR,2)))/(6.*pow(beta0,3)*pow(-1. + y,3))
    + 2.*C1ggn*((y)/(1.-y));

  }

  double B1qbar = B1q + A1q*LQ;
  double B2qbar = B2q + A2q*LQ;
  double B3qbar = B3q + A3q*LQ;

  //double LQR = 0.;
  double LQR = LR-LQ;
  
  bool term = true;  //switch off one term which gives divergent Sudakov
  complex <double> G4 = (B3qbar*(-2. + y)*y)/(2.*beta0*pow(1. - y,2)) + (A4q*(-3. + y)*pow(y,2))/(6.*pow(beta0,2)*pow(1. - y,3)) - (A3q*beta1*y*(6. - 15.*y + 5.*pow(y,2)))/(12.*pow(beta0,3)*pow(1. - y,3)) + (A3q*beta1*(-1. + 3.*y)*log1y)/(2.*pow(beta0,3)*pow(1. - y,3)) + (B2qbar*beta1*((2. - y)*y + 2.*log1y))/(2.*pow(beta0,2)*pow(1. - y,2)) + (B1qbar*((pow(beta1,2) - beta0*beta2)*pow(y,2) - pow(beta1,2)*pow(log1y,2)))/(2.*pow(beta0,3)*pow(1. - y,2)) - (A2q*(y*(8.*beta0*beta2*pow(y,2) + pow(beta1,2)*(-6. + (9. - 11.*y)*y)) + 6.*pow(beta1,2)*log1y*(-1. + y + (-1. + 3.*y)*log1y)))/(12.*pow(beta0,4)*pow(1. - y,3)) - (A1q*(y*(2.*pow(beta1,3)*pow(y,2) + pow(beta0,2)*beta3*(-6. + (15. - 7.*y)*y) + beta0*beta1*beta2*(6. + 5.*(-3. + y)*y)) + 2.*log1y*(3.*(-(pow(beta0,2)*beta3*pow(1. - y,3)) + pow(beta1,3)*pow(y,2)*(1. + y) + beta0*beta1*beta2*(1. + y*(-3. + 2.*(1. - y)*y))) + pow(beta1,3)*log1y*(3.*y + (1. - 3.*y)*log1y))))/(12.*pow(beta0,5)*pow(1. - y,3))
    + ((LQR)*(3.*y*(-2.*B2qbar*pow(beta0,3)*(-2. + y)*(-1. + y) - A2q*beta0*beta1*(-2. + y)*(-1. + y) + y*(A3q*pow(beta0,2)*(-3. + y) + A1q*(pow(beta1,2) - beta0*beta2)*(1. + y))) + 3.*beta1*log1y*(2.*(-(B1qbar*pow(beta0,2)*(-1. + y)) + A1q*beta1*y + A2q*beta0*(-1. + 3.*y)) + term*A1q*beta1*(1. - 3.*y)*log1y) + 3.*pow(beta0,2)*(y*(A1q*beta1 + B1qbar*pow(beta0,2)*(-2. + y)*(-1. + y) - A2q*beta0*(-3. + y)*y) + A1q*beta1*(1. - 3.*y)*log1y)*(LQR) + A1q*pow(beta0,4)*(-3. + y)*pow(y,2)*pow(LQR,2)))/(6.*pow(beta0,3)*pow(-1. + y,3));

  //double bound = log(2.)/resint::aass/resint::aass;
  ////if (opts.gbound)
  //if (real(G4) > bound)
  //  G4 = complex <double> (bound, imag(G4));

  return G4;
  
  //return 
  //  +(B3qbar*(-2. + ry)*ry)/(2.*beta0*pow(1. - y,2))
  //  +(A4q*(-3. + ry)*pow(ry,2))/(6.*pow(beta0,2)*pow(1. - y,3))
  //  -(A3q*beta1*ry*(6. - 15.*ry + 5.*pow(ry,2)))/(12.*pow(beta0,3)*pow(1. - y,3))
  //  +(A3q*beta1*(-1. + 3.*ry)*log1y)/(2.*pow(beta0,3)*pow(1. - y,3))    
  //  +(B2qbar*beta1*((2. - ry)*ry + 2.*log1y))/(2.*pow(beta0,2)*pow(1. - y,2))
  //  +(B1qbar*((pow(beta1,2) - beta0*beta2)*pow(ry,2) - pow(beta1,2)*pow(log1y,2)))/(2.*pow(beta0,3)*pow(1. - y,2))
  //  -(A2q*(ry*(8.*beta0*beta2*pow(ry,2) + pow(beta1,2)*(-6. + (9. - 11.*ry)*ry)) + 6.*pow(beta1,2)*log1y*(-1. + y + (-1. + 3.*ry)*log1y)))/(12.*pow(beta0,4)*pow(1. - y,3))
  //  -(A1q*(ry*(2.*pow(beta1,3)*pow(ry,2) + pow(beta0,2)*beta3*(-6. + (15. - 7.*ry)*ry) + beta0*beta1*beta2*(6. + 5.*(-3. + ry)*ry)) + 2.*log1y*(3.*(-(pow(beta0,2)*beta3*pow(1. - y,3)) + pow(beta1,3)*pow(ry,2)*(1. + ry) + beta0*beta1*beta2*(1. + ry*(-3. + 2.*(1. - y)*ry))) + pow(beta1,3)*log1y*(3.*ry + (1. - 3.*ry)*log1y))))/(12.*pow(beta0,5)*pow(1. - y,3))
  //
  //  + (
  //     (LQR)*(
  //	      +3.*ry*(-2.*B2qbar*pow(beta0,3)*(-2. + ry)*(-1. + y)
  //		    - A2q*beta0*beta1*(-2. + ry)*(-1. + y)
  //		    + ry*(A3q*pow(beta0,2)*(-3. + ry)
  //			 + A1q*(pow(beta1,2) - beta0*beta2)*(1. + y)))
  //	      + 3.*beta1*log1y*(
  //				+2.*(0.
  //				     -(B1qbar*pow(beta0,2)*(-1. + y))
  //				     + A1q*beta1*ry
  //				     + A2q*beta0*(-1. + 3.*ry)
  //				    )
  //				+ term*A1q*beta1*(1. - 3.*ry)*log1y
  //				)
  //	      + 3.*pow(beta0,2)*(ry*(A1q*beta1 + B1qbar*pow(beta0,2)*(-2. + ry)*(-1. + y)
  //				    - A2q*beta0*(-3. + ry)*ry)
  //				 + A1q*beta1*(1. - 3.*ry)*log1y)*(LQR)
  //	      + A1q*pow(beta0,4)*(-3. + ry)*pow(ry,2)*pow(LQR,2)
  //	      )
  //     )/(6.*pow(beta0,3)*pow(-1. + y,3))
  //  ;
}


complex <double> sudakov::g5(complex <double> y)
{
  double B1qbar = B1q + A1q*LQ;
  double B2qbar = B2q + A2q*LQ;
  double B3qbar = B3q + A3q*LQ;
  double B4qbar = B4q + A4q*LQ;

  //double LQR = 0.;
  double LQR = LR-LQ;

  complex <double> G5 =

    -(B4qbar*y*(3. + (-3. + y)*y))/(3.*beta0*pow(1. - y,3)) - (A5q*pow(y,2)*(6. - 4.*y + pow(y,2)))/(12.*pow(beta0,2)*pow(1. - y,4)) + (A4q*beta1*y*(-12. + 42.*y - 28.*pow(y,2) + 7.*pow(y,3)))/(36.*pow(beta0,3)*pow(1. - y,4)) + (A4q*beta1*(-1. + 4.*y)*log1y)/(3.*pow(beta0,3)*pow(1. - y,4)) + (B3qbar*beta1*(y*(3. + (-3. + y)*y) + 3.*log1y))/(3.*pow(beta0,2)*pow(1. - y,3))
    //- (B2qbar*pow(beta1,2)*(2.*y*(3. + (-3. + y)*y) + 6.*log1y + 9.*pow(log1y,2)))/(27.*pow(beta0,3)*pow(1. - y,3))
    + (A3q*(9.*beta0*beta2*(-4. + y)*pow(y,3) + pow(beta1,2)*y*(12. - 24.*y + 52.*pow(y,2) - 13.*pow(y,3)) + 12.*pow(beta1,2)*(1. - y)*log1y + 18.*pow(beta1,2)*(1. - 4.*y)*pow(log1y,2)))/(36.*pow(beta0,4)*pow(1. - y,4)) + (A2q*(3.*pow(beta0,2)*beta3*pow(y,3)*(-8. + 5.*y) + 2.*beta0*beta1*beta2*y*(6. - 21.*y + 44.*pow(y,2) - 20.*pow(y,3)) + pow(beta1,3)*y*(-12. + 42.*y - 64.*pow(y,2) + 25.*pow(y,3)) - 12.*beta1*(pow(beta1,2) - beta0*beta2)*(1. - 4.*y + 9.*pow(y,2))*log1y - 36.*pow(beta1,3)*y*pow(log1y,2) + 12.*pow(beta1,3)*(-1. + 4.*y)*pow(log1y,3)))/(36.*pow(beta0,5)*pow(1. - y,4)) + (A1q*(y*(-3.*pow(beta1,4)*pow(y,2)*(2. + 3.*y) + pow(beta0,3)*beta4*(12. - 42.*y + 40.*pow(y,2) - 13.*pow(y,3)) + beta0*pow(beta1,2)*beta2*(12. - 42.*y + 52.*pow(y,2) + 5.*pow(y,3)) + pow(beta0,2)*(2.*beta1*beta3*(-6. + 21.*y - 17.*pow(y,2) + 5.*pow(y,3)) + pow(beta2,2)*(-12. + 42.*y - 52.*pow(y,2) + 7.*pow(y,3)))) + 6.*(2.*pow(beta0,3)*beta4*pow(1. - y,4) + pow(beta1,4)*pow(y,2)*(-3. + 2.*y - 2.*pow(y,2)) + 2.*beta0*pow(beta1,2)*beta2*(1. - 4.*y + 6.*pow(y,2) - 6.*pow(y,3) + 3.*pow(y,4)) - pow(beta0,2)*(2.*pow(beta2,2)*pow(1. - y,4) + beta1*beta3*(2. - 8.*y + 9.*pow(y,2) - 10.*pow(y,3) + 4.*pow(y,4))))*log1y + 18.*pow(beta1,2)*y*(-3.*beta0*beta2*y + pow(beta1,2)*(-1. + 3.*y))*pow(log1y,2) + 6.*pow(beta1,4)*(-1. + 6.*y)*pow(log1y,3) + 3.*pow(beta1,4)*(1. - 4.*y)*pow(log1y,4)))/(36.*pow(beta0,6)*pow(1. - y,4)) - (B1qbar*((pow(beta1,3) - 2.*beta0*beta1*beta2 + pow(beta0,2)*beta3)*(3. - 2.*y)*pow(y,2) + beta1*log1y*(6.*(pow(beta1,2) - beta0*beta2)*y + pow(beta1,2)*(3. - 2.*log1y)*log1y)))/(6.*pow(beta0,4)*pow(1. - y,3))

    //+(B2qbar*(y*(-9.*beta0*beta2*(-3. + y)*y + pow(beta1,2)*(-6. + 7.*(-3. + y)*y)) +  6.*pow(beta1,2)*log1y*(-1. + 3.*log1y)))/(27.*pow(beta0,3)*pow(-1. + y,3))    

    +(B2qbar*((pow(beta1,2) - beta0*beta2)*(-3. + y)*pow(y,2) + 3.*pow(beta1,2)*pow(log1y,2)))/ (3.*pow(beta0,3)*pow(-1. + y,3))
    
    -((LQR)*(4.*A1q*pow(beta1,3)*(-1. + 4.*y)*pow(log1y,3) - 6.*pow(beta1,2)*pow(log1y,2)*(2.*B1qbar*pow(beta0,2) - A1q*beta1 - 2.*B1qbar*pow(beta0,2)*y + 6.*A1q*beta1*y + 2.*A2q*beta0*(-1. + 4.*y) + A1q*pow(beta0,2)*(1. - 4.*y)*(LQR)) + 2.*beta1*log1y*(6.*(B1qbar*pow(beta0,2)*beta1 - 2.*B2qbar*pow(beta0,3)*(-1. + y) + 2.*A2q*beta0*beta1*y - B1qbar*pow(beta0,2)*beta1*y + A1q*pow(beta1,2)*y - 3.*A1q*pow(beta1,2)*pow(y,2) + 3.*A1q*beta0*beta2*pow(y,2) + A3q*pow(beta0,2)*(-1. + 4.*y)) + 3.*pow(beta0,2)*(A2q*beta0*(2. - 8.*y) + A1q*beta1*(1. - 6.*y) + 2.*B1qbar*pow(beta0,2)*(-1. + y))*(LQR) + 2.*A1q*pow(beta0,4)*(-1. + 4.*y)*pow(LQR,2)) + y*(2.*(6.*B1qbar*pow(beta0,2)*pow(beta1,2) - 6.*B1qbar*pow(beta0,3)*beta2 - 12.*A4q*pow(beta0,3)*y + 6.*A2q*beta0*pow(beta1,2)*y - 6.*B1qbar*pow(beta0,2)*pow(beta1,2)*y + 3.*A1q*pow(beta1,3)*y - 6.*A2q*pow(beta0,2)*beta2*y + 6.*B1qbar*pow(beta0,3)*beta2*y - 3.*A1q*pow(beta0,2)*beta3*y + 8.*A4q*pow(beta0,3)*pow(y,2) + 8.*A2q*beta0*pow(beta1,2)*pow(y,2) - 2.*A1q*pow(beta1,3)*pow(y,2) - 8.*A2q*pow(beta0,2)*beta2*pow(y,2) + 4.*A1q*beta0*beta1*beta2*pow(y,2) - 2.*A1q*pow(beta0,2)*beta3*pow(y,2) - 2.*A4q*pow(beta0,3)*pow(y,3) - 2.*A2q*beta0*pow(beta1,2)*pow(y,3) + 2.*A1q*pow(beta1,3)*pow(y,3) + 2.*A2q*pow(beta0,2)*beta2*pow(y,3) - 4.*A1q*beta0*beta1*beta2*pow(y,3) + 2.*A1q*pow(beta0,2)*beta3*pow(y,3) + 6.*B3qbar*pow(beta0,4)*(-3. + 6.*y - 4.*pow(y,2) + pow(y,3)) + 2.*A3q*pow(beta0,2)*beta1*(-3. + 6.*y - 4.*pow(y,2) + pow(y,3))) - 6.*pow(beta0,2)*(-2.*A2q*beta0*beta1 - 3.*B1qbar*pow(beta0,2)*beta1 - A1q*pow(beta1,2) - 6.*A3q*pow(beta0,2)*y + 6.*B1qbar*pow(beta0,2)*beta1*y + 3.*A1q*pow(beta1,2)*y - 3.*A1q*beta0*beta2*y + 4.*A3q*pow(beta0,2)*pow(y,2) - 4.*B1qbar*pow(beta0,2)*beta1*pow(y,2) - A3q*pow(beta0,2)*pow(y,3) + B1qbar*pow(beta0,2)*beta1*pow(y,3) + 2.*B2qbar*pow(beta0,3)*(-3. + 6.*y - 4.*pow(y,2) + pow(y,3)))*(LQR) + 2.*pow(beta0,4)*(-2.*A2q*beta0*y*(6. - 4.*y + pow(y,2)) + 2.*B1qbar*pow(beta0,2)*(-3. + 6.*y - 4.*pow(y,2) + pow(y,3)) - A1q*beta1*(2. + 6.*y - 4.*pow(y,2) + pow(y,3)))*pow(LQR,2) + A1q*pow(beta0,6)*y*(6. - 4.*y + pow(y,2))*pow(LQR,3))))/(12.*pow(beta0,4)*pow(-1. + y,4))
    ;

  //double bound = log(2.)/resint::aass/resint::aass/resint::aass;
  ////if (opts.gbound)
  //if (real(G5) > bound)
  //  G5 = complex <double> (bound, imag(G5));

  return G5;
}



    

//Sudakov form factor QCD ONLY
complex <double> sudakov::sff(complex <double> b)
{
  //Return 0 if b is close to the Landau pole
  double b_L;// = a_param_.b0p_*(1./resint::_m)*exp(1./(2.*resint::aass*resconst::beta0));
  if (opts.modlog)
    //b_L = b0/scales::res * (exp(1./(2.*aass*beta0))-1.);
    //I think the correct formula is:
    b_L = b0/scales::res * sqrt(exp(1./(aass*beta0))-1.);
  else
    b_L = b0/scales::res * exp(1./(2.*aass*beta0));

  if (abs(b-b_L) < 1e-8)
    return 0;

  //In reading these formulas, notice that blog is L = log((q*bstar/b0p)^2) = log[(q/a_param)*bstar/b0] = log[Q^2 * bstar^2/b0^2], according to Eq. (13) and (17) of hep-ph/0508068.
      
  //choose bstar (b) for real axis (complex plane) integration
  //mass dependence in blim --> mass dependence in bstar
  //complex <double> blim = cx(blimit_.cblim_); //blim=(1/q)*exp(1/(2*aass*beta0))
  double blim = blim::sudakov;
  complex <double> bstar;
  if (opts.bprescription == 0 || opts.bprescription == 4 || opts.bstar_sudakov)
    //bstar = real(b)/sqrt(1.+pow(real(b)/blim,2));
    bstar = bstar::calc(b,blim);
  else //if (opts.bprescription == 1 || opts.bprescription == 2 || opts.bprescription == 3)
    bstar = b;

  //cout << "bstar " << bstar << endl;
  
  //mass dependence (q) in blog
  //complex <double> blog;
  //if (opts.modlog)
  //  blog = log(pow(scaleh_.q_*bstar/a_param_.b0p_,2)+1.); //modified sudakov
  //else
  //  blog = log(pow(scaleh_.q_*bstar/a_param_.b0p_,2));    //normal sudakov
  
  double Q = scales::res;
  complex <double> blog;
  if (!opts.modlog)
    blog = log(pow(Q*bstar/b0,2));   //normal sudakov
  else if (opts.p == 1)
    blog = log(pow(Q*bstar/b0,2) + 1.); //modified sudakov
  else
    blog = 1./opts.p*log(pow(Q*bstar/b0,2*opts.p) + 1.); //modified sudakov with exponent p
    //blog = log(pow(pow(Q*bstar/b0,2*opts.p) + 1.,1./opts.p)); //modified sudakov with exponent p
    //blog = log(pow(Q*bstar/b0,2) + double(opts.p)); //modified sudakov

  
  //if (opts.modlog)
  //  blog = log(pow(Q*bstar/b0,2) + 1.); //modified sudakov
  //else
  //  blog = log(pow(Q*bstar/b0,2));   //normal sudakov

  //double rblog;
  //if (opts.modlog)
  //  rblog = log(pow(Q*real(bstar)/b0,2) + 1.); //modified sudakov
  //else
  //  rblog = log(pow(Q*real(bstar)/b0,2));   //normal sudakov
  
  //cout << modified_.imod_ << endl;
  //cout << "blog " << blog << endl;

  complex <double> y = beta0*resint::aass*blog;

  //Numerical integration of the Sudakov form factor
  if (opts.numsud)
    logS = gint::logS;

  //Analytical solution for the Sudakov form factor
  else
    {
      if (abs(y) < 1e-16)
	{
	  S = 1.;
	  logS = 0.;
	  return S;
	}
      
      //mass dependence in f0(y), f1(y), f2(y)
      //FFN beta
      //complex <double> y = resconst::beta0*aass_.aass_*blog;
      //VFN beta
      log1y = log(1.-y);
      csud_.log1y_ = fcx(log(1.-y));

      //ry = beta0*resint::aass*rblog;
      //ry = abs(y);
      //ry = beta0*resint::aass*real(blog);
      
      fcomplex fy = fcx(y);

      logS = 0.;


      //                                blog*f0(beta0*as*blog) + f1(beta0*as*blog) + ...
      //      (1d0+cc/(q*bstar/b0p)**2)*blog*f0(beta0*as*blog) + f1(beta0*as*blog) + ...

      double cc = 0.5;
      if (opts.order_sudak >= 0)
	//if (opts.powcorr)
	//  logS += (1.+cc/pow(Q*bstar/b0,2))*blog*g1(y);
	//else
	logS += blog*g1(y);
      if (opts.order_sudak >= 1)
	logS += g2(y);
      if (opts.order_sudak >= 2)
	logS += aass*g3(y);
      if (opts.order_sudak >= 3)
	logS += pow(aass,2)*g4(y);
      if (opts.order_sudak >= 4)
	logS += pow(aass,3)*g5(y);

      /*
      if (opts.order_sudak == 0)
	//S = exp(blog*cx(f0_(fy)));
	//S = exp(blog*f0(y));
	S = exp(blog*g1(y));
      else if (opts.order_sudak == 1)
	//S = exp(blog*cx(f0_(fy))+cx(f1_(fy)));
	//S = exp(blog*f0(y)+f1(y));
	S = exp(blog*g1(y)+g2(y));
      else if (opts.order_sudak == 2)
	//S = exp(blog*cx(f0_(fy))+cx(f1_(fy))+aass_.aass_*cx(f2_(fy)));
	//S = exp(blog*f0(y)+f1(y)+resint::aass*f2(y));
	//S = exp(blog*g1(y)+g2(y)+resint::aass*(g3(y)+2.*C1qqn*((y)/(1.-y))));
	S = exp(blog*g1(y)+g2(y)+resint::aass*g3(y));
      else if (opts.order_sudak == 3)
	//S = exp(blog*g1(y)+g2(y)+resint::aass*g3(y));
	S = exp(blog*g1(y)+g2(y)+resint::aass*g3(y)+pow(resint::aass,2)*g4(y));
      else if (opts.order_sudak == 4)
	S = exp(blog*g1(y)+g2(y)+resint::aass*g3(y)+pow(resint::aass,2)*g4(y)+pow(resint::aass,3)*g5(y));
      */
      
      /*
	if (opts.order_sudak == 0)
	S = exp(blog*g1(y));
	else if (opts.order_sudak == 1)
	S = exp(blog*g1(y)+g2(y));
	else if (opts.order_sudak == 2)
	S = exp(blog*g1(y)+g2(y)+resint::aass*g3(y));
	else if (opts.order_sudak == 3)
	S = exp(blog*g1(y)+g2(y)+resint::aass*g3(y)+pow(resint::aass,2)*g4(y));
      */

      //cout << " blog " << blog << " g1 " << g1(y) << " g2 " << g2(y) << endl;
      //  cout << " y " << y << "  " << cx(fy) << endl;

      //cout << setprecision(16);
      //cout << "f0(y) " << cx(f0_(fy)) << "  " << f0(y) << endl;
      //cout << "f1(y) " << cx(f1_(fy)) << "  " << f1(y) << endl;
      //cout << "f2(y) " << cx(f2_(fy)) << "  " << f2(y) << endl;

      //cout << setprecision(16);
      //cout << "f0(y) " << f0(y) << " g1 " << g1(y) << endl;
      //cout << "f1(y) " << f1(y) << " g2 " << g2(y) << endl;
      //cout << "f2(y) " << f2(y) << " g3 " << g3(y) << endl;
      //cout << "f3(y) " << 0     << " g4 " << g4(y) << endl;

      //compare analytical and numerical integrations
      //cout << setprecision(16);
      //cout << " b " << b << endl;
      //cout << " logS " << logS << endl;
      //cout << " gint " << gint::logS << endl;
      //cout << endl;
    }

  if (opts.deltafsmb || opts.deltafsmc)
    logS += g1deltafs(bstar);
      
  if (opts.deltaismb)
    S_bb = exp(blog*g1deltaismb(y, bstar));

  if (opts.deltaismc)
    S_cc = exp(blog*g1deltaismc(y, bstar));

  /*
  ofstream bf("bline.C");
  bf << "{" << endl;
  bf << "TGraph *gb = new TGraph();" << endl;
  int nb = 1000;
  double bm = 3.;
  for (int ib = 0; ib < nb; ib++)
    {
      double bbb = ib * bm/nb;
      complex <double> blog = btilde::blog(bbb,Q);
      complex <double> y = beta0*resint::aass*blog;
      //cout << bbb << "  " << g1deltafs(bbb) << endl;
      bf << "gb->SetPoint(gb->GetN(), " << bbb << ", " << real(g1deltafs(bbb)) << ");" << endl;
      //bf << "gb->SetPoint(gb->GetN(), " << bbb << ", " << real(g1deltaismb(y, bbb)) << ");" << endl;
      //bf << "gb->SetPoint(gb->GetN(), " << bbb << ", " << real(g1deltaismc(y, bbb)) << ");" << endl;
    }
  bf << "gb->Draw();" << endl;
  bf << "}" << endl;
  exit(0);
  */
  
  S = exp(logS);

  if (opts.bprescription == 0 && abs(S) > 100.)
    {
      fcomplex fbb = fcx(b);
      //cout << "Warning! Large Sudakov, S(b) = " << S <<"; for bstar = " << bstar << " fortran S(b) = " << cx(s_(fbb)) << endl;
      //S = 0.;
      //cout << "C++ " << b << "  " << bstar << "  " << blog << "  " << y << "  " << blog*f0(y) << "  " << f1(y) << "  " << resint::aass*f2(y) << "  " <<  S << endl;
      //cout << "C++ " << b << "  " << bstar << "  " << blog << "  " << y << "  " << blog*g1(y) << "  " << g2(y) << "  " << resint::aass*g3(y) << "  " << pow(resint::aass,2)*g4(y) << "  " << S << endl;
    }
  if (isnan_ofast(real(S)) || isnan_ofast(imag(S)))
    {
      fcomplex fbb = fcx(b);
      cout << "Warning! Sudakov is nan, S(b) = " << S << " log(S) = " << logS << "; for b = " << b << " bstar = " << bstar << " Landau pole " << b_L << endl;// << " fortran S(b) = " << cx(s_(fbb)) << endl;
      cout << "C++ " << b << "  " << bstar << "  " << blog << "  " << y << "  " << blog*g1(y) << "  "
	   << g2(y) << "  " << resint::aass*g3(y) << "  " << pow(resint::aass,2)*g4(y) << "  " << pow(resint::aass,3)*g5(y)
	   << "  " <<  S << endl;
      //cout << "C++ " << b << "  " << bstar << "  " << blog << "  " << y << "  " << blog*f0(y) << "  " << f1(y) << "  " << resint::aass*f2(y) << "  " <<  S << endl;
      //should distinguish between positive divergence (-> set S = 1) and negative divergence (-> set S = 0)
      S = 0;
      //S = 1;
    }
  //cout << "C++ " << b << "  " << bstar << "  " << y << "  " << f0(y) << "  " << S << endl;
  
  //S = S*exp(-opts.g_param*pow(b,2));
  //S = S * npff::S(b,resint::_m,resint::x1,resint::x2);

  if (opts.sumlogs)
    return 1.;
  else
    return S;
}







//Sudakov form factor QCD+QED (Z boson)
complex <double> sudakov::sff(complex <double> b, double eq2)
{
  //Return 0 if b is close to the Landau pole
  double b_L;// = a_param_.b0p_*(1./resint::_m)*exp(1./(2.*resint::aass*resconst::beta0));
  if (opts.modlog)
    //b_L = b0/scales::res * (exp(1./(2.*aass*beta0))-1.);
    //I think the correct formula is:
    b_L = b0/scales::res * sqrt(exp(1./(aass*beta0))-1.);
  else
    b_L = b0/scales::res * exp(1./(2.*aass*beta0));

  if (abs(b-b_L) < 1e-8)
    return 0;

  //In reading these formulas, notice that blog is L = log((q*bstar/b0p)^2) = log[(q/a_param)*bstar/b0] = log[Q^2 * bstar^2/b0^2], according to Eq. (13) and (17) of hep-ph/0508068.
      
  //choose bstar (b) for real axis (complex plane) integration
  //mass dependence in blim --> mass dependence in bstar
  //complex <double> blim = cx(blimit_.cblim_); //blim=(1/q)*exp(1/(2*aass*beta0))
  double blim = blim::sudakov;
  complex <double> bstar;
  if (opts.bprescription == 0 || opts.bprescription == 4 || opts.bstar_sudakov)
    //bstar = real(b)/sqrt(1.+pow(real(b)/blim,2));
    bstar = bstar::calc(b,blim);
  else //if (opts.bprescription == 1 || opts.bprescription == 2 || opts.bprescription == 3)
    bstar = b;

  //cout << "bstar " << bstar << endl;
  
  //mass dependence (q) in blog
  //complex <double> blog;
  //if (opts.modlog)
  //  blog = log(pow(scaleh_.q_*bstar/a_param_.b0p_,2)+1.); //modified sudakov
  //else
  //  blog = log(pow(scaleh_.q_*bstar/a_param_.b0p_,2));    //normal sudakov
  
  double Q = scales::res;

  /*
  complex <double> blog;
  if (!opts.modlog)
    blog = log(pow(Q*bstar/b0,2));   //normal sudakov
  else if (opts.p == 1)
    blog = log(pow(Q*bstar/b0,2) + 1.); //modified sudakov
  else
    blog = 1./opts.p*log(pow(Q*bstar/b0,2*opts.p) + 1.); //modified sudakov with exponent p
    //blog = log(pow(pow(Q*bstar/b0,2*opts.p) + 1.,1./opts.p)); //modified sudakov with exponent p
    //blog = log(pow(Q*bstar/b0,2) + double(opts.p)); //modified sudakov
    */
  
  //if (opts.modlog)
  //  blog = log(pow(Q*bstar/b0,2) + 1.); //modified sudakov
  //else
  //  blog = log(pow(Q*bstar/b0,2));   //normal sudakov

  //double rblog;
  //if (opts.modlog)
  //  rblog = log(pow(Q*real(bstar)/b0,2) + 1.); //modified sudakov
  //else
  //  rblog = log(pow(Q*real(bstar)/b0,2));   //normal sudakov
  
  //cout << modified_.imod_ << endl;
  //cout << "blog " << blog << endl;

  complex <double> blog = btilde::blog(bstar,Q);
  complex <double> y = beta0*resint::aass*blog;
  complex <double> yp = beta0p*resint::aaqed*blog;

  //  cout << " aass : " << resint::aass << endl;
  //  cout << " resconst::aaqed  " << resint::aaqed << endl;

  //Numerical integration of the Sudakov form factor
  if (opts.numsud)
    { 
      logS = gint::logS;
      S = exp(logS);
    }
  //Analytical solution for the Sudakov form factor
  else
    {
      if (abs(y) < 1e-16)
	{
	  S = 1.;
	  logS = 0.;
	  return S;
	}
      
      //mass dependence in f0(y), f1(y), f2(y)
      //FFN beta
      //complex <double> y = resconst::beta0*aass_.aass_*blog;
      //VFN beta
      log1y = log(1.-y);
      csud_.log1y_ = fcx(log(1.-y));
      log1yp = log(1.-yp);
      csud_.log1yp_ = fcx(log(1.-yp));

      //ry = beta0*resint::aass*rblog;
      //ry = abs(y);
      //ry = beta0*resint::aass*real(blog);
      
      fcomplex fy = fcx(y);
      fcomplex fyp = fcx(yp);

      logS = 0.;


      //                                blog*f0(beta0*as*blog) + f1(beta0*as*blog) + ...
      //      (1d0+cc/(q*bstar/b0p)**2)*blog*f0(beta0*as*blog) + f1(beta0*as*blog) + ...
      //   cout << "andrea " << " opts.qedorder " << opts.qedorder << endl;
      double cc = 0.5;
// Andrea consideration: it is useless to put here anothere if opts.QED; this overloading of the function is caled only if we pass the electric charge!       
      if (opts.order_sudak >= 0)
	{
	  //if (opts.powcorr)
	  //  logS += (1.+cc/pow(Q*bstar/b0,2))*blog*g1(y);
	  //else
	  if (opts.qed) logS += blog*(g1(y)+eq2*g1p(yp)); 
	  else logS += blog*g1(y);
	}
      if (opts.order_sudak >= 1)
	{
	  if (opts.qed && (opts.qedorder_sudak >= 1)) 
	  {
	  	logS += g2(y)+eq2*g2p(yp)+g11(y,yp,eq2);//Since g11 is pure mixed the charge is passed in the argument. In the case of pure QED (as for g2p) the charge completle factorise.
	  }	
	  else  logS += g2(y);
	}
      if (opts.order_sudak >= 2)
	logS += aass*g3(y);
      if (opts.order_sudak >= 3)
	logS += pow(aass,2)*g4(y);
      if (opts.order_sudak >= 4)
	logS += pow(aass,3)*g5(y);

      S = exp(logS);

      /*
      if (opts.order_sudak == 0)
	//S = exp(blog*cx(f0_(fy)));
	//S = exp(blog*f0(y));
	S = exp(blog*g1(y));
      else if (opts.order_sudak == 1)
	//S = exp(blog*cx(f0_(fy))+cx(f1_(fy)));
	//S = exp(blog*f0(y)+f1(y));
	S = exp(blog*g1(y)+g2(y));
      else if (opts.order_sudak == 2)
	//S = exp(blog*cx(f0_(fy))+cx(f1_(fy))+aass_.aass_*cx(f2_(fy)));
	//S = exp(blog*f0(y)+f1(y)+resint::aass*f2(y));
	//S = exp(blog*g1(y)+g2(y)+resint::aass*(g3(y)+2.*C1qqn*((y)/(1.-y))));
	S = exp(blog*g1(y)+g2(y)+resint::aass*g3(y));
      else if (opts.order_sudak == 3)
	//S = exp(blog*g1(y)+g2(y)+resint::aass*g3(y));
	S = exp(blog*g1(y)+g2(y)+resint::aass*g3(y)+pow(resint::aass,2)*g4(y));
      else if (opts.order_sudak == 4)
	S = exp(blog*g1(y)+g2(y)+resint::aass*g3(y)+pow(resint::aass,2)*g4(y)+pow(resint::aass,3)*g5(y));
      */
      
      /*
	if (opts.order_sudak == 0)
	S = exp(blog*g1(y));
	else if (opts.order_sudak == 1)
	S = exp(blog*g1(y)+g2(y));
	else if (opts.order_sudak == 2)
	S = exp(blog*g1(y)+g2(y)+resint::aass*g3(y));
	else if (opts.order_sudak == 3)
	S = exp(blog*g1(y)+g2(y)+resint::aass*g3(y)+pow(resint::aass,2)*g4(y));
      */

      //cout << " blog " << blog << " g1 " << g1(y) << " g2 " << g2(y) << endl;
      //  cout << " y " << y << "  " << cx(fy) << endl;

      //cout << setprecision(16);
      //cout << "f0(y) " << cx(f0_(fy)) << "  " << f0(y) << endl;
      //cout << "f1(y) " << cx(f1_(fy)) << "  " << f1(y) << endl;
      //cout << "f2(y) " << cx(f2_(fy)) << "  " << f2(y) << endl;

      //cout << setprecision(16);
      //cout << "f0(y) " << f0(y) << " g1 " << g1(y) << endl;
      //cout << "f1(y) " << f1(y) << " g2 " << g2(y) << endl;
      //cout << "f2(y) " << f2(y) << " g3 " << g3(y) << endl;
      //cout << "f3(y) " << 0     << " g4 " << g4(y) << endl;

      //compare analytical and numerical integrations
      //cout << setprecision(16);
      //cout << " b " << b << endl;
      //cout << " logS " << logS << endl;
      //cout << " gint " << gint::logS << endl;
      //cout << endl;
    }

  if (opts.bprescription == 0 && abs(S) > 100.)
    {
      fcomplex fbb = fcx(b);
      cout << "Warning! Large Sudakov, S(b) = " << S <<"; for bstar = " << bstar << " fortran S(b) = " << cx(s_(fbb)) << endl;
      //S = 0.;
      //cout << "C++ " << b << "  " << bstar << "  " << blog << "  " << y << "  " << blog*f0(y) << "  " << f1(y) << "  " << resint::aass*f2(y) << "  " <<  S << endl;
      //cout << "C++ " << b << "  " << bstar << "  " << blog << "  " << y << "  " << blog*g1(y) << "  " << g2(y) << "  " << resint::aass*g3(y) << "  " << pow(resint::aass,2)*g4(y) << "  " << S << endl;
    }
  if (isnan_ofast(real(S)) || isnan_ofast(imag(S)))
    {
      fcomplex fbb = fcx(b);
      cout << "Warning! Sudakov is nan, S(b) = " << S << " log(S) = " << logS << "; for b = " << b << " bstar = " << bstar << " Landau pole " << b_L << endl;// << " fortran S(b) = " << cx(s_(fbb)) << endl;
      cout << "C++ " << b << "  " << bstar << "  " << blog << "  " << y << "  " << blog*g1(y) << "  "
	   << g2(y) << "  " << resint::aass*g3(y) << "  " << pow(resint::aass,2)*g4(y) << "  " << pow(resint::aass,3)*g5(y)
	   << "  " <<  S << endl;
      //cout << "C++ " << b << "  " << bstar << "  " << blog << "  " << y << "  " << blog*f0(y) << "  " << f1(y) << "  " << resint::aass*f2(y) << "  " <<  S << endl;
      //should distinguish between positive divergence (-> set S = 1) and negative divergence (-> set S = 0)
      S = 0;
      //S = 1;
    }
  //cout << "C++ " << b << "  " << bstar << "  " << y << "  " << f0(y) << "  " << S << endl;
  
  //S = S*exp(-opts.g_param*pow(b,2));
  //S = S * npff::S(b,resint::_m,resint::x1,resint::x2);

  if (opts.sumlogs)
    return 1.;
  else
    return S;
}

// Andrea : trial , beginnning, for W boson production  
complex <double> sudakov::sff(complex <double> b, double eq2, double eqbarp2)
{
  //Return 0 if b is close to the Landau pole
  double b_L;// = a_param_.b0p_*(1./resint::_m)*exp(1./(2.*resint::aass*resconst::beta0));
  if (opts.modlog)
    //b_L = b0/scales::res * (exp(1./(2.*aass*beta0))-1.);
    //I think the correct formula is:
    b_L = b0/scales::res * sqrt(exp(1./(aass*beta0))-1.);
  else
    b_L = b0/scales::res * exp(1./(2.*aass*beta0));

  if (abs(b-b_L) < 1e-8)
    return 0;

  //In reading these formulas, notice that blog is L = log((q*bstar/b0p)^2) = log[(q/a_param)*bstar/b0] = log[Q^2 * bstar^2/b0^2], according to Eq. (13) and (17) of hep-ph/0508068.
      
  //choose bstar (b) for real axis (complex plane) integration
  //mass dependence in blim --> mass dependence in bstar
  //complex <double> blim = cx(blimit_.cblim_); //blim=(1/q)*exp(1/(2*aass*beta0))
  double blim = blim::sudakov;
  complex <double> bstar;
  if (opts.bprescription == 0 || opts.bprescription == 4 || opts.bstar_sudakov)
    //bstar = real(b)/sqrt(1.+pow(real(b)/blim,2));
    bstar = bstar::calc(b,blim);
  else //if (opts.bprescription == 1 || opts.bprescription == 2 || opts.bprescription == 3)
    bstar = b;

  //cout << "bstar " << bstar << endl;
  
  //mass dependence (q) in blog
  //complex <double> blog;
  //if (opts.modlog)
  //  blog = log(pow(scaleh_.q_*bstar/a_param_.b0p_,2)+1.); //modified sudakov
  //else
  //  blog = log(pow(scaleh_.q_*bstar/a_param_.b0p_,2));    //normal sudakov
  
  double Q = scales::res;

  /*
  complex <double> blog;
  if (!opts.modlog)
    blog = log(pow(Q*bstar/b0,2));   //normal sudakov
  else if (opts.p == 1)
    blog = log(pow(Q*bstar/b0,2) + 1.); //modified sudakov
  else
    blog = 1./opts.p*log(pow(Q*bstar/b0,2*opts.p) + 1.); //modified sudakov with exponent p
    //blog = log(pow(pow(Q*bstar/b0,2*opts.p) + 1.,1./opts.p)); //modified sudakov with exponent p
    //blog = log(pow(Q*bstar/b0,2) + double(opts.p)); //modified sudakov
    */
  
  //if (opts.modlog)
  //  blog = log(pow(Q*bstar/b0,2) + 1.); //modified sudakov
  //else
  //  blog = log(pow(Q*bstar/b0,2));   //normal sudakov

  //double rblog;
  //if (opts.modlog)
  //  rblog = log(pow(Q*real(bstar)/b0,2) + 1.); //modified sudakov
  //else
  //  rblog = log(pow(Q*real(bstar)/b0,2));   //normal sudakov
  
  //cout << modified_.imod_ << endl;
  //cout << "blog " << blog << endl;

  complex <double> blog = btilde::blog(bstar,Q);
  complex <double> y = beta0*resint::aass*blog;
  complex <double> yp = beta0p*resint::aaqed*blog;

  //  cout << " aass : " << resint::aass << endl;
  //  cout << " resconst::aaqed  " << resint::aaqed << endl;

  //Numerical integration of the Sudakov form factor
  if (opts.numsud)
    { 
      logS = gint::logS;
      S = exp(logS);
    }
  //Analytical solution for the Sudakov form factor
  else
    {
      if (abs(y) < 1e-16)
	{
	  S = 1.;
	  logS = 0.;
	  return S;
	}
      
      //mass dependence in f0(y), f1(y), f2(y)
      //FFN beta
      //complex <double> y = resconst::beta0*aass_.aass_*blog;
      //VFN beta
      log1y = log(1.-y);
      csud_.log1y_ = fcx(log(1.-y));
      log1yp = log(1.-yp);
      csud_.log1yp_ = fcx(log(1.-yp));

      //ry = beta0*resint::aass*rblog;
      //ry = abs(y);
      //ry = beta0*resint::aass*real(blog);
      
      fcomplex fy = fcx(y);
      fcomplex fyp = fcx(yp);

      logS = 0.;


      //                                blog*f0(beta0*as*blog) + f1(beta0*as*blog) + ...
      //      (1d0+cc/(q*bstar/b0p)**2)*blog*f0(beta0*as*blog) + f1(beta0*as*blog) + ...
      //   cout << "andrea " << " opts.qedorder " << opts.qedorder << endl;
      double cc = 0.5;
// Andrea consideration: it is useless to put here anothere if opts.QED; this overloading of the function is caled only if we pass the electric charge!       
      if (opts.order_sudak >= 0)
	{
	  //if (opts.powcorr)
	  //  logS += (1.+cc/pow(Q*bstar/b0,2))*blog*g1(y);
	  //else
	  if (opts.qed) logS += blog*(g1(y)+eq2*g1p(yp)); 
	  else logS += blog*g1(y);
	}
      if (opts.order_sudak >= 1)
	{
	  if (opts.qed && (opts.qedorder_sudak >= 1)) 
	  // Andrea proposal: insertion of deltaWQED and changing of charge factors  
	  {
	  	logS += g2(y)+(eq2+eqbarp2)/2.*g2p(yp)+deltaWQED(yp)+0.5*(g11(y,yp,eq2)+g11(y,yp,eqbarp2));//Since g11 is pure mixed the charge is passed in the argument. In the case of pure QED (as for g2p) the charge completle factorise.
	  }	
	  else  logS += g2(y);
	}
      if (opts.order_sudak >= 2)
	logS += aass*g3(y);
      if (opts.order_sudak >= 3)
	logS += pow(aass,2)*g4(y);
      if (opts.order_sudak >= 4)
	logS += pow(aass,3)*g5(y);

      S = exp(logS);

      /*
      if (opts.order_sudak == 0)
	//S = exp(blog*cx(f0_(fy)));
	//S = exp(blog*f0(y));
	S = exp(blog*g1(y));
      else if (opts.order_sudak == 1)
	//S = exp(blog*cx(f0_(fy))+cx(f1_(fy)));
	//S = exp(blog*f0(y)+f1(y));
	S = exp(blog*g1(y)+g2(y));
      else if (opts.order_sudak == 2)
	//S = exp(blog*cx(f0_(fy))+cx(f1_(fy))+aass_.aass_*cx(f2_(fy)));
	//S = exp(blog*f0(y)+f1(y)+resint::aass*f2(y));
	//S = exp(blog*g1(y)+g2(y)+resint::aass*(g3(y)+2.*C1qqn*((y)/(1.-y))));
	S = exp(blog*g1(y)+g2(y)+resint::aass*g3(y));
      else if (opts.order_sudak == 3)
	//S = exp(blog*g1(y)+g2(y)+resint::aass*g3(y));
	S = exp(blog*g1(y)+g2(y)+resint::aass*g3(y)+pow(resint::aass,2)*g4(y));
      else if (opts.order_sudak == 4)
	S = exp(blog*g1(y)+g2(y)+resint::aass*g3(y)+pow(resint::aass,2)*g4(y)+pow(resint::aass,3)*g5(y));
      */
      
      /*
	if (opts.order_sudak == 0)
	S = exp(blog*g1(y));
	else if (opts.order_sudak == 1)
	S = exp(blog*g1(y)+g2(y));
	else if (opts.order_sudak == 2)
	S = exp(blog*g1(y)+g2(y)+resint::aass*g3(y));
	else if (opts.order_sudak == 3)
	S = exp(blog*g1(y)+g2(y)+resint::aass*g3(y)+pow(resint::aass,2)*g4(y));
      */

      //cout << " blog " << blog << " g1 " << g1(y) << " g2 " << g2(y) << endl;
      //  cout << " y " << y << "  " << cx(fy) << endl;

      //cout << setprecision(16);
      //cout << "f0(y) " << cx(f0_(fy)) << "  " << f0(y) << endl;
      //cout << "f1(y) " << cx(f1_(fy)) << "  " << f1(y) << endl;
      //cout << "f2(y) " << cx(f2_(fy)) << "  " << f2(y) << endl;

      //cout << setprecision(16);
      //cout << "f0(y) " << f0(y) << " g1 " << g1(y) << endl;
      //cout << "f1(y) " << f1(y) << " g2 " << g2(y) << endl;
      //cout << "f2(y) " << f2(y) << " g3 " << g3(y) << endl;
      //cout << "f3(y) " << 0     << " g4 " << g4(y) << endl;

      //compare analytical and numerical integrations
      //cout << setprecision(16);
      //cout << " b " << b << endl;
      //cout << " logS " << logS << endl;
      //cout << " gint " << gint::logS << endl;
      //cout << endl;
    }

  if (opts.bprescription == 0 && abs(S) > 100.)
    {
      fcomplex fbb = fcx(b);
      cout << "Warning! Large Sudakov, S(b) = " << S <<"; for bstar = " << bstar << " fortran S(b) = " << cx(s_(fbb)) << endl;
      //S = 0.;
      //cout << "C++ " << b << "  " << bstar << "  " << blog << "  " << y << "  " << blog*f0(y) << "  " << f1(y) << "  " << resint::aass*f2(y) << "  " <<  S << endl;
      //cout << "C++ " << b << "  " << bstar << "  " << blog << "  " << y << "  " << blog*g1(y) << "  " << g2(y) << "  " << resint::aass*g3(y) << "  " << pow(resint::aass,2)*g4(y) << "  " << S << endl;
    }
  if (isnan_ofast(real(S)) || isnan_ofast(imag(S)))
    {
      fcomplex fbb = fcx(b);
      cout << "Warning! Sudakov is nan, S(b) = " << S << " log(S) = " << logS << "; for b = " << b << " bstar = " << bstar << " Landau pole " << b_L << endl;// << " fortran S(b) = " << cx(s_(fbb)) << endl;
      cout << "C++ " << b << "  " << bstar << "  " << blog << "  " << y << "  " << blog*g1(y) << "  "
	   << g2(y) << "  " << resint::aass*g3(y) << "  " << pow(resint::aass,2)*g4(y) << "  " << pow(resint::aass,3)*g5(y)
	   << "  " <<  S << endl;
      //cout << "C++ " << b << "  " << bstar << "  " << blog << "  " << y << "  " << blog*f0(y) << "  " << f1(y) << "  " << resint::aass*f2(y) << "  " <<  S << endl;
      //should distinguish between positive divergence (-> set S = 1) and negative divergence (-> set S = 0)
      S = 0;
      //S = 1;
    }
  //cout << "C++ " << b << "  " << bstar << "  " << y << "  " << f0(y) << "  " << S << endl;
  
  //S = S*exp(-opts.g_param*pow(b,2));
  //S = S * npff::S(b,resint::_m,resint::x1,resint::x2);

  if (opts.sumlogs)
    return 1.;
  else
    return S;
}
