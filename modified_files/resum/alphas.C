#include "alphas.h"
#include "resconst.h"
#include "resint.h"
#include "coupling.h"
#include "pdf.h"
#include "scales.h"
#include "settings.h"
#include "constants.h"
#include "phasespace.h"
#include "pegasus.h"
#include "crundec_interface.h"

using namespace resconst;
using namespace constants;
using namespace resint;

const int alphas::nastps = 10;

complex <double> alphas::as;
double alphas::as0;
complex <double> alphas::asLO;
complex <double> alphas::asNLO;
complex <double> alphas::asNNLO;
complex <double> alphas::asNNNLO;
complex <double> alphas::asNNNNLO;
complex <double> alphas::as1_1l;
complex <double> alphas::as1_2l;
complex <double> alphas::as1_3l;
complex <double> alphas::as1_4l;
complex <double> alphas::as1_5l;
complex <double> alphas::as2_2l;
complex <double> alphas::as2_3l;
complex <double> alphas::as2_4l;
complex <double> alphas::as2_5l;
complex <double> alphas::as3_3l;
complex <double> alphas::as3_4l;
complex <double> alphas::as3_5l;
complex <double> alphas::as4_4l;
complex <double> alphas::as4_5l;
complex <double> alphas::as5_5l;

double alphas::bet0[NFMAX-NFMIN+1];
double alphas::bet1[NFMAX-NFMIN+1];
double alphas::bet2[NFMAX-NFMIN+1];
double alphas::bet3[NFMAX-NFMIN+1];
double alphas::bet4[NFMAX-NFMIN+1];

double alphas::beta0_msbar;
double alphas::beta1_msbar;
double alphas::beta2_msbar;
double alphas::beta3_msbar;
double alphas::beta4_msbar;

double alphas::brnf[NFMAX-NFMIN+1][6][nmax];

double alphas::b1csi;
double alphas::b2csi;
double alphas::b3csi;
double alphas::b4csi;
double alphas::b5csi;
double alphas::b6csi;

//double alphas::beta0;
//double alphas::beta1;
//double alphas::beta2;
//double alphas::beta3;
//double alphas::beta4;

complex <double> (*alphas::fbetabarptr)(complex <double> a) = 0;

void alphas::init()
{
  for (int nf = NFMIN; nf <= NFMAX; nf++)
    {
      bet0[nf-NFMIN] = (33.-2.*nf)/12.;
      bet1[nf-NFMIN] = (153.-19.*nf)/24.;
      bet2[nf-NFMIN] = 2857./128.-5033.*nf/1152.+325.*nf*nf/3456.;
      bet3[nf-NFMIN] = (149753./6. + 3564.*zeta3
			 + nf*(-1078361./162.-6508./27.*zeta3)
			 + nf*nf*(50065./162.+6472./81.*zeta3)
			 + nf*nf*nf*1093./729.)/256.; //from https://arxiv.org/pdf/1701.01404.pdf Eq.(3.6), with a 4^4=256 normalisation factor
      bet4[nf-NFMIN] = (8157455./16. + 621885./2.*zeta3 - 88209./2.*zeta4 - 288090.*zeta5
	     +nf* (-336460813./1944. - 4811164./81.*zeta3 + 33935./6.*zeta4 + 1358995./27.*zeta5)
	     +nf*nf*(25960913./1944. + 698531./81.*zeta3 - 10526./9.*zeta4 - 381760./81.*zeta5)
	     +nf*nf*nf*(-630559./5832. - 48722./243.*zeta3 + 1618./27.*zeta4 + 460./9.*zeta5)
	     +nf*nf*nf*nf*(1205./2916. - 152./81.*zeta3))/1024.;  //from https://arxiv.org/pdf/1701.01404.pdf Eq.(3.7), with a 4^5=1024 normalisation factor
    }

  beta0_msbar = bet0[resconst::NF-NFMIN];
  beta1_msbar = bet1[resconst::NF-NFMIN];
  beta2_msbar = bet2[resconst::NF-NFMIN];
  beta3_msbar = bet3[resconst::NF-NFMIN];
  beta4_msbar = bet4[resconst::NF-NFMIN];

  //Cache flavour dependent brnf
  for (int nf = NFMIN; nf <= NFMAX; nf++)
    {
      double b[6][nmax] = {0.};
      for (int nl = 0; nl < 6; nl++)
	{
	  b[nl][0] = 1.;
	  if (nl >= 1) b[nl][1] = bet1[nf-NFMIN]/bet0[nf-NFMIN]/M_PI;
	  if (nl >= 2) b[nl][2] = bet2[nf-NFMIN]/bet0[nf-NFMIN]/pi2;
	  if (nl >= 3) b[nl][3] = bet3[nf-NFMIN]/bet0[nf-NFMIN]/pi3;
	  if (nl >= 4) b[nl][4] = bet4[nf-NFMIN]/bet0[nf-NFMIN]/pi4;
	}

      for (int nl = 0; nl < 6; nl++)
	{
	  brnf[nf-NFMIN][nl][0] =  b[nl][0];
	  brnf[nf-NFMIN][nl][1] = -b[nl][1];
	  for (int k = 2; k < nmax; k++)
	    for (int kp = 0; kp < k; kp++)
	      brnf[nf-NFMIN][nl][k] += - brnf[nf-NFMIN][nl][kp]*b[nl][k-kp];
	  for (int k = 2; k < nmax; k++)
	    brnf[nf-NFMIN][nl][k] /= double(k-1);
	}
    }
}

//Truncated solution for the running of alphas up to N4LO
complex <double> alphas::calc(complex <double> q, int nloop, bool commute)
{
  //cout << "alphas evolution, q = " << q << endl;
  
  //double q0 = phasespace::m;
  //double q0 = scales::ren;
  //double q0 = coupling::zmass;
  //as0 = pdf::alphas(q0);
  
  double q0 = scales::res;

  //if (opts.alphaslha)
  //  as0 = pdf::extalphas(scales::ren);
  //else
  //  as0 = pdf::rgktalphas(scales::ren);
  as0 = pdf::alphas(scales::ren);

  //cout << endl;
  //rgkt(q, q0, as0);
  //cout << "as Runge-Kutta C++:     " << asLO << "  " << asNLO << "  " << asNNLO << "  " << asNNNLO << "  " << asNNNNLO << endl;
  //expan(q, q0, as0);
  //cout << "as expanded       :     " << asLO << "  " << asNLO << "  " << asNNLO << "  " << asNNNLO << "  " << asNNNNLO << endl;
  //iter(q, q0, as0);
  //cout << "as iterative      :     " << asLO << "  " << asNLO << "  " << asNNLO << "  " << asNNNLO << "  " << asNNNNLO << endl;
  //anal(q, q0, as0);
  //cout << "as analytic       :     " << asLO << "  " << asNLO << "  " << asNNLO << "  " << asNNNLO << "  " << asNNNNLO << endl;
  
  /*
  //Runge Kutta evolution
  //  if (opts.asrgkt)
  //    {
  //q = q0+complex <double>(0.,q0);
  //q = complex <double>(0.1,10.);
  rgkt(q, q0, as0);
  cout << "as Runge-Kutta C++:     " << asLO << "  " << asNLO << "  " << asNNLO << "  " << asNNNLO << "  " << asNNNNLO << endl;
      
  //double R20  = pow(q0,2);
  //double R2   = pow(q,2);
  //double AS0 = as0/4./M_PI;
  //int nf = 5;
  //aspar_.naord_ = 0;    asLO = as_(R2, R20, AS0, nf)*4*M_PI;
  //aspar_.naord_ = 1;   asNLO = as_(R2, R20, AS0, nf)*4*M_PI;
  //aspar_.naord_ = 2;  asNNLO = as_(R2, R20, AS0, nf)*4*M_PI;
  //aspar_.naord_ = 3; asNNNLO = as_(R2, R20, AS0, nf)*4*M_PI;
  //aspar_.naord_  = order_.npord_;
  //cout << "as Runge-Kutta fortran: " << asLO << "  " << asNLO << "  " << asNNLO << "  " << asNNNLO << endl;
  //    }
  */
  
  //compute asN4LO numerically
  if (opts.order >= 4)
  rgkt(q, q0, as0);

  complex <double> qlog = log(pow(q/q0,2));

  complex <double> lambda = as0*beta0/M_PI*qlog; //-as0*beta0/M_PI*qlog;
  complex <double> oneplambda = 1. + lambda;
  complex <double> log1plambda = log(oneplambda);
  //complex <double> log1plambda = log1p(lambda); //check if this is faster

  /*
    if (nloop >= 0) //0-loop
    as = as0;
    if (nloop >= 1) //1-loop (LO)
    as = as0/(1.+lambda);
    if (nloop >= 2) //2-loop (NLO)
    as += -(pow(as0,2)*beta1*log1plambda)/(beta0*pow(oneplambda,2)*M_PI);
    if (nloop >= 3) //3-loop (NNLO)
    as += (pow(as0,3)*((pow(beta1,2) - beta0*beta2)*lambda + pow(beta1,2)*(-1.+log1plambda)*log1plambda))/(pow(beta0,2)*pow(oneplambda,3)*pi2);
    if (nloop >= 4) //4-loop (NNNLO)
    as += (pow(as0,4)*(-(lambda*(pow(beta1,3)*lambda - 2.*beta0*beta1*beta2*(oneplambda) + pow(beta0,2)*beta3*(2.+lambda))) + beta1*log1plambda*(-4.*pow(beta1,2)*lambda + 2.*beta0*beta2*(-1.+2.*lambda) + pow(beta1,2)*(5. - 2.*log1plambda)*log1plambda)))/(2.*pow(beta0,3)*pow(oneplambda,4)*pi2*M_PI);
  */
      
  asLO    = as0/(1.+lambda); //Taylor series: asLO = as0*(1.-lambda+pow(lambda,2)-pow(lambda,3)+pow(lambda,4)-pow(lambda,5)+pow(lambda,6));
  asNLO   = asLO   - (pow(as0,2)*beta1*log1plambda)/(beta0*pow(oneplambda,2)*M_PI);
  asNNLO  = asNLO  + (pow(as0,3)*((pow(beta1,2) - beta0*beta2)*lambda + pow(beta1,2)*(-1.+log1plambda)*log1plambda))/(pow(beta0,2)*pow(oneplambda,3)*pi2);
  asNNNLO = asNNLO + (pow(as0,4)*(-(lambda*(pow(beta1,3)*lambda - 2.*beta0*beta1*beta2*(oneplambda) + pow(beta0,2)*beta3*(2.+lambda)))
				  + beta1*log1plambda*(-4.*pow(beta1,2)*lambda + 2.*beta0*beta2*(-1.+2.*lambda)
						       + pow(beta1,2)*(5. - 2.*log1plambda)*log1plambda)))/(2.*pow(beta0,3)*pow(oneplambda,4)*pi2*M_PI);
  asNNNNLO = asNNNLO + (pow(as0,5)*(lambda*(-6.*beta0*pow(beta1,2)*beta2*lambda*(4. + lambda) + pow(beta1,4)*lambda*(9. + 2.*lambda) - 2.*pow(beta0,3)*beta4*(3. + 3.*lambda + pow(lambda,2))
					    + pow(beta0,2)*(2.*pow(beta2,2)*lambda*(6. + lambda) + beta1*beta3*(6. + 9.*lambda + 4.*pow(lambda,2))))
				    + 6.*beta1*(beta0*beta1*beta2*(1. - 2.*lambda)*lambda + pow(beta1,3)*(-3. + lambda)*lambda + pow(beta0,2)*beta3*(-1. + 2.*lambda + pow(lambda,2)))*log1plambda
				    + 9.*(-2.*beta0*pow(beta1,2)*beta2*(-1. + lambda) + pow(beta1,4)*(1. + 2.*lambda))*pow(log1plambda,2)
				    - 26.*pow(beta1,4)*pow(log1plambda,3) + 6.*pow(beta1,4)*pow(log1plambda,4)))/(6.*pow(beta0,4)*pow(oneplambda,5)*pi4);



  //Approximate the exact numeric alphas with higher-order iterative analytic solutions setting beta terms to zero
  //complex <double> asiter;
  //asiter   = as0/(1.+lambda); //Taylor series: asLO = as0*(1.-lambda+pow(lambda,2)-pow(lambda,3)+pow(lambda,4)-pow(lambda,5)+pow(lambda,6));
  //asiter   += - (pow(as0,2)*beta1*log1plambda)/(beta0*pow(oneplambda,2)*M_PI);
  //asiter   += + (pow(as0,3)*((pow(beta1,2))*lambda + pow(beta1,2)*(-1.+log1plambda)*log1plambda))/(pow(beta0,2)*pow(oneplambda,3)*pi2);
  //asiter   += + (pow(as0,4)*(-(lambda*(pow(beta1,3)*lambda))
  //				  + beta1*log1plambda*(-4.*pow(beta1,2)*lambda
  //						       + pow(beta1,2)*(5. - 2.*log1plambda)*log1plambda)))/(2.*pow(beta0,3)*pow(oneplambda,4)*pi2*M_PI);
  //cout << "as iterative             " << asLO << "  " << asiter << endl;
  
  complex <double> dasNLO   = asNLO-asLO;
  complex <double> dasNNLO  = asNNLO-asNLO;
  complex <double> dasNNNLO = asNNNLO-asNNLO;

  //Loop-truncated powers before scale variations
  //To calculate loop-truncate powers write asN^kLO = asLO + dasNLO + dasNNLO + ...
  //Then expand the expression pow(asN^kLO,p), and take the terms which contribute up to k loop,
  //where asLO contributes 1 loop, dasNLO contributes 2 loops, etc...
  //alphas truncated powers of 2
  as2_2l = pow(asLO,2);
  as2_3l = pow(asLO,2)+2.*dasNLO*asLO;
  as2_4l = pow(asLO,2)+pow(dasNLO,2)+2.*asLO*dasNLO+2.*asLO*dasNNLO;
  as2_5l = pow(asLO,2)+pow(dasNLO,2)+2.*asLO*dasNLO+2.*asLO*dasNNLO+2.*dasNLO*dasNNLO+2.*asLO*dasNNNLO;
  //alphas truncated powers of 3
  as3_3l = pow(asLO,3);
  as3_4l = pow(asLO,3)+3.*pow(asLO,2)*dasNLO;
  as3_5l = pow(asLO,3)+3.*pow(asLO,2)*dasNLO+3.*pow(asLO,2)*dasNNLO+3.*asLO*pow(dasNLO,2);
  //alphas truncated powers of 4
  as4_4l = pow(asLO,4);
  as4_5l = pow(asLO,4)+4.*pow(asLO,3)*dasNLO;
  //alphas truncated powers of 5
  as5_5l = pow(asLO,5);

  //QCD coupling scale dependence

  //double LQR = log(pow(q0/scales::ren,2));
  double LQR = log(pow(scales::res/scales::ren,2));

  //cout << commute << "  " << LQR << endl;
  
  //if (nloop >= 2) as += -pow(as,2)*beta0/M_PI*LQR;
  //if (nloop >= 3) as += -pow(as,3)*LQR*(beta1 - pow(beta0,2)*LQR)/pi2;
  //if (nloop >= 4) as += -pow(as,4)*LQR*(beta2 - 5./2.*beta0*beta1*LQR + pow(beta0,3)*pow(LQR,2))/(pi2*M_PI);
  //if (nloop >= 5) as += -pow(as,5)*LQR*(beta3 - 3./2.*pow(beta1,2)*LQR - 3.*beta0*beta2*LQR + 13./3.*pow(beta0,2)*beta1*pow(LQR,2) - pow(beta0,4)*pow(LQR,3))/pi4;

  //Runge Kutta evolution
  if (opts.asrgkt)
    rgkt(q, q0, as0);

  //Truncate alphas at the exact order
  asNLO    += -as2_2l * beta0/M_PI*LQR;
  asNNLO   += -as2_3l * beta0/M_PI*LQR;
  asNNLO   += -as3_3l * LQR*(beta1 - pow(beta0,2)*LQR)/pi2;
  asNNNLO  += -as2_4l * beta0/M_PI*LQR;
  asNNNLO  += -as3_4l * LQR*(beta1 - pow(beta0,2)*LQR)/pi2;
  asNNNLO  += -as4_4l * LQR*(beta2 - 5./2.*beta0*beta1*LQR + pow(beta0,3)*pow(LQR,2))/(pi2*M_PI);
  asNNNNLO += -as2_5l * beta0/M_PI*LQR;
  asNNNNLO += -as3_5l * LQR*(beta1 - pow(beta0,2)*LQR)/pi2;
  asNNNNLO += -as4_5l * LQR*(beta2 - 5./2.*beta0*beta1*LQR + pow(beta0,3)*pow(LQR,2))/(pi2*M_PI);
  asNNNNLO += -as5_5l * LQR*(beta3 - 3./2.*pow(beta1,2)*LQR - 3.*beta0*beta2*LQR + 13./3.*pow(beta0,2)*beta1*pow(LQR,2) - pow(beta0,4)*pow(LQR,3))/pi4;
  
//  //Do not truncate alphas
//  asNLO   += -pow(asNLO,2)   * beta0/M_PI*LQR;
//  asNNLO  += -pow(asNNLO,2)  * beta0/M_PI*LQR;
//  asNNLO  += -pow(asNNLO,3)  * LQR*(beta1 - pow(beta0,2)*LQR)/pi2;
//  asNNNLO += -pow(asNNNLO,2) * beta0/M_PI*LQR;
//  asNNNLO += -pow(asNNNLO,3) * LQR*(beta1 - pow(beta0,2)*LQR)/pi2;
//  asNNNLO += -pow(asNNNLO,4) * LQR*(beta2 - 5./2.*beta0*beta1*LQR + pow(beta0,3)*pow(LQR,2))/(pi2*M_PI);

//  //Do not evolve alphas  
//  asNLO   += -pow(as0,2) * beta0/M_PI*LQR;
//  asNNLO  += -pow(as0,2) * beta0/M_PI*LQR;
//  asNNLO  += -pow(as0,3) * LQR*(beta1 - pow(beta0,2)*LQR)/pi2;
//  asNNNLO += -pow(as0,2) * beta0/M_PI*LQR;
//  asNNNLO += -pow(as0,3) * LQR*(beta1 - pow(beta0,2)*LQR)/pi2;
//  asNNNLO += -pow(as0,4) * LQR*(beta2 - 5./2.*beta0*beta1*LQR + pow(beta0,3)*pow(LQR,2))/(pi2*M_PI);
  
  
  //1 over pi normalisation
  as0      /= M_PI;
  asLO     /= M_PI;
  asNLO    /= M_PI;
  asNNLO   /= M_PI;
  asNNNLO  /= M_PI;
  asNNNNLO /= M_PI;
    
  //Recompute truncated powers including scale variations and 1 over pi normalisation
  dasNLO   = asNLO-asLO;
  dasNNLO  = asNNLO-asNLO;
  dasNNNLO = asNNNLO-asNNLO;

  as1_1l = asLO;
  as1_2l = asNLO;
  as1_3l = asNNLO;
  as1_4l = asNNNLO;
  as1_5l = asNNNNLO;
  //alphas truncated powers of 2
  as2_2l = pow(asLO,2);
  as2_3l = pow(asLO,2)+2.*dasNLO*asLO;
  as2_4l = pow(asLO,2)+pow(dasNLO,2)+2.*asLO*dasNLO+2.*asLO*dasNNLO;
  as2_5l = pow(asLO,2)+pow(dasNLO,2)+2.*asLO*dasNLO+2.*asLO*dasNNLO+2.*dasNLO*dasNNLO+2.*asLO*dasNNNLO;
  //alphas truncated powers of 3
  as3_3l = pow(asLO,3);
  as3_4l = pow(asLO,3)+3.*pow(asLO,2)*dasNLO;
  as3_5l = pow(asLO,3)+3.*pow(asLO,2)*dasNLO+3.*pow(asLO,2)*dasNNLO+3.*asLO*pow(dasNLO,2);
  //alphas truncated powers of 4
  as4_4l = pow(asLO,4);
  as4_5l = pow(asLO,4)+4.*pow(asLO,3)*dasNLO;
  //alphas truncated powers of 5
  as5_5l = pow(asLO,5);

  if      (nloop == 0) as = as0;
  else if (nloop == 1) as = asLO;
  else if (nloop == 2) as = asNLO;
  else if (nloop == 3) as = asNNLO;
  else if (nloop == 4) as = asNNNLO;
  else if (nloop == 5) as = asNNNNLO;

  //do not truncate alphas
  if (opts.asrgkt)
    {
      as1_1l = as1_2l = as1_3l = as1_4l = as1_5l = as;
      as2_2l = as2_3l = as2_4l = as2_5l = pow(as,2);
      as3_3l = as3_4l = as3_5l =  pow(as,3);
      as4_4l = as4_5l =  pow(as,4);
      as5_5l =  pow(as,5);
    }
 //cout << "as analytic             " << as << "    " "asL0           " << as0 <<"    " << endl;
  return as;
}

void alphas::fbetabar_cache(double csi, int nloop)
{
  double lk = csi;
  double lk2 = lk * lk;
  double lk3 = lk * lk2;
  double lk4 = lk * lk3;
  double lk5 = lk * lk4;
  double lk6 = lk * lk5;
  double lk7 = lk * lk6;

  double b0 = beta0;
  double b1 = beta1;
  double b2 = beta2;
  double b3 = beta3;
  double b4 = beta4;
  double b5 = beta5;
  double b6 = 0.;
  
  b1csi = 2.*pow(b0,2)*lk;
  b2csi = 5.*b0*b1*lk + 7.*pow(b0,3)*lk2;
  b3csi = 3.*pow(b1,2)*lk + 6.*b0*b2*lk + 29.*pow(b0,2)*b1*lk2 + 32.*pow(b0,4)*lk3;
  b4csi = 7.*b1*b2*lk + 7.*b0*b3*lk + 77./2.*b0*pow(b1,2)*lk2 + 38.*pow(b0,2)*b2*lk2 + 563./3.*pow(b0,3)*b1*lk3 + 517./3.*pow(b0,5)*lk4;
  b5csi = 4.*pow(b2,2)*lk + 8.*b1*b3*lk + 33./2.*pow(b1,3)*lk2 + 8.*b0*b4*lk + 98.*b0*b1*b2*lk2 + 48.*pow(b0,2)*b3*lk2 + 1193./3.*pow(b0,2)*pow(b1,2)*lk3 + 260*pow(b0,3)*b2*lk3 + 2633./2.*pow(b0,4)*b1*lk4 + 3121./3.*pow(b0,6)*lk5;
  b6csi = 9.*b2*b3*lk + 9.*b1*b4*lk + 123./2.*pow(b1,2)*b2*lk2 + 9.*b0*b5*lk + 61.*b0*pow(b2,2)*lk2 + 121.*b0*b1*b3*lk2 + 363.*b0*pow(b1,3)*lk3 + 59.*pow(b0,2)*b4*lk2 + 1071.*pow(b0,2)*b1*b2*lk3 + 346.*pow(b0,3)*b3*lk3 + 11653./3.*pow(b0,3)*pow(b1,2)*lk4 + 1899.*pow(b0,4)*b2*lk4 + 147526./15.*pow(b0,5)*b1*lk5 + 34224./5.*pow(b0,7)*lk6;
  //b7csi = 5*pow(b3,2)*lk + 10*b2*b4*lk + 10*b1*b5*lk + 75*b1*pow(b2,2)*lk2 + 149./2.*pow(b1,2)*b3*lk2 + 121*pow(b1,4)*lk3 + 10*b0*b6*lk + 148*b0*b2*b3*lk2 + 146*b0*b1*b4*lk2 + 2449./3.*b0*pow(b1,2)*b2*lk3 + 615*b0*pow(b1,2)*b2*lk*lk2 + 71*pow(b0,2)*b5*lk2 + 286./3.*pow(b0,2)*pow(b2,2)*lk3 + 610*pow(b0,2)*pow(b2,2)*lk*lk2 + 183*pow(b0,2)*b1*b3*lk3 + 1210*pow(b0,2)*b1*b3*lk*lk2 + 11587./6.*pow(b0,2)*pow(b1,3)*lk4 + 3630*pow(b0,2)*pow(b1,3)*lk*lk3 - 430./3.*pow(b0,3)*b4*lk3 + 590*pow(b0,3)*b4*lk*lk2 + 383./2.*pow(b0,3)*b1*b2*lk4 + 10710*pow(b0,3)*b1*b2*lk*lk3 - 2497./3.*pow(b0,4)*b3*lk4 + 3460*pow(b0,4)*b3*lk*lk3 - 1366*pow(b0,4)*pow(b1,2)*lk5 + 116530./3.*pow(b0,4)*pow(b1,2)*lk*lk4 - 13045./3.*pow(b0,5)*b2*lk5 + 18990*pow(b0,5)*b2*lk*lk4 - 945019./45.*pow(b0,6)*b1*lk6 + 295052./3.*pow(b0,6)*b1*lk*lk5 - 909818./45.*pow(b0,8)*lk7 + 68448*pow(b0,8)*lk*lk6;

  switch (nloop)
    {
    case 1: fbetabarptr = fbetabar1; break;
    case 2: fbetabarptr = fbetabar2; break;
    case 3: fbetabarptr = fbetabar3; break;
    case 4: fbetabarptr = fbetabar4; break;
    case 5: fbetabarptr = fbetabar5; break;
    case 6: fbetabarptr = fbetabar6; break;
    }
  
}

complex <double> alphas::fbetabar(complex <double> a, double csi, int n)
{
  double lk = csi;
  double lk2 = lk * lk;
  double lk3 = lk * lk2;
  double lk4 = lk * lk3;
  double lk5 = lk * lk4;
  double lk6 = lk * lk5;
  double lk7 = lk * lk6;
  complex <double> a2  = pow(a,2);
  complex <double> a3  = pow(a,3);
  complex <double> a4  = pow(a,4);
  complex <double> a5  = pow(a,5);
  complex <double> a6  = pow(a,6);
  complex <double> a7  = pow(a,7);
  complex <double> a8  = pow(a,8);
  complex <double> a9  = pow(a,9);
  /*
  double bt0 = beta0;
  double bt1 = beta1 + 2.*pow(beta0,2)*lk;
  double bt2 = beta2 + 5.*beta1*beta0*lk + 7.*pow(beta0,3)*lk2;
  double bt3 = beta3 + (3.*pow(beta1,2)+6.*beta0*beta2)*lk + 29.*pow(beta0,2)*beta1*lk2 + 32.*pow(beta0,4)*lk3;
  double bt4 = beta4 + (7.*beta1*beta2+7.*beta0*beta3)*lk + (77./2.*beta0*pow(beta1,2)+38.*pow(beta0,2)*beta2)*lk2 + 563./3.*pow(beta0,3)*beta1*lk3 - 517./3.*pow(beta0,5)*lk4;
  double b0 = - bt0;
  double b1 = (n > 1) ? -bt1 : 0.;
  double b2 = (n > 2) ? -bt2 : 0.;
  double b3 = (n > 3) ? -bt3 : 0.;
  double b4 = (n > 4) ? -bt4 : 0.;
  //beta function and derivatives
  complex <double> beta0p = a2 * b0 + a3 * b1 + a4 * b2 + a5 * b3 + a6 * b4;
  complex <double> beta1p = 2. * a * b0 + 3. * a2 * b1 + 4. * a3 * b2 + 5. * a4 * b3 + 6. * a5 * b4;
  complex <double> beta2p = 2. * b0 + 6. * a * b1 + 12. * a2 * b2 + 20. * a3 * b3 + 30. * a4 * b4;
  complex <double> beta3p = 6. * b1 + 24. * a * b2 + 60. * a2 * b3 + 120. * a3 * b4;
  complex <double> beta4p = 24. * b2 + 120. * a * b3 + 360. * a2 * b4;
  return beta0p
    + beta0p * beta1p * lk
    + (beta0p * pow(beta1p,2) + pow(beta0p,2) * beta2p ) * lk2 / 2.
    + (beta0p * pow(beta1p,3) + 4.*pow(beta0p,2) * beta1p * beta2p + pow(beta0p,3) * beta3p ) * lk3 / 6.;
  //+ (beta0p * pow(beta1p,4) + 11.*pow(beta0p,2)*pow(beta1p,2)*beta2p + 4.*pow(beta0p,3)*pow(beta2p,2) + 7.*pow(beta0p,3)*beta1p*beta3p + pow(beta0p,4)*beta4p) * lk4 / 24.;
  */
  
  double b0 = beta0;
  double b1 = beta1;
  double b2 = beta2;
  double b3 = beta3;
  double b4 = beta4;
  double b5 = beta5;
  double b6 = 0.;
  //if (n == 1) return -a2*b0;
  //if (n == 2) return -a2*b0 - a3*b1 + 5.*a4*b0*b1*lk + 10.*a4*pow(b0,3)*lk2;
  //if (n == 3) return -a2*b0 - a3*b1 - a4*b2 + 3.*a5*pow(b1,2)*lk + 6.*a5*b0*b2*lk + 29.*a5*pow(b0,2)*b1*lk2 + 28.*a5*pow(b0,4)*lk3;
  //if (n == 4) return -a2*b0 - a3*b1 - a4*b2 - a5*b3 + 7.*a6*b1*b2*lk + 7.*a6*b0*b3*lk + 77./2.*a6*b0*pow(b1,2)*lk2 + 38.*a6*pow(b0,2)*b2*lk2 + 563./3.*a6*pow(b0,3)*b1*lk3 + 532./3.*a6*pow(b0,5)*lk4;
  //if (n == 5) return -a2*b0 - a3*b1 - a4*b2 - a5*b3 - a6*b4 + 4.*pow(b2,2)*a7*lk + 8.*b1*b3*a7*lk + 33./2.*pow(b1,3)*a7*lk2 + 8.*b0*b4*a7*lk + 98.*b0*b1*b2*a7*lk2 + 48.*pow(b0,2)*b3*a7*lk2 + 1193./3.*pow(b0,2)*pow(b1,2)*a7*lk3 + 260*pow(b0,3)*b2*a7*lk3 + 2633./2.*pow(b0,4)*b1*a7*lk4;

  if (n == 1) return -a2*b0 + 2.*a3*pow(b0,2)*lk;
  if (n == 2) return -a2*b0 - a3*b1 + 5.*a4*b0*b1*lk + 7.*a4*pow(b0,3)*lk2;
  if (n == 3) return -a2*b0 - a3*b1 - a4*b2 + 3.*a5*pow(b1,2)*lk + 6.*a5*b0*b2*lk + 29.*a5*pow(b0,2)*b1*lk2 + 32.*a5*pow(b0,4)*lk3;
  if (n == 4) return -a2*b0 - a3*b1 - a4*b2 - a5*b3 + 7.*a6*b1*b2*lk + 7.*a6*b0*b3*lk + 77./2.*a6*b0*pow(b1,2)*lk2 + 38.*a6*pow(b0,2)*b2*lk2 + 563./3.*a6*pow(b0,3)*b1*lk3 + 517./3.*a6*pow(b0,5)*lk4;
  if (n == 5) return -a2*b0 - a3*b1 - a4*b2 - a5*b3 - a6*b4 + 4.*pow(b2,2)*a7*lk + 8.*b1*b3*a7*lk + 33./2.*pow(b1,3)*a7*lk2 + 8.*b0*b4*a7*lk + 98.*b0*b1*b2*a7*lk2 + 48.*pow(b0,2)*b3*a7*lk2 + 1193./3.*pow(b0,2)*pow(b1,2)*a7*lk3 + 260*pow(b0,3)*b2*a7*lk3 + 2633./2.*pow(b0,4)*b1*a7*lk4 + 3121./3.*pow(b0,6)*a7*lk5;
  if (n == 6) return -a2*b0 - a3*b1 - a4*b2 - a5*b3 - a6*b4 - a7*b5 + 9.*b2*b3*a8*lk + 9.*b1*b4*a8*lk + 123./2.*pow(b1,2)*b2*a8*lk2 + 9.*b0*b5*a8*lk + 61.*b0*pow(b2,2)*a8*lk2 + 121.*b0*b1*b3*a8*lk2 + 363.*b0*pow(b1,3)*a8*lk3 + 59.*pow(b0,2)*b4*a8*lk2 + 1071.*pow(b0,2)*b1*b2*a8*lk3 + 346.*pow(b0,3)*b3*a8*lk3 + 11653./3.*pow(b0,3)*pow(b1,2)*a8*lk4 + 1899.*pow(b0,4)*b2*a8*lk4 + 147526./15.*pow(b0,5)*b1*a8*lk5 + 34224./5.*pow(b0,7)*a8*lk6;
  if (n == 7) return -a2*b0 - a3*b1 - a4*b2 - a5*b3 - a6*b4 - a7*b5 - b6*a8 + 5*pow(b3,2)*a9*lk + 10*b2*b4*a9*lk + 10*b1*b5*a9*lk + 75*b1*pow(b2,2)*a9*lk2 + 149./2.*pow(b1,2)*b3*a9*lk2 + 121*pow(b1,4)*a9*lk3 + 10*b0*b6*a9*lk + 148*b0*b2*b3*a9*lk2 + 146*b0*b1*b4*a9*lk2 + 2449./3.*b0*pow(b1,2)*b2*a9*lk3 + 615*b0*pow(b1,2)*b2*a9*lk*lk2 + 71*pow(b0,2)*b5*a9*lk2 + 286./3.*pow(b0,2)*pow(b2,2)*a9*lk3 + 610*pow(b0,2)*pow(b2,2)*a9*lk*lk2 + 183*pow(b0,2)*b1*b3*a9*lk3 + 1210*pow(b0,2)*b1*b3*a9*lk*lk2 + 11587./6.*pow(b0,2)*pow(b1,3)*a9*lk4 + 3630*pow(b0,2)*pow(b1,3)*a9*lk*lk3 - 430./3.*pow(b0,3)*b4*a9*lk3 + 590*pow(b0,3)*b4*a9*lk*lk2 + 383./2.*pow(b0,3)*b1*b2*a9*lk4 + 10710*pow(b0,3)*b1*b2*a9*lk*lk3 - 2497./3.*pow(b0,4)*b3*a9*lk4 + 3460*pow(b0,4)*b3*a9*lk*lk3 - 1366*pow(b0,4)*pow(b1,2)*a9*lk5 + 116530./3.*pow(b0,4)*pow(b1,2)*a9*lk*lk4 - 13045./3.*pow(b0,5)*b2*a9*lk5 + 18990*pow(b0,5)*b2*a9*lk*lk4 - 945019./45.*pow(b0,6)*b1*a9*lk6 + 295052./3.*pow(b0,6)*b1*a9*lk*lk5 - 909818./45.*pow(b0,8)*a9*lk7 + 68448*pow(b0,8)*a9*lk*lk6;
  
  return 0;
}


//Betabar function as in Eq.(8) of https://arxiv.org/pdf/2202.03380.pdf
complex <double> alphas::betabar(complex <double> q, int nloop, double csi)
{
  double q0 = scales::ren;
  as0 = pdf::alphas(scales::ren);

  //Solve the evolution equation with a fourth-order Runge-Kutta.
  complex <double> lrrat = log(pow(q/q0,2));
  complex <double> dlr   = lrrat / double(nastps);
  complex <double> as    = as0/M_PI;
  complex <double> xk0, xk1, xk2, xk3;

  //if (nloop == 1) fbetabarptr = fbetabar1;
  //else if (nloop == 2) fbetabarptr = fbetabar2;
  //else if (nloop == 3) fbetabarptr = fbetabar3;
  //else if (nloop == 4) fbetabarptr = fbetabar4;
  //else if (nloop == 5) fbetabarptr = fbetabar5;
  
  for (int k1 = 1; k1 <= nastps; k1++)
    {
      xk0 = dlr * fbetabar(as,csi,nloop);
      xk1 = dlr * fbetabar(as+0.5*xk0,csi,nloop);
      xk2 = dlr * fbetabar(as+0.5*xk1,csi,nloop);
      xk3 = dlr * fbetabar(as+xk2,csi,nloop);

      //xk0 = dlr * fbetabarptr(as);
      //xk1 = dlr * fbetabarptr(as+0.5*xk0);
      //xk2 = dlr * fbetabarptr(as+0.5*xk1);
      //xk3 = dlr * fbetabarptr(as+xk2);
      

      as += (xk0+2.*xk1+2.*xk2+xk3)/6.;
    }

  return as;
}

//interface to the VFN alphas of pegasus
double alphas::vfn(double q)
{
  double ASI;
  int NF;  
  double R2 = pow(q,2);
  double M2 = pow(q,2);
  double ASF = 0.;
  if (R2 > asfthr_.m2t_)
    {
      NF = 6;
      double R2T = asfthr_.m2t_;
      ASI = asfthr_.ast_;
      ASF = as_(R2, R2T, asfthr_.ast_, NF);
    }
  else if (M2 > asfthr_.m2b_)
    {
      NF = 5;
      double R2B = asfthr_.m2b_;
      ASI = asfthr_.asb_;
      ASF =  as_(R2, R2B, asfthr_.asb_, NF);
    }
  else if (M2 > asfthr_.m2c_)
    {
      NF = 4;
      double R2C = asfthr_.m2c_;
      ASI = asfthr_.asc_;
      ASF =  as_(R2, R2C, asfthr_.asc_, NF);
    }
  else
    {
      NF = 3;
      double R20 = asinp_.m20_;
      ASI = asinp_.as0_;
      ASF =  as_(R2, R20, asinp_.as0_, NF);
    }
  
  //ASF = LHAPDF::alphasPDF(sqrt(M2)) / (4.* M_PI);
  return ASF*4.*M_PI;
}

// code for Runge-Kutta solutions from Pegasus asrgkt.f extended to five-loop with beta4
//
// ..File: asrgkt.f      (requires a previous call of BETAFCT)
//
//
// ..The running coupling of QCD,  
//
//         AS  =  a_s  =  alpha_s(mu_r^2)/(4 pi),
//
//    obtained by integrating the evolution equation for a fixed number
//    of massless flavours  NF.  Except at leading order (LO),  AS  is 
//    obtained using a fourth-order Runge-Kutta integration. 
//
// ..The initial and final scales  R20  and  R2,  the value  AS0  at
//    R20, and  NF  are passed as function arguments.  The coefficients 
//    of the beta function up to  a_s^5 (N^3LO)  are provided by the 
//    common-block  BETA.  The order of the expansion  NAORD  (defined 
//    as the 'n' in N^nLO) and the number of steps  NASTPS  for the 
//    integration beyond LO are given by the common-block  ASPAR.
//
// =====================================================================

//5 flavour evolution with scheme dependent beta (msbar or cmw)
void alphas::rgkt5(complex <double> q, double q0, double aass0)
{
  complex <double> lrrat = log(pow(q/q0,2));
  complex <double> dlr = lrrat / double(nastps);

  asLO     = aass0/M_PI / (1.+ beta0*aass0/M_PI*lrrat);
  asNLO    = aass0/M_PI;
  asNNLO   = aass0/M_PI;
  asNNNLO  = aass0/M_PI;
  asNNNNLO = aass0/M_PI;

  //Solve the evolution equation with a fourth-order Runge-Kutta.
  complex <double> xk0, xk1, xk2, xk3;  
  for (int k1 = 1; k1 <= nastps; k1++)
    {
      //2-loop
      xk0 = dlr * fbeta1(asNLO);
      xk1 = dlr * fbeta1(asNLO+0.5*xk0);
      xk2 = dlr * fbeta1(asNLO+0.5*xk1);
      xk3 = dlr * fbeta1(asNLO+xk2);
      asNLO += (xk0+2.*xk1+2.*xk2+xk3)/6.;

      //3-loop
      xk0 = dlr * fbeta2(asNNLO);
      xk1 = dlr * fbeta2(asNNLO + 0.5 * xk0);
      xk2 = dlr * fbeta2(asNNLO + 0.5 * xk1);
      xk3 = dlr * fbeta2(asNNLO + xk2);
      asNNLO = asNNLO + (xk0 + 2.* xk1 + 2.* xk2 + xk3)/6.;

      //4-loop
      xk0 = dlr * fbeta3(asNNNLO);
      xk1 = dlr * fbeta3(asNNNLO + 0.5 * xk0);
      xk2 = dlr * fbeta3(asNNNLO + 0.5 * xk1);
      xk3 = dlr * fbeta3(asNNNLO + xk2);
      asNNNLO = asNNNLO + (xk0 + 2.* xk1 + 2.* xk2 + xk3)/6.;

      //5-loop
      xk0 = dlr * fbeta4(asNNNNLO);
      xk1 = dlr * fbeta4(asNNNNLO + 0.5 * xk0);
      xk2 = dlr * fbeta4(asNNNNLO + 0.5 * xk1);
      xk3 = dlr * fbeta4(asNNNNLO + xk2);
      asNNNNLO = asNNNNLO + (xk0 + 2.* xk1 + 2.* xk2 + xk3)/6.;
    }

  asLO     *= M_PI;
  asNLO    *= M_PI;
  asNNLO   *= M_PI;
  asNNNLO  *= M_PI;
  asNNNNLO *= M_PI;
}

//VFN evolution in MSbar
void alphas::rgkt(complex <double> q, double q0, double aass0, int nf)
{
  complex <double> lrrat = log(pow(q/q0,2));
  complex <double> dlr = lrrat / double(nastps);

  asLO     = aass0/M_PI / (1.+ bet0[nf-NFMIN]*aass0/M_PI*lrrat);
  asNLO    = aass0/M_PI;
  asNNLO   = aass0/M_PI;
  asNNNLO  = aass0/M_PI;
  asNNNNLO = aass0/M_PI;

  //Solve the evolution equation with a fourth-order Runge-Kutta.
  complex <double> xk0, xk1, xk2, xk3;  
  for (int k1 = 1; k1 <= nastps; k1++)
    {
      //2-loop
      xk0 = dlr * fbeta1(asNLO,nf);
      xk1 = dlr * fbeta1(asNLO+0.5*xk0,nf);
      xk2 = dlr * fbeta1(asNLO+0.5*xk1,nf);
      xk3 = dlr * fbeta1(asNLO+xk2,nf);
      asNLO += (xk0+2.*xk1+2.*xk2+xk3)/6.;

      //3-loop
      xk0 = dlr * fbeta2(asNNLO,nf);
      xk1 = dlr * fbeta2(asNNLO + 0.5 * xk0,nf);
      xk2 = dlr * fbeta2(asNNLO + 0.5 * xk1,nf);
      xk3 = dlr * fbeta2(asNNLO + xk2,nf);
      asNNLO = asNNLO + (xk0 + 2.* xk1 + 2.* xk2 + xk3)/6.;

      //4-loop
      xk0 = dlr * fbeta3(asNNNLO,nf);
      xk1 = dlr * fbeta3(asNNNLO + 0.5 * xk0,nf);
      xk2 = dlr * fbeta3(asNNNLO + 0.5 * xk1,nf);
      xk3 = dlr * fbeta3(asNNNLO + xk2,nf);
      asNNNLO = asNNNLO + (xk0 + 2.* xk1 + 2.* xk2 + xk3)/6.;

      //5-loop
      xk0 = dlr * fbeta4(asNNNNLO,nf);
      xk1 = dlr * fbeta4(asNNNNLO + 0.5 * xk0,nf);
      xk2 = dlr * fbeta4(asNNNNLO + 0.5 * xk1,nf);
      xk3 = dlr * fbeta4(asNNNNLO + xk2,nf);
      asNNNNLO = asNNNNLO + (xk0 + 2.* xk1 + 2.* xk2 + xk3)/6.;
    }

  asLO     *= M_PI;
  asNLO    *= M_PI;
  asNNLO   *= M_PI;
  asNNNLO  *= M_PI;
  asNNNNLO *= M_PI;
  //cout <<" asLO " << asLO << endl;
}

void alphas::expan(complex <double> q, double q0, double aass0, int nf) //used in pegasus.C to evolve below Q0
{
  complex <double> qlog = log(pow(q/q0,2));

  complex <double> lambda = aass0*bet0[nf-NFMIN]/M_PI*qlog;
  complex <double> oneplambda = 1. + lambda;
  complex <double> log1plambda = log(oneplambda);

  asLO     = aass0/(1.+lambda);
  asNLO    = asLO    - (pow(aass0,2)*bet1[nf-NFMIN]*log1plambda)/(bet0[nf-NFMIN]*pow(oneplambda,2)*M_PI);
  asNNLO   = asNLO   + (pow(aass0,3)*((pow(bet1[nf-NFMIN],2) - bet0[nf-NFMIN]*bet2[nf-NFMIN])*lambda + pow(bet1[nf-NFMIN],2)*(-1.+log1plambda)*log1plambda))/(pow(bet0[nf-NFMIN],2)*pow(oneplambda,3)*pi2);
  asNNNLO  = asNNLO  + (pow(aass0,4)*(-(lambda*(pow(bet1[nf-NFMIN],3)*lambda - 2.*bet0[nf-NFMIN]*bet1[nf-NFMIN]*bet2[nf-NFMIN]*(oneplambda) + pow(bet0[nf-NFMIN],2)*bet3[nf-NFMIN]*(2.+lambda))) + bet1[nf-NFMIN]*log1plambda*(-4.*pow(bet1[nf-NFMIN],2)*lambda + 2.*bet0[nf-NFMIN]*bet2[nf-NFMIN]*(-1.+2.*lambda) + pow(bet1[nf-NFMIN],2)*(5. - 2.*log1plambda)*log1plambda)))/(2.*pow(bet0[nf-NFMIN],3)*pow(oneplambda,4)*pi2*M_PI);
  asNNNNLO = asNNNLO + (pow(aass0,5)*(lambda*(-6.*bet0[nf-NFMIN]*pow(bet1[nf-NFMIN],2)*bet2[nf-NFMIN]*lambda*(4. + lambda) + pow(bet1[nf-NFMIN],4)*lambda*(9. + 2.*lambda) - 2.*pow(bet0[nf-NFMIN],3)*bet4[nf-NFMIN]*(3. + 3.*lambda + pow(lambda,2)) + pow(bet0[nf-NFMIN],2)*(2.*pow(bet2[nf-NFMIN],2)*lambda*(6. + lambda) + bet1[nf-NFMIN]*bet3[nf-NFMIN]*(6. + 9.*lambda + 4.*pow(lambda,2)))) + 6.*bet1[nf-NFMIN]*(bet0[nf-NFMIN]*bet1[nf-NFMIN]*bet2[nf-NFMIN]*(1. - 2.*lambda)*lambda + pow(bet1[nf-NFMIN],3)*(-3. + lambda)*lambda + pow(bet0[nf-NFMIN],2)*bet3[nf-NFMIN]*(-1. + 2.*lambda + pow(lambda,2)))*log1plambda + 9.*(-2.*bet0[nf-NFMIN]*pow(bet1[nf-NFMIN],2)*bet2[nf-NFMIN]*(-1. + lambda) + pow(bet1[nf-NFMIN],4)*(1. + 2.*lambda))*pow(log1plambda,2) - 26.*pow(bet1[nf-NFMIN],4)*pow(log1plambda,3) + 6.*pow(bet1[nf-NFMIN],4)*pow(log1plambda,4)))/(6.*pow(bet0[nf-NFMIN],4)*pow(oneplambda,5)*pi4);
}


void alphas::iter(complex <double> q, double q0, double aass0, int nf) //Eq. A.1 of https://arxiv.org/pdf/1907.02971.pdf
{
  complex <double> qlog = log(pow(q/q0,2));

  complex <double> lambda = aass0*bet0[nf-NFMIN]/M_PI*qlog;
  complex <double> oneplambda = 1. + lambda;
  complex <double> log1plambda = log(oneplambda);

  asLO    = aass0/(1.+lambda);
  asNLO   = asLO   - (pow(aass0,2)*bet1[nf-NFMIN]*log1plambda)/(bet0[nf-NFMIN]*pow(oneplambda,2)*M_PI);
  asNNLO  = asNLO  + (pow(aass0,3)*((pow(bet1[nf-NFMIN],2) - bet0[nf-NFMIN]*bet2[nf-NFMIN])*lambda + pow(bet1[nf-NFMIN],2)*(-1.+log1plambda)*log1plambda))/(pow(bet0[nf-NFMIN],2)*pow(oneplambda,3)*pi2);
  asNNNLO = asNNLO + (pow(aass0,4)*(-(lambda*(pow(bet1[nf-NFMIN],3)*lambda - 2.*bet0[nf-NFMIN]*bet1[nf-NFMIN]*bet2[nf-NFMIN]*(oneplambda) + pow(bet0[nf-NFMIN],2)*bet3[nf-NFMIN]*(2.+lambda)))
  + bet1[nf-NFMIN]*log1plambda*(-4.*pow(bet1[nf-NFMIN],2)*lambda + 2.*bet0[nf-NFMIN]*bet2[nf-NFMIN]*(-1.+2.*lambda)
  + pow(bet1[nf-NFMIN],2)*(5. - 2.*log1plambda)*log1plambda)))/(2.*pow(bet0[nf-NFMIN],3)*pow(oneplambda,4)*pi2*M_PI);
}

void alphas::anal(complex <double> q, double q0, double aass0, int nf) //Analytic running coupling https://arxiv.org/pdf/hep-ph/0211339.pdf
{
  double lambdaqcd = q0/exp(M_PI/(2.*aass0*bet0[nf-NFMIN]));

  asLO = M_PI/bet0[nf-NFMIN]*(1./log(q*q/(lambdaqcd*lambdaqcd))+(lambdaqcd*lambdaqcd)/(lambdaqcd*lambdaqcd-q*q));
}

//exact iterative solution of alphas (see formulas in Appendix C of https://raw.githubusercontent.com/vbertone/apfelxx/master/docs/latex/pdf/TMDEvolution.pdf)
complex <double> alphas::exact(complex <double> q, double q0, double as0, int nloop, int nf)
{
  double eps = 1e-8; //required absolute accuracy

  complex <double> asL[nmax];
  asL[0] = as0/(1.+as0*bet0[nf-NFMIN]/M_PI*log(pow(q/q0,2)));
  asL[1] = asL[0]/(1. + asL[0]*(brnf[nf-NFMIN][nloop-1][1]*log(asL[0]/as0)));
  asL[2] = asL[0]/(1. + asL[0]*(brnf[nf-NFMIN][nloop-1][1]*log(asL[1]/as0) + brnf[nf-NFMIN][nloop-1][2]*(asL[0]-as0)));

  double as0p[nmax];
  as0p[0] = 1.;
  as0p[1] = as0;
  as0p[2] = as0*as0p[1];

  complex <double> asLkp[nmax];
  asLkp[0] = pow(asL[0],2);
  asLkp[1] = asL[1];

  for (int k = 3; k < nmax; k++)
    {
      as0p[k] = as0*as0p[k-1];

      complex <double> sumkp = 0.;
      for (int kp = 1; kp < k; kp++)
	sumkp += brnf[nf-NFMIN][nloop-1][kp+1]*(asLkp[k-kp-1]-as0p[kp]); //asLkp[k-kp-1] = pow(asL[k-kp-1],kp)

      asL[k] = asL[0]/(1. + asL[0]*(brnf[nf-NFMIN][nloop-1][1]*log(asL[k-1]/as0) + sumkp));
      //cout << k << "  " << asL[k] << endl;
      if (abs(asL[k]-asL[k-1]) < eps) return asL[k];

      asLkp[k-1] = asL[k-1];
      for (int kp = k-2; kp > -1; kp--)
	asLkp[kp] = asLkp[kp]*asL[kp];
    }
  return asL[nmax-1];
}

double alphas::lo(double q, double q0, double as0, int nf)
{
  double lrrat = log(pow(q/q0,2));
  double aslo = as0/M_PI / (1.+ bet0[nf-NFMIN]*as0/M_PI*lrrat);
  aslo *= M_PI;
  return aslo;
}
