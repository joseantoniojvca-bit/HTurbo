#include "cevol.h"

#include "settings.h"
#include "mellinint.h"
#include "resint.h"
#include "blim.h"
#include "scales.h"
#include "mesq.h"
#include "ccoeff.h"
#include "resconst.h"
#include "anomalous.h"
#include "alphas.h"
#include "gint.h"
#include "evoltrunc.h"
#include "bstar.h"
#include "btilde.h"
//#include "hell-x.hh"
#include "dytmatrix.h"

using namespace resconst;

int cevol::dim;

complex <double> *cevol::C0qq;

complex <double> *cevol::C1qg;
complex <double> *cevol::C1qq;
complex <double> *cevol::C1qqb;
complex <double> *cevol::C1qqp;
complex <double> *cevol::C1qqbp;

complex <double> *cevol::C2qg;
complex <double> *cevol::C2qq;
complex <double> *cevol::C2qqb;
complex <double> *cevol::C2qqp;
complex <double> *cevol::C2qqbp;

complex <double> *cevol::C3qg;
complex <double> *cevol::C3qq;
complex <double> *cevol::C3qqb;
complex <double> *cevol::C3qqp;
complex <double> *cevol::C3qqbp;

complex <double> *cevol::C4qg;
complex <double> *cevol::C4qq;
complex <double> *cevol::C4qqb;
complex <double> *cevol::C4qqp;
complex <double> *cevol::C4qqbp;

complex <double> *cevol::C0qq_1;
complex <double> *cevol::C1qg_1;
complex <double> *cevol::C1qq_1;
complex <double> *cevol::C1qqb_1;
complex <double> *cevol::C1qqp_1;
complex <double> *cevol::C1qqbp_1;
complex <double> *cevol::C2qg_1;
complex <double> *cevol::C2qq_1;
complex <double> *cevol::C2qqb_1;
complex <double> *cevol::C2qqp_1;
complex <double> *cevol::C2qqbp_1;
complex <double> *cevol::C3qg_1;
complex <double> *cevol::C3qq_1;
complex <double> *cevol::C3qqb_1;
complex <double> *cevol::C3qqp_1;
complex <double> *cevol::C3qqbp_1;
complex <double> *cevol::C4qg_1;
complex <double> *cevol::C4qq_1;
complex <double> *cevol::C4qqb_1;
complex <double> *cevol::C4qqp_1;
complex <double> *cevol::C4qqbp_1;

complex <double> *cevol::C0qq_2;
complex <double> *cevol::C1qg_2;
complex <double> *cevol::C1qq_2;
complex <double> *cevol::C1qqb_2;
complex <double> *cevol::C1qqp_2;
complex <double> *cevol::C1qqbp_2;
complex <double> *cevol::C2qg_2;
complex <double> *cevol::C2qq_2;
complex <double> *cevol::C2qqb_2;
complex <double> *cevol::C2qqp_2;
complex <double> *cevol::C2qqbp_2;
complex <double> *cevol::C3qg_2;
complex <double> *cevol::C3qq_2;
complex <double> *cevol::C3qqb_2;
complex <double> *cevol::C3qqp_2;
complex <double> *cevol::C3qqbp_2;
complex <double> *cevol::C4qg_2;
complex <double> *cevol::C4qq_2;
complex <double> *cevol::C4qqb_2;
complex <double> *cevol::C4qqp_2;
complex <double> *cevol::C4qqbp_2;

//HELLx::sqmatrix<complex<double>> I = HELLx::sqmatrix<complex<double>>(1.,0.,0.,1.);

void cevol::init()
{
  if (opts.mellin1d)
    dim = mellinint::mdim;
  else
    dim = 2*mellinint::mdim;
}

//cevol::update()
//{
//
//}

void cevol::allocate()
{
  if (opts.mellin1d)
    {
      C0qq = new complex <double>[mellinint::mdim*2];

      C1qg = new complex <double>[mellinint::mdim*2];
      C1qq = new complex <double>[mellinint::mdim*2];
      C1qqb = new complex <double>[mellinint::mdim*2];
      C1qqp = new complex <double>[mellinint::mdim*2];
      C1qqbp = new complex <double>[mellinint::mdim*2];

      C2qg = new complex <double>[mellinint::mdim*2];
      C2qq = new complex <double>[mellinint::mdim*2];
      C2qqb = new complex <double>[mellinint::mdim*2];
      C2qqp = new complex <double>[mellinint::mdim*2];
      C2qqbp = new complex <double>[mellinint::mdim*2];

      C3qg = new complex <double>[mellinint::mdim*2];
      C3qq = new complex <double>[mellinint::mdim*2];
      C3qqb = new complex <double>[mellinint::mdim*2];
      C3qqp = new complex <double>[mellinint::mdim*2];
      C3qqbp = new complex <double>[mellinint::mdim*2];

      C4qg = new complex <double>[mellinint::mdim*2];
      C4qq = new complex <double>[mellinint::mdim*2];
      C4qqb = new complex <double>[mellinint::mdim*2];
      C4qqp = new complex <double>[mellinint::mdim*2];
      C4qqbp = new complex <double>[mellinint::mdim*2];
    }
  else
    {
      C0qq_1 = new complex <double>[mellinint::mdim*2];
      C1qg_1 = new complex <double>[mellinint::mdim*2];
      C1qq_1 = new complex <double>[mellinint::mdim*2];
      C1qqb_1 = new complex <double>[mellinint::mdim*2];
      C1qqp_1 = new complex <double>[mellinint::mdim*2];
      C1qqbp_1 = new complex <double>[mellinint::mdim*2];
      C2qg_1 = new complex <double>[mellinint::mdim*2];
      C2qq_1 = new complex <double>[mellinint::mdim*2];
      C2qqb_1 = new complex <double>[mellinint::mdim*2];
      C2qqp_1 = new complex <double>[mellinint::mdim*2];
      C2qqbp_1 = new complex <double>[mellinint::mdim*2];
      C3qg_1 = new complex <double>[mellinint::mdim*2];
      C3qq_1 = new complex <double>[mellinint::mdim*2];
      C3qqb_1 = new complex <double>[mellinint::mdim*2];
      C3qqp_1 = new complex <double>[mellinint::mdim*2];
      C3qqbp_1 = new complex <double>[mellinint::mdim*2];
      C4qg_1 = new complex <double>[mellinint::mdim*2];
      C4qq_1 = new complex <double>[mellinint::mdim*2];
      C4qqb_1 = new complex <double>[mellinint::mdim*2];
      C4qqp_1 = new complex <double>[mellinint::mdim*2];
      C4qqbp_1 = new complex <double>[mellinint::mdim*2];
      C0qq_2 = new complex <double>[mellinint::mdim*2];
      C1qg_2 = new complex <double>[mellinint::mdim*2];
      C1qq_2 = new complex <double>[mellinint::mdim*2];
      C1qqb_2 = new complex <double>[mellinint::mdim*2];
      C1qqp_2 = new complex <double>[mellinint::mdim*2];
      C1qqbp_2 = new complex <double>[mellinint::mdim*2];
      C2qg_2 = new complex <double>[mellinint::mdim*2];
      C2qq_2 = new complex <double>[mellinint::mdim*2];
      C2qqb_2 = new complex <double>[mellinint::mdim*2];
      C2qqp_2 = new complex <double>[mellinint::mdim*2];
      C2qqbp_2 = new complex <double>[mellinint::mdim*2];
      C3qg_2 = new complex <double>[mellinint::mdim*2];
      C3qq_2 = new complex <double>[mellinint::mdim*2];
      C3qqb_2 = new complex <double>[mellinint::mdim*2];
      C3qqp_2 = new complex <double>[mellinint::mdim*2];
      C3qqbp_2 = new complex <double>[mellinint::mdim*2];
      C4qg_2 = new complex <double>[mellinint::mdim*2];
      C4qq_2 = new complex <double>[mellinint::mdim*2];
      C4qqb_2 = new complex <double>[mellinint::mdim*2];
      C4qqp_2 = new complex <double>[mellinint::mdim*2];
      C4qqbp_2 = new complex <double>[mellinint::mdim*2];
    }
}

void cevol::free()
{
  if (opts.mellin1d)
    {
      delete[] C0qq;

      delete[] C1qg;
      delete[] C1qq;
      delete[] C1qqb;
      delete[] C1qqp;
      delete[] C1qqbp;

      delete[] C2qg;
      delete[] C2qq;
      delete[] C2qqb;
      delete[] C2qqp;
      delete[] C2qqbp;

      delete[] C3qg;
      delete[] C3qq;
      delete[] C3qqb;
      delete[] C3qqp;
      delete[] C3qqbp;

      delete[] C4qg;
      delete[] C4qq;
      delete[] C4qqb;
      delete[] C4qqp;
      delete[] C4qqbp;
    }
  else
    {
      delete[] C0qq_1;
      delete[] C1qg_1;
      delete[] C1qq_1;
      delete[] C1qqb_1;
      delete[] C1qqp_1;
      delete[] C1qqbp_1;
      delete[] C2qg_1;
      delete[] C2qq_1;
      delete[] C2qqb_1;
      delete[] C2qqp_1;
      delete[] C2qqbp_1;
      delete[] C3qg_1;
      delete[] C3qq_1;
      delete[] C3qqb_1;
      delete[] C3qqp_1;
      delete[] C3qqbp_1;
      delete[] C4qg_1;
      delete[] C4qq_1;
      delete[] C4qqb_1;
      delete[] C4qqp_1;
      delete[] C4qqbp_1;

      delete[] C0qq_2;
      delete[] C1qg_2;
      delete[] C1qq_2;
      delete[] C1qqb_2;
      delete[] C1qqp_2;
      delete[] C1qqbp_2;
      delete[] C2qg_2;
      delete[] C2qq_2;
      delete[] C2qqb_2;
      delete[] C2qqp_2;
      delete[] C2qqbp_2;
      delete[] C3qg_2;
      delete[] C3qq_2;
      delete[] C3qqb_2;
      delete[] C3qqp_2;
      delete[] C3qqbp_2;
      delete[] C4qg_2;
      delete[] C4qq_2;
      delete[] C4qqb_2;
      delete[] C4qqp_2;
      delete[] C4qqbp_2;
    }
}
void cevol::calc(complex <double> b)
{
  //LL (no evolution)
  if (opts.order_expc == 0)
    return;

  //cout << "cevol::calc " << endl;
  
  //Mellin arrays of evolution matrices
  dytmatrix::sqmatrix<complex<double>> usg[dim*2];
  complex <double> unsp[dim*2];
  complex <double> unsm[dim*2];
  complex <double> unsv[dim*2];
  fill(usg,usg+dim*2, dytmatrix::sqmatrix<complex<double>>(1.,0.,0.,1.));
  fill(unsp,unsp+dim*2, 1.);
  fill(unsm,unsm+dim*2, 1.);
  fill(unsv,unsv+dim*2, 1.);

  int n,nn,sign;
  
  dytmatrix::sqmatrix<complex<double>> I = dytmatrix::sqmatrix<complex<double>>(1.,0.,0.,1.);
  dytmatrix::sqmatrix<complex<double>> N = dytmatrix::sqmatrix<complex<double>>(0.,0.,0.,0.);

  //alphas for the evolution
  //double asi = resint::alpqres*4.*M_PI;
  double asi = resint::alpqres*4.; //alphas/pi
  //double asf = resint::alpqres*4.;

  double blim = blim::expc;  //better use blim::sudakov?
  
  //Set b according to bstar or other prescriptions
  complex <double> bstar;
  if (opts.bprescription == 0 || opts.bprescription == 4 || opts.bstar_expc)
    //bstar = real(b)/sqrt(1.+pow(real(b)/blim,2));
    bstar = bstar::calc(b,blim);
  else
    bstar = b;

  double Q = scales::res;
  double as = resint::aass;
  double as2 = pow(as,2);
  double as3 = pow(as,3);
  
  complex <double> blog;
  if (!opts.modlog)
	  blog = log(pow(Q*bstar/b0,2));   //normal sudakov
	else if (opts.p == 1)
		blog = log(pow(Q*bstar/b0,2) + 1.); //modified sudakov
	else
		blog = 1./btilde::p*log(pow(Q*bstar/b0,2*btilde::p) + 1.); //modified sudakov with exponent p

  //lambda as defined in Eq. (25) of hep-ph/0508068.
  complex <double> xlambda = beta0*as*blog;
  complex <double> log1xlambda = log(1.-xlambda);

  complex <double> logasl;

  //LL (no evolution)
  logasl = 0.;
  
  //NLL (LO evolution)
  if (opts.order_expc >= 1)
    logasl += log1xlambda;

  //NNLL (NLO evolution)
  if (opts.order_expc >= 2)
    logasl += as* beta1/beta0*log1xlambda/(1.-xlambda);

  //NNNLL (NNLO evolution)
  if (opts.order_expc >= 3)
    logasl += as2* ((pow(beta1/beta0,2)-beta2/beta0) *xlambda/pow(1.-xlambda,2)
		    + pow(beta1/beta0,2)             *log1xlambda/pow(1.-xlambda,2)
		    - pow(beta1/beta0,2)             *pow(log1xlambda,2)/(2.*pow(1.-xlambda,2)));

  //NNNNLL (NNNLO evolution)
  if (opts.order_expc >= 4)
    logasl += as3*
      (3.*xlambda*(pow(beta0,2)*beta3*(-2. + xlambda)
		   - 2.*beta0*beta1*beta2*(-1. + xlambda)
		   + pow(beta1,3)*xlambda)
       + 6.*(-(pow(beta1,3)*xlambda)
	     + beta0*beta1*beta2*(1. + xlambda))*log1xlambda
       - 9.*pow(beta1,3)*pow(log1xlambda,2)
       + 2.*pow(beta1,3)*pow(log1xlambda,3)
       )/(6.*pow(beta0,3)*pow(1. - xlambda,3));
  
  //QCD coupling scale dependence
  double LQR = resint::LR-resint::LQ;
  if (opts.order_expc >= 2)
    logasl += as*LQR*beta0*xlambda/(1.-xlambda);

  if (opts.order_expc >= 3)
    logasl += as2*(+LQR*beta1                   *(xlambda-log1xlambda)/pow(1.-xlambda,2)
		   +LQR*beta1                   *xlambda/(1.-xlambda)                      //missing piece
		   +0.5*pow(LQR,2)*pow(beta0,2) *xlambda*(xlambda-2.)/pow(1.-xlambda,2));  //missing piece

  if (opts.order_expc >= 4)
    logasl += as3*
      (LQR*(6.*pow(beta1,2)*pow(log1xlambda,2)
	    + 6.*beta1*log1xlambda*(-3.*beta1 + pow(beta0,2)*LQR)
	    + xlambda*(-6.*pow(beta1,2) + 6.*beta0*beta2*(4. - 3.*xlambda + pow(xlambda,2))
		       - 9.*pow(beta0,2)*beta1*(3. - 3.*xlambda + pow(xlambda,2))*LQR
		       + 2.*pow(beta0,4)*(3. - 3.*xlambda + pow(xlambda,2))*pow(LQR,2))))/(6.*beta0*pow(1. - xlambda,3));
  
  
  //logasl *= -1;
  complex <double> asf = exp(-logasl)*asi;
  //complex <double> asi = exp(logasl)*asf;
  //cout << "asi " << asi*M_PI << " asf " << asf*M_PI << endl;

  for (int i = 0; i < dim*2; i++)
    {
      n    = (i<dim)? i : i-dim;
      sign = (i<dim)? mesq::positive : mesq::negative;
      if (opts.mellin1d)
	{
	  nn = n;
	  int idx = anomalous::index(nn,sign);
	  
	  /*************  All this part can be cached *******************/
	  complex <double> C1gg = 0.;
	  complex <double> C2gg = 0.;
	  complex <double> C3gg = 0.;
	  complex <double> C1gq = 0.;
	  complex <double> C2gq = 0.;
	  complex <double> C3gq = 0.;
	  complex <double> C1qg = ccoeff::C1qg[idx];
	  complex <double> C2qg = ccoeff::C2qg[idx];
	  complex <double> C3qg = ccoeff::C3qg[idx];
	  dytmatrix::sqmatrix<complex<double>> C1sg = dytmatrix::sqmatrix<complex<double>> (C1gg, C1gq, 2.*resconst::NF*C1qg, ccoeff::C1qq[idx]);
	  dytmatrix::sqmatrix<complex<double>> C2sg = dytmatrix::sqmatrix<complex<double>> (C2gg, C2gq, 2.*resconst::NF*C2qg, ccoeff::C2qq[idx]+(resconst::NF-1.)*ccoeff::C2qqp[idx]+ccoeff::C2qqb[idx]+(resconst::NF-1.)*ccoeff::C2qqbp[idx]);
	  dytmatrix::sqmatrix<complex<double>> C3sg = dytmatrix::sqmatrix<complex<double>> (C3gg, C3gq, 2.*resconst::NF*C3qg, ccoeff::C3qq[idx]+(resconst::NF-1.)*ccoeff::C3qqp[idx]+ccoeff::C3qqb[idx]+(resconst::NF-1.)*ccoeff::C3qqbp[idx]);

	  //Compute the Rk matrices with Eq. (2.21) of https://arxiv.org/pdf/hep-ph/0408244.pdf
	  dytmatrix::sqmatrix<complex<double>> R0 = N;
	  dytmatrix::sqmatrix<complex<double>> R1 = C1sg;
	  dytmatrix::sqmatrix<complex<double>> R2 = C2sg-C1sg*C1sg/2.;
	  dytmatrix::sqmatrix<complex<double>> R3 = C3sg-(C1sg*C2sg+C2sg*C1sg)/2.+C1sg*C1sg*C1sg/3.;
	  complex<double> rp = 0.;
	  complex<double> rm = 0.;
	  dytmatrix::sqmatrix<complex<double>> Ep = dytmatrix::sqmatrix<complex<double>>(1.,0.,0.,0.);
	  dytmatrix::sqmatrix<complex<double>> Em = dytmatrix::sqmatrix<complex<double>>(0.,0.,0.,1.);

	  //dytmatrix::sqmatrix<complex<double>> R0 = C1sg;
	  //dytmatrix::sqmatrix<complex<double>> R1 = C2sg-C1sg*C1sg/2.;
	  //dytmatrix::sqmatrix<complex<double>> R2 = C3sg-(C1sg*C2sg+C2sg*C1sg)/2.+C1sg*C1sg*C1sg/3.;
	  //dytmatrix::sqmatrix<complex<double>> R3 = dytmatrix::sqmatrix<complex<double>>(0.,0.,0.,0.);
	  //
	  ////Compute eigenvalues of the R0 matrix (See Eq. (2.27) of https://arxiv.org/pdf/hep-ph/0408244.pdf)
	  //complex<double> htr = R0.trace()/2.;
	  //complex<double> d = R0.det();
	  //complex<double> rp = htr + sqrt(htr*htr - d); 
	  //complex<double> rm = htr - sqrt(htr*htr - d);
	  //
	  ////Compute eigenvector projections (See Eq. (2.28) of https://arxiv.org/pdf/hep-ph/0408244.pdf)
	  //dytmatrix::sqmatrix<complex<double>> Ep = 1./(rp-rm)*(R0-rm*I);
	  //dytmatrix::sqmatrix<complex<double>> Em = 1./(rm-rp)*(R0-rp*I);

	  
	  //Compute Rktilde and Uk matrices (Eq. (2.25) and (2.31) of https://arxiv.org/pdf/hep-ph/0408244.pdf)
	  dytmatrix::sqmatrix<complex<double>> R1TL = R1; 
	  dytmatrix::sqmatrix<complex<double>> U1 = -Em*R1TL*Em/1. - Ep*R1TL*Ep/1. + Ep*R1TL*Em/(rm-rp-1.) + Em*R1TL*Ep/(rp-rm-1.); //Eq. (2.31) of https://arxiv.org/pdf/hep-ph/0408244.pdf
	  dytmatrix::sqmatrix<complex<double>> R2TL = R2+R1*U1;
	  dytmatrix::sqmatrix<complex<double>> U2 = -Em*R2TL*Em/2. - Ep*R2TL*Ep/2. + Ep*R2TL*Em/(rm-rp-2.) + Em*R2TL*Ep/(rp-rm-2.); //Eq. (2.31) of https://arxiv.org/pdf/hep-ph/0408244.pdf
	  dytmatrix::sqmatrix<complex<double>> R3TL = R3+R2*U1+R1*U2;
	  dytmatrix::sqmatrix<complex<double>> U3 = -Em*R3TL*Em/3. - Ep*R3TL*Ep/3. + Ep*R3TL*Em/(rm-rp-3.) + Em*R3TL*Ep/(rp-rm-3.); //Eq. (2.31) of https://arxiv.org/pdf/hep-ph/0408244.pdf
	  /*************  End caching ************************************/

	  //LO evolution matrix
	  //dytmatrix::sqmatrix<complex<double>> ULO = Em*exp(-logasl*rm) + Ep*exp(-logasl*rp); //Eq. (2.29) of https://arxiv.org/pdf/hep-ph/0408244.pdf
	  dytmatrix::sqmatrix<complex<double>> ULO = I;
	    
	  //Eq. (2.24) of https://arxiv.org/pdf/hep-ph/0408244.pdf
	  //LO evolution
	  if (opts.order_expc >= 1)
	    usg[i] = ULO;

	  //NLO evolution
	  if (opts.order_expc >= 2)
	    usg[i] += asf*U1*ULO - asi*ULO*U1;

	  //NNLO evolution
	  if (opts.order_expc >= 3)
	    usg[i] += asf*asf*U2*ULO - asf*asi*U1*ULO*U1 + asi*asi*ULO*(U1*U1-U2);

	  //NNNLO evolution
	  if (opts.order_expc >= 4)
	    usg[i] += asf*asf*asf*U3*ULO - asf*asf*asi*U2*ULO*U1 + asf*asi*asi*U1*ULO*(U1*U1-U2)
	      - asi*asi*asi*ULO*(U1*U1*U1-U1*U2-U2*U1+U3); //typo in Eq. (2.24) of https://arxiv.org/pdf/hep-ph/0408244.pdf, see Eq. (34) of https://doi.org/10.1016/0550-3213(94)90076-0

	  /*************  All this part can be cached *******************/
	  complex <double> C1nsp = ccoeff::C1qq[idx];
	  complex <double> C1nsm = ccoeff::C1qq[idx];
	  complex <double> C1nsv = ccoeff::C1qq[idx];

	  complex <double> C2nsp = ccoeff::C2qq[idx] + ccoeff::C2qqb[idx] - ccoeff::C2qqp[idx] - ccoeff::C2qqbp[idx];
	  complex <double> C2nsm = ccoeff::C2qq[idx] - ccoeff::C2qqb[idx] - ccoeff::C2qqp[idx] + ccoeff::C2qqbp[idx];
	  complex <double> C2nsv = ccoeff::C2qq[idx]+(resconst::NF-1.)*ccoeff::C2qqp[idx]-ccoeff::C2qqb[idx]+(resconst::NF-1.)*ccoeff::C2qqbp[idx];

	  complex <double> C3nsp = ccoeff::C3qq[idx] + ccoeff::C3qqb[idx] - ccoeff::C3qqp[idx] - ccoeff::C3qqbp[idx];
	  complex <double> C3nsm = ccoeff::C3qq[idx] - ccoeff::C3qqb[idx] - ccoeff::C3qqp[idx] + ccoeff::C3qqbp[idx];
	  complex <double> C3nsv = ccoeff::C3qq[idx]+(resconst::NF-1.)*ccoeff::C3qqp[idx]-ccoeff::C3qqb[idx]+(resconst::NF-1.)*ccoeff::C3qqbp[idx];
  
	  complex <double> R0NSp = 0.;
	  complex <double> R0NSm = 0.;
	  complex <double> R0NSv = 0.;
	  complex <double> R1NSp = C1nsp;
	  complex <double> R1NSm = C1nsm;
	  complex <double> R1NSv = C1nsv;
	  complex <double> R2NSp = C2nsp-C1nsp*C1nsp/2.;
	  complex <double> R2NSm = C2nsm-C1nsm*C1nsm/2.;
	  complex <double> R2NSv = C2nsv-C1nsv*C1nsv/2.;
	  complex <double> R3NSp = C3nsp-C1nsp*C2nsp+C1nsp*C1nsp*C1nsp/3.;
	  complex <double> R3NSm = C3nsm-C1nsm*C2nsm+C1nsm*C1nsm*C1nsm/3.;
	  complex <double> R3NSv = C3nsv-C1nsv*C2nsv+C1nsv*C1nsv*C1nsv/3.;

	  //complex <double> R0NSp = C1nsp;
	  //complex <double> R0NSm = C1nsm;
	  //complex <double> R0NSv = C1nsv;
	  //complex <double> R1NSp = C2nsp-C1nsp*C1nsp/2.;
	  //complex <double> R1NSm = C2nsm-C1nsm*C1nsm/2.;
	  //complex <double> R1NSv = C2nsv-C1nsv*C1nsv/2.;
	  //complex <double> R2NSp = C3nsp-C1nsp*C2nsp+C1nsp*C1nsp*C1nsp/3.;
	  //complex <double> R2NSm = C3nsm-C1nsm*C2nsm+C1nsm*C1nsm*C1nsm/3.;
	  //complex <double> R2NSv = C3nsv-C1nsv*C2nsv+C1nsv*C1nsv*C1nsv/3.;
	  //complex <double> R3NSp = 0.;
	  //complex <double> R3NSm = 0.;
	  //complex <double> R3NSv = 0.;
	  
	  complex <double> U1NSp = -R1NSp;
	  complex <double> U1NSm = -R1NSm;
	  complex <double> U1NSv = -R1NSv;

	  complex<double> R2NSTLp = R2NSp+R1NSp*U1NSp;
	  complex<double> R2NSTLm = R2NSm+R1NSm*U1NSm;
	  complex<double> R2NSTLv = R2NSv+R1NSv*U1NSv;
      
	  complex <double> U2NSp = -(R2NSTLp)/2.;
	  complex <double> U2NSm = -(R2NSTLm)/2.;
	  complex <double> U2NSv = -(R2NSTLv)/2.;

	  complex<double> R3NSTLp = R3NSp+R2NSp*U1NSp+R1NSp*U2NSp;
	  complex<double> R3NSTLm = R3NSm+R2NSm*U1NSm+R1NSm*U2NSm;
	  complex<double> R3NSTLv = R3NSv+R2NSv*U1NSv+R1NSv*U2NSv;

	  complex <double> U3NSp = -(R3NSTLp)/3.;
	  complex <double> U3NSm = -(R3NSTLm)/3.;
	  complex <double> U3NSv = -(R3NSTLv)/3.;
	  /*************  End caching ************************************/

	  complex<double> UNSLO = exp(-R0NSp*logasl);
	  if (opts.order_expc >= 1)
	    {
	      unsp[i] = UNSLO;
	      unsm[i] = UNSLO;
	      unsv[i] = UNSLO;
	    }
      
	  //NLO evolution
	  if (opts.order_expc >= 2)
	    {
	      unsp[i] += (asf-asi)*U1NSp*UNSLO;
	      unsm[i] += (asf-asi)*U1NSm*UNSLO;
	      unsv[i] += (asf-asi)*U1NSv*UNSLO;
	    }

	  //NNLO evolution
	  if (opts.order_expc >= 3)
	    {
	      unsp[i] += asf*asf*U2NSp*UNSLO - asf*asi*U1NSp*U1NSp*UNSLO + asi*asi*(U1NSp*U1NSp-U2NSp)*UNSLO;
	      unsm[i] += asf*asf*U2NSm*UNSLO - asf*asi*U1NSm*U1NSm*UNSLO + asi*asi*(U1NSm*U1NSm-U2NSm)*UNSLO;
	      unsv[i] += asf*asf*U2NSv*UNSLO - asf*asi*U1NSv*U1NSv*UNSLO + asi*asi*(U1NSv*U1NSv-U2NSv)*UNSLO;
	    }

	  //NNNLO evolution
	  if (opts.order_expc >= 4)
	    {
	      unsp[i] += asf*asf*asf*U3NSp*UNSLO - asf*asf*asi*U2NSp*U1NSp*UNSLO + asf*asi*asi*U1NSp*UNSLO*(U1NSp*U1NSp-U2NSp) - asi*asi*asi*UNSLO*(U1NSp*U1NSp*U1NSp-2.*U1NSp*U2NSp+U3NSp);
	      unsm[i] += asf*asf*asf*U3NSm*UNSLO - asf*asf*asi*U2NSm*U1NSm*UNSLO + asf*asi*asi*U1NSm*UNSLO*(U1NSm*U1NSm-U2NSm) - asi*asi*asi*UNSLO*(U1NSm*U1NSm*U1NSm-2.*U1NSm*U2NSm+U3NSm);
	      unsv[i] += asf*asf*asf*U3NSv*UNSLO - asf*asf*asi*U2NSv*U1NSv*UNSLO + asf*asi*asi*U1NSv*UNSLO*(U1NSv*U1NSv-U2NSv) - asi*asi*asi*UNSLO*(U1NSv*U1NSv*U1NSv-2.*U1NSv*U2NSv+U3NSv);
	    }
	}
    }

  for (int i = 0; i < dim*2; i++)
    {
      n    = (i<dim)? i : i-dim;
      sign = (i<dim)? mesq::positive : mesq::negative;
      if (opts.mellin1d)
	{
	  nn = n;
	  int idx = anomalous::index(nn,sign);
	  complex <double> C1gg = 0.;
	  complex <double> C2gg = 0.;
	  complex <double> C3gg = 0.;
	  complex <double> C1gq = 0.;
	  complex <double> C2gq = 0.;
	  complex <double> C3gq = 0.;
	  dytmatrix::sqmatrix<complex<double>> C1sg = dytmatrix::sqmatrix<complex<double>> (C1gg, C1gq, 2.*resconst::NF*ccoeff::C1qg[idx], ccoeff::C1qq[idx]);
	  dytmatrix::sqmatrix<complex<double>> C2sg = dytmatrix::sqmatrix<complex<double>> (C2gg, C2gq, 2.*resconst::NF*ccoeff::C2qg[idx], ccoeff::C2qq[idx]+(resconst::NF-1.)*ccoeff::C2qqp[idx]+ccoeff::C2qqb[idx]+(resconst::NF-1.)*ccoeff::C2qqbp[idx]);
	  dytmatrix::sqmatrix<complex<double>> C3sg = dytmatrix::sqmatrix<complex<double>> (C3gg, C3gq, 2.*resconst::NF*ccoeff::C3qg[idx], ccoeff::C3qq[idx]+(resconst::NF-1.)*ccoeff::C3qqp[idx]+ccoeff::C3qqb[idx]+(resconst::NF-1.)*ccoeff::C3qqbp[idx]);

	  C1sg = C1sg * usg[i];
	  C2sg = C2sg * usg[i];
	  C3sg = C3sg * usg[i];

	  complex <double> C1nsp = ccoeff::C1qq[idx];
	  complex <double> C1nsm = ccoeff::C1qq[idx];
	  complex <double> C1nsv = ccoeff::C1qq[idx];

	  complex <double> C2nsp = ccoeff::C2qq[idx] + ccoeff::C2qqb[idx] - ccoeff::C2qqp[idx] - ccoeff::C2qqbp[idx];
	  complex <double> C2nsm = ccoeff::C2qq[idx] - ccoeff::C2qqb[idx] - ccoeff::C2qqp[idx] + ccoeff::C2qqbp[idx];
	  complex <double> C2nsv = ccoeff::C2qq[idx]+(resconst::NF-1.)*ccoeff::C2qqp[idx]-ccoeff::C2qqb[idx]+(resconst::NF-1.)*ccoeff::C2qqbp[idx];

	  complex <double> C3nsp = ccoeff::C3qq[idx] + ccoeff::C3qqb[idx] - ccoeff::C3qqp[idx] - ccoeff::C3qqbp[idx];
	  complex <double> C3nsm = ccoeff::C3qq[idx] - ccoeff::C3qqb[idx] - ccoeff::C3qqp[idx] + ccoeff::C3qqbp[idx];
	  complex <double> C3nsv = ccoeff::C3qq[idx]+(resconst::NF-1.)*ccoeff::C3qqp[idx]-ccoeff::C3qqb[idx]+(resconst::NF-1.)*ccoeff::C3qqbp[idx];

	  C1nsp = C1nsp * unsp[i];
	  C1nsm = C1nsm * unsm[i];
	  C1nsv = C1nsv * unsv[i];

	  C2nsp = C2nsp * unsp[i];
	  C2nsm = C2nsm * unsm[i];
	  C2nsv = C2nsv * unsv[i];

	  C3nsp = C3nsp * unsp[i];
	  C3nsm = C3nsm * unsm[i];
	  C3nsv = C3nsv * unsv[i];

	  //Change of basis
	  double nf = double(resconst::NF);
	  complex <double> csg_qq,csg_qg,csg_gq,csg_gg,cns_p,cns_m,cns_v;
	  complex <double> cqqS_p_cqqbS,cqqS_m_cqqbS,cqqS,cqqbS,cqqV,cqqbV;
	  
	  csg_qq = C1sg.qq();
	  csg_qg = C1sg.qg();
	  csg_gq = C1sg.gq();
	  csg_gg = C1sg.gg();
	  cns_p = C1nsp;
	  cns_m = C1nsm;
	  cns_v = C1nsv;

	  cqqS_p_cqqbS = (csg_qq-cns_p)/nf;
	  cqqS_m_cqqbS = (cns_v-cns_m)/nf;
	  cqqS         = (cqqS_p_cqqbS+cqqS_m_cqqbS)/2.;
	  cqqbS        = (cqqS_p_cqqbS-cqqS_m_cqqbS)/2.;
	  cqqV         = (cns_p+cns_m)/2.;
	  cqqbV        = (cns_p-cns_m)/2.;

	  C1qq[idx]   = cqqV +cqqS;
	  C1qqb[idx]  = cqqbV+cqqbS;
	  C1qqp[idx]  = cqqS;
	  C1qqbp[idx] = cqqbS;
	  C1qg[idx]   = csg_qg/2./nf;
	  //C1gq[idx]   = csg_gq;
	  //C1gg[idx]   = csg_gg;

	  csg_qq = C2sg.qq();
	  csg_qg = C2sg.qg();
	  csg_gq = C2sg.gq();
	  csg_gg = C2sg.gg();
	  cns_p = C2nsp;
	  cns_m = C2nsm;
	  cns_v = C2nsv;

	  cqqS_p_cqqbS = (csg_qq-cns_p)/nf;
	  cqqS_m_cqqbS = (cns_v-cns_m)/nf;
	  cqqS         = (cqqS_p_cqqbS+cqqS_m_cqqbS)/2.;
	  cqqbS        = (cqqS_p_cqqbS-cqqS_m_cqqbS)/2.;
	  cqqV         = (cns_p+cns_m)/2.;
	  cqqbV        = (cns_p-cns_m)/2.;

	  C2qq[idx]   = cqqV +cqqS;
	  C2qqb[idx]  = cqqbV+cqqbS;
	  C2qqp[idx]  = cqqS;
	  C2qqbp[idx] = cqqbS;
	  C2qg[idx]   = csg_qg/2./nf;
	  //C2gq[idx]   = csg_gq;
	  //C2gg[idx]   = csg_gg;

	  //cout << "cevol " << C1qq[idx]	 << endl;
	}
	  
    }  
}

void cevol::css(complex <double> b)
{
  double blim = blim::expc;
  complex <double> bstar;
  if (opts.bprescription == 0 || opts.bprescription == 4 || opts.bstar_expc)
    //bstar = real(b)/sqrt(1.+pow(real(b)/blim,2));
    bstar = bstar::calc(b,blim);
  else
    bstar = b;

  double Q = scales::res;
  complex <double> mubstar = resconst::b0/bstar;
  complex <double> mubstartilde;
  if (!opts.modlog)
	  mubstartilde = mubstar;
	else if (opts.p == 1)
		mubstartilde = mubstar * Q / sqrt((pow(mubstar,2) + pow(Q,2)));
	//mubstartilde = mubstar / sqrt(1.+pow(mubstar/Q,2));
	else
		mubstartilde = mubstar * Q / pow((pow(mubstar,2*opts.p) + pow(Q,2*opts.p)),0.5/opts.p);

  //iterative alphas
  complex <double> asmub = alphas::calc(mubstartilde, opts.order_expc+1);

  complex <double> as1 = asmub/resint::aass;
  complex <double> as2 = pow(as1,2);
  complex <double> as3 = pow(as1,3);
  complex <double> as4 = pow(as1,4);
  
  int n,nn,sign;
  for (int i = 0; i < dim*2; i++)
    {
      n    = (i<dim)? i : i-dim;
      sign = (i<dim)? mesq::positive : mesq::negative;
      if (opts.mellin1d)
	{
	  nn = n;
	  int idx = anomalous::index(nn,sign);

	  C0qq[idx]   = 1.;
	  
	  C1qq[idx]   = as1*ccoeff::C1qq[idx];
	  C1qqb[idx]  = as1*ccoeff::C1qqb[idx];
	  C1qqp[idx]  = as1*ccoeff::C1qqp[idx];
	  C1qqbp[idx] = as1*ccoeff::C1qqbp[idx];
	  C1qg[idx]   = as1*ccoeff::C1qg[idx];

	  C2qq[idx]   = as2*ccoeff::C2qq[idx];
	  C2qqb[idx]  = as2*ccoeff::C2qqb[idx];
	  C2qqp[idx]  = as2*ccoeff::C2qqp[idx];
	  C2qqbp[idx] = as2*ccoeff::C2qqbp[idx];
	  C2qg[idx]   = as2*ccoeff::C2qg[idx];

	  C3qq[idx]   = as3*ccoeff::C3qq[idx];
	  C3qqb[idx]  = as3*ccoeff::C3qqb[idx];
	  C3qqp[idx]  = as3*ccoeff::C3qqp[idx];
	  C3qqbp[idx] = as3*ccoeff::C3qqbp[idx];
	  C3qg[idx]   = as3*ccoeff::C3qg[idx];

	  C4qq[idx]   = as4*ccoeff::C4qq[idx];
	  C4qqb[idx]  = as4*ccoeff::C4qqb[idx];
	  C4qqp[idx]  = as4*ccoeff::C4qqp[idx];
	  C4qqbp[idx] = as4*ccoeff::C4qqbp[idx];
	  C4qg[idx]   = as4*ccoeff::C4qg[idx];
	}
      else
	{
	  if (n < mellinint::mdim)
	    nn = n;
	  else
	    nn = n-mellinint::mdim;
	  int idx = anomalous::index(nn,sign);

	  //C0qq_1[idx] = 1.;
	  C1qq_1[idx]   = as1*ccoeff::C1qq_1[idx];
	  C1qqb_1[idx]  = as1*ccoeff::C1qqb_1[idx];
	  C1qqp_1[idx]  = as1*ccoeff::C1qqp_1[idx];
	  C1qqbp_1[idx] = as1*ccoeff::C1qqbp_1[idx];
	  C1qg_1[idx]   = as1*ccoeff::C1qg_1[idx];
	  C2qq_1[idx]   = as2*ccoeff::C2qq_1[idx];
	  C2qqb_1[idx]  = as2*ccoeff::C2qqb_1[idx];
	  C2qqp_1[idx]  = as2*ccoeff::C2qqp_1[idx];
	  C2qqbp_1[idx] = as2*ccoeff::C2qqbp_1[idx];
	  C2qg_1[idx]   = as2*ccoeff::C2qg_1[idx];
	  C3qq_1[idx]   = as3*ccoeff::C3qq_1[idx];
	  C3qqb_1[idx]  = as3*ccoeff::C3qqb_1[idx];
	  C3qqp_1[idx]  = as3*ccoeff::C3qqp_1[idx];
	  C3qqbp_1[idx] = as3*ccoeff::C3qqbp_1[idx];
	  C3qg_1[idx]   = as3*ccoeff::C3qg_1[idx];
	  C4qq_1[idx]   = as4*ccoeff::C4qq_1[idx];
	  C4qqb_1[idx]  = as4*ccoeff::C4qqb_1[idx];
	  C4qqp_1[idx]  = as4*ccoeff::C4qqp_1[idx];
	  C4qqbp_1[idx] = as4*ccoeff::C4qqbp_1[idx];
	  C4qg_1[idx]   = as4*ccoeff::C4qg_1[idx];

	  //C0qq_2[idx] = 1.;
	  C1qq_2[idx]   = as1*ccoeff::C1qq_2[idx];
	  C1qqb_2[idx]  = as1*ccoeff::C1qqb_2[idx];
	  C1qqp_2[idx]  = as1*ccoeff::C1qqp_2[idx];
	  C1qqbp_2[idx] = as1*ccoeff::C1qqbp_2[idx];
	  C1qg_2[idx]   = as1*ccoeff::C1qg_2[idx];
	  C2qq_2[idx]   = as2*ccoeff::C2qq_2[idx];
	  C2qqb_2[idx]  = as2*ccoeff::C2qqb_2[idx];
	  C2qqp_2[idx]  = as2*ccoeff::C2qqp_2[idx];
	  C2qqbp_2[idx] = as2*ccoeff::C2qqbp_2[idx];
	  C2qg_2[idx]   = as2*ccoeff::C2qg_2[idx];
	  C3qq_2[idx]   = as3*ccoeff::C3qq_2[idx];
	  C3qqb_2[idx]  = as3*ccoeff::C3qqb_2[idx];
	  C3qqp_2[idx]  = as3*ccoeff::C3qqp_2[idx];
	  C3qqbp_2[idx] = as3*ccoeff::C3qqbp_2[idx];
	  C3qg_2[idx]   = as3*ccoeff::C3qg_2[idx];
	  C4qq_2[idx]   = as4*ccoeff::C4qq_2[idx];
	  C4qqb_2[idx]  = as4*ccoeff::C4qqb_2[idx];
	  C4qqp_2[idx]  = as4*ccoeff::C4qqp_2[idx];
	  C4qqbp_2[idx] = as4*ccoeff::C4qqbp_2[idx];
	  C4qg_2[idx]   = as4*ccoeff::C4qg_2[idx];
	}
    }
}

void cevol::ctilde(complex <double> b)
{
  if (opts.order_expc == 0)
    return;

  double blim = blim::expc;
  
  //Set b according to bstar or other prescriptions
  complex <double> bstar;
  if (opts.bprescription == 0 || opts.bprescription == 4 || opts.bstar_expc)
    //bstar = real(b)/sqrt(1.+pow(real(b)/blim,2));
    bstar = bstar::calc(b,blim);
  else
    bstar = b;

  double Q = scales::res;
  double aass = resint::aass;
  double aass2 = pow(aass,2);
  double aass3 = pow(aass,3);
  
  complex <double> blog;
  if (!opts.modlog)
	  blog = log(pow(Q*bstar/b0,2));   //normal sudakov
	  else if (opts.p == 1)
		  blog = log(pow(Q*bstar/b0,2) + 1.); //modified sudakov
		  else
			  blog = 1./btilde::p*log(pow(Q*bstar/b0,2*btilde::p) + 1.); //modified sudakov with exponent p

  //lambda as defined in Eq. (25) of hep-ph/0508068.
  complex <double> xlambda = beta0*aass*blog;
  complex <double> log1xlambda = log(1.-xlambda);

  complex <double> logasl;

  //LL (no evolution)
  logasl = 0.;
  
  //NLL (LO evolution)
  if (opts.order_expc >= 1)
    logasl += log1xlambda;

  //NNLL (NLO evolution)
  if (opts.order_expc >= 2)
    logasl += aass* beta1/beta0*log1xlambda/(1.-xlambda);

  //NNNLL (NNLO evolution)
  if (opts.order_expc >= 3)
    logasl += aass2* ((pow(beta1/beta0,2)-beta2/beta0) *xlambda/pow(1.-xlambda,2)
		      + pow(beta1/beta0,2)             *log1xlambda/pow(1.-xlambda,2)
		      - pow(beta1/beta0,2)             *pow(log1xlambda,2)/(2.*pow(1.-xlambda,2)));

  //NNNNLL (NNNLO evolution)
  if (opts.order_expc >= 4)
    logasl += aass3*
      (3.*xlambda*(pow(beta0,2)*beta3*(-2. + xlambda)
		   - 2.*beta0*beta1*beta2*(-1. + xlambda)
		   + pow(beta1,3)*xlambda)
       + 6.*(-(pow(beta1,3)*xlambda)
	     + beta0*beta1*beta2*(1. + xlambda))*log1xlambda
       - 9.*pow(beta1,3)*pow(log1xlambda,2)
       + 2.*pow(beta1,3)*pow(log1xlambda,3)
       )/(6.*pow(beta0,3)*pow(1. - xlambda,3));
  
  double LQR = resint::LR-resint::LQ;
  //QCD coupling scale dependence
  if (opts.order_expc >= 2)
    logasl += aass*LQR*beta0*xlambda/(1.-xlambda);

  //from c3new.txt
  if (opts.order_expc >= 3)
    logasl += aass2*(+LQR*beta1                   *(xlambda-log1xlambda)/pow(1.-xlambda,2)
		     +LQR*beta1                   *xlambda/(1.-xlambda)                      //missing piece
    		     +0.5*pow(LQR,2)*pow(beta0,2) *xlambda*(xlambda-2.)/pow(1.-xlambda,2));  //missing piece

  if (opts.order_expc >= 4)
    logasl += aass3*
      (LQR*(6.*pow(beta1,2)*pow(log1xlambda,2)
	    + 6.*beta1*log1xlambda*(-3.*beta1 + pow(beta0,2)*LQR)
	    + xlambda*(-6.*pow(beta1,2) + 6.*beta0*beta2*(4. - 3.*xlambda + pow(xlambda,2))
		       - 9.*pow(beta0,2)*beta1*(3. - 3.*xlambda + pow(xlambda,2))*LQR
		       + 2.*pow(beta0,4)*(3. - 3.*xlambda + pow(xlambda,2))*pow(LQR,2))))/(6.*beta0*pow(1. - xlambda,3));
  
  if (opts.numexpc)
    logasl = gint::logasl_expc;
  
  complex <double> aexp = exp(-logasl); //aexp is approximately alphas(Q^2)/alphas(mub^2)
  complex <double> aexp2 = pow(aexp,2);

  //cout << "cevol aexp " << aexp << endl;
  
  //Multiply lamB, lamC, lamD by the corresponding power of alphas to avoid overflow/uunderflow problems in the exponential
  complex <double> lamB = aass*xlambda/(1.-xlambda);
  complex <double> lamC = aass2*xlambda*(xlambda-2.)/pow(1.-xlambda,2);
  complex <double> lamD = aass2*log1xlambda/pow(1.-xlambda,2);
  complex <double> lamE = aass3*xlambda/pow(1.-xlambda,3);
  complex <double> lamF = aass3*xlambda*(3.-3.*xlambda+pow(xlambda,2))/pow(1.-xlambda,3);
  complex <double> lamG = aass3*log1xlambda/pow(1.-xlambda,3);
  complex <double> lamH = aass3*pow(log1xlambda,2)/pow(1.-xlambda,3);
  
  complex <double> aexpB = exp(lamB);
  complex <double> aexpC = exp(lamC);
  complex <double> aexpD = exp(lamD);
  complex <double> aexpE = exp(lamE);
  complex <double> aexpF = exp(lamF);
  complex <double> aexpG = exp(lamG);
  complex <double> aexpH = exp(lamH);


  int n,nn,sign;
  for (int i = 0; i < dim*2; i++)
    {
      n    = (i<dim)? i : i-dim;
      sign = (i<dim)? mesq::positive : mesq::negative;
      if (opts.mellin1d)
	{
	  nn = n;
	  int idx = anomalous::index(nn,sign);

	  /*
	  complex <double> aexpqq  = 1.;
	  complex <double> aexpqg  = 1.;
	  complex <double> aexpqqb = 1.;
	  complex <double> aexpqqp = 1.;
	  complex <double> aexpqqbp = 1.;

	  if (opts.order_expc == 2)
	    {
	      aexpqq   = pow(aexpB,ccoeff::C1qq[idx]);
	      aexpqg   = pow(aexpB,ccoeff::C2qg[idx]/ccoeff::C1qg[idx]);
	      aexpqqb  = pow(aexpB,ccoeff::C3qqb[idx]/ccoeff::C2qqb[idx]);
	      aexpqqp  = pow(aexpB,ccoeff::C3qqp[idx]/ccoeff::C2qqp[idx]);
	      aexpqqbp = pow(aexpB,ccoeff::C3qqbp[idx]/ccoeff::C2qqbp[idx]);
	    }

	  aexpqg   *= aexp;
	  aexpqqb  *= aexp2;
	  aexpqqp  *= aexp2;
	  aexpqqbp *= aexp2;

	  C0qq[idx]   = aexpqq;
	  
	  C1qq[idx]   = aexpqq   * ccoeff::C1qq[idx];    
	  C1qg[idx]   = aexpqg   * ccoeff::C1qg[idx];    
	  C1qqb[idx]  = aexpqqb  * ccoeff::C1qqb[idx];   
	  C1qqp[idx]  = aexpqqp  * ccoeff::C1qqp[idx];   
	  C1qqbp[idx] = aexpqqbp * ccoeff::C1qqbp[idx];  

	  C2qq[idx]   = aexpqq   * ccoeff::C2qq[idx];    
	  C2qg[idx]   = aexpqg   * ccoeff::C2qg[idx];    
	  C2qqb[idx]  = aexpqqb  * ccoeff::C2qqb[idx];   
	  C2qqp[idx]  = aexpqqp  * ccoeff::C2qqp[idx];   
	  C2qqbp[idx] = aexpqqbp * ccoeff::C2qqbp[idx];  

	  C3qq[idx]   = aexpqq   * ccoeff::C3qq[idx];    
	  C3qg[idx]   = aexpqg   * ccoeff::C3qg[idx];    
	  C3qqb[idx]  = aexpqqb  * ccoeff::C3qqb[idx];   
	  C3qqp[idx]  = aexpqqp  * ccoeff::C3qqp[idx];   
	  C3qqbp[idx] = aexpqqbp * ccoeff::C3qqbp[idx];
	  */
	  
	  //change of basis to SG NS representation
	  complex <double> C1sgqq = ccoeff::C1qq[idx];
	  complex <double> C1sgqg = 2.*resconst::NF*ccoeff::C1qg[idx];

	  complex <double> C2sgqq = ccoeff::C2qq[idx]+(resconst::NF-1.)*ccoeff::C2qqp[idx]+ccoeff::C2qqb[idx]+(resconst::NF-1.)*ccoeff::C2qqbp[idx];
	  complex <double> C2sgqg = 2.*resconst::NF*ccoeff::C2qg[idx];

	  complex <double> C3sgqq = ccoeff::C3qq[idx]+(resconst::NF-1.)*ccoeff::C3qqp[idx]+ccoeff::C3qqb[idx]+(resconst::NF-1.)*ccoeff::C3qqbp[idx];
	  complex <double> C3sgqg = 2.*resconst::NF*ccoeff::C3qg[idx];

	  complex <double> C1nsp = ccoeff::C1qq[idx];
	  complex <double> C1nsm = ccoeff::C1qq[idx];
	  complex <double> C1nsv = ccoeff::C1qq[idx];

	  complex <double> C2nsp = ccoeff::C2qq[idx] + ccoeff::C2qqb[idx] - ccoeff::C2qqp[idx] - ccoeff::C2qqbp[idx];
	  complex <double> C2nsm = ccoeff::C2qq[idx] - ccoeff::C2qqb[idx] - ccoeff::C2qqp[idx] + ccoeff::C2qqbp[idx];
	  complex <double> C2nsv = ccoeff::C2qq[idx]+(resconst::NF-1.)*ccoeff::C2qqp[idx]-ccoeff::C2qqb[idx]+(resconst::NF-1.)*ccoeff::C2qqbp[idx];

	  complex <double> C3nsp = ccoeff::C3qq[idx] + ccoeff::C3qqb[idx] - ccoeff::C3qqp[idx] - ccoeff::C3qqbp[idx];
	  complex <double> C3nsm = ccoeff::C3qq[idx] - ccoeff::C3qqb[idx] - ccoeff::C3qqp[idx] + ccoeff::C3qqbp[idx];
	  complex <double> C3nsv = ccoeff::C3qq[idx]+(resconst::NF-1.)*ccoeff::C3qqp[idx]-ccoeff::C3qqb[idx]+(resconst::NF-1.)*ccoeff::C3qqbp[idx];

	  //Ctilde = C*V;
	  //complex <double> CTL1sgqq = C1sgqq + evoltrunc::U1[i].qq();
	  //complex <double> CTL1sgqg = C1sgqg + evoltrunc::U1[i].qg();
	  //complex <double> CTL2sgqq = C2sgqq + evoltrunc::U2[i].qq() + C1sgqq*evoltrunc::U1[i].qq() + C1sgqg*evoltrunc::U1[i].gq();
	  //complex <double> CTL2sgqg = C2sgqg + evoltrunc::U2[i].qg() + C1sgqq*evoltrunc::U1[i].qg() + C1sgqg*evoltrunc::U1[i].gg();
	  //complex <double> CTL3sgqq = C3sgqq + evoltrunc::U3[i].qq() + C1sgqq*evoltrunc::U2[i].qq() + C1sgqg*evoltrunc::U2[i].gq() + C2sgqq*evoltrunc::U1[i].qq() + C2sgqg*evoltrunc::U1[i].gq();
	  //complex <double> CTL3sgqg = C3sgqg + evoltrunc::U3[i].qg() + C1sgqq*evoltrunc::U2[i].qg() + C1sgqg*evoltrunc::U2[i].gg() + C2sgqq*evoltrunc::U1[i].qg() + C2sgqg*evoltrunc::U1[i].gg();
	  //complex <double> CTL1nsp = C1nsp + evoltrunc::U1NSp[i];
	  //complex <double> CTL1nsm = C1nsm + evoltrunc::U1NSm[i];
	  //complex <double> CTL1nsv = C1nsv + evoltrunc::U1NSv[i];
	  //complex <double> CTL2nsp = C2nsp + evoltrunc::U2NSp[i] + C1nsp*evoltrunc::U1NSp[i];
	  //complex <double> CTL2nsm = C2nsm + evoltrunc::U2NSm[i] + C1nsm*evoltrunc::U1NSm[i];
	  //complex <double> CTL2nsv = C2nsv + evoltrunc::U2NSv[i] + C1nsv*evoltrunc::U1NSv[i];
	  //complex <double> CTL3nsp = C3nsp + evoltrunc::U3NSp[i] + C1nsp*evoltrunc::U2NSp[i] + C2nsp*evoltrunc::U1NSp[i];
	  //complex <double> CTL3nsm = C3nsm + evoltrunc::U3NSm[i] + C1nsm*evoltrunc::U2NSm[i] + C2nsm*evoltrunc::U1NSm[i];
	  //complex <double> CTL3nsv = C3nsv + evoltrunc::U3NSv[i] + C1nsv*evoltrunc::U2NSv[i] + C2nsv*evoltrunc::U1NSv[i];

	  complex <double> CTL1sgqq = C1sgqq;
	  complex <double> CTL1sgqg = C1sgqg;
	  complex <double> CTL2sgqq = C2sgqq;
	  complex <double> CTL2sgqg = C2sgqg;
	  complex <double> CTL3sgqq = C3sgqq;
	  complex <double> CTL3sgqg = C3sgqg;
	  complex <double> CTL1nsp  = C1nsp ;
	  complex <double> CTL1nsm  = C1nsm ;
	  complex <double> CTL1nsv  = C1nsv ;
	  complex <double> CTL2nsp  = C2nsp ;
	  complex <double> CTL2nsm  = C2nsm ;
	  complex <double> CTL2nsv  = C2nsv ;
	  complex <double> CTL3nsp  = C3nsp ;
	  complex <double> CTL3nsm  = C3nsm ;
	  complex <double> CTL3nsv  = C3nsv ;
	  
	  //exponentiate here....
	  complex <double> aexpqq  = 1.;
	  complex <double> aexpqg  = 1.;
	  complex <double> aexpnsp = 1.;
	  complex <double> aexpnsm = 1.;
	  complex <double> aexpnsv = 1.;
  
	  if (opts.order_expc == 1)
	      aexpqg  = aexp;
	  else if (opts.order_expc == 2)
	    {
	      aexpqq  = pow(aexpB,CTL1sgqq);
	      aexpqg  = aexp*pow(aexpB,CTL2sgqg/CTL1sgqg);
	      aexpnsp = pow(aexpB,CTL1nsp);
	      aexpnsm = pow(aexpB,CTL1nsm);
	      aexpnsv = pow(aexpB,CTL1nsv);
	    }
	  else if (opts.order_expc == 3)
	    {
	      complex <double> C1,C2;
	      C1 = CTL1sgqq; C2 = CTL2sgqq;
	      aexpqq = exp(lamB*(C1)+lamC*((0.5*pow(C1,2) - C2))+lamD*(-beta1/beta0*C1)+lamC*(beta0*C1*(LQR)));
	      C1 = CTL2sgqg/CTL1sgqg; C2 = CTL3sgqg/CTL1sgqg;
	      aexpqg = aexp*exp(lamB*(C1)+lamC*((0.5*pow(C1,2) - C2))+lamD*(-beta1/beta0*C1)+lamC*(beta0*C1*(LQR)));
	      C1 = CTL1nsp; C2 = CTL2nsp;
	      aexpnsp = exp(lamB*(C1)+lamC*((0.5*pow(C1,2) - C2))+lamD*(-beta1/beta0*C1)+lamC*(beta0*C1*(LQR)));
	      C1 = CTL1nsm; C2 = CTL2nsm;
	      aexpnsm = exp(lamB*(C1)+lamC*((0.5*pow(C1,2) - C2))+lamD*(-beta1/beta0*C1)+lamC*(beta0*C1*(LQR)));
	      C1 = CTL1nsv; C2 = CTL2nsv;
	      aexpnsv = exp(lamB*(C1)+lamC*((0.5*pow(C1,2) - C2))+lamD*(-beta1/beta0*C1)+lamC*(beta0*C1*(LQR)));
	    }
	  else if (opts.order_expc == 4)
	    {
	      complex <double> C1,C2,C3;	      
	      C1 = CTL1sgqq; C2 = CTL2sgqq; C3 = CTL3sgqq;
	      aexpqq = exp(+lamB*(C1)+lamC*((0.5*pow(C1,2) - C2))+lamD*(-beta1/beta0*C1)+lamC*(beta0*C1*(LQR))+lamE*(-pow(beta1/beta0,2)*C1)+lamE*(beta2/beta0*C1)+lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))+lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))+lamH*(pow(beta1/beta0,2)*C1)+lamG*(2.*beta1*C1*LQR)+lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)+lamF*(pow(beta0,2)*C1*pow(LQR,2)));
	      C1 = CTL2sgqg/CTL1sgqg; C2 = CTL3sgqg/CTL1sgqg; C3 = - H1q/2.*C2  - (H2q/2.-pow(H1q,2)/8.)*C1  - (H3q/2.+pow(H1q,3)/16.);
	      aexpqg = aexp*exp(+lamB*(C1)+lamC*((0.5*pow(C1,2) - C2))+lamD*(-beta1/beta0*C1)+lamC*(beta0*C1*(LQR))+lamE*(-pow(beta1/beta0,2)*C1)+lamE*(beta2/beta0*C1)+lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))+lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))+lamH*(pow(beta1/beta0,2)*C1)+lamG*(2.*beta1*C1*LQR)+lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)+lamF*(pow(beta0,2)*C1*pow(LQR,2)));
	      C1 = CTL1nsp; C2 = CTL2nsp; C3 = CTL3nsp;
	      aexpnsp = exp(+lamB*(C1)+lamC*((0.5*pow(C1,2) - C2))+lamD*(-beta1/beta0*C1)+lamC*(beta0*C1*(LQR))+lamE*(-pow(beta1/beta0,2)*C1)+lamE*(beta2/beta0*C1)+lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))+lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))+lamH*(pow(beta1/beta0,2)*C1)+lamG*(2.*beta1*C1*LQR)+lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)+lamF*(pow(beta0,2)*C1*pow(LQR,2)));
	      C1 = CTL1nsm; C2 = CTL2nsm; C3 = CTL3nsm;
	      aexpnsm = exp(+lamB*(C1)+lamC*((0.5*pow(C1,2) - C2))+lamD*(-beta1/beta0*C1)+lamC*(beta0*C1*(LQR))+lamE*(-pow(beta1/beta0,2)*C1)+lamE*(beta2/beta0*C1)+lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))+lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))+lamH*(pow(beta1/beta0,2)*C1)+lamG*(2.*beta1*C1*LQR)+lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)+lamF*(pow(beta0,2)*C1*pow(LQR,2)));
	      C1 = CTL1nsv; C2 = CTL2nsv; C3 = CTL3nsv;
	      aexpnsv = exp(+lamB*(C1)+lamC*((0.5*pow(C1,2) - C2))+lamD*(-beta1/beta0*C1)+lamC*(beta0*C1*(LQR))+lamE*(-pow(beta1/beta0,2)*C1)+lamE*(beta2/beta0*C1)+lamF*(1./3.*(pow(C1,3)-3.*C1*C2+3.*C3))+lamG*(-beta1/pow(beta0,2)*(beta1*C1 - beta0*pow(C1,2) + 2.*beta0*C2))+lamH*(pow(beta1/beta0,2)*C1)+lamG*(2.*beta1*C1*LQR)+lamF*((-beta1*C1+beta0*pow(C1,2)-2.*beta0*C2)*LQR)+lamF*(pow(beta0,2)*C1*pow(LQR,2)));
	    }

	  C0qq[idx] = aexpqq;
	  CTL1sgqq *= aexpqq;
	  CTL2sgqq *= aexpqq;
	  CTL3sgqq *= aexpqq;
	  CTL1sgqg *= aexpqg;
	  CTL2sgqg *= aexpqg;
	  CTL3sgqg *= aexpqg;
	  CTL1nsp  *= aexpnsp;
	  CTL2nsp  *= aexpnsp;
	  CTL3nsp  *= aexpnsp;
	  CTL1nsm  *= aexpnsm;
	  CTL2nsm  *= aexpnsm;
	  CTL3nsm  *= aexpnsm;
	  CTL1nsv  *= aexpnsv;
	  CTL2nsv  *= aexpnsv;
	  CTL3nsv  *= aexpnsv;

	  
	  //Change of basis
	  double nf = double(resconst::NF);
	  complex <double> csg_qq,csg_qg,csg_gq,csg_gg,cns_p,cns_m,cns_v;
	  complex <double> cqqS_p_cqqbS,cqqS_m_cqqbS,cqqS,cqqbS,cqqV,cqqbV;
	  
	  csg_qq = CTL1sgqq;
	  csg_qg = CTL1sgqg;
	  cns_p = CTL1nsp;
	  cns_m = CTL1nsm;
	  cns_v = CTL1nsv;

	  cqqS_p_cqqbS = (csg_qq-cns_p)/nf;
	  cqqS_m_cqqbS = (cns_v-cns_m)/nf;
	  cqqS         = (cqqS_p_cqqbS+cqqS_m_cqqbS)/2.;
	  cqqbS        = (cqqS_p_cqqbS-cqqS_m_cqqbS)/2.;
	  cqqV         = (cns_p+cns_m)/2.;
	  cqqbV        = (cns_p-cns_m)/2.;

	  C1qq[idx]   = cqqV +cqqS;
	  C1qqb[idx]  = cqqbV+cqqbS;
	  C1qqp[idx]  = cqqS;
	  C1qqbp[idx] = cqqbS;
	  C1qg[idx]   = csg_qg/2./nf;

	  csg_qq = CTL2sgqq;
	  csg_qg = CTL2sgqg;
	  cns_p = CTL2nsp;
	  cns_m = CTL2nsm;
	  cns_v = CTL2nsv;

	  cqqS_p_cqqbS = (csg_qq-cns_p)/nf;
	  cqqS_m_cqqbS = (cns_v-cns_m)/nf;
	  cqqS         = (cqqS_p_cqqbS+cqqS_m_cqqbS)/2.;
	  cqqbS        = (cqqS_p_cqqbS-cqqS_m_cqqbS)/2.;
	  cqqV         = (cns_p+cns_m)/2.;
	  cqqbV        = (cns_p-cns_m)/2.;

	  C2qq[idx]   = cqqV +cqqS;
	  C2qqb[idx]  = cqqbV+cqqbS;
	  C2qqp[idx]  = cqqS;
	  C2qqbp[idx] = cqqbS;
	  C2qg[idx]   = csg_qg/2./nf;

	  csg_qq = CTL3sgqq;
	  csg_qg = CTL3sgqg;
	  cns_p = CTL3nsp;
	  cns_m = CTL3nsm;
	  cns_v = CTL3nsv;

	  cqqS_p_cqqbS = (csg_qq-cns_p)/nf;
	  cqqS_m_cqqbS = (cns_v-cns_m)/nf;
	  cqqS         = (cqqS_p_cqqbS+cqqS_m_cqqbS)/2.;
	  cqqbS        = (cqqS_p_cqqbS-cqqS_m_cqqbS)/2.;
	  cqqV         = (cns_p+cns_m)/2.;
	  cqqbV        = (cns_p-cns_m)/2.;

	  C3qq[idx]   = cqqV +cqqS;
	  C3qqb[idx]  = cqqbV+cqqbS;
	  C3qqp[idx]  = cqqS;
	  C3qqbp[idx] = cqqbS;
	  C3qg[idx]   = csg_qg/2./nf;
	}
    }  


  /*

  //exponentiate C
  for (int i = 0; i < dim*2; i++)
    {
      n    = (i<dim)? i : i-dim;
      sign = (i<dim)? mesq::positive : mesq::negative;
      if (opts.mellin1d)
	{
	  nn = n;
	  int idx = anomalous::index(nn,sign);

	  C1qq[idx]   = aexpqq[idx];
	  C1qg[idx]   = aexpqg[idx];

	  C2qq[idx]   = aexpqq[idx];
	  C2qqb[idx]  = aexpqqb[idx];
	  C2qqp[idx]  = aexpqqp[idx];
	  C2qqbp[idx] = aexpqqbp[idx];
	  C2qg[idx]   = aexpqg[idx];

	  C3qq[idx]   = aexpqq[idx];
	  C3qqb[idx]  = aexpqqb[idx];
	  C3qqp[idx]  = aexpqqp[idx];
	  C3qqbp[idx] = aexpqqbp[idx];
	  C3qg[idx]   = aexpqg[idx];
	}
    }
  */
  
}
