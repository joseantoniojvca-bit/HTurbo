#include "ctintmc.h"

#include "ctint.h"
#include "mesq.h"
#include "parton.h"
#include "settings.h"
#include "omegaintegr.h"
#include "pdf.h"
#include "scales.h"
#include "phasespace.h"
#include "resconst.h"
#include "constants.h"
#include "dyres_interface.h"
#include "dynnlo_interface.h"
#include "qt_interface.h"
#include "hpoly_interface.h"
#include "KinematicCuts.h"
#include "switch.h"
#include "itilde.h"
#include "isnan.h"
#include "alphaqed.h"

#include <iostream>
#include <string.h>
//#include <math.h>

using namespace std;
using namespace resconst;
using namespace constants;

const double ctlimit = 1e10;

void ctintmc::calc(double r[8], double f[], int order, bool fpc)
{
  for (int npdf = 0; npdf < opts.totpdf; npdf++)
    f[npdf] = 0.;

  //Jacobian of the change of variables from the unitary hypercube x[6] to the m, y, qt, phi, costh, philep boundaries
  double jac = 1.;
  bool status = true;

  double r1[2] = {r[0], r[1]};
  status = phasespace::gen_my(r1, jac, !opts.fixedorder, !opts.fixedorder);  //qtcut = !opts.fixedorder, qtswitching = !opts.fixedorder
  if (!status)
    {
      f[0] = 0.;
      return;
    }

  double qt;
  if (fpc) //for the FPC term integrate from qtfpc to qtcut
    {
      phasespace::gen_qt_fpc(r[2], jac);
      qt = phasespace::qt;
    }
  else if (opts.fixedorder) //In the fixed order calculation, integrate from qtcut to infinity
    {
      phasespace::gen_qt_ctfo(r[2], jac);
      qt = phasespace::qt;
      phasespace::set_qt(0.); //In the fixed order calculation evaluate kinematic cuts with qt=0
    }
  else
    {
      //Generate the boson transverse momentum between the integration boundaries
      double qtcut = max(opts.qtcut,opts.xqtcut*phasespace::m);
      double qtmn = max(qtcut, phasespace::qtmin);
      //phasespace::calcexpy();
      //double cosh2y=pow((phasespace::exppy+phasespace::expmy)*0.5,2);
      //double kinqtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+phasespace::m2,2)/(4*pow(opts.sroot,2)*cosh2y)-phasespace::m2)); //introduced max to avoid neqative argument of sqrt when y=ymax
      //double kinqtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+phasespace::m2,2)/(4*pow(opts.sroot,2))-phasespace::m2)); //introduced max to avoid negative argument of sqrt when y=ymax
      double kinqtlim = ctlimit; //There should not be any kinematic limit on qt, since the counterterm is evaluated with born level kinematic
      double switchqtlim = switching::qtlimit(phasespace::m);
      double qtlim = min(kinqtlim, switchqtlim);
      double qtmx = min(qtlim, phasespace::qtmax);
      if (qtmn >= qtmx)
	{
	  f[0] = 0.;
	  return;
	}
      status = phasespace::gen_qt(r[2], jac, qtlim, true);
      if (!status)
	{
	  f[0] = 0.;
	  return;
	}
      qt = phasespace::qt;
    }
  jac = jac *2.*qt;

  //generate boson 4-momentum, with m, qt, y, and phi
  phasespace::set_phiV(-M_PI+2.*M_PI*r[5]);

  double r2[2] = {r[3], r[4]};
  double fpc_weight = 1.;
  if (fpc)
    {
      double jj;
      phasespace::set_qt(0.);
      omegaintegr::genV4p();
      phasespace::gen_costhphi(r2, jj);
      phasespace::calcphilep();
      phasespace::genl4p();
      bool ctqt0 = Kinematics::Cuts::KeepThisEvent(phasespace::p3, phasespace::p4);

      phasespace::set_qt(qt);
      omegaintegr::genV4p();
      phasespace::gen_costhphi(r2, jac);
      phasespace::calcphilep();
      phasespace::genl4p();
      bool ctrec = Kinematics::Cuts::KeepThisEvent(phasespace::p3, phasespace::p4);

      if (ctqt0 == ctrec)
	  {
	    f[0]=0.;
	    return;
	  }
      fpc_weight = ctqt0 - ctrec;
    }
  else
    {
      //This function set the RF axes according to the selected qt-recoil prescription
      omegaintegr::genV4p();

      phasespace::gen_costhphi(r2, jac);
      phasespace::calcphilep();
      phasespace::genl4p();

      //apply lepton cuts
      if (opts.makecuts)
	if (!Kinematics::Cuts::KeepThisEvent(phasespace::p3, phasespace::p4))
	  {
	    f[0]=0.;
	    return;
	  }
    }
  
  double m = phasespace::m;
  double m2 = m*m;
  double y = phasespace::y;

  jac = jac/2./M_PI;

  //Set scales
  scales::set(m);
  scales::mcfm();
  scales::dyres(m);
  double muf = scales::fac;
  double mur = scales::ren;
  double Q   = scales::res;
  
  double LR, LF, LQ;
  if (order >= 2)
    LR = log(m2/pow(mur,2));
  LF = log(m2/pow(muf,2));
  LQ = log(m2/pow(Q,2));

  //double LR = 2.*log(phasespace::m/scales::ren);
  //double LF = 2.*log(phasespace::m/scales::fac);
  //double LQ = 2.*log(phasespace::m/scales::res);
  
  double LQ2 = pow(LQ,2);
  double LQ3 = pow(LQ,3);
  double LQ4 = pow(LQ,4);
  double LQ5 = pow(LQ,5);
  double LR2 = pow(LR,2);
  double LF2 = pow(LF,2);

  //   LL1,LL2,LL3,LL4: large log (squared) corresponding to eq. (136) 
  //   In this way normalization is fixed to dsigma/dqt2


  //double LL1 = Itilde(1)/pow(scales::res,2);
  //double LL2 = Itilde(2)/pow(scales::res,2);
  //double LL3 = Itilde(3)/pow(scales::res,2);
  //double LL4 = Itilde(4)/pow(scales::res,2);

  double x = qt/scales::res;
  double LL1 = itilde::besselk(x,1)/pow(scales::res,2)/m2;
  double LL2 = itilde::besselk(x,2)/pow(scales::res,2)/m2;
  double LL3 = itilde::besselk(x,3)/pow(scales::res,2)/m2;
  double LL4 = itilde::besselk(x,4)/pow(scales::res,2)/m2;
  double LL5 = itilde::besselk(x,5)/pow(scales::res,2)/m2;
  double LL6 = itilde::besselk(x,6)/pow(scales::res,2)/m2;
      
  
  //Scaled momentum fractions

  double cut=1e-7;
   
  double alfa = cut+(1.-cut)*r[6];
  double beta = cut+(1.-cut)*r[7];


  double qt2 = pow(phasespace::qt,2);
  double exppy = exp(y);
  double expmy = 1./exppy;
  double tau = sqrt(m2/pow(opts.sroot,2));

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
  
  double z1 = pow(x1,beta);
  double z2 = pow(x2,alfa);

  mesq::setpropagators(phasespace::m);
  double cthmom0 = 1.;
  double cthmom1 = phasespace::costh;
  double cthmom2 = pow(phasespace::costh,2);
  mesq::setmesq(cthmom0, cthmom1, cthmom2);
  
  dysetpdf_(opts.LHAPDFmember);

  //PDFs
  double fx1[2*MAXNF+1],fx2[2*MAXNF+1];
  //fdist_(opts.ih1,x1,scales::fac,fx1);
  //fdist_(opts.ih2,x2,scales::fac,fx2);
  pdf::fdist(1,x1,scales::fac,fx1);
  pdf::fdist(2,x2,scales::fac,fx2);

  //scaled PDFs
  double fx1p[2*MAXNF+1],fx2p[2*MAXNF+1];
  double f1Ap,f2Ap;
  double xx1 = pow(x1,(1.-beta));
  double xx2 = pow(x2,(1.-alfa));
  //fdist_(opts.ih1,xx1,scales::fac,fx1p);
  //fdist_(opts.ih2,xx2,scales::fac,fx2p);
  pdf::fdist(1,xx1,scales::fac,fx1p);
  pdf::fdist(2,xx2,scales::fac,fx2p);
  pdf::fdistphot(1,xx1,scales::fac,f1Ap);
  pdf::fdistphot(2,xx2,scales::fac,f2Ap);

  double sumfx1p = 0.;
  double sumfx2p = 0.;
  if (order >= 2)
    for (int f = 0; f < MAXNF; f++)
      {
	sumfx1p += fx1p[f]+fx1p[parton::charge_conj(parton::pdgid(f))];
	sumfx2p += fx2p[f]+fx2p[parton::charge_conj(parton::pdgid(f))];
      }

  double fx1pafx1pmanoti[2*MAXNF+1];
  double fx2pafx2pmanotj[2*MAXNF+1];
  memset(fx1pafx1pmanoti, 0, sizeof(fx1pafx1pmanoti));    
  memset(fx2pafx2pmanotj, 0, sizeof(fx2pafx2pmanotj));    
  if (order >= 2)
    for (int f = 0; f < MAXNF; f++)
      {
	fx1pafx1pmanoti[f] = fx1pafx1pmanoti[parton::charge_conj(parton::pdgid(f))] = sumfx1p - (fx1p[f]+fx1p[parton::charge_conj(parton::pdgid(f))]);
	fx2pafx2pmanotj[f] = fx2pafx2pmanotj[parton::charge_conj(parton::pdgid(f))] = sumfx2p - (fx2p[f]+fx2p[parton::charge_conj(parton::pdgid(f))]);
      }

  double sumposfx1p = 0.;
  double sumnegfx1p = 0.;
  double sumposfx2p = 0.;
  double sumnegfx2p = 0.;
  if (order >= 3)
    for (int f = 0; f < MAXNF; f++)
      {
	sumposfx1p += fx1p[f];
	sumnegfx1p += fx1p[parton::charge_conj(parton::pdgid(f))];
	sumposfx2p += fx2p[f];
	sumnegfx2p += fx2p[parton::charge_conj(parton::pdgid(f))];
      }
  

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


  double cz1a;
  double cz2b;
  double oz1a;
  double oz2b;
  double lz1a;
  double lz2b;
  double l1z1a;
  double l1z2b;
  double lzoz1a;
  double lzoz2b;
  //Preliminary caching
  if (order >= 1)
    {
      cz1a = z1;
      cz2b = z2;
      oz1a = 1./(1.-cz1a);
      oz2b = 1./(1.-cz2b);
      lz1a = log(cz1a);
      lz2b = log(cz2b);
      l1z1a = log(1.-cz1a);
      l1z2b = log(1.-cz2b);
      lzoz1a = log(1.-cz1a)/(1.-cz1a);
      lzoz2b = log(1.-cz2b)/(1.-cz2b);
    }

  //tc2ga1
  double p1qqc2qqhsreg1a;
  double p1qqc2qqhsd01a;
  double p1qqc2qqhsd11a;
  double p1qqc2qqhsdelta1a;
  double p1qqc2qqb1a;
  double p1qqc2qqp1a;
  double p1gqc2qg1a;
  double p1qgc2qq1a;
  double p1qgc2qqb1a;
  double p1qgc2qqp1a;
  double c2qqreghschm1a;
  double c2qqb1a;
  double c2qqp1a;
  double c2qghschm1a;
  double p1ggc2qg1a;

  double p1qqc2qqhsreg2b;
  double p1qqc2qqhsd02b;
  double p1qqc2qqhsd12b;
  double p1qqc2qqhsdelta2b;
  double p1qqc2qqb2b;
  double p1qqc2qqp2b;
  double p1gqc2qg2b;
  double p1qgc2qq2b;
  double p1qgc2qqb2b;
  double p1qgc2qqp2b;
  double c2qqreghschm2b;
  double c2qqb2b;
  double c2qqp2b;
  double c2qghschm2b;
  double p1ggc2qg2b;
  
  //tgamma3
  double p3qqv1a;           
  double p3qqvd01a;
  double p3qqvdelta1a;
  double p3qqbv1a;
  double p3sqq1a;
  double p3sqqb1a;
  double p3qg1a;
  double p3qqv2b;
  double p3qqvd02b;
  double p3qqvdelta2b;
  double p3qqbv2b;
  double p3sqq2b;
  double p3sqqb2b;
  double p3qg2b;
  if (order >= 3)
    {
      //if (cz1a >= 1) continue;
      cacheuv_(cz1a);
      p1qqc2qqhsreg1a   =  p1qqc2qqhsreg_(cz1a);   
      p1qqc2qqhsd01a    =  p1qqc2qqhsd0_(cz1a);	  
      p1qqc2qqhsd11a    =  p1qqc2qqhsd1_(cz1a);	  
      p1qqc2qqhsdelta1a =  p1qqc2qqhsdelta_(cz1a); 
      p1qqc2qqb1a	     =  p1qqc2qqb_(cz1a);	  
      p1qqc2qqp1a	     =  p1qqc2qqp_(cz1a);	  
      p1gqc2qg1a	     =  p1gqc2qg_(cz1a);	  
      p1qgc2qq1a	     =  p1qgc2qq_(cz1a);	  
      p1qgc2qqb1a	     =  p1qgc2qqb_(cz1a);	  
      p1qgc2qqp1a	     =  p1qgc2qqp_(cz1a);	  
      c2qqreghschm1a    =  c2qqreghschm_(cz1a);	  
      c2qqb1a	     =  c2qqb_(cz1a);		  
      c2qqp1a	     =  c2qqp_(cz1a);		  
      c2qghschm1a       =  c2qghschm_(cz1a);       
      p1ggc2qg1a        =  p1ggc2qg_(cz1a);
      p3qqv1a      = p3qqv_(cz1a);
      p3qqvd01a	= p3qqvd0_(cz1a);
      p3qqvdelta1a	= p3qqvdelta_(cz1a);
      p3qqbv1a	= p3qqbv_(cz1a);
      p3sqq1a	= p3sqq_(cz1a);
      p3sqqb1a	= p3sqqb_(cz1a);
      p3qg1a	= p3qg_(cz1a);
    }
  if (order >= 3)
    {
      //if (cz2b >= 1) continue;
      cacheuv_(cz2b);
      p1qqc2qqhsreg2b   =  p1qqc2qqhsreg_(cz2b);   
      p1qqc2qqhsd02b    =  p1qqc2qqhsd0_(cz2b);	  
      p1qqc2qqhsd12b    =  p1qqc2qqhsd1_(cz2b);	  
      p1qqc2qqhsdelta2b =  p1qqc2qqhsdelta_(cz2b); 
      p1qqc2qqb2b	     =  p1qqc2qqb_(cz2b);	  
      p1qqc2qqp2b	     =  p1qqc2qqp_(cz2b);	  
      p1gqc2qg2b	     =  p1gqc2qg_(cz2b);	  
      p1qgc2qq2b	     =  p1qgc2qq_(cz2b);	  
      p1qgc2qqb2b	     =  p1qgc2qqb_(cz2b);	  
      p1qgc2qqp2b	     =  p1qgc2qqp_(cz2b);	  
      c2qqreghschm2b    =  c2qqreghschm_(cz2b);	  
      c2qqb2b	     =  c2qqb_(cz2b);		  
      c2qqp2b	     =  c2qqp_(cz2b);		  
      c2qghschm2b       =  c2qghschm_(cz2b);       
      p1ggc2qg2b        =  p1ggc2qg_(cz2b);
      p3qqv2b      = p3qqv_(cz2b);
      p3qqvd02b	= p3qqvd0_(cz2b);
      p3qqvdelta2b	= p3qqvdelta_(cz2b);
      p3qqbv2b	= p3qqbv_(cz2b);
      p3sqq2b	= p3sqq_(cz2b);
      p3sqqb2b	= p3sqqb_(cz2b);
      p3qg2b	= p3qg_(cz2b);
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
      if (cz1a >= 1) continue;
	  
      //H1st non delta terms
      th1st += (fx1p[i]*cqq_(cz1a)+fx1p[g]*cqg_(cz1a))*(-lx1)*fx2[j];

      //H1st muf dependence, gammaqq and gammaqg:
      th1stF += (-lx1*((fx1p[i]-fx1[i]*cz1a)*pqq_(cz1a)+fx1p[g]*dypqg_(cz1a)))*fx2[j];
      th1stF_qed += (-lx1*((fx1p[i]-fx1[i]*cz1a)*3./4.*pqq_(cz1a)))*fx2[j] *(parton::chsq(i)+parton::chsq(j))/2.;
      th1stF_qed += (-lx1*( +f1Ap*6.*dypqg_(cz1a)   ))*fx2[j] *parton::chsq(j) * facqedis; //photon PDF

      th1stF += -pqqintx1*tdelta;
      th1stF_qed += -3./4.*pqqintx1*tdelta*(parton::chsq(i)+parton::chsq(j))/2.;
      
      //beta loop (second leg)
      if (cz2b >= 1) continue;

      //H1st non delta terms
      th1st += (fx2p[j]*cqq_(cz2b)+fx2p[g]*cqg_(cz2b))*(-lx2)*fx1[i];
      
      //H1st muf dependence, gammaqq and gammaqg:
      th1stF += (-lx2*((fx2p[j]-fx2[j]*cz2b)*pqq_(cz2b)+fx2p[g]*dypqg_(cz2b)))*fx1[i];
      th1stF_qed += (-lx2*((fx2p[j]-fx2[j]*cz2b)*3./4.*pqq_(cz2b)))*fx1[i] *(parton::chsq(i)+parton::chsq(j))/2.;
      th1stF_qed += (-lx2*(f2Ap*6.*dypqg_(cz2b) ))*fx1[i] *parton::chsq(i) * facqedis; //photon PDF

      th1stF += -pqqintx2*tdelta;
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
      
      if (cz1a >= 1) continue;
               
      //(gamma+gamma)*(gamma+gamma) term

      //First part: one gamma for each leg
      diffg1f += (-lx1*(fx1p[i]-fx1[i]*cz1a)*pqq_(cz1a) - pqqintx1*fx1[i]) ;
      diffg10 += -lx1*fx1p[g]*dypqg_(cz1a) ;
      
      //Second part: gamma*gamma terms
      //Pij * Pjk = D1ijjk (log(1-z)/(1-z))_+ + D0ijjk/(1-z)_+ 
      //          + Pijjk(z) + Deltaijjk delta(1-z)
      //First leg
      diff1 += (-lx1*((fx1p[i]-fx1[i]*cz1a)
		      *(resconst::D0qqqq*oz1a+resconst::D1qqqq*lzoz1a)
		      +fx1p[i]*pqqqq_(cz1a)+fx1p[g]*(pqqqg_(cz1a)+pqggg_(cz1a)))
		+(resconst::Deltaqqqq-resconst::D0qqqq*d0intx1-resconst::D1qqqq*d1intx1)
		*fx1[i]) ;

      //Include Pqggq
      diff1 += -lx1*sumfx1p*pqggq_(cz1a) ;
      //End of (gamma+gamma)*(gamma+gamma) term

      //Start  (C+C)*(gamma+gamma) term
      //C first leg, gamma second leg
      diffc1f += (-lx1*fx1p[i]*cqq_(cz1a)+resconst::C1qqdelta*fx1[i]) ;
      diffc10 += -lx1*fx1p[g]*cqg_(cz1a) ;
      
      //C*gamma: first leg (ignore delta term in Cqq: taken into account with th1stF)
      tcga += (fx1p[i]*cqqpqq_(cz1a)+fx1p[g]*(cqqpqg_(cz1a)+cqgpgg_(cz1a)))*(-lx1)*fx2[j] ;

      //Add Cqg*Pgq contribution
      tcga += sumfx1p*cqgpgq_(cz1a)*(-lx1)*fx2[j] ;

      //Start 2-loop AP
      // Gluon + pure singlet
      //f == gluon piece
      tgamma2 += fx1p[g]*p2qg_(cz1a)*(-lx1)*fx2[j] ;

      //f != gluon piece
      tgamma2 += sumfx1p*p2qqs_(cz1a)*(-lx1)*fx2[j] ;

      //P2qq non-singlet: regular part
      tgamma2 += fx1p[i]*p2qqv_(cz1a)*(-lx1)*fx2[j] ;

      //P2qq non-singlet: 1/(1-z)_+
      tgamma2 += 2./3.*resconst::Kappa*(-lx1*(fx1p[i]-fx1[i]*cz1a)*oz1a-d0intx1*fx1[i])*fx2[j] ;

      //P2qqb non singlet
      tgamma2 += fx1p[im]*p2qqbv_(cz1a)*(-lx1)*fx2[j] ;

      //beta loop
      double diffg2f = 0;
      double diffg20 = 0;
      double diff2 = 0;
      double diffc2f = 0;
      double diffc20 = 0;
      if (cz2b >= 1) continue;
	  
      //(gamma+gamma)*(gamma+gamma) term

      //First part: one gamma for each leg
      diffg2f += (-lx2*(fx2p[j]-fx2[j]*cz2b)*pqq_(cz2b) - pqqintx2*fx2[j]) ;
      diffg20 += -lx2*fx2p[g]*dypqg_(cz2b) ;
	  
      //Second part: gamma*gamma terms
      //Pij * Pjk = D1ijjk (log(1-z)/(1-z))_+ + D0ijjk/(1-z)_+ 
      //          + Pijjk(z) + Deltaijjk delta(1-z)
      //Second leg
      diff2 += (-lx2*((fx2p[j]-fx2[j]*cz2b)
		      *(resconst::D0qqqq*oz2b+resconst::D1qqqq*lzoz2b)
		      +fx2p[j]*pqqqq_(cz2b)+fx2p[g]*(pqqqg_(cz2b)+pqggg_(cz2b)))
		+(resconst::Deltaqqqq-resconst::D0qqqq*d0intx2-resconst::D1qqqq*d1intx2)
		*fx2[j]);
      //Include Pqggq
      diff2 += -lx2*sumfx2p*pqggq_(cz2b) ;
      //End of (gamma+gamma)*(gamma+gamma) term

      //Start  (C+C)*(gamma+gamma) term
      //gamma first leg, C second leg
      diffc2f += (-lx2*fx2p[j]*cqq_(cz2b)+resconst::C1qqdelta*fx2[j]) ;
      diffc20 += -lx2*fx2p[g]*cqg_(cz2b) ;

      //C*gamma: second leg (ignore delta term in Cqq: taken into account with th1stF)
      tcga += (fx2p[j]*cqqpqq_(cz2b)+fx2p[g]*(cqqpqg_(cz2b)+cqgpgg_(cz2b)))*(-lx2)*fx1[i] ;

      //Add Cqg*Pgq contribution
      tcga += sumfx2p*cqgpgq_(cz2b)*(-lx2)*fx1[i] ;

      //Start 2-loop AP
      // Gluon + pure singlet
      //f == gluon piece
      tgamma2 += fx2p[g]*p2qg_(cz2b)*(-lx2)*fx1[i] ;

      //f != gluon piece
      tgamma2 += sumfx2p*p2qqs_(cz2b)*(-lx2)*fx1[i] ;

      //P2qq non-singlet: regular part
      tgamma2 += fx2p[j]*p2qqv_(cz2b)*(-lx2)*fx1[i] ;

      //P2qq non-singlet: 1/(1-z)_+
      tgamma2 += 2./3.*resconst::Kappa*(-lx2*(fx2p[j]-fx2[j]*cz2b)*oz2b-d0intx2*fx2[j])*fx1[i] ;

      //P2qqb non singlet
      tgamma2 += fx2p[jm]*p2qqbv_(cz2b)*(-lx2)*fx1[i] ;

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

      //cz1[a] -> cz1a
      //fx1p[a] -> fx1p
      //* abint::abw[a] ->
      //z1[a] -> z1a
      //1[a] -> 1a
      
      //cz2[b] -> cz2b
      //fx2p[b] -> fx2p
      //* abint::abw[b] ->
      //z2[b] -> z2b
      //2[b] -> 2b

	  
      //
      // ----------------------------------------------------------------------------------------------------------------------------------------
      // tgagaga
      // ----------------------------------------------------------------------------------------------------------------------------------------	   
      //Leg 1
      cacheuv_(cz1a);
      //if (cz1a >= 1) continue;
      // qqb Channel
      tgagaga  += ctint::intfuncqt(fx1p[i],fx1[i],x1,cz1a,p1qqp1qqp1qqreg_(cz1a),p1qqp1qqp1qqd0_(cz1a),p1qqp1qqp1qqd1_(cz1a),p1qqp1qqp1qqd2_(cz1a),p1qqp1qqp1qqdelta_(cz1a))*fx2[j] ; //PqqPqqPqq
      tgagaga  += - 2.*lx1*(sumfx1p)*p1gqp1qqp1qg_(cz1a)*fx2[j] ; //PqqPqgPgq
      // PqqPqgPgq The Contribution is the same as the previous one ->  we simply put a factor 2 
      tgagaga  += - lx1*(sumfx1p)*p1qgp1ggp1gq_(cz1a)*fx2[j] ;		                                                                                                            //PqgPggPgq
      // qg Channel
      tgagaga += - lx1*fx1p[g]*p1qqp1qqp1qg_(cz1a)*fx2[j]; //PqqPqqPqg
      tgagaga += - lx1*fx1p[g]*p1qqp1qgp1gg_(cz1a)*fx2[j]; //PqqPqgPgg
      tgagaga += - lx1*fx1p[g]*p1qgp1gqp1qg_(cz1a)*fx2[j]*2.*NF;//PqgPgqPqg ! Candidate for Nf factor?!!!!
      tgagaga += - lx1*fx1p[g]*p1qgp1ggp1gg_(cz1a)*fx2[j]; // PqgPggPgg 
      // Start contributions used for simple products contributions
      diffqq1  += (-lx1*(fx1p[i]-fx1[i]*cz1a)*pqq_(cz1a) - pqqintx1*fx1[i]); // porque no se puede usar la de NNLO? si se usa fuera del loop?
      tgagaqL1 += (-lx1*((fx1p[i]-fx1[i]*cz1a) *(D0qqqq/(1.-cz1a)+D1qqqq*log(1.-cz1a)/(1.-cz1a))
			 +fx1p[i]*pqqqq_(cz1a))+(Deltaqqqq-D0qqqq*d0intx1-D1qqqq*d1intx1)*fx1[i]);
      tgagaqL1 += -lx1*(sumfx1p)*pqggq_(cz1a);
      diffg10  += -lx1*fx1p[g]*dypqg_(cz1a);
      tgagagL1 += -lx1*(fx1p[g]*(pqqqg_(cz1a)+pqggg_(cz1a)));
	      
      //Leg 2
      cacheuv_(cz2b);
      //if (cz2b >= 1) continue;
      // qqb Channel
      tgagaga += ctint::intfuncqt(fx2p[j],fx2[j],x2,cz2b,p1qqp1qqp1qqreg_(cz2b),p1qqp1qqp1qqd0_(cz2b),p1qqp1qqp1qqd1_(cz2b),p1qqp1qqp1qqd2_(cz2b),p1qqp1qqp1qqdelta_(cz2b))*fx1[i] ; //PqqPqqPqq
      tgagaga += - 2.*lx2*(sumfx2p)*p1gqp1qqp1qg_(cz2b)*fx1[i]; //PqqPqgPgq
      // PqqPqgPgq The Contribution is the same as the previous one ->  we simply put a factor 2  
      tgagaga += - lx2*(sumfx2p)*p1qgp1ggp1gq_(cz2b)*fx1[i] ; //PqgPggPgq
      // qg Channel
      tgagaga += - lx2*fx2p[g]*p1qqp1qqp1qg_(cz2b)*fx1[i]; //PqqPqqPqg 
      tgagaga += - lx2*fx2p[g]*p1qqp1qgp1gg_(cz2b)*fx1[i]; //PqqPqgPgg 
      tgagaga += - lx2*fx2p[g]*p1qgp1gqp1qg_(cz2b)*fx1[i]*2.*NF;//PqgPgqPqg  ! Candidate for Nf factor?!!!!
      tgagaga += - lx2*fx2p[g]*p1qgp1ggp1gg_(cz2b)*fx1[i]; //PqgPggPgg 
      // Start contributions for simple products contributions ! Checked 09.05.2020
      diffqq2  += (-lx2*(fx2p[j]-fx2[j]*cz2b)*pqq_(cz2b)- pqqint_(x2)*fx2[j]);
      diffg20  += -lx2*fx2p[g]*dypqg_(cz2b);
      tgagaqL2 += (-lx2*((fx2p[j]-fx2[j]*cz2b) *(D0qqqq/(1.-cz2b)+D1qqqq*log(1.-cz2b)/(1.-cz2b))
			 +fx2p[j]*pqqqq_(cz2b))+(Deltaqqqq-D0qqqq*d0intx2-D1qqqq*d1intx2)*fx2[j]);
      tgagaqL2 += -lx2*(sumfx2p)*pqggq_(cz2b);		
      tgagagL2 += -lx2*(fx2p[g]*(pqqqg_(cz2b)+pqggg_(cz2b)));		

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
      cacheuv_(cz1a);
      //if (cz1a >= 1) continue;
      // qqb Channel
      // CqqPqqPqq Leg 1 Typo solved on 09.05.2020
      tcgaga += - lx1*fx1p[i]*c1qqp1qqp1qq_(cz1a)*fx2[j];
      // CqqPqgPgq Leg 1 Checked 09.05.2020
      tcgaga += - lx1*(sumfx1p)*c1qqp1qgp1gq_(cz1a)*fx2[j];
      // CqgPgqPqq Leg 1 Checked 09.05.2020
      tcgaga += - lx1*(sumfx1p)*c1qgp1gqp1qq_(cz1a)*fx2[j];
      // CqgPggPgq Leg 1 Checked 09.05.2020
      tcgaga += - lx1*(sumfx1p)*c1qgp1ggp1gq_(cz1a)*fx2[j];
      // End qqb Channel
      //       ! qg Channel Checked 09.05.2020
      // CqqPqqPqg Leg 1
      tcgaga += - lx1*fx1p[g]*c1qqp1qqp1qg_(cz1a)*fx2[j];
      // CqqPqgPgg Leg 1
      tcgaga += - lx1*fx1p[g]*c1qqp1qgp1gg_(cz1a)*fx2[j];
      // CqgPgqPqg Leg 1
      tcgaga += - lx1*fx1p[g]*c1qgp1gqp1qg_(cz1a)*fx2[j]*2.*NF; // Candidate for NF factor?!!!!
      // CqgPggPgg Leg 1
      tcgaga += - lx1*fx1p[g]*c1qgp1ggp1gg_(cz1a)*fx2[j];		
      // Start contributions for simple products contributions Checked 09.05.2020		
      tcgahqL1 += (fx1p[i]*cqqpqq_(cz1a))*(-lx1);
      tcgahqL1 += + sumfx1p*cqgpgq_(cz1a)*(-lx1);
      tcgahgL1 += fx1p[g]*(cqqpqg_(cz1a)+cqgpgg_(cz1a))*(-lx1);
      diffc1f += (-lx1*fx1p[i]*cqq_(cz1a)+0.0*fx1[i]);

      cacheuv_(cz2b);
      //if (cz2b >= 1) continue;
      // qqb Channel		
      // CqqPqqPqq Leg 2 Typo solved on 09.05.2020
      tcgaga += - lx2*fx2p[j]*c1qqp1qqp1qq_(cz2b)*fx1[i];
      // CqqPqgPgq Leg 2 Checked 09.05.2020
      tcgaga += - lx2*(sumfx2p)*c1qqp1qgp1gq_(cz2b)*fx1[i];
      // CqgPgqPqq Leg 2 Checked 09.05.2020
      tcgaga += - lx2*(sumfx2p)*c1qgp1gqp1qq_(cz2b)*fx1[i];
      // CqgPggPgq Leg 2 Checked 09.05.2020
      tcgaga += - lx2*(sumfx2p)*c1qgp1ggp1gq_(cz2b)*fx1[i];
      //       ! qg Channel Checked 09.05.2020
      // CqqPqqPqg Leg 2 Typo solved 09.05.2020
      tcgaga += - lx2*fx2p[g]*c1qqp1qqp1qg_(cz2b)*fx1[i];
      // CqqPqgPgg Leg 2 Typo solved 09.05.2020
      tcgaga += - lx2*fx2p[g]*c1qqp1qgp1gg_(cz2b)*fx1[i];
      // CqgPgqPqg Leg 2 Typo solved 09.05.2020
      tcgaga += - lx2*fx2p[g]*c1qgp1gqp1qg_(cz2b)*fx1[i]*2.*NF; // Candidate for NF factor?!!!!
      // CqgPggPgg Leg 2 Typo solved 09.05.2020
      tcgaga += - lx2*fx2p[g]*c1qgp1ggp1gg_(cz2b)*fx1[i];
      // Start contributions for simple products contributions Checked 09.05.2020
      tcgahqL2 += (fx2p[j]*cqqpqq_(cz2b))*(-lx2);
      tcgahqL2 += + sumfx2p*cqgpgq_(cz2b)*(-lx2);
      tcgahgL2 += fx2p[g]*(cqqpqg_(cz2b)+cqgpgg_(cz2b))*(-lx2);
      diffc2f += (-lx2*fx2p[j]*cqq_(cz2b)+0.0*fx2[j]);

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
      //if (cz1a >= 1) continue;
      // qqb Channel
      tc2ga1 += + ctint::intfuncqt(fx1p[i],fx1[i],x1,cz1a,p1qqc2qqhsreg1a,p1qqc2qqhsd01a,p1qqc2qqhsd11a,0.,p1qqc2qqhsdelta1a)*fx2[j]; // P1qqC2qq
      // qq
      tc2ga1 += + fx1p[im]*p1qqc2qqb1a*(-lx1)*fx2[j]; // P1qqC2qqb
      // qqp
      tc2ga1 += + (-lx1)*(fx1pafx1pmanoti[i])*p1qqc2qqp1a*fx2[j]; // P1qqC2qqp
      tc2ga1 += + (-lx1)*(sumfx1p)*p1gqc2qg1a*fx2[j]; // P1gqC2qg
      // qg Channel
      // General C2qq*P1qg
      tc2ga1 += + (-lx1)*fx1p[g]*p1qgc2qq1a*fx2[j];
      tc2ga1 += + (-lx1)*fx1p[g]*p1qgc2qqb1a*fx2[j];
      tc2ga1 += + (-lx1)*fx1p[g]*p1qgc2qqp1a*fx2[j]*2.*(NF-1.);//OK
      tc2ga1 += + (-lx1)*fx1p[g]*p1ggc2qg1a*fx2[j]; // P1ggC2qg
      // Start contributions for simple products contributions
      diff1  = -lx1*(fx1p[i]-fx1[i]*cz1a)*H2qqD0/(1.-cz1a);
      C2qqL1 +=  0.5*diff1;
      C2qqL1 += -0.5*H2qqD0*d0intx1*fx1[i];
      C2qqL1 += + fx1p[i]*c2qqreghschm1a*(-lx1);
      C2qqL1 += + fx1p[im]*c2qqb1a*(-lx1);
      C2qqL1 += + fx1pafx1pmanoti[i]*c2qqp1a*(-lx1);
      C2qgL1 += fx1p[g]*c2qghschm1a*(-lx1);

      //Leg 2
      //if (cz2b >= 1) continue;
      // qqb Channel
      tc2ga1 += + ctint::intfuncqt(fx2p[j],fx2[j],x2,cz2b,p1qqc2qqhsreg2b,p1qqc2qqhsd02b,p1qqc2qqhsd12b,0.,p1qqc2qqhsdelta2b)*fx1[i]; // P1qqC2qq
      // qq
      tc2ga1 += + fx2p[jm]*p1qqc2qqb2b*(-lx2)*fx1[i]; // P1qqC2qqb
      // qqp
      tc2ga1 += + (-lx2)*(fx2pafx2pmanotj[j])*p1qqc2qqp2b*fx1[i]; // P1qqC2qqp
      tc2ga1 += + (-lx2)*(sumfx2p)*p1gqc2qg2b*fx1[i]; // P1gqC2qg
      // qg Channel
      // General C2qq*P1qg
      tc2ga1 += + (-lx2)*fx2p[g]*p1qgc2qq2b*fx1[i];
      tc2ga1 += + (-lx2)*fx2p[g]*p1qgc2qqb2b*fx1[i];
      tc2ga1 += + (-lx2)*fx2p[g]*p1qgc2qqp2b*fx1[i]*2.*(NF-1.); // OK
      tc2ga1 += + (-lx2)*fx2p[g]*p1ggc2qg2b*fx1[i]; // P1ggC2qg
      // Start contributions for simple products contributions
      diff2  = -lx2*(fx2p[j]-fx2[j]*cz2b)*H2qqD0/(1.-cz2b);
      C2qqL2 += 0.5*diff2;
      C2qqL2 += -0.5*H2qqD0*d0intx2*fx2[j];
      C2qqL2 += + fx2p[j]*c2qqreghschm2b*(-lx2);
      C2qqL2 += + fx2p[jm]*c2qqb2b*(-lx2);
      C2qqL2 += + fx2pafx2pmanotj[j]*c2qqp2b*(-lx2);
      C2qgL2 += fx2p[g]*c2qghschm2b*(-lx2);

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
      cacheuv_(cz1a);
      //if (cz1a >= 1) continue;
      // tc1ga2 Checked 10.05.2020
      // qq Channel
      // Singlet case C1qqP2qqS
      //  Leg 1    
      tc1ga2 += + sumfx1p * c1qqp2qqs_(cz1a)*(-lx1)*fx2[j];
      // c------ Non-singlet case C1qqP2qqV
      // // Leg 1
      tc1ga2 += + fx1p[i]*c1qqp2qqv_(cz1a)*(-lx1)*fx2[j];
      // // Non-singlet case C1qqP2qqbV
      // // Leg 1
      tc1ga2 += + fx1p[im]*c1qqp2qqbv_(cz1a)*(-lx1)*fx2[j];		
      // C1qgP2gq
      // Leg 1
      tc1ga2 += + (-lx1)*(sumfx1p)*c1qgp2gq_(cz1a)*fx2[j];
      // qg Channel     
      // Leg 1
      tc1ga2 += + fx1p[g]*(c1qqp2qg_(cz1a)+c1qgp2gg_(cz1a))*(-lx1)*fx2[j];
      // Start contributions for simple products contributions
      gamma2qqL1 += sumfx1p*p2qqs_(cz1a)*(-lx1);
      gamma2qqL1 += + fx1p[i]*p2qqv_(cz1a)*(-lx1);
      gamma2qqL1 += fx1[i]*resconst::Delta2qq;
      gamma2qqL1 += + 2./3.*Kappa*(-lx1*(fx1p[i]-fx1[i]*cz1a)/(1.-cz1a)-d0intx1*fx1[i]);
      gamma2qqL1 += + fx1p[im]*p2qqbv_(cz1a)*(-lx1);
      gamma2qgL1 += fx1p[g]*p2qg_(cz1a)*(-lx1);	 

      cacheuv_(cz2b);
      //if (cz2b >= 1) continue;
      // qqb Channel
      // Leg 2
      tc1ga2 += + sumfx2p * c1qqp2qqs_(cz2b)*(-lx2)*fx1[i];
      // c------ Non-singlet case C1qqP2qqV
      // Leg 2
      tc1ga2 += + fx2p[j]*c1qqp2qqv_(cz2b)*(-lx2)*fx1[i];
      // Non-singlet case C1qqP2qqbV
      // Leg 2
      tc1ga2 += + fx2p[jm]*c1qqp2qqbv_(cz2b)*(-lx2)*fx1[i];
      // C1qgP2gq
      // Leg 2
      tc1ga2 += + (-lx2)*(sumfx2p)*c1qgp2gq_(cz2b)*fx1[i];
      // qg Channel     
      // Leg 2
      tc1ga2 += + fx2p[g]*(c1qqp2qg_(cz2b)+c1qgp2gg_(cz2b))*(-lx2)*fx1[i];
      // Start contributions for simple products contributions
      gamma2qqL2 +=  sumfx2p*p2qqs_(cz2b)*(-lx2);
      gamma2qqL2 += + fx2p[j]*p2qqv_(cz2b)*(-lx2);
      gamma2qqL2 +=   fx2[j]*resconst::Delta2qq;
      gamma2qqL2 += + 2./3.*Kappa*(-lx2*(fx2p[j]-fx2[j]*cz2b)/(1.-cz2b)-d0intx2*fx2[j]);
      gamma2qqL2 += + fx2p[jm]*p2qqbv_(cz2b)*(-lx2);
      gamma2qgL2 += fx2p[g]*p2qg_(cz2b)*(-lx2);

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
      cacheuv_(cz1a);
      //if (cz1a >= 1) continue;
      // tga1ga2 Checked 10.05.2020
      // Usage IntFuncqT(fxip,fxi0,xxi0,z,Preg,PD0,PD1,PD2,Pdelta)
      // qq Channel
      // Leg 1
      tga1ga2 += + 2.*0.5*ctint::intfuncqt(fx1p[i],fx1[i],x1,cz1a,p1qqp2qqvreg_(cz1a),p1qqp2qqvd0_(cz1a),p1qqp2qqvd1_(cz1a),0.,p1qqp2qqvdelta_(cz1a))*fx2[j];
      // Singlet case P2qqS
      // Leg 1
      tga1ga2 += + 2.*0.5*sumfx1p *p1qqp2qqs_(cz1a)*(-lx1)*fx2[j];
      // P2qqb non singlet P2qqbV
      // Leg 1
      tga1ga2 += + 2.*0.5*fx1p[im]*p1qqp2qqbv_(cz1a)*(-lx1)*fx2[j];
      // P1qgP2gq +P1gqP2qg
      // Leg 1
      tga1ga2 += + 0.5*(-lx1)*(sumfx1p)*(p1qgp2gq_(cz1a)+p1gqp2qg_(cz1a))*fx2[j];  
      // End qqb Channel
      // qg Channel
      // P1qqP2qg + P1qgP2gg + P1ggP2qg
      // Leg 1
      tga1ga2 += + 0.5*(-lx1)*fx1p[g]*(p1qqp2qg_(cz1a)+p1qgp2gg_(cz1a)+p1ggp2qg_(cz1a))*fx2[j];
      // P1qgP2qq (all contributions inside)
      // Leg 1
      tga1ga2 += + 0.5*(-lx1)*fx1p[g]*p1qgp2qq_(cz1a)*fx2[j];

      cacheuv_(cz2b);
      //if (cz2b >= 1) continue;
      // qqb Channel
      // Leg 2
      tga1ga2 += + 2.*0.5*ctint::intfuncqt(fx2p[j],fx2[j],x2,cz2b,p1qqp2qqvreg_(cz2b),p1qqp2qqvd0_(cz2b),p1qqp2qqvd1_(cz2b),0.,p1qqp2qqvdelta_(cz2b))*fx1[i];
      // Singlet case P2qqS
      // Leg 2
      tga1ga2 += + 2.*0.5*sumfx2p *p1qqp2qqs_(cz2b)*(-lx2)*fx1[i];
      // P2qqb non singlet P2qqbV
      // Leg 2
      tga1ga2 += + 2.*0.5*fx2p[jm]*p1qqp2qqbv_(cz2b)*(-lx2)*fx1[i];
      // P1qgP2gq +P1gqP2qg
      // Leg 2
      tga1ga2 += + 0.5*(-lx2)*(sumfx2p)*(p1qgp2gq_(cz2b)+p1gqp2qg_(cz2b))*fx1[i];
      // End qqb Channel
      // qg Channel
      // P1qqP2qg + P1qgP2gg + P1ggP2qg
      // Leg 2
      tga1ga2 += + 0.5*(-lx2)*fx2p[g]*(p1qqp2qg_(cz2b)+p1qgp2gg_(cz2b)+p1ggp2qg_(cz2b))*fx1[i];
      // P1qgP2qq (all contributions inside)
      // Leg 2
      tga1ga2 += + 0.5*(-lx2)*fx2p[g]*p1qgp2qq_(cz2b)*fx1[i];

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
      cacheuv_(cz1a);
      //if (cz1a >= 1) continue;
      //qq Channel 
      tgamma3 += + ctint::intfuncqt(fx1p[i],fx1[i],x1,cz1a,p3qqv1a,p3qqvd01a,0.,0.,p3qqvdelta1a)*fx2[j]; //non singlet P3qqV      
      tgamma3 += + fx1p[im]*p3qqbv1a*(-lx1)*fx2[j];		                                           //non singlet P3qqbV
      //singlet
      if (j > g)
	{
	  tgamma3 += + (sumnegfx1p)*p3sqq1a*(-lx1)*fx2[j];
	  tgamma3 += + (sumposfx1p)*p3sqqb1a*(-lx1)*fx2[j];
	}
      else
	{
	  tgamma3 += + (sumposfx1p)*p3sqq1a*(-lx1)*fx2[j];
	  tgamma3 += + (sumnegfx1p)*p3sqqb1a*(-lx1)*fx2[j];
	}
      //qg Channel 
      tgamma3+=+fx1p[g]*p3qg1a*(-lx1)*fx2[j];

      //Leg 2
      cacheuv_(cz2b);
      //if (cz2b >= 1) continue;
      //qq Channel
      tgamma3 += + ctint::intfuncqt(fx2p[j],fx2[j],x2,cz2b,p3qqv2b,p3qqvd02b,0.,0.,p3qqvdelta2b)*fx1[i]; //non singlet P3qqV
      tgamma3 += + fx2p[jm]*p3qqbv2b*(-lx2)*fx1[i];                                                     //non singlet P3qqbV
      //singlet
      if (i > g)
	{
	  tgamma3 += + (sumnegfx2p)*p3sqq2b*(-lx2)*fx1[i];
	  tgamma3 += + (sumposfx2p)*p3sqqb2b*(-lx2)*fx1[i];
	}
      else
	{
	  tgamma3 += + (sumposfx2p)*p3sqq2b*(-lx2)*fx1[i];
	  tgamma3 += + (sumnegfx2p)*p3sqqb2b*(-lx2)*fx1[i];
	}
      //qg Channel
      tgamma3+=+fx2p[g]*p3qg2b*(-lx2)*fx1[i];

      // ----------------------------------------------------------------------------------------------------------------------------------------
      //sig21[sp] = tgamma3;
      //if (order == 2) continue;	  
      // ----------------------------------------------------------------------------------------------------------------------------------------
      // C1caC1cbarb and tH2st
      // ----------------------------------------------------------------------------------------------------------------------------------------	   				
      // Leg 1
      cacheuv_(cz1a);
      if (cz1a >= 1) continue;
      C1cbarbpiuC1ca += +(fx1p[i]*cqq_(cz1a)+fx1p[g]*cqg_(cz1a))*(-lx1)*fx2[j];
      tH2st += +fx1p[i]*cqq_(cz1a)*(-lx1)*fx2[j]*C1qqdelta*0.0;// Since in the Hard Scheme C1qqdelta is ZERO
      tH2st += +fx1p[g]*cqg_(cz1a)*(-lx1)*fx2[j]*C1qqdelta*0.0;// Since in the Hard Scheme C1qqdelta is ZERO
      tH2st += +0.5*(-lx1*(fx1p[i]-fx1[i]*cz1a)*H2qqD0/(1.-cz1a))*fx2[j];
      tH2st += +fx1p[i]*c2qqreghschm_(cz1a)*(-lx1)*fx2[j] ;
      tH2st += +fx1p[g]*c2qghschm_(cz1a)*(-lx1)*fx2[j] ;
      tH2st += +fx1p[im]*c2qqb_(cz1a)*(-lx1)*fx2[j] ;
      tH2st += +fx1pafx1pmanoti[i]*c2qqp_(cz1a)*(-lx1)*fx2[j] ;
      tH2st += -0.5*H2qqD0*d0intx1*fx1[i]*fx2[j];
      // Leg 2
      cacheuv_(cz2b);
      if (cz2b >= 1) continue;
      C1cbarbpiuC1ca += +(fx2p[j]*cqq_(cz2b)+fx2p[g]*cqg_(cz2b))*(-lx2)*fx1[i];
      tH2st += +fx2p[j]*cqq_(cz2b)*(-lx2)*fx1[i]*C1qqdelta*0.0 ;// Since in the Hard Scheme C1qqdelta is ZERO
      tH2st += +fx2p[g]*cqg_(cz2b)*(-lx2)*fx1[i]*C1qqdelta*0.0 ;// Since in the Hard Scheme C1qqdelta is ZERO
      tH2st += +0.5*(-lx2*(fx2p[j]-fx2[j]*cz2b)*H2qqD0/(1.-cz2b))*fx1[i];
      tH2st += +fx2p[j]*c2qqreghschm_(cz2b)*(-lx2)*fx1[i]; 
      tH2st += +fx2p[g]*c2qghschm_(cz2b)*(-lx2)*fx1[i] ;
      tH2st += +fx2p[jm]*c2qqb_(cz2b)*(-lx2)*fx1[i] ;
      tH2st += +fx2pafx2pmanotj[j]*c2qqp_(cz2b)*(-lx2)*fx1[i];
      tH2st += -0.5*H2qqD0*d0intx2*fx1[i]*fx2[j];

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
      if (cz1a >= 1) continue;
      tcgah += ( fx1p[i] *cqqpqq_(cz1a)+ fx1p[g] *(cqqpqg_(cz1a)+cqgpgg_(cz1a)))*(-lx1)*fx2[j];
      tcgah += ( sumfx1p )*cqgpgq_(cz1a)*(-lx1)*fx2[j];
	      
      if (cz2b >= 1) continue;
      tcgah += ( fx2p[j] *cqqpqq_(cz2b)+ fx2p[g] *(cqqpqg_(cz2b)+cqgpgg_(cz2b)))*(-lx2)*fx1[i];
      tcgah += ( sumfx2p )*cqgpgq_(cz2b)*(-lx2)*fx1[i];

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
      double sig1 = (sig12[sp]*LL2*real(mesq::mesqij[sp])+sig11[sp]*LL1*real(mesq::mesqij[sp]))*asopi + (sig12_qed[sp]*LL2*real(mesq::mesqij[sp])+sig11_qed[sp]*LL1*real(mesq::mesqij[sp]))*resconst::a0qed/M_PI;
      if (order < 3 || fpc) xmsq += -sig1;
      //if (isnan_ofast(xmsq))
      //cout << "mesq " << mesq::mesqij[sp] << " sig12 " << sig12[sp] << " LL2 " << LL2 << endl;
      if (order == 1) continue;
      
      //(as/pi)^2 factor
      double sig2 = ((sig24[sp]*LL4*real(mesq::mesqij[sp])+sig23[sp]*LL3*real(mesq::mesqij[sp])+sig22[sp]*LL2*real(mesq::mesqij[sp]))+sig21[sp]*LL1*real(mesq::mesqij[sp]))*pow(asopi,2);
      // sig2 =  sig21[sp]*qtint::LL3_real(mesq::mesqij[sp])*pow(asopi,3);
      if (order < 3 || fpc) xmsq += -sig2;
      if (order == 2) continue;
      
      //(as/pi)^3 factor
      double sig3 = ((sig36[sp]*LL6*real(mesq::mesqij[sp])+sig35[sp]*LL5*real(mesq::mesqij[sp])+sig34[sp]*LL4*real(mesq::mesqij[sp])+sig33[sp]*LL3*real(mesq::mesqij[sp]))+sig32[sp]*LL2*real(mesq::mesqij[sp])+(sig31[sp]*LL1*real(mesq::mesqij[sp])))*pow(asopi,3);
      xmsq += -sig3;
    }	  

  // Flux for Born cross section
  //double fluxborn = fbGeV2/(2*q2);
  //double Vol = 1.;
  //double shad = pow(opts.sqrts,2);

  
  double shad = pow(opts.sroot,2);
  xmsq = xmsq/shad * 3./8.;
  xmsq = xmsq * m2;
  if (isnan_ofast(xmsq))
    cout << m << " " << y << " " << phasespace::costh << "  " << xmsq << endl;
  
  //apply resummation switching
  double swtch = switching::swtch(qt, m);
  //if (swtch < 0.01) return 0.;// do not apply this cut to avoid discontinuities. Instead the phase space is limited to qt and m switching limits
  xmsq = xmsq*swtch;

  xmsq = xmsq*jac;//jacobian for the change of variable qt=qtp/sqrt(1-qtp^2/m^2)

  xmsq *= fpc_weight;
  
  //f[0] = xmsq*w; //vegas importance sampling weight
  f[0] = xmsq; //vegas importance sampling weight
}
