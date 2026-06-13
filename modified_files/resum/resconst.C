#include "resconst.h"
#include "constants.h"
#include "settings.h"
#include "accel.h"
#include "mesq.h"
#include "phasespace.h"
#include <math.h>
#include <iostream>

using namespace std;
using namespace constants;
using namespace mesq;
using namespace phasespace;

// constants.f
const int resconst::NF = 5;
//number of quarks up, down and leptons
const int resconst::nups = 2;
const int resconst::ndowns = 3;
const int resconst::nleps = 3;
//square of the quark charge
const double resconst::equ2 = pow(eqU,2);
const double resconst::eqd2 = pow(eqD,2);
// Number of flavors weighed with charges
const double resconst::NFq2 = 3.*(resconst::nups*resconst::equ2+resconst::ndowns*resconst::eqd2)+resconst::nleps;
const double resconst::NFq3 = 0.;// Quitamos la contribucion, es 0 por simetria F-antiF
const double resconst::NFq4 = 3.*(resconst::nups*pow(resconst::equ2,2)+resconst::ndowns*pow(resconst::eqd2,2))+resconst::nleps;
const double resconst::NFq6 = 3.*(resconst::nups*pow(resconst::equ2,3)+resconst::ndowns*pow(resconst::eqd2,3))+resconst::nleps;  
// Value of alpha_qed(0)
//const double resconst::aaqed = 1./137.035999084; // PDG 2022 value
const double resconst::a0qed = 7.2973525693e-3;  // CODATA 2018 value (https://inspirehep.net/literature/1898737)
// Value of pcF coefficient (Hst) DY QED
const double resconst::pcF_QED = 1.;

//DYRes precision
//const double resconst::Euler = 0.57721566;
//const double resconst::Z2 = 1.644934;
//const double resconst::Z3 = 1.202057;

//higher precision
const double resconst::Z2 = zeta2; //M_PI*M_PI/6.;
const double resconst::Z3 = zeta3; //1.20205690315959428540;
const double resconst::Z4 = zeta4; //1.08232323371113819152;
const double resconst::Z5 = zeta5; //1.03692775514336992633;
const double resconst::b0 = 2.*exp(-euler);

double resconst::beta0, resconst::beta1, resconst::beta2, resconst::beta3, resconst::beta4, resconst::beta5;
double resconst::beta0p, resconst::beta1p, resconst::beta01_QED, resconst::beta01_QCD;
double resconst::Kappa;
double resconst::A1g, resconst::A2g, resconst::A3g, resconst::A4g, resconst::B1g, resconst::B2g, resconst::B3g, resconst::C1ggn;
double resconst::A1q, resconst::A2q, resconst::A3q, resconst:: A4q, resconst:: A5q, resconst:: A6q;
double resconst::A1qp, resconst::A2qp;
double resconst::B1q, resconst::B2q, resconst::B3q, resconst::B4q, resconst::B5q;
double resconst::B1qp, resconst::B2qp, resconst::D1w_QED;
double resconst::B1qHSCHM, resconst::B2qHSCHM, resconst::B3qHSCHM;
double resconst::H1q, resconst::H2q, resconst::H3q, resconst::H4q, resconst::H1g, resconst::H2g;
double resconst::H1qp;
double resconst::K1, resconst::K2, resconst::K3, resconst::K4;
double resconst::C1ggdelta;
double resconst::C1qqn, resconst::C2qqn;
double resconst::C1qqnp;
double resconst::C1qqdelta, resconst::Delta2qq;
double resconst::D0qqqq, resconst::D1qqqq, resconst::D0gggg, resconst::D1gggg;
double resconst::Deltaqqqq, resconst::Deltagggg;
double resconst::H2qqdelta;
double resconst::H2qqD0;
double resconst::s1,resconst::s2,resconst::s3, resconst::s4;
double resconst::k1,resconst::k2,resconst::k3,resconst::k4,resconst::k3_nfv,resconst::k4_nfv;
double resconst::Delta2gg;

void resconst::init()
{
  double nf2 = NF*NF;
  double nf3 = nf2*NF;

  //Resummation coefficients
  beta0=(33.-2.*NF)/12.;
  beta1=(153.-19.*NF)/24.;
  beta2=2857./128.-5033.*NF/1152.+325.*nf2/3456.;
  beta3=(149753./6. + 3564.*Z3 + NF*(-1078361./162.-6508./27.*Z3)+nf2*(50065./162.+6472./81.*Z3)+1093./729.*nf3)/256.; //from https://arxiv.org/pdf/1701.01404.pdf Eq.(3.6), with a 4^4=256 normalisation factor
  beta4=(8157455./16. + 621885./2.*Z3 - 88209./2.*Z4 - 288090.*Z5
	 +NF* (-336460813./1944. - 4811164./81.*Z3 + 33935./6.*Z4 + 1358995./27.*Z5)
	 +NF*NF*(25960913./1944. + 698531./81.*Z3 - 10526./9.*Z4 - 381760./81.*Z5)
	 +NF*NF*NF*(-630559./5832. - 48722./243.*Z3 + 1618./27.*Z4 + 460./9.*Z5)
	 +NF*NF*NF*NF*(1205./2916. - 152./81.*Z3))/1024.;  //from https://arxiv.org/pdf/1701.01404.pdf Eq.(3.7), with a 4^5=1024 normalisation factor
  beta5 = 0.;
  
  //cout << beta0/M_PI << "  " << beta1/pow(M_PI,2) << "  " << beta2/pow(M_PI,3) << "  " << beta3/pow(M_PI,4) << "  " << beta4/pow(M_PI,5) << endl;

  Kappa=67./6.-pi2/2.-5./9.*NF;

  //gluon coefficients
  A1g=CA;
  A2g=CA/2.*(67./6.-pi2/2.-5./9.*NF); //CA/2.*Kappa;
  bool NEWA3 = true;
  if (NEWA3)
    A3g = -14.544163962735272;
  else
  A3g=CA*(13.81-2.15*NF-nf2/108.);
  B1g=-(11.*CA-2.*NF)/6.;
  B2g=CA*CA*(23./24.+(11.*pi2)/18.-3.*Z3/2.)+CF*NF/2.-CA*NF*(1./12.+pi2/9.)-11./8.*CF*CA;
  C1ggn=(pi2/2.+11./2.+pi2)/2.;

  //quark coefficients
  A1q=CF; //4./3;
  A2q=CF/2.*Kappa; //(67./6.-pi2/2.-5./9.*NF);
  //A3q=CF*(13.81-2.15*NF-nf2/108.)+CF*(CA*(29.9259-28.*Z3)-8.2963*NF/2.)*2.*(beta0*4.)/64.; //A3 from Becher & Neubert
  //cout << "A3q orig " << A3q << endl;
  //Better approximation for A3q:
  //A3q=CF*(13.768339123-2.15*NF-nf2/108.)+CF*(CA*(29.9259-28.*Z3)-8.2963*NF/2.)*2.*(beta0*4.)/64.;
  A3q=CF*( (2*CA*CA*(245./24.-67./9.*Z2+11./6.*Z3+11./5.*Z2*Z2)/8.)
	  +(2*(CF*(-55./24.+2*Z3)+CA*(-209./108.+10./9.*Z2-7./3.*Z3))/8.)*NF-nf2/108.)
    +CF*(CA*((808./27.)-28.*Z3)-(224./54.*2.)*NF/2.)*2.*(beta0*4.)/64.; //A3 from Becher & Neubert (https://arxiv.org/pdf/1007.4005.pdf Eq.(74))
  //cout << "A3q fix " << A3q << endl;

  //A4q from https://arxiv.org/pdf/1705.09127.pdf
  //A4q=(CA*CA*CA*CF*(121./3.*Z3*Z2 - 8789./162.*Z2 - 19093./54.*Z3 - 847./24.*Z4 + 132.*Z5 + 3761815./11664.)
  //     + CA*CA*CF*NF*(-22./3.*Z3*Z2 + 2731./162.*Z2 + 4955./54.*Z3 + 11./6.*Z4 -24*Z5 - 31186./243.)
  //     + CA*CF*CF*NF*(272./9.*Z3 + 11.*Z4 - 7351./144.)
  //     + CA*CF*nf2*(-103./81.*Z2 - 47./27.*Z3 +5./6.*Z4 + 13819./972.)
  //     + CF*CF*nf2*(-38./9.*Z3 - 2*Z4 + 215./24.)
  //     + CF*NF*NF*NF*(-4./9.*Z3-232./729.))/16.;
  //A4q=(-1189.19)/16.; //simplied formula for NF=5

  //A4q from https://arxiv.org/pdf/1805.05916.pdf Eq. (B.11)
  //A4q=(CA*CA*CA*CF*(121./3.*Z3*Z2 - 8789./162.*Z2 - 19093./54.*Z3 - 847./24.*Z4 + 132.*Z5 + 3761815./11664.)
  //     + CA*CA*CF*NF*(-22./3.*Z3*Z2 + 2731./162.*Z2 + 4955./54.*Z3 + 11./6.*Z4 -24*Z5 - 31186./243.)
  //     + CA*CF*CF*NF*(272./9.*Z3 + 11.*Z4 - 7351./144.)
  //     + CA*CF*nf2*(-103./81.*Z2 - 47./27.*Z3 +5./4.*Z4 + 13819./972.) //-->likely a typo 5./4.*Z4 -> 5./6.*Z4
  //     + CF*CF*nf2*(-38./9.*Z3 - 2*Z4 + 215./24.)
  //     + CF*NF*NF*NF*(-4./9.*Z3-232./729.))/16.+A4q_th;
  
  //from https://arxiv.org/pdf/1604.01404.pdf
  double gamma1r = CF*CA*(28*Z3-808./27.) + 112.*CF*NF/27.;
  double gamma2r = CF*CA*CA*(-176./3*Z3*Z2 + 6392./81.*Z2 + 12328./27.*Z3 +154./3.*Z4 - 192.*Z5 - 297029./729.)
    + CF*CA*NF* ( -824./81.*Z2 - 904./27.*Z3 +20./3.*Z4 + 62626./729.)
    + CF*nf2* ( -32./9.*Z3 - 1856./729.)
    + CF*CF*NF*(-304./9.*Z3-16.*Z4 + 1711./27.);

  double D2q = -gamma1r;
  double D3q = -gamma2r;

  double A4q_th;
  //A4q_th = 1553./256. //Pade approximant from https://arxiv.org/pdf/hep-ph/0506288.pdf
  //A4q_th = 0.42441*0.0133*pow(M_PI,4); //from Eq. (4.3) of https://arxiv.org/pdf/1805.09638.pdf

  //Exact result for A4q_th from Eq.(3.6) of https://arxiv.org/pdf/1912.12920.pdf (also in https://arxiv.org/pdf/2002.04617.pdf Eq.(6))
  double DfRAnc = 5./2.;
  double DfRRnc = 5./36.;
  A4q_th =  (CF*CA*CA*CA*(84278./81. - 88400./81.*zeta2 + 20944./27.*zeta3 + 1804.*zeta4 - 352./3.*zeta2*zeta3 - 3608./9.*zeta5 - 16.*pow(zeta3,2) - 2504./3.*zeta6)
	     + DfRAnc*(- 128.*zeta2 + 128./3.*zeta3 + 3520./3.*zeta5 - 384.*pow(zeta3,2) - 992.*zeta6)
	     + CF*CF*CF*NF*(572./9. + 592./3.*zeta3 - 320.*zeta5)
	     + CF*CF*CA*NF*(- 34066./81. + 440./3.*zeta2 + 3712./9.*zeta3 - 176.*zeta4 - 128.*zeta2*zeta3 + 160.*zeta5)
	     + CF*CA*CA*NF*(- 24137./81. + 20320./81.*zeta2 - 23104./27.*zeta3 - 176./3.*zeta4 + 448./3.*zeta2*zeta3 + 2096./9.*zeta5)
	     + NF*DfRRnc*(256.*zeta2 - 256./3.*zeta3 - 1280./3.*zeta5)
	     + CF*CF*NF*NF*(2392./81. - 640./9.*zeta3 + 32.*zeta4)
	     + CF*CA*NF*NF*(923./81. - 608./81.*zeta2 + 2240./27.*zeta3 - 112./3.*zeta4)
	     - CF*NF*NF*NF*(32./81. - 64./27.*zeta3))/pow(4.,4);
  
  A4q = A4q_th
    + (3*beta0*4*D3q + 2*beta1*16*D2q)/256.; //formula from Eq. (71) of https://arxiv.org/pdf/1007.4005.pdf with 4^4 normalisation.

  //A5q = 0.;
  if (opts.a5lc)
    //A5 with 5-loop cusp in leading colour approximation from Eq. (17) of https://arxiv.org/abs/1812.11818 
    A5q =  (- 8.29032530100E-1 + opts.A5*0.1*0.42441) * pow(M_PI,5);
  else
    //A5 with 5-loop cusp from Eq. (13) of https://arxiv.org/abs/1812.11818 
    A5q =  (- 7.44150530100E-1 + opts.A5*0.4*0.42441) * pow(M_PI,5);

  A6q = 0.;
  
  B1q=-(3.*CF)/2.;                                                                                //-2;
  B2q=CF*CF*(pi2/4.-3./16.-3.*Z3)+CA*CF*(11.*pi2/36.-193./48.+3.*Z3/2.)+CF*NF*(17./24.-pi2/18.);  //B2q=4./9.*(pi2-3./4.-12.*Z3)+(11./9.*pi2-193./12.+6.*Z3)+NF/6.*(17./3.-4./9.*pi2);
  //B2q=CF*CF*(3./2.*Z2-3./16.-3.*Z3)+CA*CF*(11.*Z2/6.-193./48.+3.*Z3/2.)+CF*NF*(17./24.-Z2/3.);
  //B2q = (13.3447 + 3.4138*NF)/16.; from https://arxiv.org/pdf/1604.01404.pdf, Eq.(13)
  //B3q=-116.685/8. // from https://arxiv.org/pdf/1705.09127.pdf, but does not look OK
  B3q=(7358.86 -721.516*NF + 20.5951* nf2)/64.; // from https://arxiv.org/pdf/1604.01404.pdf, Eq.(13)

  double B3g = B3q;

  //B4q = 0.;
  //B4q = B3q*M_PI;
  //B4q = 6.36907E+0*pow(M_PI,4);
  //B4q = (6.3690704E+0 + opts.B4*2.7500000e-05) *pow(M_PI,4);
  //updated value using results from https://arxiv.org/pdf/2505.09381
  //B4q = (6.36907201841E+0 + opts.B4*1.6040600e-06) *pow(M_PI,4);

  //updated value using results from https://arxiv.org/pdf/2604.09534
  B4q = 6.36907166043E+0 *pow(M_PI,4);

  //NF dependent result
  //B4q = (2.44039557744E+1 +opts.B4*1.604e-06 - 4.91009901731E+0*NF + 2.86101205171E-1*pow(NF,2) - 5.09535039039E-3*pow(NF,3))*pow(M_PI,4);
  
  B5q = 0.;
  
  C1qqn=CF/2.*(pi2/2.-4.);          // Only delta(1-z) part, i.e. N independent part

  // Only delta(1-z) part, i.e. N independent part
//  C2qqn = 1./2.* (CA*CF*(59.*Z3/18. - 1535./192. + 215.*pi2/216. - pi4/240.)
//		  + CF*NF*(192.*Z3 + 1143. - 152.*pi2)/864.
//		  + CF*CF*(-15.*Z3 + 511./16. - 67.*pi2/12. + 17.*pi4/45.)/4. - CF*CF* pow(pi2-8.,2)/16.);

  C2qqn = 1./2.* (CA*CF*(59.*Z3/18. - 1535./192. + 215*pi2/216. - pi4/240.)
		  + CF*NF*(192.*Z3 + 1143. - 152.*pi2)/864.
		  + CF*CF*(-15.*Z3 + 511./16. - 67.*pi2/12. + 17.*pi4/45.)/4. - CF*CF*pow(pi2-8.,2)/16.);
  
  // Delta term in c1qq coefficient
  C1qqdelta=(pi2-8.)/3.;
  C1ggdelta=(11+3*pi2)/4.;


  // Beta functions QED
  beta0p = -4.*NFq2/12.; // in DYqt it was defined as beta0QED
  beta1p = -4.*NFq4/16.; // in DYqt it was defined as beta1QED
  // Beta Mixed functions: New!!!                                    
  beta01_QCD = -1./8.*(nups*equ2+ndowns*eqd2); // in DYqt it was defined as beta01QCD
  beta01_QED = -(nups*equ2+ndowns*eqd2);// in DYqt it was defined as beta01QED
  // QED resummation coefficients
  // Coeficientes usados para la resumacion en QED: Los coeficienten tienen la carga factorizada
  // hay que multiplicar por eq^2 para recuperar la expresion correcta.
  A1qp = 1.;
  A2qp = 1./2.*(-10./9.*NFq2);
  B1qp = -3./2.;
  // Only for W or massive charged final states
  D1w_QED  = -1./2.;
  C1qqnp = 1./2.*(pi2/2.-4.); // Only delta(1-z) part, i.e. N independent part                                                                                                               


  // Delta term in P2qq splitting function (as/M_PI normalization)
  Delta2qq=16./9.*(3./8.-pi2/2.+6.*Z3)+4.*(17./24.+11.*pi2/18.-3.*Z3)-2./3.*NF*(1./6.+2.*pi2/9.);
  Delta2qq=Delta2qq/4.;

  // Delta term in P2gg splitting function (as/pi normalization)
  Delta2gg=(9.*(8./3.+3.*Z3)-2./3.*NF-2*NF)/4.;

  // Coefficients of D0 and D1 in P*P (as/M_PI normalization)
  D0qqqq=8./3.;
  D1qqqq=32./9.;
  D0gggg=6*beta0;
  D1gggg=18.;

  // Coefficients of delta(1-z) in P*P
  Deltaqqqq=4./9.*(9./4.-2.*pi2/3.);
  Deltagggg=pow(beta0, 2)-3./2.*pi2;

  // H2qq contribution: coefficient of delta(1-z)
  H2qqdelta=-2561./144.+127.*NF/72.+3.*pi2/2.-19.*NF*pi2/81.+49.*pi4/324.+58.*Z3/9.+8.*NF*Z3/27.;

  // H2qq contribution: coefficient of D0(z)
  H2qqD0=-404./27.+(56.*NF)/81.+14.*Z3;
  //(23) of https://arxiv.org/pdf/1209.0158.pdf (1/1-z)+)

  //DY scheme (opts.resscheme == 0)
  H1q = 0.;
  H1qp = 0.;   // QED
  H2q = 0.;
  H3q = 0.;
  H4q = 0.;
  H1g = 0.;
  H2g = 0.;

  //to define resummation schemes, see: https://arxiv.org/pdf/hep-ph/0008184.pdf, https://arxiv.org/pdf/1209.0158.pdf

  int NF2 = NF*NF;
  //sn are the constants of the threshold soft function, see Eq. (16S) of https://arxiv.org/pdf/1604.01404.pdf (-> not true)
  s1 = -CF*pi2/12.;
  s2 = CF*(9.*CF*pi4 + CA*(4856. - 603.*pi2 + 18.*pi4 - 2772.*zeta3) + NF*(-656. + 90.*pi2 + 504.*zeta3))/2592.;
  s3 = CF*(pow(CA,2)*(689661*pi4 - 55548*pi6 + 105*pi2*(-297481 + 89100*zeta3) + 35*(5211949 - 8161128*zeta3 + 1353024*pow(zeta3,2) + 2630232*zeta5))
		  - 7*(405*pow(CF,2)*pi6 + 4*NF2*(640 + 13770*pi2 + 297*pi4 + 37800*zeta3) - 27*CF*NF*(6*pi4 - 5*pi2*(-3131 + 2664*zeta3) + 5*(-42727 + 20928*zeta3 + 12096*zeta5)))
		  - 7*CA*(135*CF*pi2*(4856 - 603*pi2 + 18*pi4 - 2772*zeta3) + 2*NF*(5616*pi4 - 75*pi2*(7453 + 324*zeta3) + 5*(412765 - 660312*zeta3 + 402408*zeta5))))/29393280.;

  double ss1 = - 1./2.*CF*zeta2;
  double ss2 = - 41./162.*CF*NF + 7./36.*CF*NF*zeta3 + 5./24.*CF*NF*zeta2 + 607./324.*CA*CF + 5./8.*CA*CF*zeta4 - 77./72.*CA*CF*zeta3 - 67./48.*CA*CF*zeta2;
  double ss3 = - 4./6561.*CF*pow(NF,2) - 11./432.*CF*pow(NF,2)*zeta4 - 35./972.*CF*pow(NF,2)*zeta3 - 17./216.*CF*pow(NF,2)*zeta2 - 42727./31104.*pow(CF,2)*NF + 7./18.*pow(CF,2)*NF*zeta5 + 19./72.*pow(CF,2)*NF*zeta4 + 109./162.*pow(CF,2)*NF*zeta3 + 275./576.*pow(CF,2)*NF*zeta2 - 5./12.*pow(CF,2)*NF*zeta2*zeta3 - 412765./419904.*CA*CF*NF - 23./24.*CA*CF*NF*zeta5 - 13./54.*CA*CF*NF*zeta4 + 1019./648.*CA*CF*NF*zeta3 + 37265./23328.*CA*CF*NF*zeta2 + 5./72.*CA*CF*NF*zeta2*zeta3 + 5211949./839808.*pow(CA,2)*CF + 451./144.*pow(CA,2)*CF*zeta5 + 3649./1728.*pow(CA,2)*CF*zeta4 - 37783./3888.*pow(CA,2)*CF*zeta3 + 29./18.*pow(CA,2)*CF*pow(zeta3,2) - 297481./46656.*pow(CA,2)*CF*zeta2 - 1543./1512.*pow(CA,2)*CF*zeta2*zeta4 + 275./144.*pow(CA,2)*CF*zeta2*zeta3;
  double ss4 = (1.+opts.H4)*accel::levin(s1,s2,s3);
 
  s1 = ss1;
  s2 = ss2+ss1*ss1/2.;
  s3 = ss3+ss1*ss2+pow(ss1,3)/6.;
  s4 = ss4+ss1*ss3+pow(ss2,2)/2.+pow(ss1,2)*ss2/2.+pow(ss1,4)/24.;
  
  //kn are the enq of Eq.(72) of https://arxiv.org/pdf/1007.4005.pdf
  k1 = CF*(-4 + 7*pi2/12.);
  k2 = CF*(27.*CF*(7665. - 1660.*pi2 + 134.*pi4 - 3600.*zeta3) + 10.*NF*(4085. - 546.*pi2 + 72.*zeta3) + CA*(-255785. + 31830.*pi2 - 288.*pi4 + 112680.*zeta3))/25920.;
  k3 = CF*(NF2*(-190931./419904. + (403*pi2)/3888. + (43*pi4)/38880. - (13*zeta3)/243.)
		  + CA*CF*(824281./20736. - (406507*pi2)/62208. + (92237*pi4)/155520. - (739*pi6)/54432. - (13141*zeta3)/432. + (845*pi2*zeta3)/288. + (37*pow(zeta3,2))/12. - (689*zeta5)/72.)
		  + pow(CA,2)*(-51082685./1679616. + (596513*pi2)/139968. - (4303*pi4)/311040. + (299*pi6)/102060. + (505087*zeta3)/15552. - (73*pi2*zeta3)/36. - (71*pow(zeta3,2))/18. - (217*zeta5)/144.)
		  + CF*NF*(-56963./31104. + (13705*pi2)/15552. - (1463*pi4)/15552. +  (815*zeta3)/162. - (37*pi2*zeta3)/144. - (13*zeta5)/9.)
		  + CA*NF*(1700171./209952. - (201749*pi2)/139968. -  (35*pi4)/15552. - (134*zeta3)/27. + (37*pi2*zeta3)/144. -  zeta5/24.)
		  + pow(CF,2)*(-5599/384. + (4339*pi2)/2304. - (173*pi4)/480. +  (27403*pi6)/1088640. - (115*zeta3)/16. - (35*pi2*zeta3)/48. + pow(zeta3,2)/2. + (83*zeta5)/4.));
  k3_nfv = ( (-4 + pow(NC,2)) *(1./8. + (5*pi2)/96. - pi4/2880. + (7*zeta3)/48. - (5*zeta5)/6.) ) /NC; 
  k4 = 321482899./26873856.+41601./(32768.*pow(NC,4))+21812995./(2654208.*pow(NC,2))+(583291805.*pow(NC,2))/71663616.-(6361535363.*pow(NC,4))/214990848.-(5528861.*NF)/(3981312.*pow(NC,3))-(677109931.*NF)/(214990848.*NC)-(3301359761.*NC*NF)/214990848.+(237612677.*pow(NC,3)*NF)/11943936.+(818632639.*pow(NF,2))/214990848.-(8087615.*pow(NF,2))/(26873856.*pow(NC,2))-(9307799.*pow(NC,2)*pow(NF,2))/2654208.-(1865531.*pow(NF,3))/(13436928.*NC)+(1865531.*NC*pow(NF,3))/13436928.+(73.*0.)/(72.*pow(NC,3))+(617.*0.)/(288.*NC)-(2309.*NC*0.)/576.+(491.*pow(NC,3)*0.)/576.+(245.*NF*0.)/288.-(49.*NF*0.)/(72.*pow(NC,2))-49./288.*pow(NC,2)*NF*0.-(282994049.*pi2)/80621568.+(1411.*pi2)/(49152.*pow(NC,4))-(190811.*pi2)/(497664.*pow(NC,2))-(151186117.*pow(NC,2)*pi2)/161243136.+(1548736387.*pow(NC,4)*pi2)/322486272.+(2613235.*NF*pi2)/(5971968.*pow(NC,3))+(255090349.*NF*pi2)/(322486272.*NC)+(373605629.*NC*NF*pi2)/161243136.-(1143416297.*pow(NC,3)*NF*pi2)/322486272.-(112557385.*pow(NF,2)*pi2)/161243136.-(1064527.*pow(NF,2)*pi2)/(40310784.*pow(NC,2))+(116815493.*pow(NC,2)*pow(NF,2)*pi2)/161243136.+(19847.*pow(NF,3)*pi2)/(559872.*NC)-(19847.*NC*pow(NF,3)*pi2)/559872.+(25.*0.*pi2)/(96.*pow(NC,3))+(2893.*0.*pi2)/(3456.*NC)-(9595.*NC*0.*pi2)/6912.+(2009.*pow(NC,3)*0.*pi2)/6912.+(695.*NF*0.*pi2)/3456.-(139.*NF*0.*pi2)/(864.*pow(NC,2))-(139.*pow(NC,2)*NF*0.*pi2)/3456.-(5021281.*pi4)/107495424.+(77933.*pi4)/(4423680.*pow(NC,4))+(656111.*pi4)/(4976640.*pow(NC,2))-(15813319.*pow(NC,2)*pi4)/26873856.+(104417575.*pow(NC,4)*pi4)/214990848.-(457717.*NF*pi4)/(9953280.*pow(NC,3))-(14829269.*NF*pi4)/(268738560.*NC)+(91615573.*NC*NF*pi4)/268738560.-(12885589.*pow(NC,3)*NF*pi4)/53747712.-(29591.*pow(NF,2)*pi4)/746496.+(4783.*pow(NF,2)*pi4)/(466560.*pow(NC,2))+(109691.*pow(NC,2)*pow(NF,2)*pi4)/3732480.+(949.*pow(NF,3)*pi4)/(933120.*NC)-(949.*NC*pow(NF,3)*pi4)/933120.-(41.*0.*pi4)/(1920.*pow(NC,3))-(2.*0.*pi4)/(405.*NC)+(28367.*NC*0.*pi4)/829440.-(6559.*pow(NC,3)*0.*pi4)/829440.-(95.*NF*0.*pi4)/20736.+(19.*NF*0.*pi4)/(5184.*pow(NC,2))+(19.*pow(NC,2)*NF*0.*pi4)/20736.+(1099393.*pi6)/29859840.-(22679.*pi6)/(8709120.*pow(NC,4))-(722717.*pi6)/(59719680.*pow(NC,2))-(11045201.*pow(NC,2)*pi6)/418037760.+(25733.*pow(NC,4)*pi6)/5971968.+(90193.*NF*pi6)/(52254720.*pow(NC,3))-(75449.*NF*pi6)/(11612160.*NC)+(253727.*NC*NF*pi6)/52254720.-(419.*pow(NC,3)*NF*pi6)/4976640.-(1751.*pow(NF,2)*pi6)/26127360.+(1217.*pow(NF,2)*pi6)/(8709120.*pow(NC,2))-(95.*pow(NC,2)*pow(NF,2)*pi6)/1306368.+(943.*0.*pi6)/(725760.*pow(NC,3))-(41449.*0.*pi6)/(5225472.*NC)+(428999.*NC*0.*pi6)/52254720.-(16481.*pow(NC,3)*0.*pi6)/10450944.-(11.*NF*0.*pi6)/11664.+(11.*NF*0.*pi6)/(14580.*pow(NC,2))+(11.*pow(NC,2)*NF*0.*pi6)/58320.+(6530779.*pi8)/6270566400.+(5926681.*pi8)/(20901888000.*pow(NC,4))-(743663.*pi8)/(870912000.*pow(NC,2))-(19452931.*pow(NC,2)*pi8)/31352832000.+(267479.*pow(NC,4)*pi8)/1791590400.+23./256.*zeta53-(69.*zeta53)/(320.*pow(NC,4))+(1277.*zeta53)/(1280.*pow(NC,2))-61./80.*pow(NC,2)*zeta53-7./64.*pow(NC,4)*zeta53-(211844861.*zeta3)/2985984.-(23965.*zeta3)/(6144.*pow(NC,4))+(65653.*zeta3)/(6144.*pow(NC,2))+(11207041.*pow(NC,2)*zeta3)/373248.+(33976055.*pow(NC,4)*zeta3)/995328.-(820819.*NF*zeta3)/(331776.*pow(NC,3))+(3797689.*NF*zeta3)/(186624.*NC)-(17374699.*NC*NF*zeta3)/2985984.-(2000053.*pow(NC,3)*NF*zeta3)/165888.+(35681.*pow(NF,2)*zeta3)/46656.-(1005107.*pow(NF,2)*zeta3)/(746496.*pow(NC,2))+(144737.*pow(NC,2)*pow(NF,2)*zeta3)/248832.-(13.*pow(NF,3)*zeta3)/(7776.*NC)+(13.*NC*pow(NF,3)*zeta3)/7776.+(19.*0.*zeta3)/(16.*pow(NC,3))+(529.*0.*zeta3)/(96.*NC)-6487./768.*NC*0.*zeta3+1343./768.*pow(NC,3)*0.*zeta3+35./64.*NF*0.*zeta3-(7.*NF*0.*zeta3)/(16.*pow(NC,2))-7./64.*pow(NC,2)*NF*0.*zeta3+(1177831.*pi2*zeta3)/746496.-(3895.*pi2*zeta3)/(6144.*pow(NC,4))+(1058855.*pi2*zeta3)/(331776.*pow(NC,2))-(1343051.*pow(NC,2)*pi2*zeta3)/2985984.-(5502499.*pow(NC,4)*pi2*zeta3)/1492992.-(132493.*NF*pi2*zeta3)/(497664.*pow(NC,3))-(567893.*NF*pi2*zeta3)/(497664.*NC)-(51719.*NC*NF*pi2*zeta3)/497664.+(752105.*pow(NC,3)*NF*pi2*zeta3)/497664.+(19.*pow(NF,2)*pi2*zeta3)/1458.+(21509.*pow(NF,2)*pi2*zeta3)/(186624.*pow(NC,2))-(23941.*pow(NC,2)*pow(NF,2)*pi2*zeta3)/186624.-(pow(NF,3)*pi2*zeta3)/(3456.*NC)+(NC*pow(NF,3)*pi2*zeta3)/3456.+(19.*0.*pi2*zeta3)/(288.*pow(NC,3))+(91.*0.*pi2*zeta3)/(144.*NC)-(4039.*NC*0.*pi2*zeta3)/4608.+(823.*pow(NC,3)*0.*pi2*zeta3)/4608.+35./576.*NF*0.*pi2*zeta3-(7.*NF*0.*pi2*zeta3)/(144.*pow(NC,2))-7./576.*pow(NC,2)*NF*0.*pi2*zeta3+(82177.*pi4*zeta3)/829440.+(863.*pi4*zeta3)/(46080.*pow(NC,4))-(61981.*pi4*zeta3)/(414720.*pow(NC,2))+71./576.*pow(NC,2)*pi4*zeta3-(75989.*pow(NC,4)*pi4*zeta3)/829440.+(6083.*NF*pi4*zeta3)/(207360.*pow(NC,3))-(12167.*NF*pi4*zeta3)/(414720.*NC)-(1403.*NC*NF*pi4*zeta3)/138240.+(421.*pow(NC,3)*NF*pi4*zeta3)/41472.-(11.*0.*pi4*zeta3)/(1440.*pow(NC,3))+(5.*0.*pi4*zeta3)/(576.*NC)-(19.*NC*0.*pi4*zeta3)/23040.-(pow(NC,3)*0.*pi4*zeta3)/4608.+(746603.*pow(zeta3,2))/82944.+(2697.*pow(zeta3,2))/(1024.*pow(NC,4))-(677.*pow(zeta3,2))/(1536.*pow(NC,2))-(70517.*pow(NC,2)*pow(zeta3,2))/41472.-(196867.*pow(NC,4)*pow(zeta3,2))/20736.+(113.*NF*pow(zeta3,2))/(128.*pow(NC,3))-(26489.*NF*pow(zeta3,2))/(10368.*NC)+(9517.*NC*NF*pow(zeta3,2))/41472.+(59827.*pow(NC,3)*NF*pow(zeta3,2))/41472.-(673.*pow(NF,2)*pow(zeta3,2))/5184.+(2279.*pow(NF,2)*pow(zeta3,2))/(20736.*pow(NC,2))+(413.*pow(NC,2)*pow(NF,2)*pow(zeta3,2))/20736.+(21.*0.*pow(zeta3,2))/(32.*pow(NC,3))-(2755.*0.*pow(zeta3,2))/(384.*NC)+3113./384.*NC*0.*pow(zeta3,2)-305./192.*pow(NC,3)*0.*pow(zeta3,2)-55./48.*NF*0.*pow(zeta3,2)+(11.*NF*0.*pow(zeta3,2))/(12.*pow(NC,2))+11./48.*pow(NC,2)*NF*0.*pow(zeta3,2)+(1289.*pi2*pow(zeta3,2))/6912.+(131.*pi2*pow(zeta3,2))/(1536.*pow(NC,4))-(25.*pi2*pow(zeta3,2))/(288.*pow(NC,2))+(1387.*pow(NC,2)*pi2*pow(zeta3,2))/13824.-(493.*pow(NC,4)*pi2*pow(zeta3,2))/1728.+(1092277.*zeta5)/18432.+(29041.*zeta5)/(6144.*pow(NC,4))-(1110893.*zeta5)/(36864.*pow(NC,2))-(15801697.*pow(NC,2)*zeta5)/552960.-(729227.*pow(NC,4)*zeta5)/138240.+(3461.*NF*zeta5)/(18432.*pow(NC,3))-(235601.*NF*zeta5)/(55296.*NC)+(2457541.*NC*NF*zeta5)/276480.-(49313.*pow(NC,3)*NF*zeta5)/10240.-(20927.*pow(NF,2)*zeta5)/17280.+(949.*pow(NF,2)*zeta5)/(3456.*pow(NC,2))+899./960.*pow(NC,2)*pow(NF,2)*zeta5+(53.*pow(NF,3)*zeta5)/(17280.*NC)-(53.*NC*pow(NF,3)*zeta5)/17280.-(955.*0.*zeta5)/(192.*pow(NC,3))+(19915.*0.*zeta5)/(2304.*NC)-(19705.*NC*0.*zeta5)/4608.+(2795.*pow(NC,3)*0.*zeta5)/4608.-275./288.*NF*0.*zeta5+(55.*NF*0.*zeta5)/(72.*pow(NC,2))+55./288.*pow(NC,2)*NF*0.*zeta5+35./432.*pi2*zeta5+(371.*pi2*zeta5)/(768.*pow(NC,4))-(37715.*pi2*zeta5)/(55296.*pow(NC,2))-(4821.*pow(NC,2)*pi2*zeta5)/2048.+(68345.*pow(NC,4)*pi2*zeta5)/27648.+(517.*NF*pi2*zeta5)/(13824.*pow(NC,3))-(NF*pi2*zeta5)/(9.*NC)+(857.*NC*NF*pi2*zeta5)/1728.-(5837.*pow(NC,3)*NF*pi2*zeta5)/13824.+(55.*0.*pi2*zeta5)/(144.*pow(NC,3))-(35.*0.*pi2*zeta5)/(24.*NC)+1015./768.*NC*0.*pi2*zeta5-(565.*pow(NC,3)*0.*pi2*zeta5)/2304.+149./256.*zeta3*zeta5-(7.*zeta3*zeta5)/(8.*pow(NC,4))-(7.*zeta3*zeta5)/(64.*pow(NC,2))-1429./256.*pow(NC,2)*zeta3*zeta5+383./64.*pow(NC,4)*zeta3*zeta5-(40251.*zeta7)/2048.-(29919.*zeta7)/(4096.*pow(NC,4))+(243851.*zeta7)/(24576.*pow(NC,2))+(1329235.*pow(NC,2)*zeta7)/73728.-(36605.*pow(NC,4)*zeta7)/36864.+(1707.*NF*zeta7)/(2048.*pow(NC,3))-(35.*NF*zeta7)/(512.*NC)-(19265.*NC*NF*zeta7)/9216.+(24427.*pow(NC,3)*NF*zeta7)/18432.-(721.*0.*zeta7)/(192.*pow(NC,3))-(497.*0.*zeta7)/(768.*NC)+(8813.*NC*0.*zeta7)/1536.-(2051.*pow(NC,3)*0.*zeta7)/1536.;
  k4_nfv = (73.*1.)/(72.*pow(NC,3))+(617.*1.)/(288.*NC)-(2309.*NC*1.)/576.+(491.*pow(NC,3)*1.)/576.+(245.*NF*1.)/288.-(49.*NF*1.)/(72.*pow(NC,2))-49./288.*pow(NC,2)*NF*1.+(25.*1.*pi2)/(96.*pow(NC,3))+(2893.*1.*pi2)/(3456.*NC)-(9595.*NC*1.*pi2)/6912.+(2009.*pow(NC,3)*1.*pi2)/6912.+(695.*NF*1.*pi2)/3456.-(139.*NF*1.*pi2)/(864.*pow(NC,2))-(139.*pow(NC,2)*NF*1.*pi2)/3456.-(41.*1.*pi4)/(1920.*pow(NC,3))-(2.*1.*pi4)/(405.*NC)+(28367.*NC*1.*pi4)/829440.-(6559.*pow(NC,3)*1.*pi4)/829440.-(95.*NF*1.*pi4)/20736.+(19.*NF*1.*pi4)/(5184.*pow(NC,2))+(19.*pow(NC,2)*NF*1.*pi4)/20736.+(943.*1.*pi6)/(725760.*pow(NC,3))-(41449.*1.*pi6)/(5225472.*NC)+(428999.*NC*1.*pi6)/52254720.-(16481.*pow(NC,3)*1.*pi6)/10450944.-(11.*NF*1.*pi6)/11664.+(11.*NF*1.*pi6)/(14580.*pow(NC,2))+(11.*pow(NC,2)*NF*1.*pi6)/58320.+(19.*1.*zeta3)/(16.*pow(NC,3))+(529.*1.*zeta3)/(96.*NC)-6487./768.*NC*1.*zeta3+1343./768.*pow(NC,3)*1.*zeta3+35./64.*NF*1.*zeta3-(7.*NF*1.*zeta3)/(16.*pow(NC,2))-7./64.*pow(NC,2)*NF*1.*zeta3+(19.*1.*pi2*zeta3)/(288.*pow(NC,3))+(91.*1.*pi2*zeta3)/(144.*NC)-(4039.*NC*1.*pi2*zeta3)/4608.+(823.*pow(NC,3)*1.*pi2*zeta3)/4608.+35./576.*NF*1.*pi2*zeta3-(7.*NF*1.*pi2*zeta3)/(144.*pow(NC,2))-7./576.*pow(NC,2)*NF*1.*pi2*zeta3-(11.*1.*pi4*zeta3)/(1440.*pow(NC,3))+(5.*1.*pi4*zeta3)/(576.*NC)-(19.*NC*1.*pi4*zeta3)/23040.-(pow(NC,3)*1.*pi4*zeta3)/4608.+(21.*1.*pow(zeta3,2))/(32.*pow(NC,3))-(2755.*1.*pow(zeta3,2))/(384.*NC)+3113./384.*NC*1.*pow(zeta3,2)-305./192.*pow(NC,3)*1.*pow(zeta3,2)-55./48.*NF*1.*pow(zeta3,2)+(11.*NF*1.*pow(zeta3,2))/(12.*pow(NC,2))+11./48.*pow(NC,2)*NF*1.*pow(zeta3,2)-(955.*1.*zeta5)/(192.*pow(NC,3))+(19915.*1.*zeta5)/(2304.*NC)-(19705.*NC*1.*zeta5)/4608.+(2795.*pow(NC,3)*1.*zeta5)/4608.-275./288.*NF*1.*zeta5+(55.*NF*1.*zeta5)/(72.*pow(NC,2))+55./288.*pow(NC,2)*NF*1.*zeta5+(55.*1.*pi2*zeta5)/(144.*pow(NC,3))-(35.*1.*pi2*zeta5)/(24.*NC)+1015./768.*NC*1.*pi2*zeta5-(565.*pow(NC,3)*1.*pi2*zeta5)/2304.-(721.*1.*zeta7)/(192.*pow(NC,3))-(497.*1.*zeta7)/(768.*NC)+(8813.*NC*1.*zeta7)/1536.-(2051.*pow(NC,3)*1.*zeta7)/1536.;
  
  //Hard scheme
  if (opts.resscheme == 1)
    {
      H1q = k1+s1;             //H1qdelta  = 2*C1qq_delta = 2.*(pi2-8.)/3.
      H1qp = H1q*(3./4.);       // H1qp(QED) = H1q/CF
      H2q = k2+s2+k1*s1;       //H2qqdelta = 2*C2qq_delta + C1qq_delta^2 = -2561./144.+ 127.*NF/72. + 3.*pi2/2. -19.*NF*pi2/81.+ 49.*pi4/324. + 58.*zeta3/9. + 8.*NF*zeta3/27.
      H3q = k3+s3+k2*s1+k1*s2; //H3qqdelta = 2*C3qq_delta + 2*C1qq_delta*C2qq_delta;
      H4q = k4+s4+k2*s2+k1*s3+k3*s1;
      H1g = CA*pi2/2. + 11./2.; // c_H in the limit m_t \to \infty is soupposed to be 11/2: 	arXiv:1311.1654 eq.86
      double costh1 = phasespace::costh;
      H2g = pow(CA,2)*(3187./288. + 7./8.*mesq::LQ(costh1) + 157./72.*pi2 + 13./144.*pi4 - 55./18.*zeta3) + CA*CF*(-145./24. - 11./8.*mesq::LQ(costh1) - 3./4.*pi2) + 9./4.*pow(CF,2) - 5./96.*CA - 1./12.*CF - CA*NF*(287./144. + 5./36.*pi2 + 4./9.*zeta3) + CF*NF*(-41./24. + 1./2.*mesq::LQ(costh1) + zeta3);
    }

  //Semi Hard scheme
  if (opts.resscheme == 2)
    {
      H1q = (k1+s1)/2;              //C1qq_delta
      H1qp = H1q*(3./4.);           // H1qp(QED) = H1q/CF 
      H2q = (K2+s2+k1*s1)/2.;       //C2qq_delta + C1qq_delta^2/2
      H3q = (k3+s3+k2*s1+k1*s2)/2.; //C3qq_delta + C1qq_delta*C2qq_delta;
      H4q = (k4+s4+k2*s2+k1*s3+k3*s1)/2.;
    }
  //CSS scheme (see https://inspirehep.net/files/b825e7a2984b74530665e61ec61c588d 3.69 3.70)
  if (opts.resscheme == 3)
    {
      H1q = CF*(pi2/2.-4.); //Is there a typo in Eq. 3.69 of https://inspirehep.net/files/b825e7a2984b74530665e61ec61c588d? pi -> pi^2
      H1qp = H1q*(3./4.);   //H1qp(QED) = H1q/CF 
      H2q = CF*CA*(59./18.*zeta3-1535./192.+215./216.*pi2-pi4/240.)
	+1./4.*CF*CF*(-15.*zeta3+511./16.-67./12.*pi2+17./45.*pi4)
	+1./864.*CF*NF*(192.*zeta3+1143.-152.*pi2);
      H3q = 0.;
    }
  //NS scheme
  if (opts.resscheme == 4)
    {
      H1q = CF*(pi2/2.-7./2.); // Eq. (35) of https://arxiv.org/pdf/hep-ph/0008184.pdf
      H1qp = (pi2/2.-7./2.); // H1qp(QED) = H1q/CF
      H2q = 0.;
      H3q = 0.;
    }

  //TMD scheme
  if (opts.resscheme == 5)
    {
      double c1qd = (-CF*pi2/6.)/4.; // From Eq. (7.2) of https://arxiv.org/pdf/1604.07869.pdf
      H1q = 2*(C1qqn-c1qd); //CF*(pi2*7./12.-4.);
      H1qp = H1q*(3./4.);   //H1qp(QED) = H1q/CF 
      //H1q = CF*(pi2*7./12.-4.);
      double c2qd = (CF*CF*pi4/72. + CF*CA*(1214./81.-67.*pi2/36.-77./9.*zeta3+pi4/18.)+CF*TF*NF*(-328/81.+pi2*5./9.+28./9.*zeta3))/16.;
      //c2qd = C2qqdelta - H1q/2.*(-CF*pi2/24) - (H2q/2.-pow(H1q,2)/8.);
      H2q  = 2.*(C2qqn - c2qd - H1q/2.*c1qd +pow(H1q,2)/8.);
      H3q = 0.;
    }
  //cout << H1q << endl;
  //cout << H2q << endl;
  //cout << H3q << endl;

  //see 3.75 3.76 of https://inspirehep.net/files/b825e7a2984b74530665e61ec61c588d (https://inspirehep.net/literature/1839718)
  B2q = B2q - beta0*H1q;
  B3q = B3q - beta1*H1q - beta0*(2.*H2q - pow(H1q,2));
  B4q = B4q - beta2*H1q - beta1*(2.*H2q - pow(H1q,2)) - beta0*(3.*H3q - 3.*H2q*H1q + pow(H1q,3));

  B1qHSCHM = B1q;  
  B2qHSCHM = (1.0/48.0)*CF*(2.0*NF*(1.0 + 4.0*zeta2) + 9.0*CF*(-1.0 + 8.0*zeta2 - 16.0*zeta3) + 3.0*(-17.0 - 44.0*zeta2 + 72.0*zeta3));
  B3qHSCHM = (-4623*CF)/32. 
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
  
  //CMW K coefficients
  if (opts.cmw || opts.cmwfo)
    {
      K1 = A2q/A1q;
      K2 = A3q/A1q;
      K3 = A4q/A1q;
      K4 = A5q/A1q;
    }
  else
    {
      K1 = 0.;
      K2 = 0.;
      K3 = 0.;
      K4 = 0.;
    }
  //A1q=0; //<--
  //A2q=0; //<--
  //A3q=0;
  //A4q=0;
  //A5q=0;
  //B1q=0; //<--
  //B2q=0;
  //B3q=0;
  //B4q=0;
  //C1qqn=0;
}

void resconst::initcmw()
{
  //cout
  //  << "B msbar  " << B1q/pi
  //  << "  " << B2q/pi2
  //  << "  " << B3q/pi3
  //  << "  " << B4q/pi4
  //  << endl;
  init();
  if (opts.cmw)
    {
      double A2qCMW = A2q - K1*A1q;
      double A3qCMW = A3q - K2*A1q;
      double A4qCMW = A4q - K3*A1q;
      double A5qCMW = A5q - K4*A1q;

      A2q = A2qCMW;
      A3q = A3qCMW;
      A4q = A4qCMW;
      A4q = A4qCMW;
      A5q = A4qCMW;
      
      //double B2qCMW = B2q - K1*B1q;
      //double B3qCMW = B3q - K2*B1q;
      //double B4qCMW = B4q - K3*B1q;

      double B2qCMW = B2q - K1*B1q;
      double B3qCMW = B3q - 2.*K1*B2q + 2.*pow(K1,2)*B1q - K2*B1q;
      double B4qCMW = B4q - 3.*B3q*K1 + 5.*B2q*pow(K1,2) - 5.*B1q*pow(K1,3) - 2.*B2q*K2 + 5.*B1q*K1*K2 - B1q*K3;

      B2q = B2qCMW;
      B3q = B3qCMW;
      B4q = B4qCMW;
      
      double beta2CMW = beta2 + beta0*(K2-pow(K1,2))-beta1*K1;
      double beta3CMW = beta3 + 4*beta0*pow(K1,3)-6*beta0*K1*K2+2*beta0*K3+beta1*pow(K1,2)-2*beta2*K1;
      double beta4CMW = beta4 - 14*beta0*pow(K1,4) + 28*beta0*pow(K1,2)*K2 - 12*beta0*K1*K3 - 5*beta0*pow(K2,2) +
	3*beta0*K4 - 2*beta1*K1*K2 + beta1*K3 + 4*beta2*pow(K1,2) - beta2*K2 -
	3*beta3*K1;

      beta2 = beta2CMW;
      beta3 = beta3CMW;
      beta4 = beta4CMW;

      double H2qCMW = H2q - K1*H1q;
      double H3qCMW = H3q - 2.*K1*H2q + 2.*pow(K1,2)*H1q - K2*H1q;

      H2q = H2qCMW;
      H3q = H3qCMW;
    }
  //cout
  //  << "B CMW  " << B1q/pi 
  //  << "  " << B2q/pi2
  //  << "  " << B3q/pi3
  //  << "  " << B4q/pi4
  //  << endl;
}
