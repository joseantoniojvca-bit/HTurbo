#ifndef phasespace_h
#define phasespace_h

#include "settings.h"

#include <math.h>
#include <algorithm>

//bool naiveRF = false;
//bool CSRF = true;
//bool kt0RF = false;


using namespace std;

//fortran interface
extern "C" {
  void sety_(double &yy);
  void setqt_(double &qtt);
  int binner_(double p3[4], double p4[4]);

  //tools
  void dyboost_(double& gamma, double beta[3], double pin[4], double pout[4]);
  void boostv_(double& m,double p[4],double& gamma,double beta[3]);
  void rotate_(double vin[3], double& c, double& s, double ax[3], double vout[3]);
  void genp_(double& costh, double& phi, double& m, double p[4]);
  void qtweight_(double& x,double& qtmin,double& qtmax,double& qt,double& jac);
  void qtweight_lo_(double& x,double& qtmin,double& qtmax,double& q2,double& qt,double& jac);
  void qtweight_res_(double& x,double& qtmin,double& qtmax,double& qt,double& jac);
  void qt2weight_(double& x,double& qt2min,double& qt2max,double& qt2,double& jac);
  void qtweight_flat_(double& x,double& qtmin,double& qtmax,double& qt,double& jac);
  void mweight_breitw_(double& x,double& mmin2,double& mmax2,double& rmass,double& rwidth,double& m2,double& jac);
  void mweight_flat_(double& x,double& mmin,double& mmax,double& m,double& jac);
}

/** @brief Interface between calculation and user input.
 *
 * This namespace is used to provide values of kinematic parameters to integral
 * boundaries and integrands.
 *
 */
namespace phasespace
{
  enum restframeid {CS=0, naive=1, kt0=2};

  extern restframeid CurrentRF;
#pragma omp threadprivate(CurrentRF)
  
  //  extern restframeid RF;
  //  inline void setRF(restframeid rfid) {RF = rfid;};
  
  //boson integration boundaries
  extern double mmin;
  extern double mmax;
  extern double qtmin;
  extern double qtmax;
  extern double ymin;
  extern double ymax;
  extern double xfmin;
  extern double xfmax;
  extern void setbounds(double m1, double m2, double qt1, double qt2, double y1, double y2);

  extern double cthmin;
  extern double cthmax;
  extern void setcthbounds(double cth1, double cth2);

  //lepton integration boundaries
  extern double lpt2min;
  extern double lpt2max;
  extern double lexp2ymin;
  extern double lexp2ymax;
  extern void setlepbounds(double pt1, double pt2, double y1, double y2);

  extern bool lepbinner(double p3[4], double p4[4]);
  
  //point in phase space
  extern double m, qt, y, phiV;
#pragma omp threadprivate(m, qt, y, phiV)
  extern double costh, phi_lep;
#pragma omp threadprivate(costh,phi_lep)
  extern double x1, x2;
#pragma omp threadprivate(x1,x2)

  //cached values of cos and sin of phi_lep
  extern double cosphi_lep, sinphi_lep;
#pragma omp threadprivate(cosphi_lep,sinphi_lep)
  
  //sign flip the costh boundaries according to y
  //inline double getcthmin() {return cthmin;};
  //inline double getcthmax() {return cthmax;};
  inline double getcthmin()
  {
    double val = min(cthmin * (y < 0 ? -1 : 1),cthmax * (y < 0 ? -1 : 1));
    ////In W production, an mT cut corresponds to an upper cut on |costh|
    //if (opts.makecuts && opts.nproc != 3 && opts.mtcut > 0)
    //  {
    //	double cthlim = sqrt(1.-pow(opts.mtcut,2)/pow(m,2));
    //	val = max(val,-cthlim);
    //	val = min(val,cthlim);
    //  }
    return val;
  };
  inline double getcthmax()
  {
    double val = max(cthmin * (y < 0 ? -1 : 1),cthmax * (y < 0 ? -1 : 1));
    ////In W production, an mT cut corresponds to an upper cut on |costh|
    //if (opts.makecuts && opts.nproc != 3 && opts.mtcut > 0)
    //  {
    //	double cthlim = sqrt(1.-pow(opts.mtcut,2)/pow(m,2));
    //	val = max(val,-cthlim);
    //	val = min(val,cthlim);
    //  }
    return val;
  };
  
  //auxiliary useful quantities
  extern double m2, qt2;
  extern double mt, mt2;
  extern double exppy, expmy;
#pragma omp threadprivate(m2, qt2, mt, mt2, exppy, expmy)
    
  extern void set_mqtyphi(double M, double Qt, double Y, double PhiV = 0.);
  extern void set_m(double M);
  extern void set_qt(double Qt);
  extern void set_y(double Y);
  extern void set_phiV(double PhiV);
  
  inline void set_cth(double Costh) {costh = Costh;};
  inline void set_philep(double Phi_lep) {phi_lep = Phi_lep;};

  inline void calcexpy() {exppy = exp(phasespace::y); expmy=1./exppy;};
  inline void calcmt() {mt2 = m2+qt2; mt = sqrt(mt2);};

  inline void calcphilep() {cosphi_lep = cos(phi_lep-M_PI); sinphi_lep = sin(phi_lep-M_PI);};
  //inline void calcphilep() {cosphi_lep = cos(phi_lep-M_PI); sinphi_lep = sqrt(max(0.,1.-pow(cosphi_lep,2)))*((phi_lep-M_PI)>0 ? 1 : -1);}; //-> this is not necessarily faster, as cos and sin are probably calculated together with optimisation flags

  extern void set_xf_bounds(double &ymn, double &ymx);
  
  //generation of phase space variables from unitary (hyper)cubes
  extern bool gen_m(double x, double& jac, double mlim, bool qtcut = false, bool qtswitching = false);
  extern bool gen_ms(double x, double xs, double& jac, double mlim, bool qtcut = false, bool qtswitching = false);
  extern bool gen_y(double x, double& jac, double ylim);
  extern bool gen_qt(double x, double& jac, double qtlim, bool qtcut = false);
  extern bool gen_qt_ctfo(double x, double& jac);
  extern bool gen_qt_fpc(double x, double& jac);
  
  extern bool gen_mqty(const double x[3], double& jac, bool qtcut = false, bool qtswitching = false);
  extern bool gen_myqt(const double x[3], double& jac, bool qtcut = false, bool qtswitching = false);
  extern bool gen_myqts(const double x[4], double& jac, bool qtcut = false, bool qtswitching = false);
  extern bool gen_mqt(const double x[2], double& jac, bool qtcut = false, bool qtswitching = false);
  extern bool gen_my(const double x[2], double& jac, bool qtcut = false, bool qtswitching = false);
  extern bool gen_mys(const double x[3], double& jac, bool qtcut = false, bool qtswitching = false);
  extern void gen_costhphi(const double x[2], double& jac);
  extern void gen_costh(const double x, double& jac);
  extern void gen_phi(const double x, double& jac);
  extern bool gen_x2(const double x, double& jac);

  //generation of particles
  //Vector boson 4-momentum and boost
  extern double pV[4];
  extern double gam;
  extern double beta[3];
#pragma omp threadprivate(pV,gam,beta)

  extern void genV4p();
  
  //leptons 4-momenta
  extern double p4[4];
  extern double p3[4];
#pragma omp threadprivate(p4,p3)

  inline void getp3(double p[4])
  {
    p[0] = p3[0];
    p[1] = p3[1];
    p[2] = p3[2];
    p[3] = p3[3];
  }
  inline void getp4(double p[4])
  {
    p[0] = p4[0];
    p[1] = p4[1];
    p[2] = p4[2];
    p[3] = p4[3];
  }

  //rest frame axes
  extern double kap1[4];
  extern double xax[3];
  extern double yax[3];
  extern double zax[3];
#pragma omp threadprivate(kap1,xax,yax,zax)

  extern void genRFaxes(restframeid RF);
  
  extern void genl4p();

  extern double p1[4];
  extern double p2[4];
#pragma omp threadprivate(p1,p2)
  extern void genp12();

  extern double p5[4];
#pragma omp threadprivate(p5)
  extern void genp5();

  inline void rotate(double vin[3], double c, double s, double ax[3], double vout[])
  {
    //Rotate the 3-vector "vin" around the axis "ax"
    //by and angle phi with "c" = cos(phi) and "s" = sin(phi)
    //"vout" is the result
    vout[0]=(c+ax[0]*ax[0]*(1-c))      *vin[0] + (ax[0]*ax[1]*(1-c)-ax[2]*s)*vin[1] + (ax[0]*ax[2]*(1-c)+ax[1]*s)*vin[2];
    vout[1]=(ax[1]*ax[0]*(1-c)+ax[2]*s)*vin[0] + (c+ax[1]*ax[1]*(1-c))      *vin[1] + (ax[1]*ax[2]*(1-c)-ax[0]*s)*vin[2];
    vout[2]=(ax[2]*ax[0]*(1-c)-ax[1]*s)*vin[0] + (ax[2]*ax[1]*(1-c)+ax[0]*s)*vin[1] + (c+ax[2]*ax[2]*(1-c))      *vin[2];
  }

  inline void boost(double gamma, double beta[3], double pin[4], double pout[])
  {
    //boost 4-momentum pin into 4-momentum pout
    double bdotp=pin[0]*beta[0]+pin[1]*beta[1]+pin[2]*beta[2];
    pout[3]=gamma*(pin[3]-bdotp);
    pout[0]=pin[0]+gamma*beta[0]*(gamma/(gamma+1.)*bdotp-pin[3]);
    pout[1]=pin[1]+gamma*beta[1]*(gamma/(gamma+1.)*bdotp-pin[3]);
    pout[2]=pin[2]+gamma*beta[2]*(gamma/(gamma+1.)*bdotp-pin[3]);
  }

  //cached phi_lep rotation
  extern double rot1[3];
  extern double r2ax[3];
#pragma omp threadprivate(rot1,r2ax)
  
  //perform the phi_lep rotation
  inline void genl4p_phirot()
  {
    //Start from the x axis, and rotate by angle phi_lep with respect to z axis
    double c, s;
    c = phasespace::cosphi_lep;
    s = phasespace::sinphi_lep;
    phasespace::rotate(xax, c, s, zax, rot1);
    
    //define the axis r2ax = rot1 x zax for the second rotation
    r2ax[0] = rot1[1]*zax[2] - rot1[2]*zax[1];
    r2ax[1] = rot1[2]*zax[0] - rot1[0]*zax[2];
    r2ax[2] = rot1[0]*zax[1] - rot1[1]*zax[0];
  }
  
  //version of genl4p in which the phi_lep rotation is cached
  inline void genl4p_phifix()
  {
    //Generate the 4-momentum of the first lepton in the boson rest frame, using costh and phi_lep
    double p3cm[4];
  
    if (qt == 0. || CurrentRF == naive)
      //if (opts.fixedorder || CurrentRF == naive) //can use opts.fixedorder instead of qt == 0.?
      genp_(costh, phi_lep, m, p3cm);
    else
      {
	//rotate by angle theta with respect to r2ax axis
	double rot2[3];
	double c, s;
	s = phasespace::costh;
	c = sqrt(max(0.,1.-pow(s,2)));
	phasespace::rotate(rot1, c, s, r2ax, rot2);

	p3cm[3]=phasespace::m/2.; //E
	p3cm[0]=p3cm[3]*rot2[0];  //px
	p3cm[1]=p3cm[3]*rot2[1];  //py
	p3cm[2]=p3cm[3]*rot2[2];  //pz
      }
    //Boost to go in the lab frame
    //dyboost_(gam, beta, p3cm, p3);
    phasespace::boost(gam, beta, p3cm, p3);
  
    //momentum of the second lepton
    p4[0]=pV[0]-p3[0];      //px
    p4[1]=pV[1]-p3[1];      //py
    p4[2]=pV[2]-p3[2];      //pz
    p4[3]=pV[3]-p3[3];      //E
    //ensure energy is calculated with enough precision
    //p4[3]=sqrt(p4[0]*p4[0]+p4[1]*p4[1]+p4[2]*p4[2]);
  }
}

#endif
