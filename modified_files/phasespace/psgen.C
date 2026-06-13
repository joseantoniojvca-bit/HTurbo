#include "phasespace.h"
#include "settings.h"
#include "switch.h"
#include "isnan.h"
#include "cubacall.h"

#include <iostream>
#include <math.h>

//additional useful variables
double phasespace::m2;
double phasespace::qt2;
double phasespace::mt2;
double phasespace::mt;
double phasespace::exppy;
double phasespace::expmy;

//Bjorken-x
double phasespace::x1;
double phasespace::x2;

//Vector boson 4-momentum and boost
double phasespace::pV[4];
double phasespace::gam;
double phasespace::beta[3];

//leptons 4-momenta
double phasespace::p4[4];
double phasespace::p3[4];
double phasespace::p1[4];
double phasespace::p2[4];
double phasespace::p5[4];

//exponent for the change of variable in x2
const double esp=8.;

bool phasespace::gen_m(double x, double& jac, double mlim, bool qtcut, bool qtswitching)
{
  //Generate the boson invariant mass between the integration boundaries,
  //jac gets multiplied by the Jacobian of the change of variable from the unitary interval to the m boundaries
  double qtcutlim = (qtcut && opts.xqtcut != 0.) ? phasespace::qtmax/opts.xqtcut : phasespace::mmax;
  double switchlim = (qtswitching) ? switching::mlimit(phasespace::qtmin) : phasespace::mmin;
  double wsqmin = pow(max(phasespace::mmin,switchlim),2);
  double wsqmax = pow(min(mlim,min(phasespace::mmax,qtcutlim)),2);
  if (wsqmin >= wsqmax)
    return false;
  //If the mass is generated at the phase space limit, there is no room for integration in y (i.e. ymx-ymn = 0)
  if (wsqmax == pow(mlim,2) && x == 1)
    return false;
  mweight_breitw_(x,wsqmin,wsqmax,opts.rmass,opts.rwidth,m2,jac);
  m = sqrt(m2);
  return true;
}

bool phasespace::gen_ms(double x, double xs, double& jac, double mlim, bool qtcut, bool qtswitching)
{
  //Generate the boson invariant mass between the integration boundaries,
  //jac gets multiplied by the Jacobian of the change of variable from the unitary interval to the m boundaries
  double qtcutlim = (qtcut && opts.xqtcut != 0.) ? phasespace::qtmax/opts.xqtcut : phasespace::mmax;
  double switchlim = (qtswitching) ? switching::mlimit(phasespace::qtmin) : phasespace::mmin;
  double wsqmin = pow(max(phasespace::mmin,switchlim),2);
  double wsqmax = pow(min(mlim,min(phasespace::mmax,qtcutlim)),2);
  if (wsqmin >= wsqmax)
    return false;
  //If the mass is generated at the phase space limit, there is no room for integration in y (i.e. ymx-ymn = 0)
  if (wsqmax == pow(mlim,2) && x == 1)
    return false;

  if(opts.calc_type != "OLD")
  { // New implementation LEAN -- two photon decay
    mweight_flat_(x,wsqmin,wsqmax,m2,jac);
    m = sqrt(m2);
    //C----------------------------------------------------------------------C
    //C Suavizado gaussiano
    //c     xx(8) --> smearing
    double     SIGMA = opts.SigmaSmearing; // Reccomended value 1.7GeV
    double     Qsolo=m;
    double     GaussMin = atan(-Qsolo/SIGMA);
    double     GaussMax = atan(Qsolo/SIGMA);
    double     Gauss = GaussMin +xs *(GaussMax-GaussMin);
    double     JacGauss = GaussMax-GaussMin;
    double     Qmedido = Qsolo + SIGMA*tan(Gauss);
    JacGauss=JacGauss*SIGMA/(cos(Gauss)*cos(Gauss));
    // if(m < 125.)
    // 	{

    // 	  cout << " Q = " << m << endl;
    // 	  cout << " Qmedido = " << Qmedido << endl;
    // 	  cout << " JacGauss = " << JacGauss << endl;
    // 	}

    double     FactorM=exp(-(Qmedido-Qsolo)*(Qmedido-Qsolo)/(2.0*SIGMA*SIGMA))/(SIGMA*sqrt(2.0*M_PI))*JacGauss;
    m = Qmedido;
    jac *=FactorM;
    //C----------------------------------------------------------------------C
  }

  if(opts.calc_type == "OLD")
  {
    // CJ FO Breit-Wigner - implement NWA for Higgs
    mweight_breitw_(x,wsqmin,wsqmax,opts.rmass,opts.rwidth,m2,jac);
    m = sqrt(m2);
  }
  return true;
}

bool phasespace::gen_mys(const double x[3], double& jac, bool qtcut, bool qtswitching)
{
  //generate phase space as m, y
  //jac gets multiplied by the Jacobian of the change of variables from the unitary square x[2] to the m, y boundaries
  //Current version works for born kinematics, i.e. qt=0 (which is, probably, the correct implementation for born, resummed and counterterm)
  bool status = true;

  //Generate the boson invariant mass between the integration boundaries
  //kinematic limit from the relation x(1,2) < 1 ==> exp(fabs(y)) < sqrt(s)/m
  double ymn;
  if (ymin*ymax <= 0.)
    ymn = 0.;
  else
    ymn = min(fabs(ymin),fabs(ymax));
  double mlim = opts.sroot/exp(ymn);
  status = gen_ms(x[0],x[2], jac, mlim, qtcut, qtswitching);
  if (!status)
    return false;

  //Generate the boson rapidity between the integration boundaries
  //Here there is a delicate point, should consider the value of qt != 0 when setting the kinematical limit for y?
  //For counterterm and resummed piece, the cross section is calculated at born level with pt=0, and the qt is generated independently
  //However, this will determine some inconsistencies at large rapidities, because of nonzero cross sections above the kinematical limit for y

  //double mtmin = sqrt(m2+pow(phasespace::qtmin,2));
  //double tmpx=(m2+pow(opts.sroot,2))/opts.sroot/mtmin;
  //double ylim=log((tmpx+sqrt(pow(tmpx,2)-4.))/2.); //Limit y boundaries to the kinematic limit in y

  double ylim = 0.5*log(pow(opts.sroot,2)/m2); //Limit y boundaries to the kinematic limit in y

  status = gen_y(x[1], jac, ylim);
  if (!status)
    return false;

  return status;
}

bool phasespace::gen_myqts(const double x[4], double& jac, bool qtcut, bool qtswitching)
{
  //With respect to the previous routine, the order in which y and qt are generated is reversed
  //generate phase space as m, y, qt
  //jac gets multiplied by the Jacobian of the change of variables from the unitary cube x[3] to the m, y, qt boundaries
  bool status = true;

  //Generate the boson invariant mass between the integration boundaries
  //kinematic limit from the relation x(1,2) < 1 ==> exp(fabs(y)) < sqrt(s)/m
  //  cout << "ymax = " << ymax << " ; " << " ymin = " << ymin << endl;
  double ymn;
  if (ymin*ymax <= 0.)
    ymn = 0.;
  else
    ymn = min(fabs(ymin),fabs(ymax));
  double mlim = opts.sroot/exp(ymn);
  status = gen_ms(x[0],x[3], jac, mlim, qtcut, qtswitching);
  if (!status)
    return false;

  //Generate the boson rapidity between the integration boundaries
  double ylim;
  if (qtswitching)
    ylim = 0.5*log(pow(opts.sroot,2)/m2); //Limit y boundaries to the kinematic limit in y
    else
    {
      //Limit y boundaries to the kinematic limit in y, accounting for phasespace::qtmin
      double qtmn = qtcut ?
      max(max(opts.qtcut,opts.xqtcut*m), phasespace::qtmin)
      : phasespace::qtmin;
      double mtmin = sqrt(m2+pow(qtmn,2));
      double tmpx=(m2+pow(opts.sroot,2))/opts.sroot/mtmin;
      ylim=log((tmpx+sqrt(max(0.,pow(tmpx,2)-4.)))/2.); //introduced max to avoid neqative argument of sqrt
    }
    status = gen_y(x[1], jac, ylim);
    if (!status)
      return false;
  //  //fodyqt limits
  //  //.....kinematical limits on qt
  //  double z=q2/pow(energy_.sroot_,2);
  //  double xr=pow((1-z),2)-4*z*pow(qt/q,2);
  //  if (xr < 0)
  //    {
  //      f[0] = 0.;
  //      return 0;
  //    }
  //

  //Generate the boson transverse momentum between the integration boundaries
  double qtlim;
  if (qtswitching)
  {
    //in switching mode (ct and res) do not account for y in the qt limit (actually the qtlimit should be 1e10, i.e. unphysical)
    double kinqtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+m2,2)/(4*pow(opts.sroot,2))-m2)); //introduced max to avoid neqative argument of sqrt
    double switchqtlim = switching::qtlimit(m);
    // qtlim = min(kinqtlim, switchqtlim);
    qtlim = switchqtlim; //--> in CT and RES should not have any phase space qt limit
  }
  else
  {
    calcexpy();
    double cosh2y=pow((exppy+expmy)*0.5,2);
    qtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+m2,2)/(4*pow(opts.sroot,2)*cosh2y)-m2)); //introduced max to avoid neqative argument of sqrt when y=ymax
  }
  status = gen_qt(x[2], jac, qtlim, qtcut);
  if (!status)
    return false;

  return status;
}

void phasespace::set_xf_bounds(double &ymn, double &ymx)
{
  //cout << "input set_xf_bounds " << ymn << "  " << ymx << endl;

  double yxfmn, yxfmx;
  double mm;
  if (opts.exactkin)
    mm = sqrt(m*m+qt*qt); //exact kinematics including qt
  else
    mm = m; //BORN kinematic at qt = 0
  
  //xF cuts
  yxfmn = asinh(opts.sroot/mm/2.*opts.xfmin);
  yxfmx = asinh(opts.sroot/mm/2.*opts.xfmax);

  //Equivalent formulas as in Eq. 87 of https://arxiv.org/pdf/0901.0002
  //double Hmin, Hmax;
  //Hmin = opts.sroot/mm/2.*opts.xfmin;
  //Hmax = opts.sroot/mm/2.*opts.xfmax;
  //yxfmn = log(Hmin+sqrt(Hmin*Hmin+1));
  //yxfmx = log(Hmax+sqrt(Hmax*Hmax+1));

  //cout << "yxf set_xf_bounds " << yxfmn << "  " << yxfmx << endl;

  ymn = max(ymn,yxfmn);
  ymx = min(ymx,yxfmx);
  
  //  if (opts.m3dsdm)
  //  {
  //    //additional ycuts from <m>
  //    double ycmn, ycmx;
  //    double mavg = (phasespace::mmin + phasespace::mmax)/2.;
  //    double ptavg = (phasespace::qtmin + phasespace::qtmax)/2.;
  //    double mt = sqrt(mavg*mavg*ptavg*ptavg);
  //    ycmn = asinh(opts.sroot/mavg/2.*opts.xfmin);
  //    ycmx = asinh(opts.sroot/mavg/2.*opts.xfmax);
  //    ymn = max(ymn,ycmn);
  //    ymx = min(ymx,ycmx);
  //  }

  //cout << "output set_xf_bounds " << ymn << "  " << ymx << endl;
}

bool phasespace::gen_y(double x, double& jac, double ylim)
{
  //Generate the boson rapidity between the integration boundaries,
  //jac gets multiplied by the Jacobian of the change of variable from the unitary interval to the y boundaries

  double ymn = min(max(-ylim, phasespace::ymin),ylim);
  double ymx = max(min(ylim, phasespace::ymax),-ylim);

  //xF cuts
  if (opts.dsdxf)
    set_xf_bounds(ymn, ymx);
  
  if (ymn >= ymx)
    return false;

  //Uniform weight
  y = ymn+(ymx-ymn)*x;
  jac = jac*(ymx-ymn);

  ////integration in xf
  //if (opts.dsdxf && !cubacall::isvegas)
  //  {
  //    jac = jac/(opts.xfmax-opts.xfmin);
  //    jac = jac * 2*phasespace::m/opts.sroot * cosh(y);
  //  }

  return true;

  //Unweight assuming all PDFs have the form x^a*(1-x)^b, with a = 1 and b = 1
  // --> s(y) = tau^2*(1+tau^2-2*tau*cosh(y)) with tau = m/sqrt(s)
  // change of variable x = int_y0^y s(y) / int_y0^y1 s(y)
  // --> x = ((y*(1+tau^2)-2*tau*sinh(y))-a)/(b-a)
  // with a = y0*(1+tau^2)-2*tau*sinh(y0);
  // and  b = y1*(1+tau^2)-2*tau*sinh(y1);
  double tau = phasespace::m/opts.sroot;
  double tau2 = pow(tau,2);
  double a = ymn*(1+tau2)-2*tau*sinh(ymn);
  double b = ymx*(1+tau2)-2*tau*sinh(ymx);

  //need to find y by newton method: x_n+1 =  x_n - f(x_n)/f'(x_n)
  //f(y) = ((y*(1+tau^2)-2*tau*sinh(y))-a)/(b-a) - x = 0
  //f'(y) = ((1+tau^2)-2*tau*cosh(y))/(b-a)

  //cout << x << "  ";
  double yn = (x*(b-a)+a)/(1+tau2);
  for (int n = 0; n < 5; n++)
    {
      //cout << yn << "  ";
      double ynp1 = yn - ((yn*(1+tau2)-2*tau*sinh(yn)-a)/(b-a)-x) / (((1+tau2)-2*tau*cosh(yn))/(b-a));
      //ynp1 = min(max(ymn, ynp1),ymx);
      //ynp1 = max(min(ymx, ynp1),ymn);
      yn = ynp1;
    }
  //cout << yn << endl;
  y = yn;
  
  jac = jac*((b-a)/(1+tau2-2*tau*cosh(y)));
  return true;
  
}

bool phasespace::gen_qt(double x, double& jac, double qtlim, bool qtcut)
{
  //Generate the boson rapidity between the integration boundaries,
  //jac gets multiplied by the Jacobian of the change of variable from the unitary interval to the qt boundaries
  double qtmn = qtcut ?
    max(max(opts.qtcut,opts.xqtcut*m), phasespace::qtmin)
    : max(max(opts.qtcutoff,phasespace::qtmin),opts.xqtcutoff*m);
  double qtmx = min(qtlim, phasespace::qtmax);
  //if (opts.makecuts)
  //  if (phasespace::m/2 < opts.lptcut)
  //    qtmn = max(qtmn,2.*opts.lptcut);
  if (qtmn >= qtmx)
    return false;
  if (qtcut) //phase space generation for counterterm and V+j fixed order
    qtweight_(x,qtmn,qtmx,qt,jac);
    //qtweight_lo_(x,qtmn,qtmx,m2,qt,jac);
    //qtweight_flat_(x,qtmn,qtmx,qt,jac);
  else       //phase space generation for resummed cross section
    if (qtmx > 5)
      qtweight_res_(x,qtmn,qtmx,qt,jac);
    else
      qtweight_flat_(x,qtmn,qtmx,qt,jac);
  qt2 = qt*qt;
  return true;
}

bool phasespace::gen_qt_ctfo(double x, double& jac)
{
  //Generate the boson rapidity between the integration boundaries,
  //jac gets multiplied by the Jacobian of the change of variable from the unitary interval to the qt boundaries
  double qtmn = max(opts.qtcut,opts.xqtcut*m);
  double qtmx = 1e10;
  qtweight_(x,qtmn,qtmx,qt,jac);
  //qtweight_flat_(xqt,qtmn,qtmx,qt,jac);
  qt2 = qt*qt;
  return true;
}

bool phasespace::gen_qt_fpc(double x, double& jac)
{
  //Generate the boson pt between the integration boundaries,
  //jac gets multiplied by the Jacobian of the change of variable from the unitary interval to the qt boundaries
  double qtmn = m*opts.qtfpc;
  double qtmx = max(opts.qtcut,opts.xqtcut*phasespace::m);
  qtweight_(x,qtmn,qtmx,qt,jac);
  //qtweight_flat_(xqt,qtmn,qtmx,qt,jac);
  qt2 = qt*qt;
  return true;
}

bool phasespace::gen_mqty(const double x[3], double& jac, bool qtcut, bool qtswitching)
{
  //generate phase space as m, qt, y
  //jac gets multiplied by the Jacobian of the change of variables from the unitary cube x[3] to the m, qt, y boundaries
  bool status = true;

  //Generate the boson invariant mass between the integration boundaries
  //kinematic limit from the relation x(1,2) < 1 ==> exp(fabs(y)) < sqrt(s)/m
  double ymn;
  if (ymin*ymax <= 0.)
    ymn = 0.;
  else
    ymn = min(fabs(ymin),fabs(ymax));
  double mlim = opts.sroot/exp(ymn);
  status = gen_m(x[0], jac, mlim, qtcut, qtswitching);
  if (!status)
    return false;

  //Generate the boson transverse momentum between the integration boundaries
  //double qtlim = sqrt(max(0.,pow(pow(opts.sroot,2)-m2,2))/(4*pow(opts.sroot,2))); //introduced max to avoid negative argument of sqrt when  y=ymax
  //Account for ymin-ymax
  double exppylim = exp(ymn);
  double expmylim = 1./exppylim;
  double cosh2y=pow((exppylim+expmylim)*0.5,2);

  //kinematic limit from the relation E < sqrt(s) ==> (mt + qt)*cosh(y) < sqrt(s)
  double qtlim;
  if (qtswitching)
    {
      double switchqtlim = switching::qtlimit(m);
      //in switching mode (ct and res) do not account for y in the qt limit (actually the qtlimit should be 1e10, i.e. unphysical)
      //double kinqtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+m2,2)/(4*pow(opts.sroot,2))-m2)); //introduced max to avoid neqative argument of sqrt
      double kinqtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+m2,2)/(4*pow(opts.sroot,2)*cosh2y)-m2)); //kinematic limit including rapidity
      if (opts.qtlim) //--> introduce an option to have the kinematic limit in qt applied in the resummed and CT terms
	qtlim = min(kinqtlim, switchqtlim);
      else
	qtlim = switchqtlim; //--> in CT and RES should not have any phase space qt limit
    }
  else
    //qtlim = sqrt(max(0.,pow(pow(opts.sroot,2)-m2*cosh2y,2)/(4*pow(opts.sroot,2)*cosh2y))); //introduced max to avoid negative argument of sqrt when y=ymax --> wrong formula calculated by me
    qtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+m2,2)/(4*pow(opts.sroot,2)*cosh2y)-m2)); //introduced max to avoid neqative argument of sqrt when y=ymax --> apparently correct formula, cannot check
  
  
  status = gen_qt(x[1], jac, qtlim, qtcut);
  if (!status)
    return false;
  
  //Generate the boson rapidity between the integration boundaries
  calcmt();
  double tmpx=(m2+pow(opts.sroot,2))/opts.sroot/mt;
  double ylim=log((tmpx+sqrt(max(0.,pow(tmpx,2)-4.)))/2.); //Limit y boundaries to the kinematic limit in y --> introduced max to avoid neqative argument of sqrt
  status = gen_y(x[2], jac, ylim);
  if (!status)
    return false;

  return status;
}

bool phasespace::gen_myqt(const double x[3], double& jac, bool qtcut, bool qtswitching)
{
  //With respect to the previous routine, the order in which y and qt are generated is reversed
  //generate phase space as m, y, qt
  //jac gets multiplied by the Jacobian of the change of variables from the unitary cube x[3] to the m, y, qt boundaries
  bool status = true;

  //Generate the boson invariant mass between the integration boundaries
  //kinematic limit from the relation x(1,2) < 1 ==> exp(fabs(y)) < sqrt(s)/m
  double ymn;
  if (ymin*ymax <= 0.)
    ymn = 0.;
  else
    ymn = min(fabs(ymin),fabs(ymax));
  double mlim = opts.sroot/exp(ymn);
  status = gen_m(x[0], jac, mlim, qtcut, qtswitching);
  if (!status)
    return false;

  //Generate the boson rapidity between the integration boundaries
  double ylim;
  if (qtswitching)
    ylim = 0.5*log(pow(opts.sroot,2)/m2); //Limit y boundaries to the kinematic limit in y
  else
    {
      //Limit y boundaries to the kinematic limit in y, accounting for phasespace::qtmin  
      double qtmn = qtcut ?
	max(max(opts.qtcut,opts.xqtcut*m), phasespace::qtmin)
	: phasespace::qtmin;
      double mtmin = sqrt(m2+pow(qtmn,2));
      double tmpx=(m2+pow(opts.sroot,2))/opts.sroot/mtmin;
      ylim=log((tmpx+sqrt(max(0.,pow(tmpx,2)-4.)))/2.); //introduced max to avoid neqative argument of sqrt
    }
  status = gen_y(x[1], jac, ylim);
  if (!status)
    return false;
//  //fodyqt limits
//  //.....kinematical limits on qt
//  double z=q2/pow(energy_.sroot_,2);
//  double xr=pow((1-z),2)-4*z*pow(qt/q,2);
//  if (xr < 0)
//    {
//      f[0] = 0.;
//      return 0;
//    }
//  

  //Generate the boson transverse momentum between the integration boundaries
  double qtlim;
  if (qtswitching)
    {
      //in switching mode (ct and res) do not account for y in the qt limit (actually the qtlimit should be 1e10, i.e. unphysical)
      double switchqtlim = switching::qtlimit(m);
      //double kinqtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+m2,2)/(4*pow(opts.sroot,2))-m2)); //introduced max to avoid neqative argument of sqrt
      double cosh2y = pow((exppy+expmy)*0.5,2);
      double kinqtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+m2,2)/(4*pow(opts.sroot,2)*cosh2y)-m2)); //kinematic limit including rapidity
      if (opts.qtlim) //--> introduce an option to have the kinematic limit in qt applied in the resummed and CT terms
	qtlim = min(kinqtlim, switchqtlim);
      else
	qtlim = switchqtlim; //--> in CT and RES should not have any phase space qt limit
    }
  else
    {
      calcexpy();
      double cosh2y = pow((exppy+expmy)*0.5,2);
      qtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+m2,2)/(4*pow(opts.sroot,2)*cosh2y)-m2)); //introduced max to avoid neqative argument of sqrt when y=ymax
    }
  status = gen_qt(x[2], jac, qtlim, qtcut);
  if (!status)
    return false;

  return status;
}

bool phasespace::gen_mqt(const double x[2], double& jac, bool qtcut, bool qtswitching)
{
  //generate phase space as m, qt (for resummed part)
  //jac gets multiplied by the Jacobian of the change of variables from the unitary square x[2] to the m, qt boundaries
  //Current version works for born kinematics, i.e. qt=0 (which is, probably, the correct implementation for born, resummed and counterterm)
  bool status = true;

  //Generate the boson invariant mass between the integration boundaries
  //kinematic limit from the relation x(1,2) < 1 ==> exp(fabs(y)) < sqrt(s)/m
  double ymn;
  if (ymin*ymax <= 0.)
    ymn = 0.;
  else
    ymn = min(fabs(ymin),fabs(ymax));
  double mlim = opts.sroot/exp(ymn);
  status = gen_m(x[0], jac, mlim, qtcut, qtswitching);
  if (!status)
    return false;

  //Generate the boson pt between the integration boundaries
  //Here there is a delicate point, should consider the value of y != 0 when setting the kinematical limit for qt?
  //For counterterm and resummed piece, the cross section is calculated at born level with pt=0, and the qt is generated independently
  //However, this will determine some inconsistencies at large rapidities, because of nonzero cross sections above the kinematical limit for qt
  
  double switchqtlim = switching::qtlimit(m);

  double qtlim;
  //double kinqtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+m2,2)/(4*pow(opts.sroot,2))-m2)); //introduced max to avoid neqative argument of sqrt
  double cosh2y = pow((exppy+expmy)*0.5,2);
  double kinqtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+m2,2)/(4*pow(opts.sroot,2)*cosh2y)-m2)); //kinematic limit including rapidity
  if (opts.qtlim) //--> introduce an option to have the kinematic limit in qt applied in the resummed and CT terms
    qtlim = min(kinqtlim, switchqtlim);
  else
    qtlim = switchqtlim; //--> in CT and RES should not have any phase space qt limit
  
  status = gen_qt(x[1], jac, qtlim, qtcut);
  if (!status)
    return false;

  return status;
}  

bool phasespace::gen_my(const double x[2], double& jac, bool qtcut, bool qtswitching)
{
  //generate phase space as m, y
  //jac gets multiplied by the Jacobian of the change of variables from the unitary square x[2] to the m, y boundaries
  //Current version works for born kinematics, i.e. qt=0 (which is, probably, the correct implementation for born, resummed and counterterm)
  bool status = true;

  //Generate the boson invariant mass between the integration boundaries
  //kinematic limit from the relation x(1,2) < 1 ==> exp(fabs(y)) < sqrt(s)/m
  double ymn;
  if (ymin*ymax <= 0.)
    ymn = 0.;
  else
    ymn = min(fabs(ymin),fabs(ymax));
  double mlim = opts.sroot/exp(ymn);
  status = gen_m(x[0], jac, mlim, qtcut, qtswitching);
  if (!status)
    return false;

  //Generate the boson rapidity between the integration boundaries
  //Here there is a delicate point, should consider the value of qt != 0 when setting the kinematical limit for y?
  //For counterterm and resummed piece, the cross section is calculated at born level with pt=0, and the qt is generated independently
  //However, this will determine some inconsistencies at large rapidities, because of nonzero cross sections above the kinematical limit for y
  
  //double mtmin = sqrt(m2+pow(phasespace::qtmin,2));
  //double tmpx=(m2+pow(opts.sroot,2))/opts.sroot/mtmin;
  //double ylim=log((tmpx+sqrt(pow(tmpx,2)-4.))/2.); //Limit y boundaries to the kinematic limit in y

  double ylim = 0.5*log(pow(opts.sroot,2)/m2); //Limit y boundaries to the kinematic limit in y

  status = gen_y(x[1], jac, ylim);
  if (!status)
    return false;

  return status;
}  

void phasespace::gen_costhphi(const double x[2], double& jac)
{
  //generate phase space as costh, phi
  //jac gets multiplied by the Jacobian of the change of variables from the unitary square x[2] to the costh, phi_lep boundaries
  gen_costh(x[0], jac);
  gen_phi(x[1], jac);  
}

void phasespace::gen_costh(const double x, double& jac)
{
  //generate phase space as costh
  //cosine of the polar angle of the lepton in the boson rest frame
  costh=-1.+2.*x;
  jac=jac*2.;
}

void phasespace::gen_phi(const double x, double& jac)
{
  //generate phase space as phi
  //lepton azimuthal angle in the boson rest frame
  phi_lep = -M_PI+2.*M_PI*x; // --> important to generate phi in [-pi,pi] to use fast trigonometric relations as sin(phi) = sqrt(1 - cos2(phi)) * sign(phi)
  jac=jac*2.*M_PI;
}

bool phasespace::gen_x2(const double x, double& jac)
{
  //Calculate Bjorken x1 and x2
  //Phase space is generated for integration performed in dx2
  double ss = opts.sroot;
  double s = ss*ss;
  double x2min=(m2-ss*mt*expmy)/(ss*mt*exppy-s);
  double lx2min = log(x2min);
    
  if (x2min > 1. || x2min < 0.)
    {
      //cout << "error in x2min " << x2min << " m " << m << " pt " << qt  << " y " << y << endl;
      jac=0.;
      return false;
    }
  
  x2=exp((1.-pow(x,esp))*lx2min);
  jac=jac*fabs(x2*lx2min*esp*pow(x,(esp-1.)));
      
  //define x1 as a function of x2 from energy conservation
  x1=(ss*mt*exppy*x2-m2)/(x2*s-ss*mt*expmy);

  if (x1 > 1. || x2 > 1.)
    {
      jac=0.;
      return false;
    }

  //Jacobian from the change of variable in delta(p1+p2-p3-p4-p5)
  double xjj=fabs(x2*s-ss*mt*expmy);
  jac = jac/xjj;
  return true;
}
