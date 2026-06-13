//   Standard includes
#include <algorithm>
//declaracion a archivos lean*.h

const double Nc  = 3.0;
const double Na  = 8.0;
const double Nc2 = 9.0;
const double Na2 = 64.0;
const double Cf  = 1.33333333333333333333;
const double Ca  = 3.0;
const double Gq  = 2.0;

#include <cmath>
#include <quadmath.h>
//necesarios para setw
#include <iomanip>
//   nlo includes
#include "psg-phasespace_n2i2f0.h"
#include "Inter_params.h"
// DYturbo includes
#include "settings.h"
#include "coupling.h"
#include "phasespace.h"
#include "switch.h"

#define __DEB_1_8PI 0.03978873577297383395

double Lw1=1.0;
double Lphi1=1.0;
double Lz=1.0;
double Lw=1.0;
double Lphi=1.0;

double Llogzmin=1.0;
double Lzm=1.0;
double Lwm=1.0;
double Lw1m=1.0;
double Llogz=1.0;
double Llogzm=1.0;
double FACNew=1.0;

std::complex<double> LSpa[6][6];

std::complex<double> LSpb[6][6];

double LS[6][6];

typedef hadronic_event<lorentzvector<double>, hadronic_event_traits<2U,2U,0U> > _Event220;

void basic_phasespace<_Event220>::phasespace_cuts(double q2min, double q2max) 
{
  if(q2min <= 0.0)  throw "basic_phasespace<...> : Q2min must be greater than zero";
	if(q2min > q2max) throw "basic_phasespace<...> : Q2min must be less than Q2max";
	if(q2min > _M_s)  throw "basic_phasespace<...> : Q2min must be less than s";
	if(q2max > _M_s) q2max = _M_s;
	
	_M_q2min = q2min; _M_q2max = q2max;
}

void basic_phasespace<_Event220>::_S_safety_cut(const event_type& p)
{
  double s = p[-1]*p[0]; 
  int i, j, up = p.upper();
  for(i = -1; i < up; i++)  
    for(j = (i < 1 ? 1 : i+1); j <= up; j++)
        if(p[i]*p[j] < 1.0e-09*s) throw numeric_error();
}



//-------------------------------------------------------------------------------------------
// -- Interference Real qg NLO
//-------------------------------------------------------------------------------------------    
double basic_phasespace<_Event220>::operator()(event_type& p,double &smearing,const double x[])
{
  double _M_s_     = opts.sroot*opts.sroot;
  double _M_q2min_ = phasespace::mmin*phasespace::mmin;
  double _M_q2max_ = phasespace::mmax*phasespace::mmax;
  double mh = coupling::hmass;
  //----- generate the variable Q2 -----    
  double xa, xb, weight = _M_q2max_-_M_q2min_;
  double Q2 = weight*x[0] + _M_q2min_;
  double MMmin =  phasespace::mmin;
  double MMmax =  phasespace::mmax;
  
  double Q = mh - (-mh + MMmin)*pow(-1. + x[0],7) + (-mh + MMmax)*pow(x[0],7);
  Q2 = Q*Q;
  weight = fabs(2.*Q*(-7.*(-mh + MMmin)*pow(-1. + x[0],6) + 7.*(-mh + MMmax)*pow(x[0],6)));
// -- Smearing
  double     SIGMA = opts.SigmaSmearing; // Reccomended value 1.7GeV
  double     Qsolo=sqrt(Q2);
  double     GaussMin = atan(-Qsolo/SIGMA);
  double     GaussMax = atan(Qsolo/SIGMA);
  double     Gauss = GaussMin + x[7] *(GaussMax-GaussMin);
  double     JacGauss = GaussMax-GaussMin;
  double     Qmedido = Qsolo + SIGMA*tan(Gauss);
  JacGauss=JacGauss*SIGMA/(cos(Gauss)*cos(Gauss));
  double     FactorM=exp(-(Qmedido-Qsolo)*(Qmedido-Qsolo)/(2.0*SIGMA*SIGMA))/(SIGMA*sqrt(2.0*__DEB_PI))*JacGauss;
  if (opts.withsmearing) 
    {
      Q2 = Qmedido*Qmedido;
      weight *=FactorM;
    }
  //cout << "Q from PS = " << sqrt(Q2) << endl;
  //----- incoming hadrons -----
  double Eh = 0.5*std::sqrt(_M_s_);
  p[hadron(-1)] = lorentzvector<double>(0.0, 0.0,  Eh, Eh);
  p[hadron( 0)] = lorentzvector<double>(0.0, 0.0, -Eh, Eh);
  double taumin = Q2/_M_s_;
  double ltaumin = log(taumin); //ESTO CAPAZ QUE TIENE QUE SER DEFINIDO POSITIVO
  double tau = pow(taumin,x[1]);
  double ltau = log(tau);
  double Ecm = sqrt(_M_s_*tau);
  double xi = 1.0-Q2/(Ecm*Ecm);
  double tmp7,tmp2;
  double y=0.5*ltau*(1.0-2.0*(tmp7=x[2]));
  //----- incoming momentum fractions -----
  xa = sqrt(tau)*exp(y);
  xb = sqrt(tau)*exp(-y);
  double xabjac = ltau*tau*ltaumin;
  double costh = 1.0 - 2.0*(tmp2=x[3]);
  double sinth = sqrt(1.0-costh*costh); 
  double xphi = __DEB_2PI*x[4]; 
  double Eq = Ecm*xi*0.5;
  lorentzvector<double> q34;
  p[1]=lorentzvector<double>(Eq*cos(xphi)*sinth,Eq*sin(xphi)*sinth,Eq*costh, Eq);
  q34=lorentzvector<double>(-Eq*cos(xphi)*sinth,-Eq*sin(xphi)*sinth,-Eq*costh,sqrt(Eq*Eq+Q2));
  //Generamos a los difotones en su cm, con masa invariante Q2
  double Ef=sqrt(Q2)*0.5;
  double cosx4 = 1.0-2.0*x[5];
  double sinx4 = sqrt(1.0-cosx4*cosx4);	  
  double phi = __DEB_2PI*x[6];	
  p[-2]=lorentzvector<double>(Ef*cos(phi)*sinx4, Ef*sin(phi)*sinx4, Ef*cosx4, Ef);
  p[-3]=lorentzvector<double>(-Ef*cos(phi)*sinx4,-Ef*sin(phi)*sinx4,-Ef*cosx4, Ef);
  //Boost difotones; frame of incommings
  threevector<double> Bvec;
  threevector<double> bvec(q34); 
  double E34=q34.T();
  Bvec = (1.0/E34)*bvec; 
  p[-3].boost(Bvec); p[-2].boost(Bvec);
  double bz = (xa - xb)/(xa + xb);
  p[-3].boost(0, 0, bz);
  p[-2].boost(0, 0, bz);
  p[1].boost(0,0,bz);  
  weight *= xabjac*xi*__DEB_1_8PI*__DEB_1_8PI*(1.0/__DEB_2PI);
  //----- incoming partons -----
  p[-1] = xa*p[hadron(-1)];
  p[ 0] = xb*p[hadron( 0)];
  double s12 = 2.*p[-1]*p[0];
  //  cout<< p[-1]+p[ 0]-p[ 1]-p[ -2]-p[ -3]<<endl;
  // -- FLUX FACTOR added here
  return weight/2./s12;
}

//-------------------------------------------------------------------------------------------
// -- Interference CT gg NLO
//-------------------------------------------------------------------------------------------  
double basic_phasespace<_Event220>::operator()(event_type& p,double &Q2,double &smearing, bool flagrmed,double &z1, double &z2,const double x[])
{
  double _M_s_     = opts.sroot*opts.sroot;
  double _M_q2min_ = phasespace::mmin*phasespace::mmin;
  double _M_q2max_ = phasespace::mmax*phasespace::mmax;
  double  weight   = _M_q2max_-_M_q2min_;
  double deltaQ2   = _M_q2max_-_M_q2min_;
  Q2 = weight*x[0] + _M_q2min_;
  // -- Smearing
  double     SIGMA = opts.SigmaSmearing; // Reccomended value 1.7GeV
  double     Qsolo=sqrt(Q2);
  double     GaussMin = atan(-Qsolo/SIGMA);
  double     GaussMax = atan(Qsolo/SIGMA);
  double     Gauss = GaussMin + x[7] *(GaussMax-GaussMin);
  double     JacGauss = GaussMax-GaussMin;
  double     Qmedido = Qsolo + SIGMA*tan(Gauss);
  JacGauss=JacGauss*SIGMA/(cos(Gauss)*cos(Gauss));
  double     FactorM=exp(-(Qmedido-Qsolo)*(Qmedido-Qsolo)/(2.0*SIGMA*SIGMA))/(SIGMA*sqrt(2.0*__DEB_PI))*JacGauss;
  if (opts.withsmearing && flagrmed) Q2 = Qmedido*Qmedido;
  // -----
  // ----- incoming hadrons -----
  double Eh = 0.5*std::sqrt(_M_s_);
  p[hadron(-1)] = lorentzvector<double>(0.0, 0.0,  Eh, Eh);
  p[hadron( 0)] = lorentzvector<double>(0.0, 0.0, -Eh, Eh);
  double taumin = Q2/_M_s_;
  double ltaumin = log(taumin);  
  double wpsBORN = 0.0;

  // double  xx = taumin + (1.0-taumin)*(x[1]);
  // double  deltataumin = (1.0-taumin);
  // double  x1B = xx;
  // double  x2B = taumin/xx;
  // double GsdPs=deltataumin/_M_s_/xx;

  double y = 0.5*ltaumin*(1.0-2.0*(x[1]));
  double x1B = sqrt(taumin)*exp(y);
  double x2B = sqrt(taumin)*exp(-y); 
  double GsdPs = fabs(ltaumin)/_M_s_;

  double mh= opts.hmass;
  double Gammah = opts.hwidth; // Ok
  //double Ds = ((-mh*mh + Q2)*(-mh*mh + Q2) + mh*mh*Gammah*Gammah);  // In Born mesq.C
  p[-1] = x1B*p[hadron(-1)];
  p[0]  = x2B*p[hadron( 0)];
  
  double costh,sinth,philep;
  double bzB = (x1B - x2B)/(x1B + x2B);

  // -- qT generation -----------------------------------------------------
  // EXTRA TERM Espacio de fases
  // ----------------------------------------------------------------------
  double tmp,pTmax,pTJac,pT2;
  bool qtswitching = !(opts.fixedorder);
  double qtlim;
  double exppy = exp(y);
  double expmy = 1./exppy;
  if(qtswitching)
    {
      //in switching mode (ct and res) do not account for y in the qt limit (actually the qtlimit should be 1e10, i.e. unphysical)                                
      double kinqtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+Q2,2)/(4*pow(opts.sroot,2))-Q2)); //introduced max to avoid neqative argument of sqrt                  
      double switchqtlim = switching::qtlimit(sqrt(Q2));
      qtlim = switchqtlim; //--> in CT and RES should not have any phase space qt limit                                           
    }
  else
    {
      double cosh2y=pow((exppy+expmy)*0.5,2);
      qtlim = sqrt(max(0.,pow(pow(opts.sroot,2)+Q2,2)/(4*pow(opts.sroot,2)*cosh2y)-Q2)); //introduced max to avoid neqative argument of sqrt when y=ymax                
    }
  pTmax = sqrt(_M_s_)/2.;
  double qtmn = (opts.fixedorder) ?
    max(max(opts.qtcut,opts.xqtcut*sqrt(Q2)), phasespace::qtmin)
    : max(opts.qtcutoff,phasespace::qtmin);
  double qtmx = min(qtlim, phasespace::qtmax);
  // cout << " qtmx = " << qtmx << endl;
  // cout << " qtmn = " << qtmn << endl;
  if(!(opts.fixedorder)) pTmax = qtmx;
  // cout << " opts.fixedorder = " << opts.fixedorder << endl;
  // cout << "  pTmax = " << pTmax << endl; 
  pT2 = pTmax*pTmax*pow((tmp=x[6]),3.0) + qtmn*qtmn;
  pTJac = 3.0*pTmax*pTmax*pow(tmp,2.0);
  p[1] = lorentzvector<double>(0.0,sqrt(pT2),0.0,sqrt(pT2));
  //-----------------------------------------------------------------------

  // costh and philep as defined in DYTurbo
  costh = 1.0-2.0*(x[2]);
  sinth = sqrt(1.0-costh*costh);	  
  philep = __DEB_2PI*(x[3]);


  if(opts.fixedorder)
    { // Genera el par de fotones con qT = 0.
      // Generamos a los difotones en su cm, con masa invariante Q2	  
      double Ef=sqrt(Q2)*0.5;
      
      p[-2]= lorentzvector<double>(Ef*cos(philep)*sinth, Ef*sin(philep)*sinth, Ef*costh, Ef);
      p[-3]=lorentzvector<double>(-Ef*cos(philep)*sinth,-Ef*sin(philep)*sinth,-Ef*costh, Ef);
      //Boostear los difotones al cm de las incomings
      //---- boost to the original frame ----
      
      p[-2].boost(0.0,0.0,bzB);
      p[-3].boost(0.0,0.0,bzB);

      // double umio = -2.*(p[-1]*p[-2]);
      // double tmio = -2.*(p[-1]*p[-3]);
      // double t = -(1. + 1.*cosx4)/2.*Q2;
      // double u = -Q2 - t;
      // cout << " t ratio  = " << t/tmio << endl;
      // cout << " u ratio  = " << u/umio << endl;

    }
  else
    { // Genera el par de fotones con momento transverso distinto de 0.
      // Genera el recoil segun la prescripcion elegida           
      //------------------------------------------------------------------
      // what is done in phasespace::set_phiV(-M_PI+2.*M_PI*x[5]);
      double phiV = ( -__DEB_PI+__DEB_2PI*x[8]);
      double cosphiV = cos(phiV);
      double sinphiV = sqrt(max(0.,1.-pow(cosphiV,2)))*(phiV>0 ? 1 : -1);     
      // -----------------------------------------------------------------
      // what is done in omegaintegr::genV4p();
      lorentzvector<double> pV; 
      double mt2 = Q2+pT2;
      double mt = sqrt(mt2);
      double qt = sqrt(pT2);
      pV = lorentzvector<double>(qt*cosphiV,qt*sinphiV,0.5*mt*(exppy-expmy),0.5*mt*(exppy+expmy));
      double beta[3];
      double gam = pV.T()/sqrt(Q2);
      beta[0] = -pV.X()/pV.T();
      beta[1] = -pV.Y()/pV.T();
      beta[2] = -pV.Z()/pV.T();
      // -----------------------------------------------------------------
      // Initiates void phasespace::genRFaxes(restframeid RF)
      //recoil prescriptions are defined by the values of kt1 and kt2     
      double kt1, kt2;
      //CS frame prescription                                                                                                                                                                        
      if (opts.qtrec_cs)
	{
	  kt1 = pV.X()/2.;
	  kt2 = pV.Y()/2.;
	}
      //naive prescription (also called MY prescription in DYRES)                                                                                                                                    
      if (opts.qtrec_naive)
	{
	  kt1=(1.+pV.Z()/(sqrt(Q2)+pV.T()))*pV.X()/2.; //this prescription should be still symmetric, because for kt_1(y)=kt_2(-y) and kt_2(y)=kt_1(-y)              
	  kt2=(1.+pV.Z()/(sqrt(Q2)+pV.T()))*pV.Y()/2.;
	}
      //alternative k1t = 0 prescription                                                                                                                                                             
      if (opts.qtrec_kt0)
	{
	  kt1 = 0; //this prescription is unphysical. It should be at least: kt1 = pV[0] * (phasespace::y < 0 ? 1 : 0); kt2 = pV[1] * (phasespace::y < 0 ? 1 : 0)                                    
	  kt2 = 0;
	}
      //zeta1 as in Eq.(26) of arXiv:1507.06937                                                                                                                                                      
      double ktdotpV = pV.X()*kt1+pV.Y()*kt2;
      double zeta1 = 1./Q2/2.*(Q2+2.*(ktdotpV)+sqrt(pow((Q2+2.*(ktdotpV)),2)-4.*mt2*(pow(kt1,2)+pow(kt2,2))));
      double qP1 = (pV.T()-pV.Z())*opts.sroot/2.;
      //kap1 is the colliding parton a1 after the lorentz transformation from the boson rest frame to the laboratory frame                                                                           
      double kap1[4];
      kap1[3] = opts.sroot/2.*(zeta1*Q2/2./qP1+(pow(kt1,2)+pow(kt2,2))/zeta1*qP1/Q2/pow(opts.sroot,2)*2.);
      kap1[0] = kt1;
      kap1[1] = kt2;
      kap1[2] = opts.sroot/2.*(zeta1*Q2/2./qP1-(pow(kt1,2)+pow(kt2,2))/zeta1*qP1/Q2/pow(opts.sroot,2)*2.);
      //Determine the x,y,z axes of the boson rest frame for any general kt1 prescription                                                                                                            
      //Calculate z axis by boosting back kap1 from the laboratory frame to the boson rest frame       
      double bt[3];
      bt[0]=-beta[0];
      bt[1]=-beta[1];
      bt[2]=-beta[2];
      //----------------------------------------------------------------------------
      double bdotk1=kap1[0]*bt[0]+kap1[1]*bt[1]+kap1[2]*bt[2];
      double zax[3];
      double yax[3];
      double xax[3];
      zax[0]=kap1[0]+gam*bt[0]*(gam/(gam+1)*bdotk1-kap1[3]);
      zax[1]=kap1[1]+gam*bt[1]*(gam/(gam+1)*bdotk1-kap1[3]);
      zax[2]=kap1[2]+gam*bt[2]*(gam/(gam+1)*bdotk1-kap1[3]);

      //normalise z axis to unity                                                                                                               
      double zaxabs = sqrt(pow(zax[0],2)+pow(zax[1],2)+pow(zax[2],2));
      zax[0] = zax[0] / zaxabs;
      zax[1] = zax[1] / zaxabs;
      zax[2] = zax[2] / zaxabs;

      //define y axis ortogonal to z and V direction                                                                                              
      yax[0] = zax[1]*pV.Z() - zax[2]*pV.Y();
      yax[1] = zax[2]*pV.X() - zax[0]*pV.Z();
      yax[2] = zax[0]*pV.Y() - zax[1]*pV.X();

      //normalise y axis to unity                                                                                                                                    
      double yaxabs = sqrt(pow(yax[0],2)+pow(yax[1],2)+pow(yax[2],2));
      yax[0] = yax[0] / yaxabs;
      yax[1] = yax[1] / yaxabs;
      yax[2] = yax[2] / yaxabs;

     //define x axis ortogonal to y and z direction                                                                                                                                    
      xax[0] = yax[1]*zax[2] - yax[2]*zax[1];
      xax[1] = yax[2]*zax[0] - yax[0]*zax[2];
      xax[2] = yax[0]*zax[1] - yax[1]*zax[0];

      //normalise x axis to unity                                                                                                                                         
      double xaxabs = sqrt(pow(yax[0],2)+pow(yax[1],2)+pow(yax[2],2));
      xax[0] = xax[0] / xaxabs;
      xax[1] = xax[1] / xaxabs;
      xax[2] = xax[2] / xaxabs;
      //----------------------------------------------------------------------------
      // what is done in void phasespace::genl4p()
      double rot1[3];
      double c, s;
      c = costh;
      s = sqrt(max(0.,1.-pow(c,2)));
      phasespace::rotate(zax, c, s, yax, rot1);
      //rotate by angle phi_lep with respect to z axis
      double rot2[3];
      c = cos(philep);
      s = sin(philep);
      phasespace::rotate(rot1, c, s, zax, rot2);
      double p3cm[4],p3[4],p4[4];                                                                                                                                                     
      p3cm[3]=sqrt(Q2)/2.; //E                                                                                                                                                                  
      p3cm[0]=p3cm[3]*rot2[0];  //px     
      p3cm[1]=p3cm[3]*rot2[1];  //py 
      p3cm[2]=p3cm[3]*rot2[2];  //pz                                                                                                                                       
      //Boost to go in the lab frame                                                                                                                    
      //dyboost_(gam, beta, p3cm, p3);                                                                                                                                              
      phasespace::boost(gam, beta, p3cm, p3);
      //momentum of the second lepton                                                                                                                                              
      p4[0]=pV.X()-p3[0];      //px                                                                                                                                                   
      p4[1]=pV.Y()-p3[1];      //py                                                                                                                                                  
      p4[2]=pV.Z()-p3[2];      //pz                                                                                                                              
      p4[3]=pV.T()-p3[3];      //E                                                                                                                               
      //----------------------------------------------------------------------------
      
      // cout << " p4[0] = " << p4[0] << endl;
      // cout << " opts.qtrec_cs " << opts.qtrec_cs << endl;
      // cout << " opts.qtrec_naive " << opts.qtrec_naive << endl;
      // cout << " pV.X() = " << pV.X() << endl; 
      // cout << " kap1[0] = " << kap1[0] << " ; " << kt1 << endl;
      // cout << " x = " << p4[0]+p3[0]-pV.X() << endl;
      
      p[-3] = lorentzvector<double>(p3[0],p3[1],p3[2],p3[3]);
      p[-2] = lorentzvector<double>(p4[0],p4[1],p4[2],p4[3]);
      p[1]  = lorentzvector<double>(pV.X(),pV.Y(),pV.Z(),pV.T()); 

      // cout << " ratio mass = " <<  2.0*((p[-2])*(p[-3]))/Q2 << endl; 
      //      cout << " ratio = " << qt/((p[-2]+p[-3]).perp()) << endl;
    }
  





  //Genera la fracción de momento z//Scaled momentum fraction
    double cut = 0.0000001;
  double beta=cut + (1-cut)*(x[4]);
  double alfa=cut + (1-cut)*(x[5]);
  
  z1= pow(x1B,beta);
  z2= pow(x2B,alfa);       
  
  //C----------------------------------------------------------------------C
  //C Suavizado gaussiano
  //c     xx(5) --> smearing
  
  //C----------------------------------------------------------------------C
  
  //-- -------------------------------------------------------------------           
  //-- Implementacion de la función de Martin como está en el paper.
  //-- -------------------------------------------------------------------           
  
           
  
  //   cout <<  p <<endl;
  
  
  
  /*
  // Espacio de Phases de Martin que da bien con el Born
  wpsBORN = GsdPs/pow(2.0,8.0)/__DEB_PI/Q2/Ds*(deltaQ2*2.0); //*(FactorM);
  //
  wpsBORN =  wpsBORN/2.0 ; // divided by two for identical particles
  //
  */
  
  wpsBORN = GsdPs/pow(2.0,8.0)/__DEB_PI*(deltaQ2*2.0); //*(FactorM);
  wpsBORN = wpsBORN/2.0 ; 
  wpsBORN *= pTJac;
  
  // -- Este ultimo factor es para que sirva para el CTNLO y CTNNLO qq
  // -- quizás sea algun promedio sobre los colores y spines que 
  // wpsBORN =  wpsBORN*16.0; // Cosi lo spazio delle fase va bene per il CTNLO(qqbar)
  
  // wpsBORN = 1.0/8.0/__DEB_PI*GsdPs; 
  if(opts.withsmearing) wpsBORN *= FactorM;
  return wpsBORN;	   
}
  

//-------------------------------------------------------------------------------------------
// -- Interference Real gg NLO
//-------------------------------------------------------------------------------------------
double basic_phasespace<_Event220>::operator()(event_type& p, bool flagrmed,double& FacNew,const double x[])
{ 
  double _M_s_ = opts.sroot*opts.sroot;
  double _M_q2min_ = phasespace::mmin*phasespace::mmin;
  double _M_q2max_ = phasespace::mmax*phasespace::mmax;
  
  //    const random_generator& __rng = *_M_rng;
 	  
  //  double x[8];
  
  // x[0] =0.0;
  // x[1] =0.0;
  // x[2] =0.0;
  // x[3] =0.0;
  // x[4] =0.0;
  // x[5] =0.0;
  // x[6] =0.0;
  // x[7] =0.0;
  
  // x[1] = __rng();
  // x[2] = __rng();
  // x[3] = __rng();
  // x[4] = __rng();
  // x[5] = __rng();
  // x[6] = __rng();
  // x[7] = __rng();
  
	  
  double xa, xb, weight =sqrt( _M_q2max_)-sqrt(_M_q2min_);
  
  // double deltaQ2= _M_q2max_-_M_q2min_; Q=Q0+Qbin*(x[7]-0.5);
  
  double Qbin = sqrt( _M_q2max_) - sqrt(_M_q2min_);
  double Q0   = sqrt(_M_q2min_) + 0.5*Qbin;
  
  double Q=Q0+Qbin*(x[7]-0.5);
	  
  //  Q=Q0+Qbin*(x[7]-0.5);
  double Qsq=Q*Q;
  double Q2=Qsq;
  double tau= Q2/_M_s_;
  
  //C----------------------------------------------------------------------C
  //C Suavizado gaussiano
  //c     xx(8) --> smearing    
  double     SIGMA = opts.SigmaSmearing; // Reccomended value 1.7GeV
  double     Qsolo=sqrt(Q2);
  double     GaussMin = atan(-Qsolo/SIGMA);
  double     GaussMax = atan(Qsolo/SIGMA);
  double     Gauss = GaussMin +x[0] *(GaussMax-GaussMin);
  double     JacGauss = GaussMax-GaussMin;
  double     Qmedido = Qsolo + SIGMA*tan(Gauss);
  JacGauss=JacGauss*SIGMA/(cos(Gauss)*cos(Gauss));
  double     FactorM=exp(-(Qmedido-Qsolo)*(Qmedido-Qsolo)/(2.0*SIGMA*SIGMA))/(SIGMA*sqrt(2.0*__DEB_PI))*JacGauss;
  
  if (opts.withsmearing && flagrmed) Q2 = Qmedido*Qmedido;
  
  
  //  setSigma0();
  
  double xA;
  double xA0=xA=pow(tau,x[1]);
  Lw1=x[2];
  Lphi1=2.0*PI*x[3];
  Lz=pow(tau/xA0,x[4]);
  Lw=x[5];
  Lphi=Lphi1+2.0*PI*x[6];
  double xB0=tau/xA;
  double xB=xB0/Lz;
  
  Llogzmin=log(1.0-xB0);
  double jacob0=(2.0*sqrt(Q2))*tau*(-log(tau));
  double jacob=jacob0*(-log(tau/xA));
  Lzm=1.0-Lz;
  Lwm=1.0-Lw;
  Lw1m=1.0-Lw1;
  Llogz=log(Lz);
  Llogzm=log(Lzm);
  
  double Ecm = sqrt(_M_s*xA*xB);
  double xi =  1.0-Q2/(Ecm*Ecm);
  
  particle k10,k20,k11,k1A,k1B,k21,k2A,k2B,p1;
  k10=k20=photon;
  double cosTheta1=1.0-2.0*Lw1;
  
  k10.setMom(sqrt(Q2)/2.0,cosTheta1,Lphi1);
  k20.setMom(sqrt(Q2)/2.0,-cosTheta1,Lphi1+PI);
  k11=k1A=k1B=k10;
  k21=k2A=k2B=k20;
 
  double b=(xA0-xB0)/(xA0+xB0);
  k10.boostZ(b);
  k20.boostZ(b);
  
  b=(xA*Lz-xB)/(xA*Lz+xB);
  k1A.boostZ(b);
  k2A.boostZ(b);
  b=(xA-xB*Lz)/(xA+xB*Lz);
  k1B.boostZ(b);
  k2B.boostZ(b);
  double cosTheta=1.0-2.0*Lw;
  double rshat=sqrt(_M_s_*xA*xB);
  double energy=(1.0-Lz)*rshat/2.0;
  p1=particle(0.0,0.0,     " parton   ");
  p1.setMom(energy,cosTheta,Lphi);
  
  double theta=acos(cosTheta);
  k11.invBoost(0.0,theta,Lphi);
  k21.invBoost(0.0,theta,Lphi);
  b=(Lz-1.0)/(1.0+Lz);
  k11.boost(b,theta,Lphi);
  k21.boost(b,theta,Lphi);
  b=(xA-xB)/(xA+xB);
  p1.boostZ(b);
  k11.boostZ(b);
  k21.boostZ(b);	 
  
  /*
    dixon pete;
    pete.z = z;
    pete.zm = zm;
    pete.w = w;
    pete.wm = wm;
    pete.w1 = w1;
    pete.w1m = w1m;
    pete.logz= logz;
    pete.logzm= logzm;
  */
  //----- incoming hadrons -----
  double Eh = 0.5*std::sqrt(_M_s_);
  p[hadron(-1)] = lorentzvector<double>(0.0, 0.0,  Eh, Eh);
  p[hadron( 0)] = lorentzvector<double>(0.0, 0.0, -Eh, Eh);
  
  p[1]=lorentzvector<double>(p1.x(),p1.y(),p1.z(),p1.e());
  p[-2]=lorentzvector<double>(k11.x(),k11.y(),k11.z(),k11.e());
  p[-3]=lorentzvector<double>(k21.x(),k21.y(),k21.z(),k21.e());
  
  p[-1] = (xA)*p[hadron(-1)];
  p[ 0] = (xB)*p[hadron( 0)];
  
  //   cout<< p[-1]+p[ 0]+p[ 1]+p[ -2]+p[ -3]<<endl;
  
  //  cout << k11.z() << endl;
  double zr=sqrt(Lz);
  double zmr=sqrt(Lzm);
  double wr=sqrt(Lw);
  double wmr=sqrt(Lwm);
  double w1r=sqrt(Lw1);
  double w1mr=sqrt(Lw1m);
  complex<double> eplus=exp(I*Lphi/2.0);
  complex<double> eminus=exp(-I*Lphi/2.0);
  complex<double> eplus1=exp(I*(Lphi-Lphi1)/2.0);
  complex<double> eminus1=exp(-I*(Lphi-Lphi1)/2.0);
  
  //  1=pA, 2=pB, 3=p1, 4=k1, 5=k2  (4 and 5 are photons)
  
  //  Note that we should be able to remove the overall eplus
  //  and eminus phases by changing the phases of |1> and |2>.
  
  // The Spa[i][j] and Spb[i][j] are normalized by rshat=sqrt(ss*xA*xB).
  // The S[i][j] are normalized by shat=ss*xA*xB.
  
  LSpa[1][2]=1.0;
  LSpa[1][3]=zmr*wr*eplus;
  LSpa[2][3]=-zmr*wmr*eminus;
  LSpa[1][4]=eplus*(wr*wmr*w1mr*(zr-1.0)*eplus1+(Lw*zr+Lwm)*w1r*eminus1);
  LSpa[2][4]=-eminus*(wr*wmr*w1r*(zr-1.0)*eminus1+(Lwm*zr+Lw)*w1mr*eplus1);
  LSpa[3][4]=zmr*(-w1mr*wr*eplus1+wmr*w1r*eminus1);
  LSpa[1][5]=eplus*(-wr*wmr*w1r*(zr-1.0)*eplus1+(Lw*zr+Lwm)*w1mr*eminus1);
  LSpa[2][5]=-eminus*(wr*wmr*w1mr*(zr-1.0)*eminus1-(Lwm*zr+Lw)*w1r*eplus1);
  LSpa[3][5]=zmr*(w1r*wr*eplus1+wmr*w1mr*eminus1);
  LSpa[4][5]=zr;
  
  for(int i=1;i<6;i++)
    for (int j=i+1;j<6;j++) {
      LSpa[j][i]=-LSpa[i][j];
      if ((i<3)&&(j!=2)) {
	LSpb[j][i]=-conj(LSpa[i][j]);
      } else {
	LSpb[j][i]=conj(LSpa[i][j]);
      }            
      LSpb[i][j]=-LSpb[j][i];
      LS[i][j]=LS[j][i]=real(LSpa[i][j]*LSpb[j][i]);
    }

  
  double FLUX = 4.0*(p[-1]*p[0]);
  if(opts.withsmearing) FLUX /= FactorM;
  FacNew = 1.0/(__DEB_1_8PI*__DEB_1_8PI*(1.0/__DEB_2PI))/xi/Q2/Q2*Lz;
  FACNew=FLUX*FacNew;
  //return Qbin/Q2/Q2*jacob; // Espacio de fases original como esta en el codigo de Dixon
  return Qbin/Q2/Q2*jacob/FACNew; // Este es el espacio de fases Correcto ya dividido por el Flujo
}

//-------------------------------------------------------------------------------------------
// -- Interference qqbar Real LO
//-------------------------------------------------------------------------------------------
double basic_phasespace<_Event220>::operator()(event_type& p,double &Q2,double &smearing, int & flag, bool flagrmed, const double x[])
{
  // if(x[0] < 0.00001) return 0.;

    double _M_q2min_ = phasespace::mmin*phasespace::mmin;
    double _M_q2max_ = phasespace::mmax*phasespace::mmax;
     
    // double xx[9]
    // xx[0] = 0.14198;
    // xx[1] = 0.112612;
    // xx[2] = 0.0985368;
    // xx[3] = 0.39425;
    // xx[4] = 0.165576;
    // xx[5] = 0.413769;
    // xx[6] = 0.652282;
    // xx[7] = 0.182755;

    // cout << " flag : " << flag << endl;
    // cout << " flagr : " << flagr << endl;
    
    flag  = 0;

    double wpsreal=0.0;
    //----- incoming hadrons -----
    double _M_s_ = opts.sroot*opts.sroot;
    double Eh = 0.5*std::sqrt(_M_s_);
    p[hadron(-1)] = lorentzvector<double>(0.0, 0.0,  Eh, Eh);
    p[hadron( 0)] = lorentzvector<double>(0.0, 0.0, -Eh, Eh);
    
    //C----------------------------------------------------------------------C
    //C Inicio de variables de integracion y calculo del jacobiano
    //C correspondiente
    double sh =_M_s_;
    double MH = opts.hmass;
    double GH = opts.hwidth; // Ok;     
    //c     xx(1) --> tau
    double tau=MH*MH/sh+x[0] *(1.0-MH*MH/sh);
    double taujac=(1.0-MH*MH/sh);
    //c     xx(2) --> y_cm
    double  ymax= -log(tau)/2.0;
    double  ymin=  log(tau)/2.0;
    double  ycm = ymin +x[1] *(ymax-ymin);
    double  ycmjac=(ymax-ymin);
    double  s=sh*tau;
    //c randoms para la creacion de los dos (q/g y Higgs)      
    //c     xx(3) --> yj
    double   xx3 = x[2] ;
    double   tita5=__DEB_PI*xx3*xx3;
    double   tita5jac=2.0*__DEB_PI*xx3;
    double   yj5=cos(tita5);
    tita5jac=tita5jac*sin(tita5);
    //c     xx(4) --> phij
    double    phi5=2.0*__DEB_PI*x[3] ;
    double    phi5jac=2.0*__DEB_PI;
    //c randoms para el decaimiento del higgs a dos fotones (p3 y p4)
    //c     xx(5) --> yj
    double   xx5 = x[4] ;
    double   tita3=__DEB_PI*xx5*xx5;
    double   tita3jac= 2.0*__DEB_PI*xx5;
    double   yj3=cos(tita3);
    tita3jac=tita3jac*sin(tita3);
    //c     xx(6) --> phij
    double    phi3=2.0*__DEB_PI*x[5] ;
    double    phi3jac=2.0*__DEB_PI;
    double    tita4=__DEB_PI-tita3;
    double    phi4=__DEB_PI+phi3;    
    ////////////////////////////////////////////////////////////////////////  
    //c random para Q2
    //c     xx(7) --> Q2
    ////////////////////////////////////////////////////////////////////////
    double RhoMin = atan(-MH/GH);
    double RhoMax;
    if(s < 1000000.0){
      RhoMax = atan((s-MH*MH)/MH/GH); }
    else if(s > 1000000.0){
      RhoMax = atan((1000000.0-MH*MH)/MH/GH);
    };
    
    double Rho = RhoMin +x[7]*(RhoMax-RhoMin);
    double JacRho = MH*GH*(RhoMax-RhoMin)/(cos(Rho)*cos(Rho));
    Q2 = MH*MH + MH*GH*tan(Rho);

    // // // Standarizando la forma en Q2
    // double weight = _M_q2max_-_M_q2min_;
    // double deltaQ2= _M_q2max_-_M_q2min_;
    // Q2 = weight*x[7] + _M_q2min_;
    // double JacRho =deltaQ2;

    // Standarizando la forma en Q
    // double weight = phasespace::mmax-phasespace::mmin;
    // double deltaQ = phasespace::mmax-phasespace::mmin;
    // double Q = weight*x[7] + phasespace::mmin;
    // double JacRho =deltaQ;
    // Q2 = Q*Q;
    double Q2prev = Q2;
    //C----------------------------------------------------------------------C
    //C Suavizado gaussiano
    //c     xx(8) --> smearing    
    double     SIGMA = opts.SigmaSmearing; // Reccomended value 1.7GeV
    double     Qsolo=sqrt(Q2);
    double     GaussMin = atan(-Qsolo/SIGMA);
    double     GaussMax = atan(Qsolo/SIGMA);
    double     Gauss = GaussMin +x[6] *(GaussMax-GaussMin);
    double     JacGauss = GaussMax-GaussMin;
    double     Qmedido = Qsolo + SIGMA*tan(Gauss);
    JacGauss = JacGauss*SIGMA/(cos(Gauss)*cos(Gauss));
    double     FactorM=exp(-(Qmedido-Qsolo)*(Qmedido-Qsolo)/(2.0*SIGMA*SIGMA))/(SIGMA*sqrt(2.0*__DEB_PI))*JacGauss;
    
    //C----------------------------------------------------------------------C
    //C Armado del jacobiano
    if (opts.withsmearing && flagrmed) Q2 = Qmedido*Qmedido;
    double   jacob=taujac*ycmjac*tita5jac*phi5jac*tita3jac*phi3jac*JacRho;
    double   x1B=sqrt(tau) * exp(ycm);
    double   x2B=tau/x1B;
    double   mod5=(s-Q2)/2.0/sqrt(s);
    double   modQ=mod5;
    double   titaQ=__DEB_PI-tita5;
    double   phiQ=__DEB_PI+phi5;
    //C----------------------------------------------------------------------C    
    // generaciond e los momentos fisicos    
    double p3[4],p4[4],q3[4],q4[4];    
    
    //    if (opts.withsmearing && flagrmed) Q2 = Qmedido*Qmedido;

    // double   mod5=(s-Q2)/2.0/sqrt(s);
    // double   modQ=mod5;
    
    double      p51=sqrt(mod5*mod5);
    double      p52=mod5*sin(tita5)*cos(phi5);
    double      p53=mod5*sin(tita5)*sin(phi5);
    double      p54=mod5*cos(tita5);
    
    double      Q1l=sqrt(Q2+modQ*modQ);
    double      Q2l=modQ*sin(titaQ)*cos(phiQ);
    double      Q3l=modQ*sin(titaQ)*sin(phiQ);
    double      Q4l=modQ*cos(titaQ);
    double      p11=0.5*sqrt(s);
    double      p12=0.0;
    double      p13=0.0;
    double      p14=0.5*sqrt(s);
    double      p21=p11;
    double      p22=0.0;
    double      p23=0.0;
    double      p24=-p14;
    
    q3[0]=sqrt(Q2)/2.0;
    q3[1]=sqrt(Q2)/2.0*sin(tita3)*cos(phi3);
    q3[2]=sqrt(Q2)/2.0*sin(tita3)*sin(phi3);
    q3[3]=sqrt(Q2)/2.0*cos(tita3);
    q4[0]=sqrt(Q2)/2.0;
    q4[1]=sqrt(Q2)/2.0*sin(tita4)*cos(phi4);
    q4[2]=sqrt(Q2)/2.0*sin(tita4)*sin(phi4);
    q4[3]=sqrt(Q2)/2.0*cos(tita4);
    
    //c  mgamma es k1(1)=sqrt(s/4)
    
    double      gama=Q1l/sqrt(Q2)   ;
    double      beta11=p52/p51   ;
    double      beta12=p53/p51 ;
    double      beta13=p54/p51;
    double      beta21=Q2l/Q1l;
    double      beta22=Q3l/Q1l;
    double      beta23=Q4l/Q1l;
    
    double lor1[5][5];
    double lor2[5][5];
    
    lor1[1][1]=gama;
    
    lor1[1][2]=beta11*gama;
    lor1[1][3]=beta12*gama;
    lor1[1][4]=beta13*gama;
    lor1[2][1]=beta11*gama;
    lor1[2][2]=1.0+(gama-1.0)*beta11*beta11/(beta11*beta11+beta12*beta12+beta13*beta13);
    
    lor1[2][3]=(gama-1.0)*beta11*beta12/(beta11*beta11+beta12*beta12+beta13*beta13);
    
    lor1[2][4]=(gama-1.0)*beta11*beta13/(beta11*beta11+beta12*beta12+beta13*beta13);
    
    lor1[3][1]=beta12*gama;
    
    lor1[3][2]=(gama-1.0)*beta12*beta11/(beta11*beta11+beta12*beta12+beta13*beta13);
    
    lor1[3][3]=1.0+(gama-1.0)*beta12*beta12/(beta11*beta11+beta12*beta12+beta13*beta13);
    
    lor1[3][4]=(gama-1.0)*beta12*beta13/(beta11*beta11+beta12*beta12+beta13*beta13);
    
    lor1[4][1]=beta13*gama;
    
    lor1[4][2]=(gama-1.0)*beta13*beta11/(beta11*beta11+beta12*beta12+beta13*beta13);
    
    lor1[4][3]=(gama-1.0)*beta13*beta12/(beta11*beta11+beta12*beta12+beta13*beta13);
    
    lor1[4][4]=1.0+(gama-1.0)*beta13*beta13/(beta11*beta11+beta12*beta12+beta13*beta13);
    
    lor2[1][1]=gama;
    lor2[1][2]=beta21*gama;
    lor2[1][3]=beta22*gama;
    lor2[1][4]=beta23*gama;
    lor2[2][1]=beta21*gama;
    lor2[2][2]=1.0+(gama-1.0)*beta21*beta21/(beta21*beta21+beta22*beta22+beta23*beta23);
    
    
    lor2[2][3]=(gama-1.0)*beta21*beta22/(beta21*beta21+   beta22*beta22+beta23*beta23);
    
    lor2[2][4]=(gama-1.0)*beta21*beta23/(beta21*beta21+ beta22*beta22+beta23*beta23);
    
    lor2[3][1]=beta22*gama;
    
    lor2[3][2]=(gama-1.0)*beta22*beta21/(beta21*beta21+ beta22*beta22+beta23*beta23);
    
    lor2[3][3]=1.0+(gama-1.0)*beta22*beta22/(beta21*beta21+ beta22*beta22+beta23*beta23);
    
    lor2[3][4]=(gama-1.0)*beta22*beta23/(beta21*beta21+ beta22*beta22+beta23*beta23);
    
    lor2[4][1]=beta23*gama;
    
    lor2[4][2]=(gama-1.0)*beta23*beta21/(beta21*beta21+ beta22*beta22+beta23*beta23);
    
    lor2[4][3]=(gama-1.0)*beta23*beta22/(beta21*beta21+ beta22*beta22+beta23*beta23);
    
    lor2[4][4]=1.0+(gama-1.0)*beta23*beta23/(beta21*beta21+ beta22*beta22+beta23*beta23);
    
    
    p3[0] =  (lor2[1][1]*q3[0] + lor2[1][2]*q3[1]+lor2[1][3]*q3[2]+lor2[1][4]*q3[3]);
    p3[1] = -(lor2[2][1]*q3[0] + lor2[2][2]*q3[1]+lor2[2][3]*q3[2]+lor2[2][4]*q3[3]);
    p3[2] = -(lor2[3][1]*q3[0] + lor2[3][2]*q3[1]+lor2[3][3]*q3[2]+lor2[3][4]*q3[3]);
    p3[3] = -(lor2[4][1]*q3[0] + lor2[4][2]*q3[1]+lor2[4][3]*q3[2]+lor2[4][4]*q3[3]);
    
    p4[0] =  (lor2[1][1]*q4[0] + lor2[1][2]*q4[1]+lor2[1][3]*q4[2]+lor2[1][4]*q4[3]);
    p4[1] = -(lor2[2][1]*q4[0] + lor2[2][2]*q4[1]+lor2[2][3]*q4[2]+lor2[2][4]*q4[3]);
    p4[2] = -(lor2[3][1]*q4[0] + lor2[3][2]*q4[1]+lor2[3][3]*q4[2]+lor2[3][4]*q4[3]);
    p4[3] = -(lor2[4][1]*q4[0] + lor2[4][2]*q4[1]+lor2[4][3]*q4[2]+lor2[4][4]*q4[3]);
    
    double p5[4];
    
    p5[0]=p51;       
    p5[1]=-p52;
    p5[2]=-p53;
    p5[3]=-p54;
    
    // These are in the center of mass frame of the incoming partons
    p[1]  = lorentzvector<double>( p5[1], p5[2], p5[3], -p5[0]);// -- la energia va con menos
    p[-2] = lorentzvector<double>( p3[1], p3[2], p3[3], -p3[0]);// -- la energia va con menos
    p[-3] = lorentzvector<double>( p4[1], p4[2], p4[3], -p4[0]);// -- la energia va con menos
        
    p[-1] = lorentzvector<double>( p12,p13, p14,p11);
    p[0]  = lorentzvector<double>( p22, p23,p24,p21);
    
    // cout << " (p[-1]*p[0])/(p[-1]*p[-2])  :" << (p[-1]*p[0])/(p[-1]*p[-2]) << endl;
    /*
      cout << "++++++++++++++++++++++++++++++++++++++++++++++++++ "<< endl;
      cout << "Antes del Boost p[0]*p[-2]   :" << p[0]*p[-2] << endl;
    */
    // if(opts.withsmearing && flagrmed)
    //   {

    //incoming partons
    p[-1] = x1B*p[hadron(-1)];
    p[0]  = x2B*p[hadron( 0)];
    // Boost to the lab frame
    
    double bz = (x1B - x2B)/(x1B + x2B);
    p[-3].boost(0., 0., bz);
    p[-2].boost(0., 0., bz);
    p[1].boost(0.,0.,bz);
      // }
    double flux = 2.0*s; 
    // all factors --- hc2 only missing factor to be included later     
    wpsreal = jacob/flux/(16.0*__DEB_PI)/2.0/__DEB_PI*sqrt((Q2prev-s)*(Q2prev-s)/4.0/s)/(16.0*__DEB_PI*__DEB_PI)/sqrt(s)/2.0/__DEB_PI; // (*FactorM)
    wpsreal = wpsreal/2.0*2.0;   // Particulas Identicas :P ! 2*Re
    //    cout << " wpsreal before = " << wpsreal << endl;
    if(opts.withsmearing) 
      {
	wpsreal *= FactorM;
      }
 
    //   if(flagrmed)   cout << " sqrt(2.*p[-2]*p[-3])/Qmedido = " << sqrt(2.*p[-2]*p[-3])/Qmedido  << endl; 

    // cout << " sqrt(Q2) = " << sqrt(Q2) << endl;
    // cout << " Qmedido = " << Qmedido << endl;
    // cout << " opts.withsmearing : " << opts.withsmearing << endl;
    // cout << " opts.opts.SigmaSmearing = " << opts.SigmaSmearing << endl;
    // cout << " wpsreal smearing = " << wpsreal << endl;

    //    cout << " p[-1]+  p[0] -( p[1] + p[-2] + p[-3]) = " << p[-1]+  p[0] -( p[1] + p[-2] + p[-3]) << endl;

    // if ( (p[-1]+  p[0] + p[1] + p[-2] + p[-3]).Z() > 1.1)
    //   {
    // 	cout <<" Controlando cinematica " << endl;
    // 	//	cout <<" p : " << p << endl;
    // 	cout << " conserva  :" <<     p[-1] + p[0] + p[1] + p[-2] + p[-3] << endl;     
    //     cout << " wpsreal : " << wpsreal << endl;
    // 	cout << " x[0] : " << x[0] << endl;
    // 	cout << " x[1] : " << x[1] << endl;
    // 	cout << " x[2] : " << x[2] << endl;
    // 	cout << " x[3] : " << x[3] << endl;
    // 	cout << " x[4] : " << x[4] << endl;
    // 	cout << " x[5] : " << x[5] << endl;
    // 	cout << " x[6] : " << x[6] << endl;
    // 	cout << " x[7] : " << x[7] << endl;
    // 	cout << " tau : " << tau << endl;
    // 	cout << " ymax : " << ymax << endl;
    // 	cout << " ymin : " << ymin << endl;	
    //   }
 
    //    cout << " wpsreal : " << wpsreal << endl;
   
    return wpsreal;    
  }  























//-------------------------------------------------------------------------------------------
// -- Interference qg alternative Real LO
//-------------------------------------------------------------------------------------------
double basic_phasespace<_Event220>::operator()(event_type& p,double &Q2,bool flagrmed, const double x[])
{
  // if(x[0] < 0.00001) return 0.;

    double _M_q2min_ = phasespace::mmin*phasespace::mmin;
    double _M_q2max_ = phasespace::mmax*phasespace::mmax;
     
    // double xx[9]
    // xx[0] = 0.14198;
    // xx[1] = 0.112612;
    // xx[2] = 0.0985368;
    // xx[3] = 0.39425;
    // xx[4] = 0.165576;
    // xx[5] = 0.413769;
    // xx[6] = 0.652282;
    // xx[7] = 0.182755;

    // cout << " flag : " << flag << endl;
    // cout << " flagr : " << flagr << endl;
    
    double wpsreal=0.0;
    //----- incoming hadrons -----
    double _M_s_ = opts.sroot*opts.sroot;
    double Eh = 0.5*std::sqrt(_M_s_);
    p[hadron(-1)] = lorentzvector<double>(0.0, 0.0,  Eh, Eh);
    p[hadron( 0)] = lorentzvector<double>(0.0, 0.0, -Eh, Eh);
    
    //C----------------------------------------------------------------------C
    //C Inicio de variables de integracion y calculo del jacobiano
    //C correspondiente
    double sh =_M_s_;
    double MH = coupling::hmass;
    double GH = coupling::hwidth; // Ok;     
    //c     xx(1) --> tau
    double tau=MH*MH/sh+x[0] *(1.0-MH*MH/sh);
    double taujac=(1.0-MH*MH/sh);
    //c     xx(2) --> y_cm
    double  ymax= -log(tau)/2.0;
    double  ymin=  log(tau)/2.0;
    double  ycm = ymin +x[1] *(ymax-ymin);
    double  ycmjac=(ymax-ymin);
    double  s=sh*tau;
    //c randoms para la creacion de los dos (q/g y Higgs)      
    //c     xx(3) --> yj
    double   xx3 = x[2] ;
    double   tita5=__DEB_PI*xx3*xx3;
    double   tita5jac=2.0*__DEB_PI*xx3;
    double   yj5=cos(tita5);
    tita5jac=tita5jac*sin(tita5);
    //c     xx(4) --> phij
    double    phi5=2.0*__DEB_PI*x[3] ;
    double    phi5jac=2.0*__DEB_PI;
    //c randoms para el decaimiento del higgs a dos fotones (p3 y p4)
    //c     xx(5) --> yj
    double   xx5 = x[4] ;
    double   tita3=__DEB_PI*xx5*xx5;
    double   tita3jac= 2.0*__DEB_PI*xx5;
    double   yj3=cos(tita3);
    tita3jac=tita3jac*sin(tita3);
    //c     xx(6) --> phij
    double    phi3=2.0*__DEB_PI*x[5] ;
    double    phi3jac=2.0*__DEB_PI;
    double    tita4=__DEB_PI-tita3;
    double    phi4=__DEB_PI+phi3;    
    ////////////////////////////////////////////////////////////////////////  
    //c random para Q2
    //c     xx(7) --> Q2
    ////////////////////////////////////////////////////////////////////////
    // double RhoMin = atan(-MH/GH);
    // double RhoMax;
    // if (s < 1000000.0)
    //   RhoMax = atan((s-MH*MH)/MH/GH); 
    // else
    //   RhoMax = atan((1000000.0-MH*MH)/MH/GH);
    
    // double Rho = RhoMin +x[6]*(RhoMax-RhoMin);
    // double JacRho = MH*GH*(RhoMax-RhoMin)/(cos(Rho)*cos(Rho));
    // Q2 = MH*MH + MH*GH*tan(Rho);
    
    // Standarizando la forma en Q2
     double weight = _M_q2max_-_M_q2min_;
     double deltaQ2= _M_q2max_-_M_q2min_;
     Q2 = weight*x[6] + _M_q2min_;
     double JacRho =deltaQ2;
    //C----------------------------------------------------------------------C
    //C Suavizado gaussiano
    //c     xx(8) --> smearing    
    double     SIGMA = opts.SigmaSmearing; // Reccomended value 1.7GeV
    double     Qsolo=sqrt(Q2);
    double     GaussMin = atan(-Qsolo/SIGMA);
    double     GaussMax = atan(Qsolo/SIGMA);
    double     Gauss = GaussMin +x[7] *(GaussMax-GaussMin);
    double     JacGauss = GaussMax-GaussMin;
    double     Qmedido = Qsolo + SIGMA*tan(Gauss);
    JacGauss=JacGauss*SIGMA/(cos(Gauss)*cos(Gauss));
    double     FactorM=exp(-(Qmedido-Qsolo)*(Qmedido-Qsolo)/(2.0*SIGMA*SIGMA))/(SIGMA*sqrt(2.0*__DEB_PI))*JacGauss;
    if (opts.withsmearing && flagrmed) Q2 = Qmedido*Qmedido;
    //C----------------------------------------------------------------------C
    //C Armado del jacobiano
    double   jacob=taujac*ycmjac*tita5jac*phi5jac*tita3jac*phi3jac*JacRho;
    double   x1B=sqrt(tau) * exp(ycm);
    double   x2B=tau/x1B;
    double   mod5=(s-Q2)/2.0/sqrt(s);
    double   modQ=mod5;
    double   titaQ=__DEB_PI-tita5;
    double   phiQ=__DEB_PI+phi5;
    //C----------------------------------------------------------------------C    
    // generaciond e los momentos fisicos    
    double p3[4],p4[4],q3[4],q4[4];    
    
    double      p51=sqrt(mod5*mod5);
    double      p52=mod5*sin(tita5)*cos(phi5);
    double      p53=mod5*sin(tita5)*sin(phi5);
    double      p54=mod5*cos(tita5);
    
    double      Q1l=sqrt(Q2+modQ*modQ);
    double      Q2l=modQ*sin(titaQ)*cos(phiQ);
    double      Q3l=modQ*sin(titaQ)*sin(phiQ);
    double      Q4l=modQ*cos(titaQ);
    double      p11=0.5*sqrt(s);
    double      p12=0.0;
    double      p13=0.0;
    double      p14=0.5*sqrt(s);
    double      p21=p11;
    double      p22=0.0;
    double      p23=0.0;
    double      p24=-p14;
    
    q3[0]=sqrt(Q2)/2.0;
    q3[1]=sqrt(Q2)/2.0*sin(tita3)*cos(phi3);
    q3[2]=sqrt(Q2)/2.0*sin(tita3)*sin(phi3);
    q3[3]=sqrt(Q2)/2.0*cos(tita3);
    q4[0]=sqrt(Q2)/2.0;
    q4[1]=sqrt(Q2)/2.0*sin(tita4)*cos(phi4);
    q4[2]=sqrt(Q2)/2.0*sin(tita4)*sin(phi4);
    q4[3]=sqrt(Q2)/2.0*cos(tita4);
    
    //c  mgamma es k1(1)=sqrt(s/4)
    
    double      gama=Q1l/sqrt(Q2)   ;
    double      beta11=p52/p51   ;
    double      beta12=p53/p51 ;
    double      beta13=p54/p51;
    double      beta21=Q2l/Q1l;
    double      beta22=Q3l/Q1l;
    double      beta23=Q4l/Q1l;
    
    double lor1[5][5];
    double lor2[5][5];
    
    lor1[1][1]=gama;
    
    lor1[1][2]=beta11*gama;
    lor1[1][3]=beta12*gama;
    lor1[1][4]=beta13*gama;
    lor1[2][1]=beta11*gama;
    lor1[2][2]=1.0+(gama-1.0)*beta11*beta11/(beta11*beta11+beta12*beta12+beta13*beta13);
    
    lor1[2][3]=(gama-1.0)*beta11*beta12/(beta11*beta11+beta12*beta12+beta13*beta13);
    
    lor1[2][4]=(gama-1.0)*beta11*beta13/(beta11*beta11+beta12*beta12+beta13*beta13);
    
    lor1[3][1]=beta12*gama;
    
    lor1[3][2]=(gama-1.0)*beta12*beta11/(beta11*beta11+beta12*beta12+beta13*beta13);
    
    lor1[3][3]=1.0+(gama-1.0)*beta12*beta12/(beta11*beta11+beta12*beta12+beta13*beta13);
    
    lor1[3][4]=(gama-1.0)*beta12*beta13/(beta11*beta11+beta12*beta12+beta13*beta13);
    
    lor1[4][1]=beta13*gama;
    
    lor1[4][2]=(gama-1.0)*beta13*beta11/(beta11*beta11+beta12*beta12+beta13*beta13);
    
    lor1[4][3]=(gama-1.0)*beta13*beta12/(beta11*beta11+beta12*beta12+beta13*beta13);
    
    lor1[4][4]=1.0+(gama-1.0)*beta13*beta13/(beta11*beta11+beta12*beta12+beta13*beta13);
    
    lor2[1][1]=gama;
    lor2[1][2]=beta21*gama;
    lor2[1][3]=beta22*gama;
    lor2[1][4]=beta23*gama;
    lor2[2][1]=beta21*gama;
    lor2[2][2]=1.0+(gama-1.0)*beta21*beta21/(beta21*beta21+beta22*beta22+beta23*beta23);
    
    
    
    lor2[2][3]=(gama-1.0)*beta21*beta22/(beta21*beta21+   beta22*beta22+beta23*beta23);
    
    lor2[2][4]=(gama-1.0)*beta21*beta23/(beta21*beta21+ beta22*beta22+beta23*beta23);
    
    lor2[3][1]=beta22*gama;
    
    lor2[3][2]=(gama-1.0)*beta22*beta21/(beta21*beta21+ beta22*beta22+beta23*beta23);
    
    lor2[3][3]=1.0+(gama-1.0)*beta22*beta22/(beta21*beta21+ beta22*beta22+beta23*beta23);
    
    lor2[3][4]=(gama-1.0)*beta22*beta23/(beta21*beta21+ beta22*beta22+beta23*beta23);
    
    lor2[4][1]=beta23*gama;
    
    lor2[4][2]=(gama-1.0)*beta23*beta21/(beta21*beta21+ beta22*beta22+beta23*beta23);
    
    lor2[4][3]=(gama-1.0)*beta23*beta22/(beta21*beta21+ beta22*beta22+beta23*beta23);
    
    lor2[4][4]=1.0+(gama-1.0)*beta23*beta23/(beta21*beta21+ beta22*beta22+beta23*beta23);
    
    
    p3[0] =  (lor2[1][1]*q3[0] + lor2[1][2]*q3[1]+lor2[1][3]*q3[2]+lor2[1][4]*q3[3]);
    p3[1] = -(lor2[2][1]*q3[0] + lor2[2][2]*q3[1]+lor2[2][3]*q3[2]+lor2[2][4]*q3[3]);
    p3[2] = -(lor2[3][1]*q3[0] + lor2[3][2]*q3[1]+lor2[3][3]*q3[2]+lor2[3][4]*q3[3]);
    p3[3] = -(lor2[4][1]*q3[0] + lor2[4][2]*q3[1]+lor2[4][3]*q3[2]+lor2[4][4]*q3[3]);
    
    p4[0] =  (lor2[1][1]*q4[0] + lor2[1][2]*q4[1]+lor2[1][3]*q4[2]+lor2[1][4]*q4[3]);
    p4[1] = -(lor2[2][1]*q4[0] + lor2[2][2]*q4[1]+lor2[2][3]*q4[2]+lor2[2][4]*q4[3]);
    p4[2] = -(lor2[3][1]*q4[0] + lor2[3][2]*q4[1]+lor2[3][3]*q4[2]+lor2[3][4]*q4[3]);
    p4[3] = -(lor2[4][1]*q4[0] + lor2[4][2]*q4[1]+lor2[4][3]*q4[2]+lor2[4][4]*q4[3]);
    
    double p5[4];
    
    p5[0]=p51;       
    p5[1]=-p52;
    p5[2]=-p53;
    p5[3]=-p54;
    
    
    // These are in the center of mass frame of the incoming partons
    p[1] = lorentzvector<double>( -p5[1], -p5[2], -p5[3], p5[0]);
    p[-2]= lorentzvector<double>( -p3[1], -p3[2], -p3[3], p3[0]);
    p[-3]= lorentzvector<double>( -p4[1], -p4[2], -p4[3], p4[0]);
    
    
    //  p[-1]=lorentzvector<double>( p12,p13, p14,p11);
    //  p[0]= lorentzvector<double>( p22, p23,p24,p21);
    
    // cout << " (p[-1]*p[0])/(p[-1]*p[-2])  :" << (p[-1]*p[0])/(p[-1]*p[-2]) << endl;
    /*
      cout << "++++++++++++++++++++++++++++++++++++++++++++++++++ "<< endl;
      cout << "Antes del Boost p[0]*p[-2]   :" << p[0]*p[-2] << endl;
    */
    //incoming partons
    p[-1] = x1B*p[hadron(-1)];
    p[0]  = x2B*p[hadron( 0)];
    // Boost to the lab frame
    double bz = (x1B - x2B)/(x1B + x2B);
    p[-3].boost(0., 0., bz);
    p[-2].boost(0., 0., bz);
    p[1].boost(0.,0.,bz);
    
    double flux = 2.0*s; 
    // all factors --- hc2 only missing factor to be included later     
    wpsreal=jacob/flux/(16.0*__DEB_PI)/2.0/__DEB_PI*sqrt((Q2-s)*(Q2-s)/4.0/s)/(16.0*__DEB_PI*__DEB_PI)/sqrt(s)/2.0/__DEB_PI; // (*FactorM)
    // wpsreal = wpsreal/2.0*2.0;   // Particulas Identicas :P ! 2*Re -- Already taken into account in Ampi.cc
    //    cout << " wpsreal before = " << wpsreal << endl;
    if(opts.withsmearing) wpsreal *= FactorM;

    // cout << " sqrt(2.*p[-2]*p[-3])/Qmedido = " << sqrt(2.*p[-2]*p[-3])/Qmedido  << endl; 
    // cout << " sqrt(Q2) = " << sqrt(Q2) << endl;
    // cout << " Qmedido = " << Qmedido << endl;
    // cout << " opts.withsmearing : " << opts.withsmearing << endl;
    // cout << " opts.opts.SigmaSmearing = " << opts.SigmaSmearing << endl;
    // cout << " wpsreal smearing = " << wpsreal << endl;


    //    cout <<" s/p[-1]*p[0] = " << s/(2.*p[-1]*p[0]) << endl;;

    if ( (p[-1]+  p[0] -( p[1] + p[-2] + p[-3])).Z() > 10000.1)
      {
	cout <<" Controlando cinematica " << endl;
	//	cout <<" p : " << p << endl;
	cout << " conserva  :" <<     p[-1] + p[0] -( p[1] + p[-2] + p[-3]) << endl;     
        cout << " wpsreal : " << wpsreal << endl;
	cout << " x[0] : " << x[0] << endl;
	cout << " x[1] : " << x[1] << endl;
	cout << " x[2] : " << x[2] << endl;
	cout << " x[3] : " << x[3] << endl;
	cout << " x[4] : " << x[4] << endl;
	cout << " x[5] : " << x[5] << endl;
	cout << " x[6] : " << x[6] << endl;
	cout << " x[7] : " << x[7] << endl;
	cout << " tau : " << tau << endl;
	cout << " ymax : " << ymax << endl;
	cout << " ymin : " << ymin << endl;
	
      }
 
    //    cout << " wpsreal : " << wpsreal << endl;
   
    return wpsreal;    
  }  





//-------------------------------------------------------------------------------------------
// -- Interference Real qg NLO
//-------------------------------------------------------------------------------------------    
double basic_phasespace<_Event220>::operator()(event_type& p,bool flagrmed,const double x[])
{
  double _M_s_     = opts.sroot*opts.sroot;
  double _M_q2min_ = phasespace::mmin*phasespace::mmin;
  double _M_q2max_ = phasespace::mmax*phasespace::mmax;
  double mh = coupling::hmass;
  //----- generate the variable Q2 -----    
  double xa, xb, weight = _M_q2max_-_M_q2min_;
  double Q2 = weight*x[0] + _M_q2min_;
  double MMmin =  phasespace::mmin;
  double MMmax =  phasespace::mmax;
  
  double Q = mh - (-mh + MMmin)*pow(-1. + x[0],7) + (-mh + MMmax)*pow(x[0],7);
  Q2 = Q*Q;
  weight = fabs(2.*Q*(-7.*(-mh + MMmin)*pow(-1. + x[0],6) + 7.*(-mh + MMmax)*pow(x[0],6)));
// -- Smearing
  double     SIGMA = opts.SigmaSmearing; // Reccomended value 1.7GeV
  double     Qsolo=sqrt(Q2);
  double     GaussMin = atan(-Qsolo/SIGMA);
  double     GaussMax = atan(Qsolo/SIGMA);
  double     Gauss = GaussMin + x[7] *(GaussMax-GaussMin);
  double     JacGauss = GaussMax-GaussMin;
  double     Qmedido = Qsolo + SIGMA*tan(Gauss);
  JacGauss=JacGauss*SIGMA/(cos(Gauss)*cos(Gauss));
  double     FactorM=exp(-(Qmedido-Qsolo)*(Qmedido-Qsolo)/(2.0*SIGMA*SIGMA))/(SIGMA*sqrt(2.0*__DEB_PI))*JacGauss;
  if (opts.withsmearing && flagrmed) Q2 = Qmedido*Qmedido;
  //cout << "Q from PS = " << sqrt(Q2) << endl;
  //----- incoming hadrons -----
  double Eh = 0.5*std::sqrt(_M_s_);
  p[hadron(-1)] = lorentzvector<double>(0.0, 0.0,  Eh, Eh);
  p[hadron( 0)] = lorentzvector<double>(0.0, 0.0, -Eh, Eh);
  double taumin = Q2/_M_s_;
  double ltaumin = log(taumin); //ESTO CAPAZ QUE TIENE QUE SER DEFINIDO POSITIVO
  double tau = pow(taumin,x[1]);
  double ltau = log(tau);
  double Ecm = sqrt(_M_s_*tau);
  double xi = 1.0-Q2/(Ecm*Ecm);
  double tmp7,tmp2;
  double y=0.5*ltau*(1.0-2.0*(tmp7=x[2]));
  //----- incoming momentum fractions -----
  xa = sqrt(tau)*exp(y);
  xb = sqrt(tau)*exp(-y);
  double xabjac = ltau*tau*ltaumin;
  double costh = 1.0 - 2.0*(tmp2=x[3]);
  double sinth = sqrt(1.0-costh*costh); 
  double xphi = __DEB_2PI*x[4]; 
  double Eq = Ecm*xi*0.5;
  lorentzvector<double> q34;
  p[1]=lorentzvector<double>(Eq*cos(xphi)*sinth,Eq*sin(xphi)*sinth,Eq*costh, Eq);
  q34=lorentzvector<double>(-Eq*cos(xphi)*sinth,-Eq*sin(xphi)*sinth,-Eq*costh,sqrt(Eq*Eq+Q2));
  //Generamos a los difotones en su cm, con masa invariante Q2
  double Ef=sqrt(Q2)*0.5;
  double cosx4 = 1.0-2.0*x[5];
  double sinx4 = sqrt(1.0-cosx4*cosx4);	  
  double phi = __DEB_2PI*x[6];	
  p[-2]=lorentzvector<double>(Ef*cos(phi)*sinx4, Ef*sin(phi)*sinx4, Ef*cosx4, Ef);
  p[-3]=lorentzvector<double>(-Ef*cos(phi)*sinx4,-Ef*sin(phi)*sinx4,-Ef*cosx4, Ef);
  //Boost difotones; frame of incommings
  threevector<double> Bvec;
  threevector<double> bvec(q34); 
  double E34=q34.T();
  Bvec = (1.0/E34)*bvec; 
  p[-3].boost(Bvec); p[-2].boost(Bvec);
  double bz = (xa - xb)/(xa + xb);
  p[-3].boost(0, 0, bz);
  p[-2].boost(0, 0, bz);
  p[1].boost(0,0,bz);  
  weight *= xabjac*xi*__DEB_1_8PI*__DEB_1_8PI*(1.0/__DEB_2PI);
  //----- incoming partons -----
  p[-1] = xa*p[hadron(-1)];
  p[ 0] = xb*p[hadron( 0)];
  double s12 = 2.*p[-1]*p[0];
  //  cout<< p[-1]+p[ 0]-p[ 1]-p[ -2]-p[ -3]<<endl;
  // -- FLUX FACTOR added here
  if(opts.withsmearing) return weight/2./s12*FactorM;
  return weight/2./s12;
}








































  
  double basic_phasespace<_Event220>::operator()(event_type& p)
  {
	  double _M_s_ = opts.sroot*opts.sroot;
    double _M_q2min_ = phasespace::mmin*phasespace::mmin;
    double _M_q2max_ = phasespace::mmax*phasespace::mmax;
	  
  //  unsigned int up = p.upper();
 //   if(up < 0)  throw "unable to generate hadron-hadron -> V+other event";
    
    //----- generate the variable Q2 -----
  
    const random_generator& __rng = *_M_rng;
    
        //   cout << "random desde phasespace  :" << __rng() << endl;
    
    double xa, xb, weight = _M_q2max_-_M_q2min_;
 //   double deltaQ2= _M_q2max_-_M_q2min_;
    double Q2 = weight*__rng() + _M_q2min_;
   // double tau = Q2/_M_s_;
   // double ltau = -std::log(tau);
    
   
    //----- incoming hadrons -----
    double Eh = 0.5*std::sqrt(_M_s_);
    p[hadron(-1)] = lorentzvector<double>(0.0, 0.0,  Eh, Eh);
    p[hadron( 0)] = lorentzvector<double>(0.0, 0.0, -Eh, Eh);
    
    //----- outgoing non-qcd particles (in their c.m. frame) -----
    //event_type::iterator pbegin = p.begin();
    //event_type::iterator pend   = p.end();
    //weight *= _M_psgen -> operator()(Q2, pbegin, pbegin+2);
    
    double taumin = Q2/_M_s_;
    double ltaumin = log(taumin); //ESTO CAPAZ QUE TIENE QUE SER DEFINIDO POSITIVO
    double tau = pow(taumin,__rng());
    double ltau = log(tau);
    double Ecm = sqrt(_M_s_*tau);
    double xi = 1.0-Q2/(Ecm*Ecm);
    double tmp7,tmp2;
    
     double y=0.5*ltau*(1.0-2.0*(tmp7=__rng()));
    
     //----- incoming momentum fractions -----
      xa = sqrt(tau)*exp(y);
      xb = sqrt(tau)*exp(-y);
    
    
     double xabjac = ltau*tau*ltaumin;
    
       double costh = 1.0 - 2.0*(tmp2=__rng());
	  double sinth = sqrt(1.0-costh*costh); 
	  double xphi = __DEB_2PI*__rng(); 
	  double Eq = Ecm*xi*0.5;
	lorentzvector<double> q34;
	p[1]=lorentzvector<double>(Eq*cos(xphi)*sinth,Eq*sin(xphi)*sinth,Eq*costh, Eq);
	q34=lorentzvector<double>(-Eq*cos(xphi)*sinth,-Eq*sin(xphi)*sinth,-Eq*costh,sqrt(Eq*Eq+Q2));
	
	// double bpar=Eq/sqrt(Q2);
	  
	  //Generamos a los difotones en su cm, con masa invariante Q2
	  
	  double Ef=sqrt(Q2)*0.5;
	  double cosx4 = 1.0-2.0*__rng();
	  double sinx4 = sqrt(1.0-cosx4*cosx4);	  
	  double phi = __DEB_2PI*__rng();
	
	  p[-2]=lorentzvector<double>(Ef*cos(phi)*sinx4, Ef*sin(phi)*sinx4, Ef*cosx4, Ef);
	  p[-3]=lorentzvector<double>(-Ef*cos(phi)*sinx4,-Ef*sin(phi)*sinx4,-Ef*cosx4, Ef);
	  
	//  leanlorentz qf10(Ef*cos(phi)*sinx4, Ef*sin(phi)*sinx4, Ef*cosx4, Ef);
	//  leanlorentz qf20(-Ef*cos(phi)*sinx4,-Ef*sin(phi)*sinx4,-Ef*cosx4, Ef);
	
    //Boost difotones; frame of incommings
	  threevector<double> Bvec;
	  threevector<double> bvec(q34); 
 //     bvec= assigm(q34);
      double E34=q34.T();
      Bvec = (1.0/E34)*bvec; 
	  
	  p[-3].boost(Bvec); p[-2].boost(Bvec);
	  
    
    
    double bz = (xa - xb)/(xa + xb);
        p[-3].boost(0, 0, bz);
        p[-2].boost(0, 0, bz);
		p[1].boost(0,0,bz);
		
		
//		double pT2 = p[-2].perp();
//		double pT1 = p[-3].perp();
		
//		if(pT2 > pT1) swap(p[-2],p[-3]); 
		
 //     double Ejet=p[1].T();
 //     double pTJet=p[1].perp();
 //     double phiJet=p[1].phi();
        
      //  pB1 = leanlorentz(pTJet*cos(phiJet),pTJet*sin(phiJet),sqrt(Ejet*Ejet-pTJet*pTJet),Ejet);
        
        
     //   cout << p[1] <<endl;
     //   coutvec(pB1);
        
	   
    // --Chequeamos que pasen los cuts para generar
    
    
   weight *= xabjac*xi*__DEB_1_8PI*__DEB_1_8PI*(1.0/__DEB_2PI);
    
    
    //----- incoming partons -----
    p[-1] = xa*p[hadron(-1)];
    p[ 0] = xb*p[hadron( 0)];
   // cout<< p[-1]+p[ 0]-p[ 1]-p[ -2]-p[ -3]<<endl;
    
    
   
   
  
  //  cout << " weight desde born :" << weight<< endl;
    
       return weight;
  }




  double basic_phasespace<_Event220>::operator()(const event_type& q, event_type& p)
  {double _M_s_ = opts.sroot*opts.sroot;
    double _M_q2min_ = phasespace::mmin*phasespace::mmin;
    double _M_q2max_ = phasespace::mmax*phasespace::mmax;

    unsigned int up = p.upper(), uq = q.upper();
  //  if(uq < 0 || up < uq) throw "unable to generate the event";
     
    if( uq != 1 || up != 2) {
    cout << "uq  :" << uq <<endl;
    cout << "up  :" << up <<endl;
}
    
   // cout << " q NLO  :" << q << endl;
   // cout << " p NLO :" << p << endl;
 //         cout << "estoy pasando por la parte NNLO" << endl;
    //----- copy from q to p -----
    event_type::iterator pi;
    event_type::const_iterator qi;
    
   // cout << "pi  :" << pi <<endl;
   // cout << "qi  :" << qi <<endl;
    
    p[hadron(-1)] = q[hadron(-1)];
    p[hadron( 0)] = q[hadron( 0)];
    for(pi = p.begin(), qi = q.begin(); qi < q.end(); pi++, qi++) 
      *pi = *qi; 
    
    // Son siempre distintos up y uq
    
    if(up == uq) return 1.0;
    
    int emit, omit, thd;
    double weight = 1.0;
    const random_generator& __rng = *_M_rng;

    event_type::const_iterator p1 = p.begin()+4;
    event_type::const_iterator pn = p1+1;
    lorentzvector<double> pab(p[hadron(0)] + p[hadron(-1)]);
    pab /= p[hadron(0)]*p[hadron(-1)];
    
 //   cout << "weight  :" << weight <<endl;
    
	//---- generate the dipole emmissions ----
	__fe_clear_exception();

    if(uq == 0) {
		
	cout <<	" uq es igual a cero   " << endl;
      dipole_emission::boost_ifi boost;
      int emit = -1, omit = (int) (2*__rng()) - 1;
      if(emit == omit) emit = 0;
      
      dipole_emission::gendip_ifi(__rng, _M_beta, _M_eps, pab*p[emit], p[emit], p[1], p[omit], boost);
      p[-2] = boost(p[-2]); p[-3] = boost(p[-3]);
      
      //----- calculate the weight -----
      weight /= dipole_emission::jacobi_ifi(_M_beta, _M_eps, pab*p[0], p[0], p[-1], p1, pn)
	+       dipole_emission::jacobi_ifi(_M_beta, _M_eps, pab*p[-1], p[-1], p[0], p1, pn);
      
      __fe_throw_exception();
      
      weight *= 2.0; uq++; 
      if(up == 1) return weight;
    }

    double jac;
    unsigned int n;
    
  //  cout << "weight despues del if :" << weight <<endl;
  //  int ipete=1;
    for(n = uq+1, pn = p1+uq+1; n <= up; n++, pn++) {
		
		
	//	cout << " ipete  :" << ipete <<endl;
	//	cout << " pn iterado  :" << pn ;
		
		
		
      //---- choose which one emit, omit, third ----    
      emit = (int) ((n-1)*__rng()) + 1;if(emit != 1) { emit=1; cout << " emit  :" << emit <<endl;}
      omit = (int) (n*__rng()) - 1; if(omit > 0) { omit=0; cout << " omit  :" << omit <<endl;}
      thd  = (int) (n*__rng()) + 1;
      if(emit == omit) omit = n-1;
      
  //    cout << " omit    :" << omit << endl;
      
  //    cout << " _M_eps  antes :" << _M_eps << endl;
      
      
      
      
 //     cout << " p   antes:" << p <<endl;
      
   
 //     cout << " thd  :" << thd <<endl;
       
      
      //---- generate the dipole emmission ----
      if(omit > 0) dipole_emission::gendip_fff(__rng, _M_beta, _M_eps, p[emit], p[n], p[omit]);
      else dipole_emission::gendip_ffi(__rng, _M_beta, _M_eps, pab*p[omit], p[emit], p[n], p[omit]);
      
      if(thd != (int) n) std::swap(p[thd], p[n]);
      
 // if(thd != (int) n)  {  cout<< "thd   :"  << thd << endl;
//  cout<< "n   :"  << n<< endl;}
  
      
   //       cout << " p[n]   despues:  " << p[n] <<endl;    
      
 //     cout << " omit despues :" << omit << endl;
      
  //      cout << " p   antes:" << p <<endl;
 
 //     cout << " p1   despues:" << *p1 <<endl;
   
   //         cout << " pn   antes:" << *pn <<endl;
      
      //----- calculate the weight -----

      jac = dipole_emission::jacobi_fff(_M_beta, _M_eps, p1, pn) 
	+   dipole_emission::jacobi_ffi(_M_beta, _M_eps, pab*p[ 0], p[ 0], p1, pn)
	+   dipole_emission::jacobi_ffi(_M_beta, _M_eps, pab*p[-1], p[-1], p1, pn);
      
      weight *= n*n*(n-1)/jac;
      __fe_throw_exception();
    //  ipete++;
    }

 

  //  cout << " pn   despues:" << *pn <<endl;
      
 //  cout << " p   despues:" << p <<endl;
    //  cout << -p[-2]-p[-3]-p[1]-p[2]+p[-1]+p[0] << endl;
    
 if(std::isinf( weight) || std::isnan( weight)) return 0.0;
  
    return weight;
  }



double basic_phasespace<_Event220>::operator()(event_type& p,double &Q2,double &z1, double &z2)
  {
	  double _M_s_ = opts.sroot*opts.sroot;
    double _M_q2min_ = phasespace::mmin*phasespace::mmin;
    double _M_q2max_ = phasespace::mmax*phasespace::mmax;
	  
 //   unsigned int up = p.upper();
  //  if(up < 0)  throw "unable to generate hadron-hadron -> V+other event";
    
    //----- generate the variable Q2 -----
  
    const random_generator& __rng = *_M_rng;
    
        //   cout << "random desde phasespace  :" << __rng() << endl;
    
  //  double xa, xb;
    double  weight = _M_q2max_-_M_q2min_;
    double deltaQ2= _M_q2max_-_M_q2min_;
     Q2 = weight*__rng() + _M_q2min_;
   
    
    //----- incoming hadrons -----
    double Eh = 0.5*std::sqrt(_M_s_);
    p[hadron(-1)] = lorentzvector<double>(0.0, 0.0,  Eh, Eh);
    p[hadron( 0)] = lorentzvector<double>(0.0, 0.0, -Eh, Eh);
    
    
    
    
    
    double taumin = Q2/_M_s_;
    double ltaumin = log(taumin);  
    double wpsBORN=0.0;
  //VARIABLE REUTILIZADA XX[0]
    double yo=0.5*ltaumin*(1.0-2.0*(__rng() ));
    double  x1B=sqrt(taumin)*exp(yo);
    double  x2B=sqrt(taumin)*exp(-yo);
    
  //Fake para el momento del Jet no presente a este orden   
 //   lorentzvector pJet_fake;
 //   pJet_fake = lorentzvector(0.0,0.0,0.0,0.0);
    //incoming partons
         p[-1] = x1B*p[hadron(-1)];
         p[0] = x2B*p[hadron( 0)];
         
         //Generamos a los difotones en su cm, con masa invariante Q2
	  
	  double Ef=sqrt(Q2)*0.5;
	  double cosx4 = 1.0-2.0*(__rng() );
	  double sinx4 = sqrt(1.0-cosx4*cosx4);	  
	  double phi = __DEB_2PI*(__rng() );
	  
	  lorentzvector<double> qf1B,qf2B;
	  p[-2]= lorentzvector<double>(Ef*cos(phi)*sinx4, Ef*sin(phi)*sinx4, Ef*cosx4, Ef);
	  p[-3]=lorentzvector<double>(-Ef*cos(phi)*sinx4,-Ef*sin(phi)*sinx4,-Ef*cosx4, Ef);
	 
	  
	  
         
    //Boostear los difotones al cm de las incomings
      	  //---- boost to the original frame ----
        double bzB = (x1B - x2B)/(x1B + x2B);
       
		
        p[-2].boost(0.0,0.0,bzB);
        p[-3].boost(0.0,0.0,bzB);
        
  //Genera la fracción de momento z//Scaled momentum fractions
  
   double cut = 0.0000001;
   double beta=cut + (1-cut)*(__rng() );
   double alfa=cut + (1-cut)*(__rng() );
   
    z1= pow(x1B,beta);
    z2= pow(x2B,alfa);       
           
          
           wpsBORN = deltaQ2*fabs(ltaumin)*__DEB_1_8PI/_M_s_;
           
           return wpsBORN;
    
    
    
    
    
    
    
  }


double basic_phasespace<_Event220>::operator()(event_type& p,double Q2,double &z1, double &z2,double cosx4,double phi)
  {
	  double _M_s_ = opts.sroot*opts.sroot;
    double _M_q2min_ = phasespace::mmin*phasespace::mmin;
    double _M_q2max_ = phasespace::mmax*phasespace::mmax;

	  //   unsigned int up = p.upper();
  //  if(up < 0)  throw "unable to generate hadron-hadron -> V+other event";
    
    //----- generate the variable Q2 -----
  
    const random_generator& __rng = *_M_rng;
    
        //   cout << "random desde phasespace  :" << __rng() << endl;
    
  //  double xa, xb, weight = _M_q2max_-_M_q2min_;
    double deltaQ2= _M_q2max_-_M_q2min_;
//     Q2 = weight*__rng() + _M_q2min_;
   
    
    //----- incoming hadrons -----
    double Eh = 0.5*std::sqrt(_M_s_);
    p[hadron(-1)] = lorentzvector<double>(0.0, 0.0,  Eh, Eh);
    p[hadron( 0)] = lorentzvector<double>(0.0, 0.0, -Eh, Eh);
    
    
   
    
    
    double taumin = Q2/_M_s_;
    double ltaumin = log(taumin);  
    double wpsBORN=0.0;
  //VARIABLE REUTILIZADA XX[0]
    double yo=0.5*ltaumin*(1.0-2.0*(__rng() ));
  double  x1B=sqrt(taumin)*exp(yo);
  double  x2B=sqrt(taumin)*exp(-yo);
  
 
    
  //Fake para el momento del Jet no presente a este orden   
 //   lorentzvector pJet_fake;
 //   pJet_fake = lorentzvector(0.0,0.0,0.0,0.0);
    //incoming partons
         p[-1] = x1B*p[hadron(-1)];
         p[0] = x2B*p[hadron( 0)];
         
         //Generamos a los difotones en su cm, con masa invariante Q2
	  
	  double Ef=sqrt(Q2)*0.5;
//	  double cosx4 = 1.0-2.0*(__rng() );
	  double sinx4 = sqrt(1.0-cosx4*cosx4);	  
//	  double phi = __DEB_2PI*(__rng() );
	  
	  lorentzvector<double> qf1B,qf2B;
	  p[-2]= lorentzvector<double>(Ef*cos(phi)*sinx4, Ef*sin(phi)*sinx4, Ef*cosx4, Ef);
	  p[-3]=lorentzvector<double>(-Ef*cos(phi)*sinx4,-Ef*sin(phi)*sinx4,-Ef*cosx4, Ef);
	 
	  
	  
	  //// Contra termino A NLO Y NNLO
 // EXTRA TERM Espacio de fases
  double tmp;
  double pTmax = sqrt( _M_s_ )/2.0;
  double pT2 = pTmax*pTmax*pow((tmp=( __rng() )),3.0);
  double pTJac = 3.0*pTmax*pTmax*pow(tmp,2.0);
  p[1] = lorentzvector<double>(0.0,sqrt(pT2),0.0,sqrt(pT2));
 // p[1] es muleto no respeta cinematica...La cinematica es BORN
	  
         
    //Boostear los difotones al cm de las incomings
      	  //---- boost to the original frame ----
        double bzB = (x1B - x2B)/(x1B + x2B);
       
		
        p[-2].boost(0.0,0.0,bzB);
        p[-3].boost(0.0,0.0,bzB);
        
  //Genera la fracción de momento z//Scaled momentum fractions
  
   double cut = 0.0000001;
   double beta=cut + (1-cut)*(__rng() );
   double alfa=cut + (1-cut)*(__rng() );
   
    z1= pow(x1B,beta);
    z2= pow(x2B,alfa);       
           
          
           wpsBORN = deltaQ2*fabs(ltaumin)*__DEB_1_8PI/_M_s_*Q2*pTJac;
           
           return wpsBORN;
    
    
       
    
  }
  
  
  
  
  double basic_phasespace<_Event220>::operator()(event_type& p,double& cosx4, double& phi)
  {
	  
	  double _M_s_ = opts.sroot*opts.sroot;
    double _M_q2min_ = phasespace::mmin*phasespace::mmin;
    double _M_q2max_ = phasespace::mmax*phasespace::mmax;

  //  unsigned int up = p.upper();
   // if(up < 0)  throw "unable to generate hadron-hadron -> V+other event";
    
    //----- generate the variable Q2 -----
  
    const random_generator& __rng = *_M_rng;
    
        //   cout << "random desde phasespace  :" << __rng() << endl;
    
    double xa, xb, weight = _M_q2max_-_M_q2min_;
 //   double deltaQ2= _M_q2max_-_M_q2min_;
    double Q2 = weight*__rng() + _M_q2min_;
   // double tau = Q2/_M_s_;
   // double ltau = -std::log(tau);
    
    
    //----- incoming hadrons -----
    double Eh = 0.5*std::sqrt(_M_s_);
    p[hadron(-1)] = lorentzvector<double>(0.0, 0.0,  Eh, Eh);
    p[hadron( 0)] = lorentzvector<double>(0.0, 0.0, -Eh, Eh);
    
    //----- outgoing non-qcd particles (in their c.m. frame) -----
    //event_type::iterator pbegin = p.begin();
    //event_type::iterator pend   = p.end();
    //weight *= _M_psgen -> operator()(Q2, pbegin, pbegin+2);
    
    double taumin = Q2/_M_s_;
    double ltaumin = log(taumin); //ESTO CAPAZ QUE TIENE QUE SER DEFINIDO POSITIVO
    double tau = pow(taumin,__rng());
    double ltau = log(tau);
    double Ecm = sqrt(_M_s_*tau);
    double xi = 1.0-Q2/(Ecm*Ecm);
    double tmp7,tmp2;
    
     double y=0.5*ltau*(1.0-2.0*(tmp7=__rng()));
    
     //----- incoming momentum fractions -----
      xa = sqrt(tau)*exp(y);
      xb = sqrt(tau)*exp(-y);
    
    
     double xabjac = ltau*tau*ltaumin;
    
       double costh = 1.0 - 2.0*(tmp2=__rng());
	  double sinth = sqrt(1.0-costh*costh); 
	  double xphi = __DEB_2PI*__rng(); 
	  double Eq = Ecm*xi*0.5;
	lorentzvector<double> q34;
	p[1]=lorentzvector<double>(Eq*cos(xphi)*sinth,Eq*sin(xphi)*sinth,Eq*costh, Eq);
	q34=lorentzvector<double>(-Eq*cos(xphi)*sinth,-Eq*sin(xphi)*sinth,-Eq*costh,sqrt(Eq*Eq+Q2));
	
//	double bpar=Eq/sqrt(Q2);
	  
	  //Generamos a los difotones en su cm, con masa invariante Q2
	  
	  double Ef=sqrt(Q2)*0.5;
	   cosx4 = 1.0-2.0*__rng();
	  double sinx4 = sqrt(1.0-cosx4*cosx4);	  
	   phi = __DEB_2PI*__rng();
	
	  p[-2]=lorentzvector<double>(Ef*cos(phi)*sinx4, Ef*sin(phi)*sinx4, Ef*cosx4, Ef);
	  p[-3]=lorentzvector<double>(-Ef*cos(phi)*sinx4,-Ef*sin(phi)*sinx4,-Ef*cosx4, Ef);
	  
	//  leanlorentz qf10(Ef*cos(phi)*sinx4, Ef*sin(phi)*sinx4, Ef*cosx4, Ef);
	//  leanlorentz qf20(-Ef*cos(phi)*sinx4,-Ef*sin(phi)*sinx4,-Ef*cosx4, Ef);
	
    //Boost difotones; frame of incommings
	  threevector<double> Bvec;
	  threevector<double> bvec(q34); 
 //     bvec= assigm(q34);
      double E34=q34.T();
      Bvec = (1.0/E34)*bvec; 
	  
	  p[-3].boost(Bvec); p[-2].boost(Bvec);
	  
    
    
    double bz = (xa - xb)/(xa + xb);
        p[-3].boost(0, 0, bz);
        p[-2].boost(0, 0, bz);
		p[1].boost(0,0,bz);
		
		
//		double pT2 = p[-2].perp();
//		double pT1 = p[-3].perp();
		
//		if(pT2 > pT1) swap(p[-2],p[-3]); 
		
 //     double Ejet=p[1].T();
 //     double pTJet=p[1].perp();
  //    double phiJet=p[1].phi();
        
      //  pB1 = leanlorentz(pTJet*cos(phiJet),pTJet*sin(phiJet),sqrt(Ejet*Ejet-pTJet*pTJet),Ejet);
        
        
     //   cout << p[1] <<endl;
     //   coutvec(pB1);
        
	   
    // --Chequeamos que pasen los cuts para generar
    
    
   weight *= xabjac*xi*__DEB_1_8PI*__DEB_1_8PI*(1.0/__DEB_2PI);
    
    
    //----- incoming partons -----
    p[-1] = xa*p[hadron(-1)];
    p[ 0] = xb*p[hadron( 0)];
   // cout<< p[-1]+p[ 0]-p[ 1]-p[ -2]-p[ -3]<<endl;
    
    
   
   
  
  //  cout << " weight desde born :" << weight<< endl;
    
       return weight;
  }
  
  
  
  
  //-----------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------
  //----- Interferencia Espacio de fases LO
  //-----------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------
  
  
  double basic_phasespace<_Event220>::operator()(event_type& p,double &Q2,double &smearing, int & flag)
  {
    double _M_q2min_ = phasespace::mmin*phasespace::mmin;
    double _M_q2max_ = phasespace::mmax*phasespace::mmax;
	  
	  double _M_s_ = opts.sroot*opts.sroot;
 //   unsigned int up = p.upper();
  //  if(up < 0)  throw "unable to generate hadron-hadron -> V+other event";
    
    //----- generate the variable Q2 -----
  
    const random_generator& __rng = *_M_rng;
    
        //   cout << "random desde phasespace  :" << __rng() << endl;
    
  //  double xa, xb;
    double  weight = _M_q2max_-_M_q2min_;
    double deltaQ2= _M_q2max_-_M_q2min_;
     Q2 = weight*__rng() + _M_q2min_;
   
    
    //----- incoming hadrons -----
    double Eh = 0.5*std::sqrt(_M_s_);
    p[hadron(-1)] = lorentzvector<double>(0.0, 0.0,  Eh, Eh);
    p[hadron( 0)] = lorentzvector<double>(0.0, 0.0, -Eh, Eh);
    
    
    
   
    
    double taumin = Q2/_M_s_;
    double ltaumin = log(taumin);  
    double wpsBORN=0.0;


//--- Implementacion Base Normal Phase-Space
  //VARIABLE REUTILIZADA XX[0]
  // descomentar si se va a usar
  /*
    double yo=0.5*ltaumin*(1.0-2.0*(__rng() ));
  double  x1B=sqrt(taumin)*exp(yo);
  double  x2B=sqrt(taumin)*exp(-yo);
*/

  double  xx=taumin + (1.0-taumin)*(__rng() );
 
  double  deltataumin=(1.0-taumin);
  double  x1B=xx;
  double  x2B=taumin/xx;

   double GsdPs=deltataumin/_M_s_/xx;
   double mh=_MH_;
   double Gammah = _GammaH_; // Ok
   double Ds = ((-mh*mh + Q2)*(-mh*mh + Q2) + mh*mh*Gammah*Gammah);  
     

    
  //Fake para el momento del Jet no presente a este orden   
 //   lorentzvector pJet_fake;
 //   pJet_fake = lorentzvector(0.0,0.0,0.0,0.0);
    //incoming partons
         p[-1] = x1B*p[hadron(-1)];
         p[0]  = x2B*p[hadron( 0)];
         
         //Generamos a los difotones en su cm, con masa invariante Q2
	  
	  double Ef=sqrt(Q2)*0.5;
	  double cosx4 = 1.0-2.0*(__rng() );
	  double sinx4 = sqrt(1.0-cosx4*cosx4);	  
	  double phi = __DEB_2PI*(__rng() );
	  
	  lorentzvector<double> qf1B,qf2B;
	  p[-2]= lorentzvector<double>(Ef*cos(phi)*sinx4, Ef*sin(phi)*sinx4, Ef*cosx4, Ef);
	  p[-3]=lorentzvector<double>(-Ef*cos(phi)*sinx4,-Ef*sin(phi)*sinx4,-Ef*cosx4, Ef);
	 
	  
	  
         
    //Boostear los difotones al cm de las incomings
      	  //---- boost to the original frame ----
        double bzB = (x1B - x2B)/(x1B + x2B);
       
		
        p[-2].boost(0.0,0.0,bzB);
        p[-3].boost(0.0,0.0,bzB);
        
  //Genera la fracción de momento z//Scaled momentum fractions
  
  // double cut = 0.0000001;
  // double beta=cut + (1-cut)*(__rng() );
  // double alfa=cut + (1-cut)*(__rng() );
   
 //   z1= pow(x1B,beta);
 //   z2= pow(x2B,alfa);       
           
//C----------------------------------------------------------------------C
//C Suavizado gaussiano
//c     xx(5) --> smearing
/*
 double     SIGMA = SigmaSmearing; // Reccomended value 1.7GeV
 double     Qsolo=sqrt(Q2);
 double     GaussMin = atan(-Qsolo/SIGMA);
 double     GaussMax = atan(Qsolo/SIGMA);
 double     Gauss = GaussMin +__rng() *(GaussMax-GaussMin);
 double     JacGauss = GaussMax-GaussMin;
 double     Qmedido = Qsolo + SIGMA*tan(Gauss);
            JacGauss=JacGauss*SIGMA/(cos(Gauss)*cos(Gauss));
 double     FactorM=exp(-(Qmedido-Qsolo)*(Qmedido-Qsolo)/(2.0*SIGMA*SIGMA))/(SIGMA*sqrt(2.0*__DEB_PI))*JacGauss;
 */
//C----------------------------------------------------------------------C

//-- -------------------------------------------------------------------           
//-- Implementacion de la función de Martin como está en el paper.
//-- -------------------------------------------------------------------           
           
           
           
           
//	cout << "Born PS " << endl;
           
           
           
           
           
           
          // Normal Phsae-Space descomentar para la implementacion
     //      wpsBORN = deltaQ2*fabs(ltaumin)*__DEB_1_8PI/_M_s_;//*(FactorM);
           
           // Factor di smearing ya contiene su Jacobiano
           
           // Espacio de Phases de Martin
           wpsBORN = GsdPs/pow(2.0,8.0)/__DEB_PI/Q2/Ds*(deltaQ2*2.0); //*(FactorM);
           
           wpsBORN =  wpsBORN/2.0 ; // divided by two for identical particles
	   /// std::cout<< "1.0/_M_s_/xx/pow(2.0,8.0)/__DEB_PI/Q2/Ds   : " << 1.0/_M_s_/xx/pow(2.0,8.0)/__DEB_PI/Q2/Ds << std::endl;
           
          // Checked!!!
           
           return wpsBORN;
    
    
    
    
    
    
    
  }
  
  
  
  // Phase space for 2 -> 3 processes
  // All particles outgoing
  // p1+p2+p3+p4+p5
  
  

  double basic_phasespace<_Event220>::operator()(event_type& p, bool &tete)
  {
    double _M_q2min_ = phasespace::mmin*phasespace::mmin;
    double _M_q2max_ = phasespace::mmax*phasespace::mmax;
	  
	  double _M_s_ = opts.sroot*opts.sroot;
  //  unsigned int up = p.upper();
 //   if(up < 0)  throw "unable to generate hadron-hadron -> V+other event";
    
    //----- generate the variable Q2 -----
  
    const random_generator& __rng = *_M_rng;
    
        //   cout << "random desde phasespace  :" << __rng() << endl;
    
    double xa, xb, weight = _M_q2max_-_M_q2min_;
 //   double deltaQ2= _M_q2max_-_M_q2min_;
    double Q2 = weight*__rng() + _M_q2min_;
   // double tau = Q2/_M_s_;
   // double ltau = -std::log(tau);
    
   
    //----- incoming hadrons -----
    double Eh = 0.5*std::sqrt(_M_s_);
    p[hadron(-1)] = lorentzvector<double>(0.0, 0.0,  Eh, Eh);
    p[hadron( 0)] = lorentzvector<double>(0.0, 0.0, -Eh, Eh);
    
    //----- outgoing non-qcd particles (in their c.m. frame) -----
    //event_type::iterator pbegin = p.begin();
    //event_type::iterator pend   = p.end();
    //weight *= _M_psgen -> operator()(Q2, pbegin, pbegin+2);
    
    double taumin = Q2/_M_s_;
    double ltaumin = log(taumin); //ESTO CAPAZ QUE TIENE QUE SER DEFINIDO POSITIVO
    double tau = pow(taumin,__rng());
    double ltau = log(tau);
    double Ecm = sqrt(_M_s_*tau);
    double xi = 1.0-Q2/(Ecm*Ecm);
    double tmp7,tmp2;
    
     double y=0.5*ltau*(1.0-2.0*(tmp7=__rng()));
    
     //----- incoming momentum fractions -----
      xa = sqrt(tau)*exp(y);
      xb = sqrt(tau)*exp(-y);
    
    
     double xabjac = ltau*tau*ltaumin;
    
       double costh = 1.0 - 2.0*(tmp2=__rng());
	  double sinth = sqrt(1.0-costh*costh); 
	  double xphi = __DEB_2PI*__rng(); 
	  double Eq = Ecm*xi*0.5;
	lorentzvector<double> q34;
	p[1]=lorentzvector<double>(Eq*cos(xphi)*sinth,Eq*sin(xphi)*sinth,Eq*costh, Eq);
	q34=lorentzvector<double>(-Eq*cos(xphi)*sinth,-Eq*sin(xphi)*sinth,-Eq*costh,sqrt(Eq*Eq+Q2));
	
	// double bpar=Eq/sqrt(Q2);
	  
	  //Generamos a los difotones en su cm, con masa invariante Q2
	  
	  double Ef=sqrt(Q2)*0.5;
	  double cosx4 = 1.0-2.0*__rng();
	  double sinx4 = sqrt(1.0-cosx4*cosx4);	  
	  double phi = __DEB_2PI*__rng();
	
	  p[-2]=lorentzvector<double>(Ef*cos(phi)*sinx4, Ef*sin(phi)*sinx4, Ef*cosx4, Ef);
	  p[-3]=lorentzvector<double>(-Ef*cos(phi)*sinx4,-Ef*sin(phi)*sinx4,-Ef*cosx4, Ef);
	  
	//  leanlorentz qf10(Ef*cos(phi)*sinx4, Ef*sin(phi)*sinx4, Ef*cosx4, Ef);
	//  leanlorentz qf20(-Ef*cos(phi)*sinx4,-Ef*sin(phi)*sinx4,-Ef*cosx4, Ef);
	
    //Boost difotones; frame of incommings
	  threevector<double> Bvec;
	  threevector<double> bvec(q34); 
 //     bvec= assigm(q34);
      double E34=q34.T();
      Bvec = (1.0/E34)*bvec; 
	  
	  p[-3].boost(Bvec); p[-2].boost(Bvec);
	  
    
    
    double bz = (xa - xb)/(xa + xb);
        p[-3].boost(0, 0, bz);
        p[-2].boost(0, 0, bz);
		p[1].boost(0,0,bz);
		
		
//		double pT2 = p[-2].perp();
//		double pT1 = p[-3].perp();
		
//		if(pT2 > pT1) swap(p[-2],p[-3]); 
		
 //     double Ejet=p[1].T();
 //     double pTJet=p[1].perp();
 //     double phiJet=p[1].phi();
        
      //  pB1 = leanlorentz(pTJet*cos(phiJet),pTJet*sin(phiJet),sqrt(Ejet*Ejet-pTJet*pTJet),Ejet);
        
        
     //   cout << p[1] <<endl;
     //   coutvec(pB1);
        
	   
    // --Chequeamos que pasen los cuts para generar
    
    
   weight *= xabjac*xi*__DEB_1_8PI*__DEB_1_8PI*(1.0/__DEB_2PI);
    
    
    //----- incoming partons -----
    p[-1] = -xa*p[hadron(-1)];
    p[ 0] = -xb*p[hadron( 0)];
   // cout<< p[-1]+p[ 0]-p[ 1]-p[ -2]-p[ -3]<<endl;
    
    
   
   
  
  //  cout << " weight desde born :" << weight<< endl;
    
       return weight;
  }




finalStates::finalStates( const int &size ) 
{
   nParticles = size;
   p = new particle[nParticles];
   for (int i=0; i<nParticles; i++)
      p[i] = particle(0.0);
}


finalStates::finalStates(const finalStates &fs) 
{
   nParticles=fs.nParticles;
   p=new particle[nParticles];
   for (int i=0;i<nParticles;i++)
      p[i]=fs.p[i];
}


finalStates& finalStates::operator=(const finalStates &fs) 
{
   if (nParticles!=fs.nParticles) {
      delete [] p;  
      nParticles=fs.nParticles;
      p=new particle[nParticles];
   }
   for (int i=0;i<nParticles;i++)
      p[i]=fs.p[i];
   return *this;
}


void finalStates::increment(particle &k)
{
   nParticles++;
   particle *ptemp=p;
   p=new particle[nParticles];
   for (int i=0;i<nParticles-1;i++) p[i]=ptemp[i];
   p[nParticles-1]=k;
   delete [] ptemp;
}


