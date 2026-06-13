//      ampBORN.cpp
//      
//      Copyright 2010 Leandro <leandro@ubuntu>
//      
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//      
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//      
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.


#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


#include "nlo-process.h"
#include "H1.h"
#include "Inter_params.h"
#include "nlo-color.h"


// --- From DYTurbo
#include "mesq.h"
#include "coupling.h"
#include "settings.h"

//fortran interface                                                                             
extern "C" {
// FLpppp
// FSLmmpp
// FLmmpp
// FSLpppp
  std::complex<double> flpppp_ (double & x1,double & y,double & lx1,double & ly);
  std::complex<double> fslmmpp_(double & x1,double & y,double & lx1,double & ly);
  std::complex<double> flmmpp_ (double & x1,double & y,double & lx1,double & ly);
  std::complex<double> fslpppp_(double & x1,double & y,double & lx1,double & ly);
}

double  ampBORN( event_hhc2ph p)
{
  double s13 =-2.0*(p[-1]*p[-3]), s14 =-2.0*(p[-1]*p[-2]);
  //tiene la suma sobre colores el flujo y el commonfac(1),
  // Promediado sobre colores y spines  
  return Na*Nc*(s13/s14+s14/s13)*__DEB_4PI2/Nc2/8.0;  
}
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
// -------------- Interferencia LO --- LEAN
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  

//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
//---------- Fixed order LO for signal and interference  LEAN ---------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Leandro Cieri /////////////////////////////////////////////////////
 
std::complex<double> LO_MggHB( event_hhc2ph p)
{ double s12 = 2.0*(p[-1]*p[0]);
  return mesq::LO_MggH(s12);
}

std::complex<double> LO_MHyyB( event_hhc2ph p)
{ //Este elemento de matriz calcula el decaimiento del Higgs H->yy
  // Falta hacer el modulo cuadrado que se hace en el momento de multiplicar por las PDFs  
  // In the infinit top mass limit
  // Por el momento solo la contribucion top, tau y W  
  double s12 = 2.0*(p[-1]*p[0]);

  return  mesq::LO_MHyy(s12) ;
 }

 double LO_MggyyB( event_hhc2ph p, double & spass) 
 { //Este elemento de matriz calcula solo la contribuciones de helicidad que entran en juego en la interferencia
   // La contribucion del top no está tenida en cuenta
   // Suma sobre el resto de las cargas de los quarks livianos
   double u =-2.0*(p[-1]*p[-3]), t =-2.0*(p[-1]*p[-2]); 
   double s12 = 2.0*(p[-1]*p[0]);
   spass = s12;
   return  mesq::LO_Mggyy(s12,t,u);
 }

 double _amp_MggHyyB( event_hhc2ph p)
 {
   double u =-2.0*(p[-1]*p[-3]), t =-2.0*(p[-1]*p[-2]); 
   double s = 2.0*(p[-1]*p[0]);
   return mesq::_amp_MggHyy( s,  t, u);
 }

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// -- H1 for the interference in the Hard Scheme 08/07/2024
//----------------------------------------------------------------------
//----------------------------------------------------------------------

double H1intHSchm(double costh, double s)
{
  double FormFactor = 3.0*M_PI*M_PI/2.0 +11.0/2.0; //H1g Hard Scheme
  double mt = opts.mtop;
  double _VEV_ = sqrt(1./opts.Gf/sqrt(2.0));
  double vmartin = _VEV_/sqrt(2.0);
  double et =  (2.0/3.0);
  double Mh = opts.hmass;
  double Gh = opts.hwidth; // Ok
  double u  = -1./2.*s*(1.-costh);
  double t  = -1./2.*s*(1.+costh);
  double x1 = t/s;
  double y  = u/s;
  double lx1 = log(-x1);
  double ly  = log(-y);
  if (x1 > 0. || y > 0.) cout << "wrong sign x1 or y in ampBORN.cpp " << endl;
  double mur2 = s; // There is no scale dependence in the Hard functions
  double average = 1./(2.*2.*8.*8.);
  double identical = 1./2.;
  double color = 8.;
  double Sum_eq2 = sqrt(121./81.);
  double factors = (4.*Sum_eq2)*(4.*Sum_eq2); 
  std::complex<double>  AggH  = mesq::LO_MggH( s);
  std::complex<double>  AyyH  = mesq::LO_MHyy( s)/opts.aemmz; // Divided by alpha as in resolve
  std::complex<double>  AyyHt = -(s)/(4.*sqrt(2.)*M_PI*vmartin)*(3.0*et*et*mesq::F1on2(4.0*mt*mt/s));
  std::complex<double> Msig   = -(AggH*AyyH)/(s-Mh*Mh+ std::complex<double> (0.,1.)*Mh*Gh);
  std::complex<double> MsigAA = -(AggH*AyyH)*(s-Mh*Mh) /((s-Mh*Mh)*(s-Mh*Mh) + (Mh*Gh)*(Mh*Gh));
  std::complex<double> MsigBB = -(AggH*AyyH)*std::complex<double> (0.,-1.)*Mh*Gh/((s-Mh*Mh)*(s-Mh*Mh) + (Mh*Gh)*(Mh*Gh));
  std::complex<double> M1sig  = FormFactor*Msig-2.*( -(AggH*AyyHt)/(s-Mh*Mh+std::complex<double> (0.,1.)*Mh*Gh));
  // cout << " C AggH = " << AggH << endl;
  // cout << " C AyyH = " << AyyH << endl;
  // cout << " C  AyyHt = " <<  AyyHt << endl;
  // cout << " C  Msig = " <<  Msig << endl;
  // cout << " C MsigAA = " << MsigAA << endl;
  // cout << " C MsigBB = " << MsigBB << endl;
  // cout << " C M1sig = " << M1sig << endl;
  double  Poles = 1.0/6.0*(9.0*M_PI*M_PI+(33.0-2.0*5.)*log(mur2/s));
  std::complex<double> Mqback_pppp = std::complex<double>(-1. + costh*log((1.+costh)/(1.-costh)) -( 1. + costh*costh )/4.*(pow((log((1.+costh)/(1.-costh))),2.)+M_PI*M_PI),0.);
  std::complex<double> Mqback_mmpp = std::complex<double>(1.,0.);
  std::complex<double> Mqback_mmmm = Mqback_pppp;
  std::complex<double> Mqback_ppmm = Mqback_mmpp;
  std::complex<double> M1back_pppp = Poles*Mqback_pppp+3.*flpppp_(x1,y,lx1,ly)-1.0/3.0*fslpppp_(x1,y,lx1,ly);
  std::complex<double> M1back_mmpp = Poles*Mqback_mmpp+3.*flmmpp_(x1,y,lx1,ly)-1.0/3.0*fslmmpp_(x1,y,lx1,ly);
  std::complex<double> M0back_pppp = Mqback_pppp;
  std::complex<double> M0back_mmpp = Mqback_mmpp;        
  // cout << " Mqback_pppp = " <<  Mqback_pppp << endl;
  // cout << " Mqback_mmpp = " << Mqback_mmpp << endl;
  // cout << " Mqback_mmmm = " << Mqback_mmmm << endl;
  // cout << "  Mqback_ppmm = " <<  Mqback_ppmm << endl;
  // cout << " M1back_pppp = " << M1back_pppp << endl;
  // cout << " M1back_mmpp = " << M1back_mmpp << endl;
  // cout << " M0back_pppp = " << M0back_pppp << endl;
  // cout << "  M0back_mmpp = " <<  M0back_mmpp << endl;
  std::complex<double> M1back0sig_pppp = 2.*(conj(M1back_pppp)*Msig);
  std::complex<double> M1back0sig_mmpp = 2.*(conj(M1back_mmpp)*Msig);
  std::complex<double> M1back0sig = 2. * (M1back0sig_pppp + M1back0sig_mmpp); // times 2 because pppp = mmmm and mmpp = ppmm
  M1back0sig = M1back0sig*average*identical*color*sqrt(factors)*4.*M_PI*M_PI;
  // cout << " M1back0sig_pppp = " << M1back0sig_pppp << endl;
  // cout << " M1back0sig_mmpp = " <<  M1back0sig_mmpp << endl;
  // cout << " M1back0sig = " << M1back0sig << endl;
  // cout << "M1back0sig = " << M1back0sig << endl;
  std::complex<double> M0back1sig_pppp = 2.*(conj(M1sig)*M0back_pppp);
  std::complex<double> M0back1sig_mmpp = 2.*(conj(M1sig)*M0back_mmpp);
  std::complex<double> M0back1sig = 2. * (M0back1sig_pppp + M0back1sig_mmpp); // times because pppp = mmmm and mmpp = ppmm                                                                
  M0back1sig = M0back1sig*average*identical*color*sqrt(factors)*4.*M_PI*M_PI;
  std::complex<double> Minterf = 2.*(conj(Msig)*(Mqback_pppp+Mqback_mmmm+Mqback_ppmm+Mqback_mmpp));
  Minterf = Minterf*average*identical*color*sqrt(factors)*4.*M_PI*M_PI;

  
  std::complex<double> H1gYYint_1back0sig = 1.0/2.0*M1back0sig/real(Minterf);
  std::complex<double> H1gYYint_0back1sig = 1.0/2.0*M0back1sig/real(Minterf);
  std::complex<double> H1gYYintDY = 4./3.*(M_PI*M_PI/2. - 4.);
  std::complex<double> H1gYYint = H1gYYint_1back0sig + H1gYYint_0back1sig-0.*H1gYYintDY; // El H1DY no tiene nada que ver con el canal gg
  // cout << " M1back0sig = " << M1back0sig << endl;
  // cout << " M0back1sig_pppp = " << M0back1sig_pppp << endl;
  // cout << "  M0back1sig_mmpp = " << M0back1sig_mmpp<< endl;
  // cout << "   M0back1sig = " << M0back1sig<< endl;
  // cout << " Minterf = " << Minterf << endl;
  // cout << " H1gYYint_1back0sig = " << H1gYYint_1back0sig<< endl;
  // cout << "  H1gYYint_0back1sig = " << H1gYYint_0back1sig<< endl;
  // cout << "   H1gYYint = " << H1gYYint << endl;
  //  H1gYYint =H1gYYint; // In order to adjust expansion from DYTurbo see the file where it is implemented
  return real(H1gYYint); // Cross checked with resolve implementation -- Lean -- 09/07/2024
}
