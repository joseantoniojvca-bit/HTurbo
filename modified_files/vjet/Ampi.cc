//   nlo includes                                                                                                                                                                                     
#include "nlo-color.h"
#include "ampq2g1p2.h"
#include "defmacros.h"
#include <math.h>
#include "Inter_params.h"

// DYTurbo includes                                                                                                                                                                                   
#include "settings.h"
#include "coupling.h"
#include "mesq.h"
#include "../include/dyturbo/Li2.hpp"

using namespace polylogarithm;

std::complex<double> cclog(double z);

double ampq2g1p2::AmpIntqgNLO(int p1, int p2, int p3,int p4,int p5) const
{
  // double s12 = S(1,2), s13 = S(1,3), s14 = S(1,4);
  // double s23 = S(2,3), s24 = S(2,4);
  // // cout << "s12+s13+s14+s23+s24 : " << s12+s13+s14+s23+s24 << endl; 
  double s12t =  -S(1,2);
  double s23t =  S(2,3);
  double s45t =  -S(4,5),  s51t = S(5,1);
  double s34t =  -S(3,4);
  // Notice that Prasanna's variables s34 and s45 are not positive.
  // This is in agreement with change sub4 in Amp.m

  // cout << " s45t : " << S(4,5) << endl;
  // cout << " s23t : " << S(2,3) << endl;
  // double s13t = S(1,3); 
  // double s14t = S(1,4);
  // double s24t = S(2,4);
  // cout << "s12t-s13t-s23t : " << s12t-s13t-s23t << endl;
  
  //--- For check with Prasanna and Federico's results
  // s12t = -7.882008462695891;
  // s13t = -1.318636638450604;
  // s23t = -6.243661401410687;
  // double   s14t = -2.734409678428643; 
  // double   s24t = -0.8951075101637909;
  // s51t = -3.828962145816647; 
  // s45t = s12t-s13t-s23t;
  //s34t = (-s12t+s13t+s14t+s23t+s24t);
  
    // s12 =  7.882008462695891;
    // s13 = -1.318636638450604;
    // s23 = -6.243661401410687;
    // s14 = -2.734409678428643; 
    // s24 = -0.8951075101637909;
    // double s34 = S(3,4);
    
    //  cout << " s34 = " <<  S(3,4) << endl;
    //  s34 = -(-s12t+s13t+s14t+s23t+s24t);
    //  s51 = -3.828962145816647; 
    //  s45 = s12-s13-s23;
    //--------------------------------------------------

  double mh = opts.hmass;
  double Gammah = opts.hwidth; // Ok     
  double _VEV_ = sqrt(1./opts.Gf/sqrt(2.0));
  double vmartin = _VEV_/sqrt(2.0);
  double SumQq2 = 11./9.;
  // -- pre factors --
  double e2fac = 4.*__DEB_PI;
  double gs6fac = pow(4.*__DEB_PI,3);
  double Martinfac = 1./(__DEB_PI*vmartin)/(__DEB_PI*vmartin);
  double idenpart = 1./2.;
  double twotimesRealpart = 2.;
  double twogNNLOalphasfac = pow(2.*__DEB_PI,3);
  double factorfromgg = (sqrt(2.)/pow(2.,5)/pow(2.,2));
  double Ave = 1./3./8./2./2.;
  // -----------------
  //std::complex<double> Deno = 1./(std::complex<double>(0,1.)*Gammah*mh + mh*mh + s12 + s13 + s14 + s23 + s24); // Original form of the Deno
  //std::complex<double> Deno = std::complex<double>(1,0.)*(std::complex<double>(0,-1.)*Gammah*mh + mh*mh + s12 + s13 + s14 + s23 + s24)/(Gammah*mh*Gammah*mh+pow(mh*mh + s12 + s13 + s14 + s23 + s24,2)); // -- qT-behaviour OK respect to the CT-- Invariant mass phase inverted respect to the CT -- Form equivalent to precedent line
 //  std::complex<double> Deno = std::complex<double>(0.,1.)*(std::complex<double>(0,-1.)*Gammah*mh + mh*mh - s34)/(Gammah*mh*Gammah*mh+pow(mh*mh - s34,2));
 // Deno = std::complex<double>(0.,1.)*(std::complex<double>(0,1.)*Gammah*mh + mh*mh -s34)/(-Gammah*mh*Gammah*mh+pow(mh*mh -s34,2));
 // std::complex<double> prefac = conj((mesq::CorcheteMartin(s34)))*conj((mesq::SumaMartin(s34)))*((-1. + pow(Nc,2))*SumQq2)*Deno;
 
  double s34tprop = S(3,4);
  std::complex<double> Denot = ((std::complex<double>(0,-1.)*Gammah*mh - mh*mh + s34tprop)/(Gammah*mh*Gammah*mh+pow(-mh*mh + s34tprop,2))); // -- The conjugate is taken over the Higgs part and not over the background
  std::complex<double> prefact = conj((mesq::CorcheteMartin(s34tprop)))*conj((mesq::SumaMartin(s34tprop)))*((-1. + pow(Nc,2))*SumQq2)*Denot;
  prefact = conj((mesq::CorcheteMartin(s34tprop))*(mesq::SumaMartin(s34tprop)))*((-1. + pow(Nc,2))*SumQq2*Denot);
  
   // cout << " Denot : " << Denot << endl;
   // cout << " prefact : " << prefact << endl;
   // cout << " Deno : " << Deno << endl;
   // cout << " prefac : " << prefac << endl;

 double EulerGamma = 0.5772156649015329;
 double mu = 1.;
 
 // std::complex<double> IntqgNLOSempli = (((std::complex<double>(0,-72.)*(s12 + s13 + s14 + s23 + s24)*
 // 					 (pow(s12,5)*(s23 + s24) + pow(s13,3)*s24*pow(s23 + s24,2) + 
 // 					  pow(s12,4)*(s23 + s24)*(3*s13 + 3*s14 + 2*(s23 + s24)) + 
 // 					  s13*s14*s23*(pow(s23 - s24,2)*s24 + s14*pow(s23 + s24,2)) + 
 // 					  pow(s12,3)*(s23 + s24)*
 // 					  (3*pow(s13,2) + 3*pow(s14,2) + 6*s14*s23 + pow(s23,2) + 
 // 					   4*s14*s24 + 4*s23*s24 + pow(s24,2) + 
 // 					   s13*(6*s14 + 4*s23 + 6*s24)) + 
 // 					  pow(s14,2)*s23*(s14*pow(s23 + s24,2) + 
 // 							  s24*(2*pow(s23,2) + s23*s24 + pow(s24,2))) + 
 // 					  pow(s13,2)*s24*(s14*pow(s23 + s24,2) + 
 // 							  s23*(pow(s23,2) + s23*s24 + 2*pow(s24,2))) + 
 // 					  pow(s12,2)*(pow(s13,3)*(s23 + s24) + 
 // 						      pow(s14,3)*(s23 + s24) + 
 // 						      pow(s13,2)*(s23 + s24)*(3*s14 + 2*s23 + 6*s24) + 
 // 						      2*s23*s24*(pow(s23,2) + s23*s24 + pow(s24,2)) + 
 // 						      2*pow(s14,2)*(3*pow(s23,2) + 4*s23*s24 + pow(s24,2)) + 
 // 						      s14*(3*pow(s23,3) + 12.*pow(s23,2)*s24 + 
 // 							   8*s23*pow(s24,2) + pow(s24,3)) + 
 // 						      s13*(pow(s23,3) + 8.*pow(s23,2)*s24 + 
 // 							   12.*s23*pow(s24,2) + 3*pow(s24,3) + 
 // 							   3*pow(s14,2)*(s23 + s24) + 8*s14*pow(s23 + s24,2))) + 
 // 					  s12*(2*pow(s13,3)*s24*(s23 + s24) + 
 // 					       pow(s13,2)*(2*s14*
 // 							   (pow(s23,2) + 3*s23*s24 + 2*pow(s24,2)) + 
 // 							   s24*(4*pow(s23,2) + 9*s23*s24 + 3*pow(s24,2))) + 
 // 					       s14*s23*(2*pow(s14,2)*(s23 + s24) + 
 // 							s24*(4*pow(s23,2) + 3*s23*s24 + pow(s24,2)) + 
 // 							s14*(3*pow(s23,2) + 9*s23*s24 + 4*pow(s24,2))) + 
 // 					       s13*(2*pow(s14,2)*
 // 						    (2*pow(s23,2) + 3*s23*s24 + pow(s24,2)) + 
 // 						    s23*s24*(pow(s23,2) + 3*s23*s24 + 4*pow(s24,2)) + 
 // 						    2*s14*(pow(s23,3) + 5*pow(s23,2)*s24 + 
 // 							   5*s23*pow(s24,2) + pow(s24,3))))))/
 // 					(s23*(s12 + s13 + s14 + s23)*s24*(s12 + s13 + s14 + s24)) + 
 // 					(std::complex<double>(0,72)*(s12 + s13 + s14)*
 // 					 (2*pow(s12,7) + pow(s13,5)*
 // 					  (pow(s23,2) - s23*s24 + 4*pow(s24,2)) + 
 // 					  pow(s12,6)*(10*s13 + 10*s14 + 7*(s23 + s24)) + 
 // 					  pow(s12,5)*(20*pow(s13,2) + 20*pow(s14,2) + 33*s14*s23 + 
 // 						      10*pow(s23,2) + 27*s14*s24 + 16*s23*s24 + 
 // 						      10*pow(s24,2) + s13*(40*s14 + 27*s23 + 33*s24)) + 
 // 					  pow(s13,4)*(3*pow(s23,3) - 2*pow(s23,2)*s24 + 
 // 						      s23*pow(s24,2) + 6*pow(s24,3) + 
 // 						      s14*(2*pow(s23,2) - 17*s23*s24 + 11*pow(s24,2))) + 
 // 					  pow(s13,3)*(3*pow(s23,4) - pow(s23,3)*s24 - 
 // 						      9*pow(s23,2)*pow(s24,2) + 3*pow(s24,4) + 
 // 						      2*pow(s14,2)*(2*pow(s23,2) - 23*s23*s24 + 
 // 								    5*pow(s24,2)) + 
 // 						      s14*(5*pow(s23,3) - 35*pow(s23,2)*s24 - 
 // 							   29*s23*pow(s24,2) + 11*pow(s24,3))) + 
 // 					  pow(s13,2)*(pow(s23,5) - 8*pow(s23,3)*pow(s24,2) - 
 // 						      9*pow(s23,2)*pow(s24,3) - s23*pow(s24,4) + 
 // 						      pow(s24,5) + 2*pow(s14,3)*
 // 						      (5*pow(s23,2) - 23*s23*s24 + 2*pow(s24,2)) + 
 // 						      7*pow(s14,2)*(pow(s23,3) - 9*pow(s23,2)*s24 - 
 // 								    9*s23*pow(s24,2) + pow(s24,3)) + 
 // 						      s14*(4*pow(s23,4) - 24*pow(s23,3)*s24 - 
 // 							   61*pow(s23,2)*pow(s24,2) - 23*s23*pow(s24,3) + 
 // 							   4*pow(s24,4))) + 
 // 					  s14*(-2*pow(s23,2)*pow(s24,2)*pow(s23 + s24,2) + 
 // 					       pow(s14,4)*(4*pow(s23,2) - s23*s24 + pow(s24,2)) + 
 // 					       pow(s14,3)*(6*pow(s23,3) + pow(s23,2)*s24 - 
 // 							   2*s23*pow(s24,2) + 3*pow(s24,3)) + 
 // 					       pow(s14,2)*(3*pow(s23,4) - 
 // 							   9*pow(s23,2)*pow(s24,2) - s23*pow(s24,3) + 
 // 							   3*pow(s24,4)) + 
 // 					       s14*(pow(s23,5) - pow(s23,4)*s24 - 
 // 						    9*pow(s23,3)*pow(s24,2) - 
 // 						    8*pow(s23,2)*pow(s24,3) + pow(s24,5))) + 
 // 					  s13*(-2*pow(s23,2)*pow(s24,2)*pow(s23 + s24,2) + 
 // 					       pow(s14,4)*(11*pow(s23,2) - 17*s23*s24 + 
 // 							   2*pow(s24,2)) + 
 // 					       pow(s14,3)*(11*pow(s23,3) - 29*pow(s23,2)*s24 - 
 // 							   35*s23*pow(s24,2) + 5*pow(s24,3)) + 
 // 					       pow(s14,2)*(4*pow(s23,4) - 23*pow(s23,3)*s24 - 
 // 							   61*pow(s23,2)*pow(s24,2) - 24*s23*pow(s24,3) + 
 // 							   4*pow(s24,4)) + 
 // 					       s14*(pow(s23,5) - 6*pow(s23,4)*s24 - 
 // 						    27*pow(s23,3)*pow(s24,2) - 
 // 						    27*pow(s23,2)*pow(s24,3) - 6*s23*pow(s24,4) + 
 // 						    pow(s24,5))) + 
 // 					  pow(s12,4)*(20*pow(s13,3) + 20*pow(s14,3) + 
 // 						      8*pow(s23,3) + 11*pow(s23,2)*s24 + 
 // 						      11*s23*pow(s24,2) + 8*pow(s24,3) + 
 // 						      pow(s14,2)*(62*s23 + 38*s24) + 
 // 						      pow(s13,2)*(60*s14 + 38*s23 + 62*s24) + 
 // 						      s14*(43*pow(s23,2) + 53*s23*s24 + 28*pow(s24,2)) + 
 // 						      s13*(60*pow(s14,2) + 28*pow(s23,2) + 53*s23*s24 + 
 // 							   43*pow(s24,2) + 100*s14*(s23 + s24))) + 
 // 					  pow(s12,3)*(10*pow(s13,4) + 10*pow(s14,4) + 
 // 						      pow(s14,3)*(58*s23 + 22*s24) + 
 // 						      pow(s13,3)*(40*s14 + 22*s23 + 58*s24) + 
 // 						      pow(s14,2)*(73*pow(s23,2) + 62*s23*s24 + 
 // 								  25*pow(s24,2)) + 
 // 						      s14*(29*pow(s23,3) + 29*pow(s23,2)*s24 + 
 // 							   14*s23*pow(s24,2) + 16*pow(s24,3)) + 
 // 						      2*(2*pow(s23,4) + pow(s23,3)*s24 - 
 // 							 3*pow(s23,2)*pow(s24,2) + s23*pow(s24,3) + 
 // 							 2*pow(s24,4)) + 
 // 						      pow(s13,2)*(60*pow(s14,2) + 25*pow(s23,2) + 
 // 								  62*s23*s24 + 73*pow(s24,2) + 6*s14*(17*s23 + 23*s24)) \
 // 						      + s13*(40*pow(s14,3) + 16*pow(s23,3) + 14*pow(s23,2)*s24 + 
 // 							     29*s23*pow(s24,2) + 29*pow(s24,3) + 
 // 							     6*pow(s14,2)*(23*s23 + 17*s24) + 
 // 							     4*s14*(23*pow(s23,2) + 28*s23*s24 + 23*pow(s24,2)))) \
 // 					  + pow(s12,2)*(2*pow(s13,5) + 2*pow(s14,5) + pow(s23,5) - 
 // 							9*pow(s23,3)*pow(s24,2) - 9*pow(s23,2)*pow(s24,3) + 
 // 							pow(s24,5) + 3*pow(s14,4)*(9*s23 + s24) + 
 // 							pow(s14,3)*(61*pow(s23,2) + 28*s23*s24 + 
 // 								    7*pow(s24,2)) + 
 // 							pow(s13,3)*(20*pow(s14,2) + 36*s14*s23 + 
 // 								    7*pow(s23,2) + 84*s14*s24 + 28*s23*s24 + 
 // 								    61*pow(s24,2)) + 
 // 							pow(s14,2)*(40*pow(s23,3) + 26*pow(s23,2)*s24 - 
 // 								    7*s23*pow(s24,2) + 11*pow(s24,3)) + 
 // 							s14*(11*pow(s23,4) + 2*pow(s23,3)*s24 - 
 // 							     28*pow(s23,2)*pow(s24,2) - 9*s23*pow(s24,3) + 
 // 							     6*pow(s24,4)) + 
 // 							pow(s13,4)*(10*s14 + 3*(s23 + 9*s24)) + 
 // 							pow(s13,2)*(20*pow(s14,3) + 11*pow(s23,3) - 
 // 								    7*pow(s23,2)*s24 + 26*s23*pow(s24,2) + 
 // 								    40*pow(s24,3) + 90*pow(s14,2)*(s23 + s24) + 
 // 								    3*s14*(19*pow(s23,2) + 16*s23*s24 + 37*pow(s24,2))) \
 // 							+ s13*(10*pow(s14,4) + 6*pow(s23,4) - 9*pow(s23,3)*s24 - 
 // 							       28*pow(s23,2)*pow(s24,2) + 2*s23*pow(s24,3) + 
 // 							       11*pow(s24,4) + 12.*pow(s14,3)*(7*s23 + 3*s24) + 
 // 							       3*pow(s14,2)*
 // 							       (37*pow(s23,2) + 16*s23*s24 + 19*pow(s24,2)) + 
 // 							       s14*(41*pow(s23,3) - 11*pow(s23,2)*s24 - 
 // 								    11*s23*pow(s24,2) + 41*pow(s24,3)))) + 
 // 					  s12*(-(pow(s13,5)*(s23 - 5*s24)) + 
 // 					       pow(s14,5)*(5*s23 - s24) - 
 // 					       2*pow(s23,2)*pow(s24,2)*pow(s23 + s24,2) + 
 // 					       pow(s14,4)*(25*pow(s23,2) + 2*s23*s24 + pow(s24,2)) + 
 // 					       pow(s14,3)*(25*pow(s23,3) + 9*pow(s23,2)*s24 - 
 // 							   12.*s23*pow(s24,2) + 6*pow(s24,3)) + 
 // 					       pow(s14,2)*(10*pow(s23,4) - 
 // 							   31*pow(s23,2)*pow(s24,2) - 12.*s23*pow(s24,3) + 
 // 							   5*pow(s24,4)) + 
 // 					       s14*(2*pow(s23,5) - pow(s23,4)*s24 - 
 // 						    19*pow(s23,3)*pow(s24,2) - 
 // 						    20*pow(s23,2)*pow(s24,3) - 3*s23*pow(s24,4) + 
 // 						    pow(s24,5)) + 
 // 					       pow(s13,4)*(pow(s23,2) + 2*s23*s24 + 25*pow(s24,2) + 
 // 							   s14*(s23 + 19*s24)) + 
 // 					       pow(s13,3)*(6*pow(s23,3) - 12.*pow(s23,2)*s24 + 
 // 							   9*s23*pow(s24,2) + 25*pow(s24,3) + 
 // 							   2*pow(s14,2)*(7*s23 + 13*s24) + 
 // 							   2*s14*(5*pow(s23,2) - 14*s23*s24 + 29*pow(s24,2))) + 
 // 					       pow(s13,2)*(5*pow(s23,4) - 12.*pow(s23,3)*s24 - 
 // 							   31*pow(s23,2)*pow(s24,2) + 10*pow(s24,4) + 
 // 							   2*pow(s14,3)*(13*s23 + 7*s24) + 
 // 							   6*pow(s14,2)*
 // 							   (7*pow(s23,2) - 10*s23*s24 + 7*pow(s24,2)) + 
 // 							   s14*(17*pow(s23,3) - 75*pow(s23,2)*s24 - 
 // 								54*s23*pow(s24,2) + 36*pow(s24,3))) + 
 // 					       s13*(pow(s23,5) - 3*pow(s23,4)*s24 - 
 // 						    20*pow(s23,3)*pow(s24,2) - 
 // 						    19*pow(s23,2)*pow(s24,3) - s23*pow(s24,4) + 
 // 						    2*pow(s24,5) + pow(s14,4)*(19*s23 + s24) + 
 // 						    2*pow(s14,3)*
 // 						    (29*pow(s23,2) - 14*s23*s24 + 5*pow(s24,2)) + 
 // 						    pow(s14,2)*(36*pow(s23,3) - 54*pow(s23,2)*s24 - 
 // 								75*s23*pow(s24,2) + 17*pow(s24,3)) + 
 // 						    2*s14*(5*pow(s23,4) - 17*pow(s23,3)*s24 - 
 // 							   48*pow(s23,2)*pow(s24,2) - 17*s23*pow(s24,3) + 
 // 							   5*pow(s24,4)))))*log(s12 + s13 + s14))/
 // 					(pow(s12 + s13 + s14 + s23,2)*pow(s12 + s13 + s14 + s24,2)*
 // 					 (s12 + s13 + s14 + s23 + s24)) - 
 // 					(std::complex<double>(0,72)*(s12 + s13 + s14)*
 // 					 (pow(s12,5) + pow(s12,4)*(3*s13 + 3*s14 + 3*s23 + 2*s24) + 
 // 					  pow(s13,3)*(pow(s23,2) + 2*pow(s24,2)) + 
 // 					  pow(s12,3)*(3*pow(s13,2) + 3*pow(s14,2) + 8*s14*s23 + 
 // 						      4*pow(s23,2) + 3*s14*s24 + 2*s23*s24 + pow(s24,2) + 
 // 						      s13*(6*s14 + 6*s23 + 7*s24)) + 
 // 					  pow(s13,2)*(pow(s23,2)*(3*s23 + s24) + 
 // 						      s14*(pow(s23,2) - 7*s23*s24 + pow(s24,2))) + 
 // 					  pow(s12,2)*(pow(s13,3) + pow(s14,3) + 
 // 						      7*pow(s14,2)*s23 + pow(s13,2)*(3*s14 + 3*s23 + 8*s24) + 
 // 						      s23*(4*pow(s23,2) - s23*s24 - 3*pow(s24,2)) + 
 // 						      s14*(8*pow(s23,2) - s23*s24 - pow(s24,2)) + 
 // 						      s13*(3*pow(s14,2) + 10*s14*s23 + 6*pow(s23,2) + 
 // 							   8*s14*s24 + 6*s23*s24 + 4*pow(s24,2))) + 
 // 					  s23*(pow(s14,3)*(2*s23 - s24) + 
 // 					       s23*(s23 - s24)*pow(s23 + s24,2) + 
 // 					       pow(s14,2)*(2*pow(s23,2) - 3*s23*s24 - 
 // 							   3*pow(s24,2)) + 
 // 					       s14*(pow(s23,3) - 2*pow(s23,2)*s24 - 
 // 						    5*s23*pow(s24,2) - 2*pow(s24,3))) + 
 // 					  s13*(pow(s14,2)*(2*pow(s23,2) - 8*s23*s24 - 
 // 							   pow(s24,2)) + 
 // 					       s23*(3*pow(s23,3) + 2*pow(s23,2)*s24 - 
 // 						    3*s23*pow(s24,2) - 2*pow(s24,3)) + 
 // 					       s14*(2*pow(s23,3) - 9*pow(s23,2)*s24 - 
 // 						    8*s23*pow(s24,2) - pow(s24,3))) + 
 // 					  s12*(pow(s14,3)*(2*s23 - s24) + 3*pow(s13,3)*s24 + 
 // 					       pow(s14,2)*(6*pow(s23,2) - 4*s23*s24 - 2*pow(s24,2)) + 
 // 					       pow(s13,2)*(2*s14*s23 + 3*pow(s23,2) + 5*s14*s24 + 
 // 							   4*s23*s24 + 5*pow(s24,2)) + 
 // 					       s23*(3*pow(s23,3) - 5*s23*pow(s24,2) - 2*pow(s24,3)) + 
 // 					       s14*(4*pow(s23,3) - 6*pow(s23,2)*s24 - 
 // 						    9*s23*pow(s24,2) - pow(s24,3)) + 
 // 					       s13*(6*s14*s23*(s23 - s24) + pow(s14,2)*(4*s23 + s24) + 
 // 						    s23*(6*pow(s23,2) + s23*s24 - 3*pow(s24,2)))))*
 // 					 log(-s23))/
 // 					(pow(s12 + s13 + s14 + s23,2)*(s12 + s13 + s14 + s23 + s24)) - 
 // 					(std::complex<double>(0,72)*(s12 + s13 + s14)*
 // 					 (pow(s12,5) + 2*pow(s14,3)*pow(s23,2) - 
 // 					  2*s14*pow(s23,3)*s24 + pow(s14,3)*pow(s24,2) + 
 // 					  pow(s14,2)*s23*pow(s24,2) - 
 // 					  3*s14*pow(s23,2)*pow(s24,2) - pow(s23,3)*pow(s24,2) + 
 // 					  3*pow(s14,2)*pow(s24,3) + 2*s14*s23*pow(s24,3) - 
 // 					  pow(s23,2)*pow(s24,3) + 3*s14*pow(s24,4) + 
 // 					  s23*pow(s24,4) + pow(s24,5) + 
 // 					  pow(s13,3)*s24*(-s23 + 2*s24) + 
 // 					  pow(s12,4)*(3*s13 + 3*s14 + 2*s23 + 3*s24) + 
 // 					  pow(s12,3)*(3*pow(s13,2) + 3*pow(s14,2) + 7*s14*s23 + 
 // 						      pow(s23,2) + 6*s14*s24 + 2*s23*s24 + 4*pow(s24,2) + 
 // 						      s13*(6*s14 + 3*s23 + 8*s24)) - 
 // 					  pow(s13,2)*(s24*(3*pow(s23,2) + 3*s23*s24 - 
 // 							   2*pow(s24,2)) + 
 // 						      s14*(pow(s23,2) + 8*s23*s24 - 2*pow(s24,2))) + 
 // 					  pow(s12,2)*(pow(s13,3) + pow(s14,3) + 
 // 						      pow(s14,2)*(8*s23 + 3*s24) + 
 // 						      pow(s13,2)*(3*s14 + 7*s24) + 
 // 						      s24*(-3*pow(s23,2) - s23*s24 + 4*pow(s24,2)) + 
 // 						      s14*(4*pow(s23,2) + 6*s23*s24 + 6*pow(s24,2)) + 
 // 						      s13*(3*pow(s14,2) + 8*s14*s23 - pow(s23,2) + 
 // 							   10*s14*s24 - s23*s24 + 8*pow(s24,2))) + 
 // 					  s13*(pow(s14,2)*(pow(s23,2) - 7*s23*s24 + pow(s24,2)) - 
 // 					       s14*(pow(s23,3) + 8*pow(s23,2)*s24 + 
 // 						    9*s23*pow(s24,2) - 2*pow(s24,3)) + 
 // 					       s24*(-2*pow(s23,3) - 5*pow(s23,2)*s24 - 
 // 						    2*s23*pow(s24,2) + pow(s24,3))) + 
 // 					  s12*(3*pow(s14,3)*s23 - pow(s13,3)*(s23 - 2*s24) + 
 // 					       pow(s14,2)*(5*pow(s23,2) + 4*s23*s24 + 3*pow(s24,2)) + 
 // 					       pow(s13,2)*(s14*s23 - 2*pow(s23,2) + 4*s14*s24 - 
 // 							   4*s23*s24 + 6*pow(s24,2)) + 
 // 					       s14*s24*(-3*pow(s23,2) + s23*s24 + 6*pow(s24,2)) + 
 // 					       s24*(-2*pow(s23,3) - 5*pow(s23,2)*s24 + 
 // 						    3*pow(s24,3)) + 
 // 					       s13*(-pow(s23,3) - 9*pow(s23,2)*s24 - 
 // 						    6*s23*pow(s24,2) + 4*pow(s24,3) + 
 // 						    6*s14*s24*(-s23 + s24) + pow(s14,2)*(5*s23 + 2*s24))))*
 // 					 log(-s24))/
 // 					(pow(s12 + s13 + s14 + s24,2)*(s12 + s13 + s14 + s23 + s24)) + 
 // 					(12.*(pow(__DEB_PI,2)*(96*(6*s12 + 6*s13 + 6*s14 - s23 - s24)*
 // 							      pow(s12 + s13 + s14 + s23 + s24,2) + 
 // 							      std::complex<double>(0,1)*(s12 + s13 + s14)*
 // 							      (2*pow(s12,4) + 2*pow(s14,2)*pow(s23,2) + 
 // 							       s14*pow(s23,3) + pow(s23,4) + pow(s14,2)*s23*s24 + 
 // 							       2*s14*pow(s23,2)*s24 + 2*pow(s23,3)*s24 + 
 // 							       pow(s14,2)*pow(s24,2) + 3*s14*s23*pow(s24,2) + 
 // 							       2*pow(s23,2)*pow(s24,2) + 2*s14*pow(s24,3) + 
 // 							       2*s23*pow(s24,3) + pow(s24,4) + 
 // 							       pow(s13,2)*(pow(s23,2) + s23*s24 + 2*pow(s24,2)) + 
 // 							       s13*(2*pow(s23,3) + s14*pow(s23 - s24,2) + 
 // 								    3*pow(s23,2)*s24 + 2*s23*pow(s24,2) + pow(s24,3)) \
 // 							       + pow(s12,3)*(4*s13 + 4*s14 + 5*(s23 + s24)) + 
 // 							       pow(s12,2)*(2*pow(s13,2) + 2*pow(s14,2) + 
 // 									   8*s14*s23 + 5*pow(s23,2) + 6*s14*s24 + 8*s23*s24 + 
 // 									   5*pow(s24,2) + s13*(4*s14 + 6*s23 + 8*s24)) + 
 // 							       s12*(3*pow(s23,3) + 5*pow(s23,2)*s24 + 
 // 								    5*s23*pow(s24,2) + 3*pow(s24,3) + 
 // 								    pow(s14,2)*(3*s23 + s24) + 
 // 								    pow(s13,2)*(s23 + 3*s24) + 
 // 								    s14*(5*pow(s23,2) + 7*s23*s24 + 4*pow(s24,2)) + 
 // 								    s13*(4*pow(s23,2) + 7*s23*s24 + 5*pow(s24,2) + 
 // 									 4*s14*(s23 + s24))))) + 
 // 					     std::complex<double>(0,3)*(s12 + s13 + s14)*
 // 					     (2*pow(s12,4) + 2*pow(s14,2)*pow(s23,2) + 
 // 					      s14*pow(s23,3) + pow(s23,4) + pow(s14,2)*s23*s24 + 
 // 					      2*s14*pow(s23,2)*s24 + 2*pow(s23,3)*s24 + 
 // 					      pow(s14,2)*pow(s24,2) + 3*s14*s23*pow(s24,2) + 
 // 					      2*pow(s23,2)*pow(s24,2) + 2*s14*pow(s24,3) + 
 // 					      2*s23*pow(s24,3) + pow(s24,4) + 
 // 					      pow(s13,2)*(pow(s23,2) + s23*s24 + 2*pow(s24,2)) + 
 // 					      s13*(2*pow(s23,3) + s14*pow(s23 - s24,2) + 
 // 						   3*pow(s23,2)*s24 + 2*s23*pow(s24,2) + pow(s24,3)) + 
 // 					      pow(s12,3)*(4*s13 + 4*s14 + 5*(s23 + s24)) + 
 // 					      pow(s12,2)*(2*pow(s13,2) + 2*pow(s14,2) + 8*s14*s23 + 
 // 							  5*pow(s23,2) + 6*s14*s24 + 8*s23*s24 + 
 // 							  5*pow(s24,2) + s13*(4*s14 + 6*s23 + 8*s24)) + 
 // 					      s12*(3*pow(s23,3) + 5*pow(s23,2)*s24 + 
 // 						   5*s23*pow(s24,2) + 3*pow(s24,3) + 
 // 						   pow(s14,2)*(3*s23 + s24) + pow(s13,2)*(s23 + 3*s24) + 
 // 						   s14*(5*pow(s23,2) + 7*s23*s24 + 4*pow(s24,2)) + 
 // 						   s13*(4*pow(s23,2) + 7*s23*s24 + 5*pow(s24,2) + 
 // 							4*s14*(s23 + s24))))*pow(log(-s23),2) - 
 // 					     std::complex<double>(0,6)*(s12 + s13 + s14)*
 // 					     (2*pow(s12,4) + 2*pow(s14,2)*pow(s23,2) + 
 // 					      s14*pow(s23,3) + pow(s23,4) + pow(s14,2)*s23*s24 + 
 // 					      2*s14*pow(s23,2)*s24 + 2*pow(s23,3)*s24 + 
 // 					      pow(s14,2)*pow(s24,2) + 3*s14*s23*pow(s24,2) + 
 // 					      2*pow(s23,2)*pow(s24,2) + 2*s14*pow(s24,3) + 
 // 					      2*s23*pow(s24,3) + pow(s24,4) + 
 // 					      pow(s13,2)*(pow(s23,2) + s23*s24 + 2*pow(s24,2)) + 
 // 					      s13*(2*pow(s23,3) + s14*pow(s23 - s24,2) + 
 // 						   3*pow(s23,2)*s24 + 2*s23*pow(s24,2) + pow(s24,3)) + 
 // 					      pow(s12,3)*(4*s13 + 4*s14 + 5*(s23 + s24)) + 
 // 					      pow(s12,2)*(2*pow(s13,2) + 2*pow(s14,2) + 8*s14*s23 + 
 // 							  5*pow(s23,2) + 6*s14*s24 + 8*s23*s24 + 
 // 							  5*pow(s24,2) + s13*(4*s14 + 6*s23 + 8*s24)) + 
 // 					      s12*(3*pow(s23,3) + 5*pow(s23,2)*s24 + 
 // 						   5*s23*pow(s24,2) + 3*pow(s24,3) + 
 // 						   pow(s14,2)*(3*s23 + s24) + pow(s13,2)*(s23 + 3*s24) + 
 // 						   s14*(5*pow(s23,2) + 7*s23*s24 + 4*pow(s24,2)) + 
 // 						   s13*(4*pow(s23,2) + 7*s23*s24 + 5*pow(s24,2) + 
 // 							4*s14*(s23 + s24))))*log(-s23)*log(-s24) + 
 // 					     std::complex<double>(0,3)*(s12 + s13 + s14)*
 // 					     (2*pow(s12,4) + 2*pow(s14,2)*pow(s23,2) + 
 // 					      s14*pow(s23,3) + pow(s23,4) + pow(s14,2)*s23*s24 + 
 // 					      2*s14*pow(s23,2)*s24 + 2*pow(s23,3)*s24 + 
 // 					      pow(s14,2)*pow(s24,2) + 3*s14*s23*pow(s24,2) + 
 // 					      2*pow(s23,2)*pow(s24,2) + 2*s14*pow(s24,3) + 
 // 					      2*s23*pow(s24,3) + pow(s24,4) + 
 // 					      pow(s13,2)*(pow(s23,2) + s23*s24 + 2*pow(s24,2)) + 
 // 					      s13*(2*pow(s23,3) + s14*pow(s23 - s24,2) + 
 // 						   3*pow(s23,2)*s24 + 2*s23*pow(s24,2) + pow(s24,3)) + 
 // 					      pow(s12,3)*(4*s13 + 4*s14 + 5*(s23 + s24)) + 
 // 					      pow(s12,2)*(2*pow(s13,2) + 2*pow(s14,2) + 8*s14*s23 + 
 // 							  5*pow(s23,2) + 6*s14*s24 + 8*s23*s24 + 
 // 							  5*pow(s24,2) + s13*(4*s14 + 6*s23 + 8*s24)) + 
 // 					      s12*(3*pow(s23,3) + 5*pow(s23,2)*s24 + 
 // 						   5*s23*pow(s24,2) + 3*pow(s24,3) + 
 // 						   pow(s14,2)*(3*s23 + s24) + pow(s13,2)*(s23 + 3*s24) + 
 // 						   s14*(5*pow(s23,2) + 7*s23*s24 + 4*pow(s24,2)) + 
 // 						   s13*(4*pow(s23,2) + 7*s23*s24 + 5*pow(s24,2) + 
 // 							4*s14*(s23 + s24))))*pow(log(-s24),2) + 
 // 					     std::complex<double>(0,6)*(s12 + s13 + s14)*
 // 					     (2*pow(s12,4) + 2*pow(s14,2)*pow(s23,2) + 
 // 					      s14*pow(s23,3) + pow(s23,4) + pow(s14,2)*s23*s24 + 
 // 					      2*s14*pow(s23,2)*s24 + 2*pow(s23,3)*s24 + 
 // 					      pow(s14,2)*pow(s24,2) + 3*s14*s23*pow(s24,2) + 
 // 					      2*pow(s23,2)*pow(s24,2) + 2*s14*pow(s24,3) + 
 // 					      2*s23*pow(s24,3) + pow(s24,4) + 
 // 					      pow(s13,2)*(pow(s23,2) + s23*s24 + 2*pow(s24,2)) + 
 // 					      s13*(2*pow(s23,3) + s14*pow(s23 - s24,2) + 
 // 						   3*pow(s23,2)*s24 + 2*s23*pow(s24,2) + pow(s24,3)) + 
 // 					      pow(s12,3)*(4*s13 + 4*s14 + 5*(s23 + s24)) + 
 // 					      pow(s12,2)*(2*pow(s13,2) + 2*pow(s14,2) + 8*s14*s23 + 
 // 							  5*pow(s23,2) + 6*s14*s24 + 8*s23*s24 + 
 // 							  5*pow(s24,2) + s13*(4*s14 + 6*s23 + 8*s24)) + 
 // 					      s12*(3*pow(s23,3) + 5*pow(s23,2)*s24 + 
 // 						   5*s23*pow(s24,2) + 3*pow(s24,3) + 
 // 						   pow(s14,2)*(3*s23 + s24) + pow(s13,2)*(s23 + 3*s24) + 
 // 						   s14*(5*pow(s23,2) + 7*s23*s24 + 4*pow(s24,2)) + 
 // 						   s13*(4*pow(s23,2) + 7*s23*s24 + 5*pow(s24,2) + 
 // 							4*s14*(s23 + s24))))*
 // 					     Li2((s12 + s13 + s14 + s23)/s23) + 
 // 					     std::complex<double>(0,6)*(s12 + s13 + s14)*
 // 					     (2*pow(s12,4) + 2*pow(s14,2)*pow(s23,2) + 
 // 					      s14*pow(s23,3) + pow(s23,4) + pow(s14,2)*s23*s24 + 
 // 					      2*s14*pow(s23,2)*s24 + 2*pow(s23,3)*s24 + 
 // 					      pow(s14,2)*pow(s24,2) + 3*s14*s23*pow(s24,2) + 
 // 					      2*pow(s23,2)*pow(s24,2) + 2*s14*pow(s24,3) + 
 // 					      2*s23*pow(s24,3) + pow(s24,4) + 
 // 					      pow(s13,2)*(pow(s23,2) + s23*s24 + 2*pow(s24,2)) + 
 // 					      s13*(2*pow(s23,3) + s14*pow(s23 - s24,2) + 
 // 						   3*pow(s23,2)*s24 + 2*s23*pow(s24,2) + pow(s24,3)) + 
 // 					      pow(s12,3)*(4*s13 + 4*s14 + 5*(s23 + s24)) + 
 // 					      pow(s12,2)*(2*pow(s13,2) + 2*pow(s14,2) + 8*s14*s23 + 
 // 							  5*pow(s23,2) + 6*s14*s24 + 8*s23*s24 + 5*pow(s24,2) + 
 // 							  s13*(4*s14 + 6*s23 + 8*s24)) + 
 // 					      s12*(3*pow(s23,3) + 5*pow(s23,2)*s24 + 
 // 						   5*s23*pow(s24,2) + 3*pow(s24,3) + 
 // 						   pow(s14,2)*(3*s23 + s24) + pow(s13,2)*(s23 + 3*s24) + 
 // 						   s14*(5*pow(s23,2) + 7*s23*s24 + 4*pow(s24,2)) + 
 // 						   s13*(4*pow(s23,2) + 7*s23*s24 + 5*pow(s24,2) + 
 // 							4*s14*(s23 + s24))))*
 // 					     Li2((s12 + s13 + s14 + s24)/s24)))/
 // 					 pow(s12 + s13 + s14 + s23 + s24,2))/pow(__DEB_PI,2))/(576.*pow(s12 + s13 + s14,2));;


 std::complex<double> IntqgNLOulti = 2.*cclog(s23t)*pow(s34t,-1)*pow(s23t - s51t,-2)*pow(s51t,-1)*
   (s34t*(2*s34t*s45t - 3*s34t*s51t + 2*s45t*s51t)*pow(s23t,2) + 2*pow(s23t,3)*pow(s34t,2) + 
    s23t*(s34t + s51t)*(s34t*s45t*(2*s45t - s51t) + (2*s45t + s51t)*pow(s34t,2) + 
		     2*s51t*pow(s45t,2)) + (2*s23t - s34t - 2*s51t)*pow(s12t,2)*pow(s23t - s51t,2) + 
    s12t*((5*s34t + 4*s45t)*s51t*pow(s23t,2) - 2*s34t*pow(s23t,3) + 
         (s34t + 4*s45t)*(s34t + s51t)*pow(s51t,2) - 
         s23t*(4*s34t*s51t*(s45t + s51t) + (2*s45t + s51t)*pow(s34t,2) + 8*s45t*pow(s51t,2))) + 
    s45t*(s34t*(s45t - s51t) - 2*s45t*s51t)*pow(s34t + s51t,2)) - 
   2.*cclog(-s23t + s34t + s51t)*pow(s23t - s34t,-2)*pow(s34t,-1)*pow(s51t,-1)*
   (s34t*pow(s23t,2)*(2*s34t*s45t - s34t*s51t + 2*s45t*s51t - 4*pow(s34t,2)) + 
    2*pow(s23t,3)*pow(s34t,2) + s23t*(s34t + s51t)*
    (s34t*s45t*(2*s45t + s51t) - 2*s45t*pow(s34t,2) + 2*pow(s34t,3) + 2*s51t*pow(s45t,2)) + 
    pow(s12t,2)*(-((5*s34t + 4*s51t)*pow(s23t,2)) + 2*pow(s23t,3) + 
		2*s23t*(4*s34t*s51t + 2*pow(s34t,2) + pow(s51t,2)) - 
		s34t*(4*s34t*s51t + 2*pow(s34t,2) + 3*pow(s51t,2))) + 
    s12t*(-2*s34t*pow(s23t,3) + 2*s34t*(s34t + s51t)*
	 (3*s45t*s51t + s34t*(2*s45t + s51t) + pow(s34t,2)) + 
         pow(s23t,2)*(3*s34t*s51t + 4*s45t*s51t + 4*pow(s34t,2)) - 
         s23t*(s34t*s51t*(12*s45t + s51t) + 2*(s45t + 3*s51t)*pow(s34t,2) + 2*pow(s34t,3) + 
	      4*s45t*pow(s51t,2))) - s34t*(2*s34t*s45t + pow(s34t,2) + 3*pow(s45t,2))*
    pow(s34t + s51t,2)) - 2.*cclog(s51t)*pow(s23t - s34t,-2)*pow(s34t,-1)*pow(s23t - s51t,-2)*
   pow(s51t,-1)*(-(s34t*(3*s34t - 2*s45t)*s51t*(s34t + s51t)*pow(s23t,3)) + 
		2*s51t*pow(s23t,4)*pow(s34t,2) + 
		pow(s12t,2)*(2*s23t*(s23t - s51t)*s51t + s34t*s51t*(-4*s23t + 3*s51t) + 2*s51t*pow(s34t,2) + 
			    pow(s34t,3))*pow(s23t - s51t,2) - 
		s23t*(s34t + s51t)*(s51t*(3*s45t + 4*s51t)*pow(s34t,3) + (-2*s45t + s51t)*pow(s34t,4) + 
				 2*s51t*pow(s34t,2)*pow(s45t,2) + s34t*s45t*(4*s45t + s51t)*pow(s51t,2) + 
				 2*pow(s45t,2)*pow(s51t,3)) + 
		pow(s23t,2)*(s51t*(-s45t + 7*s51t)*pow(s34t,3) + s51t*pow(s34t,4) + pow(s34t,5) + 
			    s34t*s45t*(4*s45t - s51t)*pow(s51t,2) + 
			    s51t*pow(s34t,2)*(-2*s45t*s51t + 2*pow(s45t,2) + pow(s51t,2)) + 
			    2*pow(s45t,2)*pow(s51t,3)) - 
		s12t*(2*s34t*s51t*pow(s23t,4) - s51t*pow(s23t,3)*
		     (5*s34t*s51t + 4*s45t*s51t + 3*pow(s34t,2)) + 
		     s34t*(s34t + s51t)*(2*s34t*s51t + 6*s45t*s51t + pow(s34t,2))*pow(s51t,2) + 
		     pow(s23t,2)*(s51t*(2*s45t + 9*s51t)*pow(s34t,2) + s51t*pow(s34t,3) + 2*pow(s34t,4) + 
				 2*s34t*(5*s45t + 2*s51t)*pow(s51t,2) + 8*s45t*pow(s51t,3)) + 
		     s23t*(-4*s51t*(s45t + s51t)*pow(s34t,3) + (2*s45t - 3*s51t)*pow(s34t,4) - 
			  2*(3*s45t + 4*s51t)*pow(s34t,2)*pow(s51t,2) - s34t*(16*s45t + s51t)*pow(s51t,3) - 
			  4*s45t*pow(s51t,4))) + s34t*
		(2*s34t*s45t*s51t*(-s45t + s51t) + 3*pow(s45t,2)*pow(s51t,2) + 
		 pow(s34t,2)*(-(s45t*s51t) + pow(s45t,2) + pow(s51t,2)))*pow(s34t + s51t,2)) + 
   2*s34t*pow(s23t,-1)*pow(s23t - s34t,-1)*pow(s23t - s51t,-1)*pow(s23t - s34t - s51t,-1)*
   pow(s51t,-2)*(2*pow(s23t,4)*(s34t*s51t + pow(s34t,2) + pow(s51t,2)) - 
		(s34t + s51t)*pow(s23t,3)*(-4*s34t*s45t + 3*s34t*s51t - 2*s45t*s51t + 2*pow(s34t,2) + 
					4*pow(s51t,2)) + (s23t - s51t)*pow(s12t,2)*
		(s23t*s51t*(s34t + s51t) - 2*(s34t + s51t)*pow(s23t,2) + 2*pow(s23t,3) - 
		 (s34t + s51t)*pow(s51t,2)) - s23t*(2*s45t - s51t)*(s34t*(s45t - s51t) + s45t*s51t)*
		pow(s34t + s51t,2) + pow(s23t,2)*
		(-4*s34t*s45t + 2*s34t*s51t - 3*s45t*s51t + 2*pow(s45t,2) + 2*pow(s51t,2))*
		pow(s34t + s51t,2) - s12t*((s34t + s51t)*(-4*s34t*s45t + 5*s34t*s51t + 2*s45t*s51t)*
					pow(s23t,2) + 2*(2*s34t + s51t)*pow(s23t,4) - 
					s23t*(s34t + 4*s45t - s51t)*(s34t + s51t)*pow(s51t,2) - 
					pow(s23t,3)*(-4*s34t*s45t + 9*s34t*s51t + 4*pow(s34t,2) + 3*pow(s51t,2)) + 
					2*s45t*pow(s51t,2)*pow(s34t + s51t,2)) + s51t*pow(s45t,2)*pow(s34t + s51t,3)) - 
   (pow(s34t,-2)*pow(s51t,-1)*(-(s23t*s34t*(s34t - 2*s45t)*(s34t + s51t)) + 
			     2*pow(s23t,2)*pow(s34t,2) - s12t*
			     (2*s34t*pow(s23t,2) - s23t*(3*s34t*s51t + 4*s45t*s51t + pow(s34t,2)) + 
			      (s34t + s51t)*(2*s34t*s45t + s34t*s51t + 4*s45t*s51t + 2*pow(s34t,2))) + 
			     2*pow(s12t,2)*(s34t*s51t - s23t*(s34t + 2*s51t) + pow(s23t,2) + pow(s34t,2) + 
					   pow(s51t,2)) + (s34t*s45t + pow(s34t,2) + 2*pow(s45t,2))*pow(s34t + s51t,2))*
    (6.*Li2(1 - s51t*pow(s23t,-1)) + 6.*Li2(1 - s51t*pow(-s23t + s34t + s51t,-1)) + 
     pow(__DEB_PI,2) + 3.*pow(cclog(s23t*pow(-s23t + s34t + s51t,-1)),2)))/3.;

 //-- Check Prasanna and Federico's results which are in agreement in the real part
 // cout << "IntqgNLOulti : "  << IntqgNLOulti<< endl;
 // cout << " real(IntqgNLO)/real(IntqgNLOulti) : " << real(IntqgNLOSempli*std::complex<double>(1,0)*prefac)*16*pow(__DEB_PI,2)/real(IntqgNLOulti*prefact) << endl;
 // cout << "IntqgNLOulti : "  << IntqgNLOulti<< endl;
 // cout << "IntqgNLOSempli : "  << IntqgNLOSempli*std::complex<double>(1,0)*16.*pow(__DEB_PI,2)<< endl;
 // cout << " imag(IntqgNLO)/real(IntqgNLOulti) : " << imag(IntqgNLOSempli*std::complex<double>(1,0)*prefac)*16*pow(__DEB_PI,2)/real(IntqgNLOulti*prefact) << endl;
 // cout << " real(IntqgNLO)/imag(IntqgNLOulti) : " << real(IntqgNLOSempli*std::complex<double>(1,0)*prefac)*16*pow(__DEB_PI,2)/imag(IntqgNLOulti*prefact) << endl;
 // cout << " imag(IntqgNLO)/imag(IntqgNLOulti) : " << imag(IntqgNLOSempli*std::complex<double>(1,0)*prefac)*16*pow(__DEB_PI,2)/imag(IntqgNLOulti*prefact) << endl;
 // --------------------------------------------------------------------------------

 return  twotimesRealpart*Ave*factorfromgg*e2fac*gs6fac*Martinfac*idenpart*twogNNLOalphasfac*real(IntqgNLOulti*prefact)/16./pow(__DEB_PI,2);
}

std::complex<double> cclog(double z)
{
  if (z > 0.0) {
    return { std::log(z), 0.0 };
  } else  
    return { std::log(-z), __DEB_PI };
  
  //  return { std::log(std::hypot(std::real(z), std::imag(z))), std::arg(z) };
}
