//  Copyright (C) 2024 Leandro Cieri for DYTurbo
//  
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

//   nlo includes
#include "nlo-color.h"
#include "ampq2g1p2.h"
#include "defmacros.h"
#include <math.h>
#include "Inter_params.h"


////////////////////////////////////////////////////////////////////////
// Interferencia LO QQBAR and QG channels
////////////////////////////////////////////////////////////////////////

 double ampq2g1p2::LO_MqqbarH(int p1, int p2, int p3,int p4,int p5) const
  { //Este elemento de matriz calcula la produccion de Higgs gg->H
	// Falta hacer el modulo cuadrado que se hace en el momento de multiplicar por las PDFs
	// In the infinit top mass limit
	// Por el momento solo la contribucion top
	// alpha_s se mete en la parte de las PDFs cuando se hace el cuadrado

double p1p2 =  S(1,2)/2.0;  //The dot product defined in Version B5 has no factor 2
double p1p4 = -S(1,5)/2.0;
double p1p5 = -S(1,3)/2.0;
double p2p4 = -S(2,5)/2.0;
double p2p5 = -S(2,3)/2.0;
double s45  =  S(4,5);



      double vmartin = _VEV_/sqrt(2.0); // Ok

      double Color=4.0;
      double vev= _VEV_;   //Reescribimos la amplitud para que aparezca 1/Vev**2 (usamos vHiggs=174 y vev=246)
      //  double alpha_qed= alpha_qed ;
      //double vertice=alpha_qed/3.0/__DEB_PI/(vev*vev)*3.0/2.0*sqrt(norm(CorcheteMartin(s45)));
      double vertice=-alpha_qed/3.0/__DEB_PI/(vev*vev)*3.0/2.0*(16.0/9.0-8.3442);

      double promedioqq=1.0/(2.0*2.0*3.0*3.0);
      double MW = mW ;
      double SW = sqrt(sin2W);
      double el = sqrt(alpha_qed*4.0*__DEB_PI);
      double MH = _MH_;
      double Gs = 1.0; // to be icluded at the end alpha_s
      double GH = _GammaH_;
      double EQ = 1.0; // to be included at the end UU or DD



      double amp= -promedioqq*(2.0*MW*SW)*
        (Color*el*el*EQ*EQ*Gs*Gs*Gs*Gs*(p1p5*p1p5 + p2p5*p2p5)*
        (-MH*MH + 2.0*(p1p2 + p1p5 + p2p5))*(p1p4*(p1p2 + p2p5)*
        (p1p2*(p1p2 + p1p5 + p1p4) + (p1p2 + p1p4)*p2p5) +
        (p1p2*(p1p2 + p1p5)*(p1p2 + p1p5) + (p1p2*(p1p2 + p1p5) - 2.0*p1p5*p1p4)*
        p2p5)*p2p4 + (p1p2 + p1p5)*(p1p2 + p1p5)*p2p4*p2p4)*vertice)/
        (12.0*MW*p1p2*p1p4*(p1p2 + p1p5 + p1p4)*(GH*GH*MH*MH + pow((-MH*MH
        + 2.0*(p1p2 + p1p5 + p2p5)),2.0))*p2p4*(p1p2+p2p5+p2p4)*__DEB_PI*__DEB_PI*SW);



    return amp;
 }

double ampq2g1p2::LO_MqgH(int p1, int p2, int p3,int p4,int p5) const
  { //Este elemento de matriz calcula la produccion de Higgs gg->H


    double amp2;

double p1p2 =  S(1,2)/2.0;  //The dot product defined in Version B5 has no factor 2
double p1p4 = -S(1,5)/2.0;
double p1p5 = -S(1,3)/2.0;
double p2p4 = -S(2,5)/2.0;
double p2p5 = -S(2,3)/2.0;
double s45 = abs(S(4,5));

            /*
             cout << "From amp S(1,5)/2.0  :  " <<  S(2,3)/2.0 << endl;

             cout << "***************************************************************************** "<< endl;
             cout << "***************************************************************************** "<< endl;
             cout << "***************************************************************************** "<< endl;
            */




      double vmartin = _VEV_/sqrt(2.0); // Ok

    double promedioqg=1.0/(2.0*2.0*3.0*8.0);
     double Color=4.0;
      double vev= _VEV_;     // !Reescribimos la amplitud para que aparezca 1/Vev**2 (usamos vHiggs=174 y vev=246)
        //double vertice=alpha_qed/3.0/__DEB_PI/(vev*vev)*3.0/2.0*sqrt(norm(CorcheteMartin(s45)));
        double vertice=-alpha_qed/3.0/__DEB_PI/(vev*vev)*3.0/2.0*(16.0/9.0-8.3442);
       // Checkear a que configuracion de vértice pertenece estos números

       double MW = mW ;
      double SW = sqrt(sin2W);
      double el = sqrt(alpha_qed*4.0*__DEB_PI);
      double MH = _MH_;
      double Gs = 1.0; // to be icluded at the end alpha_s
      double GH = _GammaH_;
      double EQ = 1.0; // to be included at the end UU or DD
      /*
p1p2 =  10000000.0/2.0;
p1p4 = -70000000.0/2.0;
p1p5 = -80000000.0/2.0;
p2p4 = -30000000.0/2.0;
p2p5 = -11000000.0/2.0;
      */

      amp2=promedioqg*vertice*(2.0*MW*SW)*
        (Color*el*el*EQ*EQ*Gs*Gs*Gs*Gs*(p1p2*p1p2 + p2p5*p2p5)*
        (-MH*MH + 2.0*(p1p2 + p1p5 + p2p5))*
        (p1p2*p1p2*p1p2*p1p2 + p1p4*(2.0*p1p5*p1p5*(p1p5 + p1p4) + p1p5*(3.0*p1p5 +
        2.0*p1p4)*p2p5 + (p1p5 + p1p4)*p2p5*p2p5) + p1p5*p1p5*(p1p5 +
        2.0*p1p4 + p2p5)*p2p4 + p1p5*p1p5*p2p4*p2p4 + p1p2*p1p2*p1p2*(3.0*p1p5 +
        2.0*(p1p4 + p2p5 + p2p4)) + p1p2*p1p2*(3.0*p1p5*p1p5 + p1p4*p1p4 +
        (p2p5 + p2p4)*(p2p5 + p2p4) + 2.0*p1p4*(2.0*p2p5 + p2p4) + p1p5*(5.0*p1p4 +
        4.0*p2p5 + 5.0*p2p4)) + p1p2*(p1p5*p1p5*p1p5 + 2.0*p1p4*p2p5*(p1p4 + p2p5 +
        p2p4) + p1p5*p1p5*(5.0*p1p4 + 2.0*p2p5 + 4.0*p2p4) + p1p5*(2.0*p1p4*p1p4 +
        (p2p5 + p2p4)*(p2p5 + 2.0*p2p4) + p1p4*(6.0*p2p5 + 4.0*p2p4)))))/
        (12.0*MW*p1p5*p1p4*(p1p2 + p1p5 + p1p4)*(GH*GH*MH*MH +
        (-MH*MH + 2.0*(p1p2 + p1p5 + p2p5))*(-MH*MH + 2.0*(p1p2 + p1p5 + p2p5)))*(p1p2 + p1p4 + p2p4)*
         (p1p2 + p1p5 + p1p4 + p2p5 + p2p4)*__DEB_PI*__DEB_PI*SW);

     /*

       double otra = promedioqg*vertice*(2.*MW*SW)*
        (Color*el*el*EQ*EQ*Gs*Gs*Gs*Gs*(p1p2*p1p2 + p2p5*p2p5)*
        (-MH*MH + 2*(p1p2 + p1p5 + p2p5))*
        (p1p2*p1p2*p1p2*p1p2 + p1p4*(2*p1p5*p1p5*(p1p5 + p1p4) + p1p5*(3*p1p5 +
        2*p1p4)*p2p5 + (p1p5 + p1p4)*p2p5*p2p5) + p1p5*p1p5*(p1p5 +
        2*p1p4 + p2p5)*p2p4 + p1p5*p1p5*p2p4*p2p4 + p1p2*p1p2*p1p2*(3*p1p5 +
        2*(p1p4 + p2p5 + p2p4)) + p1p2*p1p2*(3*p1p5*p1p5 + p1p4*p1p4 +
        (p2p5 + p2p4)*(p2p5 + p2p4) + 2*p1p4*(2*p2p5 + p2p4) + p1p5*(5*p1p4 +
        4*p2p5 + 5*p2p4)) + p1p2*(p1p5*p1p5*p1p5 + 2*p1p4*p2p5*(p1p4 + p2p5 +
        p2p4) + p1p5*p1p5*(5*p1p4 + 2*p2p5 + 4*p2p4) + p1p5*(2*p1p4*p1p4 +
        (p2p5 + p2p4)*(p2p5 + 2*p2p4) + p1p4*(6*p2p5 + 4*p2p4)))))/
        (12.*MW*p1p5*p1p4*(p1p2 + p1p5 + p1p4)*(GH*GH*MH*MH +
        pow((-MH*MH + 2*(p1p2 + p1p5 + p2p5)),2))*(p1p2 + p1p4 + p2p4)*
         (p1p2 + p1p5 + p1p4 + p2p5 + p2p4)*__DEB_PI*__DEB_PI*SW);

   cout << amp2/otra << endl;
*/
    //   amp2 = 	-S(2,3)/2.0/p1p2*SW/MW; //p1p2/p1p4*
     //    cout << " S(2,4)  :  " << -S(2,4)/2.0 << endl;
      // cout << " S(1,5)  :  " << S(1,5) << endl;
      // cout << " S(2,4)  :  " << S(2,4) << endl;
      // cout << " S(2,5)  :  " << S(2,5) << endl;
       //cout << SW/MW << endl;
    return amp2;

 }
 
 
 
