//      Itilde.cpp
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

#include <cstdlib>
#include "wadp.h"
#include <vector>
#include <stdio.h>
#include <math.h>

#define Egamma 0.577215664902
#define z2 1.64493406685
#define z3 1.20205690316


//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//c.....Function BK(n,z)
//c.....BK(n,z) is the n-derivative of BesselK[nu,z]
//c.....with respect to nu in nu=1

//c.....Itilde defined as in the paper

wadp  ADPINT(double  ,double, double , int  );
using namespace std;
double zBK(int ,double );
double zbesselk0(double);
double zbesselk1(double);
double zbesselk2(double);
double zbesselk3(double);


double Itilde(int m,double xmio)
{      
      double argum,Itildee;
      double b0,logx;
      
     

      
      b0=2.0*exp(-Egamma);



      argum=b0*xmio;
      logx=log(xmio);
      
      if(m == 1){
         Itildee=-zBK(0,argum)/(xmio*xmio);}
      else if (m == 2){
         Itildee=2.0/(xmio*xmio)*(zBK(0,argum)*logx-zBK(1,argum));}
      else if (m == 3){
         Itildee=-3.0/(xmio*xmio)*(zBK(0,argum)*(logx*logx-z2)
               -2.0*zBK(1,argum)*logx+zBK(2,argum));}
      else if (m == 4){
         Itildee=4.0/(xmio*xmio)*(zBK(0,argum)*(pow(logx,3.0)-3.0*z2*logx+2.0*z3)
              -3.0*zBK(1,argum)*(logx*logx-z2)+3.0*zBK(2,argum)*logx
              -zBK(3,argum));
  }



      return Itildee;
     
 }     
      
      
      
      

//C     n-derivative of the function BesselK[nu,z] 
//C     with respect to nu for nu=1
//C     NOTE: IT IS MULTIPLIED by z

double zBK(int n,double z){
      
      double zbk;
     
      double max=10.0;

//     Use approximated form only for z<1.5d0

      if(n == 0){
      zbk=zbesselk0(z);}
      else if(n == 1){
      zbk=zbesselk1(z);}
      else if(n == 2){
      zbk=zbesselk2(z);}      
      else if(n == 3){ 
      zbk=zbesselk3(z);}
      
      
      // Otra integral!!!!!!!!!!!!!! Estoy implementando ADPINT
      if(z > 1.5){
		  wadp integral=ADPINT(0.0,max,z,n);
		  
      zbk=z*(integral.WRE());
            
      }

   return zbk;
      
}

//C     Approximated forms of BesselK[nu,z] for nu=1 and
//C     derivatives with respect to nu of BesselK[nu,z] at nu=1
//C     All functions multiplied by z

double zbesselk0(double z){
      
      double zm,loz,zbesselk00;
      
      zm=z/2.0;
      loz=log(zm);

      zbesselk00=1.0+z*zm*(loz-0.5*(1.0-2.0*Egamma))
             +pow(zm,4.0)*(loz-0.5*(2.5-2.0*Egamma))
              +pow(zm,6.0)/6.0*(loz-0.5*(10.0/3.0-2.0*Egamma))
              +pow(zm,8.0)/72.0*(loz-0.5*(47.0/12.0-2.0*Egamma))
              +pow(zm,10.0)/1440.0*(loz-0.5*(131.0/30.0-2.0*Egamma));
      
      return zbesselk00;
  }


 double   zbesselk1(double z){
      
      double zbesselk11,zm,loz;
      
      zm=z/2.0;
      loz=log(zm);
      zbesselk11=-(loz+Egamma)-zm*zm*(loz-1.0+Egamma)
              -0.25*pow(zm,4.0)*(loz-1.5+Egamma)
              -pow(zm,6.0)/36.0*(loz-11.0/6.0+Egamma)
              -pow(zm,8.0)/576.0*(loz-25.0/12.0+Egamma);
      

      return zbesselk11;
  }


double zbesselk2(double z){
      double a[14];
      double zbesselk22,loz,zm;
      a[0] = 1.15443132980306572;
       a[1] = 1.97811199065594511;
       a[2] = 0.154431329803065721;
       a[3] = 4.801792651508824500;
       a[4] = 0.806235643470665767;
       a[5] =-0.672784335098467139;
       a[6] = 3.285072828402112960;
       a[7] =-1.945338757678943440;
       a[8] =-0.181575166960855634;
       a[9] = 0.694195147571435559;
       a[10]=-0.607655744858515573;
       a[11]=-0.019182189839330562;
       a[12]= 0.068894530444636532;
       a[13]=-0.070514317816328185;


      zm=z/2;
      loz=log(zm);
      
      zbesselk22=loz*loz+a[0]*loz+a[1]
           +zm*zm*(2.0*pow(loz,3.0)/3.0+a[2]*loz*loz+a[3]*loz+a[4])
            +pow(zm,4.0)*(pow(loz,3.0)/3.0+a[5]*loz*loz+a[6]*loz+a[7])
            +pow(zm,6.0)*(pow(loz,3.0)/18.0+a[8]*loz*loz+a[9]*loz+a[10])
            +pow(zm,8.0)*(pow(loz,3.0)/216.0+a[11]*loz*loz+a[12]*loz+a[13]);
      return zbesselk22;
  }


double  zbesselk3(double z){
      double b[15];
      double zbesselk33,loz,zm;

       b[0] = 1.731646994704598580;
       b[1] = 5.934335971967835330;
       b[2] = 5.444874456485317730;
       b[3] =-1.268353005295401420 ;
       b[4] = 8.471041982558638170;
       b[5] =-3.026167526073320430;
       b[6] =-0.692088251323850355 ;
       b[7] = 2.809848746963509900;
       b[8] =-2.161466255000085060;
       b[9] =-0.104676472369316706;
       b[10]= 0.381989731242156681;
       b[11]=-0.367492827636283900;
       b[12]=-0.007844362856415627;
       b[13]= 0.027796539630842606;
       b[14]=-0.029917436634978395;


      zm=z/2;
      loz=log(zm);
      
      zbesselk33=pow(loz,3.0)+b[0]*loz*loz+b[1]*loz+b[2]
             +zm*zm*(pow(loz,3.0)+b[3]*loz*loz+b[4]*loz+b[5])
             +pow(zm,4.0)*(pow(loz,3.0)/4.0+b[6]*loz*loz+b[7]*loz+b[8])
             +pow(zm,6.0)*(pow(loz,3.0)/36.0+b[9]*loz*loz+b[10]*loz+b[11])
             +pow(zm,8.0)*(pow(loz,3.0)/576.0+b[12]*loz*loz+b[13]*loz+b[14]);
      zbesselk33 *=-1.0;

      return zbesselk33;
  }

      
