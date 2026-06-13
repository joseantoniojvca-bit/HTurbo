//      qt_splitt.cc
//      
//      Copyright 2012 Leandro <leandro@leandro-K53SC>
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
//      
//      


#include <iostream>
#include "qt_splitt.h"
#define __DEB_PI    3.14159265358979323846


extern "C" {
    double myli2_(double * x);    }
 

 inline double MYLI2(double x)  {
 return  myli2_(&x);
 }


//Cqg coefficient (con factor asopi normalization)
double Cqg(double z){
	
	return 0.5*z*(1.0-z);
	
}

//Cqq coefficient (con factor asopi normalization)
double Cqq(double z){
	
	return (2.0/3.0)*(1.0-z);
	
}
//Splitting Pqg (con asopi normalization)
double Pqg(double z){
	
	return 0.25*(1.0 - 2.0*z*(1.0-z));
}

//Splitting Pqq (con asopi normalization)
double Pqq(double z){
	
	return (2.0/3.0)*(1.0 + z*z)/(1.0-z);
}

//Integral de  Pqq de 0 a z
double Pqqint(double z){
	
	return -(2.0/3.0)*(z + z*z/2.0+ 2.0*log(1.0-z));
}

// FUNCIONES PROPIAS SEGUNDO ORDEN -> NNLO

// Integral of 1/(1-x) from 0 to z

      double D0int(double z){
     
      return -log(1-z);
  }

// Integral of log(1-x)/(1-x) from 0 to z

      double D1int(double z){
      
      return -0.5*log(1-z)*log(1-z);
      }

//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//
//                P*P convolutions
//
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

// Regular part of Pqq*Pqq (checked !)

      double Pqqqq(double z){
     

      return 4.0/9*(-4*log(z)/(1-z)-2*(1-z)
     +(1+z)*(3*log(z)-4*log(1-z)-3));
      }


// Pqq*Pqg (checked !)

      double Pqqqg(double z) {
     
   
      return 1.0/3*((z*z+(1-z)*(1-z))*log((1-z)/z)
       -(z-0.5)*log(z)+z-0.25);
      }

// Pqg*Pgq (checked !)

      double Pqggq(double z){
     
    
      return 1.0/3*(2.0/3/z+(1+z)*log(z)-2.0/3*z*z-0.5*(z-1));
      }


// Full Pqg*Pgg (checked !)

      double Pqggg(double z, int nf){
      double beta0=(33-2*nf)/12.0;
     

      return 1.5*(1/3.0/z+(z*z-z+0.5)*log(1-z)
        +(2*z+0.5)*log(z)+0.25+2*z-31.0/12*z*z)+beta0*Pqg(z);

      }


//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//
//                C*P convolutions
//
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

// Cqq*Pqq (without delta term in Cqq) (checked !)

      double CqqPqq(double z){
      
      return 2.0/9*(1-z)*(4*log(1-z)-2*log(z)-1);
      }

// Cqq*Pqg (without delta term in Cqq) (checked !)

      double CqqPqg(double z){
      
     
      return (-2+z+z*z-(1+2*z)*log(z))/6.0;
      }

// Cqg*Pgq (checked !)

      double CqgPgq(double z) {
     
      return (1.0/3/z-1+2*z*z/3-z*log(z))/3.0;
      }

// Cqg*Pgg (checked !)

      double CqgPgg(double z, int nf){
      double beta0=(33-2*nf)/12.0;
  
      return 3.0/4*(2*z*(1-z)*log(1-z)-4*z*log(z)
        +1.0/3/z-1-5*z+17.0*z*z/3)+beta0/2*z*(1-z);
      }

//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//
//          Two loop AP:  pqq of ESW is my 3/2 Pqq
//
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

//    Pqq NS: Eq. (4.107) ESW (no 1/(1-x)_+ and delta term)

      double P2qqV(double x, int nf){
     
      double P2qqVv=16.0/9*(-(2*log(x)*log(1-x)+1.5*log(x))*3.0/2*Pqq(x)
          -(1.5+3.5*x)*log(x)-0.5*(1+x)*log(x)*log(x)-5*(1-x))
          +4*((0.5*log(x)*log(x)+11.0/6*log(x))*3.0/2*Pqq(x)
          -(67.0/18-__DEB_PI*__DEB_PI/6)*(1+x)
          +(1+x)*log(x)+20.0/3*(1-x))
          +2.0/3.0*nf*(-log(x)*Pqq(x)+10.0/9*(1+x)-4.0/3*(1-x));

//     Change to as/pi normalization

            return P2qqVv/4.0; 
      }


//    Pqqb NS: Eq. (4.108) ESW

      double P2qqbV(double x){
     
     double  P2qqbVv=-2.0/9*(3.0*Pqq(-x)*S2(x)+2*(1+x)*log(x)+4*(1-x));
      
//    Change to as/pi normalization

     
      return P2qqbVv/4.0 ;
      }



//    Pqg Singlet: Eq. (4.110) ESW (ESW Pqg is 4 times my Pqg)

      double P2qg(double x){
      
        
      double logx=log(x);
      double logomxsx=log((1-x)/x);

      double   P2qgv=(2.0/3*(4-9*x-(1-4*x)*logx-(1-2*x)*logx*logx+4*log(1-x)
         +(2*logomxsx*logomxsx-4*logomxsx-2.0/3*__DEB_PI*__DEB_PI+10.0)*4*Pqg(x))
         +1.5*(182.0/9+14.0/9*x+40.0/9/x+(136.0/3*x-38.0/3)*logx
         -4*log(1-x)-(2+8*x)*logx*logx+8*Pqg(-x)*S2(x)
         +(-logx*logx+44.0/3*logx-2*log(1-x)*log(1-x)+4*log(1-x)+__DEB_PI*__DEB_PI/3
         -218.0/9)*4*Pqg(x)))/4.0/2.0;

//     Change to as/pi normalization // Los ultimos factores 4.0 y 2.0
//     de la ecuacion de arriba se refieren a:

//      P2qg=P2qg/4d0
  
//     Divide by 2 to eliminate 2nf factor

//      P2qg=P2qg/2d0

      return P2qgv;
      }

//    Pqq Pure Singlet appearing in ESW Eq. (4.95)
//     PSqq=PSqqb
//    Obtained through Eq.(4.101)
//     PSqq=1/2/nf (P2qq-P2qqbV-P2qqV) (contains only CF TR=2/3)

      double P2qqS(double x){
     
      double P2qqSv=2.0/3*(20 - 18*x + 54*x*x - 56*x*x*x
         +3*x*(3 + 15*x + 8*x*x)*log(x) 
         - 9*x*(1 + x)*log(x)*log(x))/(9*x);
      
     

      return P2qqSv/4.0;
      }


//    S2: Eq. (4.114) ESW

      double S2(double x){
      
     
      return -2*MYLI2(-x)+0.5*log(x)*log(x)-2*log(x)*log(1+x)-__DEB_PI*__DEB_PI/6;
      }




/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// CT gg NLO
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

// gg splitting function (not normalized !)

      double Pgg(double x){
      
      double Pggf=1.0/(1.0-x)+1.0/x-2.0+x*(1.0-x);

      return Pggf;
      }

/*
CC gg splitting function: regular part (with asopi normalization)
c NO COINCIDE CON QUELLA MIA PERO ALLA FINE DA LO STESSO
C FOGLIO 7 MIO
C CONTROLLARE ALL'ORDINE SUCCESSIVO SE TORNA ANCHE
*/

      double Pggreg(double z){
      
  double  Pggregf=3.0*((1.0-2.0*z)/z+z*(1.0-z));
  return Pggregf;
  }
// gq splitting function (with asopi normalization)

     double  Pgq( double z){
      
      double Pgqf=2.0/3.0*(1.0+(1.0-z)*(1.0-z))/z;
      return Pgqf;
  
}



// Cgq coefficient (with asopi normalization)

      double Cgq(double z){
      
      double Cgqf=2.0/3.0*z;
      return Cgqf;
  }


/*
// Integral of 1/(1-x) from 0 to z

     double D0int(double z){
      
      double D0intf=-log(1.0-z);
      return D0intf;
     
}
// Integral of log(1-x)/(1-x) from 0 to z

      double D1int(double z){
      
      double D1intf=-0.5*log(1.0-z)*log(1.0-z);
      return D1intf;
     
}

*/































