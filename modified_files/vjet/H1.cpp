//      H1.cpp
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

#include "H1.h"


#define __DEB_PI    3.14159265358979323846
#define __DEB_2PI   6.28318530717958647692
#define __DEB_4PI2  157.91367041742973790108
#define __DEB_1_8PI 0.03978873577297383395


// CARGAS DE LOS Q'S
#define Qu  0.66666666666666666666 
#define Qd -0.33333333333333333333
#define Qu4 0.19753086419753086419
#define Qd4 0.01234567901234567901

  

using namespace std;








weight_hhc2ph GeT_pdf(double, double, double, unsigned int, unsigned int);

double  H1( event_hhc2ph p){
	
	
    double t =-2.0*(p[-1]*p[-3]), u =-2.0*(p[0]*p[-3]),
                s=2.0*(p[-1]*p[0]);
         

   double  hqq =  (Cf*((-21.0 + 2.0*__DEB_PI*__DEB_PI)*(t*t + u*u) 
            + 3.0*u*(2.0*t + 3.0*u)*log(-(t/s)) + 
            3.0*(2.0*t*t + 2.0*t*u + u*u)*log(-(t/s))*log(-(t/s)) + 
            3.0*t*(3.0*t + 2.0*u)*log(-(u/s)) + 
          3.0*(t*t+2.0*t*u+2.0*u*u)*log(-(u/s))*log(-(u/s))))/(3.0*(t*t+u*u)) ;
     
  
    return (hqq+Cf*__DEB_PI*__DEB_PI/3.0 );
  
}


weight_hhc2ph  procH1(event_hhc2ph p,double muf2,unsigned int nfl, double z1, double z2){
	double Q2,x1p,x2p;
	
	Q2 = 2.0*((p[-2])*(p[-3]));
	
	
	////////////// CALCULAR x10 y x20  ////////////////////////
	double s = p[hadron(-1)]*p[hadron(0)];
    double x10 = p[-1]*p[hadron(0)]/s;
    double x20 = p[0]*p[hadron(-1)]/s;
	
	
	
	
	x1p = x10/z1;
	x2p = x20/z2;
	double LF = log(Q2/muf2);
	
	
	
	
	// Prepara las PDF's para cada una de las patas
	weight_hhc2ph _M_pdfL1,_M_pdfL2,_M_pdf00,M_CT;
	
	weight_hhc2ph retval,retvalqqL1,retvalqqL2,retvalqqL11,retvalqqL22,retvalL1,retvalL2,retval12,retval1_4;
	
	
	//Calcula el BORN	
	 double ampB=ampBORN(p);
	 
	 	
	
	
	
	
	
	_M_pdf00=GeT_pdf(x10,x20,muf2, 2, 3);
	_M_pdfL1=GeT_pdf(x1p,x20,muf2, 2, 3);
	_M_pdfL2=GeT_pdf(x10,x2p,muf2, 2, 3);
	

	//CANAL QQB
	weight_hhc2ph  Waux00,WauxL1qq,WauxL2qq,WauxL1qq0,WauxL2qq0;
	weight_hhc2ph  WauxC1QQL1,WauxC1QQL2;
	
	// Témino de la Delta ---- H1st SOLO EL H1^F
	
	retval12[12] = Qu4*ampB*H1(p);
	retval12[14] = Qd4*ampB*H1(p);
	
	Waux00 = retval12;
	
	weight_hhc2ph  Wauxdelta;
	
	Wauxdelta = Waux00*_M_pdf00;
	
	// Témino de la Delta ---- H1st SOLO EL H1^F
	// GUARDA PORQUE ESTE TÉRMINO TIENE LOS CUTS DE LA CINEMÁTICA
	// DEL BORN ASIQUE NO LO PUEDO PONER EN LA CINEMÁTICA DEL
	// CT YA QUE TIENE LA EMISION
	
	//-----------------C1qq-------------------------

//H1st: non delta term PATA 1 //Procesos U1 Ub2
  
    double tH1C1L1;
    tH1C1L1 = Cqq(z1)*(-log(x10));

	retval12[12] =tH1C1L1*Qu4*ampB;
	//Procesos D1 Db2 
	retval12[14] =tH1C1L1*Qd4*ampB;
	
	WauxC1QQL1 = retval12;
	weight_hhc2ph  WL1Cqq;
	 WL1Cqq = WauxC1QQL1*_M_pdfL1;
//H1st: non delta term PATA 2 //Procesos U1 Ub2
  
   double tH1C1L2;
   tH1C1L2 = Cqq(z2)*(-log(x20));

	retval12[12] = tH1C1L2*Qu4*ampB;
	//Procesos D1 Db2 
	retval12[14] = tH1C1L2*Qd4*ampB;
	
	WauxC1QQL2 = retval12;
	   
	 weight_hhc2ph  WL2Cqq, WC1qq;
	 WL2Cqq = WauxC1QQL2*_M_pdfL2;
	 WC1qq =(2.0*(WL1Cqq +  WL2Cqq));
	//------------------------FIN C1qq------------------------------
	
	
	
	
	
	
	
	
	
	
	//-----Term Gammaqq coefficient--pata 1 (Esto significa que la 
	// pata 1 emite parton) La contribucion con la PDF con xi0 se pone aparte
	
	double tH1stqL1;
    tH1stqL1 = Pqq(z1)*(-log(x10));

	retvalqqL1[12] =tH1stqL1*Qu4*ampB;
	//Procesos D1 Db2 
	retvalqqL1[14] =tH1stqL1*Qd4*ampB;
	
	WauxL1qq = retvalqqL1;
	weight_hhc2ph  WL1qq;
	WL1qq = WauxL1qq*_M_pdfL1;
	
	//-----Term Gammaqq coefficient--pata 2 (Esto significa que la 
	// pata 1 emite parton) La contribucion con la PDF con xi0 se pone aparte
	
	double tH1stqL2;
    tH1stqL2 = Pqq(z2)*(-log(x20));

	retvalqqL2[12] =tH1stqL2*Qu4*ampB;
	//Procesos D1 Db2 
	retvalqqL2[14] =tH1stqL2*Qd4*ampB;
	
	WauxL2qq = retvalqqL2;
	weight_hhc2ph  WL2qq;
	WL2qq = WauxL2qq*_M_pdfL2;
	
	//PATA 1 para PDF's eveluadas en xi0
	
	
   double tH1stqL11 = -Pqq(z1)*(-log(x10))*z1-Pqqint(x10);

	retvalqqL11[12] =tH1stqL11*Qu4*ampB;
	//Procesos D1 Db2 
	retvalqqL11[14] =tH1stqL11*Qd4*ampB;
	
	WauxL1qq0 = retvalqqL11;
	weight_hhc2ph  WL1qq0;
	WL1qq0 = WauxL1qq0*_M_pdf00;
	
	//PATA 2 para PDF's eveluadas en xi0
	
	
   double tH1stqL22 = -Pqq(z2)*(-log(x20))*z2-Pqqint(x20);

	retvalqqL22[12] =tH1stqL22*Qu4*ampB;
	//Procesos D1 Db2 
	retvalqqL22[14] =tH1stqL22*Qd4*ampB;
	
	WauxL2qq0 = retvalqqL22;
	weight_hhc2ph  WL2qq0;
	WL2qq0 = WauxL2qq0*_M_pdf00;
	
	// Comenzamos a sumar las contribuciones ---Menos la Delta!!!!
	weight_hhc2ph  Wqqaux,Wqq;
	
	Wqqaux =(2.0*(WL1qq + WL2qq + WL1qq0 + WL2qq0));
	
	//El factor 2 Extra se debe a que  alpha_s/2/Pi USO!!!!!!!!!
	//Y Massi usa alpha_s/Pi
	
	
	// Fin CANAL QQBAR a Primer orden
	
	
	
	
	
	
	
	
	
	
	
	
	//CANAL QG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//Vamos con la pata 1, es decir procesos del tipo G1 * (Q2)
	//Es decir van en los canales 3 y 4
	
	//----------------------C1qg--------------------------
	
   double tH1C1qgL1 = Cqg(z1)*(-log(x10));
	
	        
	retvalL1[3]=tH1C1qgL1*Qu4*ampB;
	retvalL1[4]=tH1C1qgL1*Qd4*ampB;
	
	        
	//Vamos con la pata 2, es decir procesos del tipo Q1 * (G2)
	//Es decir van en los canales 1 y 2
	
	
	
    double tH1C1qgL2 = Cqg(z2)*(-log(x20));
	
	        		        
	retvalL2[1]=tH1C1qgL2*Qu4*ampB;
	retvalL2[2]=tH1C1qgL2*Qd4*ampB;
	
	
	weight_hhc2ph  WauxC1qgL1,WauxC1qgL2,WC1qgL1,WC1qgL2;
	
	WauxC1qgL1=retvalL1;
	WauxC1qgL2=retvalL2;
	
	WC1qgL1 =((2.0)*(_M_pdfL1*WauxC1qgL1));
	WC1qgL2 =((2.0)*(_M_pdfL2*WauxC1qgL2));
	
	weight_hhc2ph  WC1qg= WC1qgL1 + WC1qgL2 ;
	
	//-------------------   FIN    C1qg  --------------------------
	
	
	//-----Term Gammaqg coefficient--pata 1
	double tH1stL1;
	tH1stL1 = Pqg(z1)*(-log(x10));
	        
	retvalL1[3]=tH1stL1*Qu4*ampB;
	retvalL1[4]=tH1stL1*Qd4*ampB;
	
	//-----Term Gammaqg coefficient--pata 2
	double tH1stL2;
	tH1stL2 = Pqg(z2)*(-log(x20));
	        
        
	retvalL2[1]=tH1stL2*Qu4*ampB;
	retvalL2[2]=tH1stL2*Qd4*ampB;
	
	//---Pata 1 por PDF pata 1
	
	weight_hhc2ph  WauxL1,WauxL2,WTL1,WTL2;
	
	WauxL1=retvalL1;
	WauxL2=retvalL2;
	
	//El factor 2 Extra se debe a que las pdfs vienen con alpha_s/2/Pi
	//Y Massi usa alpha_s/Pi
	WTL1 =((2.0)*(_M_pdfL1*WauxL1));
	WTL2 =((2.0)*(_M_pdfL2*WauxL2));
	
	//Wauxdelta contiene el H1qq o ah1
	weight_hhc2ph  tH1stMIO = WC1qq+WC1qg;
	weight_hhc2ph  tH1stF =  Wqqaux+WTL1+WTL2;
	
	return tH1stMIO + LF*tH1stF;
	
	
}





