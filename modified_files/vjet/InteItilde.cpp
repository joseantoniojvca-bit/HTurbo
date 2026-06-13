//      Intetilde.cpp
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



int IER = 0;
const int MAXINT = 500;
double FV[MAXINT],FU[MAXINT],FWe[MAXINT]
    ,U[MAXINT],V[MAXINT],RESULT[MAXINT],ERR[MAXINT];
wadp  ADPINT(double  ,double, double , int  );
wadp ADPCAL(int ,double , double , double , double ,double ,int);
void coutwadp(const wadp& );
double Itilde(int ,double );


double F(double variable,double z, int n){
	  double t = variable;;
	  int nu,nn;
      double fbb;
      double zz=z;
      nu=1;
      nn=n;
      if(nn == 0){
         fbb =exp(-zz*cosh(t))*cosh(nu*t);}
      else if(nn == 1){
         fbb=exp(-zz*cosh(t))*t*sinh(nu*t);}
      else if(nn == 2){
         fbb=exp(-zz*cosh(t))*t*t*cosh(nu*t);}
      else if(nn == 3){
         fbb=exp(-zz*cosh(t))*t*t*t*sinh(nu*t);}
            
      return fbb;
	 
	  
 }

wadp  ADPINT(double  A,double  B,double z, int n){
      
     
//c.....Integral of F(X) from A to B, with error
//c.....less than ABS(AERR) + ABS(RERR*INTEGRAL)
//c.....Best estimate of error returned in ERREST.
//c.....Error code is IER: zero if OK, non-zero if in trouble.
      
     
      std::vector< std::vector<double> > Cadpcal(MAXINT,std::vector<double>(2));
//     Work space:
      
     
      int NUMINT = 499;
      double DX = (B-A)/ NUMINT;
      for(int i = 1; i < NUMINT+1 ; i++){
         if (i == 1) {            U[i] = A;
            FU[i] = F(U[i],z,n);}
         else {
            U[i] = V[i-1];
            FU[i] = FV[i-1];
         }
         if (i == NUMINT) {
            V[i] = B ;}
         else {
            V[i] = A + DX * i;
         }
         FV[i] = F(V[i],z,n);
      //   cout << FV[i] << endl;
         
         wadp  tmp;
         tmp = ADPCAL(i, U[i], V[i],  FU[i], FV[i],z,n);
       //   coutwadp(tmp); 
         Cadpcal[i][0]= tmp.WRE();
        // cout << Cadpcal[i][0] <<endl;
         Cadpcal[i][1]= tmp.WER();
         FWe[i] = tmp.WFW();  
         
   //  coutwadp(tmp);
   //  cout << Cadpcal[i][0] <<endl;
                              }
      

    
		 
          double  ADPINT1 = 0.0;
          double   ERREST1 = 0.0;
         for(int i= 1; i< NUMINT + 1; i++){
         ADPINT1 += Cadpcal[i][0];
         ERREST1 += Cadpcal[i][1];
                                       }
                                       
                                     //  cout << (ERREST1/ADPINT1) <<endl;
        if((ERREST1/ADPINT1) > 0.01){
	std::cout << "+++ El error en la integral de Bessel es mayor al 1%, ++++ "<<std::endl;
	std::cout << "+++ revise el archivo correspondiente e incremente MAXINT, " <<std::endl;
	std::cout << "+++++++++ o añada rutina de control de error +++++++++++++ " <<std::endl;
            }                          

             return wadp(ADPINT1,ERREST1,0.0);
}




wadp ADPCAL(int ii,double Uii, double Vii, double FUii, double FVii,double z,int n){

       
//c.....Fill in details of interval I given endpoints
      
       
     double FW1 = F( (Uii + Vii) /2.0,z,n);
     double DX = Vii - Uii;
     double RESULT1 = DX * (FUii + 4.0 * FW1 + FVii) / 6.0;
     double ERR1 = fabs(DX * (FUii - 2.0 * FW1 + FVii )/ 12.0);
     
      wadp Adpcal(RESULT1, ERR1,FW1);
  //     FW(I) = F( (U(I) + V(I)) /2.)
  //    DX = V(I) - U(I)
  //    RESULT(I) = DX * (FU(I) + 4. * FW(I) + FV(I)) / 6.
  //    ERR(I) = ABS(DX * (FU(I) - 2. * FW(I) + FV(I)) / 12.)
      
      //Hay que devolver todas estas cantidades
      
      return  Adpcal;
      
}

