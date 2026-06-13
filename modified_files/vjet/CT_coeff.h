#ifndef CT_COEFF
#define CT_COEFF 

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "nlo-process.h"
#include "Integral_class.h"
#include "mstwpdf.h"
#include "qt_splitt.h"


//DECLARACION DE CONSTANTES COMUNES
#define __DEB_PI    3.14159265358979323846
#define __DEB_2PI   6.28318530717958647692
#define __DEB_4PI2  157.91367041742973790108
#define __DEB_1_8PI 0.03978873577297383395
#define Z3   1.20205690315959429
#define Z2   1.6449340668482262
#define _CFH  1.33333333333333333333
// CARGAS DE LOS Q'S
#define Qu  0.66666666666666666666 
#define Qd -0.33333333333333333333
#define Qu4 0.19753086419753086419
#define Qd4 0.01234567901234567901
#define __alfaQED2 (1.0/18769.0)
#define Qsum4 (2*Qu4+3*Qd4)

// CONSTANTES A UTILIZAR EN LOS CONTRATERMINOS
#define A1q (4.0/3.0)
#define B1q -2.0

// COEFICIENTES DE RESUMACION
#define beta0(Nf) ((33.0-2.0*Nf)/12.0)
#define beta1(Nf) ((153.0-19.0*Nf)/24.0)
// Kappa=67/6d0-(pi**2)/2d0-5d0/9d0*nf
#define Kappa(Nf) (67.0/6.0-(__DEB_PI*__DEB_PI)/2.0-5.0/9.0*Nf)
//#define Kappa(Nf) 1.0
     
#define A2q(Nf) (0.5*A1q*Kappa(Nf))

   
// B2q en el esquema viejo
//#define B2q(Nf) (4.0/9*(__DEB_PI*__DEB_PI-3.0/4-12*Z3)+(11.0/9*__DEB_PI*__DEB_PI-193.0/12+6*Z3)+Nf/6.0*(17.0/3-4.0/9*__DEB_PI*__DEB_PI))
//#define B2q(Nf) -3.39493 + 2.0*C1qqdelta*(2.75 - 0.166667*Nf) + 0.421096*Nf + 0.666667*Z3
// Bq2 en el Hard Scheme 04/11/2015
//#define B2q(Nf) (1.0/48.0)*_CFH*(2.0*Nf*(1.0 + 4.0*Z2) +    9.0* _CFH*(-1.0 + 8.0* Z2 - 16.0*Z3) +    3.0*(-17.0 - 44.0*Z2 + 72.0*Z3))

//     Delta term in c1qq coefficient
//     C1qqdelta=(pi**2-8)/3d0   ! Drell-Yan
//#define C1qqdelta(ah1) ah1/2.0/2.0 
//      C1qqdelta= ah1/2d0  ! REVISAR  2gamma !!!!!!!!!!!!!!!!!
//     Delta term in P2qq splitting function (as/pi normalization)
#define Delta2qq(Nf) (16.0/9*(3.0/8.0-__DEB_PI*__DEB_PI/2+6*Z3)+4*(17.0/24.0+11.0*__DEB_PI*__DEB_PI/18-3*Z3)-2.0/3.0*Nf*(1.0/6.0+2*__DEB_PI*__DEB_PI/9.0))/4.0


//    Coefficients of D0 and D1 in P*P (as/pi normalization)
#define D0qqqq (8.0/3)
#define D1qqqq (32.0/9)

//    Coefficients of delta(1-z) in P*P
#define Deltaqqqq (4.0/9*(9.0/4-2*__DEB_PI*__DEB_PI/3.0))

//     H2qq contribution: coefficient of delta(1-z)
#define H2qqdelta(Nf) (-2561.0/144+127.0*Nf/72+3*__DEB_PI*__DEB_PI/2-19.0*Nf*__DEB_PI*__DEB_PI/81+49.0*pow(__DEB_PI,4.0)/324 +58.0*Z3/9 + 8.0*Nf*Z3/27)

//     H2qq contribution: coefficient of D0(z)
#define H2qqD0(Nf) (-404.0/27+(56.0*Nf)/81+14*Z3)


// Declaracion de Msqc
double S[11][11]={ {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0},
                   {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0},
                   {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0},
                   {0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},
                   {0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0},
                   {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
                   {0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},
                   {0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
                   {0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
                   {0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
                   {1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
			   };

double eq2[6] = {0.0,4.0/9.0,1.0/9.0,1.0/9.0,4.0/9.0,1.0/9.0};


double Msqc(int jj, int kk){
	
return eq2[abs(jj)]*eq2[abs(jj)]*S[jj+5][kk+5];
	
	
}
     

//DECLARACION DE FUNCIONES COMUNES

double Itilde(int ,double );



extern "C" {
    double ctnnlo_(double * q2,double * xx10,double * xx20,double * x1pc
     ,double * x2pc,double * xmurc,double * xmufc,double * xmioc,int * isetc,double * a2c,int * orderc,
     double * LL1c,double * LL2c,double * LL3c,double * LL4c,double *aah1);    }
 

 inline double CTNNLO(double q2,double xx10,double xx20,double x1pc
     ,double x2pc,double xmurc,double xmufc,double xmioc,int isetc,double a2c,int orderc,
      double LL1c,double LL2c,double LL3c,double LL4c, double aah1)  {
 return  ctnnlo_(&q2,&xx10,&xx20,&x1pc
     ,&x2pc,&xmurc,&xmufc,&xmioc,&isetc,&a2c,&orderc,&LL1c,&LL2c,&LL3c,&LL4c,&aah1);
 }

extern "C" {
    double h2_(double * q2,double * xx10,double * xx20,double * z1c
     ,double * z2c,double * xmurc,double * xmufc,int * isetc,double * a2c,int * orderc,
     double *aah1);    }
 

 inline double H2(double q2,double xx10,double xx20,double z1c
     ,double z2c,double xmurc,double xmufc,int isetc,double a2c,int orderc,
       double aah1)  {
 return  h2_(&q2,&xx10,&xx20,&z1c
     ,&z2c,&xmurc,&xmufc,&isetc,&a2c,&orderc,&aah1);
 }

//C2 contributions

extern "C" {
    double c2qg_(double * x);    }
 

 inline double C2QG(double x)  {
 return  c2qg_(&x);
 }
 
 extern "C" {
    double c2qqp_(double * x);    }
 

 inline double C2QQP(double x)  {
 return  c2qqp_(&x);
 }

extern "C" {
    double c2qqb_(double * x);    }
 

 inline double C2QQB(double x)  {
 return  c2qqb_(&x);
 }
 
 extern "C" {
    double c2qqreg_(double * x);    }
 

 inline double C2QQREG(double x)  {
 return  c2qqreg_(&x);
 }

// CT gg NLO

#define A1g 3.0
#define B1g -2.0*beta0(5)

/*
C     B2g from qt code (checked !)
c OLD SCHEME LEAN CHECKED!!!!!!!!!!!!!!!!!!!!!!!!!
*/
/*
c      B2g=9*(23/24d0+(11*pi**2)/18d0-3*Z3/2d0)+
c     /    2*nf/3d0-3*nf*(1/12d0+pi**2/9d0)-11/2d0
     
C     B2g from qt code (checked !)
c HARD SCHEME LEAN CHECKED!!!!!!!!!!!!!!!!!!!!!!!!
*/
// #define  B2g=1.d0/24.d0*(33.d0*pi**2 - 2.d0*nf*(-16.d0 + pi**2) -36.d0*(8.d0 + 9.d0*Z3))


#endif


