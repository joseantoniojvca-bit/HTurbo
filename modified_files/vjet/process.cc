//  Copyright (C) 2002 Zoltan Nagy
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
#include "hhc2ph-process.h"
#include "mstwpdf.h"
//extern c_mstwpdf *mstpdf;
// DYTurbo includes
#include "parton.h"
#include "settings.h"
#include "../dyres/dyres_interface.h"
#include "../dynnlo/dynnlo_interface.h"



#define __DEB_2PI   6.28318530717958647692





extern "C" 
{
extern struct {
 //La dimension de estas variables tienen que ser las definidas en
// el archivo fortran
          
    double fext[11];
   
}pdfsext_;
}

//extern "C" 
//{
//extern struct {
//	
//	int iset2g,ih12g,ih22g;
//	
//} parpdfs2gamma_;
//}

// extern "C" {
//     void getapdfs_(double *x,double *q,double *upv,double *dnv,double *usea,double *dsea,double *str,double *sbar,
//              double *chm,double *cbar,double *bot,double *bbar,double *glu,double *phot,int *order,int *ih);    }
 

 // inline void GeT_all_PDFS(double x,double q,double &upv,double &dnv,double &usea,double &dsea,double &str,double &sbar,
 //             double &chm,double &cbar,double &bot,double &bbar,double &glu,double &phot,int &order,int &ih)  {
 //  getapdfs_( &x, &q, &upv, &dnv, &usea, &dsea, &str, &sbar,
 //              &chm, &cbar, &bot, &bbar, &glu,&phot,&order,&ih);
 // }


 extern "C" {
    void mlmpdf_(int * isetproton1, int * ih1, double * xmuf2, double * x1,double *fmul, int * nf);    }
 

 inline void MLMPDF(int isetproton1, int ih1, double xmuf2, double  x1,double *fmul, int  nf)  {
  mlmpdf_(  &isetproton1,  &ih1,  &xmuf2,  &x1, fmul, &nf);
 }


    
  
  weight_hhc2ph  GeT_pdf(double x1, double x2, double mf2, unsigned int nu, unsigned int nd)
  {
	   unsigned int nfl=nu+nd;
	    
	     int iset,ih1,ih2;
	 iset=parpdfs2gamma_.iset2g;
	 ih1= parpdfs2gamma_.ih12g;
	 ih2= parpdfs2gamma_.ih22g;
	  
	 
//	   static double __f100[13], __f200[13];
//	   double *f100 = __f100+6, *f200 = __f200+6;
	
	  

	 // x1 = 0.1;
         // x2 = 0.5;
         // mf2 = 125.*125.;


	  
      //PDFs
      double muf = sqrt(mf2);
      double fx1[2*MAXNF+1],fx2[2*MAXNF+1];


      fdist_(opts.ih1,x1,muf,fx1);
      fdist_(opts.ih2,x2,muf,fx2);
	  
 	  
	  double fmul[13];
	  
    int ia, iq;
    weight_hhc2ph retval;
    
    static  double f1[13],f2[13];
    
    static const int iu[3] = {2,4,6}, id[3] = {1,3,5};
    static const int iub[3] = {8,10,12}, idb[3] = {7,9,11}; 
    
    //----- calculate the pdfs -----
     double upv,dnv,usea,dsea,str,sbar,chm,cbar,bot,bbar,glu,phot;
    //----- PROTON 1 
    
    	  // Update all PDF flavours.
 /* mstpdf->update(x1,sqrt(mf2));
  // Then the individual flavours are accessed from the cont structure.
 
  upv = mstpdf->cont.upv;
  dnv = mstpdf->cont.dnv;
  usea = mstpdf->cont.usea;
  dsea = mstpdf->cont.dsea;
  str = mstpdf->cont.str;
  sbar = mstpdf->cont.sbar;
  chm = mstpdf->cont.chm;
  cbar = mstpdf->cont.cbar;
  bot = mstpdf->cont.bot;
  bbar = mstpdf->cont.bbar;
  glu = mstpdf->cont.glu; */
     //  GeT_all_PDFS( x1, sqrt(mf2), upv, dnv, usea, dsea, str, sbar, chm, cbar, bot, bbar, glu,phot,iset,ih1);
             
 
     f1[0]  = fx1[parton::G ];
     f1[2]  = (fx1[parton::U ]);
     f1[1]  = (fx1[parton::D ]);
     f1[3]  = fx1[parton::S ];
     f1[4]  = fx1[parton::C ];
     f1[5]  = fx1[parton::B ];
     f1[6]  = 0.0;
     f1[11] = fx1[parton::Bb ];
     f1[10] = fx1[parton::Cb ];
     f1[9]  = fx1[parton::Sb ];
     f1[7]  = fx1[parton::Db ];
     f1[8]  = fx1[parton::Ub ];
     f1[12] = 0.0; 

      // cout << " f1[0] = " <<  f1[0] << ";"<< endl;
      // cout << " f1[2] = " <<  f1[2] << ";"<< endl;
      // cout << " f1[1] = " <<  f1[1] << ";"<<endl;
      // cout << " f1[3] = " <<  f1[3] << ";"<<endl;
      // cout << " f1[4] = " <<  f1[4] << ";"<<endl;
      // cout << " f1[5] = " <<  f1[5] << ";"<<endl;
      // cout << " f1[6] = " <<  f1[6] << ";"<<endl;
      // cout << " f1[11] = " <<  f1[11] << ";"<<endl;
      // cout << " f1[10] = " <<  f1[10] << ";"<<endl;
      // cout << " f1[9] = " <<  f1[9] << ";"<<endl;
      // cout << " f1[7] = " <<  f1[7] << ";"<<endl;
      // cout << " f1[8] = " <<  f1[8] << ";"<<endl;
      // cout << " f1[12] = " <<  f1[12] << ";"<<endl;

    //----- PROTON 2 
/*   mstpdf->update(x2,sqrt(mf2));
  // Then the individual flavours are accessed from the cont structure.
 
  upv = mstpdf->cont.upv;
  dnv = mstpdf->cont.dnv;
  usea = mstpdf->cont.usea;
  dsea = mstpdf->cont.dsea;
  str = mstpdf->cont.str;
  sbar = mstpdf->cont.sbar;
  chm = mstpdf->cont.chm;
  cbar = mstpdf->cont.cbar;
  bot = mstpdf->cont.bot;
  bbar = mstpdf->cont.bbar;
  glu = mstpdf->cont.glu;  */
//GeT_all_PDFS( x2, sqrt(mf2), upv, dnv, usea, dsea, str, sbar,chm, cbar, bot, bbar, glu,phot,iset,ih2);

     f2[0]  = fx2[parton::G ];
     f2[2]  = fx2[parton::U ];
     f2[1]  = fx2[parton::D ];
     f2[3]  = fx2[parton::S ];
     f2[4]  = fx2[parton::C ];
     f2[5]  = fx2[parton::B ];
     f2[6]  = 0.0;
     f2[11] = fx2[parton::Bb ];
     f2[10] = fx2[parton::Cb ];
     f2[9]  = fx2[parton::Sb ];
     f2[7]  = fx2[parton::Db ];
     f2[8]  = fx2[parton::Ub ];
     f2[12] = 0.0; 


      // cout << " f2[0] = " <<  f2[0] << ";"<< endl;
      // cout << " f2[2] = " <<  f2[2] << ";"<< endl;
      // cout << " f2[1] = " <<  f2[1] << ";"<<endl;
      // cout << " f2[3] = " <<  f2[3] << ";"<<endl;
      // cout << " f2[4] = " <<  f2[4] << ";"<<endl;
      // cout << " f2[5] = " <<  f2[5] << ";"<<endl;
      // cout << " f2[6] = " <<  f2[6] << ";"<<endl;
      // cout << " f2[11] = " <<  f2[11] << ";"<<endl;
      // cout << " f2[10] = " <<  f2[10] << ";"<<endl;
      // cout << " f2[9] = " <<  f2[9] << ";"<<endl;
      // cout << " f2[7] = " <<  f2[7] << ";"<<endl;
      // cout << " f2[8] = " <<  f2[8] << ";"<<endl;
      // cout << " f2[12] = " <<  f2[12] << ";"<<endl;



    // //BORRAR ETIQUETAS Y PONER TODO EN UNA MISMA LINEA	
    // 	  f1[0]=glu/x1; //GLUON
    // 	  f1[2]=upv/x1 +usea/x1;//UVAL+UB
    // 	  f1[1]=dnv/x1 + dsea/x1;//DVAL+DB
    // 	  f1[3]=str/x1;//STRANGE
    // 	  f1[4]=chm/x1;//CHRM
    // 	  f1[5]=bot/x1;//BOT
    // 	  f1[6]=0.0;//TOP 
    // 	  f1[11]=bbar/x1;//BOTB
    // 	  f1[10]=cbar/x1;//CHRMB
    // 	  f1[9]=sbar/x1;//STRGB
    // 	  f1[7]=dsea/x1;//DB
    // 	  f1[8]=usea/x1;//UB
    // 	  f1[12]=0.0;//TOPB

    
    // 	  f2[0]=glu/x2; //GLUON
    // 	  f2[2]=upv/x2 + usea/x2;//UVAL+UB
    // 	  f2[1]=dnv/x2 + dsea/x2;//DVAL+DB
    // 	  f2[3]=str/x2;//STRANGE
    // 	  f2[4]=chm/x2;//CHRM
    // 	  f2[5]=bot/x2;//BOT
    // 	  f2[6]=0.0;//TOP 
    // 	  f2[11]=bbar/x2;//BOTB
    // 	  f2[10]=cbar/x2;//CHRMB
    // 	  f2[9]= sbar/x2;//STRGB
    // 	  f2[7]= dsea/x2;//DB
    // 	  f2[8]= usea/x2;//UB
    // 	  f2[12]=0.0;//TOPB



 
    //----- gluon pdfs -----
    double G1 = f1[0];
    double G2 = f2[0];
    
    //---- up type quarks -----
    double q1, q2, a1, a2;
    double U1 = 0.0, U2 = 0.0, Ub1 = 0.0, Ub2 = 0.0, X = 0.0, Xb = 0.0; 
    double D1 = 0.0, D2 = 0.0, Db1 = 0.0, Db2 = 0.0, Y = 0.0, Yb = 0.0; 
    
    for(unsigned int u = 0; u < nu && u < 3; u++) {
      ia = iub[u];
      iq = iu[u]; 
      q1 = f1[iq]; q2 = f2[iq];
      a1 = f1[ia]; a2 = f2[ia];

      U1 += q1; Ub1 += a1; U2 += q2; Ub2 += a2;
      X += q1*q2 + a1*a2; Xb += q1*a2 + a1*q2;
    }
    
    //----- down type quarks -----
    for(unsigned int d = 0; d < nd && d < 3; d++) {
      ia = idb[d];
      iq = id[d]; 
      q1 = f1[iq]; q2 = f2[iq];
      a1 = f1[ia]; a2 = f2[ia];
      
      D1 += q1; Db1 += a1; D2 += q2; Db2 += a2;
      Y += q1*q2 + a1*a2; Yb += q1*a2 + a1*q2;
    }
        
        
    retval[0]  =G1*G2;
    retval[1]  =(U1 + Ub1)*G2;
    retval[2]  =(D1 + Db1)*G2;
    retval[3]  =G1*(U2 + Ub2);
    retval[4]  =G1*(D2 + Db2);
    retval[5]  =U1*U2 + Ub1*Ub2 - X;
    retval[6]  =X;
    retval[7]  =D1*D2 + Db1*Db2 - Y;
    retval[8]  =Y;
    retval[9]  = U1*D2 + Ub1*Db2;
    retval[10] =D1*U2 + Db1*Ub2;
    retval[11] = U1*Ub2 + Ub1*U2 - Xb;
    retval[12] =Xb;
    retval[13] = D1*Db2 + Db1*D2 - Yb;
    retval[14] =Yb;
    retval[15] =U1*Db2 + Ub1*D2;
    retval[16] = D1*Ub2 + Db1*U2;
     
   
    return retval;
  
  }

 weight_hhc2ph  GeT_pdfCT(double x1, double x2, double mf2, unsigned int nu, unsigned int nd)
  {
	   unsigned int nfl=nu+nd;
	    
	     int iset,ih1,ih2;
	 iset=parpdfs2gamma_.iset2g;
	 ih1= parpdfs2gamma_.ih12g;
	 ih2= parpdfs2gamma_.ih22g;
	  
	 
//	   static double __f100[13], __f200[13];
//	   double *f100 = __f100+6, *f200 = __f200+6;
	
	  
	  
	  
      //PDFs
	 double muf = sqrt(mf2);
      double fx1[2*MAXNF+1],fx2[2*MAXNF+1];
      fdist_(opts.ih1,x1,muf,fx1);
      fdist_(opts.ih2,x2,muf,fx2);
 
	  
	  
	  
	  
	 
	  
    int ia, iq;
    weight_hhc2ph retval;
    
    static  double f1[13],f2[13];
    
    static const int iu[3] = {2,4,6}, id[3] = {1,3,5};
    static const int iub[3] = {8,10,12}, idb[3] = {7,9,11}; 
    
    //----- calculate the pdfs -----
     double upv,dnv,usea,dsea,str,sbar,chm,cbar,bot,bbar,glu,phot;
    //----- PROTON 1 
    
     f1[0] = fx1[parton::G ];
     f1[2] = fx1[parton::U ];
     f1[1] = fx1[parton::D ];
     f1[3] = fx1[parton::S ];
     f1[4] = fx1[parton::C ];
     f1[5] = fx1[parton::B ];
     f1[6] = 0.0;
     f1[11] = fx1[parton::Bb ];
     f1[10] = fx1[parton::Cb ];
     f1[9] = fx1[parton::Sb ];
     f1[7] = fx1[parton::Db ];
     f1[8] = fx1[parton::Ub ];
     f1[12] = 0.0; 

    

    //----- PROTON 2 
/*   mstpdf->update(x2,sqrt(mf2));
  // Then the individual flavours are accessed from the cont structure.
 
  upv = mstpdf->cont.upv;
  dnv = mstpdf->cont.dnv;
  usea = mstpdf->cont.usea;
  dsea = mstpdf->cont.dsea;
  str = mstpdf->cont.str;
  sbar = mstpdf->cont.sbar;
  chm = mstpdf->cont.chm;
  cbar = mstpdf->cont.cbar;
  bot = mstpdf->cont.bot;
  bbar = mstpdf->cont.bbar;
  glu = mstpdf->cont.glu;  */
//GeT_all_PDFS( x2, sqrt(mf2), upv, dnv, usea, dsea, str, sbar,chm, cbar, bot, bbar, glu,phot,iset,ih2);

     f2[0] = fx2[parton::G ];
     f2[2] = fx2[parton::U ];
     f2[1] = fx2[parton::D ];
     f2[3] = fx2[parton::S ];
     f2[4] = fx2[parton::C ];
     f2[5] = fx2[parton::B ];
     f2[6] = 0.0;
     f2[11] = fx2[parton::Bb ];
     f2[10] = fx2[parton::Cb ];
     f2[9] = fx2[parton::Sb ];
     f2[7] = fx2[parton::Db ];
     f2[8] = fx2[parton::Ub ];
     f2[12] = 0.0; 

 
    //----- gluon pdfs -----
    double G1 = f1[0];
    double G2 = f2[0];
    
    //---- up type quarks -----
    double q1, q2, a1, a2;
    double U1 = 0.0, U2 = 0.0, Ub1 = 0.0, Ub2 = 0.0,  Xb = 0.0; 
    double D1 = 0.0, D2 = 0.0, Db1 = 0.0, Db2 = 0.0,  Yb = 0.0; 
    
    for(unsigned int u = 0; u < nu && u < 3; u++) {
      ia = iub[u];
      iq = iu[u]; 
      q1 = f1[iq]; q2 = f2[iq];
      a1 = f1[ia]; a2 = f2[ia];

      U1 += q1; Ub1 += a1; U2 += q2; Ub2 += a2;
       Xb += q1*a2 + a1*q2;
    }
    
    //----- down type quarks -----
    for(unsigned int d = 0; d < nd && d < 3; d++) {
      ia = idb[d];
      iq = id[d]; 
      q1 = f1[iq]; q2 = f2[iq];
      a1 = f1[ia]; a2 = f2[ia];
      
      D1 += q1; Db1 += a1; D2 += q2; Db2 += a2;
     Yb += q1*a2 + a1*q2;
    }
        
        
    retval[0]  =G1*G2;
    retval[1]  =(U1 + Ub1)*G2;
    retval[2]  =(D1 + Db1)*G2;
    retval[3]  =G1*(U2 + Ub2);
    retval[4]  =G1*(D2 + Db2);
    retval[5]  =0.0;//U1*U2 + Ub1*Ub2 - X;
    retval[6]  =0.0;//X;
    retval[7]  =0.0;//D1*D2 + Db1*Db2 - Y;
    retval[8]  =0.0;//Y;
    retval[9]  =0.0;// U1*D2 + Ub1*Db2;
    retval[10] =0.0;//D1*U2 + Db1*Ub2;
    retval[11] =0.0;// U1*Ub2 + Ub1*U2 - Xb;
    retval[12] =Xb;
    retval[13] =0.0;// D1*Db2 + Db1*D2 - Yb;
    retval[14] =Yb;
    retval[15] =0.0;//U1*Db2 + Ub1*D2;
    retval[16] =0.0;// D1*Ub2 + Db1*U2;
     
   
    return retval;
  
  }

weight_hhc2ph  GeT_pdfCT(double x1, double x2, double mf2, unsigned int nu, unsigned int nd,double (&f1)[13], double (&f2)[13])
  {
	   unsigned int nfl=nu+nd;
	    
	     int iset,ih1,ih2;
	 iset=parpdfs2gamma_.iset2g;
	 ih1= parpdfs2gamma_.ih12g;
	 ih2= parpdfs2gamma_.ih22g;
	  
	 
//	   static double __f100[13], __f200[13];
//	   double *f100 = __f100+6, *f200 = __f200+6;
	
	  
	  
	  
 
      //PDFs
	 double muf = sqrt(mf2);
	 double fx1[2*MAXNF+1],fx2[2*MAXNF+1];
	 fdist_(opts.ih1,x1,muf,fx1);
	 fdist_(opts.ih2,x2,muf,fx2);
	  
	  
	  
	  
	 
	  
    int ia, iq;
    weight_hhc2ph retval;
    
 //   static  double f1[13],f2[13];
    
    static const int iu[3] = {2,4,6}, id[3] = {1,3,5};
    static const int iub[3] = {8,10,12}, idb[3] = {7,9,11}; 
    
    //----- calculate the pdfs -----
     double upv,dnv,usea,dsea,str,sbar,chm,cbar,bot,bbar,glu,phot;
    //----- PROTON 1 
    
    	  // Update all PDF flavours.
 /* mstpdf->update(x1,sqrt(mf2));
  // Then the individual flavours are accessed from the cont structure.
 
     f1[0] = fx1[parton::G ];
     f1[2] = fx1[parton::U ]+fx1[parton::Ub ];
     f1[1] = fx1[parton::D ]+fx1[parton::Db ];
     f1[3] = fx1[parton::S ];
     f1[4] = fx1[parton::C ];
     f1[5] = fx1[parton::B ];
     f1[6] = 0.0;
     f1[11] = fx1[parton::Bb ];
     f1[10] = fx1[parton::Cb ];
     f1[9] = fx1[parton::Sb ];
     f1[7] = fx1[parton::Db ];
     f1[8] = fx1[parton::Ub ];
     f1[12] = 0.0; 

    

    //----- PROTON 2 
/*   mstpdf->update(x2,sqrt(mf2));
  // Then the individual flavours are accessed from the cont structure.
 
  upv = mstpdf->cont.upv;
  dnv = mstpdf->cont.dnv;
  usea = mstpdf->cont.usea;
  dsea = mstpdf->cont.dsea;
  str = mstpdf->cont.str;
  sbar = mstpdf->cont.sbar;
  chm = mstpdf->cont.chm;
  cbar = mstpdf->cont.cbar;
  bot = mstpdf->cont.bot;
  bbar = mstpdf->cont.bbar;
  glu = mstpdf->cont.glu;  */
//GeT_all_PDFS( x2, sqrt(mf2), upv, dnv, usea, dsea, str, sbar,chm, cbar, bot, bbar, glu,phot,iset,ih2);

     f1[0] = fx1[parton::G ];
     f1[2] = fx1[parton::U ];
     f1[1] = fx1[parton::D ];
     f1[3] = fx1[parton::S ];
     f1[4] = fx1[parton::C ];
     f1[5] = fx1[parton::B ];
     f1[6] = 0.0;
     f1[11] = fx1[parton::Bb ];
     f1[10] = fx1[parton::Cb ];
     f1[9] = fx1[parton::Sb ];
     f1[7] = fx1[parton::Db ];
     f1[8] = fx1[parton::Ub ];
     f1[12] = 0.0; 




     f2[0] = fx2[parton::G ];
     f2[2] = fx2[parton::U ];
     f2[1] = fx2[parton::D ];
     f2[3] = fx2[parton::S ];
     f2[4] = fx2[parton::C ];
     f2[5] = fx2[parton::B ];
     f2[6] = 0.0;
     f2[11] = fx2[parton::Bb ];
     f2[10] = fx2[parton::Cb ];
     f2[9] = fx2[parton::Sb ];
     f2[7] = fx2[parton::Db ];
     f2[8] = fx2[parton::Ub ];
     f2[12] = 0.0; 


 
    //----- gluon pdfs -----
    double G1 = f1[0];
    double G2 = f2[0];
    
    //---- up type quarks -----
    double q1, q2, a1, a2;
    double U1 = 0.0, U2 = 0.0, Ub1 = 0.0, Ub2 = 0.0,  Xb = 0.0; 
    double D1 = 0.0, D2 = 0.0, Db1 = 0.0, Db2 = 0.0,  Yb = 0.0; 
    
    for(unsigned int u = 0; u < nu && u < 3; u++) {
      ia = iub[u];
      iq = iu[u]; 
      q1 = f1[iq]; q2 = f2[iq];
      a1 = f1[ia]; a2 = f2[ia];

      U1 += q1; Ub1 += a1; U2 += q2; Ub2 += a2;
       Xb += q1*a2 + a1*q2;
    }
    
    //----- down type quarks -----
    for(unsigned int d = 0; d < nd && d < 3; d++) {
      ia = idb[d];
      iq = id[d]; 
      q1 = f1[iq]; q2 = f2[iq];
      a1 = f1[ia]; a2 = f2[ia];
      
      D1 += q1; Db1 += a1; D2 += q2; Db2 += a2;
     Yb += q1*a2 + a1*q2;
    }
        
        
    retval[0]  =G1*G2;
    retval[1]  =(U1 + Ub1)*G2;
    retval[2]  =(D1 + Db1)*G2;
    retval[3]  =G1*(U2 + Ub2);
    retval[4]  =G1*(D2 + Db2);
    retval[5]  =0.0;//U1*U2 + Ub1*Ub2 - X;
    retval[6]  =0.0;//X;
    retval[7]  =0.0;//D1*D2 + Db1*Db2 - Y;
    retval[8]  =0.0;//Y;
    retval[9]  =0.0;// U1*D2 + Ub1*Db2;
    retval[10] =0.0;//D1*U2 + Db1*Ub2;
    retval[11] =0.0;// U1*Ub2 + Ub1*U2 - Xb;
    retval[12] =Xb;
    retval[13] =0.0;// D1*Db2 + Db1*D2 - Yb;
    retval[14] =Yb;
    retval[15] =0.0;//U1*Db2 + Ub1*D2;
    retval[16] =0.0;// D1*Ub2 + Db1*U2;
     
   
    return retval;
  
  }
