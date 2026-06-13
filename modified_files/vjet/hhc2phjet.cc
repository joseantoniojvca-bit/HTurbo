//  Copyright (C) 2002 Zoltan Nagy
//  2012 Leandro Cieri
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
#include <cstring>

//   nlo includes
#include "proc-hhc2ph.h"
#include "nlo-color.h"
#include "nlo-flkern.h"

#include "ampq2g1p2.h"
#include "ampq2g2p2.h"
#include "ampq4p2.h"

#include "H1.h"
#include "CT_coeff.h"
#include <math.h>
#include "Inter_params.h"
#include "defmacros.h"

#include <quadmath.h>

// --- DYTurbo includes
#include "mesq.h"

// ---------------------------------------------------------------------------------------------------------------------------
// ---  (functions called from DYTurbo)
// ---------------------------------------------------------------------------------------------------------------------------
// -- 24/07/2024 -------------------------------------------------------
// -- Todos qg reales NLO incluidos en una unica rutina
// ---------------------------------------------------------------------

void _hhc2ph_jet_base::amp_qgINTg(ampq2g1p2 *amp, double *out)
{
  //   double tmp;
  std::memset(out, 0, 17*sizeof(double));
  double tmpqg,tmpqg2;
  double tmpgq,tmpgq2;
  
  // out[0] = (amp -> wgtHint(0,-1, 1, -2,-3));
  tmpqg  = (amp -> AmpIntqgNLO(0,-1,-3,-2,1));
  tmpqg2 = (amp -> AmpIntqgNLO(0,-1,-2,-3,1));
  tmpgq  = (amp -> AmpIntqgNLO(-1,0,-3,-2,1));
  tmpgq2 = (amp -> AmpIntqgNLO(-1,0,-2,-3,1));
  // in= up[p1], gluon[p2] ;  out= ga[p3], ga[p4], up[p5]   

  // cout << " Event beginning " << endl;
  // cout << " (amp -> AmpIntqgNLO(0,-1,-2,-3,1)) " << (amp -> AmpIntqgNLO(0,-1,-2,-3,1)) << endl;
  // cout << " (amp -> AmpIntqgNLO(0,-1,-3,-2,1)) " << (amp -> AmpIntqgNLO(0,-1,-3,-2,1)) << endl;
  // cout << " (amp -> AmpIntqgNLO(-1,0,-2,-3,1)) " << (amp -> AmpIntqgNLO(-1,0,-2,-3,1)) << endl;
  // cout << " (amp -> AmpIntqgNLO(-1,0,-3,-2,1)) " << (amp -> AmpIntqgNLO(-1,0,-3,-2,1)) << endl;
  // cout << " End event " << endl;
  
  //	  cout << tmpqg << endl;
  // qg ->	  
  // out[1] = tmpqg; 
  // out[2] = tmpqg;
  // out[3] = tmpgq; 
  // out[4] = tmpgq;

  out[1] = (tmpqg+tmpqg2)/2.; 
  out[2] = (tmpqg+tmpqg2)/2.;
  
  out[3] = (tmpgq+tmpgq2)/2.; 
  out[4] = (tmpgq+tmpgq2)/2.;
  
}
// -- 07/07/2024 -------------------------------------------------------
// -- Todos los reales NLO incluidos en una unica rutina
// ---------------------------------------------------------------------
void _hhc2ph_jet_base::amp_ggINTg(ampq2g1p2 *amp, double *out)
{
  //   double tmp;
  std::memset(out, 0, 17*sizeof(double));
  
  //   std::cout << " desde hhc2phjet.cc  amp_ggINTg " << std::endl;
  
  //  gg ->
      if(opts.calc_type=="Int") out[0] = (amp -> wgtHint(0,-1, 1, -2,-3));
      if(opts.calc_type=="Sig") out[0] = (amp -> wgtHsig(0,-1, 1, -2,-3));
      if(opts.calc_type=="SigInt") out[0] = (amp -> wgtHint(0,-1, 1, -2,-3))+(amp -> wgtHsig(0,-1, 1, -2,-3));;
  
}

// ---------------------------------------------------------------------------------------------------------------------------
// ---  CTgg
// ---------------------------------------------------------------------------------------------------------------------------

 void _hhc2ph_jet_base::amp_CTgg(ampq2g1p2 *amp, double *out, event_hhc2ph  p, double muf2,double mur2, unsigned int Nf, double z1, double z2,int order)
 {
     
  // p[-3] = lorentzvector<double> (-37.4177, 18.3014, 100.698, 108.973);
  // p[-2] = lorentzvector<double> (37.4177, -18.3014, -13.109, 43.6678) ;
  // p[-1] = lorentzvector<double> (0, 0, 120.114, 120.114) ;
  // p[0]  = lorentzvector<double> (0, 0, -32.5259, 32.5259);

  // mur2 = 125.009*125.009;
  // muf2 = 125.009*125.009;
  // z1 = 0.1;
  // z2 = 0.2;
  
   // cout <<  " muf2 = " <<  sqrt(muf2) << endl; 
   // cout <<  " mur2 = " <<  sqrt(mur2) << endl; 

   std::memset(out, 0, 17*sizeof(double));
  ////////////// CALCULAR x10 y x20  ////////////////////////
  double s = p[hadron(-1)]*p[hadron(0)];
  double x10 = p[-1]*p[hadron(0)]/s;
  double x20 = p[0]*p[hadron(-1)]/s;
  int ipdf[11] = {11,10,9,7,8,0,2,1,3,4,5};
  double Q2,pT,x1p,x2p,wctnnlo;
  double xmur=sqrt(mur2);
  double xmuf=sqrt(muf2);
  Q2 = 2.0*((p[-2])*(p[-3]));
  pT = (p[1]).perp();

  // cout << " pT = " << pT << endl;
  // cout << " (p[-2]+p[-3]).perp() = " <<  (p[-2]+p[-3]).perp()  << endl;
  // Q2 = 125.009*125.009;
  // pT = 10.;
  double xmio;// = pT/sqrt(Q2);
  if (opts.fmures > 0)
    a_param_.a_param_ = 1./opts.kmures;//for fixed resummation scale need to recompute a_param   
  else
    a_param_.a_param_ = sqrt(Q2)/scales::res;
  xmio = pT*sqrt(1./(Q2/pow(a_param_.a_param_,2)));
  //cout << " a_param_.a_param_ = " << a_param_.a_param_ << endl;
  x1p = x10/z1;
  x2p = x20/z2;
  double LR=log(Q2/mur2);
  double LF=log(Q2/muf2);
  // This parameter should be ON when resummation is ON
  double LQ =0.0;
  if(!opts.fixedorder) LQ = 2.*log(a_param_.a_param_);
  weight_hhc2ph Wctnlofinal;
  weight_hhc2ph _M_pdfL1,_M_pdfL2,_M_pdf00,_M_pdfpp;
  weight_hhc2ph retval,retvalggL1,retvalggL2;
  double mh= opts.hmass;
  double Gammah = opts.hwidth; // Ok
  // cout << "mh : " << mh << endl;
  // cout << "Gammah : " << Gammah << endl;
  double Ds = ((-mh*mh + Q2)*(-mh*mh + Q2) + mh*mh*Gammah*Gammah);  // In Born mesq.C 
  // cout << " Ds from CT : " << Ds << endl;
  // double ampB = _amp_MggHyyB(p);
  double s12 = 2.*p[-1]*p[0];double t = -2.*p[-1]*p[-3]; double u = -2.*p[-1]*p[-2];
  // cout << " sq2 = " << s12 << endl;
  double ampB = mesq::_amp_MggHyy(s12,t,u)/Ds;
  //double ampBmesq = real(mesq::amp_MggHyy(s12,t,u))*Ds;
  double LL1 = Itilde(1,xmio)/(Q2*Q2);//Chequear que sean los de Massi
  double LL2 = Itilde(2,xmio)/(Q2*Q2);
  double f10[13],f20[13],f1p[13],f2p[13]; 
  _M_pdf00=GeT_pdfCT(x10,x20,muf2, 2, 3,f10,f20);
  _M_pdfL1=GeT_pdfCT(x1p,x20,muf2, 2, 3);
  _M_pdfL2=GeT_pdfCT(x10,x2p,muf2, 2, 3);
  _M_pdfpp=GeT_pdfCT(x1p,x2p,muf2, 2, 3,f1p,f2p);
  //CANAL GG
  weight_hhc2ph Waux00;
  retval[0] = ampB;
  Waux00 = retval;

  weight_hhc2ph Wauxdelta;
  Wauxdelta = Waux00*_M_pdf00; // LEAN Ok!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  weight_hhc2ph WtH1stgammadelta,Wgammadeltaterm;
  Wgammadeltaterm[0] = 2.0*beta0(Nf)*ampB;
  WtH1stgammadelta = Wgammadeltaterm *_M_pdf00;// LEAN Ok!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  /////////////////////////////////////////////////////////////////////////////////////////////////////	
  // Scheme selection   
  /////////////////////////////////////////////////////////////////////////////////////////////////////	
  // H1gH or C1ggdelta should be put to zero. In this way you are selecting the scheme.   
  /////////////////////////////////////////////////////////////////////////////////////////////////////	   
  // 
  // HARD SCHEME
  //
  /////////////////////////////////////////////////////////////////////////////////////////////////////	
  //   H1gH=2.d0*(11+3*pi**2)/4d0
  double   H1gH=2.0*(11.0+3.0*__DEB_PI*__DEB_PI)/4.0;
  H1gH = 0.0;
  /////////////////////////////////////////////////////////////////////////////////////////////////////	
  weight_hhc2ph WtH1sthardm,WtH1stdeltahardm,WtH1sthard,WtH1stdeltahard;
  double tH1sthard = H1gH*ampB;//Checked uno por cada pata LEAN
  double tH1stdeltahard = H1gH*ampB;//
  WtH1sthardm[0] = tH1sthard;
  WtH1stdeltahardm[0] = tH1stdeltahard;   
  WtH1sthard = WtH1sthardm*_M_pdf00;
  WtH1stdeltahard = WtH1stdeltahardm*_M_pdf00;
  /////////////////////////////////////////////////////////////////////////////////////////////////////	
  // Old Scheme
  /////////////////////////////////////////////////////////////////////////////////////////////////////	
  double C1ggdelta=(11.0+3.0*__DEB_PI*__DEB_PI)/4.0;
  /////////////////////////////////////////////////////////////////////////////////////////////////////	
  weight_hhc2ph WtH1stoldm,WtH1stold;
  double tH1stold = 2.0*C1ggdelta*ampB;//Checked uno por cada pata LEAN
  WtH1stoldm[0] = tH1stold;
  WtH1stold = WtH1stoldm*_M_pdf00;
  /////////////////////////////////////////////////////////////////////////////////////////////////////	
  WtH1sthard += WtH1stold;
  /////////////////////////////////////////////////////////////////////////////////////////////////////		
  /////////////////////////////////////////////////////////////////////////////////////////////////////	   
  //-----Term Gammaqq coefficient--pata 1 (Esto significa que la
  // pata 1 emite parton) La contribucion con la PDF con xi0 se pone aparte
  /////////////////////////////////////////////////////////////////////////////////////////////////////	
  /////////////////////////////////////////////////////////////////////////////////////////////////////	
  double tH1stqL1 = (3.0/(1.0-z1)+Pggreg(z1))*(-log(x10)); // LEAN Ok!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  double  tH1stqL10 = (-3.0*z1/(1.0-z1))*(-log(x10))-3.0*D0int(x10); // LEAN Ok!!!!!!!!!!!!!!!!!!!!!
  retvalggL1[0] = tH1stqL1*ampB; // LEAN Ok!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  weight_hhc2ph retvalggL10;
  retvalggL10[0] = tH1stqL10*ampB; // LEAN Ok!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  weight_hhc2ph WL1gg;
  WL1gg = retvalggL1*_M_pdfL1 +retvalggL10*_M_pdf00  ;
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////		
  /////////////////////////////////////////////////////////////////////////////////////////////////////	
  //-----Term Gammaqq coefficient--pata 2 (Esto significa que la
  // pata 2 emite parton) La contribucion con la PDF con xi0 se pone aparte
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////		
  /////////////////////////////////////////////////////////////////////////////////////////////////////	
  double  tH1stqL2 = (3.0/(1.0-z2)+Pggreg(z2))*(-log(x20)); // LEAN Ok!!!!!!!!!!!!!!!!!!!!!!!!
  double  tH1stqL20 = (-3.0*z2/(1.0-z2))*(-log(x20))-3.0*D0int(x20); // LEAN Ok!!!!!!!!!!!!!!!
  retvalggL2[0] = tH1stqL2*ampB; // LEAN Ok!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  weight_hhc2ph retvalggL20;
  retvalggL20[0] = tH1stqL20*ampB; // LEAN Ok!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  weight_hhc2ph WL2gg;
  WL2gg = retvalggL2*_M_pdfL2 +retvalggL20*_M_pdf00  ;
  /////////////////////////////////////////////////////////////////////////////////////////////////////	
  /////////////////////////////////////////////////////////////////////////////////////////////////////		
  /////////////////////////////////////////////////////////////////////////////////////////////////////		
  weight_hhc2ph WtH1stm,WtH1stFm;
  double tH1stmL1=0.0,tH1stmL2=0.0,tH1stFmL1=0.0,tH1stFmL2=0.0;
  for(int l=1;l<6;l++)
    {
      // first leg
      tH1stmL1 += -log(x10)*(f1p[ipdf[l+5]]+f1p[ipdf[-l+5]])*Cgq(z1)*ampB;
      // second leg
      tH1stmL2 += -log(x20)*(f2p[ipdf[l+5]]+f2p[ipdf[-l+5]])*Cgq(z2)*ampB;
      // H1st : muf dependence Pgq, first leg
      tH1stFmL1 += -log(x10)*(f1p[ipdf[l+5]]+f1p[ipdf[-l+5]])*Pgq(z1)*ampB;
      // H1st : muf dependence Pgq, second leg
      tH1stFmL2 += -log(x20)*(f2p[ipdf[l+5]]+f2p[ipdf[-l+5]])*Pgq(z2)*ampB;	
    }
  
  WtH1stm[4] = tH1stmL1*f20[0];
  WtH1stm[3] = tH1stmL2*f10[0];
  WtH1stFm[4] = tH1stFmL1*f20[0];
  WtH1stFm[3] = tH1stFmL2*f10[0];	
  /////////////////////////////////////////////////////////////////////////////////////////////////////	
  /////////////////////////////////////////////////////////////////////////////////////////////////////		
  /////////////////////////////////////////////////////////////////////////////////////////////////////	
  weight_hhc2ph  WtH1stF;
  // tH1stF
  WtH1stF = WtH1stFm +  WtH1stgammadelta + WL1gg + WL2gg;
  weight_hhc2ph  WtH1st;
  // tH1st            
  WtH1st= WtH1stm+WtH1sthard;  
  weight_hhc2ph sig12,sig11,sig1;
  //   First order
  sig12=-0.5*A1g*Wauxdelta;
  sig11=(-1.0*B1g*Wauxdelta)+(-1.0*WtH1stF)-A1g*Wauxdelta*LQ;
  
  sig1=sig12*LL2+sig11*LL1;
  Wctnlofinal = ((-1.0)*(sig1));
  for(int i=0;i<17;i++) out[i] = Wctnlofinal[i]; 
  //////////////////////////////////////////////////////////////////////////////// 
  ////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////// 
  // ONLY GG CHANNEL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////// 	
  if(opts.gg_only)
    {
      out[1]=0.0;out[2]=0.0;out[3]=0.0;out[4]=0.0;
    }
  // -- Only qg Channel
  if(opts.qg_only) out[0] = 0.0;

  // cout << p << endl;
  // cout << " s = " << s << endl;
  // cout << " x10 = " << x10 << endl;
  // cout << " x20 = " << x20 << endl;
  // cout << " xmur = " << xmur << endl;
  // cout << " xmuf = " << xmuf << endl;
  // cout << " sqrt(Q2) = " << sqrt(Q2) << endl;
  // cout << "  pT = " << pT << endl;
  // cout << " x1p = " << x1p << endl;
  // cout << " x2p = " << x2p << endl;
  // cout << " ampB = " << ampB << endl;
  // cout << " LL1 = " << LL1 << endl;
  // cout << " LL2 = " << LL2 << endl;
  // cout << " Nf : " << Nf << endl;
  // cout << " Waux00 : " << Waux00 << endl;
  // cout << "_M_pdf00 : " << _M_pdf00 << endl; 
  // cout << " WtH1stgammadelta ; " << WtH1stgammadelta << endl;
  // cout << " out : " << *out << endl;

}



// ---------------------------------------------------------------------------------------------------------------------------
// QQBAR and QG CHANNELS
// ---------------------------------------------------------------------------------------------------------------------------
 void _hhc2ph_jet_base::amp_MqqbarHyy(ampq2g1p2 *amp, double *out)
  {
    //   double tmp;
    std::memset(out, 0, 17*sizeof(double));
    double tmpqgqg,tmpqggq,tmpqqbar;
    tmpqgqg  = (amp -> LO_MqgH(-1,0, 1,-2,-3));
    tmpqggq  = (amp -> LO_MqgH(0,-1, 1,-2,-3));
    tmpqqbar = (amp -> LO_MqqbarH(-1, 0, 1,-2,-3));
    double UU,DD,UG,DG,GU,GD;
    UU = (2.0/3.0)*(2.0/3.0);
    DD = (1.0/3.0)*(1.0/3.0);
    UG = (2.0/3.0)*(2.0/3.0);
    DG = (1.0/3.0)*(1.0/3.0);
    GU = (2.0/3.0)*(2.0/3.0);
    GD = (1.0/3.0)*(1.0/3.0);

    if(opts.LO_qg_ON)
      {    
	out[1] =UG* tmpqgqg; //original implementation from Beta5 and Beta4 version
	out[2] =DG* tmpqgqg;
	out[3] =GU* tmpqggq;
	out[4] =GD* tmpqggq;

	// out[1] =UG* tmpqggq;
	// out[2] =DG* tmpqggq;
	// out[3] =GU* tmpqgqg;
	// out[4] =GD* tmpqgqg;
      }
    
    if(opts.LO_qqbar_ON)
      {
	out[12] =UU* tmpqqbar;
	out[14] =DD* tmpqqbar;
      }
    // No contributions from other channels    
  }
// ---------------------------------------------------------------------------------------------------------------------------













// ---------------------------------------------------------------------------------------------------------------------------
// ---  (Not used functions to be called from DYTurbo)
// ---------------------------------------------------------------------------------------------------------------------------

void _hhc2ph_jet_base::amp_LO(ampq2g1p2 *amp, double *out)
  {
    double tmp;
    std::memset(out, 0, 17*sizeof(double));

    //  gg ->

    //  ug -> u    (1/Nc * 1/Na)
    //  dg -> d    (1/Nc * 1/Na)


    //  gu -> u    (1/Na * 1/Nc)
    //  gd -> d    (1/Na * 1/Nc)


    //  uU ->
    //  uu ->
    //  dD ->
    //  dd ->
    //  ud ->
    //  du ->
    //  uUb ->

    //  uub -> g   (1/Nc * 1/Nc)
    out[12] = Qu4*(tmp = (amp -> LO_tree(-1,0, 1, -2,-3))/Nc2);

    //  dDb ->
    //  ddb -> g   (1/Nc * 1/Nc)
    out[14] = Qd4*tmp;

  //for(int i=0;i<17;i++) out[i]=1.0;

    //  udb ->
    //  dub ->
  }


void _hhc2ph_jet_base::amp_CTNLO(ampq2g1p2 *amp, double *out, event_hhc2ph  p, double muf2,double mur2, unsigned int Nf,
 double z1, double z2,int order){
// double tmp;
    std::memset(out, 0, 17*sizeof(double));

//for(int i=0;i<17;i++) out[i] = 1.0;
//	return;
////////////// CALCULAR x10 y x20  ////////////////////////
	double s = p[hadron(-1)]*p[hadron(0)];
    double x10 = p[-1]*p[hadron(0)]/s;
    double x20 = p[0]*p[hadron(-1)]/s;
	int ipdf[11] = {11,10,9,7,8,0,2,1,3,4,5};





double Q2,pT,x1p,x2p,wctnnlo;
	double xmur=sqrt(mur2);
	double xmuf=sqrt(muf2);
	Q2 = 2.0*((p[-2])*(p[-3]));
	pT = (p[1]).perp();
	double xmio = pT/sqrt(Q2);
	x1p = x10/z1;
	x2p = x20/z2;
	double LR=log(Q2/mur2);
    double LF=log(Q2/muf2);

	//Peso del contratérmino a segundo orden (NNNLO) canales qq y qg




	weight_hhc2ph Wctnlofinal,Wqq,WTL1,WTL2;






//El order siempre debe estar en 2 pues el orden 1 se calcula directamente en c++
//El peso de orden 2 (NNLO) se pasa por referencia al main_c++

	// Prepara las PDF's para cada una de las patas


	weight_hhc2ph _M_pdfL1,_M_pdfL2,_M_pdf00,_M_pdfpp;

	weight_hhc2ph retval,retvalqqL1,retvalqqL2,retvalqqL11,retvalqqL22,retval34,retval12;
	weight_hhc2ph M_CT;

	//Calcula el BORN



	 double ampB=ampBORN(p);
	 double ah1 = H1(p)/2;
//	 ah1 =0.0;
	// El factor se explica mas abajo
	// es porque en el paper 0508068 se desarrolla en potencias de alpha_s/Pi
	 double C1qqdelta;// = H1(p)/2.0/2.0;
	 C1qqdelta = 0.0;
//	 double  C1qqdeltaDY=(__DEB_PI*__DEB_PI-8.0)/3.0;
//	 double B2q_gg= (4.0/9*(__DEB_PI*__DEB_PI-3.0/4-12*Z3)+(11.0/9*__DEB_PI*__DEB_PI-193.0/12+6*Z3)+Nf/6.0*(17.0/3-4.0/9*__DEB_PI*__DEB_PI))+(11.0/3.0-Nf/3)*(-(__DEB_PI*__DEB_PI*2.0/3.0-8.0)/4.0*4.0/3.0  + (ah1/2.0/2.0 -__DEB_PI*__DEB_PI/9.0 ));

      //Schema viejo
//	 double B2qNf= -3.39493 + 2.0*C1qqdelta*(2.75 - 0.166667*Nf) + 0.421096*Nf + 0.666667*Z3;
        // Hard Scheme
     double B2qNf=(1.0/48.0)*_CFH*(2.0*Nf*(1.0 + 4.0*Z2) +    9.0* _CFH*(-1.0 + 8.0* Z2 - 16.0*Z3) +    3.0*(-17.0 - 44.0*Z2 + 72.0*Z3));

	// 28/10/2015
// Cambio en el B2q. La suma de B2q con C1qq (que incluye el C1qqdelta en el esquema correspondiente)
// debe ser independiente del proceso.
//

	double LL1 = Itilde(1,xmio)/(Q2*Q2);//Chequear que sean los de Massi
	double LL2 = Itilde(2,xmio)/(Q2*Q2);
	double LL3 = Itilde(3,xmio)/(Q2*Q2);//Chequear que sean los de Massi
	double LL4 = Itilde(4,xmio)/(Q2*Q2);



//	  if(order==2){
//	wctnnlo=CTNNLO(Q2,x10, x20, x1p
 //   ,x2p, xmur, xmuf, xmio, 83, ampB*__alfaQED2 , order,
  //    LL1, LL2, LL3, LL4,ah1);


 //     for(int i=0;i<17;i++) out[i]=0.0;
 //      out[0] = wctnnlo;
 //     return;

  //    }


//	if(order==1){


		double f10[13],f20[13],f1p[13],f2p[13];

	_M_pdf00=GeT_pdfCT(x10,x20,muf2, 2, 3,f10,f20);
	_M_pdfL1=GeT_pdfCT(x1p,x20,muf2, 2, 3);
	_M_pdfL2=GeT_pdfCT(x10,x2p,muf2, 2, 3);
	_M_pdfpp=GeT_pdfCT(x1p,x2p,muf2, 2, 3,f1p,f2p);

	// vuelve a la notacion antigua de primero u y despues d
	// No hace falta, se hace con los indices directamente
//	swap(f1[1],f1[2]);
//	swap(f2[1],f2[2]);
//	swap(f1[7],f1[8]);
//	swap(f2[7],f2[8]);



	//CANAL QQB
	weight_hhc2ph Waux00,WauxL1qq,WauxL2qq,WauxL1qq0,WauxL2qq0;


	// Témino de la Delta

	retval[12] = Qu4*ampB;
	retval[14] = Qd4*ampB;

	Waux00 = retval;
//Terminos sin convoluciones
	weight_hhc2ph Wauxdelta;

	Wauxdelta = Waux00*_M_pdf00;





    //-----Term Gammaqq coefficient--pata 1 (Esto significa que la
	// pata 1 emite parton) La contribucion con la PDF con xi0 se pone aparte

	double tH1stqL1;
    tH1stqL1 = Pqq(z1)*(-log(x10));

	retvalqqL1[12] =tH1stqL1*Qu4*ampB;
	//Procesos D1 Db2
	retvalqqL1[14] =tH1stqL1*Qd4*ampB;

	WauxL1qq = retvalqqL1;
	weight_hhc2ph WL1qq;
	WL1qq = WauxL1qq*_M_pdfL1;

	//-----Term Gammaqq coefficient--pata 2 (Esto significa que la
	// pata 1 emite parton) La contribucion con la PDF con xi0 se pone aparte

	double tH1stqL2;
    tH1stqL2 = Pqq(z2)*(-log(x20));

	retvalqqL2[12] =tH1stqL2*Qu4*ampB;
	//Procesos D1 Db2
	retvalqqL2[14] =tH1stqL2*Qd4*ampB;

	WauxL2qq =retvalqqL2;
	weight_hhc2ph WL2qq;
	WL2qq = WauxL2qq*_M_pdfL2;

	//PATA 1 para PDF's eveluadas en xi0


   double tH1stqL11 = -Pqq(z1)*(-log(x10))*z1-Pqqint(x10);

	retvalqqL11[12] =tH1stqL11*Qu4*ampB;
	//Procesos D1 Db2
	retvalqqL11[14] =tH1stqL11*Qd4*ampB;

	WauxL1qq0 = retvalqqL11;
	weight_hhc2ph WL1qq0;
	WL1qq0 = WauxL1qq0*_M_pdf00;

	//PATA 2 para PDF's eveluadas en xi0


   double tH1stqL22 = -Pqq(z2)*(-log(x20))*z2-Pqqint(x20);

	retvalqqL22[12] =tH1stqL22*Qu4*ampB;
	//Procesos D1 Db2
	retvalqqL22[14] =tH1stqL22*Qd4*ampB;

	WauxL2qq0 =retvalqqL22;
	weight_hhc2ph WL2qq0;
	WL2qq0 = WauxL2qq0*_M_pdf00;

	// Comenzamos a sumar las contribuciones ---Menos la Delta!!!!
	weight_hhc2ph Wqqaux;

	Wqqaux = WL1qq + WL2qq + WL1qq0 + WL2qq0;
	//El factor 2 Extra se debe a que  alpha_s/2/Pi USO!!!!!!!!!
	//Y Massi usa alpha_s/Pi
	Wqq = ((2.0*LL1*B1q)*Wauxdelta) + ((2.0*LL1)*Wqqaux)
	          +((2.0*0.5*A1q*LL2)*Wauxdelta);
//	Wqq =  realweightprod((2.0*0.5*A1q*LL2),Wauxdelta);
//	cout << LL2 << endl;

	// Fin CANAL QQBAR a Primer orden












	//CANAL QG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//Vamos con la pata 1, es decir procesos del tipo G1 * (Q2)
	//Es decir van en los canales 3 y 4
	double tH1stL1;
	//H1st: non delta terms, pata 1:

//	tH1stL1 = Cqg(z1)*(-log(x10))*0.0;


//	retval[3]=tH1stL1*Qu4*ampB;
//	retval[4]=tH1stL1*Qd4*ampB;


	//Vamos con la pata 2, es decir procesos del tipo Q1 * (G2)
	//Es decir van en los canales 1 y 2
	double tH1stL2;
	//H1st: non delta terms, pata 1:

//	tH1stL2 = Cqg(z2)*(-log(x20))*0.0;


//	retval[1]=tH1stL2*Qu4*ampB;
//	retval[2]=tH1stL2*Qd4*ampB;




	//-----Term Gammaqg coefficient--pata 1
	tH1stL1 = Pqg(z1)*(-log(x10));

	retval34[3]=tH1stL1*Qu4*ampB*LL1;
	retval34[4]=tH1stL1*Qd4*ampB*LL1;

	//-----Term Gammaqg coefficient--pata 2
	tH1stL2 = Pqg(z2)*(-log(x20));


	retval12[1]=tH1stL2*Qu4*ampB*LL1;
	retval12[2]=tH1stL2*Qd4*ampB*LL1;

	//---Pata 1 por PDF pata 1

	weight_hhc2ph WauxL1,WauxL2;

	WauxL1=retval34;
	WauxL2=retval12;

	//El factor 2 Extra se debe a que las pdfs vienen con alpha_s/2/Pi
	//Y Massi usa alpha_s/Pi
	WTL1 =((2.0)* (_M_pdfL1*WauxL1));
	WTL2 =((2.0)*(_M_pdfL2*WauxL2));

	Wctnlofinal=Wqq+WTL1+WTL2;

// Parte a NNLO

 if(order==2){

// H1st delta term
weight_hhc2ph tH1st00,tH1stL1p0,tH1stL20p,tH1st,tH1stsust;
double tmpst;



//------------------Eschema vecchio ----------------
//------------------


//tH1st00[12] = (tmpst=2.0*C1qqdelta*ampB)*Qu4;
//tH1st00[14] = tmpst*Qd4;
//tH1st00 = tH1st00*_M_pdf00;


// 6 de Noviembre de 2015
// En el hard scheme C1qqdelta =0 pero debemos incluir el valor de
// H1qgamgam que es este ah1
// ah1 esta escrito en el hard scheme
// ah1/2 = H1notagian
// en este CTNNLO las cantidades se desarrollan en potencias de alpha_s/Pi
// por lo tanto hay que dividir por dos

tH1st00[12] = (tmpst=ah1*ampB)*Qu4;
tH1st00[14] = tmpst*Qd4;
tH1st00 = tH1st00*_M_pdf00;
tH1stsust = tH1st00;



	//H1st non delta terms
	// NO entra a orden 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//Lo dejo escrito para más adelante
	//Pata 1!!!!
tH1stL1p0[12]=(tmpst=Cqq(z1)*(-log(x10))*ampB)*Qu4;
tH1stL1p0[14]=tmpst*Qd4;

tH1stL1p0[3]=(tmpst=Cqg(z1)*(-log(x10))*ampB)*Qu4;
tH1stL1p0[4]=tmpst*Qd4;
tH1stL1p0 = tH1stL1p0*_M_pdfL1;
    // Pata 2!!!
tH1stL20p[12]=(tmpst=Cqq(z2)*(-log(x20))*ampB)*Qu4;
tH1stL20p[14]=tmpst*Qd4;

tH1stL20p[1]=(tmpst=Cqg(z2)*(-log(x20))*ampB)*Qu4;
tH1stL20p[2]=tmpst*Qd4;
tH1stL20p = tH1stL20p*_M_pdfL2;

tH1st = tH1stL20p + tH1stL1p0 + tH1st00;

// Now (gamma+gamma)*(gamma+gamma) term: to be used later
double diffg1fp,diffg1f0,diffg1gp,diffg2fp,diffg2f0,diffg2gp;
weight_hhc2ph tgagapp,tgaga00,tgagap0,tgaga0p,tgaga;

// First part: one gamma for each leg:

diffg1fp = (-log(x10))*Pqq(z1);
diffg1gp = (-log(x10))*Pqg(z1);
diffg1f0 = (-log(x10))*(-z1*Pqq(z1))-Pqqint(x10);

diffg2fp = (-log(x20))*Pqq(z2);
diffg2gp = (-log(x20))*Pqg(z2);
diffg2f0 = (-log(x20))*(-z2*Pqq(z2))-Pqqint(x20);





/*
// Terminos pp
for(int ljj=-5;ljj<5+1;ljj++){
	for(int lkk=-5;lkk<5+1;lkk++){

	if(Msqc( ljj,  lkk)!=0.0){
tgagapp[0] += diffg1gp*diffg2gp*Msqc( ljj,  lkk)*ampB;
}
}
}
*/


tgagapp[0] = diffg1gp*diffg2gp*ampB*Qsum4*2.0; //Checkeado consistencia, falta el vero!!!!
tgagapp[1] = (tmpst=diffg1fp*diffg2gp*ampB)*Qu4;
tgagapp[2] = tmpst*Qd4;
tgagapp[3] = (tmpst=diffg1gp*diffg2fp*ampB)*Qu4;
tgagapp[4] = tmpst*Qd4;
tgagapp[12] = (tmpst=diffg1fp*diffg2fp*ampB)*Qu4;
tgagapp[14] = tmpst*Qd4;
tgagapp = tgagapp*_M_pdfpp;

// Terminos 0p
tgaga0p[1]= (tmpst=diffg1f0*diffg2gp*ampB)*Qu4;
tgaga0p[2]= tmpst*Qd4;
tgaga0p[12]=(tmpst=diffg1f0*diffg2fp*ampB)*Qu4;
tgaga0p[14]=tmpst*Qd4;
tgaga0p = tgaga0p*_M_pdfL2;

// Terminos p0
tgagap0[3]= (tmpst=diffg2f0*diffg1gp*ampB)*Qu4;
tgagap0[4]= tmpst*Qd4;
tgagap0[12]=(tmpst=diffg1fp*diffg2f0*ampB)*Qu4;
tgagap0[14]=tmpst*Qd4;
tgagap0 = tgagap0*_M_pdfL1;

// Terminos 00
tgaga00[12]=(tmpst=diffg1f0*diffg2f0*ampB)*Qu4;
tgaga00[14]=tmpst*Qd4;
tgaga00 = tgaga00*_M_pdf00;


// Radiogaga fuera del loop Checkeado!!!
tgaga = 2.0*(tgaga00 + tgaga0p +tgagap0 + tgagapp);

//     Second part: gamma*gamma terms
//c     Pij * Pjk = D1ijjk (log(1-z)/(1-z))_+ + D0ijjk/(1-z)_+
//c              + Pijjk(z) + Deltaijjk delta(1-z)
//C     First leg
double diff100,diff1p0,diff1p0g,diff200,diff20p,diff20pg;
weight_hhc2ph diff1loop00,diff2loop00,diff1loopp0,diff2loop0p;

diff100=(-log(x10))*(-z1*(D0qqqq/(1-z1) + D1qqqq*log(1-z1)/(1-z1)))+(Deltaqqqq -D0qqqq*D0int(x10) - D1qqqq*D1int(x10));
diff1loop00[12] = (tmpst=diff100*ampB)*Qu4;
diff1loop00[14] = tmpst*Qd4;
diff1loop00 = diff1loop00*_M_pdf00;

diff1p0 = (-log(x10))*((D0qqqq/(1-z1) + D1qqqq*log(1-z1)/(1-z1)) + Pqqqq(z1));
diff1p0g = (-log(x10))* (Pqqqg(z1)+Pqggg(z1,Nf));

diff1loopp0[3] = (tmpst=diff1p0g*ampB)*Qu4;
diff1loopp0[4] = tmpst*Qd4;
diff1loopp0[12] = (tmpst=diff1p0*ampB)*Qu4;
diff1loopp0[14] = tmpst*Qd4;
diff1loopp0 = diff1loopp0*_M_pdfL1;

//Second leg
diff200=(-log(x20))*(-z2*(D0qqqq/(1-z2) + D1qqqq*log(1-z2)/(1-z2)))+(Deltaqqqq -D0qqqq*D0int(x20) - D1qqqq*D1int(x20));
diff2loop00[12] = (tmpst=diff200*ampB)*Qu4;
diff2loop00[14] = (tmpst)*Qd4;
diff2loop00 = diff2loop00*_M_pdf00;

diff20p = (-log(x20))*((D0qqqq/(1-z2) + D1qqqq*log(1-z2)/(1-z2)) + Pqqqq(z2));
diff20pg = (-log(x20))* (Pqqqg(z2)+Pqggg(z2,Nf));

diff2loop0p[1] = (tmpst=diff20pg*ampB)*Qu4;
diff2loop0p[2] = tmpst*Qd4;
diff2loop0p[12] = (tmpst=diff20p*ampB)*Qu4;
diff2loop0p[14] = tmpst*Qd4;
diff2loop0p = diff2loop0p*_M_pdfL2;



//     Include Pqggq
tgaga = tgaga + diff2loop0p + diff2loop00 + diff1loopp0 + diff1loop00;

// Tambien esta parte de Radio gaga checkeada!!!!

//    End of (gamma+gamma)*(gamma+gamma) term: FLGQ non trivial here ! DONE

//    Start  (C+C)*(gamma+gamma) term

//    gamma first leg, C second leg
double diffc2f0,diffc2fp,diffc20gp,diffc1f0,diffc1fp,diffc10gp;
weight_hhc2ph tcgapp,tcga00,tcgap0,tcga0p,tcga;

// First part: one gamma for each leg:

diffc2fp = (-log(x20))*Cqq(z2);
diffc2f0 = C1qqdelta;
diffc20gp = (-log(x20))*Cqg(z2);

diffc1fp = (-log(x10))*Cqq(z1);
diffc1f0 = C1qqdelta;
diffc10gp = (-log(x10))*Cqg(z1);

// Terminos pp
tcgapp[0] = diffg1gp*diffc20gp*ampB*Qsum4*2.0;
tcgapp[1] = (tmpst=diffg1fp*diffc20gp*ampB)*Qu4;
tcgapp[2] = tmpst*Qd4;
tcgapp[3] = (tmpst=diffg1gp*diffc2fp*ampB)*Qu4;
tcgapp[4] = tmpst*Qd4;
tcgapp[12] =(tmpst=diffg1fp*diffc2fp*ampB)*Qu4;
tcgapp[14] =tmpst*Qd4;
//tcgapp = tcgapp*_M_pdfpp;

// Terminos 0p
tcga0p[1]= (tmpst=diffg1f0*diffc20gp*ampB)*Qu4;
tcga0p[2]= tmpst*Qd4;
tcga0p[12]=(tmpst=diffg1f0*diffc2fp*ampB)*Qu4;
tcga0p[14]=tmpst*Qd4;
//tcga0p = tcga0p*_M_pdfL2;

// Terminos p0
tcgap0[3]= (tmpst=diffc2f0*diffg1gp*ampB)*Qu4;
tcgap0[4]= tmpst*Qd4;
tcgap0[12]=(tmpst=diffg1fp*diffc2f0*ampB)*Qu4;
tcgap0[14]=tmpst*Qd4;
//tcgap0 = tcgap0*_M_pdfL1;

// Terminos 00
tcga00[12]=(tmpst=diffg1f0*diffc2f0*ampB)*Qu4;
tcga00[14]=tmpst*Qd4;
//tcga00 = tcga00*_M_pdf00;
//tcga = tcga00 + tcga0p +tcgap0 + tcgapp;

// First leg, gamma second leg
/*
diffc1fp = (-log(x10))*Cqq(z1);
diffc1f0 = C1qqdelta;
diffc10gp = (-log(x10))*Cqg(z1);
*/
// Terminos pp
tcgapp[0] += diffc10gp*diffg2gp*ampB*Qsum4*2.0;
tcgapp[1] += (tmpst=diffc1fp*diffg2gp*ampB)*Qu4;
tcgapp[2] += tmpst*Qd4;
tcgapp[3] += (tmpst=diffc10gp*diffg2fp*ampB)*Qu4;
tcgapp[4] += tmpst*Qd4;
tcgapp[12] +=(tmpst=diffc1fp*diffg2fp*ampB)*Qu4;
tcgapp[14] +=tmpst*Qd4;
tcgapp = tcgapp*_M_pdfpp;

// Terminos 0p
tcga0p[1]+= (tmpst=diffc1f0*diffg2gp*ampB)*Qu4;
tcga0p[2]+= tmpst*Qd4;
tcga0p[12]+=(tmpst=diffc1f0*diffg2fp*ampB)*Qu4;
tcga0p[14]+=tmpst*Qd4;
tcga0p = tcga0p*_M_pdfL2;

// Terminos p0
tcgap0[3]+= (tmpst=diffg2f0*diffc10gp*ampB)*Qu4;
tcgap0[4]+= tmpst*Qd4;
tcgap0[12]+=(tmpst=diffc1fp*diffg2f0*ampB)*Qu4;
tcgap0[14]+=tmpst*Qd4;
tcgap0 = tcgap0*_M_pdfL1;

// Terminos 00
tcga00[12]+=(tmpst=diffc1f0*diffg2f0*ampB)*Qu4;
tcga00[14]+=tmpst*Qd4;
tcga00 = tcga00*_M_pdf00;
tcga = tcga00 + tcga0p +tcgap0 + tcgapp;


//   C*gamma: first leg (ignore delta term in Cqq: taken into account with tH1stF)

weight_hhc2ph tcgamulp0,tcgamul0p;


tcgamulp0[3]= (tmpst=(-log(x10))*(CqqPqg(z1)+CqgPgg(z1,Nf))*ampB)*Qu4;
tcgamulp0[4]= tmpst*Qd4;
tcgamulp0[12]= (tmpst=(-log(x10))*CqqPqq(z1)*ampB)*Qu4;
tcgamulp0[14]= tmpst*Qd4;
tcgamulp0 = tcgamulp0*_M_pdfL1;

// c    C*gamma: second leg (ignore delta term in Cqq: taken into account with tH1stF)

tcgamul0p[1]= (tmpst=(-log(x20))*(CqqPqg(z2)+CqgPgg(z2,Nf))*ampB)*Qu4;
tcgamul0p[2]= tmpst*Qd4;
tcgamul0p[12]= (tmpst=(-log(x20))*CqqPqq(z2)*ampB)*Qu4;
tcgamul0p[14]= tmpst*Qd4;
tcgamul0p = tcgamul0p*_M_pdfL2;


tcga += tcgamulp0 + tcgamul0p;


//double Kappa= 67.0/6.0-(__DEB_PI*__DEB_PI)/2.0-5.0/9.0*5.0;

// tcga chequeado!!!!

//  Start 2-loop AP
weight_hhc2ph tgamma2p0,tgamma20p,tgamma2,tgamma00;

tgamma2p0[3] = (tmpst=(-log(x10))*P2qg(z1)*ampB)*Qu4;
tgamma2p0[4] = tmpst*Qd4;

tgamma2p0[12] = (tmpst=(-log(x10))*P2qqV(z1,Nf)*ampB)*Qu4;
tgamma2p0[14] = tmpst*Qd4;
tgamma2p0[12] += (tmpst=(-log(x10))/(1-z1)*2.0/3.0*Kappa(Nf)*ampB)*Qu4;
tgamma2p0[14] += tmpst*Qd4;

tgamma2p0 = tgamma2p0*_M_pdfL1;

tgamma20p[1] =(tmpst=(-log(x20))*P2qg(z2)*ampB)*Qu4;
tgamma20p[2] =tmpst*Qd4;

tgamma20p[12] =(tmpst=(-log(x20))*P2qqV(z2,Nf)*ampB)*Qu4;
tgamma20p[14] =tmpst*Qd4;
tgamma20p[12] += (tmpst=(-log(x20))/(1-z2)*2.0/3.0*Kappa(Nf)*ampB)*Qu4;
tgamma20p[14] += tmpst*Qd4;

tgamma20p = tgamma20p*_M_pdfL2;

tgamma00[12] = (tmpst=((-log(x10))*(-z1)/(1-z1)-D0int(x10))*2.0/3.0*Kappa(Nf)*ampB)*Qu4;
tgamma00[14] = tmpst*Qd4;

tgamma00[12] += (tmpst=((-log(x20))*(-z2)/(1-z2)-D0int(x20))*2.0/3.0*Kappa(Nf)*ampB)*Qu4;
tgamma00[14] += tmpst*Qd4;
tgamma00 = tgamma00*_M_pdf00;

//cout << Kappa(Nf) << endl;

tgamma2 = tgamma20p + tgamma2p0 +tgamma00;


// tgamma2 Chequeado!!! mucha dispersion, en ambos codigos tarda en converger y los errores
// estan subestimados

// Ahora metemos los canales qq,qq',etc que van a ir en el canal 10 al final
// Todos estos términos están calculados dentro del born

 /* static const int jdo[11] = {-5,-4,-3,-2,-1,0,1,2,3,4,5};
 static const int kdo[11] = {-5,-4,-3,-2,-1,0,1,2,3,4,5};
*/

double carga=0.0;
double mdiff1=0.0,mdiff2=0.0,mtgaga=0.0,mtcga=0.0,mtgamma2=0.0;
double eq2[6] = {0.0,4.0/9.0,1.0/9.0,1.0/9.0,4.0/9.0,1.0/9.0};
for(int jj=-5;jj<6;jj++){
	for(int kk=-5;kk<6;kk++){

	if(Msqc( jj,  kk)!=0.0){

//	cout << "Msqc("<<jj<<","<<kk<<")  ="<<Msqc( jj,  kk) << endl;


//	cout << "Msqc("<<jj<<","<<kk<<")  ="<<Msqc( jj,  kk) <<"  jj   : "<< jj <<"    ipdf[jj+5]  :  " <<ipdf[jj+5]<< "     f10[ipdf[jj+5]]    :" << f10[ipdf[jj+5]] << endl;

//	carga += eq2[abs(jj)]*Msqc( jj,  kk)*ampB*f10[ipdf[jj+5]]*f20[ipdf[kk+5]] ;

	//*f20[ipdf[kk+5]]*f10[ipdf[jj+5]];

	mdiff1=mdiff2=0.0;

	// Include Pqggq !Diff contiene todo lo anterior tambien!!!
	for(int l=1;l<6;l++){
	mdiff1 += - log(x10)*(f1p[ipdf[l+5]]+f1p[ipdf[-l+5]])*Pqggq(z1);
	mdiff2 += - log(x20)*(f2p[ipdf[l+5]]+f2p[ipdf[-l+5]])*Pqggq(z2);
}
    mtgaga += mdiff1*f20[ipdf[kk+5]]*Msqc( jj,  kk)*ampB;
    mtgaga += mdiff2*f10[ipdf[jj+5]]*Msqc( jj,  kk)*ampB;


	for(int ll=1;ll<6;ll++){
		mtcga=mtcga+(f1p[ipdf[ll+5]]+f1p[ipdf[-ll+5]])*CqgPgq(z1)*(-log(x10))*f20[ipdf[kk+5]]*Msqc( jj,  kk)*ampB;
	    mtcga=mtcga+(f2p[ipdf[ll+5]]+f2p[ipdf[-ll+5]])*CqgPgq(z2)*(-log(x20))*f10[ipdf[jj+5]]*Msqc( jj,  kk)*ampB;
}


   for(int lll=-5;lll<6;lll++){
  if(lll!=0){

  mtgamma2 = mtgamma2 + f1p[ipdf[lll+5]]*P2qqS(z1)*(-log(x10))*f20[ipdf[kk+5]]*Msqc( jj,  kk)*ampB;
  mtgamma2 = mtgamma2 + f2p[ipdf[lll+5]]*P2qqS(z2)*(-log(x20))*f10[ipdf[jj+5]]*Msqc( jj,  kk)*ampB;

}
}

mtgamma2 = mtgamma2 + f1p[ipdf[-jj+5]]*P2qqbV(z1)*(-log(x10))*f20[ipdf[kk+5]]*Msqc( jj,  kk)*ampB;
mtgamma2 = mtgamma2 + f2p[ipdf[-kk+5]]*P2qqbV(z2)*(-log(x20))*f10[ipdf[jj+5]]*Msqc( jj,  kk)*ampB;
}
}
}

//cout << "carga   :" << carga << endl;

tgaga[10] = mtgaga;
tcga[10] = mtcga;
tgamma2[10] = mtgamma2;


//weight_hhc2ph wcarga;

//wcarga[1] =mtgaga;

//cout << "wcarga   :" << wcarga << endl;

weight_hhc2ph sig12,sig11,sig24,sig23,sig22,sig21,sig1,sig2,tH1stF;
tH1stF = (Wqqaux) + (WTL1 + WTL2)/(2.0*LL1);

//   First order
      sig12=-0.5*A1q*Wauxdelta;
      sig11=(-1.0*B1q*Wauxdelta)+(-1.0*tH1stF);


//   Second order
      sig24=(A1q)*(A1q)/8*Wauxdelta;// Checkeado!!!

      sig23=((-beta0(Nf)*A1q/3)*Wauxdelta)+((-0.5*A1q)*sig11);// Checkeado!!!
      sig22=(0.5*(beta0(Nf)*A1q*LR-A2q(Nf)))*Wauxdelta+
          (-0.5*A1q*(tH1st+LF*tH1stF))+
          (-0.5*(B1q-beta0(Nf))*sig11)
          +0.5*B1q*tH1stF
          +0.5*tgaga;// Checkeado!!!
      sig21=((-1.0*beta0(Nf)*LR)*sig11)+(-1.0*B1q*(tH1st+(LF*tH1stF)))+
          ((-1.0*LF)*tgaga)+((-1.0*B2qNf)*Wauxdelta)+(beta0(Nf)*(tH1st-tH1stsust))+((-1.0)*tcga)+((-1.0)*tgamma2); // Checkeado!!!
//     Include missing delta term from C*gamma (no factor 2 here !)
//      sig21=sig21+((-1.0*C1qqdelta)*tH1stF); // Checkeado!!! Va en el schema viejo!!! activar si se pasa al viejo
      sig21=sig21+((-1.0*ah1)*tH1stF);
//     Include missing term from contact term in 2 loop AP
      sig21=sig21+((-2*Delta2qq(Nf))*Wauxdelta); //Checkeado!!!






      sig1=sig12*LL2+sig11*LL1;
      sig2=sig24*LL4+sig23*LL3+sig22*LL2+sig21*LL1;



// tH1st // testeado!!

//sig1 = ((-2.0)*sig1); // testeado!

//sig24 = ((-2.0*LL4)*sig24); // testeado!

//sig23 = ((-2.0*LL3)*sig23); //testeado! pego en el palo, chequear de nuevo

//Wctnlofinal = ((-4.0)*(sig2+sigctgg)); // El factor a NLO es 2 a NNLO es 4 y así...
Wctnlofinal = ((-4.0)*(sig2));
}


//	for(int i=0;i<17;i++) out[i] = Wctnlofinal[i];
	for(int i=0;i<17;i++) out[i] = Wctnlofinal[i];






















}







void _hhc2ph_jet_base::amp_H1(ampq2g1p2 *amp, double *out, event_hhc2ph  p, double muf2,double mur2, unsigned int Nf,
 double z1, double z2,int order){


// double tmp;
    std::memset(out, 0, 17*sizeof(double));

//for(int i=0;i<17;i++) out[i] = 1.0;
//	return;
////////////// CALCULAR x10 y x20  ////////////////////////
	double s = p[hadron(-1)]*p[hadron(0)];
    double x10 = p[-1]*p[hadron(0)]/s;
    double x20 = p[0]*p[hadron(-1)]/s;
	int ipdf[11] = {11,10,9,7,8,0,2,1,3,4,5};





    double Q2,pT,x1p,x2p,wctnnlo;
	double xmur=sqrt(mur2);
	double xmuf=sqrt(muf2);
	Q2 = 2.0*((p[-2])*(p[-3]));
	pT = (p[1]).perp();
	double xmio = pT/sqrt(Q2);
	x1p = x10/z1;
	x2p = x20/z2;
	double LR=log(Q2/mur2);
    double LF=log(Q2/muf2);

	//Peso del contratérmino a segundo orden (NNNLO) canales qq y qg




	weight_hhc2ph Wctnlofinal,Wqq,WTL1,WTL2;






//El order siempre debe estar en 2 pues el orden 1 se calcula directamente en c++
//El peso de orden 2 (NNLO) se pasa por referencia al main_c++

	// Prepara las PDF's para cada una de las patas


	weight_hhc2ph _M_pdfL1,_M_pdfL2,_M_pdf00,_M_pdfpp;

	weight_hhc2ph retval,retvalqqL1,retvalqqL2,retvalqqL11,retvalqqL22,retval34,retval12;
	weight_hhc2ph M_CT;

	//Calcula el BORN



	 double ampB=ampBORN(p);
	 double ah1 = H1(p)/2;
//	 ah1 =0.0;
	// El factor se explica mas abajo
	// es porque en el paper 0508068 se desarrolla en potencias de alpha_s/Pi
	 double C1qqdelta;// = H1(p)/2.0/2.0;
	 C1qqdelta = 0.0;
//	 double  C1qqdeltaDY=(__DEB_PI*__DEB_PI-8.0)/3.0;
//	 double B2q_gg= (4.0/9*(__DEB_PI*__DEB_PI-3.0/4-12*Z3)+(11.0/9*__DEB_PI*__DEB_PI-193.0/12+6*Z3)+Nf/6.0*(17.0/3-4.0/9*__DEB_PI*__DEB_PI))+(11.0/3.0-Nf/3)*(-(__DEB_PI*__DEB_PI*2.0/3.0-8.0)/4.0*4.0/3.0  + (ah1/2.0/2.0 -__DEB_PI*__DEB_PI/9.0 ));

      //Schema viejo
//	 double B2qNf= -3.39493 + 2.0*C1qqdelta*(2.75 - 0.166667*Nf) + 0.421096*Nf + 0.666667*Z3;
        // Hard Scheme
     double B2qNf=(1.0/48.0)*_CFH*(2.0*Nf*(1.0 + 4.0*Z2) +    9.0* _CFH*(-1.0 + 8.0* Z2 - 16.0*Z3) +    3.0*(-17.0 - 44.0*Z2 + 72.0*Z3));

	// 28/10/2015
// Cambio en el B2q. La suma de B2q con C1qq (que incluye el C1qqdelta en el esquema correspondiente)
// debe ser independiente del proceso.
//



//	  if(order==2){
//	wctnnlo=CTNNLO(Q2,x10, x20, x1p
 //   ,x2p, xmur, xmuf, xmio, 83, ampB*__alfaQED2 , order,
  //    LL1, LL2, LL3, LL4,ah1);


 //     for(int i=0;i<17;i++) out[i]=0.0;
 //      out[0] = wctnnlo;
 //     return;

  //    }


//	if(order==1){


		double f10[13],f20[13],f1p[13],f2p[13];

	_M_pdf00=GeT_pdfCT(x10,x20,muf2, 2, 3,f10,f20);
	_M_pdfL1=GeT_pdfCT(x1p,x20,muf2, 2, 3);
	_M_pdfL2=GeT_pdfCT(x10,x2p,muf2, 2, 3);
	_M_pdfpp=GeT_pdfCT(x1p,x2p,muf2, 2, 3,f1p,f2p);

	// vuelve a la notacion antigua de primero u y despues d
	// No hace falta, se hace con los indices directamente
//	swap(f1[1],f1[2]);
//	swap(f2[1],f2[2]);
//	swap(f1[7],f1[8]);
//	swap(f2[7],f2[8]);



	//CANAL QQB
	weight_hhc2ph Waux00,WauxL1qq,WauxL2qq,WauxL1qq0,WauxL2qq0;


	// Témino de la Delta

	retval[12] = Qu4*ampB;
	retval[14] = Qd4*ampB;

	Waux00 = retval;
//Terminos sin convoluciones
	weight_hhc2ph Wauxdelta;

	Wauxdelta = Waux00*_M_pdf00;





    //-----Term Gammaqq coefficient--pata 1 (Esto significa que la
	// pata 1 emite parton) La contribucion con la PDF con xi0 se pone aparte

	double tH1stqL1;
    tH1stqL1 = Pqq(z1)*(-log(x10));

	retvalqqL1[12] =tH1stqL1*Qu4*ampB;
	//Procesos D1 Db2
	retvalqqL1[14] =tH1stqL1*Qd4*ampB;

	WauxL1qq = retvalqqL1;
	weight_hhc2ph WL1qq;
	WL1qq = WauxL1qq*_M_pdfL1;

	//-----Term Gammaqq coefficient--pata 2 (Esto significa que la
	// pata 1 emite parton) La contribucion con la PDF con xi0 se pone aparte

	double tH1stqL2;
    tH1stqL2 = Pqq(z2)*(-log(x20));

	retvalqqL2[12] =tH1stqL2*Qu4*ampB;
	//Procesos D1 Db2
	retvalqqL2[14] =tH1stqL2*Qd4*ampB;

	WauxL2qq =retvalqqL2;
	weight_hhc2ph WL2qq;
	WL2qq = WauxL2qq*_M_pdfL2;

	//PATA 1 para PDF's eveluadas en xi0


   double tH1stqL11 = -Pqq(z1)*(-log(x10))*z1-Pqqint(x10);

	retvalqqL11[12] =tH1stqL11*Qu4*ampB;
	//Procesos D1 Db2
	retvalqqL11[14] =tH1stqL11*Qd4*ampB;

	WauxL1qq0 = retvalqqL11;
	weight_hhc2ph WL1qq0;
	WL1qq0 = WauxL1qq0*_M_pdf00;

	//PATA 2 para PDF's eveluadas en xi0


   double tH1stqL22 = -Pqq(z2)*(-log(x20))*z2-Pqqint(x20);

	retvalqqL22[12] =tH1stqL22*Qu4*ampB;
	//Procesos D1 Db2
	retvalqqL22[14] =tH1stqL22*Qd4*ampB;

	WauxL2qq0 =retvalqqL22;
	weight_hhc2ph WL2qq0;
	WL2qq0 = WauxL2qq0*_M_pdf00;

	// Comenzamos a sumar las contribuciones ---Menos la Delta!!!!
	weight_hhc2ph Wqqaux;

	Wqqaux = WL1qq + WL2qq + WL1qq0 + WL2qq0;
	//El factor 2 Extra se debe a que  alpha_s/2/Pi USO!!!!!!!!!
	//Y Massi usa alpha_s/Pi
	Wqq = ((2.0*B1q)*Wauxdelta) + ((2.0)*Wqqaux)
	          +((2.0*0.5*A1q)*Wauxdelta);
//	Wqq =  realweightprod((2.0*0.5*A1q*LL2),Wauxdelta);
//	cout << LL2 << endl;

	// Fin CANAL QQBAR a Primer orden












	//CANAL QG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//Vamos con la pata 1, es decir procesos del tipo G1 * (Q2)
	//Es decir van en los canales 3 y 4
	double tH1stL1;
	//H1st: non delta terms, pata 1:

//	tH1stL1 = Cqg(z1)*(-log(x10))*0.0;


//	retval[3]=tH1stL1*Qu4*ampB;
//	retval[4]=tH1stL1*Qd4*ampB;


	//Vamos con la pata 2, es decir procesos del tipo Q1 * (G2)
	//Es decir van en los canales 1 y 2
	double tH1stL2;
	//H1st: non delta terms, pata 1:

//	tH1stL2 = Cqg(z2)*(-log(x20))*0.0;


//	retval[1]=tH1stL2*Qu4*ampB;
//	retval[2]=tH1stL2*Qd4*ampB;




	//-----Term Gammaqg coefficient--pata 1
	tH1stL1 = Pqg(z1)*(-log(x10));

	retval34[3]=tH1stL1*Qu4*ampB;
	retval34[4]=tH1stL1*Qd4*ampB;

	//-----Term Gammaqg coefficient--pata 2
	tH1stL2 = Pqg(z2)*(-log(x20));


	retval12[1]=tH1stL2*Qu4*ampB;
	retval12[2]=tH1stL2*Qd4*ampB;

	//---Pata 1 por PDF pata 1

	weight_hhc2ph WauxL1,WauxL2;

	WauxL1=retval34;
	WauxL2=retval12;

	//El factor 2 Extra se debe a que las pdfs vienen con alpha_s/2/Pi
	//Y Massi usa alpha_s/Pi
	WTL1 =( (_M_pdfL1*WauxL1));
	WTL2 =((_M_pdfL2*WauxL2));



weight_hhc2ph sig1,sig2,tH1stF;
tH1stF = (Wqqaux) + (WTL1 + WTL2);

// H1st delta term
weight_hhc2ph tH1st00,tH1stL1p0,tH1stL20p,tH1st,tH1stsust;
double tmpst;



//------------------Eschema vecchio ----------------
//------------------


//tH1st00[12] = (tmpst=2.0*C1qqdelta*ampB)*Qu4;
//tH1st00[14] = tmpst*Qd4;
//tH1st00 = tH1st00*_M_pdf00;


// 6 de Noviembre de 2015
// En el hard scheme C1qqdelta =0 pero debemos incluir el valor de
// H1qgamgam que es este ah1
// ah1 esta escrito en el hard scheme
// ah1/2 = H1notagian
// en este CTNNLO las cantidades se desarrollan en potencias de alpha_s/Pi
// por lo tanto hay que dividir por dos

tH1st00[12] = (tmpst=ah1*ampB)*Qu4;
tH1st00[14] = tmpst*Qd4;
tH1st00 = tH1st00*_M_pdf00;
tH1stsust = tH1st00;


//H1st non delta terms
// NO entra a orden 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//Lo dejo escrito para más adelante
//Pata 1!!!!
tH1stL1p0[12]=(tmpst=Cqq(z1)*(-log(x10))*ampB)*Qu4;
tH1stL1p0[14]=tmpst*Qd4;

tH1stL1p0[3]=(tmpst=Cqg(z1)*(-log(x10))*ampB)*Qu4;
tH1stL1p0[4]=tmpst*Qd4;
tH1stL1p0 = tH1stL1p0*_M_pdfL1;
    // Pata 2!!!
tH1stL20p[12]=(tmpst=Cqq(z2)*(-log(x20))*ampB)*Qu4;
tH1stL20p[14]=tmpst*Qd4;

tH1stL20p[1]=(tmpst=Cqg(z2)*(-log(x20))*ampB)*Qu4;
tH1stL20p[2]=tmpst*Qd4;
tH1stL20p = tH1stL20p*_M_pdfL2;

tH1st = tH1stL20p + tH1stL1p0 + tH1st00;

 sig1 = (tH1st+LF*tH1stF);

 Wctnlofinal=2.0*sig1;

// Parte a NNLO

 if(order==2){



// Now (gamma+gamma)*(gamma+gamma) term: to be used later
double diffg1fp,diffg1f0,diffg1gp,diffg2fp,diffg2f0,diffg2gp;
weight_hhc2ph tgagapp,tgaga00,tgagap0,tgaga0p,tgaga;

// First part: one gamma for each leg:

diffg1fp = (-log(x10))*Pqq(z1);
diffg1gp = (-log(x10))*Pqg(z1);
diffg1f0 = (-log(x10))*(-z1*Pqq(z1))-Pqqint(x10);

diffg2fp = (-log(x20))*Pqq(z2);
diffg2gp = (-log(x20))*Pqg(z2);
diffg2f0 = (-log(x20))*(-z2*Pqq(z2))-Pqqint(x20);



 // --------  H2 contribution
 weight_hhc2ph tH2stpp,tH2stp0,tH2st0p,tH2st;
 // ------- H2st gg contribution
 tH2stpp[0] = Cqg(z1)*(-log(x10))*(-log(x20))*Cqg(z2)*ampB*Qsum4*2.0;
 // ------- H2stqqbar from C1*C1 (without delta term)
 tH2stpp[12] = (tmpst=Cqq(z1)*(-log(x10))*(-log(x20))*Cqq(z2)*ampB)*Qu4;
 tH2stpp[14] = tmpst*Qd4;

 //---------- regular-delta
 tH2stp0[12] =  (tmpst=Cqq(z1)*(-log(x10))*(C1qqdelta)*ampB)*Qu4;
 tH2stp0[14] =  tmpst*Qd4;
 tH2st0p[12] =  (tmpst=Cqq(z2)*(-log(x20))*(C1qqdelta)*ampB)*Qu4;
 tH2st0p[14] =  tmpst*Qd4;

 // ------------ H2st qg contribution from C1*C1
 tH2stpp[3] = (tmpst=Cqg(z1)*Cqq(z2)*(-log(x10))*(-log(x20))*ampB)*Qu4;
 tH2stpp[4] = tmpst*Qd4;
 tH2stpp[1] = (tmpst=Cqg(z2)*Cqq(z1)*(-log(x10))*(-log(x20))*ampB)*Qu4;
 tH2stpp[2] = tmpst*Qd4;
 // ------------ regular-delta
 tH2stp0[3] = (tmpst=Cqg(z1)*C1qqdelta*(-log(x10))*ampB)*Qu4;
 tH2stp0[4] = tmpst*Qd4;
 tH2st0p[1] = (tmpst=Cqg(z2)*C1qqdelta*(-log(x20))*ampB)*Qu4;
 tH2st0p[2] = tmpst*Qd4;


 // ------------ C2qq regular part, first leg
 tH2stp0[12] += (tmpst=C2QQREG(z1)*(-log(x10))*ampB)*Qu4;
 tH2stp0[14] += tmpst*Qd4;
 // ------------ C2qq regular part, second leg
 tH2st0p[12] += (tmpst=C2QQREG(z2)*(-log(x20))*ampB)*Qu4;
 tH2st0p[14] += tmpst*Qd4;


 // ------------ C2qg first leg
 tH2stp0[3] += (tmpst=C2QG(z1)*(-log(x10))*ampB)*Qu4;
  tH2stp0[4] += tmpst*Qd4;
 // ------------ C2qg second leg
 tH2st0p[1] += (tmpst=C2QG(z2)*(-log(x20))*ampB)*Qu4;
 tH2st0p[2] += tmpst*Qd4;
 // ------------ C2qqbar contribution first leg
 tH2stp0[12] += (tmpst=C2QQP(z1)*(-log(x10))*ampB)*Qu4;
 tH2stp0[14] += tmpst*Qd4;
 // ------------ C2qqbar contribution second leg
 tH2st0p[12] += (tmpst=C2QQP(z2)*(-log(x20))*ampB)*Qu4;
 tH2st0p[14] += tmpst*Qd4;



 tH2stpp = tH2stpp*_M_pdfpp;
 tH2st0p = tH2st0p*_M_pdfL2;
 tH2stp0 = tH2stp0*_M_pdfL1;

 tH2st = tH2stpp + tH2st0p + tH2stp0;

// End H2st

tgagapp[0] = diffg1gp*diffg2gp*ampB*Qsum4*2.0; //Checkeado consistencia, falta el vero!!!!
tgagapp[1] = (tmpst=diffg1fp*diffg2gp*ampB)*Qu4;
tgagapp[2] = tmpst*Qd4;
tgagapp[3] = (tmpst=diffg1gp*diffg2fp*ampB)*Qu4;
tgagapp[4] = tmpst*Qd4;
tgagapp[12] = (tmpst=diffg1fp*diffg2fp*ampB)*Qu4;
tgagapp[14] = tmpst*Qd4;
tgagapp = tgagapp*_M_pdfpp;

// Terminos 0p
tgaga0p[1]= (tmpst=diffg1f0*diffg2gp*ampB)*Qu4;
tgaga0p[2]= tmpst*Qd4;
tgaga0p[12]=(tmpst=diffg1f0*diffg2fp*ampB)*Qu4;
tgaga0p[14]=tmpst*Qd4;
tgaga0p = tgaga0p*_M_pdfL2;

// Terminos p0
tgagap0[3]= (tmpst=diffg2f0*diffg1gp*ampB)*Qu4;
tgagap0[4]= tmpst*Qd4;
tgagap0[12]=(tmpst=diffg1fp*diffg2f0*ampB)*Qu4;
tgagap0[14]=tmpst*Qd4;
tgagap0 = tgagap0*_M_pdfL1;

// Terminos 00
tgaga00[12]=(tmpst=diffg1f0*diffg2f0*ampB)*Qu4;
tgaga00[14]=tmpst*Qd4;
tgaga00 = tgaga00*_M_pdf00;


// Radiogaga fuera del loop Checkeado!!!
tgaga = 2.0*(tgaga00 + tgaga0p +tgagap0 + tgagapp);

//     Second part: gamma*gamma terms
//c     Pij * Pjk = D1ijjk (log(1-z)/(1-z))_+ + D0ijjk/(1-z)_+
//c              + Pijjk(z) + Deltaijjk delta(1-z)
//C     First leg
double diff100,diff1p0,diff1p0g,diff200,diff20p,diff20pg;
weight_hhc2ph diff1loop00,diff2loop00,diff1loopp0,diff2loop0p;

diff100=(-log(x10))*(-z1*(D0qqqq/(1-z1) + D1qqqq*log(1-z1)/(1-z1)))+(Deltaqqqq -D0qqqq*D0int(x10) - D1qqqq*D1int(x10));
diff1loop00[12] = (tmpst=diff100*ampB)*Qu4;
diff1loop00[14] = tmpst*Qd4;
diff1loop00 = diff1loop00*_M_pdf00;

diff1p0 = (-log(x10))*((D0qqqq/(1-z1) + D1qqqq*log(1-z1)/(1-z1)) + Pqqqq(z1));
diff1p0g = (-log(x10))* (Pqqqg(z1)+Pqggg(z1,Nf));

diff1loopp0[3] = (tmpst=diff1p0g*ampB)*Qu4;
diff1loopp0[4] = tmpst*Qd4;
diff1loopp0[12] = (tmpst=diff1p0*ampB)*Qu4;
diff1loopp0[14] = tmpst*Qd4;
diff1loopp0 = diff1loopp0*_M_pdfL1;

//Second leg
diff200=(-log(x20))*(-z2*(D0qqqq/(1-z2) + D1qqqq*log(1-z2)/(1-z2)))+(Deltaqqqq -D0qqqq*D0int(x20) - D1qqqq*D1int(x20));
diff2loop00[12] = (tmpst=diff200*ampB)*Qu4;
diff2loop00[14] = (tmpst)*Qd4;
diff2loop00 = diff2loop00*_M_pdf00;

diff20p = (-log(x20))*((D0qqqq/(1-z2) + D1qqqq*log(1-z2)/(1-z2)) + Pqqqq(z2));
diff20pg = (-log(x20))* (Pqqqg(z2)+Pqggg(z2,Nf));

diff2loop0p[1] = (tmpst=diff20pg*ampB)*Qu4;
diff2loop0p[2] = tmpst*Qd4;
diff2loop0p[12] = (tmpst=diff20p*ampB)*Qu4;
diff2loop0p[14] = tmpst*Qd4;
diff2loop0p = diff2loop0p*_M_pdfL2;



//     Include Pqggq
tgaga = tgaga + diff2loop0p + diff2loop00 + diff1loopp0 + diff1loop00;

// Tambien esta parte de Radio gaga checkeada!!!!

//    End of (gamma+gamma)*(gamma+gamma) term: FLGQ non trivial here ! DONE

//    Start  (C+C)*(gamma+gamma) term

//    gamma first leg, C second leg
double diffc2f0,diffc2fp,diffc20gp,diffc1f0,diffc1fp,diffc10gp;
weight_hhc2ph tcgapp,tcga00,tcgap0,tcga0p,tcga;

// First part: one gamma for each leg:

diffc2fp = (-log(x20))*Cqq(z2);
diffc2f0 = C1qqdelta;
diffc20gp = (-log(x20))*Cqg(z2);

diffc1fp = (-log(x10))*Cqq(z1);
diffc1f0 = C1qqdelta;
diffc10gp = (-log(x10))*Cqg(z1);

// Terminos pp
tcgapp[0] = diffg1gp*diffc20gp*ampB*Qsum4*2.0;
tcgapp[1] = (tmpst=diffg1fp*diffc20gp*ampB)*Qu4;
tcgapp[2] = tmpst*Qd4;
tcgapp[3] = (tmpst=diffg1gp*diffc2fp*ampB)*Qu4;
tcgapp[4] = tmpst*Qd4;
tcgapp[12] =(tmpst=diffg1fp*diffc2fp*ampB)*Qu4;
tcgapp[14] =tmpst*Qd4;
//tcgapp = tcgapp*_M_pdfpp;

// Terminos 0p
tcga0p[1]= (tmpst=diffg1f0*diffc20gp*ampB)*Qu4;
tcga0p[2]= tmpst*Qd4;
tcga0p[12]=(tmpst=diffg1f0*diffc2fp*ampB)*Qu4;
tcga0p[14]=tmpst*Qd4;
//tcga0p = tcga0p*_M_pdfL2;

// Terminos p0
tcgap0[3]= (tmpst=diffc2f0*diffg1gp*ampB)*Qu4;
tcgap0[4]= tmpst*Qd4;
tcgap0[12]=(tmpst=diffg1fp*diffc2f0*ampB)*Qu4;
tcgap0[14]=tmpst*Qd4;
//tcgap0 = tcgap0*_M_pdfL1;

// Terminos 00
tcga00[12]=(tmpst=diffg1f0*diffc2f0*ampB)*Qu4;
tcga00[14]=tmpst*Qd4;
//tcga00 = tcga00*_M_pdf00;
//tcga = tcga00 + tcga0p +tcgap0 + tcgapp;

// First leg, gamma second leg
/*
diffc1fp = (-log(x10))*Cqq(z1);
diffc1f0 = C1qqdelta;
diffc10gp = (-log(x10))*Cqg(z1);
*/
// Terminos pp
tcgapp[0] += diffc10gp*diffg2gp*ampB*Qsum4*2.0;
tcgapp[1] += (tmpst=diffc1fp*diffg2gp*ampB)*Qu4;
tcgapp[2] += tmpst*Qd4;
tcgapp[3] += (tmpst=diffc10gp*diffg2fp*ampB)*Qu4;
tcgapp[4] += tmpst*Qd4;
tcgapp[12] +=(tmpst=diffc1fp*diffg2fp*ampB)*Qu4;
tcgapp[14] +=tmpst*Qd4;
tcgapp = tcgapp*_M_pdfpp;

// Terminos 0p
tcga0p[1]+= (tmpst=diffc1f0*diffg2gp*ampB)*Qu4;
tcga0p[2]+= tmpst*Qd4;
tcga0p[12]+=(tmpst=diffc1f0*diffg2fp*ampB)*Qu4;
tcga0p[14]+=tmpst*Qd4;
tcga0p = tcga0p*_M_pdfL2;

// Terminos p0
tcgap0[3]+= (tmpst=diffg2f0*diffc10gp*ampB)*Qu4;
tcgap0[4]+= tmpst*Qd4;
tcgap0[12]+=(tmpst=diffc1fp*diffg2f0*ampB)*Qu4;
tcgap0[14]+=tmpst*Qd4;
tcgap0 = tcgap0*_M_pdfL1;

// Terminos 00
tcga00[12]+=(tmpst=diffc1f0*diffg2f0*ampB)*Qu4;
tcga00[14]+=tmpst*Qd4;
tcga00 = tcga00*_M_pdf00;
tcga = tcga00 + tcga0p +tcgap0 + tcgapp;


//   C*gamma: first leg (ignore delta term in Cqq: taken into account with tH1stF)

weight_hhc2ph tcgamulp0,tcgamul0p;


tcgamulp0[3]= (tmpst=(-log(x10))*(CqqPqg(z1)+CqgPgg(z1,Nf))*ampB)*Qu4;
tcgamulp0[4]= tmpst*Qd4;
tcgamulp0[12]= (tmpst=(-log(x10))*CqqPqq(z1)*ampB)*Qu4;
tcgamulp0[14]= tmpst*Qd4;
tcgamulp0 = tcgamulp0*_M_pdfL1;

// c    C*gamma: second leg (ignore delta term in Cqq: taken into account with tH1stF)

tcgamul0p[1]= (tmpst=(-log(x20))*(CqqPqg(z2)+CqgPgg(z2,Nf))*ampB)*Qu4;
tcgamul0p[2]= tmpst*Qd4;
tcgamul0p[12]= (tmpst=(-log(x20))*CqqPqq(z2)*ampB)*Qu4;
tcgamul0p[14]= tmpst*Qd4;
tcgamul0p = tcgamul0p*_M_pdfL2;


tcga += tcgamulp0 + tcgamul0p;


//double Kappa= 67.0/6.0-(__DEB_PI*__DEB_PI)/2.0-5.0/9.0*5.0;

// tcga chequeado!!!!

//  Start 2-loop AP
weight_hhc2ph tgamma2p0,tgamma20p,tgamma2,tgamma00;

tgamma2p0[3] = (tmpst=(-log(x10))*P2qg(z1)*ampB)*Qu4;
tgamma2p0[4] = tmpst*Qd4;

tgamma2p0[12] = (tmpst=(-log(x10))*P2qqV(z1,Nf)*ampB)*Qu4;
tgamma2p0[14] = tmpst*Qd4;
tgamma2p0[12] += (tmpst=(-log(x10))/(1-z1)*2.0/3.0*Kappa(Nf)*ampB)*Qu4;
tgamma2p0[14] += tmpst*Qd4;

tgamma2p0 = tgamma2p0*_M_pdfL1;

tgamma20p[1] =(tmpst=(-log(x20))*P2qg(z2)*ampB)*Qu4;
tgamma20p[2] =tmpst*Qd4;

tgamma20p[12] =(tmpst=(-log(x20))*P2qqV(z2,Nf)*ampB)*Qu4;
tgamma20p[14] =tmpst*Qd4;
tgamma20p[12] += (tmpst=(-log(x20))/(1-z2)*2.0/3.0*Kappa(Nf)*ampB)*Qu4;
tgamma20p[14] += tmpst*Qd4;

tgamma20p = tgamma20p*_M_pdfL2;

tgamma00[12] = (tmpst=((-log(x10))*(-z1)/(1-z1)-D0int(x10))*2.0/3.0*Kappa(Nf)*ampB)*Qu4;
tgamma00[14] = tmpst*Qd4;

tgamma00[12] += (tmpst=((-log(x20))*(-z2)/(1-z2)-D0int(x20))*2.0/3.0*Kappa(Nf)*ampB)*Qu4;
tgamma00[14] += tmpst*Qd4;
tgamma00 = tgamma00*_M_pdf00;

//cout << Kappa(Nf) << endl;

tgamma2 = tgamma20p + tgamma2p0 +tgamma00;


// tgamma2 Chequeado!!! mucha dispersion, en ambos codigos tarda en converger y los errores
// estan subestimados

// Ahora metemos los canales qq,qq',etc que van a ir en el canal 10 al final
// Todos estos términos están calculados dentro del born

 /* static const int jdo[11] = {-5,-4,-3,-2,-1,0,1,2,3,4,5};
 static const int kdo[11] = {-5,-4,-3,-2,-1,0,1,2,3,4,5};
*/

double carga=0.0;
double diffH2st=0.0,tH2std=0.0;
double mdiff1=0.0,mdiff2=0.0,mtgaga=0.0,mtcga=0.0,mtgamma2=0.0;
double eq2[6] = {0.0,4.0/9.0,1.0/9.0,1.0/9.0,4.0/9.0,1.0/9.0};
for(int jj=-5;jj<6;jj++){
	for(int kk=-5;kk<6;kk++){

	if(Msqc( jj,  kk)!=0.0){

//--------------- Begining of H2st

// ----- H2st qqbar channel:  D0(z), first leg

	diffH2st =  (-log(x10))*(f1p[ipdf[jj+5]]-f10[ipdf[jj+5]]*z1)*H2qqD0(Nf)/(1-z1);
	tH2std += 0.5*diffH2st*f20[ipdf[kk+5]]*Msqc( jj,  kk)*ampB;
	tH2std += -0.5*H2qqD0(Nf)*D0int(x10)*f10[ipdf[jj+5]]*f20[ipdf[kk+5]]*Msqc( jj,  kk)*ampB;
// ----- H2st, qqbar channel: D0(z), second leg
    diffH2st =  (-log(x20))*(f2p[ipdf[kk+5]]-f20[ipdf[kk+5]]*z2)*H2qqD0(Nf)/(1-z2);
    tH2std += 0.5*diffH2st*f10[ipdf[jj+5]]*Msqc( jj,  kk)*ampB;
	tH2std += -0.5*H2qqD0(Nf)*D0int(x20)*f20[ipdf[kk+5]]*f10[ipdf[jj+5]]*Msqc( jj,  kk)*ampB;



	for(int aH2=1;aH2<6;aH2++){
		// Cqqp contribution: first leg
		if(aH2 != abs(jj)){tH2std += (f1p[ipdf[aH2+5]]+f1p[ipdf[-aH2+5]])*C2QQP(z1)*(-log(x10))*f20[ipdf[kk+5]]*Msqc( jj,  kk)*ampB; }
		if(aH2 != abs(kk)){tH2std += (f2p[ipdf[aH2+5]]+f2p[ipdf[-aH2+5]])*C2QQP(z2)*(-log(x20))*f10[ipdf[jj+5]]*Msqc( jj,  kk)*ampB; }
}



	mdiff1=mdiff2=0.0;

	// Include Pqggq !Diff contiene todo lo anterior tambien!!!
	for(int l=1;l<6;l++){
	mdiff1 += - log(x10)*(f1p[ipdf[l+5]]+f1p[ipdf[-l+5]])*Pqggq(z1);
	mdiff2 += - log(x20)*(f2p[ipdf[l+5]]+f2p[ipdf[-l+5]])*Pqggq(z2);
}
    mtgaga += mdiff1*f20[ipdf[kk+5]]*Msqc( jj,  kk)*ampB;
    mtgaga += mdiff2*f10[ipdf[jj+5]]*Msqc( jj,  kk)*ampB;


	for(int ll=1;ll<6;ll++){
		mtcga=mtcga+(f1p[ipdf[ll+5]]+f1p[ipdf[-ll+5]])*CqgPgq(z1)*(-log(x10))*f20[ipdf[kk+5]]*Msqc( jj,  kk)*ampB;
	    mtcga=mtcga+(f2p[ipdf[ll+5]]+f2p[ipdf[-ll+5]])*CqgPgq(z2)*(-log(x20))*f10[ipdf[jj+5]]*Msqc( jj,  kk)*ampB;
}


   for(int lll=-5;lll<6;lll++){
  if(lll!=0){

  mtgamma2 = mtgamma2 + f1p[ipdf[lll+5]]*P2qqS(z1)*(-log(x10))*f20[ipdf[kk+5]]*Msqc( jj,  kk)*ampB;
  mtgamma2 = mtgamma2 + f2p[ipdf[lll+5]]*P2qqS(z2)*(-log(x20))*f10[ipdf[jj+5]]*Msqc( jj,  kk)*ampB;

}
}

mtgamma2 = mtgamma2 + f1p[ipdf[-jj+5]]*P2qqbV(z1)*(-log(x10))*f20[ipdf[kk+5]]*Msqc( jj,  kk)*ampB;
mtgamma2 = mtgamma2 + f2p[ipdf[-kk+5]]*P2qqbV(z2)*(-log(x20))*f10[ipdf[jj+5]]*Msqc( jj,  kk)*ampB;
}
}
}

//cout << "carga   :" << carga << endl;

tgaga[10] = mtgaga;
tcga[10] = mtcga;
tgamma2[10] = mtgamma2;
tH2st[10] = tH2std;
weight_hhc2ph pepa;
pepa[10] = tH2std;







// ----------------------- From H2.f (old) ------------------------- //
       sig2 = (tH2st);
      // sig2 = pepa;

       sig2+=(0.5*beta0(Nf)*LF*LF*tH1stF+tgamma2*LF-beta0(Nf)*LR*(tH1st+LF*tH1stF)+LF*tcga+0.5*LF*LF*tgaga);

       sig2+=(2.0*Delta2qq(Nf)*(Wauxdelta))*LF;

//     Include missing terms from Hard Scheme -tH1stsust
       sig2+=(ah1*(tH1st-tH1stsust));
       sig2+=(ah1*LF*tH1stF);


//Wctnlofinal = ((-4.0)*(sig2+sigctgg)); // El factor a NLO es 2 a NNLO es 4 y así...
Wctnlofinal = ((4.0)*(sig2));
}


//	for(int i=0;i<17;i++) out[i] = Wctnlofinal[i];
	for(int i=0;i<17;i++) out[i] = Wctnlofinal[i];












//El order siempre debe estar en 2 pues el orden 1 se calcula directamente en c++
//El peso de orden 2 (NNLO) se pasa por referencia al main_c++
//Lo que se devuelve aqui abajo son los canales qq y qg a orden alpha_s












}

void _hhc2ph_jet_base::amp_H2f(ampq2g1p2 *amp, double *out,double mur2,unsigned int Nf)
  {
    double tmp;
    std::memset(out, 0, 17*sizeof(double));

    //  gg ->

    //  ug -> u    (1/Nc * 1/Na)
    //  dg -> d    (1/Nc * 1/Na)


    //  gu -> u    (1/Na * 1/Nc)
    //  gd -> d    (1/Na * 1/Nc)


    //  uU ->
    //  uu ->
    //  dD ->
    //  dd ->
    //  ud ->
    //  du ->
    //  uUb ->

    //  uub -> g   (1/Nc * 1/Nc)
    out[12] = Qu4*((amp -> H2f(-1,0, 1, -2,-3,mur2,Nf,sqrt(Qu4))));


      //  dDb ->
    //  ddb -> g   (1/Nc * 1/Nc)
    out[14] = Qd4*((amp -> H2f(-1,0, 1, -2,-3,mur2,Nf,sqrt(Qd4))));

  //for(int i=0;i<17;i++) out[i]=1.0;

    //  udb ->
    //  dub ->
  }



void _hhc2ph_jet_base::amp_BOX(ampq2g1p2 *amp, double *out)
  {
 //   double tmp;
    std::memset(out, 0, 17*sizeof(double));

    //  gg ->
     out[0] = (amp -> BOXMTR(0,-1, 1, -2,-3));
    //  ug -> u    (1/Nc * 1/Na)
    //  dg -> d    (1/Nc * 1/Na)


    //  gu -> u    (1/Na * 1/Nc)
    //  gd -> d    (1/Na * 1/Nc)


    //  uU ->
    //  uu ->
    //  dD ->
    //  dd ->
    //  ud ->
    //  du ->
    //  uUb ->

    //  uub -> g   (1/Nc * 1/Nc)


    //  dDb ->
    //  ddb -> g   (1/Nc * 1/Nc)


  //for(int i=0;i<17;i++) out[i]=1.0;

    //  udb ->
    //  dub ->
  }
  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
// -------------- Interferencia LeanR ---------------------------------------
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  

  

  void _hhc2ph_jet_base::amp_ggffg(ampq2g1p2 *amp, double *out)
  {
 //   double tmp;
    std::memset(out, 0, 17*sizeof(double));

 // std::cout << " desde hhc2phjet.cc  amp_ggffg " << std::endl;

    //  gg ->
    //out[0] = (amp -> wgtR(0,-1, 1, -2,-3));
     //double factors = 3.0*Lzm*Lz*__DEB_PI/16.0*pow((3.0*(1.0/9.0)+2.0*(4.0/9.0)),2);
     
     out[0] =(amp -> wgtR(0,-1, 1, -2,-3));
     //  out[0] =1.0*Lzm*Lz;
    // No contributions from other channels

  }
  void _hhc2ph_jet_base::amp_ggHg(ampq2g1p2 *amp, double *out)
  {
 //   double tmp;
    std::memset(out, 0, 17*sizeof(double));

   
    //    std::cout << " desde hhc2phjet.cc  amp_ggHg " << std::endl;

    //  gg ->
     out[0] = (amp -> wgtHsig(0,-1, 1, -2,-3));
    // No contributions from other channels

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

 void _hhc2ph_jet_base::amp_MggHyy(ampq2g1p2 *amp, double *out)
  {
   //   double tmp;
    std::memset(out, 0, 17*sizeof(double));
    double s12=0.0;
   // debemos usar solo : LO_Mggyy, LO_MHyy, LO_MggH
   // falta dividir las contribuciones en Señal, interferencia real e imaginaria
   // Contrib_Inter = "SIGNAL_LO"; // Contrib_Inter = SIGNAL_LO, Int_LO, SigInt   
    double Aggyyp  = (amp -> LO_Mggyy(-1,0, 1, -2,-3,s12));
    std::complex<double> Aggh   = ((amp -> LO_MggH(-1,0, 1, -2,-3)));
    std::complex<double> Aggyy = std::complex<double> (Aggyyp,0.0);
  //  std::complex<double> Aggh  = std::complex<double> (Agghp ,0.0);
    std::complex<double> Ahyy   = ((amp -> LO_MHyy(-1,0, 1, -2,-3)));
    std::complex<double> N_int  = Aggh*Ahyy*conj(Aggyy);
    double N_int_Re = real(N_int);
    double N_Int_Im = imag(N_int); //Questo meno viene messo nel codice Beta4 pero no nel paper di Martin?!!!!!!!!!!!!!!!!!
                                      //Pregunta!!!!!!!!!!!!! Chicos revisar!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                      // Si hay un typo en el paper de Martin!!!!!!!!!! 
    double AtotR=0.0;
      
    /*
    cout << " Aggyyp : " << Aggyyp << endl;
    cout << " Aggh : " <<Aggh <<endl;
    cout << " Aggyy : " << Aggyy << endl;
    cout << " Ahyy  : " << Ahyy  << endl;
    cout << " N_int : " << N_int << endl;
    cout << " N_int_Re : " << N_int_Re <<endl;
    cout << " N_Int_Im : " << N_Int_Im <<endl;
    */


     double mh=_MH_;
     double Gammah = _GammaH_; // Ok
     double s=s12; // Ok
    // s=(124.95)*(124.95);
     
     
          if (Contrib_Inter == "SIGNAL_LO"){
        AtotR=abs((Ahyy*Aggh)*conj(Ahyy*Aggh) );
      //std::cout << AtotR << std::endl;
	  } else if (Contrib_Inter == "Int_LO") 
      { AtotR= - 2.0*(s-mh*mh)*N_int_Re - 2.0*mh*Gammah*N_Int_Im;
		  // Real Part checked!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		  // AtotR=-2.0*(s-mh*mh)*N_int_Re ;
		  //   std::cout << AtotR << std::endl;
      } else if(Contrib_Inter == "SigInt") {
        AtotR=abs((Ahyy*Aggh)*conj(Ahyy*Aggh))- 2.0*(s-mh*mh)*N_int_Re - 2.0*mh*Gammah*N_Int_Im;
      //  AtotR=-2.0*(s-mh*mh)*N_int_Re ;
      //std::cout << Ahyy << std::endl;
      }
	  //	  cout << " AtotR : " << AtotR << endl;
//----------------------------------------------------------------------
//----------------------------------------------------------------------      
      // -- Calculo normal mio --- comparable al codigo Beta4
      /*      
      double xnorm = 1.0/(32.0* __DEB_PI   * __DEB_PI   );
      double Ds = 1.0/((-mh*mh + s)*(-mh*mh + s) + mh*mh*Gammah*Gammah);
      double M = Ds*AtotR/(2.0*s); 
      double sig0=M*xnorm/16.0; // extra 4 Pi factor squared  2__DEB_2PI  due to alpha_s expansion as alpha_s/(2 __DEB_2PI  ) Incluido en nlo-process_i2f0.h
     

//-- FactorM is included in the phase space


      */
//----------------------------------------------------------------------      
//----------------------------------------------------------------------      
      
//	  AtotR =  0.802746;
       double sig0=AtotR;
     
 //    cout << "sig0 from original:  "<< sig0 << endl;
      
/*
std::complex<double> a = std::complex<double> (3.0 , 7.0);
std::complex<double> b = std::complex<double> (5.0 , 2.0);

std::cout << "a*b   :" << a*b << std::endl;
std::cout << "a*conj(a)   :" << a*conj(a) << std::endl;
std::cout << "abs(a*conj(a))   :" << abs(a*conj(a)) << std::endl;
*/

//std::cout << "Contrib_Inter  :" << Contrib_Inter << std::endl;

    //  gg ->
     out[0] = sig0;
    // No contributions from other channels

  }

//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  




//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  
//--------------------------------------------------------------------------  


void _hhc2ph_jet_base::amp_H1f(ampq2g1p2 *amp, double *out)
  {
    double tmp;
    std::memset(out, 0, 17*sizeof(double));

    //  gg ->

    //  ug -> u    (1/Nc * 1/Na)
    //  dg -> d    (1/Nc * 1/Na)


    //  gu -> u    (1/Na * 1/Nc)
    //  gd -> d    (1/Na * 1/Nc)


    //  uU ->
    //  uu ->
    //  dD ->
    //  dd ->
    //  ud ->
    //  du ->
    //  uUb ->

    //  uub -> g   (1/Nc * 1/Nc)
    out[12] =Qu4* (tmp = (( (amp -> LO_H1f(-1,0, 1, -2,-3)))  )    );

    //  dDb ->
    //  ddb -> g   (1/Nc * 1/Nc)
    out[14] = Qd4*tmp;

  //for(int i=0;i<17;i++) out[i]=1.0;

    //  udb ->
    //  dub ->
  }


///////////////////////////////////////////////////////////////////////////



  void _hhc2ph_jet_base::amp_tree(ampq2g1p2 *amp, double *out)
  {
    double tmp;
    std::memset(out, 0, 17*sizeof(double));

    //  gg ->

    //  ug -> u    (1/Nc * 1/Na)
    //  dg -> d    (1/Nc * 1/Na)
    out[1] = Qu4*(tmp = -(amp -> su3_tree(-1,1, 0, -2,-3))/(Nc*Na));
    out[2] = Qd4*tmp;

    //  gu -> u    (1/Na * 1/Nc)
    //  gd -> d    (1/Na * 1/Nc)
    out[3] = Qu4*(tmp = -(amp -> su3_tree(0,1, -1, -2,-3))/(Nc*Na));
    out[4] = Qd4*tmp;

    //  uU ->
    //  uu ->
    //  dD ->
    //  dd ->
    //  ud ->
    //  du ->
    //  uUb ->

    //  uub -> g   (1/Nc * 1/Nc)
    out[12] = Qu4*(tmp = (amp -> su3_tree(-1,0, 1, -2,-3))/Nc2);

    //  dDb ->
    //  ddb -> g   (1/Nc * 1/Nc)
    out[14] = Qd4*tmp;

  //for(int i=0;i<17;i++) out[i]=1.0;

    //  udb ->
    //  dub ->
  }

  std::pair<double, std::complex<double> > _hhc2ph_jet_base::
  amp_ccag(ampq2g1p2 *amp, int i, int j)
  {
    _Pair cc(amp -> su3_cc(i,j, 1, -1, 0, -2,-3));
    cc.first /= -Na*Nc;
    cc.second /= -Na*Nc;
    return cc;
  }

  std::pair<double, std::complex<double> > _hhc2ph_jet_base::
  amp_ccga(ampq2g1p2 *amp, int i, int j)
  {
    _Pair cc(amp -> su3_cc(i,j, 1, 0, -1, -2,-3));
    cc.first /= -Na*Nc;
    cc.second /= -Na*Nc;
    return cc;
  }

  std::pair<double, std::complex<double> > _hhc2ph_jet_base::
  amp_ccqa(ampq2g1p2 *amp, int i, int j)
  {
    _Pair cc(amp -> su3_cc(i,j, -1, 0, 1, -2,-3));
    cc.first /= Nc2;
    cc.second /= Nc2;
    return cc;
  }

  void _hhc2ph_jet_base::amp_kp(double al, ampq2g1p2 *amp, su3_kp_i2 *out)
  {
    static su3_kp_i2 tmp;

    //  ug -> u    (1/Nc * 1/Na)
    //  dg -> d    (1/Nc * 1/Na)
    amp -> su3_kp(Nf, -1,0, -1,1, 0, -2,-3, &tmp, al);
    out[0] = (-Qu4/(Nc*Na))*tmp;
    out[1] = (-Qd4/(Nc*Na))*tmp;

    //  gu -> u    (1/Na * 1/Nc)
    //  gd -> d    (1/Na * 1/Nc)
    amp -> su3_kp(Nf, -1,0, 0,1, -1, -2,-3, &tmp, al);
    out[2] = (-Qu4/(Nc*Na))*tmp;
    out[3] = (-Qd4/(Nc*Na))*tmp;

    //  uub -> g   (1/Nc * 1/Nc)
    //  ddb -> g   (1/Nc * 1/Nc)
    amp -> su3_kp(Nf, -1,0, -1,0, 1, -2,-3, &tmp, al);
    out[4] = (Qu4/Nc2)*tmp;
    out[5] = (Qd4/Nc2)*tmp;
  }

  void _hhc2ph_jet_base::amp_1loop(ampq2g1p2 *amp, double *out)
  {
    //  ug -> u    (1/Nc * 1/Na)
    //  dg -> d    (1/Nc * 1/Na)
    amp -> su3_1loop(Nu, Nd, -1,1, 0, -2,-3, out);
    out[0] /= Nc*Na;
    out[1] /= Nc*Na;

    //  gu -> u    (1/Na * 1/Nc)
    //  gd -> d    (1/Na * 1/Nc)
    amp -> su3_1loop(Nu, Nd, 0,1, -1, -2,-3, out+2);
    out[2] /= Nc*Na;
    out[3] /= Nc*Na;

    //  uub -> g   (1/Nc * 1/Nc)
    //  ddb -> g   (1/Nc * 1/Nc)
    amp -> su3_1loop(Nu, Nd, -1,0, 1, -2,-3, out+4);
    out[4] /= Nc2;
    out[5] /= Nc2;
  }

  void _hhc2ph_jet_base::amp_1loop_mch(ampq2g1p2 *amp, double *out)
  {
    //  ug -> u    (1/Nc * 1/Na)
    //  dg -> d    (1/Nc * 1/Na)
    amp -> su3_1loop_mch(Nu, Nd, -1,1, 0, -2,-3, out);
    out[0] /= Nc*Na;
    out[1] /= Nc*Na;

    //  gu -> u    (1/Na * 1/Nc)
    //  gd -> d    (1/Na * 1/Nc)
    amp -> su3_1loop_mch(Nu, Nd, 0,1, -1, -2,-3, out+2);
    out[2] /= Nc*Na;
    out[3] /= Nc*Na;

    //  uub -> g   (1/Nc * 1/Nc)
    //  ddb -> g   (1/Nc * 1/Nc)
    amp -> su3_1loop_mch(Nu, Nd, -1,0, 1, -2,-3, out+4);
    out[4] /= Nc2;
    out[5] /= Nc2;
  }


void _hhc2ph_jet_base::amp_1loop_mch(ampq2g1p2 *amp, double *out,double rng1,double rng2)
  {
    //  ug -> u    (1/Nc * 1/Na)
    //  dg -> d    (1/Nc * 1/Na)
    amp -> su3_1loop_mch(Nu, Nd, -1,1, 0, -2,-3, out,rng1,rng2);
    out[0] /= Nc*Na;
    out[1] /= Nc*Na;

    //  gu -> u    (1/Na * 1/Nc)
    //  gd -> d    (1/Na * 1/Nc)
    amp -> su3_1loop_mch(Nu, Nd, 0,1, -1, -2,-3, out+2,rng1,rng2);
    out[2] /= Nc*Na;
    out[3] /= Nc*Na;

    //  uub -> g   (1/Nc * 1/Nc)
    //  ddb -> g   (1/Nc * 1/Nc)
    amp -> su3_1loop_mch(Nu, Nd, -1,0, 1, -2,-3, out+4,rng1,rng2);
    out[4] /= Nc2;
    out[5] /= Nc2;
  }

  void _hhc2ph_jet_base::amp_tree(ampq2g2p2 *amp1, ampq4p2 *amp2, double *out)
  {
    double tmp[2];



    //          TWO QUARKS TWO GLUONS SUBPROCESSES
    //  gg -> qqb   (1/Na * 1/Na)
    out[0] =(Nu*Qu4 + Nd*Qd4)*(amp1 -> su3_tree(2,1, -1, 0, -2,-3))/Na2;

    //  ug -> ug    (1/Nc * 1/Na)
    //  dg -> dg    (1/Nc * 1/Na)
    out[1] =Qu4*(tmp[0] = (amp1 -> su3_tree(-1,1, 0, 2, -2,-3))/(Nc*Na));
    out[2] =Qd4*tmp[0];

    //  gu -> ug    (1/Na * 1/Nc)
    //  gd -> dg    (1/Na * 1/Nc)
    out[3] =Qu4*(tmp[0] = (amp1 -> su3_tree(0,1, -1, 2, -2,-3))/(Nc*Na));
    out[4] =Qd4*tmp[0];

    //  uub -> gg   (1/Nc * 1/Nc * 1/2)
    out[12] = Qu4*(tmp[0] = 0.5*(amp1 -> su3_tree(-1,0, 1, 2, -2,-3))/Nc2);

    //  ddb -> gg   (1/Nc * 1/Nc * 1/2)
    out[14] =Qd4*tmp[0];


    //          FOUR QUARKS SUBPROCESSES
    //  uU -> uU    (1/Nc * 1/Nc)
    //  uu -> uu    (1/Nc * 1/Nc * 1/2)
    //  dD -> dD    (1/Nc * 1/Nc)
    //  dd -> dd    (1/Nc * 1/Nc * 1/2)
    amp2 -> su3_tree(Qu, -1,1, 0, 2, -2,-3, tmp);
    out[5] = tmp[0]/Nc2;
    out[6] = tmp[1]/(2.0*Nc2);
    out[7] =out[5]/16.0;
    out[8] = out[6]/16.0;

    //  ud -> ud    (1/Nc * 1/Nc)
   out[9] = (amp2 -> su3_tree(Qu, Qd, -1,1, 0,2, -2,-3))/Nc2;

    //  du -> du    (1/Nc * 1/Nc)
    out[10] =(amp2 -> su3_tree(Qd, Qu, -1,1, 0, 2, -2,-3))/Nc2;

    //  uUb -> uUb    (1/Nc * 1/Nc)
    //  uub -> uub    (1/Nc * 1/Nc)
    //  dDb -> dDb    (1/Nc * 1/Nc)
    //  ddb -> ddb    (1/Nc * 1/Nc)
    amp2 -> su3_tree(Qu, -1,1, 2, 0, -2,-3, tmp);
    out[11]  = tmp[0]/Nc2;
    out[12] += tmp[1]/Nc2;
    out[13]  = out[11]/16.0;
    out[14] += tmp[1]/(16.0*Nc2);

    //  uub -> UUb    (1/Nc * 1/Nc * (nu-1))
    //  uub -> ddb    (1/Nc * 1/Nc * nd)
   tmp[0] = (amp2 -> su3_tree(Qu, Qu, -1,0, 2, 1, -2,-3))/Nc2;
    out[12] += (Nu-1)*tmp[0]
      + Nd*(amp2 -> su3_tree(Qu, Qd, -1,0, 2, 1, -2,-3))/Nc2;

    //  ddb -> DDb    (1/Nc * 1/Nc * (nd-1))
    //  ddb -> uub    (1/Nc * 1/Nc * nu)
    out[14] += (Nd-1)*tmp[0]/16.0
      + Nu*(amp2 -> su3_tree(Qd, Qu, -1,0, 2, 1, -2,-3))/Nc2;

    //  udb -> udb    (1/Nc * 1/Nc)
    out[15] =(amp2 -> su3_tree(Qu, Qd, -1,1, 2, 0, -2,-3))/Nc2;

    //  dub -> dub    (1/Nc * 1/Nc)
    out[16] = (amp2 -> su3_tree(Qd, Qu, -1,1, 2, 0, -2,-3))/Nc2;
  }


  void _hhc2ph_jet_base::
  amp_tree_mch(ampq2g2p2 *amp1, ampq4p2 *amp2, double *out)
  {
    double tmp[2];

    //          TWO QUARKS TWO GLUONS SUBPROCESSES
    //  gg -> qqb   (1/Na * 1/Na)
    out[0] = (Nu*Qu4 + Nd*Qd4)*(amp1 -> su3_tree_mch(2,1, -1, 0, -2,-3))/Na2;

    //  ug -> ug    (1/Nc * 1/Na)
    //  dg -> dg    (1/Nc * 1/Na)
    out[1] = Qu4*(tmp[0] = (amp1 -> su3_tree_mch(-1,1, 0, 2, -2,-3))/(Nc*Na));
    out[2] = Qd4*tmp[0];

    //  gu -> ug    (1/Na * 1/Nc)
    //  gd -> dg    (1/Na * 1/Nc)
    out[3] = Qu4*(tmp[0] = (amp1 -> su3_tree_mch(0,1, -1, 2, -2,-3))/(Nc*Na));
    out[4] = Qd4*tmp[0];

    //  uub -> gg   (1/Nc * 1/Nc * 1/2)
    out[12] = Qu4*(tmp[0] = 0.5*(amp1 -> su3_tree_mch(-1,0, 1, 2, -2,-3))/Nc2);

    //  ddb -> gg   (1/Nc * 1/Nc * 1/2)
    out[14] = Qd4*tmp[0];


    //          FOUR QUARKS SUBPROCESSES
    //  uU -> uU    (1/Nc * 1/Nc)
    //  uu -> uu    (1/Nc * 1/Nc * 1/2)
    //  dD -> dD    (1/Nc * 1/Nc)
    //  dd -> dd    (1/Nc * 1/Nc * 1/2)
    amp2 -> su3_tree_mch(Qu, -1,1, 0, 2, -2,-3, tmp);
    out[5] = tmp[0]/Nc2;
    out[6] = tmp[1]/(2.0*Nc2);
    out[7] = out[5]/16.0;
    out[8] = out[6]/16.0;

    //  ud -> ud    (1/Nc * 1/Nc)
    out[9] = (amp2 -> su3_tree_mch(Qu, Qd, -1,1, 0,2, -2,-3))/Nc2;

    //  du -> du    (1/Nc * 1/Nc)
    out[10] = (amp2 -> su3_tree_mch(Qd, Qu, -1,1, 0, 2, -2,-3))/Nc2;

    //  uUb -> uUb    (1/Nc * 1/Nc)
    //  uub -> uub    (1/Nc * 1/Nc)
    //  dDb -> dDb    (1/Nc * 1/Nc)
    //  ddb -> ddb    (1/Nc * 1/Nc)
    amp2 -> su3_tree_mch(Qu, -1,1, 2, 0, -2,-3, tmp);
    out[11]  = tmp[0]/Nc2;
    out[12] += tmp[1]/Nc2;
    out[13]  = out[11]/16.0;
    out[14] += tmp[1]/(16.0*Nc2);

    //  uub -> UUb    (1/Nc * 1/Nc * (nu-1))
    //  uub -> ddb    (1/Nc * 1/Nc * nd)
    tmp[0] = (amp2 -> su3_tree_mch(Qu, Qu, -1,0, 2, 1, -2,-3))/Nc2;
    out[12] += (Nu-1)*tmp[0]
      + Nd*(amp2 -> su3_tree_mch(Qu, Qd, -1,0, 2, 1, -2,-3))/Nc2;

    //  ddb -> DDb    (1/Nc * 1/Nc * (nd-1))
    //  ddb -> uub    (1/Nc * 1/Nc * nu)
    out[14] += (Nd-1)*tmp[0]/16.0
      + Nu*(amp2 -> su3_tree_mch(Qd, Qu, -1,0, 2, 1, -2,-3))/Nc2;

    //  udb -> udb    (1/Nc * 1/Nc)
    out[15] = (amp2 -> su3_tree_mch(Qu, Qd, -1,1, 2, 0, -2,-3))/Nc2;

    //  dub -> dub    (1/Nc * 1/Nc)
    out[16] = (amp2 -> su3_tree_mch(Qd, Qu, -1,1, 2, 0, -2,-3))/Nc2;
  }


#define A(i) kp[i].tree
#define C(i) kp[i].cca
#define P(i) kp[i].pa
#define G(i) kp[i].ga


  void _hhc2ph_jet_base::
  __conv_x1(double eta, double x, double xjac, double al, const su3_kp_i2 *kp, weight_hhc2ph *S)
  {
    double lne = std::log(1.0-eta);
    double lie = lne*(lne - 2.0*std::log(eta)) - 2.0*__specfunc_li2(eta);

    //----- K term -----
    double k[4][2];
    Kgg(x, xjac, Nf, al, k[0]); Kgq(x, xjac, al, k[1]);
    Kqg(x, xjac, al, k[2]); Kqq(x, xjac, al, k[3]);
    k[0][1] += Ca*lie; k[3][1] += Cf*lie;

    S[0][0]  =                2.0*k[1][0]*(Nu*A(0) + Nd*A(1));
    S[0][1]  = k[3][0]*A(0);
    S[0][2]  = k[3][0]*A(1);
    S[0][3]  = k[0][0]*A(2) + k[1][0]*A(4);
    S[0][4]  = k[0][0]*A(3) + k[1][0]*A(5);
    S[0][12] = k[3][0]*A(4) + k[2][0]*A(2);
    S[0][14] = k[3][0]*A(5) + k[2][0]*A(3);
    S[0][5] = S[0][6] = S[0][10] = S[0][11] = S[0][16] = k[2][0]*A(2);
    S[0][7] = S[0][8] = S[0][ 9] = S[0][13] = S[0][15] = k[2][0]*A(3);

    //----- Ktilde term -----
    double tmp1, tmp2, t[4][2];
    Tgg(x, xjac, al, t[0]); Tgq(x, xjac, al, t[1]);
    Tqg(x, xjac, al, t[2]); Tqq(x, xjac, al, t[3]);
    t[0][1] += Ca*lne*lne; t[3][1] += Cf*lne*lne;

    S[0][0]  +=                2.0*t[1][0]*(Nu*C(0) + Nd*C(1));
    S[0][1]  += t[3][0]*C(0);
    S[0][2]  += t[3][0]*C(1);
    S[0][3]  += t[0][0]*C(2) + t[1][0]*C(4);
    S[0][4]  += t[0][0]*C(3) + t[1][0]*C(5);
    S[0][12] += t[3][0]*C(4) + t[2][0]*C(2);
    S[0][14] += t[3][0]*C(5) + t[2][0]*C(3);

    tmp1 = t[2][0]*C(2);
    tmp2 = t[2][0]*C(3);

    S[0][5] += tmp1; S[0][6] += tmp1; S[0][10] += tmp1; S[0][11] += tmp1;
    S[0][16] += tmp1; S[0][7] += tmp2; S[0][8] += tmp2; S[0][9] += tmp2;
    S[0][13] += tmp2; S[0][15] += tmp2;

    //----- P term -----
    double p[4][2];
    Pgg(x, xjac, Nf, p[0]); Pgq(x, xjac, p[1]);
    Pqg(x, xjac, p[2]); Pqq(x, xjac, p[3]);
    p[0][1] += 2.0*Ca*lne; p[3][1] += 2.0*Cf*lne;

    S[0][0]  +=                2.0*p[1][0]*(Nu*P(0) + Nd*P(1));
    S[0][1]  += p[3][0]*P(0);
    S[0][2]  += p[3][0]*P(1);
    S[0][3]  += p[0][0]*P(2) + p[1][0]*P(4);
    S[0][4]  += p[0][0]*P(3) + p[1][0]*P(5);
    S[0][12] += p[3][0]*P(4) + p[2][0]*P(2);
    S[0][14] += p[3][0]*P(5) + p[2][0]*P(3);

    tmp1 = p[2][0]*P(2);
    tmp2 = p[2][0]*P(3);

    S[0][5] += tmp1; S[0][6] += tmp1; S[0][10] += tmp1; S[0][11] += tmp1;
    S[0][16] += tmp1; S[0][7] += tmp2; S[0][8] += tmp2; S[0][9] += tmp2;
    S[0][13] += tmp2; S[0][15] += tmp2;

    //----- G term -----
    double g[2];
    g[0] = (x > 1.0-al ? xjac/(x-x*x) : 0.0);
    g[1] = -xjac/(1.0-x) + al-std::log(al) + lne;

    S[0][1]  += g[0]*G(0);
    S[0][2]  += g[0]*G(1);
    S[0][3]  += g[0]*G(2);
    S[0][4]  += g[0]*G(3);
    S[0][12] += g[0]*G(4);
    S[0][14] += g[0]*G(5);

    //----- dirac delta and + description terms -----
    S[2][1]  = g[1]*G(0) + k[3][1]*A(0) + t[3][1]*C(0) + p[3][1]*P(0);
    S[2][2]  = g[1]*G(1) + k[3][1]*A(1) + t[3][1]*C(1) + p[3][1]*P(1);
    S[2][3]  = g[1]*G(2) + k[0][1]*A(2) + t[0][1]*C(2) + p[0][1]*P(2);
    S[2][4]  = g[1]*G(3) + k[0][1]*A(3) + t[0][1]*C(3) + p[0][1]*P(3);
    S[2][12] = g[1]*G(4) + k[3][1]*A(4) + t[3][1]*C(4) + p[3][1]*P(4);
    S[2][14] = g[1]*G(5) + k[3][1]*A(5) + t[3][1]*C(5) + p[3][1]*P(5);
    S[2][5] = S[2][6] = S[2][10] = S[2][11] = S[2][16] = 0.0;
    S[2][7] = S[2][8] = S[2][ 9] = S[2][13] = S[2][15] = 0.0;

    //---- factorization scale dependent term ----
    S[3][0]  =               -2.0*p[1][0]*(Nu*A(0) + Nd*A(1));
    S[3][1]  = -p[3][0]*A(0);
    S[3][2]  = -p[3][0]*A(1);
    S[3][3]  = -p[0][0]*A(2) - p[1][0]*A(4);
    S[3][4]  = -p[0][0]*A(3) - p[1][0]*A(5);
    S[3][12] = -p[3][0]*A(4) - p[2][0]*A(2);
    S[3][14] = -p[3][0]*A(5) - p[2][0]*A(3);
    S[3][5] = S[3][6] = S[3][10] = S[3][11] = S[3][16] = -p[2][0]*A(2);
    S[3][7] = S[3][8] = S[3][ 9] = S[3][13] = S[3][15] = -p[2][0]*A(3);

    S[5][1]  = -p[3][1]*A(0);
    S[5][2]  = -p[3][1]*A(1);
    S[5][3]  = -p[0][1]*A(2);
    S[5][4]  = -p[0][1]*A(3);
    S[5][12] = -p[3][1]*A(4);
    S[5][14] = -p[3][1]*A(5);
    S[5][5] = S[5][6] = S[5][10] = S[5][11] = S[5][16] = 0.0;
    S[5][7] = S[5][8] = S[5][ 9] = S[5][13] = S[5][15] = 0.0;
  }

#undef C
#undef P
#undef G

#define C(i) kp[i].ccb
#define P(i) kp[i].pb
#define G(i) kp[i].gb


  void _hhc2ph_jet_base::
  __conv_x2(double eta, double x, double xjac, double al, const su3_kp_i2 *kp, weight_hhc2ph *S)
  {
    double lne = std::log(1.0-eta);
    double lie = lne*(lne - 2.0*std::log(eta)) - 2.0*__specfunc_li2(eta);

    //----- K term -----
    double k[4][2];
    Kgg(x, xjac, Nf, al, k[0]); Kgq(x, xjac, al, k[1]);
    Kqg(x, xjac, al, k[2]); Kqq(x, xjac, al, k[3]);
    k[0][1] += Ca*lie; k[3][1] += Cf*lie;

    S[1][0]  =                2.0*k[1][0]*(Nu*A(2) + Nd*A(3));
    S[1][1]  = k[0][0]*A(0) + k[1][0]*A(4);
    S[1][2]  = k[0][0]*A(1) + k[1][0]*A(5);
    S[1][3]  = k[3][0]*A(2);
    S[1][4]  = k[3][0]*A(3);
    S[1][12] = k[3][0]*A(4) + k[2][0]*A(0);
    S[1][14] = k[3][0]*A(5) + k[2][0]*A(1);
    S[1][5] = S[1][6] = S[1][ 9] = S[1][11] = S[1][15] = k[2][0]*A(0);
    S[1][7] = S[1][8] = S[1][10] = S[1][13] = S[1][16] = k[2][0]*A(1);

    //----- Ktilde term -----
    double tmp1, tmp2, t[4][2];
    Tgg(x, xjac, al, t[0]); Tgq(x, xjac, al, t[1]);
    Tqg(x, xjac, al, t[2]); Tqq(x, xjac, al, t[3]);
    t[0][1] += Ca*lne*lne; t[3][1] += Cf*lne*lne;

    S[1][0]  +=                2.0*t[1][0]*(Nu*C(2) + Nd*C(3));
    S[1][1]  += t[0][0]*C(0) + t[1][0]*C(4);
    S[1][2]  += t[0][0]*C(1) + t[1][0]*C(5);
    S[1][3]  += t[3][0]*C(2);
    S[1][4]  += t[3][0]*C(3);
    S[1][12] += t[3][0]*C(4) + t[2][0]*C(0);
    S[1][14] += t[3][0]*C(5) + t[2][0]*C(1);

    tmp1 = t[2][0]*C(0);
    tmp2 = t[2][0]*C(1);

    S[1][5] += tmp1; S[1][6] += tmp1; S[1][9] += tmp1; S[1][11] += tmp1;
    S[1][15] += tmp1; S[1][7] += tmp2; S[1][8] += tmp2; S[1][10] += tmp2;
    S[1][13] += tmp2; S[1][16] += tmp2;

    //----- P term -----
    double p[4][2];
    Pgg(x, xjac, Nf, p[0]); Pgq(x, xjac, p[1]);
    Pqg(x, xjac, p[2]); Pqq(x, xjac, p[3]);
    p[0][1] += 2.0*Ca*lne; p[3][1] += 2.0*Cf*lne;

    S[1][0]  +=                2.0*p[1][0]*(Nu*P(2) + Nd*P(3));
    S[1][1]  += p[0][0]*P(0) + p[1][0]*P(4);
    S[1][2]  += p[0][0]*P(1) + p[1][0]*P(5);
    S[1][3]  += p[3][0]*P(2);
    S[1][4]  += p[3][0]*P(3);
    S[1][12] += p[3][0]*P(4) + p[2][0]*P(0);
    S[1][14] += p[3][0]*P(5) + p[2][0]*P(1);

    tmp1 = p[2][0]*P(0);
    tmp2 = p[2][0]*P(1);

    S[1][5] += tmp1; S[1][6] += tmp1; S[1][9] += tmp1; S[1][11] += tmp1;
    S[1][15] += tmp1; S[1][7] += tmp2; S[1][8] += tmp2; S[1][10] += tmp2;
    S[1][13] += tmp2; S[1][16] += tmp2;

    //----- G term -----
    double g[2];
    g[0] = (x > 1.0-al ? xjac/(x-x*x) : 0.0);
    g[1] = -xjac/(1.0-x) + al-std::log(al) + lne;

    S[1][1]  += g[0]*G(0);
    S[1][2]  += g[0]*G(1);
    S[1][3]  += g[0]*G(2);
    S[1][4]  += g[0]*G(3);
    S[1][12] += g[0]*G(4);
    S[1][14] += g[0]*G(5);

    //----- dirac delta and + description terms -----
    S[2][1]  += g[1]*G(0) + k[0][1]*A(0) + t[0][1]*C(0) + p[0][1]*P(0);
    S[2][2]  += g[1]*G(1) + k[0][1]*A(1) + t[0][1]*C(1) + p[0][1]*P(1);
    S[2][3]  += g[1]*G(2) + k[3][1]*A(2) + t[3][1]*C(2) + p[3][1]*P(2);
    S[2][4]  += g[1]*G(3) + k[3][1]*A(3) + t[3][1]*C(3) + p[3][1]*P(3);
    S[2][12] += g[1]*G(4) + k[3][1]*A(4) + t[3][1]*C(4) + p[3][1]*P(4);
    S[2][14] += g[1]*G(5) + k[3][1]*A(5) + t[3][1]*C(5) + p[3][1]*P(5);

    //---- factorization scale dependent term ----
    S[4][0]  =               -2.0*p[1][0]*(Nu*A(2) + Nd*A(3));
    S[4][1]  = -p[0][0]*A(0) - p[1][0]*A(4);
    S[4][2]  = -p[0][0]*A(1) - p[1][0]*A(5);
    S[4][3]  = -p[3][0]*A(2);
    S[4][4]  = -p[3][0]*A(3);
    S[4][12] = -p[3][0]*A(4) - p[2][0]*A(0);
    S[4][14] = -p[3][0]*A(5) - p[2][0]*A(1);
    S[4][5] = S[1][6] = S[1][ 9] = S[1][11] = S[1][15] = -p[2][0]*A(0);
    S[4][7] = S[1][8] = S[1][10] = S[1][13] = S[1][16] = -p[2][0]*A(1);

    S[5][1]  -= p[0][1]*A(0);
    S[5][2]  -= p[0][1]*A(1);
    S[5][3]  -= p[3][1]*A(2);
    S[5][4]  -= p[3][1]*A(3);
    S[5][12] -= p[3][1]*A(4);
    S[5][14] -= p[3][1]*A(5);
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
 //////////////////////////////////////////////////////////////////////////////// 
 ////////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////// 
 // CT N^3LO gg
 ////////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////// 
  
  
  
 
//////////////////////////////////////////////////////////////////////////////// 
 ////////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////// 
 // QG AND QQBAR CHANNELS INTERFERENCE 2 -> 3
 ////////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////// 
  

