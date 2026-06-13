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

// -- 07/07/2024 -- Functions called from DYTurbo ------------------------------------------
// -----------------------------------------------------------------------------------------

// -- Elementos de matriz y funciones del Real a NLO interferencia y 

void ampq2g1p2::ggg_signal(int p1, int p2, int p3,int p4,int p5,
			   std::vector<std::complex<double> >& helamps) const
{
  double s45 =S(4,5); 
  complex<double> hprod = -3.0/4.0*mesq::SumaMartin(s45);
  double _VEV__ = sqrt(1./opts.Gf/sqrt(2.0));
  double Asq=(1.0/(3.0*__DEB_PI))*(1.0/(3.0*__DEB_PI))/ _VEV__/ _VEV__; //Ok
  double gsq=4.0*__DEB_PI; // Ok
  double xn=3.0;
  double V=8.0;
  //double fac=Asq*hdecay*gsq;// hdecay To be included directly in the Matrix Element
  double remainingFac = (2.0*__DEB_PI)*(2.0*__DEB_PI)*(2.0*__DEB_PI); // come from the definition of alpha_s which is as/2/pi
  complex<double> fac=sqrt(Asq*gsq*V*xn*remainingFac)*hprod; // Correction factor for production
  
  std::complex<double> photonsPP, photonsMM;
  complex<double> EffCoupl;
  EffCoupl = fac;
  // Hyy effective coupling constant WITH THE \alpha_QED REMOVED!!
  
  //double prefacsq=GF*pow( _MH_,4.0)/(8.0*sqrt(2)*__DEB_PI*__DEB_PI); //(esq/(4d0*pi))**2 to be included in process2ph....
  double prefacsq=opts.Gf*pow(1.0,4.0)/(8.0*sqrt(2.)*__DEB_PI*__DEB_PI); //(esq/(4d0*pi))**2 to be included in process2ph....
  
  complex<double>  msqgamgam = sqrt(prefacsq)*mesq::CorcheteMartin(s45); 
  
  complex<double> hdecayhel =  Higgs_prop(s45)*msqgamgam; // Corchete Martin directly evaluated in Mh as in MCFM fixed s45S
  // complex<double> hdecayhel = msqgamgam; // Corchete Martin directly evaluated in Mh as in MCFM

  //  photonsPP = -ICPLX*EffCoupl*Bb(4,5)*Bb(4,5)/sqrt(2.0); // Manca un B (Messo) 
  photonsPP =ICPLX*hdecayhel*EffCoupl*sqrt(S(1,2))*LSpb[4][5]*LSpb[4][5]*S(1,2);//*sqrt(S(1,2)) para ajustar a la definicion de los LSpx de Bern
  // El factor 1/sqrt(2.0) es para recuperar en la suma (Elemento de matrix al cuadrado para la señal) el factor MH^4
  //  photonsMM = -ICPLX*EffCoupl*Aa(4,5)*Aa(4,5)/sqrt(2.0); // Manca un A (Messo)
  photonsMM = ICPLX*hdecayhel*EffCoupl*sqrt(S(1,2))*LSpa[4][5]*LSpa[4][5]*S(1,2);//*sqrt(S(1,2)) para ajustar a la definicion de los LSpx de Bern
  
  //  (1+2+3+4+5+)
  
  helamps[0] = photonsPP*LS[4][5]*LS[4][5]/LSpa[1][2]/LSpa[2][3]/LSpa[3][1];  // OK
  
  //  (1-2+3+4+5+)
  
  helamps[1] = photonsPP*LSpb[3][2]*LSpb[3][2]*LSpb[3][2]/LSpb[2][1]/LSpb[1][3];   // OK
  
  //  (1+2-3+4+5+)
  
  helamps[2] = photonsPP*LSpb[1][3]*LSpb[1][3]*LSpb[1][3]/LSpb[1][2]/LSpb[2][3];// Ok signo menos viene della tracia
  
  //  (1+2+3-4+5+)
  
  helamps[3] = -photonsPP*LSpb[1][2]*LSpb[1][2]*LSpb[1][2]/LSpb[2][3]/LSpb[3][1];// Ok signo menos viene della tracia
  
  //  (1+2+3+4-5+)
  
  helamps[4] = std::complex<double>(0.,0.);
  
  //  (1+2+3+4+5-)
  
  helamps[5] = std::complex<double>(0.,0.);
  
  //  (1-2-3+4+5+)
  
  helamps[6] = photonsPP*LSpa[1][2]*LSpa[1][2]*LSpa[1][2]/LSpa[2][3]/LSpa[3][1];// Ok signo menos viene della tracia
  
  //  (1-2+3-4+5+)
  
  helamps[7] = -photonsPP*LSpa[1][3]*LSpa[1][3]*LSpa[1][3]/LSpa[1][2]/LSpa[2][3];// Ok signo menos viene della tracia
  
  //  (1-2+3+4-5+)
  
  helamps[8] = std::complex<double>(0.,0.);
  
  //  (1-2+3+4+5-)
  
  helamps[9] = std::complex<double>(0.,0.);
  
  //  (1+2-3-4+5+)
  
  helamps[10] = photonsPP*LSpa[2][3]*LSpa[2][3]*LSpa[2][3]/LSpa[2][1]/LSpa[1][3];   // OK
  
  //  (1+2-3+4-5+)
  
  helamps[11] = std::complex<double>(0.,0.);
  
  //  (1+2-3+4+5-)
  
  helamps[12] = std::complex<double>(0.,0.);
  
  //  (1+2+3-4-5+)
  
  helamps[13] = std::complex<double>(0.,0.);
  
  //  (1+2+3-4+5-)
  
  helamps[14] = std::complex<double>(0.,0.);
  
  //  (1+2+3+4-5-)
  
  helamps[15] = photonsMM*LS[4][5]*LS[4][5]/LSpa[1][2]/LSpa[2][3]/LSpa[3][1];   // OK
  
  //  (1-2-3-4-5-)
  
  helamps[16] = -photonsMM*LS[4][5]*LS[4][5]/LSpb[1][2]/LSpb[2][3]/LSpb[3][1];   // OK
  
  //  (1+2-3-4-5-)
  
  helamps[17] = photonsMM*LSpa[2][3]*LSpa[2][3]*LSpa[2][3]/LSpa[2][1]/LSpa[1][3];   // OK
  
  //  (1-2+3-4-5-)
  
  helamps[18] = -photonsMM*LSpa[1][3]*LSpa[1][3]*LSpa[1][3]/LSpa[1][2]/LSpa[2][3];// Ok signo menos viene della tracia
  
  //  (1-2-3+4-5-)
  
  helamps[19] = photonsMM*LSpa[1][2]*LSpa[1][2]*LSpa[1][2]/LSpa[2][3]/LSpa[3][1];// Ok signo menos viene della tracia
  
  //  (1-2-3-4+5-)
  
  helamps[20] = std::complex<double>(0.,0.);
  
  //  (1-2-3-4-5+)
  
  helamps[21] = std::complex<double>(0.,0.);
  
  //  (1+2+3-4-5-)
  
  helamps[22] = -photonsMM*LSpb[1][2]*LSpb[1][2]*LSpb[1][2]/LSpb[2][3]/LSpb[3][1];// Ok signo menos viene della tracia
  
  //  (1+2-3+4-5-)
  
  helamps[23] = photonsMM*LSpb[1][3]*LSpb[1][3]*LSpb[1][3]/LSpb[1][2]/LSpb[2][3];// Ok signo menos viene della tracia
  
  //  (1+2-3-4+5-)
  
  helamps[24] = std::complex<double>(0.,0.);
  
  //  (1+2-3-4-5+)
  
  helamps[25] = std::complex<double>(0.,0.);
  
  //  (1-2+3+4-5-)
  
  helamps[26] = photonsMM*LSpb[3][2]*LSpb[3][2]*LSpb[3][2]/LSpb[2][1]/LSpb[1][3];   // OK
  
  //  (1-2+3-4+5-)
  
  helamps[27] = std::complex<double>(0.,0.);
  
  //  (1-2+3-4-5+)
  
  helamps[28] = std::complex<double>(0.,0.);
  
  //  (1-2-3+4+5-)
  
  helamps[29] = std::complex<double>(0.,0.);
  
  //  (1-2-3+4-5+)
  
  helamps[30] = std::complex<double>(0.,0.);
  
  //  (1-2-3-4+5+)
  
  helamps[31] = -photonsPP*LS[4][5]*LS[4][5]/LSpb[1][2]/LSpb[2][3]/LSpb[3][1];   // OK
  
  /*
    cout << " Bb(1,3)   "  << Bb(1,3) << endl;
    cout << " Bb(3,1)   "  << Bb(3,1) << endl;
    cout << " LSpb[i3][i1]*sqrt(S(1,2))   "  << LSpb[3][1]*sqrt(S(1,2)) << endl;
    cout << " LSpb[i1][i3]*sqrt(S(1,2))   "  << LSpb[1][3]*sqrt(S(1,2)) << endl;
    cout << " LSpb[i1][i3]*sqrt(S(1,3))   "  << LSpb[1][3]*sqrt(S(1,3)) << endl;
    cout << " LSpb[i1][i2]*sqrt(S(1,2))   "  << LSpb[1][2]*sqrt(S(1,2)) << endl;
    cout << " Bb(1,2)   "  << Bb(1,2) << endl;
    cout << " LSpb[i1][i2]*LSab[2][1]   "  << LSpb[1][2]*LSpa[2][1] << endl;
    cout << " Bb(1,2)*Aa(2,1)/S(1,2)   "  << Bb(1,2)*Aa(2,1)/S(1,2) << endl;
    cout << " LSpb[i1][i3]*LSpb[3][1]   "  << LSpb[1][3]*LSpa[3][1] << endl;
    
    cout << " LS[1][3]   "  << LS[1][3] << endl;
    cout << " Bb(1,3)*Aa(3,1)/S(1,2)   "  << Bb(1,3)*Aa(3,1)/S(1,2) << endl;
    cout << " S(1,3)/S(1,2)   "  << S(1,3)/S(1,2) << endl;
    cout << " B(1,3)*A(3,1)/S(1,2)   "  << B(1,3)*A(3,1)/S(1,2) << endl;
    cout << " S(1,3)/S(1,2)   "  << S(1,3)/S(1,2) << endl;
    
    cout << " LS[2][5]   "  << LS[2][5] << endl;
    cout << " Bb(2,5)*Aa(5,2)/S(1,2)   "  << Bb(2,5)*Aa(5,2)/S(1,2) << endl;
    cout << " S(2,5)/S(1,2)   "  << S(2,5)/S(1,2) << endl;
    cout << " B(2,5)*A(5,2)/S(1,2)   "  << B(2,5)*A(5,2)/S(1,2) << endl;
    cout << " S(2,5)/S(1,2)   "  << S(2,5)/S(1,2) << endl;  
    
    cout << " LS[4][3]   "  << LS[4][3] << endl;
    cout << " Bb(4,3)*Aa(3,4)/S(1,2)   "  << Bb(4,3)*Aa(3,4)/S(1,2) << endl;
    cout << " S(4,3)/S(1,2)   "  << S(4,3)/S(1,2) << endl;
    cout << " B(4,3)*A(3,4)/S(1,2)   "  << B(4,3)*A(3,4)/S(1,2) << endl;
    cout << " S(4,3)/S(1,2)   "  << S(4,3)/S(1,2) << endl;
  */
  
  /*
    cout << " Bb(2,3)*Aa(3,2)/S(1,2)   "  << Bb(2,3)*Aa(3,2)/S(1,2) << endl;
    cout << " LSpb[2][3]*LSpa[3][2]    "  << LSpb[2][3]*LSpa[3][2]  << endl;
    cout << " LSpb[1][3]*LSpa[3][1]    "  << LSpb[1][3]*LSpa[3][1]  << endl;
    cout << " B(2,3)*A(3,2)/S(1,2)   "  << B(2,3)*A(3,2)/S(1,2) << endl;
    cout << " B(2,3)/S(1,2)   "  << B(2,3)/sqrt(S(1,2)) << endl;
    cout << " A(3,2)/S(1,2)   "  <<  A(3,2)/sqrt(S(1,2)) << endl;
    cout << " LSpb[2][3]    "  << LSpb[2][3]  << endl;
    cout << " LSpa[3][2]    "  << LSpa[3][2]  << endl;
    
  */
  
  /*
   for(int i=1;i<6;i++)
   for (int j=i+1;j<6;j++) {
   LSpa[j][i]=-LSpa[i][j];
   if ((i<3)&&(j!=2)) {
   LSpb[j][i]=-conj(LSpa[i][j]);
   } else {
   LSpb[j][i]=conj(LSpa[i][j]);
   }            
   LSpb[i][j]=-LSpb[j][i];
   LS[i][j]=LS[j][i]=real(LSpa[i][j]*LSpb[j][i]);
   }
   
  */
}


double ampq2g1p2::wgtHsig(int p1, int p2, int p3,int p4,int p5) const
{
  std::vector<std::complex<double> > helamps_sig(32,0.);
  std::vector<std::complex<double> > helamps_back(32,0.);
  double Msq=0., Msq_sig=0., Msq_int=0.;
  int ii;
  std::complex<double> H_prop;
  
  H_prop = Higgs_prop(S(4,5));
  double H_propnorm = norm(H_prop ); 
  
  double s45 =S(4,5); 
  double _MH__ = opts.hmass;
  double _GammaH__ = opts.hwidth;
  double prefacsq=(opts.Gf)*pow( _MH__,4.0)/(8.0*sqrt(2.)*__DEB_PI*__DEB_PI); //(esq/(4d0*pi))**2 to be included in process2ph....
  double msqgamgam = prefacsq*norm(mesq::CorcheteMartin(s45)); 
  
  double hdecay=msqgamgam/((s45- _MH__* _MH__)*(s45- _MH__* _MH__)+( _MH__* _GammaH__)*( _MH__* _GammaH__));// To be included directly in the Matrix Element
  double hdecayhel=msqgamgam* H_propnorm; // checked
  //c--   calculate propagators
  double ss=S(1,2);
  double tt=-S(1,3);
  double uu=-S(2,3);
  double mhsq=ss+tt+uu;
  //cout << "mhsq  :" << sqrt(mhsq) << endl;
  double avegg = 1/4.0/64.0; 
  double xn=3.0;
  double V=8.0;
  double MatrixsqPela=(pow(mhsq,4.0)+pow(ss,4.0)+pow(tt,4.0)+pow(uu,4.0))/(ss*tt*uu); // To be implemented with Msq (helamps)
     
  ggg_signal(p2,p1,p3,p4,p5,helamps_sig);

  for(ii=0;ii<32;ii++)
    {
      Msq_sig+=norm(helamps_sig[ii]);
    }
  Msq = Msq_sig/4.0; // facor 4 para que coincida Matrix Pela con Msq_sig
  double msq00=avegg*  Msq ; // checked
  return msq00;
  
}

double ampq2g1p2::wgtHint(int p1, int p2, int p3,int p4,int p5) const
{
  std::vector<std::complex<double> > helamps_sig(32,0.);
  std::vector<std::complex<double> > helamps_back(32,0.);
  double MsqRe=0.,MsqIm=0.,Msq;
  int ii;
  std::complex<double> H_prop;
  
  ggg_signal(p2,p1,p3,p4,p5,helamps_sig);
  Lggg_background(helamps_back);
  
  double s45 = S(4,5);

  for(ii=0;ii<32;ii++)
  {
  //  MsqIm+=imag(conj(helamps_back[ii])*helamps_sig[ii]);
      MsqRe+=real(conj(helamps_back[ii])*helamps_sig[ii]);
  }
  //  MsqIm*=2.0*imag( Higgs_prop(s45));
  MsqRe*= 2.0;
  //  Msq = MsqIm + MsqRe;
  Msq = MsqRe;

 double factorsDsq =3.0*Lzm*Lz*__DEB_PI/16.0*pow((3.0*(1.0/9.0)+2.0*(4.0/9.0)),2);
 double factorsD =sqrt(factorsDsq);
 double avegg =1.0/4.0/64.0;
 double CorFac = FACNew/avegg*factorsDsq;

 double CorFacrt = sqrt(CorFac);			   
 return Msq /4.0*avegg*CorFacrt;
}

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
  
  // p1p2 = 8.413e+06;
  // p1p4 = -7.22502e+06;
  // p1p5 = -817948;
  // p2p4 = -742599;
  // p2p5 = -7.58724e+06;
  // s45  = 15624.2;
  
  double _VEV__ = sqrt(1./opts.Gf/sqrt(2.0));
  
  double vmartin = _VEV__/sqrt(2.0); // Ok
  
  double Color=4.0;
  double vev= _VEV__;   //Reescribimos la amplitud para que aparezca 1/Vev**2 (usamos vHiggs=174 y vev=246)
  //  double alpha_qed= alpha_qed ;
  
  double vertice=opts.aemmz/3.0/__DEB_PI/(vev*vev)*3.0/2.0*real(mesq::CorcheteMartin(s45));
  //double vertice=-opts.aemmz/3.0/__DEB_PI/(vev*vev)*3.0/2.0*(16.0/9.0-8.3442);
  
  double promedioqq=1.0/(2.0*2.0*3.0*3.0);
  double MW = opts.wmass ;
  double SW = sqrt(opts.xw);
  //      cout << " SW : " << SW << endl;
  double el = sqrt(opts.aemmz*4.0*__DEB_PI);
  double MH = opts.hmass;
  double Gs = 1.0; // to be icluded at the end alpha_s
  double GH = opts.hwidth;
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
  double s45  = abs(S(4,5));
  
  // p1p2 = 8.413e+06;
  // p1p4 = -7.22502e+06;
  // p1p5 = -817948;
  // p2p4 = -742599;
  // p2p5 = -7.58724e+06;
  // s45  = 15624.2;
 
  /*
    cout << "From amp S(1,5)/2.0  :  " <<  S(2,3)/2.0 << endl;
    
    cout << "***************************************************************************** "<< endl;
    cout << "***************************************************************************** "<< endl;
    cout << "***************************************************************************** "<< endl;
  */
  
  double _VEV__ = sqrt(1./opts.Gf/sqrt(2.0)); 
  double vmartin = _VEV__/sqrt(2.0); // Ok
  
  double promedioqg=1.0/(2.0*2.0*3.0*8.0);
  double Color=4.0;
  double vev= _VEV__;     // !Reescribimos la amplitud para que aparezca 1/Vev**2 (usamos vHiggs=174 y vev=246)
  // double vertice=opts.aemmz/3.0/__DEB_PI/(vev*vev)*3.0/2.0*sqrt(norm(CorcheteMartin(s45)));
  double vertice=opts.aemmz/3.0/__DEB_PI/(vev*vev)*3.0/2.0*real(CorcheteMartin(s45));
  //  cout << "real(CorcheteMartin(s45)) = " << real(CorcheteMartin(s45)) << endl;
  // double vertice=-opts.aemmz/3.0/__DEB_PI/(vev*vev)*3.0/2.0*(16.0/9.0-8.3442);
  // Checkear a que configuracion de vértice pertenece estos números

  double MW = opts.wmass ;
  double SW = sqrt(opts.xw);
  double el = sqrt(opts.aemmz*4.0*__DEB_PI);
  double MH = opts.hmass;
  double Gs = 1.0; // to be icluded at the end alpha_s
  double GH = opts.hwidth;
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
      //      cout << " amp qg = " << amp2 << endl;
    return amp2;

 }
 
 
 
 
 
 
 















// extern double FACNew;


// Leandro Cieri /////////////////////////////////////////////////////
  double ampq2g1p2::LO_tree(int p1, int p2, int p3,int p4,int p5) const
  {
    double s13 = S(1,5), s14 = S(1,4);


    //tiene la suma sobre colores el flujo y el commonfac(1),
    // Promediado sobre colores y spines

    return Na*Nc*(s13/s14+s14/s13)*__DEB_4PI2/8.0;

 }

 double ampq2g1p2::BOXMTR(int p1, int p2, int p3,int p4,int p5) const
  {
    double t =S(1,5), u =S(1,4),
                s=S(1,2);


           double av=1.0/4.0/64.0;
          double xlst=log(-s/t);
          double xlsu=log(-s/u);
         double  xltu=log(t/u);
         double  rstu=(s*s+t*t)/u/u;
         double  rsut=(s*s+u*u)/t/t;
         double  rtus=(t*t+u*u)/s/s;
       double     pi2= __DEB_PI* __DEB_PI;
// missing charges (Sum Q_i^2)^2 * g^4   (g=alpha_s*4*pi)
      double      temp=av*16.0*pi2*2.0/(pi2*pi2)*(
       1.0/8.0*(( rstu*xlst*xlst+2.0*(s-t)/u*xlst )*( rstu*xlst*xlst+2.0*(s-t)/u*xlst )
      +( rsut*xlsu*xlsu+2.0*(s-u)/t*xlsu )*( rsut*xlsu*xlsu+2.0*(s-u)/t*xlsu )
      +( rtus*(xltu*xltu+pi2)+2.0*(t-u)/s*xltu )*( rtus*(xltu*xltu+pi2)+2.0*(t-u)/s*xltu ) )
      + 1.0/2.0 *
      ( rstu*xlst*xlst+2.0*(s-t)/u*xlst
        + rsut*xlsu*xlsu+2.0*(s-u)/t*xlsu
       + rtus*(xltu*xltu+pi2)+2.0*(t-u)/s*xltu )
      + pi2/2.0* ( (rstu*xlst+(s-t)/u)*(rstu*xlst+(s-t)/u)
      + (rsut*xlsu+(s-u)/t )*(rsut*xlsu+(s-u)/t ) )
       +4.0);
       //Falta dividir por 1/2 de particulas identicas!!!Listoen el 0.5!!!!!
     //  return __DEB_2PI4*4.0;

	std::complex<double> pppp = -1.0/2.0*(t*t+u*u)/(s*s)*(log(t/u)*log(t/u)+__DEB_PI*__DEB_PI)-(t-u)/s*log(t/u)-1.0;
	std::complex<double> mmmm = pppp;
	std::complex<double> ppmm = 1.0;
	std::complex<double> mmpp = ppmm;
	std::complex<double> pmpp = 1.0;
	std::complex<double> pmmp = -1.0/2.0*(t*t+s*s)/(u*u)*(log(-t/s)*log(-t/s))-(t-s)/u*log(-t/s)-1.0-ICPLX*__DEB_PI*((t*t+s*s)/(u*u)*log(-t/s)+(t-s)/u);
	std::complex<double> mppp = pmpp;
	std::complex<double> mpmp = -1.0/2.0*(u*u+s*s)/(t*t)*(log(-u/s)*log(-u/s))-(u-s)/t*log(-u/s)-1.0-ICPLX*__DEB_PI*((u*u+s*s)/(t*t)*log(-u/s)+(u-s)/t);
	std::complex<double> mmpm = pmpp;
	std::complex<double> mmmp = pmpp;
	std::complex<double> mpmm = pmpp;
	std::complex<double> mppm = pmmp;
	std::complex<double> pmmm = mppp;
	std::complex<double> pmpm = mpmp;
	std::complex<double> pppm = mmmp;
	std::complex<double> ppmp = mmpm;

	double tempNeri = real(pppp*conj(pppp)+pppm*conj(pppm)+ppmp*conj(ppmp)+ppmm*conj(ppmm)+pmpp*conj(pmpp)+pmpm*conj(pmpm)+pmmp*conj(pmmp)+pmmm*conj(pmmm)+mppp*conj(mppp)+mppm*conj(mppm)+mpmp*conj(mpmp)+mpmm*conj(mpmm)+mmpp*conj(mmpp)+mmpm*conj(mmpm)+mmmp*conj(mmmp)+mmmm*conj(mmmm));

//	return tempNeri*av*0.5*8.0*4.0*4.0*SumQ4*__DEB_2PI*__DEB_2PI;

	double BoxLean = temp*0.5*SumQ4*4.0*__DEB_2PI *__DEB_2PI*__DEB_2PI*__DEB_2PI ;

	double BoxNeri = tempNeri*av*0.5*8.0*4.0*4.0*SumQ4*__DEB_2PI*__DEB_2PI;

//	cout << BoxLean <<endl << BoxNeri <<endl << BoxLean/BoxNeri <<endl ;
// Ya está testeado, coinciden!

      return temp*0.5*SumQ4*4.0*__DEB_2PI *__DEB_2PI*__DEB_2PI*__DEB_2PI ;


// El 4 __DEB_2PI *__DEB_2PI Factores que antes venian desde el main
  //  return 4.0*__DEB_2PI *__DEB_2PI*__DEB_2PI*__DEB_2PI ; // El 4 __DEB_2PI *__DEB_2PI Factores que antes venian desde el main

 }

 double ampq2g1p2::LO_H1f(int p1, int p2, int p3,int p4,int p5) const
  {
     double t =S(1,5), u =S(2,5),
                s=S(1,2);
double s13 = S(1,5), s14 = S(1,4);


   double  hqq =  (Cf*((-21.0 + 2.0*__DEB_PI*__DEB_PI)*(t*t + u*u)
            + 3.0*u*(2.0*t + 3.0*u)*log(-(t/s)) +
            3.0*(2.0*t*t + 2.0*t*u + u*u)*log(-(t/s))*log(-(t/s)) +
            3.0*t*(3.0*t + 2.0*u)*log(-(u/s)) +
          3.0*(t*t+2.0*t*u+2.0*u*u)*log(-(u/s))*log(-(u/s))))/(3.0*(t*t+u*u)) ;


    return (hqq+Cf*__DEB_PI*__DEB_PI/3.0 )*Na*Nc*(s13/s14+s14/s13)*__DEB_4PI2/8.0/Nc2;

 }

///////////////////////////////////////////////////////////////////




  double ampq2g1p2::su3_tree(int p1, int p2, int p3,int p4,int p5) const
  {
    double s13 = S(1,3), s14 = S(1,4), s15 = S(1,5), s12 = S(1,2),
      s24 = S(2,4), s25 = S(2,5), s23 = S(2,3);

    return 8.0*Na*s12/(s23*s13*s14*s24*s15*s25)
      *(s13*s23*(s13*s13+s23*s23) + s14*s24*(s14*s14+s24*s24)
	+ s15*s25*(s15*s15+s25*s25));
  }

  std::complex<double>
  ampq2g1p2::amphtree(int p1, int p2, int p3, int p4, int p5) const
  {
    _ComplexD a12 = A(1,2), a13 = A(1,3), a14 = A(1,4), a15 = A(1,5),
      a23 = A(2,3), a24 = A(2,4), a25 = A(2,5);

    return -8.0*Na*a12*a12*(a15*a24/(a14*a25)+a14*a25/(a15*a24))/(a13*a13*a23*a23);
  }

  //  COLOR CORRALETED AMPLITUDES
#define Cond12(i,j)                                     \
(i == p1 && j == p2) || (i == p2 && j == p1)

#define Cond13(i,j)                                     \
(i == p1 && j == p3) || (i == p3 && j == p1) ||         \
(i == p3 && j == p2) || (i == p2 && j == p3)

  std::pair<double, std::complex<double> > ampq2g1p2::
  su3_cc(int pi,int pj, int p1, int p2, int p3, int p4, int p5) const
  {
    double s = 0.0, cc = su3_tree(p1, p2, p3, p4, p5);
    _ComplexD hcc = 0.0;

    if(Cond12(pi,pj))      s = 0.5/Nc;
    else if(Cond13(pi,pj)) s = -0.5*Nc;
    else throw("Error in ampq2g1p2::su3_cc(...)");

    if(pj == p3) hcc = amphtree(p1, p2, p3, p4, p5);

    return _Pair(s*cc, s*hcc);
  }

#define X(i,j) (-std::log(std::abs(Sij(i,j))))

  void ampq2g1p2::su3_kp(unsigned int nf, int pa, int pb, int p1,
			 int p2, int p3, int p4, int p5, su3_kp_i2 *res,
			 double al) const
  {
    double cc, xi, xj, xq = Gq/Cf, xg = Gg(nf)/Nc;
    int p[3] = {p1,p2,p3};

    double Ta = (pa == p1 || pa == p2 ? Cf : Ca);
    double Tb = (pb == p1 || pb == p2 ? Cf : Ca);
    double x[3] = {xq, xq, xg};

    res->tree = su3_tree(p1,p2, p3, p4,p5);
    res->pa = res->pb = res->ga = res->gb = res->loop = 0.0;

    for(unsigned i = 0; i < 3; i++)
      for(unsigned j = i+1; j < 3; j++) {
	if(Cond12(p[i], p[j]))      cc = 0.5*(res->tree)/Nc;
	else if(Cond13(p[i], p[j])) cc = -0.5*Nc*(res->tree);
	else throw("Error in ampq2g1p2::su3_kp(...)");

	//----- loop contribution -----
	if(i < 2) {
	  xi = Xq(Sij(p[i],p[j]), 1.0);
	  xj = j < 2 ? xi : Xg(Sij(p[i],p[j]), 1.0, nf);
	} else {
	  xi = Xg(Sij(p[i],p[j]), 1.0, nf);
	  xj = j < 2 ? Xq(Sij(p[i],p[j]), 1.0) : xi;
	}

	res->loop += (xi+xj)*cc;

	//----- finite contributions -----
	if(p[i] == pa || p[j] == pa) {
	  int iq = p[i] == pa ? j : i;
	  res->pa += X(pa, p[iq])*cc;
	  if(p[iq] != pb) res->ga += x[iq]*cc;
	  else res->cca = -cc/Ta;
	}

	if(p[i] == pb || p[j] == pb) {
	  int iq = p[i] == pb ? j : i;
	  res->pb += X(pb, p[iq])*cc;
	  if(p[iq] != pa) res->gb += x[iq]*cc;
	  else res->ccb = -cc/Tb;
	}
      }

    res->pa /= Ta; res->pb /= Tb;
    res->loop += (Gg(nf)+Kg(nf,al) + 2.0*(Gq+Kq(al)) - Cf)*(res->tree);
  }

#define Qu2 0.44444444444444444444
#define Qd2 0.11111111111111111111

  void ampq2g1p2::
  matrix_1loop(unsigned int nu, unsigned int nd, const _AmpPrim *A, const _AmpPrim *B, amp_1loop *res)
  {
    double TrQ2 = (nu*Qu2 + nd*Qd2);
    _ComplexD Asum(0.0), Bsum, Csum(0.0), Tsum(0.0);


    for(unsigned i = 0; i < 6; i++) {
      Tsum += A[i].Atree;
      Asum += A[i].Asusy - A[i].AL;
      Csum += A[i].Asf;
    }

    Bsum = B[0].Asusy-B[0].AL+B[0].Asf + B[1].Asusy-B[1].AL+B[1].Asf;

    res -> U0 = Qu2*Tsum;
    res -> D0 = Qd2*Tsum;

    res -> U1 = -Qu2*Bsum - Qu2*(Asum+Csum)/Nc2 - TrQ2*Csum/Nc;
    res -> D1 = -Qd2*Bsum - Qd2*(Asum+Csum)/Nc2 - TrQ2*Csum/Nc;
  }

#define __PrimitiveAmplitudes(h3,h4,h5)		\
  A1##h3##h4##h5(p1,p2, p3,p4,p5, A[0]);	\
  A1##h4##h5##h3(p1,p2, p4,p5,p3, A[1]);	\
  A1##h5##h3##h4(p1,p2, p5,p3,p4, A[2]);	\
  			       	                \
  A1##h3##h5##h4(p1,p2, p3,p5,p4, A[3]);	\
  A1##h5##h4##h3(p1,p2, p5,p4,p3, A[4]);	\
  A1##h4##h3##h5(p1,p2, p4,p3,p5, A[5]);	\
						\
						\
  A2##h3##h4##h5(p1, p3, p2, p4,p5, B[0]);	\
  A2##h3##h5##h4(p1, p3, p2, p5,p4, B[1])

  void ampq2g1p2::
  matrix_1loop_mpmpp(unsigned int nu, unsigned int nd, int p1, int p2,
		     int p3, int p4, int p5, amp_1loop *res) const
  {
    static _AmpPrim A[6], B[2];
    __PrimitiveAmplitudes(m,p,p);
    matrix_1loop(nu, nd, A, B, res);
  }

  void ampq2g1p2::
  matrix_1loop_mppmp(unsigned int nu, unsigned int nd, int p1, int p2,
		     int p3, int p4, int p5, amp_1loop *res) const
  {
    static _AmpPrim A[6], B[2];
    __PrimitiveAmplitudes(p,m,p);
    matrix_1loop(nu, nd, A, B, res);
  }

  void ampq2g1p2::
  matrix_1loop_mpppm(unsigned int nu, unsigned int nd, int p1, int p2,
		     int p3, int p4, int p5, amp_1loop *res) const
  {
    static _AmpPrim A[6], B[2];
    __PrimitiveAmplitudes(p,p,m);
    matrix_1loop(nu, nd, A, B, res);
  }

  void ampq2g1p2::
  matrix_1loop_mppmm(unsigned int nu, unsigned int nd, int p1, int p2,
		     int p3, int p4, int p5, amp_1loop *res) const
  {
    static _AmpPrim A[6], B[2];
    __PrimitiveAmplitudes(p,m,m);
    matrix_1loop(nu, nd, A, B, res);
  }

  void ampq2g1p2::
  matrix_1loop_mpmpm(unsigned int nu, unsigned int nd, int p1, int p2,
		     int p3, int p4, int p5, amp_1loop *res) const
  {
    static _AmpPrim A[6], B[2];
    __PrimitiveAmplitudes(m,p,m);
    matrix_1loop(nu, nd, A, B, res);
  }

  void ampq2g1p2::
  matrix_1loop_mpmmp(unsigned int nu, unsigned int nd, int p1, int p2,
		     int p3, int p4, int p5, amp_1loop *res) const
  {
    static _AmpPrim A[6], B[2];
    __PrimitiveAmplitudes(m,m,p);
    matrix_1loop(nu, nd, A, B, res);
  }

#define Amp1Loop				\
  out[0] += real(M.U1*conj(M.U0));		\
  out[1] += real(M.D1*conj(M.D0))


  void ampq2g1p2::
  su3_1loop(unsigned int nu, unsigned int nd, int p1, int p2,
	    int p3, int p4, int p5, double *out) const
  {
    static amp_1loop M;
    out[0] = out[1] = 0.0;

    matrix_1loop_mpppm(nu, nd, p1,p2,p3,p4,p5, &M); Amp1Loop;
    matrix_1loop_mppmp(nu, nd, p1,p2,p3,p4,p5, &M); Amp1Loop;
    matrix_1loop_mpmpp(nu, nd, p1,p2,p3,p4,p5, &M); Amp1Loop;
    matrix_1loop_mppmm(nu, nd, p1,p2,p3,p4,p5, &M); Amp1Loop;
    matrix_1loop_mpmpm(nu, nd, p1,p2,p3,p4,p5, &M); Amp1Loop;
    matrix_1loop_mpmmp(nu, nd, p1,p2,p3,p4,p5, &M); Amp1Loop;

    swap();
    matrix_1loop_mpppm(nu, nd, p1,p2,p3,p4,p5, &M); Amp1Loop;
    matrix_1loop_mppmp(nu, nd, p1,p2,p3,p4,p5, &M); Amp1Loop;
    matrix_1loop_mpmpp(nu, nd, p1,p2,p3,p4,p5, &M); Amp1Loop;
    matrix_1loop_mppmm(nu, nd, p1,p2,p3,p4,p5, &M); Amp1Loop;
    matrix_1loop_mpmpm(nu, nd, p1,p2,p3,p4,p5, &M); Amp1Loop;
    matrix_1loop_mpmmp(nu, nd, p1,p2,p3,p4,p5, &M); Amp1Loop;
    swap();

    out[0] *= 4.0*Na*Nc;
    out[1] *= 4.0*Na*Nc;
  }

  void ampq2g1p2::
  su3_1loop_mch(unsigned int nu, unsigned int nd, int p1, int p2,
		int p3, int p4, int p5, double *out) const
  {
    static amp_1loop M;
    unsigned int hpm = (unsigned int) (2*_M_rng());
    unsigned int hel = (unsigned int) (6*_M_rng());

    if(hpm == 1) swap();
    switch(hel){
    case 0: matrix_1loop_mpppm(nu, nd, p1,p2,p3,p4,p5, &M); break;
    case 1: matrix_1loop_mppmp(nu, nd, p1,p2,p3,p4,p5, &M); break;
    case 2: matrix_1loop_mpmpp(nu, nd, p1,p2,p3,p4,p5, &M); break;
    case 3: matrix_1loop_mppmm(nu, nd, p1,p2,p3,p4,p5, &M); break;
    case 4: matrix_1loop_mpmpm(nu, nd, p1,p2,p3,p4,p5, &M); break;
    case 5: matrix_1loop_mpmmp(nu, nd, p1,p2,p3,p4,p5, &M); break;
    }
    if(hpm == 1) swap();

    out[0] = 48.0*Na*Nc*real(M.U1*conj(M.U0));
    out[1] = 48.0*Na*Nc*real(M.D1*conj(M.D0));
  }

void ampq2g1p2::
  su3_1loop_mch(unsigned int nu, unsigned int nd, int p1, int p2,
		int p3, int p4, int p5, double *out,double rng1, double rng2) const
  {
    static amp_1loop M;
    unsigned int hpm = (unsigned int) (2*rng1);
    unsigned int hel = (unsigned int) (6*rng2);

    if(hpm == 1) swap();
    switch(hel){
    case 0: matrix_1loop_mpppm(nu, nd, p1,p2,p3,p4,p5, &M); break;
    case 1: matrix_1loop_mppmp(nu, nd, p1,p2,p3,p4,p5, &M); break;
    case 2: matrix_1loop_mpmpp(nu, nd, p1,p2,p3,p4,p5, &M); break;
    case 3: matrix_1loop_mppmm(nu, nd, p1,p2,p3,p4,p5, &M); break;
    case 4: matrix_1loop_mpmpm(nu, nd, p1,p2,p3,p4,p5, &M); break;
    case 5: matrix_1loop_mpmmp(nu, nd, p1,p2,p3,p4,p5, &M); break;
    }
    if(hpm == 1) swap();

    out[0] = 48.0*Na*Nc*real(M.U1*conj(M.U0));
    out[1] = 48.0*Na*Nc*real(M.D1*conj(M.D0));
  }


// Acá empiezan las reales (LeanR) //////////////////////////////////////////

double ampq2g1p2::wgtR(int p1, int p2, int p3,int p4,int p5) const
{
	
	/*
dixon pete;
	
double  z= pete.z;
double  zm = pete.zm;
double  w = pete.w;
double  wm = pete.wm;
double  w1 = pete.w1;
double  w1m = pete.w1m;
double  logz = pete.logz;
double  logzm = pete.logzm;
*/	
	
	//cout <<" Lw1 from amp  :" << Lw1 << endl;
	
	// cout <<" LSa[2][3]  :" << LSpa[2][3] << endl;

  std::vector<std::complex<double> > helamps(32,0.);
  double Msq=0.0;
  int ii;

 // ggg_background(p1,p2,p3,p4,p5,helamps);
  
  Lggg_background(helamps);

 for(ii=0;ii<32;ii++)
  {
    Msq+=norm(helamps[ii]);
  }

// This normalization is valid for the BACKGROUND but NOT FOR THE SIGNAL.
// It was cross-checked by Kora and with Dixon 0206194, and includes,
// along with the "standard" |M|^2 factors, the 1/2 for the identical photons
// and a (sort of random...) 4pi^2 which seems to be a relic of this code.
// The coupling constants alpha_s^2 alpha_QED^2 have been extracted and are to
// be restored elsewhere.


// TO BE CHECKED AGAIN!!!!!!
 // return Msq *SumQ4*3.0/(3.0*3.0-1.0)*16.*__DEB_PI*__DEB_PI*__DEB_PI;
 double factorsDsq =3.0*Lzm*Lz*__DEB_PI/16.0*pow((3.0*(1.0/9.0)+2.0*(4.0/9.0)),2);
 double factorsD =sqrt(factorsDsq);
 double avegg =1.0/4.0/64.0;
 double CorFac = FACNew/avegg*factorsDsq;
			   
			        //cout << "  FACNew desde AMPS  :" <<  FACNew << endl;
			        	//		        cout << "  Lz desde AMPS  :" <<  Lz << endl;

  // return Msq*3.0*Lzm*Lz/4.0*PI/16.0*pow((3.0*(1.0/9.0)+2.0*(4.0/9.0)),2);
     return Msq/4.0*CorFac*avegg; // Only to be clear how (and which) to include in the Matrix element

// IS THE FOLLOWING COMMENT MEANINGFUL?

// Matrix elements are multiplied by (4*pi)^2*sqrt(shat).

//  double pcrossGG=Nc*zm*z*Msq/4.0;

// sigma0=hbarCsqPB*PI/Qsq/Qsq/16.0*pow(alpha*as2pi*(3.0*(1.0/9.0)+2.0*(4.0/9.0)),2);

//  return sigma0*as2pi*pcrossGG*fgg;

// fgg=pdf(iparton,xA,muF)*pdf(iparton,xB,muF);
//

}





// Lean Checked again!!! Ok!
std::complex<double> ampq2g1p2::Higgs_prop(double ss) const
{

  return ICPLX/std::complex<double>(ss-_MH_*_MH_,_MH_*_GammaH_);

}


// Lean Checked again!!! Ok!
std::complex<double> ff( double xx)
{

if( xx >= 1.0 ) return asin( sqrt(1.0/xx) )*asin( sqrt(1.0/xx) );

	return -1.0/4.0*(  log((1.0+sqrt(1.0-xx))/(1.0-sqrt(1.0-xx))) - ICPLX* __DEB_PI)*(  log((1.0+sqrt(1.0-xx))/(1.0-sqrt(1.0-xx))) - ICPLX* __DEB_PI);
	
}

// Lean Checked again!!! Ok!
std::complex<double> F1( double xx)
{

	return 2.0 + 3.0*xx*(1.0+(2.0-xx)*ff(xx));
	
	
}

// Lean Checked again!!! Ok!
std::complex<double> F1on2( double xx)
{

 return -2.0*xx*(1.0+(1.0-xx)*ff(xx));
	
}



std::complex<double> CorcheteMartin(double ss) 
{


  return mesq::CorcheteMartin(ss);
  // with b and c
   return F1on2(4.0*mtau*mtau/ss) + F1(4.0*mW*mW/ss)+ 3.0*( (2.0/3.0)*(2.0/3.0)*F1on2(4.0*mc*mc/ss) + (2.0/3.0)*(2.0/3.0)*F1on2(4.0*mt*mt/ss) + (-1.0/3.0)*(-1.0/3.0)*F1on2(4.0*mb*mb/ss));

// without tau and c
// return  F1(4.0*mW*mW/ss)+ 3.0*( (2.0/3.0)*(2.0/3.0)*F1on2(4.0*mt*mt/ss) + (-1.0/3.0)*(-1.0/3.0)*F1on2(4.0*mb*mb/ss));

// without b and c
// return F1on2(4.0*mtau*mtau/ss)+F1(4.0*mW*mW/ss)+3.0*((2.0/3.0)*(2.0/3.0)*F1on2(4.0*mt*mt/ss));
 


}

std::complex<double> SumaMartin(double ss) 
{

// with b and c
// return F1on2(4.0*mt*mt/ss);

// without  c
// return (F1on2(4.0*mt*mt/ss)+F1on2(4.0*mb*mb/ss));

// without b and c
  return mesq::SumaMartin(ss);
 return (F1on2(4.0*mt*mt/ss)+F1on2(4.0*mb*mb/ss)+F1on2(4.0*mc*mc/ss));
 


}




std::complex<double> ampq2g1p2::A5yy(
             std::complex<double>(ampq2g1p2::*amp_base)(int,int,int,int,int)const,
                    int p1, int p2, int p3, int p4, int p5) const
{
  std::complex<double> res;

  res =
    (this->*amp_base)(p1,p2,p3,p4,p5)+(this->*amp_base)(p1,p4,p5,p2,p3)+(this->*amp_base)(p1,p2,p4,p5,p3)
   +(this->*amp_base)(p1,p2,p3,p5,p4)+(this->*amp_base)(p1,p5,p4,p2,p3)+(this->*amp_base)(p1,p2,p5,p4,p3)
   +(this->*amp_base)(p1,p2,p4,p3,p5)+(this->*amp_base)(p1,p5,p2,p3,p4)+(this->*amp_base)(p1,p4,p2,p5,p3)
   +(this->*amp_base)(p1,p2,p5,p3,p4)+(this->*amp_base)(p1,p4,p2,p3,p5)+(this->*amp_base)(p1,p5,p2,p4,p3);

   return res;

}

// Acá van las A_5^[1/2]_5;1 de la ec 4 de 9905283.
// Note that the actual expressions for the amplitudes for the various
// helicities can be actually found in 9302280 (which is ref. [9] of 9905283).
// With respect to the paper, the limit eps->0 has been taken and
// there is a different overall normalization by

/*

std::complex<double> ampq2g1p2::Mppppp(int pi1, int pi2, int pi3, int pi4, int pi5) const
{
return -(S(i1,i2)*S(i2,i3)+S(i2,i3)*S(i3,i4)+S(i3,i4)*S(i4,i5)
       +S(i4,i5)*S(i5,i1)+S(i5,i1)*S(i1,i2)
       +B(i1,i2)*A(i2,i3)*B(i3,i4)*A(i4,i1)
       -A(i1,i2)*B(i2,i3)*A(i3,i4)*B(i4,i1))
        /A(i1,i2)/A(i2,i3)/A(i3,i4)/A(i4,i5)/A(i5,i1)/6.0;
}

std::complex<double> ampq2g1p2::Mmmmmm(int pi1, int pi2, int pi3, int pi4, int pi5)   const
{
return (S(i1,i2)*S(i2,i3)+S(i2,i3)*S(i3,i4)+S(i3,i4)*S(i4,i5)
       +S(i4,i5)*S(i5,i1)+S(i5,i1)*S(i1,i2)
       +A(i1,i2)*B(i2,i3)*A(i3,i4)*B(i4,i1)
       -B(i1,i2)*A(i2,i3)*B(i3,i4)*A(i4,i1))
        /B(i1,i2)/B(i2,i3)/B(i3,i4)/B(i4,i5)/B(i5,i1)/6.0;
}

std::complex<double> ampq2g1p2::Mmpppp(int pi1, int pi2, int pi3, int pi4, int pi5)   const
{
return -((S(i2,i3)+S(i3,i4)+S(i4,i5))*B(i2,i5)*B(i2,i5)
        -B(i2,i4)*A(i4,i3)*B(i3,i5)*B(i2,i5)
        -B(i1,i2)*B(i1,i5)/A(i1,i2)/A(i1,i5)*(
            A(i1,i2)*A(i1,i2)*A(i1,i3)*A(i1,i3)
                     *B(i2,i3)/A(i2,i3)
           +A(i1,i3)*A(i1,i3)*A(i1,i4)*A(i1,i4)
                     *B(i3,i4)/A(i3,i4)
           +A(i1,i4)*A(i1,i4)*A(i1,i5)*A(i1,i5)
                     *B(i4,i5)/A(i4,i5)))
        /B(i1,i2)/A(i2,i3)/A(i3,i4)/A(i4,i5)/B(i5,i1)/3.0;
}

std::complex<double> ampq2g1p2::Mpmmmm(int pi1, int pi2, int pi3, int pi4, int pi5)   const
{
return ((S(i2,i3)+S(i3,i4)+S(i4,i5))*A(i2,i5)*A(i2,i5)
        -A(i2,i4)*B(i4,i3)*A(i3,i5)*A(i2,i5)
        -A(i1,i2)*A(i1,i5)/B(i1,i2)/B(i1,i5)*(
            B(i1,i2)*B(i1,i2)*B(i1,i3)*B(i1,i3)
                     *A(i2,i3)/B(i2,i3)
           +B(i1,i3)*B(i1,i3)*B(i1,i4)*B(i1,i4)
                     *A(i3,i4)/B(i3,i4)
           +B(i1,i4)*B(i1,i4)*B(i1,i5)*B(i1,i5)
                     *A(i4,i5)/B(i4,i5)))
        /A(i1,i2)/B(i2,i3)/B(i3,i4)/B(i4,i5)/A(i5,i1)/3.0;
}

std::complex<double> ampq2g1p2::Mmmppp(int pi1, int pi2, int pi3, int pi4, int pi5)  const
{

  //  Note: nothing should depend on mm.

  double mm=1.0;

  std::complex<double> l23,l51;

  if (S(i2,i3)<0) {
     l23=log(mm/(-S(i2,i3)));
  } else {
     l23=log(mm/S(i2,i3))+ICPLX*__DEB_PI;
  }

  if (S(i5,i1)<0) {
     l51=log(mm/(-S(i5,i1)));
  } else {
     l51=log(mm/S(i5,i1))+ICPLX*__DEB_PI;
  }


  std::complex<double> Vf=-0.5*(l23+l51);
  std::complex<double> Vs=-1.0/3.0*Vf;

  std::complex<double> tree=pow(A(i1,i2),4)/A(i1,i2)/A(i2,i3)/A(i3,i4)
                      /A(i4,i5)/A(i5,i1);

  std::complex<double> Ff=-0.5*A(i1,i2)*A(i1,i2)*(
                    A(i2,i3)*B(i3,i4)*A(i4,i1)
                   +A(i2,i4)*B(i4,i5)*A(i5,i1))
              *L0(-S(i2,i3),-S(i5,i1))/S(i5,i1)
           /A(i2,i3)/A(i3,i4)/A(i4,i5)/A(i5,i1);

  std::complex<double> Fs=-1.0/3.0*B(i3,i4)*A(i4,i1)*A(i2,i4)*B(i4,i5)*(
                    A(i2,i3)*B(i3,i4)*A(i4,i1)
                   +A(i2,i4)*B(i4,i5)*A(i5,i1))
              *L2(-S(i2,i3),-S(i5,i1))/S(i5,i1)/S(i5,i1)/S(i5,i1)
           /A(i3,i4)/A(i4,i5)
             -1.0/3.0*Ff
             -1.0/3.0*A(i3,i5)*B(i3,i5)*B(i3,i5)*B(i3,i5)
                   /B(i1,i2)/B(i2,i3)/A(i3,i4)
                      /A(i4,i5)/B(i5,i1)
             +1.0/3.0*A(i1,i2)*B(i3,i5)*B(i3,i5)
                   /B(i2,i3)/A(i3,i4)
                      /A(i4,i5)/B(i5,i1)
             +1.0/6.0*A(i1,i2)*B(i3,i4)
                      *A(i4,i1)*A(i2,i4)*B(i4,i5)
                    /S(i2,i3)/A(i3,i4)/A(i4,i5)/S(i5,i1);

return -( (Vf+Vs)*tree +Ff+Fs );
}
std::complex<double> ampq2g1p2::Mppmmm(int pi1, int pi2, int pi3, int pi4, int pi5)  const
{

  //  Note: nothing should depend on mm.

  double mm=1.0;

  std::complex<double> l23,l51;

  if (S(i2,i3)<0) {
     l23=log(mm/(-S(i2,i3)));
  } else {
     l23=log(mm/S(i2,i3))+ICPLX*__DEB_PI;
  }

  if (S(i5,i1)<0) {
     l51=log(mm/(-S(i5,i1)));
  } else {
     l51=log(mm/S(i5,i1))+ICPLX*__DEB_PI;
  }


  std::complex<double> Vf=-0.5*(l23+l51);
  std::complex<double> Vs=-1.0/3.0*Vf;

  std::complex<double> tree=pow(B(i1,i2),4)/B(i1,i2)/B(i2,i3)/B(i3,i4)
                      /B(i4,i5)/B(i5,i1);

  std::complex<double> Ff=-0.5*B(i1,i2)*B(i1,i2)*(
                    B(i2,i3)*A(i3,i4)*B(i4,i1)
                   +B(i2,i4)*A(i4,i5)*B(i5,i1))
              *L0(-S(i2,i3),-S(i5,i1))/S(i5,i1)
           /B(i2,i3)/B(i3,i4)/B(i4,i5)/B(i5,i1);

  std::complex<double> Fs=-1.0/3.0*A(i3,i4)*B(i4,i1)*B(i2,i4)*A(i4,i5)*(
                    B(i2,i3)*A(i3,i4)*B(i4,i1)
                   +B(i2,i4)*A(i4,i5)*B(i5,i1))
              *L2(-S(i2,i3),-S(i5,i1))/S(i5,i1)/S(i5,i1)/S(i5,i1)
           /B(i3,i4)/B(i4,i5)
             -1.0/3.0*Ff
             -1.0/3.0*B(i3,i5)*A(i3,i5)*A(i3,i5)*A(i3,i5)
                   /A(i1,i2)/A(i2,i3)/B(i3,i4)
                      /B(i4,i5)/A(i5,i1)
             +1.0/3.0*B(i1,i2)*A(i3,i5)*A(i3,i5)
                   /A(i2,i3)/B(i3,i4)
                      /B(i4,i5)/A(i5,i1)
             +1.0/6.0*B(i1,i2)*A(i3,i4)
                      *B(i4,i1)*B(i2,i4)*A(i4,i5)
                    /S(i2,i3)/B(i3,i4)/B(i4,i5)/S(i5,i1);

return ( (Vf+Vs)*tree +Ff+Fs );
}

std::complex<double> ampq2g1p2::Mmpmpp(int pi1, int pi2, int pi3, int pi4, int pi5)  const
{

  //  Note: nothing should depend on mm.

  double mm=1.0;

  std::complex<double> l34,l51;

  if (S(i3,i4)<0) {
     l34=log(mm/(-S(i3,i4)));
  } else {
     l34=log(mm/S(i3,i4))+ICPLX*__DEB_PI;
  }

  if (S(i5,i1)<0) {
     l51=log(mm/(-S(i5,i1)));
  } else {
     l51=log(mm/S(i5,i1))+ICPLX*__DEB_PI;
  }

  std::complex<double> Vf=-0.5*(l34+l51);
  std::complex<double> Vs=-1.0/3.0*Vf;

  std::complex<double> tree=pow(A(i1,i3),4)/A(i1,i2)/A(i2,i3)/A(i3,i4)
                      /A(i4,i5)/A(i5,i1);

  std::complex<double> Ff=-A(i1,i3)*A(i1,i3)*A(i4,i1)*B(i2,i4)*B(i2,i4)
                  *Ls1(-S(i2,i3),-S(i5,i1),-S(i3,i4),-S(i5,i1))
                   /A(i4,i5)/A(i5,i1)/S(i5,i1)/S(i5,i1)
             +A(i1,i3)*A(i1,i3)*A(i5,i3)*B(i2,i5)*B(i2,i5)
                  *Ls1(-S(i1,i2),-S(i3,i4),-S(i5,i1),-S(i3,i4))
                   /A(i3,i4)/A(i4,i5)/S(i3,i4)/S(i3,i4)
             -0.5*A(i1,i3)*A(i1,i3)*A(i1,i3)*(
                    A(i1,i5)*B(i5,i2)*A(i2,i3)
                   -A(i3,i4)*B(i4,i2)*A(i2,i1))
              *L0(-S(i3,i4),-S(i5,i1))/S(i5,i1)
           /A(i1,i2)/A(i2,i3)/A(i3,i4)/A(i4,i5)/A(i5,i1);

  std::complex<double> Fs=-A(i1,i2)*A(i2,i3)*A(i3,i4)*A(i4,i1)
                      *A(i4,i1)*B(i2,i4)*B(i2,i4)*(
                2.0*Ls1(-S(i2,i3),-S(i5,i1),-S(i3,i4),-S(i5,i1))
                   +L1(-S(i2,i3),-S(i5,i1))
                   +L1(-S(i3,i4),-S(i5,i1)))
                  /A(i4,i5)/A(i5,i1)/A(i2,i4)/A(i2,i4)
                  /S(i5,i1)/S(i5,i1)
             +A(i3,i2)*A(i2,i1)*A(i1,i5)*A(i5,i3)
                      *A(i5,i3)*B(i2,i5)*B(i2,i5)*(
                2.0*Ls1(-S(i1,i2),-S(i3,i4),-S(i5,i1),-S(i3,i4))
                   +L1(-S(i1,i2),-S(i3,i4))
                   +L1(-S(i5,i1),-S(i3,i4)))
                  /A(i5,i4)/A(i4,i3)/A(i2,i5)/A(i2,i5)
                  /S(i3,i4)/S(i3,i4)
             +2.0/3.0*A(i2,i3)*A(i2,i3)
                     *A(i4,i1)*A(i4,i1)*A(i4,i1)
                     *B(i2,i4)*B(i2,i4)*B(i2,i4)
                     *L2(-S(i2,i3),-S(i5,i1))
                     /A(i4,i5)/A(i5,i1)/A(i2,i4)
                     /S(i5,i1)/S(i5,i1)/S(i5,i1)
             -2.0/3.0*A(i2,i1)*A(i2,i1)
                     *A(i5,i3)*A(i5,i3)*A(i5,i3)
                     *B(i2,i5)*B(i2,i5)*B(i2,i5)
                     *L2(-S(i1,i2),-S(i3,i4))
                     /A(i5,i4)/A(i4,i3)/A(i2,i5)
                     /S(i3,i4)/S(i3,i4)/S(i3,i4)
             +L2(-S(i3,i4),-S(i5,i1))/S(i5,i1)/S(i5,i1)/S(i5,i1)*(
	        1.0/3.0*A(i1,i3)*B(i2,i4)*B(i2,i5)*(
                    A(i1,i5)*B(i5,i2)*A(i2,i3)
                   -A(i3,i4)*B(i4,i2)*A(i2,i1))/A(i4,i5)
               +2.0/3.0*A(i1,i2)*A(i1,i2)*A(i3,i4)*A(i3,i4)
                    *A(i4,i1)*B(i2,i4)*B(i2,i4)*B(i2,i4)
                    /A(i4,i5)/A(i5,i1)/A(i2,i4)
               -2.0/3.0*A(i3,i2)*A(i3,i2)*A(i1,i5)*A(i1,i5)
                    *A(i5,i3)*B(i2,i5)*B(i2,i5)*B(i2,i5)
                    /A(i5,i4)/A(i4,i3)/A(i2,i5))
            +1.0/6.0*A(i1,i3)*A(i1,i3)*A(i1,i3)*(
                    A(i1,i5)*B(i5,i2)*A(i2,i3)
                   -A(i3,i4)*B(i4,i2)*A(i2,i1))
              *L0(-S(i3,i4),-S(i5,i1))/S(i5,i1)
           /A(i1,i2)/A(i2,i3)/A(i3,i4)/A(i4,i5)/A(i5,i1)
             +1.0/3.0*B(i2,i4)*B(i2,i4)*B(i2,i5)*B(i2,i5)
           /B(i1,i2)/B(i2,i3)/B(i3,i4)/A(i4,i5)/B(i5,i1)
             -1.0/3.0*A(i1,i2)*A(i4,i1)*A(i4,i1)
                     *B(i2,i4)*B(i2,i4)*B(i2,i4)
                   /A(i4,i5)/A(i5,i1)/A(i2,i4)
                      /B(i2,i3)/B(i3,i4)/S(i5,i1)
             +1.0/3.0*A(i3,i2)*A(i5,i3)*A(i5,i3)
                     *B(i2,i5)*B(i2,i5)*B(i2,i5)
                   /A(i5,i4)/A(i4,i3)/A(i2,i5)
                      /B(i2,i1)/B(i1,i5)/S(i3,i4)
             +1.0/6.0*A(i1,i3)*A(i1,i3)*B(i2,i4)*B(i2,i5)
                      /S(i3,i4)/A(i4,i5)/S(i5,i1);

return -( (Vf+Vs)*tree +Ff+Fs );
}

std::complex<double> ampq2g1p2::Mpmpmm(int pi1, int pi2, int pi3, int pi4, int pi5)   const
{

  //  Note: nothing should depend on mm.

  double mm=1.0;

  std::complex<double> l34,l51;

  if (S(i3,i4)<0) {
     l34=log(mm/(-S(i3,i4)));
  } else {
     l34=log(mm/S(i3,i4))+ICPLX*__DEB_PI;
  }

  if (S(i5,i1)<0) {
     l51=log(mm/(-S(i5,i1)));
  } else {
     l51=log(mm/S(i5,i1))+ICPLX*__DEB_PI;
  }

  std::complex<double> Vf=-0.5*(l34+l51);
  std::complex<double> Vs=-1.0/3.0*Vf;

  std::complex<double> tree=pow(B(i1,i3),4)/B(i1,i2)/B(i2,i3)/B(i3,i4)
                      /B(i4,i5)/B(i5,i1);

  std::complex<double> Ff=-B(i1,i3)*B(i1,i3)*B(i4,i1)*A(i2,i4)*A(i2,i4)
                  *Ls1(-S(i2,i3),-S(i5,i1),-S(i3,i4),-S(i5,i1))
                   /B(i4,i5)/B(i5,i1)/S(i5,i1)/S(i5,i1)
             +B(i1,i3)*B(i1,i3)*B(i5,i3)*A(i2,i5)*A(i2,i5)
                  *Ls1(-S(i1,i2),-S(i3,i4),-S(i5,i1),-S(i3,i4))
                   /B(i3,i4)/B(i4,i5)/S(i3,i4)/S(i3,i4)
             -0.5*B(i1,i3)*B(i1,i3)*B(i1,i3)*(
                    B(i1,i5)*A(i5,i2)*B(i2,i3)
                   -B(i3,i4)*A(i4,i2)*B(i2,i1))
              *L0(-S(i3,i4),-S(i5,i1))/S(i5,i1)
           /B(i1,i2)/B(i2,i3)/B(i3,i4)/B(i4,i5)/B(i5,i1);

  std::complex<double> Fs=-B(i1,i2)*B(i2,i3)*B(i3,i4)*B(i4,i1)
                      *B(i4,i1)*A(i2,i4)*A(i2,i4)*(
                2.0*Ls1(-S(i2,i3),-S(i5,i1),-S(i3,i4),-S(i5,i1))
                   +L1(-S(i2,i3),-S(i5,i1))
                   +L1(-S(i3,i4),-S(i5,i1)))
                  /B(i4,i5)/B(i5,i1)/B(i2,i4)/B(i2,i4)
                  /S(i5,i1)/S(i5,i1)
             +B(i3,i2)*B(i2,i1)*B(i1,i5)*B(i5,i3)
                      *B(i5,i3)*A(i2,i5)*A(i2,i5)*(
                2.0*Ls1(-S(i1,i2),-S(i3,i4),-S(i5,i1),-S(i3,i4))
                   +L1(-S(i1,i2),-S(i3,i4))
                   +L1(-S(i5,i1),-S(i3,i4)))
                  /B(i5,i4)/B(i4,i3)/B(i2,i5)/B(i2,i5)
                  /S(i3,i4)/S(i3,i4)
             +2.0/3.0*B(i2,i3)*B(i2,i3)
                     *B(i4,i1)*B(i4,i1)*B(i4,i1)
                     *A(i2,i4)*A(i2,i4)*A(i2,i4)
                     *L2(-S(i2,i3),-S(i5,i1))
                     /B(i4,i5)/B(i5,i1)/B(i2,i4)
                     /S(i5,i1)/S(i5,i1)/S(i5,i1)
             -2.0/3.0*B(i2,i1)*B(i2,i1)
                     *B(i5,i3)*B(i5,i3)*B(i5,i3)
                     *A(i2,i5)*A(i2,i5)*A(i2,i5)
                     *L2(-S(i1,i2),-S(i3,i4))
                     /B(i5,i4)/B(i4,i3)/B(i2,i5)
                     /S(i3,i4)/S(i3,i4)/S(i3,i4)
             +L2(-S(i3,i4),-S(i5,i1))/S(i5,i1)/S(i5,i1)/S(i5,i1)*(
	        1.0/3.0*B(i1,i3)*A(i2,i4)*A(i2,i5)*(
                    B(i1,i5)*A(i5,i2)*B(i2,i3)
                   -B(i3,i4)*A(i4,i2)*B(i2,i1))/B(i4,i5)
               +2.0/3.0*B(i1,i2)*B(i1,i2)*B(i3,i4)*B(i3,i4)
                    *B(i4,i1)*A(i2,i4)*A(i2,i4)*A(i2,i4)
                    /B(i4,i5)/B(i5,i1)/B(i2,i4)
               -2.0/3.0*B(i3,i2)*B(i3,i2)*B(i1,i5)*B(i1,i5)
                    *B(i5,i3)*A(i2,i5)*A(i2,i5)*A(i2,i5)
                    /B(i5,i4)/B(i4,i3)/B(i2,i5))
            +1.0/6.0*B(i1,i3)*B(i1,i3)*B(i1,i3)*(
                    B(i1,i5)*A(i5,i2)*B(i2,i3)
                   -B(i3,i4)*A(i4,i2)*B(i2,i1))
              *L0(-S(i3,i4),-S(i5,i1))/S(i5,i1)
           /B(i1,i2)/B(i2,i3)/B(i3,i4)/B(i4,i5)/B(i5,i1)
             +1.0/3.0*A(i2,i4)*A(i2,i4)*A(i2,i5)*A(i2,i5)
           /A(i1,i2)/A(i2,i3)/A(i3,i4)/B(i4,i5)/A(i5,i1)
             -1.0/3.0*B(i1,i2)*B(i4,i1)*B(i4,i1)
                     *A(i2,i4)*A(i2,i4)*A(i2,i4)
                   /B(i4,i5)/B(i5,i1)/B(i2,i4)
                      /A(i2,i3)/A(i3,i4)/S(i5,i1)
             +1.0/3.0*B(i3,i2)*B(i5,i3)*B(i5,i3)
                     *A(i2,i5)*A(i2,i5)*A(i2,i5)
                   /B(i5,i4)/B(i4,i3)/B(i2,i5)
                      /A(i2,i1)/A(i1,i5)/S(i3,i4)
             +1.0/6.0*B(i1,i3)*B(i1,i3)*A(i2,i4)*A(i2,i5)
                      /S(i3,i4)/B(i4,i5)/S(i5,i1);

return ( (Vf+Vs)*tree +Ff+Fs );
}

*/


//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
//---------- Fixed order LO for signal and interference  LEAN ---------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Leandro Cieri /////////////////////////////////////////////////////

  std::complex<double> ampq2g1p2::LO_MggH(int p1, int p2, int p3,int p4,int p5) const
  { //Este elemento de matriz calcula la produccion de Higgs gg->H
	// Falta hacer el modulo cuadrado que se hace en el momento de multiplicar por las PDFs
	// In the infinit top mass limit
	// Por el momento solo la contribucion top
	// alpha_s se mete en la parte de las PDFs cuando se hace el cuadrado

    double s12 = S(1,2);


    cout << " ampq2g1p2::LO_MggH " << endl;

      double vmartin = _VEV_/sqrt(2.0); // Ok
      //				 s12 = 16094.8;
     //     vmartin=174.0;
     //   s12=(124.95)*(124.95);

      std::complex<double>   Aggh =-s12/(8.0*sqrt(2.0)*__DEB_PI  *vmartin)*mesq::SumaMartin(s12);

   // Checked!!!!
   //  std::cout << "Aggh  :" << Aggh  << std::endl;
   //   std::cout << "asin(sqrt(s12/(mt*mt))/2.0)  :"<< asin(sqrt(s12/(mt*mt))/2.0)  << std::endl;
   //   std::cout << "s12  :"<< sqrt(S(1,2)) << std::endl;

  //  std::cout << "Aggh  :" << Aggh  << std::endl;

    return Aggh ;

 }



  std::complex<double> ampq2g1p2::LO_MHyy(int p1, int p2, int p3,int p4,int p5) const
  { //Este elemento de matriz calcula el decaimiento del Higgs H->yy
	// Falta hacer el modulo cuadrado que se hace en el momento de multiplicar por las PDFs
	// In the infinit top mass limit
	// Por el momento solo la contribucion top, tau y W

    double s12 = S(1,2);

    cout << " ampq2g1p2::LO_MHyy " << endl;
    //				 s12 = 16094.8;

      double vmartin = _VEV_/sqrt(2.0);

   //    vmartin=174.0;
      //      cout << " alpha_qed : " << alpha_qed << endl;
      std::complex<double>   Ahyy =  -(alpha_qed*s12)/(4.0*sqrt(2.0)*__DEB_PI*vmartin)*mesq::CorcheteMartin(s12);
 // Checked!!!!
    return Ahyy ;

 }


 double ampq2g1p2::LO_Mggyy(int p1, int p2, int p3,int p4,int p5,double & spass) const
  { //Este elemento de matriz calcula solo la contribuciones de helicidad que entran en juego en la interferencia
	// La contribucion del top no está tenida en cuenta
	// Suma sobre el resto de las cargas de los quarks livianos

    double s12 = S(1,2), t = S(1,5), u = S(1,4);
			

    cout << " ampq2g1p2::LO_Mggyy " << endl;
			 /*
 s12 = 16094.8;
			 t = -10185.6;
				 u = -5909.13;



     // Checked!!!!

			 */
   //   s12=(124.95)*(124.95);
   //   u=-10000.0;
   //   t=-u-s12;
    //  double z = (t-u)/s12;
   //   std::cout << "z  :" << z << std::endl;

    spass = s12;
	

  double Aggyy =  2.0*alpha_qed*(-(((u-t)*log(u/t))/s12)
             -((t*t + u*u)*(__DEB_PI*__DEB_PI + log(u/t)*log(u/t)))/(2.0*s12*s12))*
             11.0/9.0; //alpha_s to be included at final moment
  
  //  std::cout << "Aggyy  :" << Aggyy << std::endl; 
 /*
   std::cout << "Aggyy  :" << Aggyy << std::endl;
     std::cout << "S(1,5)  :" << S(1,5) << std::endl;
      std::cout << "S(1,4)  :" << S(1,4) << std::endl;
      std::cout << "S(1,2)  :" << S(1,2) << std::endl;
     */

      // Checked!!!!
     return  Aggyy ;
 }


 
 
 
 
 
 
 
 
 
 
 
 void Lggg_background(std::vector<std::complex<double> >& helamps) 
{ 



//std::cout << " S12  " << S[1][3]  << std::endl;

//  (1+2+3+4+5+)   (Checked out versus Zvi.)

helamps[0] = Mppppp(1,2,3,4,5)+Mppppp(1,4,5,2,3)+Mppppp(1,2,4,5,3)
   +Mppppp(1,2,3,5,4)+Mppppp(1,5,4,2,3)+Mppppp(1,2,5,4,3)
   +Mppppp(1,2,4,3,5)+Mppppp(1,5,2,3,4)+Mppppp(1,4,2,5,3)
   +Mppppp(1,2,5,3,4)+Mppppp(1,4,2,3,5)+Mppppp(1,5,2,4,3);

//Msq+=norm(M);

//  (1-2+3+4+5+)   (Checked out versus Zvi.)

helamps[1]  = Mmpppp(1,2,3,4,5)+Mmpppp(1,4,5,2,3)+Mmpppp(1,2,4,5,3)
   +Mmpppp(1,2,3,5,4)+Mmpppp(1,5,4,2,3)+Mmpppp(1,2,5,4,3)
   +Mmpppp(1,2,4,3,5)+Mmpppp(1,5,2,3,4)+Mmpppp(1,4,2,5,3)
   +Mmpppp(1,2,5,3,4)+Mmpppp(1,4,2,3,5)+Mmpppp(1,5,2,4,3);

//Msq+=norm(M);

//  (1+2-3+4+5+)   (Checked out versus Zvi.)

helamps[2]  = Mmpppp(2,3,4,5,1)+Mmpppp(2,3,1,4,5)+Mmpppp(2,4,5,3,1)
   +Mmpppp(2,3,5,4,1)+Mmpppp(2,3,1,5,4)+Mmpppp(2,5,4,3,1)
   +Mmpppp(2,4,3,5,1)+Mmpppp(2,3,4,1,5)+Mmpppp(2,5,3,1,4)
   +Mmpppp(2,5,3,4,1)+Mmpppp(2,3,5,1,4)+Mmpppp(2,4,3,1,5);

// Msq+=norm(M);

//  (1+2+3-4+5+)   (Checked out versus Zvi.)

helamps[3]  = Mmpppp(3,4,5,1,2)+Mmpppp(3,1,4,5,2)+Mmpppp(3,1,2,4,5)
   +Mmpppp(3,5,4,1,2)+Mmpppp(3,1,5,4,2)+Mmpppp(3,1,2,5,4)
   +Mmpppp(3,5,1,2,4)+Mmpppp(3,4,1,5,2)+Mmpppp(3,1,4,2,5)
   +Mmpppp(3,4,1,2,5)+Mmpppp(3,5,1,4,2)+Mmpppp(3,1,5,2,4);

//Msq+=norm(M);

//  (1+2+3+4-5+)   (Checked out versus Zvi.)
// 5
helamps[4]  = Mmpppp(4,5,1,2,3)+Mmpppp(4,5,2,3,1)+Mmpppp(4,5,3,1,2)
   +Mmpppp(4,1,2,3,5)+Mmpppp(4,2,3,1,5)+Mmpppp(4,3,1,2,5)
   +Mmpppp(4,3,5,1,2)+Mmpppp(4,1,5,2,3)+Mmpppp(4,2,5,3,1)
   +Mmpppp(4,1,2,5,3)+Mmpppp(4,2,3,5,1)+Mmpppp(4,3,1,5,2);

//Msq+=norm(M);

//  (1+2+3+4+5-)   (Checked out versus Zvi.)

helamps[5]  = Mmpppp(5,1,2,3,4)+Mmpppp(5,2,3,1,4)+Mmpppp(5,3,1,2,4)
   +Mmpppp(5,4,1,2,3)+Mmpppp(5,4,2,3,1)+Mmpppp(5,4,3,1,2)
   +Mmpppp(5,1,2,4,3)+Mmpppp(5,2,3,4,1)+Mmpppp(5,3,1,4,2)
   +Mmpppp(5,3,4,1,2)+Mmpppp(5,1,4,2,3)+Mmpppp(5,2,4,3,1);

//Msq+=norm(M);

//  (1-2-3+4+5+)   (Checked out versus Zvi.)

helamps[6]  = Mmmppp(1,2,3,4,5)+Mmpmpp(2,3,1,4,5)+Mmmppp(1,2,4,5,3)
   +Mmmppp(1,2,3,5,4)+Mmpmpp(2,3,1,5,4)+Mmmppp(1,2,5,4,3)
   +Mmmppp(1,2,4,3,5)+Mmpmpp(1,5,2,3,4)+Mmpmpp(1,4,2,5,3)
   +Mmmppp(1,2,5,3,4)+Mmpmpp(1,4,2,3,5)+Mmpmpp(1,5,2,4,3);

//Msq+=norm(M);

//  (1-2+3-4+5+)   (Checked out versus Zvi.)

helamps[7]  = Mmpmpp(1,2,3,4,5)+Mmmppp(3,1,4,5,2)+Mmmppp(3,1,2,4,5)
   +Mmpmpp(1,2,3,5,4)+Mmmppp(3,1,5,4,2)+Mmmppp(3,1,2,5,4)
   +Mmpmpp(3,5,1,2,4)+Mmpmpp(3,4,1,5,2)+Mmmppp(3,1,4,2,5)
   +Mmpmpp(3,4,1,2,5)+Mmpmpp(3,5,1,4,2)+Mmmppp(3,1,5,2,4);

//Msq+=norm(M);

//  (1-2+3+4-5+)   (Checked out versus Zvi.)

helamps[8]  = Mmpmpp(4,5,1,2,3)+Mmmppp(1,4,5,2,3)+Mmpmpp(1,2,4,5,3)
   +Mmmppp(4,1,2,3,5)+Mmpmpp(1,5,4,2,3)+Mmpmpp(4,3,1,2,5)
   +Mmpmpp(1,2,4,3,5)+Mmmppp(4,1,5,2,3)+Mmmppp(1,4,2,5,3)
   +Mmmppp(4,1,2,5,3)+Mmmppp(1,4,2,3,5)+Mmpmpp(4,3,1,5,2);

//Msq+=norm(M);

//  (1-2+3+4+5-)   (Checked out versus Zvi.)

helamps[9]  = Mmmppp(5,1,2,3,4)+Mmpmpp(1,4,5,2,3)+Mmpmpp(5,3,1,2,4)
   +Mmpmpp(5,4,1,2,3)+Mmmppp(1,5,4,2,3)+Mmpmpp(1,2,5,4,3)
   +Mmmppp(5,1,2,4,3)+Mmmppp(1,5,2,3,4)+Mmpmpp(5,3,1,4,2)
   +Mmpmpp(1,2,5,3,4)+Mmmppp(5,1,4,2,3)+Mmmppp(1,5,2,4,3);

//Msq+=norm(M);

//  (1+2-3-4+5+)   (Checked out versus Zvi.)
// 11
helamps[10]  = Mmmppp(2,3,4,5,1)+Mmmppp(2,3,1,4,5)+Mmpmpp(3,1,2,4,5)
   +Mmmppp(2,3,5,4,1)+Mmmppp(2,3,1,5,4)+Mmpmpp(3,1,2,5,4)
   +Mmpmpp(2,4,3,5,1)+Mmmppp(2,3,4,1,5)+Mmpmpp(2,5,3,1,4)
   +Mmpmpp(2,5,3,4,1)+Mmmppp(2,3,5,1,4)+Mmpmpp(2,4,3,1,5);

//Msq+=norm(M);

//  (1+2-3+4-5+)   (Checked out versus Zvi.)

helamps[11]  = Mmpmpp(2,3,4,5,1)+Mmpmpp(4,5,2,3,1)+Mmmppp(2,4,5,3,1)
   +Mmpmpp(4,1,2,3,5)+Mmmppp(4,2,3,1,5)+Mmpmpp(2,5,4,3,1)
   +Mmmppp(2,4,3,5,1)+Mmpmpp(2,3,4,1,5)+Mmmppp(4,2,5,3,1)
   +Mmpmpp(4,1,2,5,3)+Mmmppp(4,2,3,5,1)+Mmmppp(2,4,3,1,5);

//Msq+=norm(M);

//  (1+2-3+4+5-)   (Checked out versus Zvi.)

helamps[12]  = Mmpmpp(5,1,2,3,4)+Mmmppp(5,2,3,1,4)+Mmpmpp(2,4,5,3,1)
   +Mmpmpp(2,3,5,4,1)+Mmpmpp(5,4,2,3,1)+Mmmppp(2,5,4,3,1)
   +Mmpmpp(5,1,2,4,3)+Mmmppp(5,2,3,4,1)+Mmmppp(2,5,3,1,4)
   +Mmmppp(2,5,3,4,1)+Mmpmpp(2,3,5,1,4)+Mmmppp(5,2,4,3,1);

//Msq+=norm(M);

//  (1+2+3-4-5+)   (Checked out versus Zvi.)

helamps[13]  = Mmmppp(3,4,5,1,2)+Mmpmpp(3,1,4,5,2)+Mmpmpp(4,5,3,1,2)
   +Mmpmpp(3,5,4,1,2)+Mmpmpp(4,2,3,1,5)+Mmmppp(4,3,1,2,5)
   +Mmmppp(4,3,5,1,2)+Mmmppp(3,4,1,5,2)+Mmpmpp(3,1,4,2,5)
   +Mmmppp(3,4,1,2,5)+Mmpmpp(4,2,3,5,1)+Mmmppp(4,3,1,5,2);

//Msq+=norm(M);

//  (1+2+3-4+5-)   (Checked out versus Zvi.)

helamps[14]  = Mmpmpp(3,4,5,1,2)+Mmpmpp(5,2,3,1,4)+Mmmppp(5,3,1,2,4)
   +Mmmppp(3,5,4,1,2)+Mmpmpp(3,1,5,4,2)+Mmpmpp(5,4,3,1,2)
   +Mmmppp(3,5,1,2,4)+Mmpmpp(5,2,3,4,1)+Mmmppp(5,3,1,4,2)
   +Mmmppp(5,3,4,1,2)+Mmmppp(3,5,1,4,2)+Mmpmpp(3,1,5,2,4);

//Msq+=norm(M);

//  (1+2+3+4-5-)   (Checked out versus Zvi.)

helamps[15]  = Mmmppp(4,5,1,2,3)+Mmmppp(4,5,2,3,1)+Mmmppp(4,5,3,1,2)
   +Mmmppp(5,4,1,2,3)+Mmmppp(5,4,2,3,1)+Mmmppp(5,4,3,1,2)
   +Mmpmpp(4,3,5,1,2)+Mmpmpp(4,1,5,2,3)+Mmpmpp(4,2,5,3,1)
   +Mmpmpp(5,3,4,1,2)+Mmpmpp(5,1,4,2,3)+Mmpmpp(5,2,4,3,1);

//Msq+=norm(M);

//  (1-2-3-4-5-)   (Checked out versus Zvi.)

helamps[16]  = Mmmmmm(1,2,3,4,5)+Mmmmmm(1,4,5,2,3)+Mmmmmm(1,2,4,5,3)
   +Mmmmmm(1,2,3,5,4)+Mmmmmm(1,5,4,2,3)+Mmmmmm(1,2,5,4,3)
   +Mmmmmm(1,2,4,3,5)+Mmmmmm(1,5,2,3,4)+Mmmmmm(1,4,2,5,3)
   +Mmmmmm(1,2,5,3,4)+Mmmmmm(1,4,2,3,5)+Mmmmmm(1,5,2,4,3);

//Msq+=norm(M);

//  (1+2-3-4-5-)   (Checked out versus Zvi.)

helamps[17]  = Mpmmmm(1,2,3,4,5)+Mpmmmm(1,4,5,2,3)+Mpmmmm(1,2,4,5,3)
   +Mpmmmm(1,2,3,5,4)+Mpmmmm(1,5,4,2,3)+Mpmmmm(1,2,5,4,3)
   +Mpmmmm(1,2,4,3,5)+Mpmmmm(1,5,2,3,4)+Mpmmmm(1,4,2,5,3)
   +Mpmmmm(1,2,5,3,4)+Mpmmmm(1,4,2,3,5)+Mpmmmm(1,5,2,4,3);

//Msq+=norm(M);

//  (1-2+3-4-5-)   (Checked out versus Zvi.)

helamps[18]  = Mpmmmm(2,3,4,5,1)+Mpmmmm(2,3,1,4,5)+Mpmmmm(2,4,5,3,1)
   +Mpmmmm(2,3,5,4,1)+Mpmmmm(2,3,1,5,4)+Mpmmmm(2,5,4,3,1)
   +Mpmmmm(2,4,3,5,1)+Mpmmmm(2,3,4,1,5)+Mpmmmm(2,5,3,1,4)
   +Mpmmmm(2,5,3,4,1)+Mpmmmm(2,3,5,1,4)+Mpmmmm(2,4,3,1,5);

//Msq+=norm(M);

//  (1-2-3+4-5-)   (Checked out versus Zvi.)

helamps[19]  = Mpmmmm(3,4,5,1,2)+Mpmmmm(3,1,4,5,2)+Mpmmmm(3,1,2,4,5)
   +Mpmmmm(3,5,4,1,2)+Mpmmmm(3,1,5,4,2)+Mpmmmm(3,1,2,5,4)
   +Mpmmmm(3,5,1,2,4)+Mpmmmm(3,4,1,5,2)+Mpmmmm(3,1,4,2,5)
   +Mpmmmm(3,4,1,2,5)+Mpmmmm(3,5,1,4,2)+Mpmmmm(3,1,5,2,4);

//Msq+=norm(M);

//  (1-2-3-4+5-)   (Checked out versus Zvi.)

helamps[20]  = Mpmmmm(4,5,1,2,3)+Mpmmmm(4,5,2,3,1)+Mpmmmm(4,5,3,1,2)
   +Mpmmmm(4,1,2,3,5)+Mpmmmm(4,2,3,1,5)+Mpmmmm(4,3,1,2,5)
   +Mpmmmm(4,3,5,1,2)+Mpmmmm(4,1,5,2,3)+Mpmmmm(4,2,5,3,1)
   +Mpmmmm(4,1,2,5,3)+Mpmmmm(4,2,3,5,1)+Mpmmmm(4,3,1,5,2);

//Msq+=norm(M);

//  (1-2-3-4-5+)   (Checked out versus Zvi.)

helamps[21]  = Mpmmmm(5,1,2,3,4)+Mpmmmm(5,2,3,1,4)+Mpmmmm(5,3,1,2,4)
   +Mpmmmm(5,4,1,2,3)+Mpmmmm(5,4,2,3,1)+Mpmmmm(5,4,3,1,2)
   +Mpmmmm(5,1,2,4,3)+Mpmmmm(5,2,3,4,1)+Mpmmmm(5,3,1,4,2)
   +Mpmmmm(5,3,4,1,2)+Mpmmmm(5,1,4,2,3)+Mpmmmm(5,2,4,3,1);

//Msq+=norm(M);

//  (1+2+3-4-5-)   (Checked out versus Zvi.)

helamps[22]  = Mppmmm(1,2,3,4,5)+Mpmpmm(2,3,1,4,5)+Mppmmm(1,2,4,5,3)
   +Mppmmm(1,2,3,5,4)+Mpmpmm(2,3,1,5,4)+Mppmmm(1,2,5,4,3)
   +Mppmmm(1,2,4,3,5)+Mpmpmm(1,5,2,3,4)+Mpmpmm(1,4,2,5,3)
   +Mppmmm(1,2,5,3,4)+Mpmpmm(1,4,2,3,5)+Mpmpmm(1,5,2,4,3);

//Msq+=norm(M);

//  (1+2-3+4-5-)   (Checked out versus Zvi.)

helamps[23]  = Mpmpmm(1,2,3,4,5)+Mppmmm(3,1,4,5,2)+Mppmmm(3,1,2,4,5)
   +Mpmpmm(1,2,3,5,4)+Mppmmm(3,1,5,4,2)+Mppmmm(3,1,2,5,4)
   +Mpmpmm(3,5,1,2,4)+Mpmpmm(3,4,1,5,2)+Mppmmm(3,1,4,2,5)
   +Mpmpmm(3,4,1,2,5)+Mpmpmm(3,5,1,4,2)+Mppmmm(3,1,5,2,4);

//Msq+=norm(M);

//  (1+2-3-4+5-)   (Checked out versus Zvi.)

helamps[24]  = Mpmpmm(4,5,1,2,3)+Mppmmm(1,4,5,2,3)+Mpmpmm(1,2,4,5,3)
   +Mppmmm(4,1,2,3,5)+Mpmpmm(1,5,4,2,3)+Mpmpmm(4,3,1,2,5)
   +Mpmpmm(1,2,4,3,5)+Mppmmm(4,1,5,2,3)+Mppmmm(1,4,2,5,3)
   +Mppmmm(4,1,2,5,3)+Mppmmm(1,4,2,3,5)+Mpmpmm(4,3,1,5,2);

//Msq+=norm(M);

//  (1+2-3-4-5+)   (Checked out versus Zvi.)

helamps[25]  = Mppmmm(5,1,2,3,4)+Mpmpmm(1,4,5,2,3)+Mpmpmm(5,3,1,2,4)
   +Mpmpmm(5,4,1,2,3)+Mppmmm(1,5,4,2,3)+Mpmpmm(1,2,5,4,3)
   +Mppmmm(5,1,2,4,3)+Mppmmm(1,5,2,3,4)+Mpmpmm(5,3,1,4,2)
   +Mpmpmm(1,2,5,3,4)+Mppmmm(5,1,4,2,3)+Mppmmm(1,5,2,4,3);

//Msq+=norm(M);

//  (1-2+3+4-5-)   (Checked out versus Zvi.)

helamps[26]  = Mppmmm(2,3,4,5,1)+Mppmmm(2,3,1,4,5)+Mpmpmm(3,1,2,4,5)
   +Mppmmm(2,3,5,4,1)+Mppmmm(2,3,1,5,4)+Mpmpmm(3,1,2,5,4)
   +Mpmpmm(2,4,3,5,1)+Mppmmm(2,3,4,1,5)+Mpmpmm(2,5,3,1,4)
   +Mpmpmm(2,5,3,4,1)+Mppmmm(2,3,5,1,4)+Mpmpmm(2,4,3,1,5);

//Msq+=norm(M);

//  (1-2+3-4+5-)   (Checked out versus Zvi.)

helamps[27]  = Mpmpmm(2,3,4,5,1)+Mpmpmm(4,5,2,3,1)+Mppmmm(2,4,5,3,1)
   +Mpmpmm(4,1,2,3,5)+Mppmmm(4,2,3,1,5)+Mpmpmm(2,5,4,3,1)
   +Mppmmm(2,4,3,5,1)+Mpmpmm(2,3,4,1,5)+Mppmmm(4,2,5,3,1)
   +Mpmpmm(4,1,2,5,3)+Mppmmm(4,2,3,5,1)+Mppmmm(2,4,3,1,5);

//Msq+=norm(M);

//  (1-2+3-4-5+)   (Checked out versus Zvi.)

helamps[28]  = Mpmpmm(5,1,2,3,4)+Mppmmm(5,2,3,1,4)+Mpmpmm(2,4,5,3,1)
   +Mpmpmm(2,3,5,4,1)+Mpmpmm(5,4,2,3,1)+Mppmmm(2,5,4,3,1)
   +Mpmpmm(5,1,2,4,3)+Mppmmm(5,2,3,4,1)+Mppmmm(2,5,3,1,4)
   +Mppmmm(2,5,3,4,1)+Mpmpmm(2,3,5,1,4)+Mppmmm(5,2,4,3,1);

//Msq+=norm(M);

//  (1-2-3+4+5-)   (Checked out versus Zvi.)

helamps[29]  = Mppmmm(3,4,5,1,2)+Mpmpmm(3,1,4,5,2)+Mpmpmm(4,5,3,1,2)
   +Mpmpmm(3,5,4,1,2)+Mpmpmm(4,2,3,1,5)+Mppmmm(4,3,1,2,5)
   +Mppmmm(4,3,5,1,2)+Mppmmm(3,4,1,5,2)+Mpmpmm(3,1,4,2,5)
   +Mppmmm(3,4,1,2,5)+Mpmpmm(4,2,3,5,1)+Mppmmm(4,3,1,5,2);

//helamps[0] Msq+=norm(M);

//  (1-2-3+4-5+)   (Checked out versus Zvi.)

helamps[30]  = Mpmpmm(3,4,5,1,2)+Mpmpmm(5,2,3,1,4)+Mppmmm(5,3,1,2,4)
   +Mppmmm(3,5,4,1,2)+Mpmpmm(3,1,5,4,2)+Mpmpmm(5,4,3,1,2)
   +Mppmmm(3,5,1,2,4)+Mpmpmm(5,2,3,4,1)+Mppmmm(5,3,1,4,2)
   +Mppmmm(5,3,4,1,2)+Mppmmm(3,5,1,4,2)+Mpmpmm(3,1,5,2,4);

//Msq+=norm(M);

//  (1-2-3-4+5+)   (Checked out versus Zvi.)
// 32
helamps[31]  = Mppmmm(4,5,1,2,3)+Mppmmm(4,5,2,3,1)+Mppmmm(4,5,3,1,2)
   +Mppmmm(5,4,1,2,3)+Mppmmm(5,4,2,3,1)+Mppmmm(5,4,3,1,2)
   +Mpmpmm(4,3,5,1,2)+Mpmpmm(4,1,5,2,3)+Mpmpmm(4,2,5,3,1)
   +Mpmpmm(5,3,4,1,2)+Mpmpmm(5,1,4,2,3)+Mpmpmm(5,2,4,3,1);

//Msq+=norm(M);

// Matrix elements are multiplied by (4*pi)^2*sqrt(shat).

 // double pcrossGG=Nc*zm*z*Msq/4.0;

 // return sigma0*as2pi*pcrossGG*fgg;
}




complex<double> Mppppp(int i1, int i2, int i3, int i4, int i5)
{
return -(LS[i1][i2]*LS[i2][i3]+LS[i2][i3]*LS[i3][i4]+LS[i3][i4]*LS[i4][i5]
       +LS[i4][i5]*LS[i5][i1]+LS[i5][i1]*LS[i1][i2]
       +LSpb[i1][i2]*LSpa[i2][i3]*LSpb[i3][i4]*LSpa[i4][i1]
       -LSpa[i1][i2]*LSpb[i2][i3]*LSpa[i3][i4]*LSpb[i4][i1])
        /LSpa[i1][i2]/LSpa[i2][i3]/LSpa[i3][i4]/LSpa[i4][i5]/LSpa[i5][i1]/6.0;
}

complex<double> Mmmmmm(int i1, int i2, int i3, int i4, int i5)
{
return (LS[i1][i2]*LS[i2][i3]+LS[i2][i3]*LS[i3][i4]+LS[i3][i4]*LS[i4][i5]
       +LS[i4][i5]*LS[i5][i1]+LS[i5][i1]*LS[i1][i2]
       +LSpa[i1][i2]*LSpb[i2][i3]*LSpa[i3][i4]*LSpb[i4][i1]
       -LSpb[i1][i2]*LSpa[i2][i3]*LSpb[i3][i4]*LSpa[i4][i1])
        /LSpb[i1][i2]/LSpb[i2][i3]/LSpb[i3][i4]/LSpb[i4][i5]/LSpb[i5][i1]/6.0;
}

complex<double> Mmpppp(int i1, int i2, int i3, int i4, int i5)
{
return -((LS[i2][i3]+LS[i3][i4]+LS[i4][i5])*LSpb[i2][i5]*LSpb[i2][i5]
        -LSpb[i2][i4]*LSpa[i4][i3]*LSpb[i3][i5]*LSpb[i2][i5]
        -LSpb[i1][i2]*LSpb[i1][i5]/LSpa[i1][i2]/LSpa[i1][i5]*(
            LSpa[i1][i2]*LSpa[i1][i2]*LSpa[i1][i3]*LSpa[i1][i3]
                     *LSpb[i2][i3]/LSpa[i2][i3]
           +LSpa[i1][i3]*LSpa[i1][i3]*LSpa[i1][i4]*LSpa[i1][i4]
                     *LSpb[i3][i4]/LSpa[i3][i4]
           +LSpa[i1][i4]*LSpa[i1][i4]*LSpa[i1][i5]*LSpa[i1][i5]
                     *LSpb[i4][i5]/LSpa[i4][i5]))
        /LSpb[i1][i2]/LSpa[i2][i3]/LSpa[i3][i4]/LSpa[i4][i5]/LSpb[i5][i1]/3.0;
}

complex<double> Mpmmmm(int i1, int i2, int i3, int i4, int i5)
{
return ((LS[i2][i3]+LS[i3][i4]+LS[i4][i5])*LSpa[i2][i5]*LSpa[i2][i5]
        -LSpa[i2][i4]*LSpb[i4][i3]*LSpa[i3][i5]*LSpa[i2][i5]
        -LSpa[i1][i2]*LSpa[i1][i5]/LSpb[i1][i2]/LSpb[i1][i5]*(
            LSpb[i1][i2]*LSpb[i1][i2]*LSpb[i1][i3]*LSpb[i1][i3]
                     *LSpa[i2][i3]/LSpb[i2][i3]
           +LSpb[i1][i3]*LSpb[i1][i3]*LSpb[i1][i4]*LSpb[i1][i4]
                     *LSpa[i3][i4]/LSpb[i3][i4]
           +LSpb[i1][i4]*LSpb[i1][i4]*LSpb[i1][i5]*LSpb[i1][i5]
                     *LSpa[i4][i5]/LSpb[i4][i5]))
        /LSpa[i1][i2]/LSpb[i2][i3]/LSpb[i3][i4]/LSpb[i4][i5]/LSpa[i5][i1]/3.0;
}

complex<double> Mmmppp(int i1, int i2, int i3, int i4, int i5)
{

  //  Note: nothing should depend on mm.

  double mm=1.0;

  complex<double> l23,l51;

  if (LS[i2][i3]<0) {
     l23=log(mm/(-LS[i2][i3]));
  } else {
     l23=log(mm/LS[i2][i3])+I*PI;
  }

  if (LS[i5][i1]<0) {
     l51=log(mm/(-LS[i5][i1]));
  } else {
     l51=log(mm/LS[i5][i1])+I*PI;
  }


  complex<double> Vf=-0.5*(l23+l51);
  complex<double> Vs=-1.0/3.0*Vf;

  complex<double> tree=pow(LSpa[i1][i2],4)/LSpa[i1][i2]/LSpa[i2][i3]/LSpa[i3][i4]
                      /LSpa[i4][i5]/LSpa[i5][i1];

  complex<double> Ff=-0.5*LSpa[i1][i2]*LSpa[i1][i2]*(
                    LSpa[i2][i3]*LSpb[i3][i4]*LSpa[i4][i1]
                   +LSpa[i2][i4]*LSpb[i4][i5]*LSpa[i5][i1])
              *LL0(-LS[i2][i3],-LS[i5][i1])/LS[i5][i1]
           /LSpa[i2][i3]/LSpa[i3][i4]/LSpa[i4][i5]/LSpa[i5][i1];

  complex<double> Fs=-1.0/3.0*LSpb[i3][i4]*LSpa[i4][i1]*LSpa[i2][i4]*LSpb[i4][i5]*(
                    LSpa[i2][i3]*LSpb[i3][i4]*LSpa[i4][i1]
                   +LSpa[i2][i4]*LSpb[i4][i5]*LSpa[i5][i1])
              *LL2(-LS[i2][i3],-LS[i5][i1])/LS[i5][i1]/LS[i5][i1]/LS[i5][i1]
           /LSpa[i3][i4]/LSpa[i4][i5]
             -1.0/3.0*Ff
             -1.0/3.0*LSpa[i3][i5]*LSpb[i3][i5]*LSpb[i3][i5]*LSpb[i3][i5]
                   /LSpb[i1][i2]/LSpb[i2][i3]/LSpa[i3][i4]
                      /LSpa[i4][i5]/LSpb[i5][i1]
             +1.0/3.0*LSpa[i1][i2]*LSpb[i3][i5]*LSpb[i3][i5]
                   /LSpb[i2][i3]/LSpa[i3][i4]
                      /LSpa[i4][i5]/LSpb[i5][i1]
             +1.0/6.0*LSpa[i1][i2]*LSpb[i3][i4]
                      *LSpa[i4][i1]*LSpa[i2][i4]*LSpb[i4][i5]
                    /LS[i2][i3]/LSpa[i3][i4]/LSpa[i4][i5]/LS[i5][i1];

return -( (Vf+Vs)*tree +Ff+Fs );
}
complex<double> Mppmmm(int i1, int i2, int i3, int i4, int i5)
{

  //  Note: nothing should depend on mm.

  double mm=1.0;

  complex<double> l23,l51;

  if (LS[i2][i3]<0) {
     l23=log(mm/(-LS[i2][i3]));
  } else {
     l23=log(mm/LS[i2][i3])+I*PI;
  }

  if (LS[i5][i1]<0) {
     l51=log(mm/(-LS[i5][i1]));
  } else {
     l51=log(mm/LS[i5][i1])+I*PI;
  }


  complex<double> Vf=-0.5*(l23+l51);
  complex<double> Vs=-1.0/3.0*Vf;

  complex<double> tree=pow(LSpb[i1][i2],4)/LSpb[i1][i2]/LSpb[i2][i3]/LSpb[i3][i4]
                      /LSpb[i4][i5]/LSpb[i5][i1];

  complex<double> Ff=-0.5*LSpb[i1][i2]*LSpb[i1][i2]*(
                    LSpb[i2][i3]*LSpa[i3][i4]*LSpb[i4][i1]
                   +LSpb[i2][i4]*LSpa[i4][i5]*LSpb[i5][i1])
              *LL0(-LS[i2][i3],-LS[i5][i1])/LS[i5][i1]
           /LSpb[i2][i3]/LSpb[i3][i4]/LSpb[i4][i5]/LSpb[i5][i1];

  complex<double> Fs=-1.0/3.0*LSpa[i3][i4]*LSpb[i4][i1]*LSpb[i2][i4]*LSpa[i4][i5]*(
                    LSpb[i2][i3]*LSpa[i3][i4]*LSpb[i4][i1]
                   +LSpb[i2][i4]*LSpa[i4][i5]*LSpb[i5][i1])
              *LL2(-LS[i2][i3],-LS[i5][i1])/LS[i5][i1]/LS[i5][i1]/LS[i5][i1]
           /LSpb[i3][i4]/LSpb[i4][i5]
             -1.0/3.0*Ff
             -1.0/3.0*LSpb[i3][i5]*LSpa[i3][i5]*LSpa[i3][i5]*LSpa[i3][i5]
                   /LSpa[i1][i2]/LSpa[i2][i3]/LSpb[i3][i4]
                      /LSpb[i4][i5]/LSpa[i5][i1]
             +1.0/3.0*LSpb[i1][i2]*LSpa[i3][i5]*LSpa[i3][i5]
                   /LSpa[i2][i3]/LSpb[i3][i4]
                      /LSpb[i4][i5]/LSpa[i5][i1]
             +1.0/6.0*LSpb[i1][i2]*LSpa[i3][i4]
                      *LSpb[i4][i1]*LSpb[i2][i4]*LSpa[i4][i5]
                    /LS[i2][i3]/LSpb[i3][i4]/LSpb[i4][i5]/LS[i5][i1];

return ( (Vf+Vs)*tree +Ff+Fs );
}

complex<double> Mmpmpp(int i1, int i2, int i3, int i4, int i5)
{

  //  Note: nothing should depend on mm.

  double mm=1.0;

  complex<double> l34,l51;

  if (LS[i3][i4]<0) {
     l34=log(mm/(-LS[i3][i4]));
  } else {
     l34=log(mm/LS[i3][i4])+I*PI;
  }

  if (LS[i5][i1]<0) {
     l51=log(mm/(-LS[i5][i1]));
  } else {
     l51=log(mm/LS[i5][i1])+I*PI;
  }

  complex<double> Vf=-0.5*(l34+l51);
  complex<double> Vs=-1.0/3.0*Vf;

  complex<double> tree=pow(LSpa[i1][i3],4)/LSpa[i1][i2]/LSpa[i2][i3]/LSpa[i3][i4]
                      /LSpa[i4][i5]/LSpa[i5][i1];

  complex<double> Ff=-LSpa[i1][i3]*LSpa[i1][i3]*LSpa[i4][i1]*LSpb[i2][i4]*LSpb[i2][i4]
                  *LLs1(-LS[i2][i3],-LS[i5][i1],-LS[i3][i4],-LS[i5][i1])
                   /LSpa[i4][i5]/LSpa[i5][i1]/LS[i5][i1]/LS[i5][i1]
             +LSpa[i1][i3]*LSpa[i1][i3]*LSpa[i5][i3]*LSpb[i2][i5]*LSpb[i2][i5]
                  *LLs1(-LS[i1][i2],-LS[i3][i4],-LS[i5][i1],-LS[i3][i4])
                   /LSpa[i3][i4]/LSpa[i4][i5]/LS[i3][i4]/LS[i3][i4]
             -0.5*LSpa[i1][i3]*LSpa[i1][i3]*LSpa[i1][i3]*(
                    LSpa[i1][i5]*LSpb[i5][i2]*LSpa[i2][i3]
                   -LSpa[i3][i4]*LSpb[i4][i2]*LSpa[i2][i1])
              *LL0(-LS[i3][i4],-LS[i5][i1])/LS[i5][i1]
           /LSpa[i1][i2]/LSpa[i2][i3]/LSpa[i3][i4]/LSpa[i4][i5]/LSpa[i5][i1];

  complex<double> Fs=-LSpa[i1][i2]*LSpa[i2][i3]*LSpa[i3][i4]*LSpa[i4][i1]
                      *LSpa[i4][i1]*LSpb[i2][i4]*LSpb[i2][i4]*(
                2.0*LLs1(-LS[i2][i3],-LS[i5][i1],-LS[i3][i4],-LS[i5][i1])
                   +LL1(-LS[i2][i3],-LS[i5][i1])
                   +LL1(-LS[i3][i4],-LS[i5][i1]))
                  /LSpa[i4][i5]/LSpa[i5][i1]/LSpa[i2][i4]/LSpa[i2][i4]
                  /LS[i5][i1]/LS[i5][i1]
             +LSpa[i3][i2]*LSpa[i2][i1]*LSpa[i1][i5]*LSpa[i5][i3]
                      *LSpa[i5][i3]*LSpb[i2][i5]*LSpb[i2][i5]*(
                2.0*LLs1(-LS[i1][i2],-LS[i3][i4],-LS[i5][i1],-LS[i3][i4])
                   +LL1(-LS[i1][i2],-LS[i3][i4])
                   +LL1(-LS[i5][i1],-LS[i3][i4]))
                  /LSpa[i5][i4]/LSpa[i4][i3]/LSpa[i2][i5]/LSpa[i2][i5]
                  /LS[i3][i4]/LS[i3][i4]
             +2.0/3.0*LSpa[i2][i3]*LSpa[i2][i3]
                     *LSpa[i4][i1]*LSpa[i4][i1]*LSpa[i4][i1]
                     *LSpb[i2][i4]*LSpb[i2][i4]*LSpb[i2][i4]
                     *LL2(-LS[i2][i3],-LS[i5][i1])
                     /LSpa[i4][i5]/LSpa[i5][i1]/LSpa[i2][i4]
                     /LS[i5][i1]/LS[i5][i1]/LS[i5][i1]
             -2.0/3.0*LSpa[i2][i1]*LSpa[i2][i1]
                     *LSpa[i5][i3]*LSpa[i5][i3]*LSpa[i5][i3]
                     *LSpb[i2][i5]*LSpb[i2][i5]*LSpb[i2][i5]
                     *LL2(-LS[i1][i2],-LS[i3][i4])
                     /LSpa[i5][i4]/LSpa[i4][i3]/LSpa[i2][i5]
                     /LS[i3][i4]/LS[i3][i4]/LS[i3][i4]
             +LL2(-LS[i3][i4],-LS[i5][i1])/LS[i5][i1]/LS[i5][i1]/LS[i5][i1]*(
	        1.0/3.0*LSpa[i1][i3]*LSpb[i2][i4]*LSpb[i2][i5]*(
                    LSpa[i1][i5]*LSpb[i5][i2]*LSpa[i2][i3]
                   -LSpa[i3][i4]*LSpb[i4][i2]*LSpa[i2][i1])/LSpa[i4][i5]
               +2.0/3.0*LSpa[i1][i2]*LSpa[i1][i2]*LSpa[i3][i4]*LSpa[i3][i4]
                    *LSpa[i4][i1]*LSpb[i2][i4]*LSpb[i2][i4]*LSpb[i2][i4]
                    /LSpa[i4][i5]/LSpa[i5][i1]/LSpa[i2][i4]
               -2.0/3.0*LSpa[i3][i2]*LSpa[i3][i2]*LSpa[i1][i5]*LSpa[i1][i5]
                    *LSpa[i5][i3]*LSpb[i2][i5]*LSpb[i2][i5]*LSpb[i2][i5]
                    /LSpa[i5][i4]/LSpa[i4][i3]/LSpa[i2][i5])
            +1.0/6.0*LSpa[i1][i3]*LSpa[i1][i3]*LSpa[i1][i3]*(
                    LSpa[i1][i5]*LSpb[i5][i2]*LSpa[i2][i3]
                   -LSpa[i3][i4]*LSpb[i4][i2]*LSpa[i2][i1])
              *LL0(-LS[i3][i4],-LS[i5][i1])/LS[i5][i1]
           /LSpa[i1][i2]/LSpa[i2][i3]/LSpa[i3][i4]/LSpa[i4][i5]/LSpa[i5][i1]
             +1.0/3.0*LSpb[i2][i4]*LSpb[i2][i4]*LSpb[i2][i5]*LSpb[i2][i5]
           /LSpb[i1][i2]/LSpb[i2][i3]/LSpb[i3][i4]/LSpa[i4][i5]/LSpb[i5][i1]
             -1.0/3.0*LSpa[i1][i2]*LSpa[i4][i1]*LSpa[i4][i1]
                     *LSpb[i2][i4]*LSpb[i2][i4]*LSpb[i2][i4]
                   /LSpa[i4][i5]/LSpa[i5][i1]/LSpa[i2][i4]
                      /LSpb[i2][i3]/LSpb[i3][i4]/LS[i5][i1]
             +1.0/3.0*LSpa[i3][i2]*LSpa[i5][i3]*LSpa[i5][i3]
                     *LSpb[i2][i5]*LSpb[i2][i5]*LSpb[i2][i5]
                   /LSpa[i5][i4]/LSpa[i4][i3]/LSpa[i2][i5]
                      /LSpb[i2][i1]/LSpb[i1][i5]/LS[i3][i4]
             +1.0/6.0*LSpa[i1][i3]*LSpa[i1][i3]*LSpb[i2][i4]*LSpb[i2][i5]
                      /LS[i3][i4]/LSpa[i4][i5]/LS[i5][i1];

return -( (Vf+Vs)*tree +Ff+Fs );
}

complex<double> Mpmpmm(int i1, int i2, int i3, int i4, int i5)
{

  //  Note: nothing should depend on mm.

  double mm=1.0;

  complex<double> l34,l51;

  if (LS[i3][i4]<0) {
     l34=log(mm/(-LS[i3][i4]));
  } else {
     l34=log(mm/LS[i3][i4])+I*PI;
  }

  if (LS[i5][i1]<0) {
     l51=log(mm/(-LS[i5][i1]));
  } else {
     l51=log(mm/LS[i5][i1])+I*PI;
  }

  complex<double> Vf=-0.5*(l34+l51);
  complex<double> Vs=-1.0/3.0*Vf;

  complex<double> tree=pow(LSpb[i1][i3],4)/LSpb[i1][i2]/LSpb[i2][i3]/LSpb[i3][i4]
                      /LSpb[i4][i5]/LSpb[i5][i1];

  complex<double> Ff=-LSpb[i1][i3]*LSpb[i1][i3]*LSpb[i4][i1]*LSpa[i2][i4]*LSpa[i2][i4]
                  *LLs1(-LS[i2][i3],-LS[i5][i1],-LS[i3][i4],-LS[i5][i1])
                   /LSpb[i4][i5]/LSpb[i5][i1]/LS[i5][i1]/LS[i5][i1]
             +LSpb[i1][i3]*LSpb[i1][i3]*LSpb[i5][i3]*LSpa[i2][i5]*LSpa[i2][i5]
                  *LLs1(-LS[i1][i2],-LS[i3][i4],-LS[i5][i1],-LS[i3][i4])
                   /LSpb[i3][i4]/LSpb[i4][i5]/LS[i3][i4]/LS[i3][i4]
             -0.5*LSpb[i1][i3]*LSpb[i1][i3]*LSpb[i1][i3]*(
                    LSpb[i1][i5]*LSpa[i5][i2]*LSpb[i2][i3]
                   -LSpb[i3][i4]*LSpa[i4][i2]*LSpb[i2][i1])
              *LL0(-LS[i3][i4],-LS[i5][i1])/LS[i5][i1]
           /LSpb[i1][i2]/LSpb[i2][i3]/LSpb[i3][i4]/LSpb[i4][i5]/LSpb[i5][i1];

  complex<double> Fs=-LSpb[i1][i2]*LSpb[i2][i3]*LSpb[i3][i4]*LSpb[i4][i1]
                      *LSpb[i4][i1]*LSpa[i2][i4]*LSpa[i2][i4]*(
                2.0*LLs1(-LS[i2][i3],-LS[i5][i1],-LS[i3][i4],-LS[i5][i1])
                   +LL1(-LS[i2][i3],-LS[i5][i1])
                   +LL1(-LS[i3][i4],-LS[i5][i1]))
                  /LSpb[i4][i5]/LSpb[i5][i1]/LSpb[i2][i4]/LSpb[i2][i4]
                  /LS[i5][i1]/LS[i5][i1]
             +LSpb[i3][i2]*LSpb[i2][i1]*LSpb[i1][i5]*LSpb[i5][i3]
                      *LSpb[i5][i3]*LSpa[i2][i5]*LSpa[i2][i5]*(
                2.0*LLs1(-LS[i1][i2],-LS[i3][i4],-LS[i5][i1],-LS[i3][i4])
                   +LL1(-LS[i1][i2],-LS[i3][i4])
                   +LL1(-LS[i5][i1],-LS[i3][i4]))
                  /LSpb[i5][i4]/LSpb[i4][i3]/LSpb[i2][i5]/LSpb[i2][i5]
                  /LS[i3][i4]/LS[i3][i4]
             +2.0/3.0*LSpb[i2][i3]*LSpb[i2][i3]
                     *LSpb[i4][i1]*LSpb[i4][i1]*LSpb[i4][i1]
                     *LSpa[i2][i4]*LSpa[i2][i4]*LSpa[i2][i4]
                     *LL2(-LS[i2][i3],-LS[i5][i1])
                     /LSpb[i4][i5]/LSpb[i5][i1]/LSpb[i2][i4]
                     /LS[i5][i1]/LS[i5][i1]/LS[i5][i1]
             -2.0/3.0*LSpb[i2][i1]*LSpb[i2][i1]
                     *LSpb[i5][i3]*LSpb[i5][i3]*LSpb[i5][i3]
                     *LSpa[i2][i5]*LSpa[i2][i5]*LSpa[i2][i5]
                     *LL2(-LS[i1][i2],-LS[i3][i4])
                     /LSpb[i5][i4]/LSpb[i4][i3]/LSpb[i2][i5]
                     /LS[i3][i4]/LS[i3][i4]/LS[i3][i4]
             +LL2(-LS[i3][i4],-LS[i5][i1])/LS[i5][i1]/LS[i5][i1]/LS[i5][i1]*(
	        1.0/3.0*LSpb[i1][i3]*LSpa[i2][i4]*LSpa[i2][i5]*(
                    LSpb[i1][i5]*LSpa[i5][i2]*LSpb[i2][i3]
                   -LSpb[i3][i4]*LSpa[i4][i2]*LSpb[i2][i1])/LSpb[i4][i5]
               +2.0/3.0*LSpb[i1][i2]*LSpb[i1][i2]*LSpb[i3][i4]*LSpb[i3][i4]
                    *LSpb[i4][i1]*LSpa[i2][i4]*LSpa[i2][i4]*LSpa[i2][i4]
                    /LSpb[i4][i5]/LSpb[i5][i1]/LSpb[i2][i4]
               -2.0/3.0*LSpb[i3][i2]*LSpb[i3][i2]*LSpb[i1][i5]*LSpb[i1][i5]
                    *LSpb[i5][i3]*LSpa[i2][i5]*LSpa[i2][i5]*LSpa[i2][i5]
                    /LSpb[i5][i4]/LSpb[i4][i3]/LSpb[i2][i5])
            +1.0/6.0*LSpb[i1][i3]*LSpb[i1][i3]*LSpb[i1][i3]*(
                    LSpb[i1][i5]*LSpa[i5][i2]*LSpb[i2][i3]
                   -LSpb[i3][i4]*LSpa[i4][i2]*LSpb[i2][i1])
              *LL0(-LS[i3][i4],-LS[i5][i1])/LS[i5][i1]
           /LSpb[i1][i2]/LSpb[i2][i3]/LSpb[i3][i4]/LSpb[i4][i5]/LSpb[i5][i1]
             +1.0/3.0*LSpa[i2][i4]*LSpa[i2][i4]*LSpa[i2][i5]*LSpa[i2][i5]
           /LSpa[i1][i2]/LSpa[i2][i3]/LSpa[i3][i4]/LSpb[i4][i5]/LSpa[i5][i1]
             -1.0/3.0*LSpb[i1][i2]*LSpb[i4][i1]*LSpb[i4][i1]
                     *LSpa[i2][i4]*LSpa[i2][i4]*LSpa[i2][i4]
                   /LSpb[i4][i5]/LSpb[i5][i1]/LSpb[i2][i4]
                      /LSpa[i2][i3]/LSpa[i3][i4]/LS[i5][i1]
             +1.0/3.0*LSpb[i3][i2]*LSpb[i5][i3]*LSpb[i5][i3]
                     *LSpa[i2][i5]*LSpa[i2][i5]*LSpa[i2][i5]
                   /LSpb[i5][i4]/LSpb[i4][i3]/LSpb[i2][i5]
                      /LSpa[i2][i1]/LSpa[i1][i5]/LS[i3][i4]
             +1.0/6.0*LSpb[i1][i3]*LSpb[i1][i3]*LSpa[i2][i4]*LSpa[i2][i5]
                      /LS[i3][i4]/LSpb[i4][i5]/LS[i5][i1];

return ( (Vf+Vs)*tree +Ff+Fs );
}




















complex<double> LL0(double r1, double r2)
{

  double r=r1/r2;

  complex<double> lr;

  if (r<0) {
      if (r1<0) {
          lr=log(-r)-I*PI;
      } else {
          lr=log(-r)+I*PI;
      }
  } else {
      lr=log(r);
  }

  return lr/(1.0-r);
}

complex<double> LL1(double r1, double r2)
{
  double r=r1/r2;

  complex<double> lr;

  if (r<0) {
      if (r1<0) {
          lr=log(-r)-I*PI;
      } else {
          lr=log(-r)+I*PI;
      }
  } else {
      lr=log(r);
  }

  double rm=1.0-r;

  return (lr+rm)/rm/rm;
}

complex<double> LL2(double r1, double r2)
{

  double r=r1/r2;
  double rm=1.0-r;

  complex<double> lr;

  if (r<0) {
      if (r1<0) {
          lr=log(-r)-I*PI;
      } else {
          lr=log(-r)+I*PI;
      }
  } else {
      lr=log(r);
  }

  return (lr-(r-1.0/r)/2.0)/rm/rm/rm;
}

complex<double> LLs1(double r1, double r2, double r3, double r4)
{
  double R1=r1/r2;
  double R2=r3/r4;

  double r1m=1.0-R1;
  double r2m=1.0-R2;
  double rm=1.0-R1-R2;

  complex<double> lr1, lr2;

  complex<double> cLi21, cLi22;

  if (R1<0) {
      if (r1<0) {
          lr1=log(-R1)-I*PI;
          cLi21=I*PI*log(r1m);
      } else {
          lr1=log(-R1)+I*PI;
          cLi21=-I*PI*log(r1m);
      }
  } else {
      lr1=log(R1);
      cLi21=0.0;
  }

  if (R2<0) {
      if (r3<0) {
          lr2=log(-R2)-I*PI;
          cLi22=I*PI*log(r2m);
      } else {
          lr2=log(-R2)+I*PI;
          cLi22=-I*PI*log(r2m);
      }
  } else {
      lr2=log(R2);
      cLi22=0.0;
  }

  return (LLi2(r1m)+cLi21+LLi2(r2m)+cLi22
      +lr1*lr2-PI*PI/6.0+rm*(LL0(r1,r2)+LL0(r3,r4)))
             /rm/rm;
}







const int NCUT = 22;

double LLi(const double &x)
{
   if (x<=0.5) {
     double L=0.0;
     double n=0.0;
     double r=1.0;
     for(int i=0;i<NCUT;i++) {
       n+=1.0;
       r*=x;
       if (r<=1.0e-10)
          return L;
       L+=r/n/n;
     }
     return L;
   } else {
     double xm = (1.0-x);          
     double L = PI*PI/6.0-log(x)*log(xm);
     double n = 0.0;
     double r = 1.0;
     for(int i=0;i<NCUT;i++) {
       n+=1.0;
       r*=xm;
       if (r<=1.0e-10) 
          return L;
       L-=r/n/n;
     }
     return L;
   }
}

double LLi2(double xx)
{
   double x=xx;
   double L=0.0;
   double b=1.0;
   if (x<0.0) {
     L+=PI*PI/6.0-log(1.0-x)*log(-x);
     b=-1.0;
     x=1.0-x;
   }
   if (x>1.0) {
     double lx=log(x);
     L+=b*(PI*PI/3.0-0.5*lx*lx);
     b=-b;
     x=1.0/x;
   }
   return L+b*LLi(x);
}

double LLi2(double rr, double ttheta)
{
   double r=rr;
   double theta=ttheta;
   if (fabs(r-1.0)<1.0e-10) {
       return pow(PI-fabs(theta),2)/4.0 - PI*PI/12.0;
   }
   double L;
   double lr;
   double b;
   if (r>1.0) {
      r=1.0/r;
      double lr=log(r);
      theta=-theta;
      b=-1.0;
      L= PI*PI/3.0-lr*lr/2.0+theta*theta/2.0-PI*fabs(theta);
   } else {
      b=1.0;
      L=0.0;
   }
   double x=r*cos(theta);
   if (x>0.5) {
      x=1.0-x;
      double y=-r*sin(theta);
      double rold=r;
      double told=theta;
      r=sqrt(x*x+y*y);
      theta=atan2(y,x);
      L+=b*(PI*PI/6.0-log(r)*log(rold)+theta*told);
      b=-b;
   }
   if (r<0.53) {
     double n=0.0;
     double xr=1.0;
     for(int i=0;i<NCUT;i++) {
       n+=1.0;
       xr*=r;
       double zr=xr*cos(n*theta);
       if (fabs(zr)<1.0e-10)
          return L;
       L+=b*zr/n/n;
     }
     return L;
   } else {
     double rm = (1.0-r);          
     double sh=2.0*sin(theta/2.0);
     L += b*(pow(PI-fabs(theta),2)/4.0-PI*PI/12.0-log(r)*log(fabs(sh)));
     double Zn = 0.0;
     double k = 0.0;
     double xr=rm;
     double Sn=1.0;
     for(int i=0;i<NCUT;i++) {
        k+=1.0;
        xr*=rm;
        Sn/=sh;
        Zn+=Sn*cos( (theta-PI)*k/2.0 )/k;
        double zr=Zn*xr;          
        if (fabs(zr)<1.0e-10) 
            return L;
        L-=b*zr/(k+1.0);
     }
     return L;
   }
}

// The following are from Lance Dixon:

const double PISQ6=PI*PI/6.0;
const double ZETA3=1.202056903;

// this version uses 't Hooft and Veltman's change of variable
// good to ~ 10^(-16)

/*
double li2(double x){
 double x_0 = -0.30;
 double x_1 = 0.25;
 double x_2 = 0.51;
 if (x == 1) return PISQ6;
 if (x <= x_0){ 
   double temp = log(abs(1.0-x));
   return -li2(-x/(1.0-x)) - temp*temp/2 ; }
 else if (x < x_1){
   double z = - log(1.0-x);
   double temp = z*(1.0-z/4.0*(1.0-z/9.0*(1.0-z*z/100.0
                  *(1.0-5.0*z*z/294.0*(1.0-7.0*z*z/360.0
                  *(1.0-5.0*z*z/242.0*(1.0-7601.0*z*z/354900.0
                  *(1.0-91.0*z*z/4146.0*(1.0-3617.0*z*z/161840.0)
                   ))))))));
   return temp; }
   else if (x < x_2) return - li2(-x) + li2(x*x)/2.0 ;
   else { return PISQ6 - li2(1.0-x) 
                  - log(abs(x))*log(abs(1.0-x)) ; }
}
*/

double LLi3(double x){
 double x_0 = -1.0;
 double x_1 = -0.85;
 double x_2 = 0.25;
 double x_3 = 0.63;
 double x_4 =  1.0;
 if (x == 1) return ZETA3;
 if (x <= x_0){ 
   double lnx = log(-x);
   return LLi3(1.0/x) - PISQ6*lnx - lnx*lnx*lnx/6.0; }
 else if (x < x_1){
   return LLi3(x*x)/4.0 - LLi3(-x); }
   else if (x < x_2){
     double z = - log(1.0-x);
     double temp = z*(1.0-3.0*z/8.0*(1.0-17.0*z/81.0*(1.0-15*z/136.0
                    *(1.0-28.0*z/1875.0*(1.0+5.0*z/8.0*(1.0-304.0*z/7203.0
                    *(1.0+945.0*z/2432.0*(1.0-44.0*z/675.0*(1.0+7.0*z/24.0
                    *(1.0-26104.0*z/307461.0*(1.0+1925.0*z/8023.0
                    *(1.0-53598548.0*z/524808375.0
                    *(1.0+22232925.0*z/107197096.0
                     )))))))))))));
     return temp; }
     else if (x < x_3){
       return LLi3(x*x)/4.0 - LLi3(-x); }
       else if (x < x_4){
         double ln1x = log(1.0-x); 
         return -LLi3(1.0-x) - LLi3(-x/(1.0-x)) + ZETA3 + PISQ6*ln1x
           - log(x)*ln1x*ln1x/2.0 + ln1x*ln1x*ln1x/6.0; }
       else { 
         double lnx = log(x);
         return LLi3(1/x) + 2.0*PISQ6*lnx - lnx*lnx*lnx/6.0; }
}

// the tetralog, Li4.
// good to ~ 1.7 * 10^(-12) (worst areas:  x = 0.9501, also x = - 0.9701)

double LLi4(double x){
 double x_0 = -1.0;
 double x_1 = -0.97;
 double x_2 = 0.25;
 double x_3 = 0.95;
 double x_4 =  1.0;
 if (x == -1) return -0.35 * PISQ6*PISQ6;
 if (x == 1) return 0.4 * PISQ6*PISQ6;
 if (x <= x_0){ 
   double lnx = log(-x);
   return - LLi4(1./x) - 0.5 * PISQ6*lnx*lnx
          - 1./24. * lnx*lnx*lnx*lnx - 0.7 * PISQ6*PISQ6; }
 else if (x < x_1){
   return LLi4(x*x)/8. - LLi4(-x); }
   else if (x < x_2){
     double z = - log(1.-x);
     double temp = z*(1.-7.*z/16.*(1.-151.*z/567.*(1.0-411.*z/2416.
                    *(1.-24986.*z/256875.*(1.-805.*z/49972.
                    *(1.+583406.*z/1159683.*(1.-7455.*z/137272.
                    *(1.+659921.*z/2444175.*(1.-251559.*z/2639684.
                    *(1.+24259894.*z/136410197.*(1.-30625595.*z/218339046.
                    *(1.+2134239258113.*z/16698772722450.
                    *(1.-1640443805715.*z/8536957032452.
                     )))))))))))));
     return temp; }
     else if (x < x_3){
       return LLi4(x*x)/8. - LLi4(-x); }
       else if (x < x_4){
        double y = 1.-x; 
        double lny = log(y);
        return 0.4*PISQ6*PISQ6 + ZETA3 * log(1.-y) + 0.5*PISQ6 * y*y 
         + (-11./36.+0.5*PISQ6+1./6.*lny) * y*y*y
         + (11./24.*PISQ6+1./4.*lny-19./48.) * y*y*y*y
         + (7./24.*lny+5./12.*PISQ6-599./1440.) * y*y*y*y*y
         + (137./360.*PISQ6+5./16.*lny-79./192.) * y*y*y*y*y*y
         + (7./20.*PISQ6-3343./8400.+29./90.*lny) * y*y*y*y*y*y*y
          + (363./1120*PISQ6-21977./57600.+469./1440.*lny) * y*y*y*y*y*y*y*y;
                        }
       else { 
         double lnx = log(x);
         return - LLi4(1/x) + PISQ6*lnx*lnx 
                - 1./24. * lnx*lnx*lnx*lnx + 0.8 * PISQ6*PISQ6; }
}
























void ampq2g1p2::ggg_background(int p1, int p2, int p3,int p4, int p5,
                        std::vector<std::complex<double> >& helamps) const
{
/*
//  (1+2+3+4+5+)

  helamps[0] = Mppppp(p1,p2,p3,p4,p5)+Mppppp(p1,p4,p5,p2,p3)+Mppppp(p1,p2,p4,p5,p3)
   +Mppppp(p1,p2,p3,p5,p4)+Mppppp(p1,p5,p4,p2,p3)+Mppppp(p1,p2,p5,p4,p3)
   +Mppppp(p1,p2,p4,p3,p5)+Mppppp(p1,p5,p2,p3,p4)+Mppppp(p1,p4,p2,p5,p3)
   +Mppppp(p1,p2,p5,p3,p4)+Mppppp(p1,p4,p2,p3,p5)+Mppppp(p1,p5,p2,p4,p3);

//  (1-2+3+4+5+)

  helamps[1] = Mmpppp(p1,p2,p3,p4,p5)+Mmpppp(p1,p4,p5,p2,p3)+Mmpppp(p1,p2,p4,p5,p3)
   +Mmpppp(p1,p2,p3,p5,p4)+Mmpppp(p1,p5,p4,p2,p3)+Mmpppp(p1,p2,p5,p4,p3)
   +Mmpppp(p1,p2,p4,p3,p5)+Mmpppp(p1,p5,p2,p3,p4)+Mmpppp(p1,p4,p2,p5,p3)
   +Mmpppp(p1,p2,p5,p3,p4)+Mmpppp(p1,p4,p2,p3,p5)+Mmpppp(p1,p5,p2,p4,p3);

//  (1+2-3+4+5+)

  helamps[2] = Mmpppp(p2,p3,p4,p5,p1)+Mmpppp(p2,p3,p1,p4,p5)+Mmpppp(p2,p4,p5,p3,p1)
   +Mmpppp(p2,p3,p5,p4,p1)+Mmpppp(p2,p3,p1,p5,p4)+Mmpppp(p2,p5,p4,p3,p1)
   +Mmpppp(p2,p4,p3,p5,p1)+Mmpppp(p2,p3,p4,p1,p5)+Mmpppp(p2,p5,p3,p1,p4)
   +Mmpppp(p2,p5,p3,p4,p1)+Mmpppp(p2,p3,p5,p1,p4)+Mmpppp(p2,p4,p3,p1,p5);

//  (1+2+3-4+5+)

  helamps[3] = Mmpppp(p3,p4,p5,p1,p2)+Mmpppp(p3,p1,p4,p5,p2)+Mmpppp(p3,p1,p2,p4,p5)
   +Mmpppp(p3,p5,p4,p1,p2)+Mmpppp(p3,p1,p5,p4,p2)+Mmpppp(p3,p1,p2,p5,p4)
   +Mmpppp(p3,p5,p1,p2,p4)+Mmpppp(p3,p4,p1,p5,p2)+Mmpppp(p3,p1,p4,p2,p5)
   +Mmpppp(p3,p4,p1,p2,p5)+Mmpppp(p3,p5,p1,p4,p2)+Mmpppp(p3,p1,p5,p2,p4);

//  (1+2+3+4-5+)

  helamps[4] = Mmpppp(p4,p5,p1,p2,p3)+Mmpppp(p4,p5,p2,p3,p1)+Mmpppp(p4,p5,p3,p1,p2)
   +Mmpppp(p4,p1,p2,p3,p5)+Mmpppp(p4,p2,p3,p1,p5)+Mmpppp(p4,p3,p1,p2,p5)
   +Mmpppp(p4,p3,p5,p1,p2)+Mmpppp(p4,p1,p5,p2,p3)+Mmpppp(p4,p2,p5,p3,p1)
   +Mmpppp(p4,p1,p2,p5,p3)+Mmpppp(p4,p2,p3,p5,p1)+Mmpppp(p4,p3,p1,p5,p2);

//  (1+2+3+4+5-)

  helamps[5] = Mmpppp(p5,p1,p2,p3,p4)+Mmpppp(p5,p2,p3,p1,p4)+Mmpppp(p5,p3,p1,p2,p4)
   +Mmpppp(p5,p4,p1,p2,p3)+Mmpppp(p5,p4,p2,p3,p1)+Mmpppp(p5,p4,p3,p1,p2)
   +Mmpppp(p5,p1,p2,p4,p3)+Mmpppp(p5,p2,p3,p4,p1)+Mmpppp(p5,p3,p1,p4,p2)
   +Mmpppp(p5,p3,p4,p1,p2)+Mmpppp(p5,p1,p4,p2,p3)+Mmpppp(p5,p2,p4,p3,p1);

//  (1-2-3+4+5+)

  helamps[6] = Mmmppp(p1,p2,p3,p4,p5)+Mmpmpp(p2,p3,p1,p4,p5)+Mmmppp(p1,p2,p4,p5,p3)
   +Mmmppp(p1,p2,p3,p5,p4)+Mmpmpp(p2,p3,p1,p5,p4)+Mmmppp(p1,p2,p5,p4,p3)
   +Mmmppp(p1,p2,p4,p3,p5)+Mmpmpp(p1,p5,p2,p3,p4)+Mmpmpp(p1,p4,p2,p5,p3)
   +Mmmppp(p1,p2,p5,p3,p4)+Mmpmpp(p1,p4,p2,p3,p5)+Mmpmpp(p1,p5,p2,p4,p3);

//  (1-2+3-4+5+)

  helamps[7] = Mmpmpp(p1,p2,p3,p4,p5)+Mmmppp(p3,p1,p4,p5,p2)+Mmmppp(p3,p1,p2,p4,p5)
   +Mmpmpp(p1,p2,p3,p5,p4)+Mmmppp(p3,p1,p5,p4,p2)+Mmmppp(p3,p1,p2,p5,p4)
   +Mmpmpp(p3,p5,p1,p2,p4)+Mmpmpp(p3,p4,p1,p5,p2)+Mmmppp(p3,p1,p4,p2,p5)
   +Mmpmpp(p3,p4,p1,p2,p5)+Mmpmpp(p3,p5,p1,p4,p2)+Mmmppp(p3,p1,p5,p2,p4);

//  (1-2+3+4-5+)

  helamps[8] = Mmpmpp(p4,p5,p1,p2,p3)+Mmmppp(p1,p4,p5,p2,p3)+Mmpmpp(p1,p2,p4,p5,p3)
   +Mmmppp(p4,p1,p2,p3,p5)+Mmpmpp(p1,p5,p4,p2,p3)+Mmpmpp(p4,p3,p1,p2,p5)
   +Mmpmpp(p1,p2,p4,p3,p5)+Mmmppp(p4,p1,p5,p2,p3)+Mmmppp(p1,p4,p2,p5,p3)
   +Mmmppp(p4,p1,p2,p5,p3)+Mmmppp(p1,p4,p2,p3,p5)+Mmpmpp(p4,p3,p1,p5,p2);

//  (1-2+3+4+5-)

  helamps[9] = Mmmppp(p5,p1,p2,p3,p4)+Mmpmpp(p1,p4,p5,p2,p3)+Mmpmpp(p5,p3,p1,p2,p4)
   +Mmpmpp(p5,p4,p1,p2,p3)+Mmmppp(p1,p5,p4,p2,p3)+Mmpmpp(p1,p2,p5,p4,p3)
   +Mmmppp(p5,p1,p2,p4,p3)+Mmmppp(p1,p5,p2,p3,p4)+Mmpmpp(p5,p3,p1,p4,p2)
   +Mmpmpp(p1,p2,p5,p3,p4)+Mmmppp(p5,p1,p4,p2,p3)+Mmmppp(p1,p5,p2,p4,p3);

//  (1+2-3-4+5+)

  helamps[10] = Mmmppp(p2,p3,p4,p5,p1)+Mmmppp(p2,p3,p1,p4,p5)+Mmpmpp(p3,p1,p2,p4,p5)
   +Mmmppp(p2,p3,p5,p4,p1)+Mmmppp(p2,p3,p1,p5,p4)+Mmpmpp(p3,p1,p2,p5,p4)
   +Mmpmpp(p2,p4,p3,p5,p1)+Mmmppp(p2,p3,p4,p1,p5)+Mmpmpp(p2,p5,p3,p1,p4)
   +Mmpmpp(p2,p5,p3,p4,p1)+Mmmppp(p2,p3,p5,p1,p4)+Mmpmpp(p2,p4,p3,p1,p5);

//  (1+2-3+4-5+)

  helamps[11] = Mmpmpp(p2,p3,p4,p5,p1)+Mmpmpp(p4,p5,p2,p3,p1)+Mmmppp(p2,p4,p5,p3,p1)
   +Mmpmpp(p4,p1,p2,p3,p5)+Mmmppp(p4,p2,p3,p1,p5)+Mmpmpp(p2,p5,p4,p3,p1)
   +Mmmppp(p2,p4,p3,p5,p1)+Mmpmpp(p2,p3,p4,p1,p5)+Mmmppp(p4,p2,p5,p3,p1)
   +Mmpmpp(p4,p1,p2,p5,p3)+Mmmppp(p4,p2,p3,p5,p1)+Mmmppp(p2,p4,p3,p1,p5);

//  (1+2-3+4+5-)

  helamps[12] = Mmpmpp(p5,p1,p2,p3,p4)+Mmmppp(p5,p2,p3,p1,p4)+Mmpmpp(p2,p4,p5,p3,p1)
   +Mmpmpp(p2,p3,p5,p4,p1)+Mmpmpp(p5,p4,p2,p3,p1)+Mmmppp(p2,p5,p4,p3,p1)
   +Mmpmpp(p5,p1,p2,p4,p3)+Mmmppp(p5,p2,p3,p4,p1)+Mmmppp(p2,p5,p3,p1,p4)
   +Mmmppp(p2,p5,p3,p4,p1)+Mmpmpp(p2,p3,p5,p1,p4)+Mmmppp(p5,p2,p4,p3,p1);

//  (1+2+3-4-5+)

  helamps[13] = Mmmppp(p3,p4,p5,p1,p2)+Mmpmpp(p3,p1,p4,p5,p2)+Mmpmpp(p4,p5,p3,p1,p2)
   +Mmpmpp(p3,p5,p4,p1,p2)+Mmpmpp(p4,p2,p3,p1,p5)+Mmmppp(p4,p3,p1,p2,p5)
   +Mmmppp(p4,p3,p5,p1,p2)+Mmmppp(p3,p4,p1,p5,p2)+Mmpmpp(p3,p1,p4,p2,p5)
   +Mmmppp(p3,p4,p1,p2,p5)+Mmpmpp(p4,p2,p3,p5,p1)+Mmmppp(p4,p3,p1,p5,p2);

//  (1+2+3-4+5-)

  helamps[14] = Mmpmpp(p3,p4,p5,p1,p2)+Mmpmpp(p5,p2,p3,p1,p4)+Mmmppp(p5,p3,p1,p2,p4)
   +Mmmppp(p3,p5,p4,p1,p2)+Mmpmpp(p3,p1,p5,p4,p2)+Mmpmpp(p5,p4,p3,p1,p2)
   +Mmmppp(p3,p5,p1,p2,p4)+Mmpmpp(p5,p2,p3,p4,p1)+Mmmppp(p5,p3,p1,p4,p2)
   +Mmmppp(p5,p3,p4,p1,p2)+Mmmppp(p3,p5,p1,p4,p2)+Mmpmpp(p3,p1,p5,p2,p4);

//  (1+2+3+4-5-)

  helamps[15] = Mmmppp(p4,p5,p1,p2,p3)+Mmmppp(p4,p5,p2,p3,p1)+Mmmppp(p4,p5,p3,p1,p2)
   +Mmmppp(p5,p4,p1,p2,p3)+Mmmppp(p5,p4,p2,p3,p1)+Mmmppp(p5,p4,p3,p1,p2)
   +Mmpmpp(p4,p3,p5,p1,p2)+Mmpmpp(p4,p1,p5,p2,p3)+Mmpmpp(p4,p2,p5,p3,p1)
   +Mmpmpp(p5,p3,p4,p1,p2)+Mmpmpp(p5,p1,p4,p2,p3)+Mmpmpp(p5,p2,p4,p3,p1);

//  (1-2-3-4-5-)

  helamps[16] = Mmmmmm(p1,p2,p3,p4,p5)+Mmmmmm(p1,p4,p5,p2,p3)+Mmmmmm(p1,p2,p4,p5,p3)
   +Mmmmmm(p1,p2,p3,p5,p4)+Mmmmmm(p1,p5,p4,p2,p3)+Mmmmmm(p1,p2,p5,p4,p3)
   +Mmmmmm(p1,p2,p4,p3,p5)+Mmmmmm(p1,p5,p2,p3,p4)+Mmmmmm(p1,p4,p2,p5,p3)
   +Mmmmmm(p1,p2,p5,p3,p4)+Mmmmmm(p1,p4,p2,p3,p5)+Mmmmmm(p1,p5,p2,p4,p3);

//  (1+2-3-4-5-)

  helamps[17] = Mpmmmm(p1,p2,p3,p4,p5)+Mpmmmm(p1,p4,p5,p2,p3)+Mpmmmm(p1,p2,p4,p5,p3)
   +Mpmmmm(p1,p2,p3,p5,p4)+Mpmmmm(p1,p5,p4,p2,p3)+Mpmmmm(p1,p2,p5,p4,p3)
   +Mpmmmm(p1,p2,p4,p3,p5)+Mpmmmm(p1,p5,p2,p3,p4)+Mpmmmm(p1,p4,p2,p5,p3)
   +Mpmmmm(p1,p2,p5,p3,p4)+Mpmmmm(p1,p4,p2,p3,p5)+Mpmmmm(p1,p5,p2,p4,p3);

//  (1-2+3-4-5-)

  helamps[18] = Mpmmmm(p2,p3,p4,p5,p1)+Mpmmmm(p2,p3,p1,p4,p5)+Mpmmmm(p2,p4,p5,p3,p1)
   +Mpmmmm(p2,p3,p5,p4,p1)+Mpmmmm(p2,p3,p1,p5,p4)+Mpmmmm(p2,p5,p4,p3,p1)
   +Mpmmmm(p2,p4,p3,p5,p1)+Mpmmmm(p2,p3,p4,p1,p5)+Mpmmmm(p2,p5,p3,p1,p4)
   +Mpmmmm(p2,p5,p3,p4,p1)+Mpmmmm(p2,p3,p5,p1,p4)+Mpmmmm(p2,p4,p3,p1,p5);

//  (1-2-3+4-5-)

  helamps[19] = Mpmmmm(p3,p4,p5,p1,p2)+Mpmmmm(p3,p1,p4,p5,p2)+Mpmmmm(p3,p1,p2,p4,p5)
   +Mpmmmm(p3,p5,p4,p1,p2)+Mpmmmm(p3,p1,p5,p4,p2)+Mpmmmm(p3,p1,p2,p5,p4)
   +Mpmmmm(p3,p5,p1,p2,p4)+Mpmmmm(p3,p4,p1,p5,p2)+Mpmmmm(p3,p1,p4,p2,p5)
   +Mpmmmm(p3,p4,p1,p2,p5)+Mpmmmm(p3,p5,p1,p4,p2)+Mpmmmm(p3,p1,p5,p2,p4);

//  (1-2-3-4+5-)

  helamps[20] = Mpmmmm(p4,p5,p1,p2,p3)+Mpmmmm(p4,p5,p2,p3,p1)+Mpmmmm(p4,p5,p3,p1,p2)
   +Mpmmmm(p4,p1,p2,p3,p5)+Mpmmmm(p4,p2,p3,p1,p5)+Mpmmmm(p4,p3,p1,p2,p5)
   +Mpmmmm(p4,p3,p5,p1,p2)+Mpmmmm(p4,p1,p5,p2,p3)+Mpmmmm(p4,p2,p5,p3,p1)
   +Mpmmmm(p4,p1,p2,p5,p3)+Mpmmmm(p4,p2,p3,p5,p1)+Mpmmmm(p4,p3,p1,p5,p2);

//  (1-2-3-4-5+)

  helamps[21] = Mpmmmm(p5,p1,p2,p3,p4)+Mpmmmm(p5,p2,p3,p1,p4)+Mpmmmm(p5,p3,p1,p2,p4)
   +Mpmmmm(p5,p4,p1,p2,p3)+Mpmmmm(p5,p4,p2,p3,p1)+Mpmmmm(p5,p4,p3,p1,p2)
   +Mpmmmm(p5,p1,p2,p4,p3)+Mpmmmm(p5,p2,p3,p4,p1)+Mpmmmm(p5,p3,p1,p4,p2)
   +Mpmmmm(p5,p3,p4,p1,p2)+Mpmmmm(p5,p1,p4,p2,p3)+Mpmmmm(p5,p2,p4,p3,p1);

//  (1+2+3-4-5-)

  helamps[22] = Mppmmm(p1,p2,p3,p4,p5)+Mpmpmm(p2,p3,p1,p4,p5)+Mppmmm(p1,p2,p4,p5,p3)
   +Mppmmm(p1,p2,p3,p5,p4)+Mpmpmm(p2,p3,p1,p5,p4)+Mppmmm(p1,p2,p5,p4,p3)
   +Mppmmm(p1,p2,p4,p3,p5)+Mpmpmm(p1,p5,p2,p3,p4)+Mpmpmm(p1,p4,p2,p5,p3)
   +Mppmmm(p1,p2,p5,p3,p4)+Mpmpmm(p1,p4,p2,p3,p5)+Mpmpmm(p1,p5,p2,p4,p3);

//  (1+2-3+4-5-)

  helamps[23] = Mpmpmm(p1,p2,p3,p4,p5)+Mppmmm(p3,p1,p4,p5,p2)+Mppmmm(p3,p1,p2,p4,p5)
   +Mpmpmm(p1,p2,p3,p5,p4)+Mppmmm(p3,p1,p5,p4,p2)+Mppmmm(p3,p1,p2,p5,p4)
   +Mpmpmm(p3,p5,p1,p2,p4)+Mpmpmm(p3,p4,p1,p5,p2)+Mppmmm(p3,p1,p4,p2,p5)
   +Mpmpmm(p3,p4,p1,p2,p5)+Mpmpmm(p3,p5,p1,p4,p2)+Mppmmm(p3,p1,p5,p2,p4);

//  (1+2-3-4+5-)

  helamps[24] = Mpmpmm(p4,p5,p1,p2,p3)+Mppmmm(p1,p4,p5,p2,p3)+Mpmpmm(p1,p2,p4,p5,p3)
   +Mppmmm(p4,p1,p2,p3,p5)+Mpmpmm(p1,p5,p4,p2,p3)+Mpmpmm(p4,p3,p1,p2,p5)
   +Mpmpmm(p1,p2,p4,p3,p5)+Mppmmm(p4,p1,p5,p2,p3)+Mppmmm(p1,p4,p2,p5,p3)
   +Mppmmm(p4,p1,p2,p5,p3)+Mppmmm(p1,p4,p2,p3,p5)+Mpmpmm(p4,p3,p1,p5,p2);

//  (1+2-3-4-5+)

  helamps[25] = Mppmmm(p5,p1,p2,p3,p4)+Mpmpmm(p1,p4,p5,p2,p3)+Mpmpmm(p5,p3,p1,p2,p4)
   +Mpmpmm(p5,p4,p1,p2,p3)+Mppmmm(p1,p5,p4,p2,p3)+Mpmpmm(p1,p2,p5,p4,p3)
   +Mppmmm(p5,p1,p2,p4,p3)+Mppmmm(p1,p5,p2,p3,p4)+Mpmpmm(p5,p3,p1,p4,p2)
   +Mpmpmm(p1,p2,p5,p3,p4)+Mppmmm(p5,p1,p4,p2,p3)+Mppmmm(p1,p5,p2,p4,p3);

//  (1-2+3+4-5-)

  helamps[26] = Mppmmm(p2,p3,p4,p5,p1)+Mppmmm(p2,p3,p1,p4,p5)+Mpmpmm(p3,p1,p2,p4,p5)
   +Mppmmm(p2,p3,p5,p4,p1)+Mppmmm(p2,p3,p1,p5,p4)+Mpmpmm(p3,p1,p2,p5,p4)
   +Mpmpmm(p2,p4,p3,p5,p1)+Mppmmm(p2,p3,p4,p1,p5)+Mpmpmm(p2,p5,p3,p1,p4)
   +Mpmpmm(p2,p5,p3,p4,p1)+Mppmmm(p2,p3,p5,p1,p4)+Mpmpmm(p2,p4,p3,p1,p5);

//  (1-2+3-4+5-)

  helamps[27] = Mpmpmm(p2,p3,p4,p5,p1)+Mpmpmm(p4,p5,p2,p3,p1)+Mppmmm(p2,p4,p5,p3,p1)
   +Mpmpmm(p4,p1,p2,p3,p5)+Mppmmm(p4,p2,p3,p1,p5)+Mpmpmm(p2,p5,p4,p3,p1)
   +Mppmmm(p2,p4,p3,p5,p1)+Mpmpmm(p2,p3,p4,p1,p5)+Mppmmm(p4,p2,p5,p3,p1)
   +Mpmpmm(p4,p1,p2,p5,p3)+Mppmmm(p4,p2,p3,p5,p1)+Mppmmm(p2,p4,p3,p1,p5);

//  (1-2+3-4-5+)

  helamps[28] = Mpmpmm(p5,p1,p2,p3,p4)+Mppmmm(p5,p2,p3,p1,p4)+Mpmpmm(p2,p4,p5,p3,p1)
   +Mpmpmm(p2,p3,p5,p4,p1)+Mpmpmm(p5,p4,p2,p3,p1)+Mppmmm(p2,p5,p4,p3,p1)
   +Mpmpmm(p5,p1,p2,p4,p3)+Mppmmm(p5,p2,p3,p4,p1)+Mppmmm(p2,p5,p3,p1,p4)
   +Mppmmm(p2,p5,p3,p4,p1)+Mpmpmm(p2,p3,p5,p1,p4)+Mppmmm(p5,p2,p4,p3,p1);

//  (1-2-3+4+5-)

  helamps[29] = Mppmmm(p3,p4,p5,p1,p2)+Mpmpmm(p3,p1,p4,p5,p2)+Mpmpmm(p4,p5,p3,p1,p2)
   +Mpmpmm(p3,p5,p4,p1,p2)+Mpmpmm(p4,p2,p3,p1,p5)+Mppmmm(p4,p3,p1,p2,p5)
   +Mppmmm(p4,p3,p5,p1,p2)+Mppmmm(p3,p4,p1,p5,p2)+Mpmpmm(p3,p1,p4,p2,p5)
   +Mppmmm(p3,p4,p1,p2,p5)+Mpmpmm(p4,p2,p3,p5,p1)+Mppmmm(p4,p3,p1,p5,p2);

//  (1-2-3+4-5+)

  helamps[30] = Mpmpmm(p3,p4,p5,p1,p2)+Mpmpmm(p5,p2,p3,p1,p4)+Mppmmm(p5,p3,p1,p2,p4)
   +Mppmmm(p3,p5,p4,p1,p2)+Mpmpmm(p3,p1,p5,p4,p2)+Mpmpmm(p5,p4,p3,p1,p2)
   +Mppmmm(p3,p5,p1,p2,p4)+Mpmpmm(p5,p2,p3,p4,p1)+Mppmmm(p5,p3,p1,p4,p2)
   +Mppmmm(p5,p3,p4,p1,p2)+Mppmmm(p3,p5,p1,p4,p2)+Mpmpmm(p3,p1,p5,p2,p4);

//  (1-2-3-4+5+)

  helamps[31] = Mppmmm(p4,p5,p1,p2,p3)+Mppmmm(p4,p5,p2,p3,p1)+Mppmmm(p4,p5,p3,p1,p2)
   +Mppmmm(p5,p4,p1,p2,p3)+Mppmmm(p5,p4,p2,p3,p1)+Mppmmm(p5,p4,p3,p1,p2)
   +Mpmpmm(p4,p3,p5,p1,p2)+Mpmpmm(p4,p1,p5,p2,p3)+Mpmpmm(p4,p2,p5,p3,p1)
   +Mpmpmm(p5,p3,p4,p1,p2)+Mpmpmm(p5,p1,p4,p2,p3)+Mpmpmm(p5,p2,p4,p3,p1);
   * 
   * 
*/
}
