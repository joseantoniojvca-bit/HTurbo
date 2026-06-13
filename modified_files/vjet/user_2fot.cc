#include "Set_Root_user_2fot.h" 
#include "params.h"
#include "Inter_params.h"

// DYTurbo includes
#include "Kinematics.h"
#include "KinematicDefinitions.h"
#include "KinematicCuts.h"
#include "settings.h"
#include "coupling.h"
#include "phasespace.h"
#include "codata.h"

#include "scales.h"


void _Evento::Set_VEGAS(_Evento & Evnt){

    Evnt.iifborn  = _iifborn; // BORN
    Evnt.iifnnlo  = _iifnnlo; //  NLO
    Evnt.iinnlo   = _iinnlo; // NNLO
    Evnt.iihyylo   = _iihyylo; // NNLO
    Evnt.iirhyy   = _iirhyy; // NNLO
    
    
    Evnt.i0   = _i0;
    Evnt.ifnn = _ifnn;
    Evnt.inn  = _inn;
    
    Evnt.n0   = _n0;
    Evnt.fnn0 = _fnn0;
    Evnt.nn0  = _nn0;
    
    Evnt.ncl0  = _ncl0;
    Evnt.fnncl = _fnncl;
    Evnt.nncl  = _nncl;
};



void _Evento::Set_ps(_Evento & Evnt){

    // Set the Ecm squared
    Evnt._M_s = opts.sroot*opts.sroot;
    
    // Limit in the invariant Mass of the diphoton pair
    Evnt.q2min = phasespace::mmin*phasespace::mmin;
    Evnt.q2max = phasespace::mmax*phasespace::mmax;
};

void Integral_class::Cuts(event_hhc2ph p, bool &_cuts,double rng1,double rng2){
	   
  //no photon belongs the beam
  double Mff = sqrt(2.0*(p[-3]*p[-2]));
  if( Mff <  phasespace::mmin || Mff >  phasespace::mmax ) return;


  //Cuts for the pt of each photon and rapidity
  double fTmin1 = opts.phpt_cut_h;
  double fTmin2 = opts.phpt_cut_s;
  if (p[-3].perp() <  p[-2].perp() ) swap(fTmin1,fTmin2);
  if (p[-3].perp() < fTmin1 || p[-2].perp() < fTmin2) return;
  if (abs(p[-3].rapidity()) > opts.phycut || abs(p[-2].rapidity()) > opts.phycut) return;
  // -------------------------------------------------
  // crack (eta)
  // -------------------------------------------------
  crack_eta_max = opts.crack_eta_max;
  crack_eta_min = opts.crack_eta_min;
  // -------------------------------------------------
  if(opts.apply_crack)
    {
      if (abs(p[-3].rapidity()) < crack_eta_max && abs(p[-3].rapidity()) > crack_eta_min) return;
      if (abs(p[-2].rapidity()) < crack_eta_max && abs(p[-2].rapidity()) > crack_eta_min) return;
    }
  //----- cuts ----- 
  double pT12min = max(opts.qtcut,phasespace::qtmin); // Internal cut, don't touch this
  double pT12max = phasespace::qtmax;
  //  cout << " pT12min : " << pT12min << endl;
  // double pT12max = 40.0;  
  //----- photon pT and eta -----
  double pT12 = (p[-3]+p[-2]).perp();
  //pT12=1.0;
  if(rng1 == 19.0 && rng2 == 19.0 )
    { 
      // EXCLUSIVE CUTS KINEMATICS NON BORN REAL Interference LO    
      if(U_photon_izolation(p, opts.ph_iso_R0, opts.ph_iso_E0, opts.n_isol) == false) return;    
      
      // The integral on the reals of the interference at LO are finite
      if(opts.Dixon_eta_jet_cut) 
	{
	  if(fabs(p[1].rapidity()) < opts.eta_jet && p[1].perp() > pT12max) return;
          if(pT12 < pT12min) return;
	}
      else
	{
	if(pT12 < pT12min || pT12 > pT12max) return;// don't touch this internal cut!
	}
    }
  
  if(rng1 != 20.0 && rng1!=2.0 && rng1 != 19.0 && rng2 != 19.0 ){ // EXCLUSIVE CUTS KINEMATICS NON BORN
      
  //Separation between jets // not for diphoton studies!!!
  //Separacion entre los dos Jets 
  //double  Rjet12 = U_radius(p[1],p[2]);
  //if (Rjet12<1.0) return;
  //if(p[2].perp() < 20.0) return;
  //if(p[1].perp() < 20.0) return;
 
  if(rng2 != 20.0)
    {
      // EXCLUSIVE CUTS KINEMATICS NON BORN
      //  photon izolation: U_photon_izolation(p, R, eps or Etmax, n)
      if(U_photon_izolation(p, opts.ph_iso_R0, opts.ph_iso_E0, opts.n_isol ) == false) return;

      // The integral on the reals of the interference at LO are finite
      if(opts.Dixon_eta_jet_cut) 
	{
	  if(fabs(p[1].rapidity()) < opts.eta_jet && p[1].perp() > pT12max) return;
          if(pT12 < pT12min) return;
	}
      else
	{
	if(pT12 < pT12min || pT12 > pT12max) return;// don't touch this internal cut!
	}
     
    }
  
  }
  
  if (rng1==2.0 && rng2 ==20.0) 
    {

      // cout << "p[1].perp() = " << p[1].perp() << endl;
      // cout << "(p[-2]+p[-3]).perp() = " << (p[-2]+p[-3]).perp() << endl;
      if(opts.fixedorder) 
	{
	  if(p[1].perp()<pT12min || p[1].perp() > pT12max) return;  // EXCLUSIVE CUTS KINEMATICS NON BORN -- COUNTERTERM
	}
      else
	{ 
	  if(pT12 < pT12min || pT12 > pT12max) return;  // EXCLUSIVE CUTS KINEMATICS NON BORN -- COUNTERTERM  
	}
      // if(p[1].perp()<pT12min || p[1].perp() > pT12max) return;  // EXCLUSIVE CUTS KINEMATICS NON BORN -- COUNTERTERM
    }
  ///////////////// RADIO ENTRE LOS DOS FOTONES /////////////////////////
  // Separation between the two photons
//   double Rgg=U_radius(p[-2],p[-3]);
//   if (Rgg < opts.Rgamgam) return;
  ///////////////////////////////////////////////////////////////////////    
 _cuts = true;
	
};

void Integral_class::userfunc(event_hhc2ph p,double &weightfinal, amplitude_hhc2ph& amp,double rng1 , double rng2)
{
  // ----------------------------------------------------------------------------------------
  double ncont, wtopd, wroot;
  // ----------------------------------------------------------------------------------------
  if(contrib ==1.0) ncont=abs(_n0);
  else if(contrib ==2.0) ncont=abs(_fnn0);
  else ncont=abs(_nn0);
  // ----------------------------------------------------------------------------------------  	
  //----- type of the contribution -----
  static const double coef = gevfb * opts.aemmz * opts.aemmz ;
  // diphoton invariant mass
  double Mff = sqrt(2.0*(p[-3]*p[-2]));  
  phasespace::m = Mff;
  //----- photon pT and eta -----
  double pT12 = (p[-3]+p[-2]).perp();
  // Set the scales to evaluate the weight
  double mr2,mf2,mu2,scale=Mff;
  
  // --- Ojo con las escalas porque no remiten al inputfile de DYTURBO

  if(Scale=="Mff"){scale=Mff;};
  if(Scale=="MH"){scale= coupling::hmass;};
  if(Scale=="MT"){scale=sqrt(Mff*Mff+pT12*pT12);};
  
      //Set scales
      scales::set(Mff);
      scales::mcfm();
      double muf = scales::fac;
      double mur = scales::ren;
      mf2 = muf*muf;
      mr2 = mur*mur;

      //      cout << " muf : " << muf << " ; " << " mur : " << mur << endl;
 
  // // mu2 = Mff*Mff; 
  // mu2 =scale*scale;
  
  // mf2 = muF*muF*mu2;
  // mr2 = muR*muR*mu2; 
  
  // cout << " sqrt(mf2) : " << sqrt(mf2) << endl;
  // cout << " sqrt(mr2) : " << sqrt(mr2) << endl;
  
  // cout << " muF : " << muF << endl;
  // cout << " muR : " << muR << endl;
  
  // apply the cuts and compute the weight
    weightfinal=0.0;
    bool flcuts=false;
    Cuts(p,flcuts,rng1,rng2);
    if (flcuts) {  
      weight_hhc2ph wt = amp(mr2, mf2, coef,rng1,rng2);
      //cout << " desde user_func " << endl;
      //cout <<  wt  << endl;
      
      for(int i=0;i<17;i++) weightfinal +=wt[i]; 
    }    
    if (std::isinf( weightfinal) || std::isnan( weightfinal)) weightfinal=0.0;
    
    // weights for the histograms
    // wtopd = weightfinal*vweght;
    // wroot = weightfinal*vweght/ncont;
    
    // if contpt == 2U you don't have to save p[1] or p[2] in the event
    // because is a fake jet. p[1] and p[2] in such events (contpt == 2U)
    // are [0,0,0,0]


    // if the weight is non-zero, compute the diphoton variables and fill the histograms
    double var[50];
    // if (weightfinal!=0) {
    // 	// compute the diphoton variables and store in var the ones that will be used to fill the histograms
    // 	FillVars(var, p);
    // 	// fill the topd histograms (booked in Set_hist.f)
    // 	if (SaveTopdHistos) {
    // 	    for (int i=0; i<nhists; i++)
    // 		fillhist(i+1, var[i], wtopd);
    // 	}
	
	// // fill the ROOT histograms/ntuple (booked in main_c++.cc)
	// if (SaveRootNtuple) T->Fill(Mff,wroot);
	// if (SaveRootHistos) {
	//     for (int i=0; i<nhists; i++)
	// 	hists[i]->Fill(var[i],wroot);
	// }
    //}
};

double U_radius(const lorentzvector<double>& p1, const lorentzvector<double>& p2)
{
  static const double pi = 3.14159265358979323846;
  static const double twopi = 6.28318530717958647692;

  double deta = p1.rapidity() - p2.rapidity();
  double dphi = p1.phi() - p2.phi();

  if(dphi >= pi) dphi = fmod(pi+dphi, twopi) - pi;
  else if(dphi < -pi) dphi = -fmod(pi-dphi, twopi) + pi;

  return sqrt(deta*deta+dphi*dphi);
}

//   izolation cone
  struct izolation_cone {
    double pT, R;
  };

  //   lists of the izolation cones
  bounded_vector<izolation_cone> cone1, cone2;

  //   order the cones
  struct radius_sort {
    bool operator()(const izolation_cone& c1, const izolation_cone& c2) const {
      return c1.R < c2.R;
    }
  };

//#define Xi(E, x) (eps*std::pow((1.0-std::cos(x))/(1.0-std::cos(delta)),  n))
#define Xi(E, x) (eps)

bool U_photon_izolation(const event_hhc2ph& p, double delta, double eps, double n)
{

     // Etmax=eps now!!!
     // eps: if you use a fraction of the transverse momentum of the photon instead a total amount of energy
     // you have to write p[i].perp() multiplying.

  unsigned int np = p.upper();
  
  if(p[np].perp() < 1e-9) np--;
  cone1.resize(1,np); cone2.resize(1,np); 
  
  for(unsigned int i = 1; i <= np; i++) {
    cone1[i].pT = cone2[i].pT = p[i].perp();
    cone1[i].R = U_radius(p[-3], p[i]);
    cone2[i].R = U_radius(p[-2], p[i]);
  }
  
  if(np > 1) {
    sort(cone1.begin(), cone1.end(), radius_sort());
    sort(cone2.begin(), cone2.end(), radius_sort());
  }
     
        
  double Esum1 = 0.0, Esum2 = 0.0;
     // Etmax=eps now!!!
     // eps: if you use a fraction of the transverse momentum of the photon instead a total amount of energy
     // you have to write p[i].perp() multiplying.
     // And uncomment the next line.
      double pT1 = p[-3].perp(), pT2 = p[-2].perp();
  for(unsigned int i = 1; i <= np; i++) {
	
    Esum1 += delta > cone1[i].R ? cone1[i].pT : 0.0;
    Esum2 += delta > cone2[i].R ? cone2[i].pT : 0.0;
    if(Xi(pT1, (delta < cone1[i].R ? delta : cone1[i].R)) < Esum1) return false;
    if(Xi(pT2, (delta < cone2[i].R ? delta : cone2[i].R)) < Esum2) return false;
  }

  return true;
}


void FillVars(double* var, const event_hhc2ph& p) {
    // diphoton invariant mass
    double Mff = sqrt(2.0*(p[-3]*p[-2]));

    // photon pT
    double pT1 = p[-3].perp(), pT2 = p[-2].perp();

    // Jet pT
    double pT = p[1].perp(); 
    
    // diphoton pt
    double pT12 = (p[-3]+p[-2]).perp();
    
    // photon rapidities
    double eta1p = p[-3].rapidity(), eta2p = p[-2].rapidity();
    
    // diphoton rapidity
    double etapar=(p[-2]+p[-3]).rapidity();
    
    // rapidity difference
    double deltay12= eta1p-eta2p;
    
    //Costheta*
    double CosTHETAStarApp =  tanh(deltay12*0.5) ;
    double CosTHETAStar =  abs(2.0*pT1*pT2*sinh(eta1p-eta2p)/(Mff*sqrt(Mff*Mff+pT12*pT12)));
    
    // azimuth difference
    static const double pi = 3.14159265358979323846;
    double PHIgg=fabs(p[-3].phi()-p[-2].phi());
    if (PHIgg>pi) {PHIgg=2.0*pi-PHIgg;}
    
    if(pT1 < pT2) {	
	std::swap(pT1, pT2);
	std::swap(eta1p, eta2p);
    }

    // for (int i=0; i<nhists; i++) {
    // 	if (histvar[i] == "Mgg") var[i] = Mff;
    // 	else if (histvar[i] == "eta") var[i] = etapar;
    // 	else if (histvar[i] == "dphi") var[i] = PHIgg;
    // 	else if (histvar[i] == "pTgg") var[i] = pT12;
    // 	else if (histvar[i] == "pT") var[i] = pT;
    // 	else if (histvar[i] == "CosP") var[i] = CosTHETAStar;
    // 	else if (histvar[i] == "CosA") var[i] = CosTHETAStarApp;
    //     else if (histvar[i] == "pT1") var[i] = pT1;
    //     else if (histvar[i] == "pT2") var[i] = pT2;

    // }
}
