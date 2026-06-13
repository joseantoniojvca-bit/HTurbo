#include "params.h"
#include "Inter_params.h"

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>


// output options
bool SaveRootHistos = true;
bool SaveRootNtuple = false;
bool SaveTopdHistos = true;

// kinematic cuts
double ph1_ptmin = 25.0; // in GeV
double ph2_ptmin = 22.0; // in GeV
double ph_eta_max = 2.37; 
double crack_eta_min = 1.37;
double crack_eta_max = 1.52;
double ph_deltaR_min = 0.4;

// inv. mass cuts. 
// mgg_min should be *at least* 1 GeV. 
// mgg_max can be even larger than 1TeV but calculation time will be increased
// Higgs
// const  double mgg_min = 100.;
// const  double mgg_max = 160.;
// SM diphoton
double mgg_min = 1.;
double mgg_max = 1000.;
double muR = 1.0;
double muF = 1.0;
//double Scale = Mff;
std::string Scale= "Mff";


// isolation cuts
double ph_iso_R0 = 0.4; 
double ph_iso_E0 = 4.0; // in GeV

// cm energy
double ecm = 7000.; // in GeV

//
// 2gNNLO options
//
// contributions to calculate
int _iifborn = 1;
int _iifnnlo = 1;
int _iinnlo  = 1;
//int _iilean  = 1;
int _iihyylo = 1;
int _iirhyy  = 1; 
// run options
int _i0   = 1;
int _ifnn = 1;
int _inn  = 1;
// iterations
// 1st time
int _n0    = -10; // BORN
int _fnn0  = -20; // NLO
int _nn0   = -50;// NNLO
// 2nd time (histos, events)
/*
   int _n0    = -20; // BORN
   int _fnn0  = -50; // NLO
   int _nn0   = -100;// NNLO
*/
// number of calls
// 1st time (only total xsection)
int _ncl0  = 200000;  // BORN
int _fnncl = 500000;  // NLO
int _nncl  = 3000000; // NNLO

// 2nd time (histos, events)
/*
 * 
  int _ncl0  = 200000; // BORN
  int _fnncl = 500000; // NLO
  int _nncl  = 10000000; // NNLO
  * 
  * 
*/

// Interference
std::string Contrib_Inter = "SigInt";
double _VEV_= 246.221;
double _MH_= 125.0;        // To adjust to most recent value
double _GammaH_= 0.0042;  // To adjust to most recent value

//#define mt 172.5         // To adjust to most recent value
double mt = 168.2;         // Paper Martin
double mb = 4.5;// To adjust to most recent value
double mc = 1.4;

double mW = 80.399;        // To adjust to most recent value

double mtau = 1.77684;    // To adjust to most recent value   
//#define mtau 1.744     // Paper Martin

double mZ = 91.1876;
int withsmearing=0;
double sin2W = 0.23116;
double _GammaW_= 2.085;
double _GammaZ_= 2.4952;
double SigmaSmearing = 1.7;
//double GF=;

void ReadParams(const char* filename) {
    
    std::ifstream infile(filename);
    if (!infile.good()) {
	std::cout << "Input file " << filename << " not found" << std::endl;
	exit(-1);
    }
    while (!infile.eof()) {
	char buf[500];
	infile.getline(buf,500,'\n');
	if (buf[0]=='#') continue;
	std::stringstream str(buf);
	std::string token;
	str >> token;
	if( !infile.good()  ) break; 
	
	if      (token=="ph1_ptmin") str >> ph1_ptmin;
	else if (token=="ph2_ptmin") str >> ph2_ptmin;
	else if (token=="ph_eta_max") str >> ph_eta_max;
	else if (token=="crack_eta_min") str >> crack_eta_min;
	else if (token=="crack_eta_max") str >> crack_eta_max;
	else if (token=="ph_deltaR_min") str >> ph_deltaR_min;
	else if (token=="mgg_min") str >> mgg_min;
	else if (token=="mgg_max") str >> mgg_max;
	else if (token=="muR") str >>muR;
    else if (token=="muF") str >>muF;
    else if (token=="Scale") str >>Scale;
	else if (token=="ph_iso_R0") str >> ph_iso_R0;
	else if (token=="ph_iso_E0") str >> ph_iso_E0;
	else if (token=="ecm") str >> ecm;
	else if (token=="SaveRootHistos") str >> SaveRootHistos;
	else if (token=="SaveRootNtuple") str >> SaveRootNtuple;
	else if (token=="SaveTopdHistos") str >> SaveTopdHistos;
	
	else if (token=="iifborn") str >> _iifborn;
	else if (token=="iifnnlo") str >> _iifnnlo;
	else if (token=="iinnlo") str >> _iinnlo;
	
    else if (token=="iihyylo") str >> _iihyylo;
	else if (token=="iirhyy") str >> _iirhyy;
	
	else if (token=="Contrib_Inter") str >> Contrib_Inter;
	
	else if (token=="VEV") str >> _VEV_;
	else if (token=="MH") str >> _MH_;
	else if (token=="GammaH") str >> _GammaH_;
	else if (token=="GammaW") str >> _GammaW_;
	else if (token=="GammaZ") str >> _GammaZ_;
	else if (token=="mt") str >> mt;
	else if (token=="mb") str >> mb;
	else if (token=="mc") str >> mc;
	else if (token=="mW") str >> mW;
	else if (token=="mZ") str >> mZ;
	else if (token=="mtau") str >> mtau;
	else if (token=="sin2W") str >> sin2W;
	else if (token=="withsmearing") str >> withsmearing;
	else if (token=="SigmaSmearing") str >> SigmaSmearing;
//	else if (token=="GF") str >> GF;
	
	
	
	else if (token=="i0") str >> _i0;
	else if (token=="ifnn") str >> _ifnn;
	else if (token=="inn") str >> _inn;
	
	else if (token=="n0") str >> _n0;
	else if (token=="fnn0") str >> _fnn0;
	else if (token=="nn0") str >> _nn0;
	
	else if (token=="ncl0") str >> _ncl0;
	else if (token=="fnncl") str >> _fnncl;
	else if (token=="nncl") str >> _nncl;
    }
}

std::string OnOff(int i) {
    return (i==0) ? "OFF" : "ON";
}

std::string Option(int i) {
    if (i==1) return "RESTART";
    else if (i==2) return "CONTINUE";
    else if (i==3) return "RESTART W/ PREVIOUS GRIDS";
    else return "UNKNOWN";
}

void DumpParams() {
    std::cout << std::endl;
    std::cout << "Settings:" << std::endl;
    std::cout << std::endl;
    std::cout << "Contributions to run and options:" << std::endl;
    std::cout << "BORN: " << OnOff(_iifborn) << " " << Option(_i0) << " " << _n0 << " iterations " << _ncl0 << " calls" << std::endl;
    std::cout << "NLO : " << OnOff(_iifnnlo) << " " << Option(_ifnn) << " " << _fnn0 << " iterations " << _fnncl << " calls" << std::endl;
    std::cout << "NNLO: " << OnOff(_iinnlo) << " " << Option(_inn) << " " << _nn0 << " iterations " << _nncl << " calls" << std::endl;
    
    
    std::cout << "Interference LO:  " << OnOff(_iihyylo) << " " << Option(_ifnn) << " " << _fnn0 << " iterations " <<_fnncl << " calls" << std::endl;
    std::cout << "Interference Reals:  " << OnOff(_iirhyy) << " " << Option(_ifnn) << " " <<_fnn0 << " iterations " <<_fnncl << " calls" << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "Interference type:  " << OnOff(1) << " " << Contrib_Inter << std::endl;
    
     
    std::cout << std::endl;
    std::cout << "Interference values" << std::endl;
    std::cout << "VEV : " <<  _VEV_ << " GeV" << std::endl;
    std::cout << "MH : " <<  _MH_ << " GeV" << std::endl;
    std::cout << "mt : " << mt << " GeV"<< std::endl;
    std::cout << "mb : " << mb << " GeV"<< std::endl;
    std::cout << "mc : " << mc << " GeV"<< std::endl;
    std::cout << "mW  : " << mW << " GeV"<< std::endl;
    std::cout << "mtau : " << mtau << " GeV"<< std::endl;
    
    
    
    std::cout << std::endl;
    std::cout << "Acceptance cuts" << std::endl;
    std::cout << "ph1_ptmin : " << ph1_ptmin << " GeV" << std::endl;
    std::cout << "ph2_ptmin : " << ph2_ptmin << " GeV" << std::endl;
    std::cout << "ph_eta_max : " << ph_eta_max << std::endl;
    std::cout << "crack_eta_min : " << crack_eta_min << std::endl;
    std::cout << "crack_eta_max : " << crack_eta_max << std::endl;
    std::cout << "ph_deltaR_min : " << ph_deltaR_min << std::endl;
    std::cout << "mgg_min : " << mgg_min << " GeV" << std::endl;
    std::cout << "mgg_max : " << mgg_max << " GeV" << std::endl;
    
    std::cout << "Scale:  " << " " << Scale << std::endl;
    std::cout << "muR : " << muR << std::endl;
    std::cout << "muF : " << muF << std::endl;
    
    std::cout << "ph_iso_R0 : " << ph_iso_R0 << std::endl;
    std::cout << "ph_iso_E0 : " << ph_iso_E0 << " GeV" << std::endl;
    std::cout << "ecm : " << ecm << " GeV" << std::endl;
    std::cout << std::endl;
}
